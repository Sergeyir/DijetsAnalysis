#include <iostream>
#include <string>
#include <chrono>

#include "Pythia8/Pythia.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include "TFile.h"
#include "TH1D.h"

#include "../lib/Box.h"
#include "../lib/ProgressBar.h"
#include "../lib/InputTool.h"

using namespace Pythia8;

struct
{
	double energy;
	double pt_min;
	std::string pdf_set;
	double nevents;
	double abs_max_y;
	std::string hadronization, isr, fsr;

	double fastjet_r_par;
	fastjet::Strategy strategy = fastjet::Best;

	std::string input_file;

	//leptons and gauge bozons id set to exclude from the jet algorithm
	std::set<int> exclude_id = {12, 14, 16, 18};
} Par;

struct FastJetVector
{
	std::vector<fastjet::PseudoJet> input, inclusive;

	void Clear()
	{
		input.clear();
		inclusive.clear();
	}
};

unsigned long GetRandomSeed()
{
	auto now = std::chrono::high_resolution_clock::now();
	auto now_ms = std::chrono::time_point_cast<std::chrono::milliseconds>(now);
	auto epoch = now_ms.time_since_epoch();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(epoch);
	return static_cast<unsigned long>(duration.count()) % 900000000;
}

void SetParameters()
{
	std::array<std::string, 9> expected_name = 
	{
		"Energy", 
		"pdf_set", 
		"pT_min", 
		"Events", 
		"abs_max_y", 
		"Hadronization", 
		"ISR", "FSR", 
		"FastJet_R"
	};
	
	std::array<std::string, 9> name;
	ReadFile("../input/" + Par.input_file + ".cmnd", 
		name[0], Par.energy,
		name[1], Par.pdf_set,
		name[2], Par.pt_min,
		name[3], Par.nevents,
		name[4], Par.abs_max_y,
		name[5], Par.hadronization,
		name[6], Par.isr,
		name[7], Par.fsr,
		name[8], Par.fastjet_r_par);

	for (int i = 0; i < 9; i++)
	{
		if (name[i] != expected_name[i]) PrintError("Input file parameter name mismatch: " + 
			name[0] + " vs " + expected_name[i]);
	}
}

void PrintParameters(unsigned long seed)
{
	Box box("Parameters");
	
	box.AddEntry("CM beams energy, TeV", Par.energy/1e3, 3);
	box.AddEntry("Minimum pT, GeV", Par.pt_min, 3);
	box.AddEntry("Pdf set", Par.pdf_set);

	box.AddEntry("|ymax|", Par.abs_max_y, 3);
	
	box.AddEntry("Hadronization", Par.hadronization);
	box.AddEntry("Initial state radiation", Par.isr);
	box.AddEntry("Final state radiation", Par.fsr);

	box.AddEntry("Number of events, 1e3", Par.nevents/1e3);
	box.AddEntry("Seed", seed);
	box.AddEntry("FastJet: R", Par.fastjet_r_par);
	
	box.Print();
}

bool IsExcludedPart(int id)
{
	id = abs(id);
	if (auto search = Par.exclude_id.find(id); search != Par.exclude_id.end()) return true;
	return false;
}

void MakeSpectra(TH1D *mult_hist)
{
	for (int i = 1; i < mult_hist->GetXaxis()->GetNbins(); i++)
	{
		const double scale = 2.*3.14159265359*mult_hist->GetXaxis()->GetBinCenter(i)*mult_hist->GetXaxis()->GetBinWidth(i);
		mult_hist->SetBinContent(i, mult_hist->GetBinContent(i)/scale);
		mult_hist->SetBinError(i, mult_hist->GetBinError(i)/scale);
	}
}

int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		PrintError("Input file name was not passed!", 0);
		Print("Input files specified in your input directory (.cmnd)");
		system("cd ../input && ls *cmnd | sed 's/\.cmnd$//'");
		PrintInfo("Pass the input file name without .cmnd extension");
		exit(0);
	}
	else if (argc > 2) PrintWarning("Only the first argument is passed. Others are allways ignored");
	
	Par.input_file = static_cast<std::string>(argv[1]);
	SetParameters();
	
	Pythia pythia;
	unsigned long seed = GetRandomSeed();
	
	//setting pythia parameters
	pythia.readString("Beams:eCM = " + to_string(Par.energy));
	pythia.readString("PDF:pSet = LHAPDF6:" + Par.pdf_set);
	pythia.readString("HardQCD:all = on");
	pythia.readString("PhaseSpace::pTHatMin = " + to_string(Par.pt_min));
	
	pythia.readString("Random:setSeed = on");
	pythia.readString("Random:seed = " + to_string(seed));

	pythia.readString("Print:quiet = on");

	pythia.readString("HadronLevel:Hadronize = " + Par.hadronization);
	pythia.readString("PartonLevel:ISR = " + Par.isr);
	pythia.readString("PartonLevel:FSR = " + Par.fsr);

	//creating structure of vectors for fastjet
	FastJetVector fjv;
	
	//setting fastjet parameters
	fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, Par.fastjet_r_par, Par.strategy);

	//initializing pythia
	pythia.init();

	//creating directory for output
	system("mkdir ../output");
	
	//printing parameters info
	PrintParameters(seed);

	//jets multiplicity vs pt
	TH1D hist_njets = TH1D("jets_multiplicity", "jets", 200, 0., 200.);
	//pair of jets multiplicity vs y
	TH1D hist_jet_pairs = TH1D("jets_pairs", "jets", 100, 
		0., static_cast<double>(ceil(Par.abs_max_y*2.)));

	//progress bar
	ProgressBar pbar("FANCY");

	//sum of weighted events
	double sum_weight = 0.;
	
	//events loop
	for (unsigned long i = 0; i < Par.nevents; i++)
	{
		if (!pythia.next()) continue;

		sum_weight += pythia.info.weight();
			
		//particles in event loop
		for (int j = 0; j < pythia.event.size(); j++)
		{
			if (!pythia.event[j].isFinal()) continue;
			
			//no neutrino or photons
			if (IsExcludedPart(pythia.event[j].id())) continue;

			fjv.input.push_back(fastjet::PseudoJet(
				pythia.event[j].px(), 
				pythia.event[j].py(), 
				pythia.event[j].pz(), 
				pythia.event[j].e()));
		}

		fastjet::ClusterSequence cluster_seq(fjv.input, jet_def);
		fjv.inclusive = cluster_seq.inclusive_jets();
		
		//Printing progress in a progressbar
		pbar.Print(static_cast<double>(i)/Par.nevents);
		
		//jets loop
		for (int j = 0; j < fjv.inclusive.size(); j++)
		{
			if (abs(fjv.inclusive[j].rap()) > Par.abs_max_y) continue;
			//loop to form pairs of jets
			
			hist_njets.Fill(fjv.inclusive[j].pt(), pythia.info.weight());
			
			if (fjv.inclusive[j].pt() < Par.pt_min) continue;
			
			for (int k = j+1; k < fjv.inclusive.size(); k++)
			{
				if (fjv.inclusive[k].pt() < Par.pt_min) continue;
				if (abs(fjv.inclusive[k].rap()) > Par.abs_max_y) continue;
				
				const double delta_y = abs(fjv.inclusive[j].rap() - fjv.inclusive[k].rap());
				
				hist_jet_pairs.Fill(delta_y, pythia.info.weight());
			}
		}
		
		//clearing vectors for the next event
		fjv.Clear();
	}
	
	const double sigma_pb = pythia.info.sigmaGen()*1e9;

	pbar.Print(1);

	std::string output_file_name = "../output/jets_" + Par.input_file + ".root";
	TFile output = TFile(output_file_name.c_str(), "RECREATE");
	
	hist_njets.Scale(sigma_pb/sum_weight);
	hist_jet_pairs.Scale(sigma_pb/sum_weight);

	TH1D *dsigma_dpt = dynamic_cast<TH1D *>(hist_njets.Clone());
	TH1D *dsigma_ddy = dynamic_cast<TH1D *>(hist_jet_pairs.Clone());
	
	dsigma_dpt->SetName("dsigma_dpt");
	dsigma_ddy->SetName("dsigma_ddy");
	
	MakeSpectra(dsigma_dpt);
	MakeSpectra(dsigma_ddy);
	
	dsigma_dpt->Write();
	dsigma_ddy->Write();
	hist_njets.Write();
	hist_jet_pairs.Write();

	output.Close();

	PrintInfo("File " + output_file_name + " was written");
	
	return 0;
}
