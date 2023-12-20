#include <iostream>
#include <string>
#include <chrono>
#include <array>

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
	const double energy = 7000.;
	const double ptmin = 25.;
	const double nevents = 1e4;
	const double abs_max_y = 4.7;
	
	const double fastjet_r_par = 0.4;
	fastjet::Strategy strategy = fastjet::Best;

	std::string pdf_set = "LHAPDF6:NNPDF31_lo_as_0118";

	//neutrinos id set to exclude from the jet algorithm
	std::set<int> exclude_id = {12, 14, 16, 18};
} Par;

struct FastJetVector
{
	std::vector<fastjet::PseudoJet> input, inclusive;
};

unsigned int GetRandomSeed()
{
	auto now = std::chrono::high_resolution_clock::now();
	auto now_ms = std::chrono::time_point_cast<std::chrono::milliseconds>(now);
	auto epoch = now_ms.time_since_epoch();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(epoch);
	return static_cast<unsigned int>(duration.count()) % 900000000;
}

void PrintParameters(unsigned int seed)
{
	Box box("Parameters");
	box.AddEntry("CM beams energy, TeV", Par.energy/1e3, 3);
	box.AddEntry("PDF set", Par.pdf_set);
	box.AddEntry("Minimum pT, GeV", Par.ptmin, 3);
	box.AddEntry("|ymax|", Par.abs_max_y, 3);
	box.AddEntry("seed", seed);
	box.Print();
}

bool IsExcludedPart(int id)
{
	id = abs(id);
	if (auto search = Par.exclude_id.find(id); search != Par.exclude_id.end()) return true;
	return false;
}

bool IsParton(int id)
{
	if (abs(id) < 9 || id == 21) return true;
	return false;
}

int main()
{
	Pythia pythia;
	unsigned int seed = GetRandomSeed();
	
	//setting pythia parameters
	pythia.readString("Beams:eCM = " + to_string(Par.energy));
	pythia.readString("HardQCD:all = on");
	pythia.readString("PDF:pSet = " + Par.pdf_set);
	pythia.readString("PhaseSpace::pTHatMin = " + to_string(Par.ptmin));
	
	pythia.readString("Random:setSeed = on");
	pythia.readString("Random:seed = " + to_string(seed));
	pythia.readString("Print:quiet = on");
	
	//creating a structure of vectors for fastjet
	
	//setting fastjet parameters
	fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, Par.fastjet_r_par, Par.strategy);

	//initializing pythia
	pythia.init();

	//creating directory for output
	system("mkdir ../output");
	
	//printing parameters info
	PrintParameters(seed);

	//partons multiplicity vs pt
	TH1D hist_part_pt = TH1D("part_mult_pt", "dsigma/dpt", 200, 0., 200.);
	//pair of partons multiplicity vs delta y
	TH1D hist_part_dy = TH1D("part_mult_dy", "dsigma/ddy", 
		200, 0., static_cast<double>(ceil(Par.abs_max_y*2.)));

	//jets multiplicity vs pt
	TH1D hist_jet_pt = TH1D("jet_mult_pt", "dsigma/dpt", 200, 0., 200.);
	//pair of jets multiplicity vs delta y
	TH1D hist_jet_dy = TH1D("jet_mult_dy", "dsigma/ddy", 
		200, 0., static_cast<double>(ceil(Par.abs_max_y*2.)));

	TH1D gen_info = TH1D("gen_info", "info", 2, 0, 2);
	
	//progress bar
	ProgressBar pbar("FANCY");

	//events loop
	for (long i = 0; i < Par.nevents; i++)
	{
		if (!pythia.next()) continue;
		FastJetVector fjv;
		
		//particles in event loop
		for (int j = 0; j < pythia.event.size(); j++)
		{
			if (!pythia.event[j].isFinal()) 
			{
				if (pythia.event[j].status() != -23) continue;
				if (!IsParton(pythia.event[j].id())) continue;
				if (abs(pythia.event[j].y()) > Par.abs_max_y) continue;
				
				hist_part_pt.Fill(pythia.event[j].pT(), pythia.info.weight());
				if (pythia.event[j].pT() < Par.ptmin) continue;
				
				for (int k = j+1; k < pythia.event.size(); k++)
				{
					if (pythia.event[k].status() != -23) continue;
					if (!IsParton(pythia.event[k].id())) continue;
					if (pythia.event[k].pT() < Par.ptmin) continue;
					if (abs(pythia.event[k].y()) > Par.abs_max_y) continue;
					
					const double delta_y = abs(pythia.event[j].y() - pythia.event[k].y());
					hist_part_dy.Fill(delta_y, pythia.info.weight());
				}

			}
			else
			{
				fjv.input.push_back(fastjet::PseudoJet(
					pythia.event[j].px(), 
					pythia.event[j].py(), 
					pythia.event[j].pz(), 
					pythia.event[j].e()));
			}
		}

		fastjet::ClusterSequence cluster_seq(fjv.input, jet_def);
		fjv.inclusive = cluster_seq.inclusive_jets(Par.ptmin);
		
		pbar.Print(static_cast<double>(i)/Par.nevents);
		
		//jets loop
		for (int j = 0; j < fjv.inclusive.size(); j++)
		{
			if (abs(fjv.inclusive[j].rap()) > Par.abs_max_y) continue;
			hist_jet_pt.Fill(fjv.inclusive[j].pt(), pythia.info.weight());
			
			if (fjv.inclusive[j].pt() < Par.ptmin) continue;
			
			//loop to form pairs of jets
			for (int k = j+1; k < fjv.inclusive.size(); k++)
			{
				if (fjv.inclusive[k].pt() < Par.ptmin) continue;
				if (abs(fjv.inclusive[k].rap()) > Par.abs_max_y) continue;
				
				const double delta_y = abs(fjv.inclusive[j].rap() - fjv.inclusive[k].rap());
				hist_jet_dy.Fill(delta_y, pythia.info.weight());
			}
		}
	}
	pbar.Print(1);

	gen_info.SetBinContent(1, pythia.info.nAccepted());
	gen_info.SetBinContent(2, pythia.info.sigmaGen()*1.e9);
	
	std::string output_file_name = "../output/gen.root";
	TFile output = TFile(output_file_name.c_str(), "RECREATE");

	//cross sections for the quick access in TFile
	TH1D *hist_part_dsigma_dpt = (TH1D *) hist_part_pt.Clone("part_dsigma_dpt");
	TH1D *hist_part_dsigma_ddy = (TH1D *) hist_part_dy.Clone("part_dsigma_ddy");
	TH1D *hist_jet_dsigma_dpt = (TH1D *) hist_jet_pt.Clone("jet_dsigma_dpt");
	TH1D *hist_jet_dsigma_ddy = (TH1D *) hist_jet_dy.Clone("jet_dsigma_ddy");

	hist_part_dsigma_dpt->Scale(pythia.info.sigmaGen()*1e9/
		(pythia.info.nAccepted()*hist_part_dsigma_dpt->GetXaxis()->GetBinWidth(1)));
	hist_part_dsigma_ddy->Scale(pythia.info.sigmaGen()*1e9/
		(pythia.info.nAccepted()*hist_part_dsigma_ddy->GetXaxis()->GetBinWidth(1)));
	hist_jet_dsigma_dpt->Scale(pythia.info.sigmaGen()*1e9/
		(pythia.info.nAccepted()*hist_jet_dsigma_dpt->GetXaxis()->GetBinWidth(1)));
	hist_jet_dsigma_ddy->Scale(pythia.info.sigmaGen()*1e9/
		(pythia.info.nAccepted()*hist_jet_dsigma_ddy->GetXaxis()->GetBinWidth(1)));
	
	gen_info.Write();
	
	hist_part_pt.Write();
	hist_part_dy.Write();
	hist_jet_pt.Write();
	hist_jet_dy.Write();
	
	hist_part_dsigma_dpt->Write();
	hist_part_dsigma_ddy->Write();
	hist_jet_dsigma_dpt->Write();
	hist_jet_dsigma_ddy->Write();
	
	output.Close();
	
	PrintInfo("File " + output_file_name + " was written");
	
	return 0;
}
