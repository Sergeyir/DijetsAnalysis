#include <iostream>/
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
	const double energy = 7000.;
	const double ptmin = 25.;
	const double nevents = 1e5;
	const double abs_max_y = 4.7;
	
	const double fastjet_r_par = 0.4;
	fastjet::Strategy strategy = fastjet::Best;

	std::string input_file = "../input/default.cmnd";

	//neutrinos id set to exclude from the jet algorithm
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

void PrintParameters(unsigned long seed)
{
	Box box("Parameters");
	box.AddEntry("CM beams energy, TeV", Par.energy/1e3, 3);
	box.AddEntry("Minimum pT, GeV", Par.ptmin, 3);
	box.AddEntry("|ymax|", Par.abs_max_y, 3);
	box.Print();

	Print("Pythia parameters imported from file", Par.input_file);
	system(("cat " + Par.input_file).c_str());
}

bool IsExcludedPart(int id)
{
	id = abs(id);
	if (auto search = Par.exclude_id.find(id); search != Par.exclude_id.end()) return true;
	return false;
}

int main()
{
	Pythia pythia;
	unsigned long seed = GetRandomSeed();
	
	//setting pythia parameters
	pythia.readString("Beams:eCM = " + to_string(Par.energy));
	pythia.readFile(Par.input_file);
	pythia.readString("PhaseSpace::pTHatMin = " + to_string(Par.ptmin));
	
	pythia.readString("Random:setSeed = on");
	pythia.readString("Random:seed = " + to_string(seed));
	pythia.readString("Print:quiet = on");
	
	//creating a structure of vectors for fastjet
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
	TH1D hist_dsigma_dpt = TH1D("dsigma_dpt", "dsigma/dpt", 40, 0., 200.);
	//pair of jets multiplicity vs y
	TH1D hist_dsigma_ddy = TH1D("dsigma_ddy", "dsigma/ddy", 
		40, 0., static_cast<double>(ceil(Par.abs_max_y*2.)));
	
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
			if (IsExcludedPart(pythia.event[j].id())) continue;
			
			fjv.input.push_back(fastjet::PseudoJet(
				pythia.event[j].px(), 
				pythia.event[j].py(), 
				pythia.event[j].pz(), 
				pythia.event[j].e()));
		}

		fastjet::ClusterSequence cluster_seq(fjv.input, jet_def);
		fjv.inclusive = cluster_seq.inclusive_jets();
		
		pbar.Print(static_cast<double>(i)/Par.nevents);
		
		//jets loop
		for (int j = 0; j < fjv.inclusive.size(); j++)
		{
			if (abs(fjv.inclusive[j].rap()) > Par.abs_max_y) continue;
			hist_dsigma_dpt.Fill(fjv.inclusive[j].pt(), pythia.info.weight());
			
			if (fjv.inclusive[j].pt() < Par.ptmin) continue;
			
			//loop to form pairs of jets
			for (int k = j+1; k < fjv.inclusive.size(); k++)
			{
				if (fjv.inclusive[k].pt() < Par.ptmin) continue;
				if (abs(fjv.inclusive[k].rap()) > Par.abs_max_y) continue;
				
				const double delta_y = abs(fjv.inclusive[j].rap() - fjv.inclusive[k].rap());
				hist_dsigma_ddy.Fill(delta_y, pythia.info.weight());
			}
		}
		//clearing vectors for the next event
		fjv.Clear();
	}
	pbar.Print(1);
	
	const double sigma_pb = pythia.info.sigmaGen()*1e9;
	
	std::string output_file_name = "../output/gen_default.root";
	TFile output = TFile(output_file_name.c_str(), "RECREATE");
	
	hist_dsigma_dpt.Scale(sigma_pb/(sum_weight*hist_dsigma_dpt.GetXaxis()->GetBinWidth(1)));
	hist_dsigma_ddy.Scale(sigma_pb*2.*pow(M_PI,2)/(sum_weight*hist_dsigma_ddy.GetXaxis()->GetBinWidth(1)));
	
	hist_dsigma_dpt.Write();
	hist_dsigma_ddy.Write();

	output.Close();
	
	PrintInfo("File " + output_file_name + " was written");
	
	return 0;
}
