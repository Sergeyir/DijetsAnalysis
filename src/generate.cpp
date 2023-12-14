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

using namespace Pythia8;

struct
{
	const double energy = 7000;
	const double pt_min = 25;
	std::string pdf_set = "NNPDF31_lo_as_0118";
	const double nevents = 1e5;
	const double abs_max_y = 4.7;

	const double fastjet_r_par = 0.4;
	fastjet::Strategy strategy = fastjet::Best;
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

void PrintParameters(const int setup, unsigned long seed)
{
	Box box("Parameters");
	
	box.AddEntry("CM beams energy, TeV", Par.energy/1e3, 3);
	box.AddEntry("Minimum pT, GeV", Par.pt_min, 3);
	box.AddEntry("Pdf set", Par.pdf_set);

	box.AddEntry("|ymax|", Par.abs_max_y, 3);
	
	if (setup == 1)
	{
		box.AddEntry("Hadronization", "off");
		box.AddEntry("Initial state radiation", "off");
		box.AddEntry("Final state radiation", "off");
	}

	box.AddEntry("Number of events, 1e3", Par.nevents/1e3);
	box.AddEntry("Seed", seed);
	box.AddEntry("FastJet: R", Par.fastjet_r_par);
	
	box.Print();
}

bool IsNeutrinoOrPhotonId(const int id)
{
	if (id == 22 || id == 12 || id == 14 || id == 16 || id == 18) return true;
	return false;
}

int main(int argc, char *argv[])
{
	if (argc = 0)
	{
		PrintError("Parameters were not passed");
	}
	
	int setup = atoi(argv[1]);
	if (setup != 0 && setup != 1)
	{
		PrintError("Unknown setup: " + to_string(setup));
	}
	
	Pythia pythia;
	unsigned long seed = GetRandomSeed();
	
	//setting pythia parameters
	pythia.readString("Beams:eCM = " + to_string(Par.energy));
	//pythia.readString("PhaseSpace::pTHatMin = " + to_string(Par.pt_min));
	pythia.readString("PDF:pSet = LHAPDF6:" + Par.pdf_set);
	pythia.readString("HardQCD:all = on");
	
	pythia.readString("Random:setSeed = on");
	pythia.readString("Random:seed = " + to_string(seed));

	pythia.readString("Print:quiet = on");

	if (setup == 1)
	{
		pythia.readString("PartonLevel:ISR = off");
		pythia.readString("PartonLevel:FSR = off");
		pythia.readString("HadronLevel:Hadronize = off");
	}

	//creating structure of vectors for fastjet
	FastJetVector fjv;
	
	//setting fastjet parameters
	fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, Par.fastjet_r_par, Par.strategy);

	//initializing pythia
	pythia.init();

	//creating directory for output
	system("mkdir ../output");
	
	//printing parameters info
	PrintParameters(setup, seed);

	//azimuthal angle hist of pair of jets
	TH1D hist_theta = TH1D("theta", "theta", 64, -1.6, 1.6);
	//jets multiplicity vs pt
	TH1D hist_njets = TH1D("jets_multiplicity", "jets", 100, 0, 100);
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
			if (IsNeutrinoOrPhotonId(pythia.event[j].id())) continue;

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
			hist_njets.Fill(fjv.inclusive[j].pt(), pythia.info.weight());

			if (abs(fjv.inclusive[j].rap()) > Par.abs_max_y) continue;
			
			//loop to form pairs of jets
			for (int k = j+1; k < fjv.inclusive.size(); k++)
			{
				if (abs(fjv.inclusive[k].rap()) > Par.abs_max_y) continue;
				
				const double delta_y = abs(fjv.inclusive[j].rap() - fjv.inclusive[k].rap());
				
				/*
				const double theta = atan(
					(fjv.inclusive[j].py() + fjv.inclusive[k].py())/
					(fjv.inclusive[j].px() + fjv.inclusive[k].px()));
				hist_theta.Fill(theta, pythia.info.weight());
				*/
				
				hist_jet_pairs.Fill(delta_y, pythia.info.weight());
			}
		}
		
		//clearing vectors for the next event
		fjv.Clear();
	}

	pbar.Print(1);

	std::string output_file_name = "../output/jets" + to_string(setup) + ".root";
	TFile output = TFile(output_file_name.c_str(), "RECREATE");

	//normalizing hists by the sum of weighted events
	hist_theta.Scale(1./sum_weight);
	hist_njets.Scale(1./sum_weight);
	hist_jet_pairs.Scale(1./sum_weight);

	hist_theta.Write();
	hist_njets.Write();
	hist_jet_pairs.Write();

	output.Close();

	PrintInfo("File " + output_file_name + " was written");
	
	return 0;
}
