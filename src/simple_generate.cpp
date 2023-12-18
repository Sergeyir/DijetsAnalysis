#include <iostream>/
#include <string>
#include <chrono>

#include "Pythia8/Pythia.h"

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
	std::string input_file = "../input/simple.cmnd";
} Par;

void PrintParameters(unsigned long seed)
{
	Box box("Parameters");
	box.AddEntry("CM beams energy, TeV", Par.energy/1e3, 3);
	box.AddEntry("Minimum pT, GeV", Par.ptmin, 3);
	box.AddEntry("|ymax|", Par.abs_max_y, 3);
	box.AddEntry("seed", seed);
	box.Print();

	Print("Pythia parameters imported from file", Par.input_file);
	system(("cat ../input/" + Par.input_file).c_str());
}

unsigned long GetRandomSeed()
{
	auto now = std::chrono::high_resolution_clock::now();
	auto now_ms = std::chrono::time_point_cast<std::chrono::milliseconds>(now);
	auto epoch = now_ms.time_since_epoch();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(epoch);
	return static_cast<unsigned long>(duration.count()) % 900000000;
}


bool IsParton(int id)
{
	if (abs(id) < 9 || id == 21) return true;
	return false;
}

void MakeSpectra(TH1D *mult_hist)
{
	for (int i = 1; i < mult_hist->GetXaxis()->GetNbins(); i++)
	{
		const double scale = 2.*M_PI*mult_hist->GetXaxis()->GetBinCenter(i)*mult_hist->GetXaxis()->GetBinWidth(i);
		mult_hist->SetBinContent(i, mult_hist->GetBinContent(i)/scale);
		mult_hist->SetBinError(i, mult_hist->GetBinError(i)/scale);
	}
}

int main()
{
	Pythia pythia;
	unsigned long seed = GetRandomSeed();
	
	//setting pythia parameters
	pythia.readString("Beams:eCM = " + to_string(Par.energy));
	pythia.readFile(Par.input_file.c_str());
	pythia.readString("PhaseSpace::pTHatMin = " + to_string(Par.ptmin));
	
	pythia.readString("Random:setSeed = " + to_string(seed));
	pythia.readString("Random:seed = " + to_string(seed));
	pythia.readString("Print:quiet = on");
	
	//initializing pythia
	pythia.init();

	//creating directory for output
	system("mkdir ../output");
	
	//printing parameters info
	PrintParameters(seed);
	
	//parton multiplicities
	TH1D hist_mult_pt = TH1D("part_mult_pt", "part", 200, 0., 200.);
	
	//parton spectra
	TH1D hist_dsigma_ddy = TH1D("dsigma_ddy", "dsigma_ddy", 
		100, 0., static_cast<double>(ceil(Par.abs_max_y*2.)));

	//progress bar
	ProgressBar pbar("FANCY");

	//sum of weighted events
	double sum_weight = 0.;
	
	//events loop
	for (unsigned long i = 0; i < Par.nevents; i++)
	{
		pbar.Print(static_cast<double>(i)/Par.nevents);
		if (!pythia.next()) continue;
		sum_weight += pythia.info.weight();
			
		//particles event loop
		for (int j = 0; j < pythia.event.size(); j++)
		{
			//only outcoming final partons
			if (!pythia.event[j].isFinal()) continue;
			if (abs(pythia.event[j].status()) != 62) continue;
			if (!IsParton(pythia.event[j].id())) continue;
			if (abs(pythia.event[j].y()) > Par.abs_max_y) continue;
			
			hist_mult_pt.Fill(pythia.event[j].pT(), pythia.info.weight());

			if (pythia.event[j].pT() > Par.ptmin) continue;
			
			for (int k = j+1; k < pythia.event.size(); k++)
			{
				if (!pythia.event[k].isFinal()) continue;
				if (abs(pythia.event[k].status()) != 62) continue;
				if (!IsParton(pythia.event[k].id())) continue;
				if (abs(pythia.event[k].y()) > Par.abs_max_y) continue;

				//pythia clones particles for the next step that was disabled
				//checking mothers of the cloned particles to get dijets
				if (pythia.event[pythia.event[j].mother1()].mother1() != 
					pythia.event[pythia.event[k].mother1()].mother1() || 
					pythia.event[pythia.event[j].mother2()].mother1() != 
					pythia.event[pythia.event[k].mother2()].mother1()) continue;
				
				if (pythia.event[k].pT() > Par.ptmin) continue;
				
				const double delta_y = abs(pythia.event[j].y() - pythia.event[k].y());
				hist_dsigma_ddy.Fill(delta_y, pythia.info.weight());
			}
		}
	}
	pbar.Print(1);
	
	const double sigma_pb = pythia.info.sigmaGen()*1e9;
	
	std::string output_file_name = "../output/gen_simple.root";
	TFile output = TFile(output_file_name.c_str(), "RECREATE");
	
	hist_mult_pt.Scale(sigma_pb/sum_weight);
	hist_dsigma_ddy.Scale(sigma_pb/(2.*M_PI*hist_dsigma_ddy.GetXaxis()->GetBinWidth(1)*sum_weight));
	
	TH1D *hist_dsigma_dpt = dynamic_cast<TH1D *>(hist_mult_pt.Clone());
	hist_dsigma_dpt->SetName("dsigma_dpt");
	
	MakeSpectra(hist_dsigma_dpt);
	
	hist_dsigma_dpt->Write();
	hist_dsigma_ddy.Write();
	
	hist_mult_pt.Write();
	
	output.Close();
	PrintInfo("File " + output_file_name + " was written");
	
	return 0;
}