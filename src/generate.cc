#include <iostream>
#include <string>

#include "Pythia8/Pythia.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include "LHAPDF/LHAPDF.h"

using namespace Pythia8;

struct
{
	std::string snn = "7";
	std::string pdf_set = "NNPDF31_lo_as_0118";
} Par;

int main()
{
	Pythia pythia;

	LHAPDF::PDF* pdf = LHAPDF::mkPDF(Par.pdf_set, 0);
	pythia.readString("Beams:eCM = 200");
	//pythia.readString("PDF:Set = " + Par.pdf_set);
	
	return 0;
}
