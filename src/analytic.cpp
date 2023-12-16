#include <iostream>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TH1.h"

#include "LHAPDF/LHAPDF.h"

#include "../lib/ProgressBar.h"
#include "../lib/Box.h"
#include "../lib/Tool.h"

using namespace LHAPDF;
using namespace Tool;

const double pi = 3.14159265359;

struct
{
	const double energy = 7000;
	std::string pdfset_name = "NNPDF31_lo_as_0118";
	const double abs_max_y = 4.7;
	const double ptmin = 25;
	const double nsteps = 10.;

	const double s = energy*energy;
	const PDF *pdf = mkPDF(pdfset_name);
} Par;

//cross sections dsigma/dOmega for different processes
//qq'->qq'
double CS_QQp_QQp(const double s, const double t, const double u)
{
	return 1./(9.*s)*(s*s + u*u)/(t*t);
}

//qq->qq
double CS_QQ_QQ(const double s, const double t, const double u)
{
	return 1./(9.*s)*(((t*t + s*s)/(u*u) + (s*s + u*u)/(t*t)) - 2.*s*s/(3.*u*t));
}

//qqbar->q'qbar'
double CS_QQbar_QpQbarp(const double s, const double t, const double u)
{
	return 1./(9.*s)*(t*t + u*u)/(s*s);
}

//qqbar->qqbar
double CS_QQbar_QQbar(const double s, const double t, const double u)
{
	return 1./(9.*s)*((t*t + u*u)/(s*s) + (s*s + u*u)/(t*t) - 2.*u*u/(3.*s*t));
}

//qqbar->gg
double CS_QQbar_GG(const double s, const double t, const double u)
{
	return 8./(27.*s)*(t*t + u*u)*(1./(t*u) - 9./(4.*s*s));
}

//gg->qqbar
double CS_GG_QQbar(const double s, const double t, const double u)
{
	return 1./(24.*s)*(t*t + u*u)*(1./(t*u) - 9./(4.*s*s));
}

//gq->gq
double CS_GQ_GQ(const double s, const double t, const double u)
{
	return 1./(9.*s)*(s*s + u*u)*(-1./(s*u) + 9./(4.*t*t));
}

//gg->gg
double CS_GG_GG(const double s, const double t, const double u)
{
	return 9./(8.*s)*(3. - (u*t)/(s*s) - (s*u)/(t*t) - (s*t)/(u*u));
}

//variables shortcuts
double CosTheta(const double s, const double pt) {return sqrt(1.-(4.*pt*pt/s));}
double T(const double s, const double cos_theta) {return -s/2.*(1.-cos_theta);}
double U(const double s, const double cos_theta) {return -s/2.*(1.+cos_theta);}

//returns dsigma/dOmega for id1+id2->X+X process
//y = y1 - y2
double CS_XX_XX(const int id1, const int id2, const double s, const double pt, const double y)
{
	double cos_theta = CosTheta(s, pt);
	if (y < 0) cos_theta *= -1.;
	
	const double t = T(s, cos_theta);
	const double u = U(s, cos_theta);

	const double alpha_s = Par.pdf->alphasQ2(pt*pt);
	
	double result = 0.;
	if (id1 == 0 && id2 == 0)
	{
		result = CS_GG_GG(s, t, u) + CS_GG_QQbar(s, t, u);
	}
	else if ((id1 != 0 && id2 == 0) || (id1 == 0 && id2 != 0))
	{
		result = CS_GQ_GQ(s, t, u);
	}
	else if (id1 == id2)
	{
		result = CS_QQ_QQ(s, t, u);
	}
	else if (id1 == -id2)
	{
		result = CS_QQbar_QQbar(s, t, u) + CS_QQbar_GG(s, t, u) + CS_QQbar_QpQbarp(s, t, u);
	}
	else result = CS_QQp_QQp(s, t, u);
	return result*alpha_s*alpha_s;
}

double X1(const double pt, const double sqrt_s, const double y1, const double y2)
{
	return 2.*pt/sqrt_s*exp((y1+y2)/2.)*cosh((y1-y2)/2.);
}

double X2(const double pt, const double sqrt_s, const double y1, const double y2)
{
	return 2.*pt/sqrt_s*exp(-(y1+y2)/2.)*cosh((y1-y2)/2.);
}

double DsigmaDpTDy1Dy2(const double pt, const double s, const double y1, const double y2, const double x1, const double x2)
{
	double result = 0.;
	
	//summation over pairs of flavours
	for (int id1 = -5; id1 <= 5; id1++)
	{
		for (int id2 = -5; id2 <= 5; id2++)
		{
			result += 8.*pi*pt*
				Par.pdf->xfxQ2(id1, x1, pt*pt)*Par.pdf->xfxQ2(id2, x2, pt*pt)*
				CS_XX_XX(id1, id2, s, pt, y1 - y2)/(x1*x2*s)*1e3;
			//1e3 is to get pb instead of fb
		}
	}
	return result;
}

//dsigma/dpT
double GetDsigmaDpT(const double pt)
{	
	//result for the integral over y1
	double int_y1_result = 0.;
	
	for (double y1 = -Par.abs_max_y; y1 <= Par.abs_max_y; 
		y1 += 2.*Par.abs_max_y/Par.nsteps)
	{
		//result for the integral over y2
		double int_y2_result = 0.;
		//normalization for integral over y1 and y2 since y1 and y2 are interchangeable in the cross section
		double int_y2_norm = 0.;
		//y range can differ from defined abs_max_y 
		//Since it may be kinematically imposible for y to have value outside current_abs_max_y
		//Only 1 variable is required since y1 and y2 are interchangeable in the cross section
		double current_abs_max_y = 0.;
		
		for (double y2 = -Par.abs_max_y; y2 <= Par.abs_max_y; 
			y2 += 2.*Par.abs_max_y/Par.nsteps)
		{
			const double ptmax = Par.energy/(2.*cosh(abs(y1-y2)/2.));
			
			if (pt > ptmax) continue;
			
			double x1 = X1(pt, Par.energy, y1, y2);
			double x2 = X2(pt, Par.energy, y1, y2);
			const double s = Par.s*x1*x2;
			
			if (x1 >= 1. || x2 >= 1. || pt*pt*4. >= s) continue;

			int_y2_result += DsigmaDpTDy1Dy2(pt, s, y1, y2, x1, x2);
			int_y2_norm += 1.;
			current_abs_max_y = Maximum(current_abs_max_y, abs(y2));
		}
		if (int_y2_norm > 0) int_y2_result *= 4.*current_abs_max_y*current_abs_max_y/(int_y2_norm*int_y2_norm);
		int_y1_result += int_y2_result;
	}
	return int_y1_result;
}

//dsigma/d dDeltay
double GetDsigmaDdy(const double delta_y)
{	
	//result and normalization for the integral over pt
	double int_pt_result = 0.;
	double int_pt_norm = 0.;
	
	const double ptmax = Par.energy/(2.*cosh(delta_y/2.));
	
	for (double pt = Par.ptmin; pt <= ptmax; pt += 1.)
	{
		double current_result = 0.;
		//y range can differ from defined abs_max_y 
		//Since it may be kinematically imposible for y to have value outside current_abs_max_y
		//Only 1 variable is required since y1 and y2 are simmetrical relative to each other
		double current_abs_max_y = 0.;

		//result and normalization over y1 and y2
		double int_y_result = 0.;
		double int_y_norm = 0.;
		
		for (double y1 = -Par.abs_max_y; y1 <= Par.abs_max_y; 
			y1 += 2.*Par.abs_max_y/Par.nsteps)
		{
			double y2 = y1 - delta_y;
			if (abs(y2) <= Par.abs_max_y) 
			{
				const double x1 = X1(pt, Par.energy, y1, y2);
				const double x2 = X2(pt, Par.energy, y1, y2);
				
				if (x1 < 1. && x2 < 1.)
				{
					int_y_result += DsigmaDpTDy1Dy2(pt, Par.s*x1*x2, y1, y2, x1, x2);
					int_y_norm += 1.;
					current_abs_max_y = Maximum(current_abs_max_y, abs(y1), abs(y2));
				}
			}
			
			//switching y2 since it can have 2 different values
			y2 = delta_y + y1;
			if (abs(y2) <= Par.abs_max_y)
			{
				const double x1 = X1(pt, Par.energy, y1, y2);
				const double x2 = X2(pt, Par.energy, y1, y2);
				
				if (x1 < 1. && x2 < 1.)
				{
					int_y_result += DsigmaDpTDy1Dy2(pt, Par.s*x1*x2, y1, y2, x1, x2);
					int_y_norm += 1.;
					current_abs_max_y = Maximum(current_abs_max_y, abs(y1), abs(y2));
				}
			}
		}

		//multiplying cross section by y ranges
		if (int_y_norm > 0.) int_y_result *= 4.*current_abs_max_y*current_abs_max_y/int_y_norm;
		int_pt_result += int_y_result;
		int_pt_norm += 1.;
	}
	if (int_pt_norm > 0) int_pt_result *= ((ptmax - Par.ptmin))/int_pt_norm;
	return int_pt_result;
}

int main()
{
	TH1D dsigma_dpt = TH1D("dsigma_dpt", "dsigma/dpT", 200, 0, 200);
	TH1D dsigma_ddy = TH1D("dsigma_ddy", "dsigma/dDeltay", 100, 0, 
		static_cast<double>(ceil(Par.abs_max_y*2)));

	ProgressBar pbar = ProgressBar("FANCY");
	pbar.SetText("dsigma/dpT");
	
	//performing monte-carlo integration for dsigma/dpT and filling the hist with the result
	for (int i = 1; i <= dsigma_dpt.GetXaxis()->GetNbins(); i++)
	{
		pbar.Print(static_cast<double>(i)/static_cast<double>(dsigma_dpt.GetXaxis()->GetNbins()));
		const double pt = dsigma_dpt.GetXaxis()->GetBinCenter(i);
		dsigma_dpt.SetBinContent(i, GetDsigmaDpT(pt));
	}
	pbar.Print(1);
	
	pbar.Reset();
	pbar.SetText("dsigma/ddy");

	//performing monte-carlo integration for dsigma/ddeltay and filling the hist with the result
	for (int i = 1; i <= dsigma_ddy.GetXaxis()->GetNbins(); i++)
	{
		pbar.Print(static_cast<double>(i)/static_cast<double>(dsigma_ddy.GetXaxis()->GetNbins()));
		const double delta_y = dsigma_ddy.GetXaxis()->GetBinCenter(i);
		dsigma_ddy.SetBinContent(i, 1.);//GetDsigmaDdy(delta_y));
	}
	pbar.Print(1);
	
	system("mkdir ../output");
	TFile output = TFile("../output/analytic.root", "RECREATE");

	dsigma_dpt.Write();
	dsigma_ddy.Write();

	output.Close();
	return 0;
}
