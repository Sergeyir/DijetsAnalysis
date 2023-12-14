#include <iostream>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TH1.h"

#include "LHAPDF/LHAPDF.h"

#include "../lib/ProgressBar.h"
#include "../lib/Box.h"

using namespace LHAPDF;

struct
{
	const double energy = 7000;
	std::string pdfset_name = "NNPDF31_lo_as_0118";
	const PDF *pdf = mkPDF(pdfset_name);

	const double s = energy*energy;

	const double abs_max_y = 4.7;

	const double ptmin = 0;
	const double ptmax = energy/2.;

	const double nsteps = 100.;
} Par;

//cross sections dsigma/dOmega for different processes
//qq'->qq'
double CS_QQp_QQp(const double s, const double t, const double u)
{
	return 1./(2.*s)*4./9.*(s*s + u*u)/(t*t);
}

//qq->qq
double CS_QQ_QQ(const double s, const double t, const double u)
{
	return 1./(4.*s)*(4./9.*((s*s + u*u)/(t*t) + (s*s + t*t)/(u*u) - 8./27.*s*s/(u*t)));
}

//qqbar->q'qbar'
double CS_QpQbarp_QpQbarp(const double s, const double t, const double u)
{
	return 1./(2.*s)*4./9.*(t*t+u*u)/(s*s);
}

//qqbar->qqbar
double CS_QQbar_QQbar(const double s, const double t, const double u)
{
	return 1./(2.*s)*(4./9.*((s*s+u*u)/(t*t) + (t*t + u*u)/(s*s)) - 8./27.*u*u/(s*t));
}

//qqbar->gg
double CS_QQbar_GG(const double s, const double t, const double u)
{
	return 1./(4.*s)*(32./27.*(t*t + u*u)/(t*u) - 8./3.*(t*t + u*u)/(s*s));
}

//gg->qqbar
double CS_GG_QQbar(const double s, const double t, const double u)
{
	return 1./(2.*s)*(1./6.*(t*t + u*u)/(t*u) - 3./8.*(t*t + u*u)/(s*s));
}

//gq->gq
double CS_GQ_GQ(const double s, const double t, const double u)
{
	return 1./(2.*s)*(-4./9.*(s*s + u*u)/(s*u) + (u*u + s*s)/(t*t));
}

//gg->gg
double CS_GG_GG(const double s, const double t, const double u)
{
	return 9./(8.*s)*(3. - (t*u)/(s*s) - (s*u)/(t*t) - (s*t)/(u*u));
}

//returns dsigma/dOmega for pid1+pid2->X+X process
double CS_XX_XX(const int pid1, const int pid2, const double s, const double t, const double u)
{
	if (pid1 == 0 && pid2 == 0)
	{
		return CS_GG_GG(s, t, u) + CS_GG_QQbar(s, t, u);
	}
	if ((pid1 != 0 && pid2 == 0) || (pid1 == 0 && pid2 != 0))
	{
		return CS_GQ_GQ(s, t, u);
	}
	if (pid1 == pid2)
	{
		return CS_QQ_QQ(s, t, u);
	}
	if (pid1 == -pid2)
	{
		return CS_QQbar_QQbar(s, t, u) + CS_QQbar_GG(s, t, u);
	}
	if ((pid1 > 0 && pid2) < 0 || (pid1 < 0 && pid2 > 0))
	{
		return CS_QpQbarp_QpQbarp(s, t, u);
	}
	return CS_QQp_QQp(s, t, u);
}

//variables shortcuts
double CosTheta(const double s, const double pt) {return sqrt(1.-(4.*pt*pt/s));}
double T(const double s, const double cos_theta) {return s/2.*(1.-cos_theta);}
double U(const double s, const double cos_theta) {return -s/2.*(1.+cos_theta);}

double X1(const double pt, const double sqrt_s, const double y1, const double y2)
{
	return 2.*pt/sqrt_s*exp((y1+y2)/2.)*cosh((y1-y2)/4.);
}

double X2(const double pt, const double sqrt_s, const double y1, const double y2)
{
	return 2.*pt/sqrt_s*exp(-(y1+y2)/2.)*cosh((y1-y2)/4.);
}

//dsigma/dpT
double dsigmadpT(const double pt)
{
	const double alpha_s = Par.pdf->alphasQ2(Par.s);
		
	double result = 0.;
	
	//all flavours loop
	for (int id1 = -5; id1 <= 5; id1++)
	{
		for (int id2 = -5; id2 <= 5; id2++)
		{
			double norm = 0.;
			double current_result = 0.;
			for (double y1 = -Par.abs_max_y; y1 <= Par.abs_max_y; 
				y1 += 2.*Par.abs_max_y/Par.nsteps)
			{
				for (double y2 = -Par.abs_max_y; y2 <= Par.abs_max_y; 
					y2 += 2.*Par.abs_max_y/Par.nsteps)
				{
					const double ptmax = Par.energy/(2.*cosh(abs(y1-y2)/2.));
					
					if (pt > ptmax) continue;
					
					double x1 = X1(pt, Par.energy, y1, y2);
					double x2 = X2(pt, Par.energy, y1, y2);
					const double s = Par.s*x1*x2;
					
					if (x1 > 1. || x2 > 1. || 4.*pt*pt > s) continue;
					
					double cos_theta = CosTheta(s, pt);
					if (y1 - y2 < 0) cos_theta *= -1.;
					
					current_result += Par.pdf->xfxQ2(id1, x1, pt*pt)*Par.pdf->xfxQ2(id2, x2, pt*pt)*
						CS_XX_XX(id1, id2, s, T(s, cos_theta), U(s, cos_theta))/x1/x2/s;
					norm += 1.;
				}
			}
			result += current_result/norm;
		}
	}
	return result*alpha_s*2.*3.14159265359/pt;
}

//dsigma/dp^2T d Deltay
double dsigmadp2TdDy(const double delta_y, const double pt)
{
	const double alpha_s = Par.pdf->alphasQ2(Par.s);
		
	double result = 0.;
	
	//all flavours loop
	for (int id1 = -5; id1 <= 5; id1++)
	{
		for (int id2 = -5; id2 <= 5; id2++)
		{
			for (double y1 = -Par.abs_max_y; y1 <= Par.abs_max_y; y1 += 2.*Par.abs_max_y/(Par.nsteps))
			{
				//delta function requires y2 to be equal to delta_y - y1 or y1 - delta_y
				double y2 = delta_y - y1;
				
				if (abs(y2) > Par.abs_max_y) continue;
				
				double x1 = X1(pt, Par.energy, y1, y2);
				double x2 = X2(pt, Par.energy, y1, y2);

				const double s = Par.s*x1*x2;

				result += Par.pdf->xfxQ2(id1, x1, Par.s)*
					Par.pdf->xfxQ2(id2, x2, Par.s)/(Par.nsteps*2.)*
					CS_XX_XX(id1, id2, s, T(s, CosTheta(s, pt)), U(s, CosTheta(s, pt)));
				
				//switching to y2 = delta_y - y2
				x1 = X1(pt, Par.energy, y1, -y2);
				x2 = X2(pt, Par.energy, y1, -y2);
				
				result += Par.pdf->xfxQ2(id1, x1, Par.s)*
					Par.pdf->xfxQ2(id2, x2, Par.s)/(Par.nsteps*2.)*
					CS_XX_XX(id1, id2, s, T(s, -CosTheta(s, pt)), U(s, -CosTheta(s, pt)));
				
				//denomitaor nsteps for all y1 steps
				//denominator 2 is for 2 possible values of y2
			}
		}
	}
	return result*alpha_s*4.*3.14159265359/Par.s;
}

int main()
{
	TH1D hist_njets = TH1D("analytic_jets_multiplicity", "jets", 200, 0, 200);

	ProgressBar pbar = ProgressBar("FANCY");

	for (int i = 1; i <= hist_njets.GetXaxis()->GetNbins(); i++)
	{
		pbar.Print(static_cast<double>(i)/static_cast<double>(hist_njets.GetXaxis()->GetNbins()));
		const double pt = hist_njets.GetXaxis()->GetBinCenter(i);
		hist_njets.SetBinContent(i, dsigmadpT(pt));
	}

	pbar.Print(1);
	
	system("mkdir ../output");
	TFile output = TFile("../output/analytic.root", "RECREATE");

	hist_njets.Write();

	output.Close();
	return 0;
}
