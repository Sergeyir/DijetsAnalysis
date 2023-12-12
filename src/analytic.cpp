#include <iostream>
#include <string>
#include <cmath>

#include "LHAPDF/LHAPDF.h"

using namespace LHAPDF;

struct
{
	const double energy = 7000;
	std::string pdfset_name = "NNPDF31_lo_as_0118";
	const PDF *pdf = mkPDF(pdfset_name);

	const double abs_max_y = 4.7;

	const double ptmin = 0;
	const double ptmax = energy/2.;

	const int nsteps = 1000;
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
double CS_QQb_GG(const double s, const double t, const double u)
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

//
double GetPDFInt(const int pid1, const int pid2)
{	
	double pdf_int = 0.;
	for (double x1 = 1e-3; x1 <= 1; x1 += 1e-3)
	{
		for (double x2 = 1e-3; x2 <= 1; x2 += 1e-3)
		{
			pdf_int += x1*x2*Par.pdf->xfxQ2(pid1, x1, Par.energy*Par.energy)*
				Par.pdf->xfxQ2(pid1, x2, Par.energy*Par.energy);
		}
	}
	return pdf_int;
}

//variables shortcuts
double CosTheta(const double pt) {return sqrt(1.-(4.*pt/Par.energy));}
double T(const double cos_theta) {return Par.energy/2.*(1.-cos_theta);}
double U(const double cos_theta) {return -1.*Par.energy/2.*(1.+cos_theta);}

double X1(const double pt, const double s, const double y1, const double y2)
{
	return 2.*pt/sqrt(s)*exp((y1+y2)/2.)*cosh((y1-y2)/4.);
}

double X2(const double pt, const double s, const double y1, const double y2)
{
	return 2.*pt/sqrt(s)*exp((y1+y2)/2.)*cosh((y1-y2)/4.);
}

double dsigmadpTdDy(const double dSdOmega, const double s, const double delta_y, double (*CS_XX_XX)(const double, const double, const double)
{
	const double alpha_s = Par.pdf->alphasQ2(Par.energy*Par.energy);
	const double cos_theta = CosTheta(pt);
	const double dsigma_domega = CS_XX_XX(s, T(cos_theta), U(cos_theta));
	
	for (double y1 = -Par.abs_max_y; y1 <= Par.abs_max_y; y1 += Par.abs_max_y/1000.)
	{
		for (double y2 = -Par.abs_max_y; y2 <= Par.abs_max_y; y2 += Par.abs_max_y/1000.)
		{
			
		}
	}

	return alpha_s*Par.energy/(4.*3.14159265359);
}

int main()
{
	return 0;
}
