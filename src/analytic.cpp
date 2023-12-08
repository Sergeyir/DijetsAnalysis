#include <iostream>
#include <string>

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

//variables shortcuts
double CosTheta(const double pt) {return sqrt(1.-(4.*pt/Par.energy));}
double T(const double cos_theta) {return Par.energy/2.*(1.-cos_theta);}
double U(const double cos_theta) {return -1.*Par.energy/2.*(1.+cos_theta);}

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

double GetCSSingle(const double pt)
{
	return 1.;
}

int main()
{
	std::cout << Par.pdf->xfxQ2(1, 0.3, Par.energy*Par.energy/200.);
	return 0;
}
