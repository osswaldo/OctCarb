#ifndef CALCULATIONS_H
#define CALCULATIONS_H

#include <cmath>
#include <complex>
#include <fstream>
#include <string>
#include <vector>

#include "enumerations.h"
#include "structs.h"

class Calculations{

private:
	float *IcomC_x; //array for interpolation
	float *IcomC_y; //array for interpolation
	int sizeOf_IcomC; //size of array

		/*
	float *fN_x; //array for interpolation
	float *fN_y; //array for interpolation
	int sizeOf_fN; //size of array

	float *fO_x; //array for interpolation
	float *fO_y; //array for interpolation
	int sizeOf_fO; //size of array
	*/

public:

	//constructor including initialization of interpolation functions
	//Calculations(std::string dataPath);
	Calculations();

	//destructor including release of memory
	~Calculations();

	//degree to rad
	double degToRad(double deg);

	//rad to degree
	double radToDeg(double rad);

	//s to Theta
	double sToTheta(double s, double wavelength);

	//T to s
	double thetaToS(double s, double wavelength);

	//hypergeometrical function 2F1
	std::complex<double> hypgeo2F1(double a, double b, double c, std::complex<double> z);

	//interpolation of the form factor of carbon
	double fCpara(double s);

	//interpolation of the form factor of carbon with respect to anisotropy
	double fCperp(double dan, double s);

	//interpolation of the form factor of nitrogen
	double fN(double s);

	//interpolation of the form factor of oxygen
	double fO(double s);

	//interpolation of the form factor of sulfor
	double fS(double s);

	//interpolation of the form factor of hydrogen
	double fH(double s);


	//incoherent scattering
	double IcomC(double s);
	double IcomH(double s);
	double IcomN(double s);
	double IcomO(double s);
	double IcomS(double s);

	double Recoil(double s, double wavelength);
	double Qabs(double s, double wavelength);
	double Dq(double s);
	double Dlam(double s, double wavelength);
	double Dlamc(double s, double wavelength);

	//attenuation function Q for compton scattering; due to secondary monochromation, pass-band b [angstr.]
	double Q(double b, double s, double wavelength);

	double IincC(bool useQ, double b, double s, double wavelength);

	//polarization factor P (P_ to avoid problems with a C macro in wstp.h)
	//double P_(double s, double wavelength);
	double P_(double s, double wavelength, bool polarized, double polarizationDegree);
	double PTheta_(double theta, bool polarized, double polarizationDegree);

	//I_intra - scattering from the layers

	double shk(int h, int k, double lcc);
	double g(double q, double lcc, int h, int k, double s);
	double g0(double q);
	//std::complex<double> F(double m, double nu, double lm, double lcc, double sig1, int h, int k, double s);

    //analytical expression for F without epsilon_1
	std::complex<double> F(double m, double nu, double lm, double lcc, double sig1, int h, int k, double s);

    //analytical expression for F with epsilon_1
    std::complex<double> F_eps1(double m, double nu, double lm, double lcc, double sig1, double eps1, int h, int k, double s);

	//analytical expression for J_hk without epsilon_1
	double Jhk(double nu, double lm, double lcc, double sig1, int h, int k, double q, double s);
	double JhkXprefactor(double nu, double lm, double lcc, double sig1, int h, int k, double q, double s); // J_hk times prefactor
	int Jhk_prefactor(int h, int k); //prefactor for J_hk

    //analytical expression for J_hk with epsilon_1
    double Jhk_eps1(double nu, double lm, double lcc, double sig1, double eps1, int h, int k, double q, double s);
    double JhkXprefactor_eps1(double nu, double lm, double lcc, double sig1, double eps1, int h, int k, double q, double s); // J_hk times prefactor

	double n0S0(double lcc);
    
    //analytical expression for Iintra with epsilon_1
	double Iintra(double nu, double lm, double lcc, double sig1, double q, double s);
    //analytical expression for Iintra with epsilon_1
    double Iintra_eps1(double nu, double lm, double lcc, double sig1, double eps1, double q, double s);


	//I_inter - scattering from the stacking of layers

	std::complex<double> H(double a3min, double da3, double sig3, double s);
	std::complex<double> HN(double mu, double Nm, double a3min, double da3, double sig3, double s);
	double I1DI0(double mu, double Nm, double a3min, double da3, double sig3, double s);
	double I0(double mu, double Nm, double a3min, double da3, double s);
	double I1D(double mu, double Nm, double a3min, double da3, double sig3, double s);
	double DD(double u3, double s);
	double Iinter(double mu, double Nm, double a3min, double da3, double sig3, double u3, double eta, double lcc, double q, double s);


	//I_coh - coherent scattering, csp = coherent scattering parameters for I_intra and I_inter, concs = concentrations of corresponding phase
	double Icoh(double cno, std::vector<csp> *csp, double cN, double cO, double cS, double cH, double dan, double s, Enumerations::radiationType radiationType);

	//I_obs - theoretical observed intensity
	double iObs(bool useA, double density, double absroptionConstantCorrection, double sampleThickness, bool transmission, bool useGradient, double g, bool useCorrAutoColl, double par_r, double par_delta, double par_l, double const1, double const2, bool useQ, double b, double k, double cno, std::vector<csp> *csp, double cN, double cO, double cS, double cH, double dan, double s, Enumerations::radiationType radiationType, double wavelength, bool polarizationCorrection, bool polarizedBeam, double polarizationDegree, bool coh, bool inc);

	//I_fit - fit function, csp = coherent scattering parameters for I_intra and I_inter
	double Ifit(bool useQ, double b, double cno, std::vector<csp> *csp, double cN, double cO, double cS, double cH, double dan, double s, Enumerations::radiationType radiationType, double wavelength, bool coh, bool inc);

	//I_incoh - incoherent scattering for all elements
	double Iinc(Enumerations::radiationType radiationType, double cno, double cN, double cO, double cS, double cH, bool useQ, double b, double s, double wavelength);

	//calculation of absorption factor A for transmission and reflection geometry. For transmission, a endless beamwidth can be used and and a infininte depth of penetration is assumed (as usual for neutron scattering). For X-ray diffraction using reflection geometry, the absorption is only non-constant for samples, which are thinner than the depth of penetration of X-rays in carbon. In general, this is only the case for thin films and not for powders.
	double A(double theta, double density, double absroptionConstantCorrection, double wavelength, Enumerations::radiationType radiationType, bool transmission, double sampleThickness);

	//exponential term for tuning of gradient/slope as function of s using parameter g
	double T(double g, double s);

	//conversion factor corrAutoColl for correction of automatic collimator
	//use to convert intensity for fixed colliminator to intensity for automatic collimator
	//needs radius of diffractometer r, divergence angle in degrees, irradiated length l, wavelength and scattering vector s
	double corrAutoColl(double r, double delta, double l, double wavelength, double s);

	//data correction

	//calculation of mut from thickness d of sample and absorption factor mue_ab
	double mut(double d, double mue_ab);

	//calculation of muR from radius of diffractometer r and absorption factor mue_ab
	double muR(double r, double mue_ab);

	//calculation of Dz as function of z = 2 theta [degrees]
	double Dz(double muR , double mut , double z);

	//calculation of intensity positions via D2thtx as function of i = position/index of intenstiy in QVector<QPointF> data
	double PositionByD2thtx(double R, double mue_ab, double t, int i);

	//calculation of opacity for intensity correction as function of z = 2 theta [degrees]
	double opacity(double d, double mue_ab, double z);

	//calculation of the pragmatic approach for intensity correction as function of z = 2 theta [degrees]
	double pragmatic(double p, double z);


	//auxiliary functions

	//hyperbolic cotangent
	double coth(double x);

	//cosecant
	double csc(double x);
};

#endif // CALCULATIONS_H
