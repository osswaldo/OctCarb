#include "pdf.h"
#include "calculations.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "numericalRecipes/nr.h"
#include "numericalRecipes/nrutil.h"
#include "numericalRecipes/svdfit.c"

#include "numericalRecipes/svdcmp.c"
#include "numericalRecipes/svbksb.c"
#include "numericalRecipes/svdvar.c"
#include "numericalRecipes/fleg.c"
#include "numericalRecipes/pythag.c"
#include "numericalRecipes/ran1.c"

void fpoly(float x, float p[], int terms) {
	int j;

	p[1] = 1.0;
	for (j = 2; j <= terms; j++) {
		p[j]=p[j-1]*x;
	}

	p[1] = 1.0;
	p[2] = 0.0;
	p[3] = pow(x, 2);
}

std::vector<double> PDF::polynomFit(std::vector<double> Q, std::vector<double> IofQ, int Qmin) {

	Qmin = 0;

	int i;
	float chisq,*x,*y,*sig,*a,*w,**cvm,**u,**v;

	std::vector<double> Qfit;
	std::vector<double> IofQfit;
	for (int i = 0; i < Q.size(); i++) {
		if (Q[i] >= Qmin) {
			Qfit.push_back(Q[i]);
			IofQfit.push_back(IofQ[i]);
		}
	}

	const int NPT = Qfit.size();
	const int NPOL = 3;
	const int NITER = 100;

	x=vector(1,NPT);
	y=vector(1,NPT);
	sig=vector(1,NPT);
	a=vector(1,NPOL);
	w=vector(1,NPOL);
	cvm=matrix(1,NPOL,1,NPOL);
	u=matrix(1,NPT,1,NPOL);
	v=matrix(1,NPT,1,NPOL);

	for (int i = 0; i < NPT; i++) {
		x[i] = Qfit[i];
		y[i] = IofQfit[i];
	}

	svdfit(x,y,sig,NPT,a,NPOL,u,v,w,&chisq,fpoly,NITER);
	svdvar(v,NPOL,w,cvm);

	std::vector<double> result;

	for (i = NPOL;i >= 1; i--) {
		result.push_back(a[i]);
		//std::cout << a[i] << std::endl;
	}

	free_matrix(v,1,NPOL,1,NPOL);
	free_matrix(u,1,NPT,1,NPOL);
	free_matrix(cvm,1,NPOL,1,NPOL);
	free_vector(w,1,NPOL);
	free_vector(a,1,NPOL);
	free_vector(sig,1,NPT);
	free_vector(y,1,NPT);
	free_vector(x,1,NPT);

	return result;
}

std::vector<std::vector<double>> PDF::getPDFData(std::vector<double> Qin, std::vector<double> IofQin, double Qmin, double Qmax, double Rstep, double Rmax, double dfactor, double scale) {

	// Qin: initial Q-vector
	// IofQin: initial intensity-Vektor

	double RCalc, PDFofRCalc;
	std::vector<double> R, PDFofR;
	std::vector<std::vector<double>> result;

	std::vector<double> Q, IofQ;

	for(int i = 0; i < Qin.size(); i++) {
		if (Qin[i] >= Qmin && Qin[i] <= Qmax) {
			Q.push_back(Qin[i]);
			IofQ.push_back(IofQin[i]);
		}
	}

	int idata,ifou;
	double factor, constant;
	int Rlen;

	double Qave,Qdiff,IofQave;
	double Qa, Qb, IofQa, IofQb, fofQa, fofQb, fofQave, x;

	factor = 1.0;    /* default */
	constant = 0.0;    /* default */

	Rlen = (int) (Rmax/Rstep + 1.5);

	for(ifou = 1; ifou <= Rlen; ifou++) {
		RCalc = (ifou - 1) * Rstep;
		PDFofRCalc = 0.0;
		for(idata = 0; idata < Q.size() - 1; idata++) {
			Qa = Q.at(idata);
			Qb = Q.at(idata + 1);
			Qave = (Qa + Qb) / 2.0;
			Qdiff = Qb - Qa;
			IofQa = IofQ.at(idata);
			IofQb = IofQ.at(idata+1);
			IofQave = (IofQa + IofQb) / 2.0;
			fofQa = (IofQa - constant) * Qa * sin(Qa * RCalc);
			fofQb = (IofQb - constant) * Qb * sin(Qa * RCalc);
			fofQave =
			  (IofQave - constant) * Qave * sin(Qa * RCalc);
			PDFofRCalc
			  += Qdiff/6.0 * ( fofQa + 4*fofQave + fofQb );
		}
		//std::cout << "Test 6" << std::endl;
		x = Q.at(0) * RCalc;
		PDFofRCalc += (IofQ.at(0) - constant) / (RCalc * RCalc) * (sin(x) - x * cos(x));
		PDFofRCalc *= 2.0/M_PI;
		PDFofRCalc *= factor;
		PDFofRCalc /= dfactor;
		PDFofRCalc /= scale;

		R.push_back(RCalc);
		PDFofR.push_back(PDFofRCalc);
	}

	result.push_back(R);
	result.push_back(PDFofR);

	return result;
}
