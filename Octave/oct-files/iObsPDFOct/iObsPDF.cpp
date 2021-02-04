#include "octave/oct.h"
#include "calculations.cpp"
#include "pdf.cpp"

#include <iostream>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>


DEFUN_DLD (iObsPDF, args, nargout, "Calculates the iObs function for Octave 5") {
	if (args.length() == 50) {
		double cno				= args(0).double_value();
		double mu				= args(1).double_value();
		double beta				= args(2).double_value();
		double a3				= args(3).double_value();
		double da3				= args(4).double_value();
		double sig3				= args(5).double_value();
		double u3				= args(6).double_value();
		double eta				= args(7).double_value();
		double nu				= args(8).double_value();
		double alpha			= args(9).double_value();
		double lcc				= args(10).double_value();
		double sig1				= args(11).double_value();
		double q				= args(12).double_value();
		double cH				= args(13).double_value();
		double cN				= args(14).double_value();
		double cO				= args(15).double_value();
		double cS				= args(16).double_value();
		double dan				= args(17).double_value();
		double k				= args(18).double_value();
		double const1			= args(19).double_value();
		double const2			= args(20).double_value();
		bool useQ				= args(21).bool_value();
		double b				= args(22).double_value();
		bool useA				= args(23).bool_value();
		double density			= args(24).double_value();
		double sampleThickness	= args(25).double_value();
		double transmission		= args(26).bool_value();
		double absroptionCorrection = args(27).double_value();
		bool useP				= args(28).bool_value();
		bool polarizedBeam		= args(29).bool_value();
		double polarizationDegree = args(30).double_value();
		bool useGradient		= args(31).bool_value();
		double g				= args(32).double_value();
		bool useCorrAutoColl	= args(33).bool_value();
		double par_r			= args(34).double_value();
		double par_delta		= args(35).double_value();
		double par_l			= args(36).double_value();
		int radiation			= args(37).int_value();
		double wavelength		= args(38).double_value();
		RowVector  S			= args(39).row_vector_value();
		bool coh				= args(40).bool_value();
		bool inc				= args(41).bool_value();
		double Smin				= args(42).double_value();
		double Smax				= args(43).double_value();
		double Sstep			= args(44).double_value();
		double Rmin				= args(45).double_value();
		double Rmax				= args(46).double_value();
		double Rstep			= args(47).double_value();
		double dfactor			= args(48).double_value();
		double scale			= args(49).double_value();
		
		std::vector<csp> *cspData = new std::vector<csp>(1);
		(cspData->data()[0]).mu    = mu;
		(cspData->data()[0]).Nm    = mu/beta;
		(cspData->data()[0]).a3min = a3-da3;
		(cspData->data()[0]).da3   = da3;
		(cspData->data()[0]).sig3  = sig3;
		(cspData->data()[0]).u3    = u3;
		(cspData->data()[0]).eta   = eta;
		(cspData->data()[0]).nu    = nu;
		(cspData->data()[0]).lm    = nu/alpha;
		(cspData->data()[0]).sig1  = sig1;
		(cspData->data()[0]).lcc   = lcc;
		(cspData->data()[0]).q     = q;
		(cspData->data()[0]).concn = 1;

		Enumerations enums;
		Enumerations::radiationType radiationType;
		switch(radiation) {
			case 0:
				radiationType = enums.X_ray;
				break;
			case 1:
				radiationType = enums.neutron;
				break;
			default:
				radiationType = enums.X_ray;
				break;
		}

		int points = (Smax-Smin)/Sstep+1.5;

		double Qmin = 2*M_PI*Smin;
		double Qmax = 2*M_PI*Smax;

		std::vector<double> Q;
		std::vector<double> IofQin;

		Calculations *calculations = new Calculations();
		PDF *pdf = new PDF();

		double s = Smin;

		for (octave_idx_type i = 0; i < points; i++) {
			IofQin.push_back(calculations->iObs(useA, density, absroptionCorrection, sampleThickness, transmission, useGradient, g, useCorrAutoColl, par_r, par_delta, par_l, const1, const2, useQ, b, k, cno, cspData, cN, cO, cS, cH, dan, s, radiationType, wavelength, useP, polarizedBeam, polarizationDegree, coh, inc));

			Q.push_back(2*M_PI*s);
			s += Sstep;
		}

		std::vector<double> Qin = Q;

		std::vector<std::vector<double>> iObsPDFOut = pdf->getPDF(Qin, IofQin, Qmin, Qmax, Rmin, Rmax, Rstep, dfactor, scale);

		double Rs = iObsPDFOut.at(0).size();

		Matrix iObsPDF = Matrix(2, Rs);
		for (int i = 0; i < Rs; i++) {
			iObsPDF(0, i) = iObsPDFOut.at(0).at(i);
			iObsPDF(1, i) = iObsPDFOut.at(1).at(i);
		}

		delete calculations;
		cspData->clear();
		delete cspData;

		return octave_value(iObsPDF);
	} else {
		return octave_value(0.);
	}
}
