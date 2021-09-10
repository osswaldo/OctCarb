#ifndef PDF_H
#define PDF_H

#include <vector>

class PDF {
	public:
		//function for making Pleczek transformation result
		std::vector<double> polynomFit(std::vector<double> Q, std::vector<double> IofQ, int Qmin);

		//function for getting PDF data
		std::vector<std::vector<double>> getPDFData(std::vector<double> Q, std::vector<double> IofQ, double Qmin, double Qmax, double Rstep, double Rmax, double dfactor, double scale);
};

#endif // iObsPDF_H
