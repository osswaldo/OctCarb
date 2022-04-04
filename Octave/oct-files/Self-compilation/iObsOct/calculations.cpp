#include <chrono>
#include <iostream>

#include "calculations.h"

#include "numericalRecipes/nrutil.h"
#include "numericalRecipes/nrcomplex.h"

#include "numericalRecipes/bessj0.c"
#include "numericalRecipes/bsstep.c"
#include "numericalRecipes/factln.c"
//#include "numericalRecipes/gammln.c"
#include "numericalRecipes/hypdrv.c"
#include "numericalRecipes/hypgeo.c"
#include "numericalRecipes/hypser.c"
#include "numericalRecipes/mmid.c"
#include "numericalRecipes/midinf.c"
#include "numericalRecipes/odeint.c"
#include "numericalRecipes/pzextr.c"
#include "numericalRecipes/spline.c"
#include "numericalRecipes/trapzd.c"
#include "numericalRecipes/gaulag.c"

#include "igam.cpp"
#include "simpson_int.h"

static double lm = 10;
static double nu = 4;
static double s_hk = 1;
static double lcc = 1.42;
static double s = 0;
static double sig1 = 0;
static double eps1 = 0;

//constructor including initialization of interpolation functions
//Calculations::Calculations(std::string dataPath){
Calculations::Calculations(){

	/*
	std::vector<std::vector<double>> *temp = new std::vector<std::vector<double>>;
	int  skipLines = 4;

	std::ifstream DataFile;
	std::string dummy;
	std::string x;
	std::string y;

	int counter = 0;
	//BEGIN IIcomC
	DataFile.open(dataPath + "/IIcomC.txt", std::ifstream::in);
	DataFile.peek();
	for (int i = 0; i < skipLines; ++i)
		getline(DataFile, dummy,  '\n');
	counter = 0;
	while(DataFile.good()) {
		getline(DataFile, x, ';');
		getline(DataFile, y, '\n');
		std::vector<double> test;
		test.push_back(stod(x));
		test.push_back(stod(y));
		temp->push_back(test);
		counter++;
		DataFile.peek();
	}
	DataFile.close();

	sizeOf_IcomC = temp->size();
	IcomC_x = new float[sizeOf_IcomC];
	IcomC_y = new float[sizeOf_IcomC];
	for (int i = 0; i < temp->size(); ++i){
		IcomC_x[i] = temp->at(i).at(0);
		IcomC_y[i] = temp->at(i).at(1);
	}
	temp->clear();
	//END IIcomC

	//BEGIN fN
	DataFile.open(dataPath + "/fN.txt", std::ifstream::in);
	DataFile.peek();
	for (int i = 0; i < skipLines; ++i)
		getline(DataFile, dummy,  '\n');
	counter = 0;
	while(DataFile.good()) {
		getline(DataFile, x, ';');
		getline(DataFile, y, '\n');
		std::vector<double> test;
		test.push_back(stod(x));
		test.push_back(stod(y));
		temp->push_back(test);
		counter++;
		DataFile.peek();
	}
	DataFile.close();

	sizeOf_fN = temp->size();
	fN_x = new float[sizeOf_fN];
	fN_y = new float[sizeOf_fN];
	for (int i = 0; i < temp->size(); ++i){
		fN_x[i] = temp->at(i).at(0);
		fN_y[i] = temp->at(i).at(1);
	}
	temp->clear();
	//END fN

	//BEGIN fO
	DataFile.open(dataPath + "/fO.txt", std::ifstream::in);
	DataFile.peek();
	for (int i = 0; i < skipLines; ++i)
		getline(DataFile, dummy,  '\n');
	counter = 0;
	while(DataFile.good()) {
		getline(DataFile, x, ';');
		getline(DataFile, y, '\n');
		std::vector<double> test;
		test.push_back(stod(x));
		test.push_back(stod(y));
		temp->push_back(test);
		counter++;
		DataFile.peek();
	}
	DataFile.close();

	sizeOf_fO = temp->size();
	fO_x = new float[sizeOf_fO];
	fO_y = new float[sizeOf_fO];
	for (int i = 0; i < temp->size(); ++i){
		fO_x[i] = temp->at(i).at(0);
		fO_y[i] = temp->at(i).at(1);
	}
	temp->clear();
	//END fO

	delete temp;
	*/
}

//destructor including release of memory
Calculations::~Calculations(){
	/*
	delete[] IcomC_x;
	delete[] IcomC_y;

	delete[] fN_x;
	delete[] fN_y;

	delete[] fO_x;
	delete[] fO_y;
	*/
}

//degree to rad
double Calculations::degToRad(double deg) {
	return deg/180*M_PI;
}

//rad to degree
double Calculations::radToDeg(double rad) {
	return rad/M_PI*360;
}

//s to Theta
double Calculations::sToTheta(double s, double wavelength) {
	return asin(s * wavelength/2);
}

//s to Theta
double Calculations::thetaToS(double theta, double wavelength) {
	return 2 * sin(theta)/wavelength;
}

//hypergeometrical function 2F1
std::complex<double> Calculations::hypgeo2F1(double a, double b, double c, std::complex<double> z){
	fcomplexNR a_ = ComplexNR(a,0.);
	fcomplexNR b_ = ComplexNR(b,0.);
	fcomplexNR c_ = ComplexNR(c,0.);
	fcomplexNR z_ = ComplexNR(z.real(), z.imag());
	//fcomplexNR erg_ = hypgeo(a_,b_,c_,z_);
	fcomplexNR erg_ = hypgeo(a_,b_,c_,z_,0);
	std::complex<double> erg (erg_.r, erg_.i);
	return erg;
}


//ITC-Vol C, chapter 6.1
double getAtomicFormFactor(double s, double a1, double b1, double a2, double b2, double a3, double b3, double a4, double b4, double c, double d1, double d2, double d3, double d4, double cut) {
	if (s < cut) {
		return a1*exp(-b1*pow(s/2, 2)) + a2*exp(-b2*pow(s/2, 2)) + a3*exp(-b3*pow(s/2, 2)) + a4*exp(-b4*pow(s/2, 2)) + c;
	} else {
		return exp(d1+d2*(s/2)+d3/10*pow((s/2), 2)+d4/100*pow((s/2), 3));
	}
}

//interpolation of the form factor of carbon
double Calculations::fCpara(double s){
	return fCperp(0, s);
}

//interpolation of the form factor of carbon with respect to anisotropy
double Calculations::fCperp(double dan, double s){
	return getAtomicFormFactor(s, 2.31, 20.84339+dan, 1.02, 10.2075, 1.5886, 0.5687, 0.865, 51.6512, 0.2156, 1.7056, -1.5676, 1.1893, -0.42715, 3.74);
}
//interpolation of the form factor of hydrogen
double Calculations::fH(double s) {
	return getAtomicFormFactor(s, 0.489918, 20.6593, 0.262003, 7.74039, 0.196767, 49.5519, 0.049879, 2.20159, 0.001305, -10, 0, 0, 0, 3.81);
}

//interpolation of the form factor of nitrogen
double Calculations::fN(double s) {
	/*
	float fN_y2[sizeOf_fN]; //sec. deriv. of interpol. funct.
	spline(fN_x, fN_y, sizeOf_fN, 1E32, 1E32, fN_y2);
	float *result = new float;
	splint(fN_x, fN_y, fN_y2, sizeOf_fN, s, result);
	double result_double = (*result);
	delete result;
	return result_double;*/

	return getAtomicFormFactor(s, 12.2126, 0.0057, 3.1322, 9.8933, 2.0125, 28.9975, 1.1663,0.5826, -11.529, 1.5494, -1.2019, 0.51064, 0.2472, 3.88);
}

//interpolation of the form factor of oxygen
double Calculations::fO(double s) {
	/*float fO_y2[sizeOf_fO]; //sec. deriv. of interpol. funct.
	spline(fO_x, fO_y, sizeOf_fO, 1E32, 1E32, fO_y2);
	float *result = new float;
	splint(fO_x, fO_y, fO_y2, sizeOf_fO, s, result);
	double result_double = (*result);
	delete result;
	return result_double;*/

	return getAtomicFormFactor(s, 3.0485, 13.2771, 2.2868, 5.7011, 1.5463, 0.3239, 0.867, 32.9089, 0.2508,1.3053, -0.83742, -0.16738, 0.475, 3.79);
}

//interpolation of the form factor of sulfor
double Calculations::fS(double s) {
	return getAtomicFormFactor(s, 6.9053, 1.4679, 5.2034, 22.2151, 1.4379, 0.2536, 1.5863, 56.172, 0.8669, 1.104, -0.40325, 0.20094, -0.26058, 5.86);
}

//incoherent scattering
double Poly(double d, double c, double b, double a, double s) {
	return d*pow(s, 3) + c*pow(s, 2) + b*s + a;
}

//ITC-Vol C, chapter 7.4.3, spline fit (5 pieces, order 4)
double Calculations::IcomC(double s){
	/*double xData[] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.6, 1.8, 2, 3, 4, 5};
	double yData[] = {0, 0.355, 1.203, 2.141, 2.914, 3.451, 3.826, 4.068, 4.238, 4.372, 4.486, 4.589, 4.686, 4.78, 4.871, 5.051, 5.208, 5.348, 5.781, 5.93, 6};
	int sizeOf_IcomC = *(&xData + 1) - xData;
	IcomC_x = new float[sizeOf_IcomC];
	IcomC_y = new float[sizeOf_IcomC];
	for (int i = 0; i < sizeOf_IcomC; i++) {
		IcomC_x[i] = xData[i];
		IcomC_y[i] = yData[i];
	}

	float IcomC_y2[sizeOf_IcomC]; //sec. deriv. of interpol. funct.
	spline(IcomC_x, IcomC_y, sizeOf_IcomC, 1E32, 1E32, IcomC_y2);
	float *result = new float;
	splint(IcomC_x, IcomC_y, IcomC_y2, sizeOf_IcomC, s, result);
	double result_double = (*result);
	delete result;
	return result_double;*/

	if (s <= 0.4) {
		return Poly(-33.464901, 26.248635, 2.168414, -0.025249, s);
	} else if (s <= 0.8) {
		return Poly(11.133906, -13.909246, 7.104169, 2.900144, s-0.4);
	} else if (s <= 1.2) {
		return Poly(0.268793, -0.548559, 1.321047, 4.228902, s-0.8);
	} else if (s <= 1.8) {
		return Poly(-0.012456, -0.226008, 1.01122, 4.686754, s-1.2);
	} else if (s <= 5) {
		return Poly(0.030782, -0.248428, 0.726559, 5.209433, s-1.8);
	} else {
		return 6;
	}
}

//ITC-Vol C, chapter 7.4.3, spline fit (5 pieces, order 4)
double Calculations::IcomH(double s) {
	if (s <= 0.4) {
		return Poly(-8.00932253, 5.64507889, 0.92450043, -0.00040318, s);
	} else if (s <= 0.8) {
		return Poly(3.37822894, -3.96610816, 1.59608872, 0.76001297, s-0.4);
	} else if (s <= 1.2) {
		return Poly(-0.20221011, 0.08776657, 0.04475209, 0.9800778, s-0.8);
	} else if (s <= 1.8) {
		return Poly(0.21696067, -0.15488556, 0.01790449, 0.99907984, s-1.2);
	} else {
		return 1;
	}
}

//ITC-Vol C, chapter 7.4.3, spline fit (5 pieces, order 4)
double Calculations::IcomN(double s) {
	if (s <= 0.4) {
		return Poly(-31.81024, 27.76358, 1.12406, -0.0001, s);
	} else if (s <= 0.8) {
		return Poly(5.8462, -10.40871, 8.06601, 2.85585, s-0.4);
	} else if (s <= 1.2) {
		return Poly(2.57319, -3.39327, 2.54522, 4.79101, s-0.8);
	} else if (s <= 1.8) {
		return Poly(0.05541, -0.30544, 1.06574, 5.43086, s-1.2);
	} else if (s <= 5.5) {
		return Poly(0.02043, -0.20571, 0.75905, 5.97232, s-1.8);
	} else {
		return 7;
	}
}

//ITC-Vol C, chapter 7.4.3, spline fit (5 pieces, order 4)
double Calculations::IcomO(double s) {
	if (s <= 0.4) {
		return Poly(-30.69434, 29.26312, 0.18316, -0.00183, s);
	} else if (s <= 0.8) {
		return Poly(2.01067, -7.57008, 8.86037, 2.7891, s-0.4);
	} else if (s <= 1.2) {
		return Poly(3.59061, -5.15727, 3.76943, 5.25072, s-0.8);
	} else if (s <= 1.8) {
		return Poly(0.37942, -0.84854, 1.3671, 6.16313, s-1.2);
	} else if (s <= 6) {
		return Poly(0.01315, -0.16558, 0.75863, 6.75987, s-1.8);
	} else {
		return 8;
	}
}

//ITC-Vol C, chapter 7.4.3, spline fit (5 pieces, order 4)
double Calculations::IcomS(double s) {
	if (s <= 0.4) {
		return Poly(-47.5631949, 37.8134763, 5.115743, -0.0579072, s);
	} else if (s <= 0.8) {
		return Poly(16.1916979, -19.2623576, 12.5361905, 4.9945018, s-0.4);
	} else if (s <= 1.2) {
		return Poly(-1.5021487, 0.1676798, 4.8983194, 7.9632694, s-0.8);
	} else if (s <= 1.8) {
		return Poly(0.3550809, -1.6348986, 4.3114319, 9.8532884, s-1.2);
	} else if (s <= 4) {
		return Poly(0.1434152, -0.995753, 2.7330409, 11.9282815, s-1.8);
	} else if (s <= 16) {
		return Poly(0.0018618, -0.0492128, 0.4341161, 14.6486118, s-4);
	} else {
		return 16;
	}
}

double Calculations::Recoil(double s, double wavelength){
	double lam = wavelength;
	return (pow((lam / (lam + 0.0485 * pow(lam,2)/4. * pow(s,2))),3));
}

double Calculations::Qabs(double s, double wavelength){
	return (1./(1. + 0.0485 * wavelength * 3./8. * pow(s,2)));
}

double Calculations::Dq(double s){
	return (3.05 * pow(s,2) / (pow(0.53,2) + pow(s,2)));
}

double Calculations::Dlam(double s, double wavelength){
	return (pow(wavelength,2) / 137. * s * Dq(s));
}

double Calculations::Dlamc(double s, double wavelength){
	return ((0.0485 / 4. * pow(s,2) - 1.5E-3 * pow(s,4) / pow((pow(0.53,2) + pow(s,2)),2)) * pow(wavelength,2));
}

//attenuation function Q for compton scattering; due to secondary monochromation, pass-band b [angstr.]
double Calculations::Q(double b, double s, double wavelength){
	return (1 / (1 + Dlam(s,wavelength)/b) / (1. + pow(M_PI,2) * pow(Dlamc(s,wavelength),2) / pow((Dlam(s,wavelength) + b),2)));
}

double Calculations::IincC(bool useQ, double b, double s, double wavelength){
	return IcomC(s) * Recoil(s,wavelength) * Qabs(s,wavelength) * (useQ ? Q(b,s,wavelength) : 1);
}

//polarization factor P (P_ to avoid problems with one C macro)
/*double Calculations::P_(double s, double wavelength) {
	return ((1 + pow(cos(2. * asin(wavelength * s / 2.)),2)) / 2.);
}*/

//polarization factor P (P_ to avoid problems with one C macro) using theta as input
/*double Calculations::P_(double s, double wavelength, bool polarizedBeam, double polarizationDegree) {
	PTheta_(sToTheta(s, wavelength), polarizedBeam, polarizationDegree);
}*/

//polarization factor P (P_ to avoid problems with one C macro) using theta as input
double Calculations::PTheta_(double theta, bool polarizedBeam, double polarizationDegree) {
	if (polarizedBeam) {
		return 1 - pow(sin(2*theta), 2) * pow(cos(degToRad(polarizationDegree)), 2);
	} else {
		return 0.5 * (1 + pow(cos(2*theta), 2));
	}
}

//I_intra - scattering from the layers
double Calculations::shk(int h, int k, double lcc){
	return (sqrt(pow(h,2) + pow(k,2) + h*k) * 2. / 3. / lcc);
}

double Calculations::g(double q, double s_hk, double s) {
    double ret = 0.0;
    if (q == 0.) {
		ret = 1;
    } else if (q < 0) {
        double qAbs = abs(q);
        if (s < s_hk) {
            ret = (sqrt(qAbs) / atanh(sqrt(qAbs))/(1. + qAbs));
        } else {
            ret = ((1. + qAbs) * sqrt(qAbs) / atanh(sqrt(qAbs)) / ( pow((1. - qAbs),2) + 4. * qAbs * pow(s_hk,2)/pow(s,2) ));
        }
	} else {
		double qAbs = abs(q);
		ret = ((1. + qAbs) * sqrt(qAbs) / atanh(sqrt(qAbs)) / ( pow((1. + qAbs),2) + 4. * qAbs - 4 * qAbs ));
    }
    return ret;
}

double Calculations::g0(double q){
    if (q == 0) {
		return 1;
    } else {
        double qAbs = abs(q);
        double ret = ((1. + qAbs) * sqrt(qAbs) / atanh(sqrt(qAbs))/ pow((1. - qAbs),2));
        if (q > 0) {
            return ret;
        } else {
            return 1/ret;
        }
    }
}

//analytical expression for F without epsilon_1
std::complex<double> Calculations::F(double m, double nu, double lm, double lcc, double sig1, double s_hk, double s){
	double al = nu / lm;
    double b = 2. * pow(M_PI, 2) * pow(s_hk,2) * pow(sig1,2) * 2./3. / lcc;
	std::complex<double> t (al+b, -2. * M_PI * s);
	double a_hypgeo = ((1+m)/2);
	double b_hypgeo = ((2+m)/2);
	double c_hypgeo = 1;
    std::complex<double> z_hypgeo ( (-1. * (4. * pow(M_PI,2) * pow(s_hk,2))),0);
	z_hypgeo /= pow(t,2);
	return ( pow(al,m) * pow(t,-1-m) * (double)(exp(gammln(1+m))) * hypgeo2F1(a_hypgeo, b_hypgeo, c_hypgeo, z_hypgeo));
}

//numerical expression for F with epsilon_1
float integrate_Simpson(float (*func)(float), float a, float b, int maxIter = 20, double eps = 1.0e-6) {
    float trapzd(float (*func)(float), float a, float b, int n);

    float s,st,ost=0.0,os=0.0;

    for (int i = 1; i <= maxIter; i++) {
        st = trapzd(func, a, b, i);
        s = (4.0 * st - ost)/3.0;
        if (i > 5) {
            if (fabs(s-os) < eps*fabs(os) || (s == 0.0 && os == 0.0)) {
                return s;
            }
        }
        os=s;
        ost=st;
    }
    return 0.0;
}

//numerical expression for F with epsilon_1
float integrate_GauLag(float (*func)(float), float start, float end, float stepwidth) {
    //int points = (end-start)/stepwidth;
    int steps = 10;
    int points = 10;

    int i;
    float alf=1.0,checkw,checkx,xx,*x,*w;

    x=vector(1, points);
    w=vector(1, points);
    for (;;) {
        gaulag(x,w,steps,alf);
        //printf("%3s %10s %14s\n","#","x(i)","w(i)");
        for (i=1;i<=steps;i++) //printf("%3d %14.6e %14.6e\n",i,x[i],w[i]);
        checkx=checkw=0.0;
        for (i=1;i<=steps;i++) {
            checkx += x[i];
            checkw += w[i];
        }
        //printf("\nCheck value: %15.7e  should be: %15.7e\n",checkx,steps*(steps+alf));
        //printf("\nCheck value: %15.7e  should be: %15.7e\n",checkw, exp(gammln(1.0+alf)));
        /* demonstrate the use of GAULAG for an integral */
        for (xx=0.0,i=1;i<=steps;i++) xx += w[i]*func(x[i]);
        //printf("\nIntegral from gaulag: %12.6f\n",xx);
        //printf("Actual value:         %12.6f\n",1.0/(2.0*sqrt(2.0)));
    }
    free_vector(w,1,points);
    free_vector(x,1,points);
    return xx;
}

double Calculations::Pl(double r, double nu, double lm, double s) {
    double alpha = nu/lm;

    int nu_round;
    if ( nu - floor(nu) == 0.5 ) {
        if (((int)(floor(nu))) % 2 == 0) {
            nu_round = floor(nu);
        } else {
            nu_round = ceil(nu);
        }
    } else {
        nu_round = floor(nu + 0.5);
    }

    double prefactor = exp(-alpha*r)/nu;

    double sum = 0;

    for (int m = 0; m <= (nu_round - 1); m++){
        sum += (nu_round - m)/((double)(exp(factln(m)))) * pow(alpha*r,m);
    }

    return prefactor*sum;
}

double Calculations::Pd(double r, double lcc, double sig1, double eps1, double s_hk) {
    double errorSig = pow(sig1, 2)*2/(3*lcc)*r;
    double errorEps = pow(eps1, 2) * pow(r, 2);
    return exp(-2*pow(M_PI, 2) * (errorSig + errorEps) * pow(s_hk, 2));
}

double Calculations::PlPd(double r, double nu, double lm, double lcc, double sig1, double eps1, double s_hk, double s) {
    //std::cout << sig1 << "\t" << Pd(r, lcc, sig1, eps1, s_hk) << std::endl;
    return Pl(r, nu, lm, s) * Pd(r, lcc, sig1, eps1, s_hk);
}

double Calculations::Integrand(double r, double nu, double lm, double lcc, double sig1, double eps1, double s_hk, double s) {
    return PlPd(r, nu, lm, lcc, sig1, eps1, s_hk, s) * bessj0(2*M_PI*r*s_hk) * sin(2*M_PI*r*s);
}

float Integrand_Ronly(float r) {
    Calculations *calculations = new Calculations();

    return calculations->Integrand(r, nu, lm, lcc, sig1, eps1, s_hk, s);
}

double Calculations::Integrand_V2(double r, double nu, double lm, double lcc, double sig1, double eps1, double shk, double s) {
    double sinus = 2*M_PI*r*(s-s_hk)+M_PI_4;

    return 1/sqrt(r) * PlPd(r, nu, lm, lcc, sig1, eps1, shk, s) * sin(sinus);
}

float Integrand_V2_Ronly(float r) {
    Calculations *calculations = new Calculations();

    return calculations->Integrand_V2(r, nu, lm, lcc, sig1, eps1, s_hk, s);
}

//analytical expression for J_hk without epsilon_1
double Calculations::Jhk(double nu, double lm, double lcc, double sig1, int h, int k, double q, double s){
    s_hk = shk(h,k,lcc);
	int nu_round;
	if ( nu - floor(nu) == 0.5 )
		if (((int)(floor(nu))) % 2 == 0)
			nu_round = floor(nu);
		else
			nu_round = ceil(nu);
	else
		nu_round = floor(nu + 0.5);
    double t = g(q,s_hk,s) * 1 / nu / s;
	std::complex<double> sum (0.,0.);
	for (int m = 0; m <= (nu_round  - 1); ++m){
        sum += ( ( (nu_round - m) / (double)(exp(factln(m))) ) * F(m,nu,lm,lcc,sig1,s_hk,s) );
	}
	return t * sum.imag();
}

//numerical expression for J_hk with epsilon_1
double Calculations::Jhk_eps1(double nu_1, double lm_1, double lcc_1, double sig1_1, double eps1_1, int h, int k, double q, double s_1){
    nu = nu_1;
    lm = lm_1;
    lcc = lcc_1;
    sig1 = sig1_1;
    eps1 = eps1_1;
    s = s_1;
    s_hk = shk(h,k,lcc_1);

    double t = g(q,s_hk,s) * 1 / s;

    double acc= 0.001;
    double lim1=acc;
    float lim2=5*lm;

    //float t;

    //printf("\nIntegral of func with 2^(n-1) points\n");
    //printf("%6s %24s\n","n","approx. integral");
    //for (int i=1;i<=iterations;i++) {
    //    t=trapzd(F_esp1_ohneInt_Ronly,lim1,lim2,i);
    //    printf("%6d %20.6f\n",i,t);
    //}

    //return integrate_Simpson(F_esp1_ohneInt_Ronly, lim1, lim2, 100, acc);
    //return trapzd(F_esp1_ohneInt_Ronly, lim1, lim2, 2);
    //return integrate_GauLag(F_esp1_ohneInt_Ronly, lim1, lim2, acc);

    double sum = 0;
    double r = lim1+acc;
    double temp = Integrand_Ronly(lim1);
    double res = 0;
    do {
        res = Integrand_Ronly(r);
        sum += (temp + res)/2*acc;
        temp = res;
        r += acc;
    } while (r <= lim2);

    return t * sum;
}

//numerical expression for J_hk with epsilon_1
double Calculations::Jhk_eps1_V2(double nu_1, double lm_1, double lcc_1, double sig1_1, double eps1_1, int h, int k, double q, double s_1){
    nu = nu_1;
    lm = lm_1;
    lcc = lcc_1;
    sig1 = sig1_1;
    eps1 = eps1_1;
    s = s_1;
    s_hk = shk(h, k, lcc);

    double t = g(q,s_hk,s) * 1 / (sqrt(2)*M_PI*s*sqrt(s+s_hk));

    double acc= 1e-3;
    float lim1=acc;
    float lim2=5*lm;

    //float t;

    //printf("\nIntegral of func with 2^(n-1) points\n");
    //printf("%6s %24s\n","n","approx. integral");
    //for (int i=1;i<=iterations;i++) {
    //    t=trapzd(F_esp1_ohneInt_Ronly,lim1,lim2,i);
    //    printf("%6d %20.6f\n",i,t);
    //}

    //return integrate_Simpson(F_esp1_ohneInt_Ronly, lim1, lim2, 100, acc);
    //return trapzd(F_esp1_ohneInt_Ronly, lim1, lim2, 2);
    //return integrate_GauLag(F_esp1_ohneInt_Ronly, lim1, lim2, acc);

    double sum = 0;
    double r = lim1+acc;
    double temp = Integrand_V2_Ronly(lim1);
    double res = 0;
    do {
        res = Integrand_V2_Ronly(r);
        sum += (temp + res)/2*acc;
        temp = res;
        r += acc;
    } while (r <= lim2);

    double prefactor = 1;
    if (s < s_hk) {
      prefactor = 2*s/(s+s_hk);
    }

    return t * prefactor * sum;
}

int Calculations::Jhk_prefactor(int h, int k){
	int multiplicity = 12;
	if ( (h == 0) && (k == 0))
		multiplicity = 1;
	else
		if ( (k == 0) || (h == k) )
			multiplicity = 6;
	int structFactSquared = 1;
	if ( ( (h-k) % 3) == 0 )
		structFactSquared = 4;
	return (multiplicity * structFactSquared);
}

//analytical expression for JhkXprefactor without epsilon_1
double Calculations::JhkXprefactor_eps1(double nu, double lm, double lcc, double sig1, double eps1, int h, int k, double q, double s){
    return (Jhk_prefactor (h, k) * Jhk_eps1(nu, lm, lcc, sig1, eps1, h, k, q, s)); // J_hk times prefactor
}

//analytical expression for JhkXprefactor with epsilon_1
double Calculations::JhkXprefactor(double nu, double lm, double lcc, double sig1, int h, int k, double q, double s){
    return (Jhk_prefactor (h, k) * Jhk(nu, lm, lcc, sig1, h, k, q, s)); // J_hk times prefactor
}

double Calculations::n0S0(double lcc){
	return ( pow(3.,(3./2.)) * pow (lcc,2));
}

//analytical expression for Iintra without epsilon_1
double Calculations::Iintra(double nu, double lm, double lcc, double sig1, double q, double s){
	/*switch (radiationType){
		case Enumerations::X_ray :{
			double Jhk_sum = 0.;
			if (wavelength <= 0.5)
				for (int h = 1; h <= 7; ++h)
					for (int k = 0; k <= h; ++k)
						Jhk_sum += JhkXprefactor(nu, lm, lcc, sig1, h, k, q, s);
			else{
				Jhk_sum += JhkXprefactor(nu, lm, lcc, sig1, 1, 0, q, s);
				Jhk_sum += JhkXprefactor(nu, lm, lcc, sig1, 1, 1, q, s);
				Jhk_sum += JhkXprefactor(nu, lm, lcc, sig1, 2, 0, q, s);
				Jhk_sum += JhkXprefactor(nu, lm, lcc, sig1, 2, 1, q, s);
				Jhk_sum += JhkXprefactor(nu, lm, lcc, sig1, 3, 0, q, s);
				if (wavelength <= 1.25){
					Jhk_sum += JhkXprefactor(nu, lm, lcc, sig1, 2, 2, q, s);
					Jhk_sum += JhkXprefactor(nu, lm, lcc, sig1, 3, 1, q, s);
					Jhk_sum += JhkXprefactor(nu, lm, lcc, sig1, 4, 0, q, s);
					if (wavelength <= 1.){
						Jhk_sum += JhkXprefactor(nu, lm, lcc, sig1, 3, 2, q, s);
						Jhk_sum += JhkXprefactor(nu, lm, lcc, sig1, 4, 1, q, s);
						Jhk_sum += JhkXprefactor(nu, lm, lcc, sig1, 5, 0, q, s);
						if (wavelength <= 0.75){
							Jhk_sum += JhkXprefactor(nu, lm, lcc, sig1, 3, 3, q, s);
							Jhk_sum += JhkXprefactor(nu, lm, lcc, sig1, 4, 2, q, s);
							Jhk_sum += JhkXprefactor(nu, lm, lcc, sig1, 5, 1, q, s);
							Jhk_sum += JhkXprefactor(nu, lm, lcc, sig1, 6, 0, q, s);
							Jhk_sum += JhkXprefactor(nu, lm, lcc, sig1, 4, 3, q, s);
							Jhk_sum += JhkXprefactor(nu, lm, lcc, sig1, 5, 2, q, s);
						}
					}
				}
			}
			return ( 1 / n0S0(lcc) * Jhk_sum );
			break;
		}
		case Enumerations::neutron :{
			double Jhk_sum = 0.;
			for (int h = 1; h <= nu; ++h)
				for (int k = 0; k <= h; ++k)
					Jhk_sum += JhkXprefactor(nu, lm, lcc, sig1, h, k, q, s);
			return ( 1 / n0S0(lcc) * Jhk_sum );
			break;
		}
		default :
			return 0.;
			break;
	}*/

    // Standard
    nu = 7;
	double Jhk_sum = 0.;
    double max = pow(nu+1, 2);
    for (int h = 1; h <= nu; h++) {
        for (int k = 0; k <= h; ++k) {
            if ((pow(h, 2) + pow(k, 2) + h * k) < max) {
                double nuTemp = nu;
                nuTemp = 4;
                if (nuTemp > 10) {
                    nuTemp = 10;
                }
                Jhk_sum += JhkXprefactor(nuTemp, lm, lcc, sig1, h, k, q, s);
            }
        }
    }

    /*
    // WANS Paper vergleich
    double Jhk_sum = 0.;
    for (int h = 1; h <= 7; ++h)
        for (int k = 0; k <= h; ++k)
            Jhk_sum += JhkXprefactor(nu, lm, lcc, sig1, h, k, q, s);
    */
	return ( 1 / n0S0(lcc) * Jhk_sum );
}

//analytical expression for Iintra with epsilon_1
double Calculations::Iintra_eps1(double nu, double lm, double lcc, double sig1, double eps1, double q, double s){
    nu = 7;

    double Jhk_sum = 0.;
    double max = pow(nu+1, 2);
    for (int h = 1; h <= nu; h++) {
        for (int k = 0; k <= h; ++k) {
            if ((pow(h, 2) + pow(k, 2) + h * k) < max) {
                //double nuTemp = nu;
                double nuTemp = 4;
                if (nuTemp > 10) {
                    nuTemp = 10;
                }
                //if (eps1 == 0) {
                //    Jhk_sum += JhkXprefactor(nuTemp, lm, lcc, sig1, h, k, q, s);
                //} else {
                    Jhk_sum += JhkXprefactor_eps1(nuTemp, lm, lcc, sig1, eps1, h, k, q, s);
                //}
            }
        }
    }

    /*
    // WANS Paper vergleich
    double Jhk_sum = 0.;
    for (int h = 1; h <= 7; ++h)
        for (int k = 0; k <= h; ++k)
            Jhk_sum += JhkXprefactor(nu, lm, lcc, sig1, h, k, q, s);
    */
    return ( 1 / n0S0(lcc) * Jhk_sum );
}

//I_inter - scattering from the stacking of layers

std::complex<double> Calculations::H(double a3min, double da3, double sig3, double s){
	std::complex<double> i (0.,1.);
	return ( exp(2. * M_PI * i * a3min * s) / pow( ( 1. - 2. * M_PI * i * pow(sig3,2) * s / da3) ,( pow(da3,2) / pow(sig3,2) ) ) );
}

std::complex<double> Calculations::HN(double mu, double Nm, double a3min, double da3, double sig3, double s){
	return (pow( (1. - Nm * log(H(a3min,da3,sig3,s)) / mu ),-mu));
}

double Calculations::I1DI0(double mu, double Nm, double a3min, double da3, double sig3, double s){
	std::complex<double> HH;
	std::complex<double> HHN;
	HH = H(a3min,da3,sig3,s);
	HHN = HN(mu,Nm,a3min,da3,sig3,s);
	std::complex<double> t;
	t = ( (1. + HH)/(1. - HH) - 2. * HH * (1. - HHN) / Nm / pow((1. - HH),2) );
	return t.real();
}

double Calculations::I0(double mu, double Nm, double a3min, double da3, double s){
	double a3m = a3min + da3;
	return ( 1 /(2 * pow(M_PI,2) * pow(s,2) * Nm * pow(a3m,2)) *
			 (1 - (cos(mu * atan(2 * M_PI * a3m * Nm * s / mu))) / (pow((1 + 4 * pow(M_PI,2) * pow (a3m,2) * pow(Nm,2) * pow(s,2) / pow(mu,2) ), (mu/2.))) ) );
}

double Calculations::I1D(double mu, double Nm, double a3min, double da3, double sig3, double s){
	return (I1DI0(mu,Nm,a3min,da3,sig3,s)-I0(mu, Nm, a3min, da3, s));
}

double Calculations::DD(double u3, double s){
	return (exp(-4. * pow(M_PI,2) * pow(u3,2) * pow(s,2)));
}

double Calculations::Iinter(double mu, double Nm, double a3min, double da3, double sig3, double u3, double eta, double lcc, double q, double s){
	return (g0(q) * 2 / M_PI / pow(s,2) / n0S0(lcc) * (1. + eta * DD(u3,s) * (I1D(mu,Nm,a3min,da3,sig3,s) -1)));
}


//I_coh - coherent scattering, csp = coherent scattering parameters for I_intra and I_inter, concs = concentrations of corresponding phase
double Calculations::Icoh(double cno, std::vector<csp> *csp, double cN, double cO, double cS, double cH, double dan, double s, Enumerations::radiationType radiationType){
	double I_inter = 0.;
	double I_intra = 0.;
	int size = csp->size();
    for (int i = 0; i < size; ++i){
        I_inter += ( (csp->at(i).concn > 0.) ? ( csp->at(i).concn * Iinter(csp->at(i).mu, csp->at(i).Nm, csp->at(i).a3min, csp->at(i).da3, csp->at(i).sig3, csp->at(i).u3, csp->at(i).eta, csp->at(i).lcc, csp->at(i).q, s) ) : (0.) );
        double eps1 = csp->at(i).eps1;
        if (eps1 == 0.0) {
            I_intra += ( (csp->at(i).concn > 0.) ? ( csp->at(i).concn * Iintra(csp->at(i).nu, csp->at(i).lm, csp->at(i).lcc, csp->at(i).sig1, csp->at(i).q, s) ) : (0.) );
        } else {
            I_intra += ( (csp->at(i).concn > 0.) ? ( csp->at(i).concn * Iintra_eps1(csp->at(i).nu, csp->at(i).lm, csp->at(i).lcc, csp->at(i).sig1, eps1, csp->at(i).q, s) ) : (0.) );
        }
    }

	switch (radiationType){
		case Enumerations::X_ray :{
			double cC = 1 - cH - cN - cO - cS;
			double cCordered = cC * (1-cno);

            return cCordered * (pow(fCperp(dan,s),2) * I_inter + pow(fCpara(s),2) * I_intra);

			/*
			return ( ( (1-cno) *
					 ( (pow(fCperp(dan,s),2) * I_inter) +
					   (pow(fCpara(s),2) * I_intra ) ) +
					 (pow(fCpara(s),2) * cno) ) +
					 (pow(fN(s),2)*cN) + (pow(fO(s),2)*cO) + (pow(fH(s),2)*cH) + (pow(fS(s),2)*cS));
					 */

			break;
		}
		case Enumerations::neutron :{ //structure factors taken from http://www.ncnr.nist.gov/resources/n-lengths/
			//https://www.ncnr.nist.gov/resources/n-lengths
			/*return ( ( (1-cno) *
					 ( (6.6460 * I_inter) +
					   (6.6460 * I_intra) ) +
					 (6.6460 * cno)) +
					 (9.36*cN) + (5.803*cO) + (-3.739*cH) + (2.847*cS) );*/

			return 6.6460 * (I_inter + I_intra);

			break;
		}
		default :
			return 0.;
			break;
	}
}


//I_obs - theoretical observed intensity
double Calculations::iObs(bool useA, double density, double absorptionConstantCorrection, double sampleThickness, bool transmission, bool useGradient, double g, bool useCorrAutoColl, double par_r, double par_delta, double par_l, double const1, double const2, bool useQ, double b, double k, double cno, std::vector<csp> *csp, double cN, double cO, double cS, double cH, double dan, double s, Enumerations::radiationType radiationType, double wavelength, bool polarizationCorrection, bool polarizedBeam, double polarizationDegree, bool coh, bool inc) {
	double theta = sToTheta(s, wavelength);
    double factA = (useA ? A(theta, density, absorptionConstantCorrection, wavelength, radiationType, transmission, sampleThickness) : 1);
	double factGradient = (useGradient ? T(g, s) : 1);
	double factPolarization = (polarizationCorrection && radiationType == Enumerations::X_ray ? (PTheta_(theta, polarizedBeam, polarizationDegree)) : 1);

	//Only Useful for XRD in Bragg-Brentano or Debye-Scherrer geometry
	double factCorrAutoColl = (useCorrAutoColl ? (1. / corrAutoColl(par_r, par_delta, par_l, wavelength, s)) : 1);

	double iFitCalc = Ifit(useQ, b, cno, csp, cN, cO, cS, cH, dan, s, radiationType, wavelength, coh, inc);

	return pow(10, (log10(k * factA * factGradient * factCorrAutoColl * factPolarization * iFitCalc + const1) + const2));
}

//I_fit - fit function, csp = coherent scattering parameters for I_intra and I_inter
double Calculations::Ifit(bool useQ, double b, double cno, std::vector<csp> *csp, double cN, double cO, double cS, double cH, double dan, double s, Enumerations::radiationType radiationType, double wavelength, bool coh, bool inc){
	switch (radiationType){
		case Enumerations::X_ray :{
			double scat = 0.;
			if (coh) {
				if (inc) {
					scat = Icoh(cno, csp, cN, cO, cS, cH, dan, s, radiationType) + Iinc(radiationType, cno, cN, cO, cS, cH, useQ, b, s, wavelength);
				} else {
					scat = Icoh(cno, csp, cN, cO, cS, cH, dan, s, radiationType);
				}
			} else if (inc) {
				scat = Iinc(radiationType, cno, cN, cO, cS, cH, useQ, b, s, wavelength);
			}
			return scat;
			break;
		}
		case Enumerations::neutron :{
			return Icoh(cno, csp, cN, cO, cS, cH, dan, s, radiationType);
			break;
		}
		default :
			return 0.;
			break;
	}
}

//I_fit - fit function, csp = coherent scattering parameters for I_intra and I_inter
double Calculations::Iinc(Enumerations::radiationType radiationType, double cno, double cN, double cO, double cS, double cH, bool useQ, double b, double s, double wavelength) {
	switch(radiationType) {
		case Enumerations::X_ray :{
            double cC = 1 - cH - cN - cO - cS;
            double cCdisordered = cC *(cno);
            //double cCordered = cC *(1-cno);

            double comptenCorr = Recoil(s,wavelength) * Qabs(s,wavelength) * (useQ ? Q(b,s,wavelength) : 1);

            double IincC = cC * IcomC(s) * comptenCorr;
            double IincForeignCom = (cH * IcomH(s) + cN * IcomN(s) + cO * IcomO(s) + cS * IcomS(s)) * comptenCorr;

            return	IincC +
                    IincForeignCom +
                    pow(fCpara(s),2) * cCdisordered +
					pow(fH(s),2) * cH +
					pow(fN(s),2) * cN +
					pow(fO(s),2) * cO +
					pow(fS(s),2) * cS;
		}
		case Enumerations::neutron :{
			return 0;
		}
	}
}

//calculation of absorption factor A for transmission and reflection geometry. For transmission, a endless beamwidth can be used and and a infininte depth of penetration is assumed (as usual for neutron scattering). For X-ray diffraction using reflection geometry, the absorption is only non-constant for samples, which are thinner than the depth of penetration of X-rays in carbon. In general, this is only the case for thin films and not for powders.
double Calculations::A(double theta, double density, double absorptionConstantCorrection, double wavelength, Enumerations::radiationType radiationType, bool transmission, double sampleThickness) {
	double l = sampleThickness;
	double mue_ab = 0;
	if (radiationType == Enumerations::X_ray) {
        mue_ab = pow(wavelength, -3.089) * pow(10, 1.081) * density * absorptionConstantCorrection;
	} else {
		double waveVector = 1/wavelength * 2 * M_PI;
		double sigAbsInc = 0.001E-24;
        double sigAbs = 4 * M_PI/waveVector * sigAbsInc * absorptionConstantCorrection;
		double n = density * 6.02214076E23/12.011;
		mue_ab = n * sigAbs;
	}

    mue_ab = mue_ab * absorptionConstantCorrection;

	if (transmission) {
		return l/cos(theta)*exp(-l*mue_ab/cos(theta));
	} else {
		return 1/(2*mue_ab)*(1-exp(-2*l*mue_ab/sin(theta)));
	}
}

//exponential term for tuning of gradient/slope as function of s using parameter g
double Calculations::T(double g, double s){
	return exp(g * s);
}

//conversion factor corrAutoColl for correction of automatic collimator
//use to convert intensity for fixed colliminator to intensity for automatic collimator
//needs radius of diffractometer r, divergence angle in degrees, irradiated length l, wavelength and scattering vector s
double Calculations::corrAutoColl(double r, double deltaDegree, double l, double wavelength, double s){
	double deltaRad = degToRad(deltaDegree); //delta in radians
	return ( (r * (wavelength * s / 2) * sin(deltaRad)) / ( l * ( pow((wavelength * s / 2), 2) - pow(sin(deltaRad/2), 2) ) ) );
}


//data correction

//calculation of mut from thickness d of sample and absorption factor mue_ab
double Calculations::mut(double d, double mue_ab){
	return (d * mue_ab);
}

//calculation of muR from radius of diffractometer r and absorption factor mue_ab
double Calculations::muR(double r, double mue_ab){
	return (r * mue_ab);
}

//calculation of Dz as function of z = 2 theta [degrees]
double Calculations::Dz(double muR_, double mut_, double z){
	double degree = (M_PI / 180.);
	return ( sin(z * degree) / 2. / muR_ * (1. - 2. * mut_ / sin(z/2. * degree) / (exp(2 * mut_ / sin(z/2 * degree)) - 1)) * 180 / M_PI);
}

//calculation of intensity positions via D2thtx as function of i = position/index of intenstiy in std::vector<QPointF> data
double Calculations::PositionByD2thtx(double R, double mue_ab, double t, int i){
	double x = (10. + ( 0.25 * i ) );
	double xx = ( (x/2.) * (M_PI / 180.) ); //x divided by 2 in radians
	double muet = mue_ab * t; //product of mue_ab and t
	double D2thtx = ( ( 1. / (mue_ab * M_PI * R) ) * 180. * cos(xx) * ( muet - muet * coth( muet * csc(xx) ) + sin(xx) ) );
	return  x + D2thtx;
}

//calculation of opacity for intensity correction as function of z = 2 theta [degrees]
double Calculations::opacity(double d, double mue_ab, double z){
	return ( exp( (2 * d * mue_ab) / sin( (z/2) * (M_PI / 180.) ) ) );
}

//calculation of the pragmatic approach for intensity correction as function of z = 2 theta [degrees]
double Calculations::pragmatic(double p, double z){
	return ( exp( -1. * p * pow(z,2) ) );
}


//auxiliary functions

//hyperbolic cotangent
double Calculations::coth(double x){
	return ( cosh(x) / sinh(x) );
}

//cosecant
double Calculations::csc(double x){
	return ( 1. / sin (x) );
}
