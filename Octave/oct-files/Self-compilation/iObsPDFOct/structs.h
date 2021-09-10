#ifndef STRUCTS
#define STRUCTS

struct wavelength{
	std::string name;
	double value;
};

struct parsMean{
	std::string  symbol;
	std::string  meaning;
};

//read Data-File
struct dataFile {
	double cno;
	double mu;
	double beta;
	double a3;
	double da3;
	double sig3;
	double u3;
	double eta;
	double nu;
	double alpha;
	double sig1;
	double lcc;
	double q;
	double cN;
	double cO;
	double dan;
	double k;
	double const1;
	double const2;
	double g;
	double mue_ab;
	double p;
	double R;
	double t;
	double La;
	double lm;
	double kapa;
	double Nm;
	double N;
	double Lc;
	double kapc;
	double a3min;
	double eps3;
	double kapr;
};

//coherent scattering parameters for I_intra and I_inter
struct csp{
	double mu;
	double Nm;
	double a3min;
	double da3;
	double sig3;
	double u3;
	double eta;
	double nu;
	double lm;
	double sig1;
	//Oliver: eps1 added
	double eps1;
	double lcc;
	double q;
	double concn; // 0 <= concentration of phase <= 1
};

struct AutoFitElement{
	bool autoFit;  //include in automatic fitting
	double min;    //minimal value for call of fit
	double max;    //maximal value for call of fit
	bool useMin;   //use minimal value
	bool useMax;   //use maximal value
};

#endif // STRUCTS

