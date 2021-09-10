#ifndef ENUMERATIONS
#define ENUMERATIONS

class Enumerations{

public:

	//enumerations for pars-array indices
	//do not change numbers!
	//1. index = parName
	//2. index = parType (0:value, 1:minimum, 2:maximum; 3:decrement; 4:increment)
	enum parName        {cno = 0, mu, beta, a3, da3, sig3, u3, eta, nu, alpha, sig1, lcc, q, cN, cO, dan, k, const1, const2, g, mue_ab, p, R, t, countOf_parName};
	enum parNameAutoFit {cnoAF = 0, muAF, betaAF, a3AF, da3AF, sig3AF, u3AF, etaAF, nuAF, alphaAF, sig1AF, lccAF, qAF, danAF, kAF, const1AF, const2AF, gAF, countOf_parNameAutoFit};
	enum parType        {val = 0, min, max, dec, inc, countOf_parType};

	enum drawingOrder {dotsInFront = 0, lineInFront};

	enum XFormat {twotheta_degree = 0 , s_oneDIVnanometer, s_oneDIVangstrom, countOf_XFormat};

	enum YFormat {linear = 0, log_10, countOf_YFormat};

	enum corrPositionMethod {Dz = 0, D2thtx, countOf_corrPositionMethod};

	enum corrIntensityMethod {LambertBeer = 0, Pragmatic, countOf_corrIntensityMethod};

	enum dataToWrite {DataCurve = 0, FitCurve, DeviationCurve, countOf_dataToWrite};

	enum errorType {NoError = 0, DialogCancelled, OtherError, countOf_errorType};

	enum radiationType {X_ray = 0, neutron, countOf_radiation};

	enum MathematicaKernelType {local = 0, remote, countOf_MathematicaKernelType};

};

#endif // ENUMERATIONS
