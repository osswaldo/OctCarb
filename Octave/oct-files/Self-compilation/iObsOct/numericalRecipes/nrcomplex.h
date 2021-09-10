#ifndef _NR_COMPLEXNR_H_
#define _NR_COMPLEXNR_H_

#ifndef _FCOMPLEXNR_DECLARE_T_
typedef struct FCOMPLEXNR {float r,i;} fcomplexNR;
#define _FCOMPLEXNR_DECLARE_T_
#endif /* _FCOMPLEXNR_DECLARE_T_ */

#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */

fcomplexNR Cadd(fcomplexNR a, fcomplexNR b);
fcomplexNR Csub(fcomplexNR a, fcomplexNR b);
fcomplexNR Cmul(fcomplexNR a, fcomplexNR b);
fcomplexNR ComplexNR(float re, float im);
fcomplexNR Conjg(fcomplexNR z);
fcomplexNR Cdiv(fcomplexNR a, fcomplexNR b);
float Cabs(fcomplexNR z);
fcomplexNR Csqrt(fcomplexNR z);
fcomplexNR RCmul(float x, fcomplexNR a);

#else /* ANSI */
/* traditional - K&R */

fcomplexNR Cadd();
fcomplexNR Csub();
fcomplexNR Cmul();
fcomplexNR ComplexNR();
fcomplexNR Conjg();
fcomplexNR Cdiv();
float Cabs();
fcomplexNR Csqrt();
fcomplexNR RCmul();

#endif /* ANSI */

#endif /* _NR_COMPLEXNR_H_ */
