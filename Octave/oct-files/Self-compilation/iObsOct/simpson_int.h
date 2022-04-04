#ifndef __SIMPSON_INT_H
#define __SIMPSON_INT_H

#ifdef __cplusplus
extern "C" {
#endif

double simpson_int(double(*f)(double), double a, double b, double eps);


#ifdef __cplusplus
}
#endif

#endif
