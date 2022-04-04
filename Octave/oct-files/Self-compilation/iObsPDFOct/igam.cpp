//Incomplete gamma integral

//https://en.wikipedia.org/wiki/Incomplete_gamma_function#Software_Implementation
//https://www.degruyter.com/document/doi/10.1515/9781501506062-023/pdf

//#include "numericalRecipes/nrutil.h"
//#include "numericalRecipes/nrcomplex.h"

#include <random>
#include <math.h>

// Continued Fraction Computation
// 6.5.31 Handbook of Mathematical Functions , page 263
// Recursive implementation
double upper_incomplete_gamma(double a, double x, int d = 0, int iterations = 100) {
   if (d == iterations) {
        if ((d % 2) == 1) {
            return 1.0; // end iterations
        }
    } else {
        double m = d/2;
        return x + (m - a);
    }

    if (d == 0) {
        return pow(x, a) * exp(-x)/upper_incomplete_gamma(a, x, d++);
    } else if ((d % 2) == 1) {
        double m = 1.0 + ((d -1.0)/2.0);
        return x + (m - a)/upper_incomplete_gamma(a, x, d++);
    } else {
        double m = d/2;
        return 1 + m/upper_incomplete_gamma(a, x, d++);
    }
}


// 6.5.31 Handbook of Mathematical Functions , page 263
// Recursive implementation
double upper_incomplete_gamma2(double a, double x, int d = 0, int iterations = 100) {
    if (d == iterations) {
        return 1.0; // end iterations
    }

    if (d == 0) {
        return (pow(x, a) * exp(-x))/upper_incomplete_gamma2(a, x, d++);
    } else if ((d % 2) == 1) {
        double m = 1.0 + ((d -1.0)/2.0);
        return x + ((m - a)/upper_incomplete_gamma2(a, x, d++));
    } else {
        double m = d*2-1;
        return (m-a) + x + (d*(a-d)/upper_incomplete_gamma2(a, x, d++));
    }
}

double lower_incomplete_gamma(double a, double x, int d = 0, int iterations = 100) {
    if (d == iterations) {
        if ((d % 2) == 1) {
            return 1.0; // end iterations
        }
    } else {
        double m = d/2;
        return x + (m - a);
    }

    if (d == 0) {
        return pow(x, a) * exp(-x)/lower_incomplete_gamma(a, x, d++);
    } else if ((d % 2) == 1) {
        double m = d - 1;
        double n = (d - 1)/2;
        return a + m - (a + n) * x/lower_incomplete_gamma(a, x, d++);
    } else {
        double m = d - 1;
        double n = d/2;
        return a + m + n * x/lower_incomplete_gamma(a, x, d++);
    }

}

double lower_incomplete_gamma2(double a, double x) {
    return exp(lgamma(a)) - upper_incomplete_gamma2(a, x);
}

double complimentary_incomplete_gamma(double a, double x) {
    return 1 - upper_incomplete_gamma(a, x);
}
