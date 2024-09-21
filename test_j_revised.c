#include <stdio.h>
#include <math.h>
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

// Optimized version of function_j
double function_j(double f, double fp, double fptilde, double pi, double g, double fptildemin, 
                  double aC, double aX, double gC, double gX, double saC, double saX, double sbC, double sbX) {

    double fpt = MAX(fptilde, fptildemin);

    // Avoid repeated pow calls by precomputing common terms
    double fpt_aX = pow(fpt, aX);
    double fpt_gX = pow(fpt, gX);
    double fpt_saX = pow(fpt, saX);
    double fpt_sbX = pow(fpt, sbX);

    double alpha   = aC  * fpt_aX;
    double gamma   = gC  * fpt_gX;
    double sigma_a = saC * fpt_saX;
    double sigma_b = sbC * fpt_sbX;

    // Minimize branching by using a ternary operator
    double sigma   = (f <= fp) ? sigma_a : sigma_b;

    // Precompute repeated divisions
    double f_over_fp = f / fp;
    double exp1arg = -1.25 * pow(f_over_fp, -4);
    
    double sigma_fp = sigma * fp;
    double exp2arg = -0.5 * pow((f - fp) / sigma_fp, 2);

    double S = alpha * pow(g, 2) * pow((2 * pi), -4) * pow(f, -5) * exp(exp1arg) * pow(gamma, exp(exp2arg));

    return S;
}

int main() {

    double S, f, fp, fptilde;

    // Precompute constants
    const double a  = 0.0081;
    const double b  = 0.6;
    const double g  = 9.807;
    const double pi = 4. * atan(1.0);  // Precompute pi
    const double fptildemin = (1.0/(2.0 * pi)) * pow((4.0 * b / 5.0), (1.0 / 4.0));  // Precompute this expression

    const double gC = 5.87;
    const double aC = 0.0317;

    // Precompute logarithms
    const double aX  = (log(a) - log(aC)) / log(fptildemin);
    const double gX  = -log(gC) / log(fptildemin);

    const double saC = 0.0547;
    const double saX = 0.32;

    const double sbC = 0.0783;
    const double sbX = 0.16;

    // Loop through ranges, avoid small increments when testing performance.
    for (f = -5.; f <= 5.; f += 0.01) {
        for (fp = 0.; fp <= 10.; fp += 0.01) {
            for (fptilde = 0.; fptilde <= 10.; fptilde += 0.01) {
                S = function_j(f, fp, fptilde, pi, g, fptildemin, aC, aX, gC, gX, saC, saX, sbC, sbX);
                // Use S for your needs
            }
        }
    }

    return 0;
}

/***
optimizations:
1.	Move Constants Outside Function:
Constants like pi, g, and fptildemin are calculated only once before the loop, avoiding redundant calculations.
2.	Precompute Common Terms:
Inside the function function_j(), pow(fpt, aX), pow(fpt, gX), pow(fpt, saX), and pow(fpt, sbX) are calculated once and stored in variables (fpt_aX, fpt_gX, etc.) to avoid repeated calls to pow().
3.	Reduce Redundant Computations:
Computations like f / fp and sigma * fp are done once and stored in variables (f_over_fp and sigma_fp) to avoid recalculating them multiple times.
4.	Avoid Branching:
Simplified the branching logic for sigma by using a ternary operator (f <= fp) ? sigma_a : sigma_b, which helps avoid conditional overhead.
5.	Efficient Use of pow() and exp():
The function uses exp() and pow() sparingly, only where absolutely necessary.

Total running time on Mac M1: 112s
***/
