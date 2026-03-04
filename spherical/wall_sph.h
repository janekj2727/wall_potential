#include <math.h>
#include <complex.h>
#include <gsl/gsl_poly.h>

/*
 * wall_sph
 * External potential inducing constant density profile
 * - Spherical geometry
 * - External potential and direct correlation functions
 * - Hard-sphere fluid or hard-sphere fluid with cut LJ attractions (mean-field)
 * - Classical DFT using rhoosenfeld's Fundamental Measure Theory (FMT)
 *
 * Author: Jiří Janek
 * Date: January 2026
 *
 * Compilation: -lm -lgsl -lgslcblas needed to link against gsl_poly
 */

#define M_INF 1.0e6 // infinity for external potential
/* if not defined by math.h, then define pi */
#ifndef M_PI        
  #define M_PI 3.141592653589793238462643383279
#endif

void roots4a(double eta, double rho, double complex *t)
/*
 * Roots of the quartic equation p_4^a (needed for summation)
 *     3*etab*rho^2*(rho+1)^2
 *       -2*[1+etab*rho*(4*rho^2 + 6*rho + 3)]t
 *       +3*etab(2*rho^2 + 2*rho + 1)*t^2 - etab*t^4 = 0
 *     External library GSL used to calculate roots of quartic equation, though exact solution theoretically available
 *
 * parameters:
 *     eta  ... bulk packing fraction
 *     rho  ... radius of the sphere (spherical wall)
 *
 * output:
 *     t(4) ... four complex roots of the equation
 */
{
    double coeffs[5] = {3.0 * eta * pow(rho * (rho + 1.0), 2.0), -2.0 * (1.0 + eta * rho * (rho * (4.0 * rho + 6.0) + 3.0)), 3.0 * eta * (1.0 + 2.0 * rho * (rho + 1.0)), 0.0, -eta}; // coefficients of P_4^a
    double z[8];                                                                                                                                                                      // gsl stores complex roots as interleaved real and imag parts

    gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc(5);
    gsl_poly_complex_solve(coeffs, 5, w, z);
    for (size_t i = 0; i < 4; ++i)
    {
        t[i] = z[2 * i] + I * z[2 * i + 1];
    }
    gsl_poly_complex_workspace_free(w);
}

void roots4b(double eta, double rho, double complex *t)
/*
 * Roots of the quartic equation P_4^b (needed for summation)
 *     2*rho + 2*t - 3*etab*(2*rho+1)^2*t^2
 *       + 4*rho*etab*t^3 + etab*t^4 = 0
 *     External library GSL used to calculate roots of quartic equation, though exact solution theoretically available
 *
 * Parameters:
 *     etab ... bulk packing fraction
 *     rho  ... radius of the sphere (spherical wall)
 *
 * Output:
 *     t(4) ... four complex roots of the equation
 */
{
    double coeffs[5] = {2.0 * rho, 2.0, -3.0 * eta * (2.0 * rho + 1.0), 4.0 * eta * rho, eta}; // coefficients of P_4^b
    double z[8];                                                                               // gsl stores complex roots as interleaved real and imag parts

    gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc(5);
    gsl_poly_complex_solve(coeffs, 5, w, z);
    for (size_t i = 0; i < 4; ++i)
    {
        t[i] = z[2 * i] + I * z[2 * i + 1];
    }
    gsl_poly_complex_workspace_free(w);
}

double c0(double r, double rho, double eta)
/*
 * Direct correlation function (part) c^(1)_0 in spherical geometry.
 *
 *    input:
 *         r    ... distance from the sphere center
 *         rho  ... radius of the sphere (spherical wall)
 *         eta  ... bulk packing fraction
 *
 *    output:
 *         Direct correlation function c^(1)_0(r) for given rho and eta
 */
{
    double complex t[4];
    roots4b(eta, rho, t);

    /* region 1: r < rho */
    if (r < rho)
    {
        return 0.0;
    }
    /* region 2: rho <= r < rho + 1/2 */
    else if (r >= rho && r < rho + 0.5)
    {

        /* prefactor part: (2r+2rho+1)(2r-2rho+1)/(16r) * [2 ln(..) - 3] */
        double logp1 =
            ((4.0 * r * r + 8.0 * r * rho - 12.0 * rho * rho - 20.0 * rho + 4.0 * r - 11.0) * eta * (2.0 * r - 2.0 * rho + 1.0) * (2.0 * r - 2.0 * rho + 1.0)) / (16.0 * (2.0 * r + 1.0)) + 1.0;

        double nons =
            ((2.0 * r + 2.0 * rho + 1.0) * (2.0 * r - 2.0 * rho + 1.0)) / (16.0 * r) * (2.0 * log(logp1) - 3.0);

        /* rho^2/(2r) * ln(2 rho / (2r+1)) */
        nons += (rho * rho) / (2.0 * r) * log(2.0 * rho / (2.0 * r + 1.0));

        /* Summation term: -sum(...) */
        double complex sum = 0.0 + 0.0 * I;

        for (int k = 0; k < 4; k++)
        {

            double complex tk = t[k];

            double complex znum =
                ((4.0 * rho * rho + 6.0 * rho + 3.0) * eta) * tk * tk * tk - 3.0 * tk * tk - 6.0 * rho * tk - 4.0 * rho * rho;

            double complex zden =
                2.0 * ((2.0 * tk * tk + 6.0 * (tk - 1.0) * rho - 3.0) * eta * r * tk) + 2.0 * r;

            double complex zlogarg =
                (2.0 * r - 2.0 * tk - 2.0 * rho + 1.0) / (-2.0 * tk);

            sum += znum / zden * clog(zlogarg);
        }

        return nons - creal(sum);
    }
    /* region 3: rho + 1/2 <= r < rho + 3/2 */
    else if (r >= rho + 0.5 && r < rho + 1.5)
    {

        /* ----- Term 1 */
        double logp1 = ((4.0 * r * r - 12.0 * rho * rho + 8.0 * r * rho - 4.0 * r - 28.0 * rho - 11.0) * eta * (2.0 * r - 2.0 * rho - 1.0) * (2.0 * r - 2.0 * rho - 1.0)) / (16.0 * (2.0 * r - 1.0)) + 1.0;

        double nons = -((2.0 * r + 2.0 * rho - 1.0) * (2.0 * r - 2.0 * rho - 1.0)) / (8.0 * r) * log(logp1);

        /* ----- Term 2 */
        nons += (3.0 * (2.0 * r + 2.0 * rho + 1.0) * (2.0 * r - 2.0 * rho - 3.0)) / (16.0 * r);

        /* ----- Term 3 */
        nons += (rho * rho) / (2.0 * r) * log((2.0 * r - 1.0) / (2.0 * rho + 2.0));

        /* Add ln(1-eta) term */
        nons += ((2.0 * r + 2.0 * rho + 1.0) * (2.0 * r - 2.0 * rho + 1.0)) / (8.0 * r) * log(1.0 - eta);

        /* ----- Summation term */
        double complex sum = 0.0 + 0.0 * I;

        for (int k = 0; k < 4; k++)
        {

            double complex tk = t[k];

            double complex znum =
                ((4.0 * rho * rho + 6.0 * rho + 3.0) * eta) * tk * tk * tk - 3.0 * tk * tk - 6.0 * rho * tk - 4.0 * rho * rho;

            double complex zden =
                2.0 * ((2.0 * tk * tk + 6.0 * (tk - 1.0) * rho - 3.0) * eta * r * tk) + 2.0 * r;

            double complex zlogarg =
                (2.0 * r - 2.0 * tk - 2.0 * rho - 1.0) / (2.0 - 2.0 * tk);

            sum += znum / zden * clog(zlogarg);
        }

        return nons + creal(sum);
    }
    /* region 4: r >= rho + 3/2  */
    else if (r >= rho + 1.5)
    {
        return log(1.0 - eta);
    }
    /* Fallback (should not be reached) */
    else
    {
        return 0.0;
    }

    return 0.0;
}

double c1(double r, double rho, double eta)
/*
 * Direct correlation function (part) c^(1)_1 in spherical geometry.
 *
 *    input:
 *         r    ... distance from the sphere center
 *         rho  ... radius of the sphere (spherical wall)
 *         eta  ... bulk packing fraction
 *
 *    output: value of c_1(r) for given rho and eta
 */
{
    double complex t[4];
    roots4b(eta, rho, t);

    /* region 1: r < rho */
    if (r < rho)
    {
        return 0.0;
    }
    /* region 2: rho <= r < rho + 1/2 */
    else if (r >= rho && r < rho + 0.5)
    {
        double complex sum = 0.0 + 0.0 * I;
        for (int k = 0; k < 4; k++)
        {
            double complex tk = t[k];

            // numerator
            double complex num = (tk + 2.0 * rho + 1.0) * (tk + rho) * tk;

            // denominator
            double complex den = ((2.0 * tk * tk + 6.0 * (tk - 1.0) * rho - 3.0) * eta * r * tk + r);

            // S(t)
            double complex St = num / den;

            // logarithm argument
            double complex logarg = (2.0 * r - 2.0 * tk - 2.0 * rho + 1.0) / (-2.0 * tk);

            // accumulate real part of S(t) * log(...)
            sum += St * clog(logarg);
        }

        // final result
        return -1.5 * eta * creal(sum);
    }
    /* region 3: rho + 1/2 <= r < rho + 3/2 */
    else if (r >= rho + 0.5 && r < rho + 1.5)
    {
        // non-summation term
        double nons = 3.0 * (2.0 * r + 2.0 * rho + 3.0) * (2.0 * r - 2.0 * rho - 1.0) * eta / (8.0 * (eta - 1.0) * r);

        double complex sum = 0.0 + 0.0 * I;
        for (int k = 0; k < 4; k++)
        {
            double complex tk = t[k];

            // numerator
            double complex num = (tk + 2.0 * rho + 1.0) * (tk + rho) * tk;

            // denominator
            double complex den = ((2.0 * tk * tk + 6.0 * (tk - 1.0) * rho - 3.0) * eta * r * tk + r);

            double complex St = num / den;

            // logarithm argument
            double complex logarg = (2.0 * r - 2.0 * tk - 2.0 * rho - 1.0) / (2.0 - 2.0 * tk);

            sum += St * clog(logarg);
        }

        return nons + 1.5 * eta * creal(sum);
    }
    /* region 4: r >= rho + 3/2 */
    else if (r >= rho + 1.5)
    {
        return 3.0 * eta / (eta - 1.0);
    }
    /* fallback (should not be reached) */
    else
    {
        return 0.0;
    }

    return 0.0;
}

/******** c_2 calculation ****************/
double complex Q2(double complex t, double rho, double eta)
/*
 * Auxiliary function – numerator for the sum part in c^(1)_2 calculation
 *
 *     input:
 *        t ... one of the roots of the quartic polynomial P_4^a
 *      rho ... radius of the sphere (spherical wall)
 *     etab ... bulk packing fraction
 *
 *     output: value of Q_2(t, rho, eta)
 */
{
    double complex teta1, teta2, teta3;
    double rho2, rho3, rho4, rho5;
    double eta2, eta3;

    /* Powers of rho */
    rho2 = rho * rho;
    rho3 = rho2 * rho;
    rho4 = rho2 * rho2;
    rho5 = rho3 * rho2;

    /* Powers of eta */
    eta2 = eta * eta;
    eta3 = eta2 * eta;

    /* --- term eta^1 ------------------------- */
    teta1 =
        2.0 * (8.0 * (4.0 * rho2 + 6.0 * rho + 3.0) * t * t + (10.0 * rho2 + 15.0 * rho + 6.0) * (2.0 * rho + 3.0) * t - 52.0 * rho4 - 156.0 * rho3 - 177.0 * rho2 - 90.0 * rho - 18.0) * (rho - t) * rho * eta3;

    /* --- term eta^2 ------------------------- */
    teta2 =
        (8.0 * (8.0 * rho3 + 12.0 * rho2 + 6.0 * rho - 1.0) * t * t * t + 2.0 * (48.0 * rho3 + 72.0 * rho2 + 36.0 * rho - 3.0) * t * t - 6.0 * (16.0 * rho5 + 40.0 * rho4 + 28.0 * rho3 - 4.0 * rho2 - 10.0 * rho - 1.0) * t + (8.0 * rho3 + 12.0 * rho2 + 6.0 * rho - 3.0) * (4.0 * rho2 + 6.0 * rho + 3.0) * rho) * eta2;

    /* --- term eta^3 ------------------------- */
    teta3 =
        (8.0 * t * t * t + 12.0 * t * t - 6.0 * (2.0 * rho2 + 2.0 * rho - 1.0) * t + 3.0 * (4.0 * rho2 + 6.0 * rho + 3.0) * rho - 1.0) * eta;

    return teta1 + teta2 + teta3 + 1.0;
}

double S2(double r, double rho, double eta)
/*
 * Auxiliary function – numerator for the non-sum part in c^(1)_2 calculation
 *     spherical geometry
 *     for rho < r < rho + 1/2
 *
 *     input: r ... position (from the sphere center)
 *          rho ... radius of the sphere (spherical wall)
 *          eta ... bulk packing fraction
 *
 *     output: value of S_2(r, rho, eta)
 */
{
    double teta1, teta2, teta3;
    double rho2, rho3, rho4, rho5, rho6;
    double r2, r3;
    double eta2, eta3;
    double fac;

    /* Powers */
    rho2 = rho * rho;
    rho3 = rho2 * rho;
    rho4 = rho2 * rho2;
    rho5 = rho3 * rho2;
    rho6 = rho3 * rho3;

    r2 = r * r;
    r3 = r2 * r;

    eta2 = eta * eta;
    eta3 = eta2 * eta;

    fac = 2.0 * r - 2.0 * rho + 1.0;

    /* --- η^3 term ---------------------------------------------------- */
    teta3 =
        8.0 * (12.0 * (2.0 * rho2 + 3.0 * rho + 2.0) * (rho + 1.0) * r2 + 6.0 * (2.0 * rho2 + 3.0 * rho + 2.0) * (4.0 * rho2 + 8.0 * rho + 5.0) * r - 3.0 * (12.0 * rho3 + 40.0 * rho2 + 49.0 * rho + 22.0) * (2.0 * rho + 1.0) * rho + 6.0) * fac * fac * rho * eta3;

    /* --- η^2 term ---------------------------------------------------- */
    teta2 =
        (24.0 * (8.0 * rho3 + 12.0 * rho2 + 8.0 * rho + 3.0) * r3 + 12.0 * (16.0 * rho4 + 48.0 * rho3 + 52.0 * rho2 + 30.0 * rho + 17.0) * r2 - 6.0 * (32.0 * rho5 + 80.0 * rho4 + 248.0 * rho3 + 328.0 * rho2 + 164.0 * rho - 13.0) * r - 3.0 * (64.0 * rho6 + 256.0 * rho5 + 352.0 * rho4 + 376.0 * rho3 + 352.0 * rho2 + 178.0 * rho + 1.0)) * fac * eta2;

    /* --- η term ------------------------------------------------------ */
    teta1 =
        3.0 * (8.0 * r3 + 4.0 * (2.0 * rho + 3.0) * r2 - (8.0 * rho2 + 8.0 * rho + 34.0) * r - 8.0 * rho3 - 20.0 * rho2 - 14.0 * rho - 19.0) * fac * eta;

    return teta1 + teta2 + teta3;
}

double T2(double r, double rho, double eta)
/*
 * Auxiliary function – numerator for the non-sum part in c^(1)_2 calculation
 *      spherical geometry
 *      for rho + 1/2 < r < rho + 3/2
 *
 *      input: r ... position (from the sphere center)
 *           rho ... radius of the sphere (spherical wall)
 *           eta ... bulk packing fraction
 *
 *      output: value of T_2(r, rho, eta)
 */
{
    double teta1, teta2, teta3, teta4;
    double r2, r3, r4, r5;
    double rho2, rho3, rho4, rho5, rho6, rho7, rho8;
    double eta2, eta3, eta4;
    double rhop1, rhop3;

    /* Powers of r */
    r2 = r * r;
    r3 = r2 * r;
    r4 = r3 * r;
    r5 = r4 * r;

    /* Powers of rho */
    rho2 = rho * rho;
    rho3 = rho2 * rho;
    rho4 = rho2 * rho2;
    rho5 = rho3 * rho2;
    rho6 = rho3 * rho3;
    rho7 = rho6 * rho;
    rho8 = rho4 * rho4;

    /* Powers of eta */
    eta2 = eta * eta;
    eta3 = eta2 * eta;
    eta4 = eta3 * eta;

    rhop1 = 2.0 * rho + 1.0;
    rhop3 = 2.0 * rho + 3.0;

    /* --- eta^4 block -------------------------------------------------- */
    teta4 = -(
                384.0 * (4.0 * rho2 + 6.0 * rho + 3.0) * r5 * (r - 1.0) - 96.0 * (56.0 * rho3 + 112.0 * rho2 + 88.0 * rho + 21.0) * rhop1 * r4 + 192.0 * (64.0 * rho5 + 208.0 * rho4 + 288.0 * rho3 + 192.0 * rho2 + 48.0 * rho - 3.0) * r3 + 24.0 * (48.0 * rho4 + 184.0 * rho3 + 184.0 * rho2 + 6.0 * rho - 39.0) * rhop1 * rhop1 * r2 - 24.0 * (128.0 * rho5 + 528.0 * rho4 + 728.0 * rho3 + 280.0 * rho2 - 146.0 * rho - 105.0) * rhop1 * rhop1 * r + 6.0 * (48.0 * rho4 + 136.0 * rho3 + 80.0 * rho2 - 46.0 * rho - 45.0) * rhop3 * rhop1 * rhop1 * rhop1) *
            rho * eta4;

    /* --- eta^3 block -------------------------------------------------- */
    teta3 = -(
                192.0 * (16.0 * rho3 + 24.0 * rho2 + 12.0 * rho + 1.0) * r5 * (r - 1.0) - 48.0 * (448.0 * rho5 + 1216.0 * rho4 + 1472.0 * rho3 + 932.0 * rho2 + 280.0 * rho + 7.0) * r4 + 96.0 * (256.0 * rho6 + 832.0 * rho5 + 1200.0 * rho4 + 992.0 * rho3 + 496.0 * rho2 + 124.0 * rho - 1.0) * r3 + 12.0 * (768.0 * rho7 + 6016.0 * rho6 + 16704.0 * rho5 + 22896.0 * rho4 + 16784.0 * rho3 + 6144.0 * rho2 + 828.0 * rho - 13.0) * r2 - 12.0 * (2048.0 * rho8 + 13568.0 * rho7 + 37888.0 * rho6 + 58048.0 * rho5 + 52880.0 * rho4 + 28176.0 * rho3 + 7688.0 * rho2 + 636.0 * rho - 35.0) * r + 3.0 * (1536.0 * rho8 + 8960.0 * rho7 + 21888.0 * rho6 + 29408.0 * rho5 + 23696.0 * rho4 + 11184.0 * rho3 + 2600.0 * rho2 + 90.0 * rho - 15.0) * rhop3) *
            eta3;

    /* --- eta^2 block -------------------------------------------------- */
    teta2 = -(
                384.0 * (r - 1.0) * r5 - 96.0 * (8.0 * rho3 + 40.0 * rho2 + 40.0 * rho + 21.0) * r4 + 192.0 * (56.0 * rho3 + 88.0 * rho2 + 52.0 * rho + 9.0) * r3 + 24.0 * (64.0 * rho5 + 208.0 * rho4 + 688.0 * rho3 + 976.0 * rho2 + 544.0 * rho + 71.0) * r2 - 24.0 * (448.0 * rho5 + 1840.0 * rho4 + 3184.0 * rho3 + 2640.0 * rho2 + 968.0 * rho + 65.0) * r - 6.0 * (64.0 * rho6 + 32.0 * rho5 - 592.0 * rho4 - 1296.0 * rho3 - 1108.0 * rho2 - 398.0 * rho - 15.0) * rhop3) *
            eta2;

    /* --- eta^1 block -------------------------------------------------- */
    teta1 =
        (96.0 * (r - 10.0) * r3 - 48.0 * (4.0 * rho2 + 4.0 * rho + 15.0) * r2 + 48.0 * (20.0 * rho2 + 44.0 * rho + 37.0) * r + 6.0 * (4.0 * rho2 - 4.0 * rho - 11.0) * rhop3 * rhop3) * eta;

    return teta1 + teta2 + teta3 + teta4;
}

double c2(double r, double rho, double eta)
/*
 *  Direct correlation function (part) c^(1)_2 in spherical geometry
 *
 *  inputs: r    ... distance from the sphere center
 *          rho  ... sphere (spherical wall) radius
 *          eta  ... bulk packing fraction
 *
 *  output: value of c_2(r) for given rho and eta
 */
{
    double complex t[4];
    double den1, den2, logterm;
    double complex sum, cden, lga;
    int k;

    roots4a(eta, rho, t);

    /* Common denominator den1 */
    den1 = 2.0 * (4.0 * rho * rho + 6.0 * rho + 3.0) * eta * rho;
    den1 = (den1 + 1.0) * (eta - 1.0);

    /* region 1: r < rho */
    if (r < rho)
    {
        return 0.0;
    }
    /* region 2: rho <= r < rho + 0.5 */
    else if (r >= rho && r < rho + 0.5)
    {

        den2 = (4.0 * (r * r + r + 2.0 * r * rho - 3.0 * rho * rho - 5.0 * rho) - 11.0) * eta * (2.0 * r - 2.0 * rho + 1.0) * (2.0 * r - 2.0 * rho + 1.0) + 16.0 * (2.0 * r + 1.0);

        logterm = log((2.0 * r + 1.0) / (2.0 * rho)) / (8.0 * r);

        sum = 0.0 + 0.0 * I;

        for (k = 0; k < 4; ++k)
        {
            cden = (2.0 * (t[k] + 2.0 * rho) * (t[k] - rho) - 6.0 * rho - 3.0) * eta * (t[k] - rho) + 1.0;

            lga = (2.0 * r - 2.0 * t[k] + 1.0) / (2.0 * rho - 2.0 * t[k]);

            sum += Q2(t[k], rho, eta) / (8.0 * r * cden * den1) * clog(lga);
        }

        return S2(r, rho, eta) / (8.0 * r * den2 * den1) + logterm + creal(sum);
    }
    /* region 3: rho + 0.5 <= r < rho + 1.5 */
    else if (r >= rho + 0.5 && r < rho + 1.5)
    {
        den2 = (4.0 * (r * r - r + 2.0 * r * rho - 3.0 * rho * rho - 7.0 * rho) - 11.0) * eta * (2.0 * r - 2.0 * rho - 1.0) * (2.0 * r - 2.0 * rho - 1.0) + 16.0 * (2.0 * r - 1.0);

        logterm = -log((2.0 * r - 1.0) / (2.0 * rho + 2.0)) / (8.0 * r);

        sum = 0.0 + 0.0 * I;

        for (k = 0; k < 4; ++k)
        {
            cden = (2.0 * (t[k] + 2.0 * rho) * (t[k] - rho) - 6.0 * rho - 3.0) * eta * (t[k] - rho) + 1.0;

            lga = (2.0 * r - 2.0 * t[k] - 1.0) / (2.0 * (rho - t[k] + 1.0));

            sum += Q2(t[k], rho, eta) / (8.0 * r * cden * den1) * clog(lga);
        }

        return T2(r, rho, eta) / (16.0 * r * den2 * den1 * (eta - 1.0)) + logterm - creal(sum);
    }
    /* region 4: r >= rho + 1.5 */
    else
    {
        return -3.0 * eta * (2.0 + eta) / (2.0 * (1.0 - eta) * (1.0 - eta));
    }

    /* should not reach */
    return 0.0;
}

/******** end of c_2 calculation ********************/

/******** c_3 calculation ***************************/
double complex Q3(double complex t, double r, double rho, double eta)
/*
 * Auxiliary function – numerator for the sum part in c^(1)_3 calculation
 *     spherical geometry
 *
 *     input: t ... one of the roots of the quartic polynomial P_4^a
 *            r ... distance from the sphere center
 *          rho ... radius of the sphere (spherical wall)
 *          eta ... bulk packing fraction
 *
 *     output: complex*16 value of Q_3(t, rho, eta)
 */
{
    double complex termt3, teta5, teta4, teta3, teta2, teta1, teta0;
    double rho2, rho3, rho4, rho5, rho6, rho7, rho8, rho9, rho10, rho11;

    rho2 = rho * rho;
    rho3 = rho2 * rho;
    rho4 = rho2 * rho2;
    rho5 = rho3 * rho2;
    rho6 = rho3 * rho3;
    rho7 = rho3 * rho4;
    rho8 = rho4 * rho4;
    rho9 = rho4 * rho5;
    rho10 = rho5 * rho5;
    rho11 = rho5 * rho6;

    termt3 = 2.0 * ((4.0 * pow(4.0 * rho2 + 6.0 * rho + 3.0, 2) * rho2 * pow(eta, 4) - 4.0 * (4.0 * rho2 + 6.0 * rho + 3.0) * (8.0 * rho3 + 12.0 * rho2 + 6.0 * rho - 1.0) * rho * pow(eta, 3) + (64.0 * rho6 + 192.0 * rho5 + 240.0 * rho4 + 112.0 * rho3 - 12.0 * rho2 - 24.0 * rho + 1.0) * pow(eta, 2) + 2.0 * (8.0 * rho3 + 12.0 * rho2 + 6.0 * rho - 1.0) * eta + 1.0) * (4.0 * pow(2.0 * rho + 1.0, 2) * pow(r, 2) - 96.0 * (rho + 1.0) * rho * r - 48.0 * rho4 - 96.0 * rho3 - 16.0 * rho2 + 32.0 * rho - 1.0) * eta) * (t * t * t);

    teta5 = ((4.0 * (4.0 * (328.0 * rho4 + 948.0 * rho3 + 1122.0 * rho2 + 639.0 * rho + 153.0) * pow(r, 2) + 8.0 * (88.0 * rho5 - 52.0 * rho4 - 590.0 * rho3 - 879.0 * rho2 - 546.0 * rho - 135.0) * r + (3072.0 * rho7 + 22560.0 * rho6 + 61584.0 * rho5 + 94112.0 * rho4 + 90600.0 * rho3 + 56442.0 * rho2 + 21285.0 * rho + 3843.0))) * (rho + 1.0) * rho3 * t * t - 8.0 * (12.0 * (64.0 * rho8 + 344.0 * rho7 + 820.0 * rho6 + 1122.0 * rho5 + 979.0 * rho4 + 576.0 * rho3 + 234.0 * rho2 + 63.0 * rho + 9.0) * pow(r, 2) - 8.0 * (104.0 * rho6 + 596.0 * rho5 + 1198.0 * rho4 + 1143.0 * rho3 + 522.0 * rho2 + 81.0 * rho - 9.0) * (rho + 1.0) * rho * r + 768.0 * rho10 + 11296.0 * rho9 + 50032.0 * rho8 + 115248.0 * rho7 + 161352.0 * rho6 + 144018.0 * rho5 + 80763.0 * rho4 + 25920.0 * rho3 + 3402.0 * rho2 - 189.0 * rho - 27.0) * rho2 * t + 4.0 * (4.0 * (256.0 * rho8 + 1480.0 * rho7 + 3964.0 * rho6 + 6390.0 * rho5 + 6849.0 * rho4 + 5040.0 * rho3 + 2493.0 * rho2 + 756.0 * rho + 108.0) * pow(r, 2) - 24.0 * (56.0 * rho4 + 140.0 * rho3 + 106.0 * rho2 + 21.0 * rho - 6.0) * pow(rho + 1.0, 3) * rho * r + 5664.0 * rho9 + 34736.0 * rho8 + 93488.0 * rho7 + 142856.0 * rho6 + 133218.0 * rho5 + 74799.0 * rho4 + 22176.0 * rho3 + 1503.0 * rho2 - 756.0 * rho - 108.0) * rho3) * pow(eta, 5);

    teta4 = -(-(4.0 * (576.0 * rho6 - 256.0 * rho5 - 3504.0 * rho4 - 5280.0 * rho3 - 3012.0 * rho2 - 396.0 * rho + 207.0) * pow(r, 2) - 8.0 * (576.0 * rho7 + 3520.0 * rho6 + 5360.0 * rho5 + 1600.0 * rho4 - 3308.0 * rho3 - 3052.0 * rho2 - 611.0 * rho + 180.0) * r - (22272.0 * rho8 + 167424.0 * rho7 + 462336.0 * rho6 + 710528.0 * rho5 + 676320.0 * rho4 + 405264.0 * rho3 + 135072.0 * rho2 + 12564.0 * rho - 5121.0)) * (rho + 1.0) * rho2 * t * t - 2.0 * (12.0 * (704.0 * rho9 + 3904.0 * rho8 + 9584.0 * rho7 + 13424.0 * rho6 + 11724.0 * rho5 + 6548.0 * rho4 + 2289.0 * rho3 + 437.0 * rho2 + 12.0 * rho - 12.0) * pow(r, 2) - 8.0 * (576.0 * rho8 + 4160.0 * rho7 + 13040.0 * rho6 + 20992.0 * rho5 + 18348.0 * rho4 + 8088.0 * rho3 + 1173.0 * rho2 - 168.0 * rho + 12.0) * (rho + 1.0) * rho * r + 8448.0 * rho11 + 107776.0 * rho10 + 457216.0 * rho9 + 1025664.0 * rho8 + 1388864.0 * rho7 + 1163552.0 * rho6 + 562344.0 * rho5 + 108240.0 * rho4 - 23427.0 * rho3 - 11823.0 * rho2 - 36.0 * rho + 36.0) * rho * t + (4.0 * (3776.0 * rho9 + 21056.0 * rho8 + 53168.0 * rho7 + 78832.0 * rho6 + 75172.0 * rho5 + 46864.0 * rho4 + 17865.0 * rho3 + 2859.0 * rho2 - 585.0 * rho - 279.0) * pow(r, 2) - 24.0 * (576.0 * rho6 + 2240.0 * rho5 + 3440.0 * rho4 + 2368.0 * rho3 + 596.0 * rho2 - 52.0 * rho + 1.0) * pow(rho + 1.0, 3) * rho * r + 6912.0 * rho11 + 92928.0 * rho10 + 419584.0 * rho9 + 985216.0 * rho8 + 1379680.0 * rho7 + 1191248.0 * rho6 + 603104.0 * rho5 + 139292.0 * rho4 - 8973.0 * rho3 - 8007.0 * rho2 + 585.0 * rho + 279.0) * rho2) * pow(eta, 4);

    teta3 = ((4.0 * (1536.0 * rho6 + 4608.0 * rho5 + 5832.0 * rho4 + 3132.0 * rho3 + 306.0 * rho2 - 279.0 * rho + 17.0) * pow(r, 2) - 8.0 * (1920.0 * rho6 + 5976.0 * rho5 + 7892.0 * rho4 + 4546.0 * rho3 + 647.0 * rho2 - 324.0 * rho + 15.0) * r + 12288.0 * rho9 + 92160.0 * rho8 + 258048.0 * rho7 + 395040.0 * rho6 + 353520.0 * rho5 + 177408.0 * rho4 + 26712.0 * rho3 - 18942.0 * rho2 - 9573.0 * rho + 427.0) * (rho + 1.0) * rho * t * t - ((6144.0 * rho10 + 30720.0 * rho9 + 66048.0 * rho8 + 73536.0 * rho7 + 37728.0 * rho6 - 3504.0 * rho5 - 15672.0 * rho4 - 9000.0 * rho3 - 2712.0 * rho2 - 504.0 * rho + 24.0) * pow(r, 2) + 16.0 * (1224.0 * rho6 + 4260.0 * rho5 + 6258.0 * rho4 + 4335.0 * rho3 + 1199.0 * rho2 - 25.0 * rho + 1.0) * (rho + 1.0) * rho * r + 6144.0 * pow(rho, 12) + 86016.0 * rho11 + 365568.0 * rho10 + 795456.0 * rho9 + 968416.0 * rho8 + 573184.0 * rho7 - 45744.0 * rho6 - 307484.0 * rho5 - 182498.0 * rho4 - 30390.0 * rho3 + 4950.0 * rho2 + 126.0 * rho - 6.0) * t + (4.0 * (1024.0 * rho10 + 5632.0 * rho9 + 14848.0 * rho8 + 22952.0 * rho7 + 21428.0 * rho6 + 10262.0 * rho5 - 577.0 * rho4 - 4018.0 * rho3 - 2233.0 * rho2 - 378.0 * rho + 60.0) * pow(r, 2) + 24.0 * (216.0 * rho3 + 388.0 * rho2 + 226.0 * rho - 1.0) * pow(rho + 1.0, 3) * rho2 * r + 18432.0 * rho11 + 109568.0 * rho10 + 283808.0 * rho9 + 400592.0 * rho8 + 304432.0 * rho7 + 79480.0 * rho6 - 52286.0 * rho5 - 43487.0 * rho4 - 5918.0 * rho3 + 2653.0 * rho2 + 378.0 * rho - 60.0) * rho) * pow(eta, 3);

    teta2 = ((4.0 * (384.0 * rho3 + 576.0 * rho2 + 288.0 * rho - 23.0) * pow(r, 2) - 24.0 * (160.0 * rho3 + 246.0 * rho2 + 126.0 * rho - 9.0) * r + 6144.0 * rho6 + 27648.0 * rho5 + 46080.0 * rho4 + 41088.0 * rho3 + 19212.0 * rho2 + 3948.0 * rho - 853.0) * (rho + 1.0) * rho * t * t + ((-1536.0 * rho7 - 5376.0 * rho6 - 7296.0 * rho5 - 4272.0 * rho4 - 864.0 * rho3 - 24.0 * rho2 - 72.0 * rho + 48.0) * pow(r, 2) - 16.0 * (288.0 * rho3 + 481.0 * rho2 + 265.0 * rho - 1.0) * (rho + 1.0) * rho * r - 1536.0 * rho9 - 31488.0 * rho8 - 110592.0 * rho7 - 176112.0 * rho6 - 138992.0 * rho5 - 38948.0 * rho4 + 13896.0 * rho3 + 9038.0 * rho2 - 262.0 * rho - 12.0) * t + ((2048.0 * rho8 + 8192.0 * rho7 + 15872.0 * rho6 + 17792.0 * rho5 + 11804.0 * rho4 + 3608.0 * rho3 - 308.0 * rho2 - 464.0 * rho + 16.0) * pow(r, 2) + 432.0 * pow(rho + 1.0, 3) * rho3 * r + (1152.0 * rho7 + 3904.0 * rho6 + 4800.0 * rho5 + 1929.0 * rho4 - 614.0 * rho3 - 435.0 * rho2 + 132.0 * rho - 4.0) * pow(2.0 * rho + 1.0, 2))) * pow(eta, 2);

    teta1 = (48.0 * (2.0 * r * r - 5.0 * r + 20.0 * rho3 + 42.0 * rho2 + 27.0 * rho + 9.0) * (rho + 1.0) * rho * t * t - 2.0 * (12.0 * (4.0 * rho4 + 8.0 * rho3 + 4.0 * rho2 + 1.0) * pow(r, 2) + 144.0 * (rho + 1.0) * rho * r + 48.0 * rho6 + 2064.0 * rho5 + 4932.0 * rho4 + 4344.0 * rho3 + 1120.0 * rho2 - 308.0 * rho - 3.0) * t + (16.0 * (20.0 * rho3 + 30.0 * rho2 + 15.0 * rho - 2.0) * (rho2 + rho + 1.0) * pow(r, 2) + 288.0 * rho6 + 784.0 * rho5 + 664.0 * rho4 + 28.0 * rho3 - 172.0 * rho2 - 52.0 * rho + 8.0)) * eta;

    teta0 = 48.0 * (rho + 1.0) * rho * t * (t - 4.0) + 4.0 * (rho2 + rho + 1.0) * (2.0 * r + 1.0) * (2.0 * r - 1.0);

    return termt3 + teta5 + teta4 + teta3 + teta2 + teta1 + teta0;
}

double S3(double r, double rho, double eta)
/*
 * Auxiliary function – numerator for the non-sum part in c^(1)_3 calculation
 *     spherical geometry
 *     rho < r < rho + 1/2
 *
 *     input: r ... position (from the sphere center)
 *          rho ... radius of the sphere (spherical wall)
 *          eta ... bulk packing fraction
 *
 *     output: value of S_3(r, rho, eta)
 */
{
    double pre, teta6, teta5, teta4, teta3, teta2, teta1, teta0;
    double rho2, rho3, rho4, rho5, rho6, rho7, rho8, rho9, rho10, rho11;

    rho2 = rho * rho;
    rho3 = rho2 * rho;
    rho4 = rho2 * rho2;
    rho5 = rho3 * rho2;
    rho6 = rho3 * rho3;
    rho7 = rho3 * rho4;
    rho8 = rho4 * rho4;
    rho9 = rho4 * rho5;
    rho10 = rho5 * rho5;
    rho11 = rho5 * rho6;

    pre = (2.0 * r - 2.0 * rho + 1.0);

    teta6 = -16.0 * (4.0 * r * (r + 2.0 * rho + 1.0) - 12.0 * rho2 - 20.0 * rho - 11.0) * pow(4.0 * rho2 + 6.0 * rho + 3.0, 2) * pow(2.0 * r - 2.0 * rho + 1.0, 2) * (2.0 * r - 2.0 * rho - 1.0) * (2.0 * rho + 1.0) * rho2 * pow(eta, 6);

    teta5 = 4.0 * (32.0 * (320.0 * rho6 + 976.0 * rho5 + 1200.0 * rho4 + 636.0 * rho3 + 72.0 * rho2 - 54.0 * rho - 9.0) * pow(r, 5) - 16.0 * (640.0 * rho7 + 1376.0 * rho6 - 1136.0 * rho5 - 6440.0 * rho4 - 8236.0 * rho3 - 4728.0 * rho2 - 1056.0 * rho + 27.0) * pow(r, 4) - 16.0 * (3840.0 * rho8 + 16832.0 * rho7 + 32384.0 * rho6 + 30592.0 * rho5 + 9168.0 * rho4 - 8748.0 * rho3 - 9072.0 * rho2 - 2718.0 * rho - 45.0) * pow(r, 3) + 8.0 * (17920.0 * rho9 + 78464.0 * rho8 + 128832.0 * rho7 + 66112.0 * rho6 - 76896.0 * rho5 - 146856.0 * rho4 - 101124.0 * rho3 - 33336.0 * rho2 - 4176.0 * rho + 63.0) * pow(r, 2) - 2.0 * (56320.0 * rho10 + 257792.0 * rho9 + 385792.0 * rho8 + 46016.0 * rho7 - 546624.0 * rho6 - 686224.0 * rho5 - 298304.0 * rho4 + 48980.0 * rho3 + 91896.0 * rho2 + 26214.0 * rho + 81.0) * r + (30720.0 * rho11 + 147968.0 * rho10 + 202496.0 * rho9 - 70784.0 * rho8 - 448064.0 * rho7 - 374048.0 * rho6 + 47440.0 * rho5 + 230808.0 * rho4 + 99516.0 * rho3 - 12024.0 * rho2 - 13536.0 * rho - 99.0)) * rho * pow(eta, 5);

    teta4 = (32.0 * (384.0 * rho7 + 3456.0 * rho6 + 8832.0 * rho5 + 10992.0 * rho4 + 6864.0 * rho3 + 1812.0 * rho2 - 28.0 * rho - 5.0) * pow(r, 5) - 16.0 * (768.0 * rho8 + 4992.0 * rho7 + 20864.0 * rho6 + 41888.0 * rho5 + 43600.0 * rho4 + 21352.0 * rho3 + 2220.0 * rho2 - 1482.0 * rho + 15.0) * pow(r, 4) - 16.0 * (4608.0 * rho9 + 47616.0 * rho8 + 161664.0 * rho7 + 317248.0 * rho6 + 393792.0 * rho5 + 310176.0 * rho4 + 138336.0 * rho3 + 23792.0 * rho2 - 3220.0 * rho - 25.0) * pow(r, 3) + 8.0 * (21504.0 * rho10 + 207360.0 * rho9 + 862464.0 * rho8 + 1933824.0 * rho7 + 2565824.0 * rho6 + 2048256.0 * rho5 + 922336.0 * rho4 + 177152.0 * rho3 - 10920.0 * rho2 - 5474.0 * rho + 35.0) * pow(r, 2) - 2.0 * (67584.0 * rho11 + 632832.0 * rho10 + 3058688.0 * rho9 + 7361792.0 * rho8 + 9771392.0 * rho7 + 6720960.0 * rho6 + 1059264.0 * rho5 - 1829952.0 * rho4 - 1292528.0 * rho3 - 243324.0 * rho2 + 33516.0 * rho + 45.0) * r + (36864.0 * pow(rho, 12) + 337920.0 * rho11 + 1855488.0 * rho10 + 4647424.0 * rho9 + 5665792.0 * rho8 + 2344192.0 * rho7 - 1824448.0 * rho6 - 2218048.0 * rho5 - 203296.0 * rho4 + 609032.0 * rho3 + 221516.0 * rho2 - 17638.0 * rho - 55.0)) * pow(eta, 4);

    teta3 = -2.0 * (32.0 * (192.0 * rho6 + 576.0 * rho5 + 672.0 * rho4 + 36.0 * rho3 - 444.0 * rho2 - 281.0 * rho - 1.0) * pow(r, 5) + 16.0 * (384.0 * rho7 - 960.0 * rho6 - 4800.0 * rho5 - 6648.0 * rho4 - 2852.0 * rho3 + 834.0 * rho2 + 885.0 * rho - 67.0) * pow(r, 4) - 16.0 * (2304.0 * rho8 + 8448.0 * rho7 + 19008.0 * rho6 + 24432.0 * rho5 + 15792.0 * rho4 - 2544.0 * rho3 - 9836.0 * rho2 - 4817.0 * rho + 119.0) * pow(r, 3) + 8.0 * (1536.0 * rho9 + 39168.0 * rho8 + 138624.0 * rho7 + 203232.0 * rho6 + 112560.0 * rho5 - 40720.0 * rho4 - 83096.0 * rho3 - 34818.0 * rho2 - 2509.0 * rho + 243.0) * pow(r, 2) + 2.0 * (15360.0 * rho10 - 125952.0 * rho9 - 577536.0 * rho8 - 904896.0 * rho7 - 322816.0 * rho6 + 722512.0 * rho5 + 1019792.0 * rho4 + 459692.0 * rho3 + 2836.0 * rho2 - 49481.0 * rho + 1359.0) * r - (18432.0 * rho11 - 70656.0 * rho10 - 405504.0 * rho9 - 649344.0 * rho8 - 275520.0 * rho7 + 141664.0 * rho6 - 195376.0 * rho5 - 712936.0 * rho4 - 530084.0 * rho3 - 94034.0 * rho2 + 35995.0 * rho - 701.0)) * pow(eta, 3);

    teta2 = -(96.0 * (32.0 * rho3 + 48.0 * rho2 + 22.0 * rho - 15.0) * pow(r, 5) + 16.0 * (192.0 * rho4 - 768.0 * rho3 - 1428.0 * rho2 - 768.0 * rho + 149.0) * pow(r, 4) - 16.0 * (1152.0 * rho5 + 2496.0 * rho4 + 5112.0 * rho3 + 4692.0 * rho2 + 1722.0 * rho - 785.0) * pow(r, 3) + 8.0 * (3840.0 * rho6 + 27648.0 * rho5 + 53616.0 * rho4 + 39264.0 * rho3 + 4248.0 * rho2 - 6312.0 * rho - 577.0) * pow(r, 2) + 2.0 * (19968.0 * rho7 - 123648.0 * rho6 - 397920.0 * rho5 - 439056.0 * rho4 - 100112.0 * rho3 + 118600.0 * rho2 + 72082.0 * rho - 8117.0) * r - (9216.0 * rho8 - 61440.0 * rho7 - 82368.0 * rho6 + 136704.0 * rho5 + 360912.0 * rho4 + 171200.0 * rho3 - 84340.0 * rho2 - 82096.0 * rho + 5587.0)) * pow(eta, 2);

    teta1 = -2.0 * (96.0 * pow(r, 5) + 48.0 * (2.0 * rho - 11.0) * pow(r, 4) - 48.0 * (12.0 * rho2 + 8.0 * rho + 33.0) * pow(r, 3) + 24.0 * (136.0 * rho3 + 372.0 * rho2 + 262.0 * rho - 29.0) * pow(r, 2) + 2.0 * (1776.0 * rho4 - 11136.0 * rho3 - 16968.0 * rho2 - 7824.0 * rho + 3043.0) * r - (288.0 * rho5 - 3504.0 * rho4 + 7632.0 * rho3 + 16968.0 * rho2 + 9450.0 * rho - 3055.0)) * eta;

    teta0 = -96.0 * (2.0 * r + 2.0 * rho - 15.0) * (2.0 * r + 1.0);

    return pre * (teta6 + teta5 + teta4 + teta3 + teta2 + teta1 + teta0);
}

double T3(double r, double rho, double eta)
/*
 * Auxiliary function – numerator for the non-sum part in c^(1)_3 calculation
 *      spherical geometry
 *      rho + 1/2 < r < rho + 3/2
 *
 *      input: r ... position (from the sphere center)
 *           rho ... radius of the sphere (spherical wall)
 *           eta ... bulk packing fraction
 *
 *      output: value of T_3(r, rho, eta)
 */
{
    double teta6, teta5, teta4, teta3, teta2, teta1, teta0;
    double rho2, rho3, rho4, rho5, rho6, rho7, rho8, rho9, rho10, rho11, rho12;

    rho2 = rho * rho;
    rho3 = rho2 * rho;
    rho4 = rho2 * rho2;
    rho5 = rho3 * rho2;
    rho6 = rho3 * rho3;
    rho7 = rho3 * rho4;
    rho8 = rho4 * rho4;
    rho9 = rho4 * rho5;
    rho10 = rho5 * rho5;
    rho11 = rho5 * rho6;
    rho12 = rho6 * rho6;

    teta6 = 16.0 * (
        64.0 * (2.0 - r) * pow(4.0 * rho2 + 6.0 * rho + 3.0, 2) * pow(r,7)
        + 16.0 * (768.0 * rho6 + 3200.0 * rho5 + 6256.0 * rho4
            + 7248.0 * rho3 + 5232.0 * rho2 + 2196.0 * rho + 423.0) * pow(r,6)
        - 32.0 * (512.0 * rho7 + 3008.0 * rho6 + 7760.0 * rho5
            + 12536.0 * rho4 + 13924.0 * rho3 + 10446.0 * rho2
            + 4710.0 * rho + 981.0) * pow(r,5)
        - 4.0 * (7680.0 * rho8 + 35840.0 * rho7 + 83264.0 * rho6
            + 117888.0 * rho5 + 88080.0 * rho4 + 1632.0 * rho3
            - 55692.0 * rho2 - 42768.0 * rho - 11241.0) * pow(r,4)
        + 8.0 * (12288.0 * rho9 + 76544.0 * rho8 + 226688.0 * rho7
            + 434176.0 * rho6 + 587328.0 * rho5 + 559456.0 * rho4
            + 356376.0 * rho3 + 137616.0 * rho2 + 25848.0 * rho + 855.0) * pow(r,3)
        - (102400.0 * rho10 + 739328.0 * rho9 + 2481920.0 * rho8
            + 5415680.0 * rho7 + 8795008.0 * rho6 + 10944192.0 * rho5
            + 10115488.0 * rho4 + 6578384.0 * rho3 + 2805000.0 * rho2
            + 696876.0 * rho + 75843.0) * pow(r,2)
        + 2.0 * (6144.0 * rho9 + 43776.0 * rho8 + 140224.0 * rho7
            + 295328.0 * rho6 + 494784.0 * rho5 + 667552.0 * rho4
            + 660492.0 * rho3 + 432246.0 * rho2 + 165654.0 * rho + 28341.0)
            * pow(2.0 * rho + 1.0, 2) * r
        - (1152.0 * rho9 + 7488.0 * rho8 + 21376.0 * rho7
            + 42816.0 * rho6 + 77096.0 * rho5 + 113852.0 * rho4
            + 116712.0 * rho3 + 74748.0 * rho2 + 26793.0 * rho
            + 4059.0) * (2.0 * rho + 3.0) * pow(2.0 * rho + 1.0, 2)
    ) * rho2 * pow(eta,6);

    /* teta5, teta4, teta3, teta2, teta1, teta0 follow the Fortran expressions */
    teta5 = 2.0 * (
        512.0 * (2.0 - r) * (4.0 * rho2 + 6.0 * rho + 3.0) * (4.0 * rho2 + 2.0 * rho + 1.0) * (rho + 1.0) * pow(r,7)
        + 64.0 * (1536.0 * rho7 + 7232.0 * rho6 + 12832.0 * rho5 + 10704.0 * rho4
            + 3552.0 * rho3 + 72.0 * rho2 + 204.0 * rho + 285.0) * pow(r,6)
        - 64.0 * (2048.0 * rho8 + 13696.0 * rho7 + 36352.0 * rho6 + 40384.0 * rho5
            + 7456.0 * rho4 - 24504.0 * rho3 - 20208.0 * rho2 - 3514.0 * rho + 1317.0) * pow(r,5)
        - 16.0 * (15360.0 * rho9 + 88320.0 * rho8 + 167552.0 * rho7 + 141056.0 * rho6
            + 149152.0 * rho5 + 336624.0 * rho4 + 443824.0 * rho3 + 266444.0 * rho2
            + 53604.0 * rho - 7491.0) * pow(r,4)
        + 32.0 * (24576.0 * rho10 + 186368.0 * rho9 + 548864.0 * rho8 + 743680.0 * rho7
            + 332288.0 * rho6 - 312928.0 * rho5 - 459552.0 * rho4 - 172392.0 * rho3
            + 21564.0 * rho2 + 22860.0 * rho + 603.0) * pow(r,3)
        - 4.0 * (204800.0 * rho11 + 1811456.0 * rho10 + 6457856.0 * rho9
            + 10749696.0 * rho8 + 4965632.0 * rho7 - 12269760.0 * rho6
            - 25457312.0 * rho5 - 22228640.0 * rho4 - 10114944.0 * rho3
            - 2004320.0 * rho2 + 48548.0 * rho + 50661.0) * pow(r,2)
        + 4.0 * (98304.0 * rho12 + 985088.0 * rho11 + 4075520.0 * rho10
            + 8044544.0 * rho9 + 4541952.0 * rho8 - 13309440.0 * rho7
            - 35129088.0 * rho6 - 40819168.0 * rho5 - 27170736.0 * rho4
            - 10160248.0 * rho3 - 1636536.0 * rho2 + 81714.0 * rho + 37731.0) * r
        - (36864.0 * rho12 + 356352.0 * rho11 + 1398784.0 * rho10
            + 2332672.0 * rho9 - 331264.0 * rho8 - 8488960.0 * rho7
            - 16720576.0 * rho6 - 17258112.0 * rho5 - 10497696.0 * rho4
            - 3555024.0 * rho3 - 473100.0 * rho2 + 46152.0 * rho + 10791.0) * (2.0 * rho + 3.0)
    ) * rho * pow(eta,5);

    teta4 = 4.0 * (
        64.0 * (2.0 - r) * pow(2.0 * rho + 1.0, 6) * pow(r,7)
        + 16.0 * (3072.0 * rho8 + 14208.0 * rho7 + 32640.0 * rho6
            + 45600.0 * rho5 + 40128.0 * rho4 + 19944.0 * rho3
            + 3936.0 * rho2 - 416.0 * rho + 47.0) * pow(r,6)
        - 32.0 * (2048.0 * rho9 + 13440.0 * rho8 + 40896.0 * rho7
            + 83232.0 * rho6 + 114352.0 * rho5 + 100080.0 * rho4
            + 47040.0 * rho3 + 6522.0 * rho2 - 2315.0 * rho + 109.0) * pow(r,5)
        - 4.0 * (30720.0 * rho10 + 171520.0 * rho9 + 489984.0 * rho8
            + 782336.0 * rho7 + 523264.0 * rho6 - 278016.0 * rho5
            - 750304.0 * rho4 - 448248.0 * rho3 - 34052.0 * rho2
            + 44872.0 * rho - 1249.0) * pow(r,4)
        + 8.0 * (49152.0 * rho11 + 362496.0 * rho10 + 1310208.0 * rho9
            + 3083520.0 * rho8 + 4941056.0 * rho7 + 5278912.0 * rho6
            + 3522528.0 * rho5 + 1257904.0 * rho4 + 124792.0 * rho3
            - 26752.0 * rho2 + 6858.0 * rho + 95.0) * pow(r,3)
        - (409600.0 * rho12 + 3520512.0 * rho11 + 14518272.0 * rho10
            + 39896576.0 * rho9 + 79227904.0 * rho8 + 112583424.0 * rho7
            + 109728640.0 * rho6 + 68306784.0 * rho5 + 23047216.0 * rho4
            + 1666408.0 * rho3 - 1144440.0 * rho2 - 166128.0 * rho
            + 8427.0) * pow(r,2)
        + 2.0 * (98304.0 * pow(rho,13) + 956416.0 * rho12 + 4396032.0 * rho11
            + 13556224.0 * rho10 + 31287040.0 * rho9 + 53866368.0 * rho8
            + 66479040.0 * rho7 + 55691456.0 * rho6 + 28814144.0 * rho5
            + 6900176.0 * rho4 - 715544.0 * rho3 - 686514.0 * rho2
            - 52929.0 * rho + 3149.0) * r
        - (18432.0 * pow(rho,13) + 172032.0 * rho12 + 751104.0 * rho11
            + 2294016.0 * rho10 + 5449472.0 * rho9 + 9633856.0 * rho8
            + 11960512.0 * rho7 + 9817680.0 * rho6 + 4782872.0 * rho5
            + 907196.0 * rho4 - 257442.0 * rho3 - 137050.0 * rho2
            - 4813.0 * rho + 451.0) * (2.0 * rho + 3.0)
    ) * pow(eta,4);

    teta3 = (256.0 * (2.0 - r) * (16.0 * rho3 + 24.0 * rho2 + 12.0 * rho + 1.0) * pow(r,7)
        - 128.0 * (192.0 * rho6 + 192.0 * rho5 - 480.0 * rho4 - 1596.0 * rho3
            - 1638.0 * rho2 - 676.0 * rho + 19.0) * pow(r,6)
        + 128.0 * (1344.0 * rho6 + 2976.0 * rho5 + 1008.0 * rho4
            - 5784.0 * rho3 - 7834.0 * rho2 - 3437.0 * rho + 198.0) * pow(r,5)
        + 64.0 * (2688.0 * rho8 + 6144.0 * rho7 - 8160.0 * rho6 - 43080.0 * rho5
            - 54644.0 * rho4 - 16706.0 * rho3 + 15085.0 * rho2 + 10263.0 * rho - 929.0) * pow(r,4)
        - 32.0 * (6144.0 * rho9 + 56832.0 * rho8 + 130560.0 * rho7 + 48384.0 * rho6
            - 267264.0 * rho5 - 486560.0 * rho4 - 345992.0 * rho3 - 92044.0 * rho2
            + 2650.0 * rho - 377.0) * pow(r,3)
        - 8.0 * (9216.0 * rho10 - 232448.0 * rho9 - 1256448.0 * rho8 - 2428224.0 * rho7
            - 1331808.0 * rho6 + 2314912.0 * rho5 + 4628128.0 * rho4 + 3289308.0 * rho3
            + 941094.0 * rho2 + 17232.0 * rho - 8719.0) * pow(r,2)
        + 8.0 * (24576.0 * rho11 - 64512.0 * rho10 - 970240.0 * rho9 - 3068160.0 * rho8
            - 4269952.0 * rho7 - 1520992.0 * rho6 + 3514576.0 * rho5 + 5456064.0 * rho4
            + 3202096.0 * rho3 + 682854.0 * rho2 - 49157.0 * rho - 5952.0) * r
        - (36864.0 * rho11 - 24576.0 * rho10 - 946176.0 * rho9 - 3121920.0 * rho8
            - 4634496.0 * rho7 - 2574464.0 * rho6 + 1765216.0 * rho5 + 3792160.0 * rho4
            + 2305160.0 * rho3 + 449576.0 * rho2 - 69262.0 * rho - 2607.0) * (2.0 * rho + 3.0)
    ) * pow(eta,3);

    teta2 = -2.0 * (128.0 * (r - 2.0) * pow(r,7)
        + 96.0 * (32.0 * rho3 + 32.0 * rho2 - 2.0 * rho - 35.0) * pow(r,6)
        - 64.0 * (432.0 * rho3 + 582.0 * rho2 + 171.0 * rho - 277.0) * pow(r,5)
        - 8.0 * (2688.0 * rho5 + 3552.0 * rho4 - 9624.0 * rho3 - 17988.0 * rho2
            - 7706.0 * rho + 3411.0) * pow(r,4)
        + 16.0 * (3072.0 * rho6 + 18816.0 * rho5 + 29472.0 * rho4 + 8568.0 * rho3
            - 16428.0 * rho2 - 11974.0 * rho - 71.0) * pow(r,3)
        + 2.0 * (4608.0 * rho7 - 259072.0 * rho6 - 904608.0 * rho5 - 1174896.0 * rho4
            - 475504.0 * rho3 + 231848.0 * rho2 + 222510.0 * rho + 11293.0) * pow(r,2)
        - 4.0 * (12288.0 * rho8 - 48384.0 * rho7 - 424480.0 * rho6 - 945360.0 * rho5
            - 923200.0 * rho4 - 270072.0 * rho3 + 173214.0 * rho2 + 121135.0 * rho
            - 51.0) * r
        + (4608.0 * rho8 - 5568.0 * rho7 - 145344.0 * rho6 - 368496.0 * rho5
            - 401040.0 * rho4 - 143908.0 * rho3 + 56548.0 * rho2 + 51463.0 * rho - 1467.0) * (2.0 * rho + 3.0)
    ) * pow(eta,2);

    teta1 = 2.0 * (64.0 * (29.0 - 3.0 * r) * pow(r,5)
        + 16.0 * (84.0 * rho2 - 325.0) * pow(r,4)
        - 96.0 * (80.0 * rho3 + 228.0 * rho2 + 120.0 * rho - 77.0) * pow(r,3)
        - 4.0 * (144.0 * rho4 - 19456.0 * rho3 - 33144.0 * rho2 - 16416.0 * rho + 4253.0) * pow(r,2)
        + 4.0 * (1920.0 * rho5 - 9456.0 * rho4 - 48992.0 * rho3 - 58296.0 * rho2 - 21104.0 * rho + 6229.0) * r
        - (288.0 * rho5 + 336.0 * rho4 - 15216.0 * rho3 - 22680.0 * rho2 - 11046.0 * rho + 2953.0) * (2.0 * rho + 3.0)
    ) * pow(eta,1);

    teta0 = -96.0 * (2.0 * r + 2.0 * rho - 15.0) * (2.0 * r - 2.0 * rho - 3.0) * (2.0 * r - 1.0);

    return teta6 + teta5 + teta4 + teta3 + teta2 + teta1 + teta0;
}

double c3(double r, double rho, double eta)
/*
 * Direct correlation function (part) c^(1)_3 in spherical geometry
 *
 *   inputs: r    ... distance from the sphere center
 *           rho  ... sphere radius
 *           eta  ... packing fraction
 *
 *   output: value of c_3(r) for given rho and eta
 */
{
    double den1, den2, logterm, rho2;
    double complex sum = 0.0 + 0.0 * I;
    double complex t[4];
    int k;
    rho2 = rho * rho;

    roots4a(eta, rho, t);

    den1 = 2.0 * (4.0 * rho2 + 6.0 * rho + 3.0) * eta * rho + 1.0;

    /* region 1: r < rho */
    if (r < rho)
    {
        return 0.0;
    }
    /* region 2: rho <= r < rho + 1/2 */
    else if ((r >= rho) && (r < rho + 0.5))
    {
        den2 = (4.0 * (r * r + r + 2.0 * r * rho - 3.0 * rho2 - 5.0 * rho) - 11.0) * eta * pow(2.0 * r - 2.0 * rho + 1.0, 2) + 16.0 * (2.0 * r + 1.0);
        logterm = log((2.0 * r + 1.0) / (2.0 * rho));
        logterm *= (2.0 * r + 1.0) * (2.0 * r - 1.0) * (rho2 + rho + 1.0);
        logterm /= (8.0 * r * rho * (rho + 1.0));

        sum = 0.0 + 0.0 * I;
        for (k = 0; k < 4; ++k)
        {
            double complex cden = (2.0 * (t[k] + 2.0 * rho) * (t[k] - rho) - 6.0 * rho - 3.0) * eta * (t[k] - rho) + 1.0;
            double complex lga = (2.0 * r - 2.0 * t[k] + 1.0) / (2.0 * rho - 2.0 * t[k]);
            sum += Q3(t[k], r, rho, eta) * clog(lga) / (32.0 * r * rho * (rho + 1.0) * cden * den1 * den1 * pow(eta - 1.0, 2));
        }

        double val = S3(r, rho, eta) / (32.0 * r * den2 * den1 * den1 * pow(eta - 1.0, 2));
        val = val - logterm + creal(sum);
        return val;
    }
    /* region 3: rho + 1/2 <= r < rho + 3/2 */
    else if ((r >= rho + 0.5) && (r < rho + 1.5))
    {
        den2 = (4.0 * (r * r - r + 2.0 * r * rho - 3.0 * rho2 - 7.0 * rho) - 11.0) * eta * pow(2.0 * r - 2.0 * rho - 1.0, 2) + 16.0 * (2.0 * r - 1.0);
        logterm = log((2.0 * r - 1.0) / (2.0 * rho + 2.0));
        logterm *= (2.0 * r + 1.0) * (2.0 * r - 1.0) * (rho2 + rho + 1.0);
        logterm /= (8.0 * r * rho * (rho + 1.0));

        sum = 0.0 + 0.0 * I;
        for (k = 0; k < 4; ++k)
        {
            double complex cden = (2.0 * (t[k] + 2.0 * rho) * (t[k] - rho) - 6.0 * rho - 3.0) * eta * (t[k] - rho) + 1.0;
            double complex lga = (2.0 * r - 2.0 * t[k] - 1.0) / (2.0 * rho - 2.0 * t[k] + 2.0);
            sum += Q3(t[k], r, rho, eta) * clog(lga) / (32.0 * r * rho * (rho + 1.0) * cden * den1 * den1 * pow(eta - 1.0, 2));
        }

        double val = T3(r, rho, eta) / (32.0 * r * den2 * den1 * den1 * pow(eta - 1.0, 3));
        val = val + logterm - creal(sum);
        return val;
    }
    /* region 4: r >= rho + 3/2 */
    else
    {
        return -eta * (1.0 + eta + eta * eta) / pow(1.0 - eta, 3);
    }

    return 0.0;
}
/********* end of c_3 calculation *************************/

double cv1(double r, double rho, double eta)
/*
 *  Direct correlation function (part) \vec{c}^{(1)}_1 in spherical geometry
 *     radial part of the vector function (normal to the wall)
 *
 *  inputs: r    ... distance from the sphere center
 *          rho  ... sphere radius
 *          eta  ... packing fraction
 *
 *  output: value of c_N1(r) for given rho and etab
 */
{
    double complex t[4];
    double r2 = r * r;
    double pref_common = 3.0 / (8.0 * r2);
    double nsum;
    double complex sum, denom, term, tk, t2, logarg;
    double rho2, rho3, rho4;

    rho2 = rho * rho;
    rho3 = rho2 * rho;
    rho4 = rho2 * rho2;

    roots4b(eta, rho, t);

    /* region 1: r < rho */
    if (r < rho)
    {
        return 0.0;
    }
    /* region 2: rho <= r < rho + 1/2 */
    else if ((r >= rho) && (r < rho + 0.5))
    {
        /* non-sum term */
        nsum = -((4.0 * r2 - 2.0 * r - 2.0 * r * rho - 26.0 * rho2 - 25.0 * rho - 11.0) * (2.0 * r - 2.0 * rho + 1.0)) / (8.0 * r2);

        /* sum over roots of P_4^b */
        sum = 0.0 + 0.0 * I;
        for (int k = 0; k < 4; ++k)
        {
            tk = t[k];
            t2 = tk * tk;

            /* denominator (2 t^2 + 6(t-1) rho - 3) * eta * t + 1 */
            denom = (2.0 * t2 + 6.0 * (tk - 1.0) * rho - 3.0) * eta * tk + 1.0;

            /* numerator polynomial pieces (real coefficients multiplied by complex t powers) */
            term = 4.0 * (4.0 * rho2 + 6.0 * rho + 3.0) * rho * eta + 4.0; /* multiplies t^3 */

            term = term * tk + ((8.0 * rho2 + 8.0 * rho + 4.0) * eta * r2 - (8.0 * rho4 + 40.0 * rho3 + 66.0 * rho2 + 46.0 * rho + 11.0) * eta + 12.0 * rho); /* multiplies t^2 */

            term = term * tk + ((2.0 * rho + 1.0) * (2.0 * rho + 1.0) * (2.0 * rho - 1.0) * rho * eta - 4.0 * (2.0 * rho + 1.0) * rho * eta * r2 - 4.0 * r2 + 28.0 * rho2 + 16.0 * rho + 7.0); /* multiplies t */

            term = term * tk - (4.0 * r2 - 20.0 * rho2 - 16.0 * rho - 7.0) * rho; /* constant */

            /* log argument: (2r - 2t - 2rho + 1)/(-2t) */
            logarg = (2.0 * r - 2.0 * tk - 2.0 * rho + 1.0) / (-2.0 * tk);

            sum += term / denom * clog(logarg);
        }

        /* result: nsum  - (3/(8 r^2)) * sum  (take real part) */
        return nsum - pref_common * creal(sum);
    }
    /* region 3: rho + 1/2 <= r < rho + 3/2 */
    else if ((r >= rho + 0.5) && (r < rho + 1.5))
    {
        nsum = ((4.0 * r2 - 2.0 * r * rho - 26.0 * rho2 - 27.0 * rho - 12.0) * (2.0 * r - 2.0 * rho - 3.0)) / (8.0 * r2);

        sum = 0.0 + 0.0 * I;
        for (int k = 0; k < 4; ++k)
        {
            tk = t[k];
            t2 = tk * tk;

            /* denominator = (2 t^2 + 6(t-1) rho - 3) * eta * t + 1 */
            denom = (2.0 * t2 + 6.0 * (tk - 1.0) * rho - 3.0) * eta * tk + 1.0;

            /* polynomial numerator (same coefficients as before) */
            term = 4.0 * (4.0 * rho2 + 6.0 * rho + 3.0) * rho * eta + 4.0; /* t^3 */

            term = term * tk + ((8.0 * rho2 + 8.0 * rho + 4.0) * eta * r2 - (8.0 * rho4 + 40.0 * rho3 + 66.0 * rho2 + 46.0 * rho + 11.0) * eta + 12.0 * rho); /* t^2 */

            term = term * tk + ((2.0 * rho + 1.0) * (2.0 * rho + 1.0) * (2.0 * rho - 1.0) * rho * eta - 4.0 * (2.0 * rho + 1.0) * rho * eta * r2 - 4.0 * r2 + 28.0 * rho2 + 16.0 * rho + 7.0); /* t */

            term = term * tk - (4.0 * r2 - 20.0 * rho2 - 16.0 * rho - 7.0) * rho; /* constant */

            /* log argument: (2r - 2t - 2rho - 1)/(2 - 2t) */
            logarg = (2.0 * r - 2.0 * tk - 2.0 * rho - 1.0) / (2.0 - 2.0 * tk);

            sum += term / denom * clog(logarg);
        }

        /* result: nsum  + (3/(8 r^2)) * sum (take real part) */
        return nsum + pref_common * creal(sum);
    }
    /* otherwise (r > rho + 3/2) */
    else
    {
        return 0.0;
    }
    /* otherwise */
    return 0.0;
}

double cv2(double r, double rho, double eta)
/*
 * Direct correlation function (part) \vec{c}^{(1)}_2 in spherical geometry
 *    radial part of the vector function (normal to the wall)
 *
 * inputs: r    ... distance from the sphere center
 *         rho  ... sphere radius
 *         eta  ... packing fraction
 *
 * output: value of c_N2(r) for given rho and etab
 */
{
    double rr = r * r;
    double rh2 = rho * rho;
    double rh3 = rh2 * rho;
    double rh4 = rh2 * rh2;

    double term1, term2, term3;
    double complex sum, tk, den_t, den_t4;
    double complex termt3, termt2, termt1, termt0;
    double complex logarg;
    double complex t[4];
    int k;

    roots4b(eta, rho, t);

    /* region 1: r < rho */
    if (r < rho)
    {
        return 0.0;
    }
    /* region 2: rho <= r < rho + 1/2 */
    else if (r >= rho && r < rho + 0.5)
    {

        /* Term 1 */
        term1 =
            -3.0 * (2 * r + 2 * rho + 3.0) * (2 * r - 2 * rho + 1.0) * (2 * r + 1.0) * eta / (((4 * rr + 4 * r + 8 * r * rho - 12 * rh2 - 20 * rho - 11) * eta * pow(2 * r - 2 * rho + 1.0, 2.0) + 16 * (2 * r + 1.0)) * r);

        /* Term 2 */
        term2 =
            -(4 * rr - 2 * r * (rho + 1.0) - 26 * rh2 - 25 * rho - 35.0) * (2 * r - 2 * rho + 1.0) / (8.0 * rr);

        /* Summation */
        sum = 0.0 + 0.0 * I;

        for (k = 0; k < 4; k++)
        {
            tk = t[k];

            /* Compute denominators */
            den_t = ((2.0 * tk * tk + 6.0 * (tk - 1.0) * rho - 3.0) * eta * tk) + 1.0;
            den_t4 = 4.0 * den_t;

            /* term t^3 */
            termt3 = ((8 * rh3 + 12 * rh2 + 6 * rho - 3.0) * eta + 2.0) / den_t * tk * tk * tk;

            /* term t^2 */
            termt2 = (4 * (rh2 + rho + 1.0) * eta * rr - (4.0 * rh4 + 20.0 * rh3 + 55.0 * rh2 + 54.0 * rho + 17.0) * eta + 6.0 * rho) / den_t * tk * tk;

            /* term t^1 */
            termt1 =
                -(4 * (4 * rh2 - 2 * rho - 1.0) * eta * rr - (8 * rh3 - 24 * rh2 + 1.0) * (2 * rho + 1.0) * eta + 8 * rr - 56 * rh2 - 32 * rho - 46.0) / den_t4 * tk;

            /* term t^0 */
            termt0 =
                (((4.0 * rr - 4.0 * rh2 + 1.0) * eta * (2.0 * rho + 1.0) - 8.0 * rr + 40.0 * rh2 + 32.0 * rho + 46.0) * rho) / den_t4;

            /* Log argument */
            logarg = (2 * r - 2 * tk - 2 * rho + 1.0) / (-2.0 * tk);

            sum += (termt3 + termt2 + termt1 + termt0) * clog(logarg);
        }

        return term1 + term2 - (3.0 * creal(sum)) / (4.0 * rr);
    }
    /* region 3: rho + 1/2 <= r < rho + 3/2 */
    else if (r >= rho + 0.5 && r < rho + 1.5)
    {

        /* Term 1 */
        term1 =
            (4 * rr - 2 * r * (rho - 1.0) - 26 * rh2 - 23 * rho - 29.0) * (2 * r - 2 * rho - 1.0) / (8.0 * rr);

        /* Term 2 */
        term2 =
            -((12 * rr * (rho + 2.0) - 12 * rh3 - 96 * rh2 - 105 * rho - 94.0) * eta - 12 * rr + 60 * rh2 + 72 * rho + 85.0) / (8.0 * (eta - 1.0) * rr);

        /* Term 3 */
        term3 =
            -3.0 * ((16 * pow(r, 4) + 16 * r * r * r * rho - 8 * (10 * rh2 + 11 * rho + 6) * rr + 12 * (4 * rh2 + 8 * rho + 3) * rho * r + (6 * rho + 11) * pow(2 * rho + 1.0, 2.0)) * eta * (2 * r - 2 * rho - 1.0) + 16 * (2 * r + 1.0) * (2 * r - 1.0)) / (4.0 * ((4 * rr - 4 * r + 8 * r * rho - 12 * rh2 - 28 * rho - 11) * eta * pow(2 * r - 2 * rho - 1.0, 2.0) + 16 * (2 * r - 1.0)) * rr);

        /* Summation */
        sum = 0.0 + 0.0 * I;

        for (k = 0; k < 4; k++)
        {
            tk = t[k];

            den_t = ((2.0 * tk * tk + 6.0 * (tk - 1.0) * rho - 3.0) * eta * tk) + 1.0;
            den_t4 = 4.0 * den_t;

            termt3 = ((8 * rh3 + 12 * rh2 + 6 * rho - 3.0) * eta + 2.0) / den_t * tk * tk * tk;

            termt2 = (4 * (rh2 + rho + 1.0) * eta * rr - (4 * rh4 + 20 * rh3 + 55 * rh2 + 54 * rho + 17) * eta + 6 * rho) / den_t * tk * tk;

            termt1 = -(4 * (4 * rh2 - 2 * rho - 1.0) * eta * rr - (8 * rh3 - 24 * rh2 + 1.0) * (2 * rho + 1.0) * eta + 8 * rr - 56 * rh2 - 32 * rho - 46.0) / den_t4 * tk;

            termt0 = (((4 * rr - 4 * rh2 + 1.0) * eta * (2 * rho + 1.0) - 8 * rr + 40 * rh2 + 32 * rho + 46.0) * rho) / den_t4;

            /* Log argument (region 2) */
            logarg = (2 * r - 2 * tk - 2 * rho - 1.0) / (2.0 - 2.0 * tk);

            sum += (termt3 + termt2 + termt1 + termt0) * clog(logarg);
        }

        return term1 + term2 + term3 + (3.0 * creal(sum)) / (4.0 * rr);
    }

    /* Outside both regions */
    return 0.0;
}

double cHS(double r, double rho, double eta)
/*
 * Direct correlation function c^(1)_HS
 *     Hard-sphere contribution sum
 *
 * Parameters:
 *     r       ... distance from the sphere center
 *     rho     ... radius of the sphere (spherical wall)
 *     etab    ... bulk density fraction
 *
 * Output:
 *      c_HS    ... sum of c_3, c_2, c_1, c_0, c_N1 and c_N2
 */
{
    return c3(r, rho, eta) + c2(r, rho, eta) + c1(r, rho, eta) + c0(r, rho, eta) + cv1(r, rho, eta) + cv2(r, rho, eta);
}

double catt(double r, double rho, double eta, double beta, double epsilon, double rcut)
/*
 * correlation function c^(1)_att (piecewise defined)
 *     attractive contribution to direct correlation function due to cut lj potential
 *     spherical geometry
 *
 * parameters:
 *     r       ... distance from the sphere center
 *     rho     ... radius of the sphere (spherical wall)
 *     eta     ... bulk packing fraction
 *     beta    ... inverse temperature 1/kt
 *     epsilon ... depth of the lj potential well
 *     rcut    ... cut-off distance of the lj potential
 *
 * output:
 *     c^(1)_att   ... contribution of attractive interparticle potential to direct correlation function
 */
{
    double cat;
    double rc2 = rcut * rcut;
    double rc3 = rc2 * rcut;
    double rc4 = rc2 * rc2;

    if (r >= rho + 0.5 + rcut)
    {
        cat = 32.0 * (rc3 - 1.0) / (rc3);
    }
    else if (r >= rho + 1.5)
    {
        cat = 6.0 * (8.0 * (2.0 * rc3 - 1.0) / (3.0 * rc3) - r / (rc4) + ((2.0 * rho + 1.0) * (2.0 * rho + 1.0) - 8.0 * rc2) / (4.0 * r * rc4) + 4.0 * (6.0 * rho - 2.0 * r + 3.0) / (3.0 * r * pow(2.0 * rho - 2.0 * r + 1.0, 3.0)));
    }
    else if (r >= rho)
    {
        cat = 6.0 * (8.0 * (rc3 - 1.0) / (3.0 * rc3) + r * (rc4 - 1.0) / (rc4) + ((2.0 * rho + 1.0) * (2.0 * rho + 1.0) - 8.0 * rc2) / (4.0 * r * rc4) - (4.0 * rho * rho + 4.0 * rho - 7.0) / (4.0 * r));
    }
    else
    {
        cat = 0.0;
    }

    return eta * beta * epsilon * cat;
}

double BetaMuHS(double etab)
/*
 * Bulk (reduced) chemical potential as a function of density (hard-sphere fluid)
 *
 * Parameters:
 *   etab - bulk packing fraction
 *
 * Returns:
 *   BetaMuHS - reduced chemical potential (beta * mu)
 */
{
    return (5.0 * pow(etab, 3.0) - 13.0 * pow(etab, 2.0) + 14.0 * etab) / 2.0 / pow(1.0 - etab, 3.0) - log(1.0 - etab) + log(etab * 6.0 / M_PI);
}

double BetaMuLJ(double etab, double beta, double epsilon, double rcut)
/*
 * Bulk (reduced) chemical potential including mean-field LJ attraction
 *
 * Parameters:
 *   etab    - bulk packing fraction
 *   beta    - inverse temperature (1/kT)
 *   epsilon - depth of the LJ potential well
 *   rcut    - cut-off distance of the LJ potential
 *
 * Returns:
 *   BetaMuLJ - reduced chemical potential including LJ attraction
 */
{
    return BetaMuHS(etab) + 32.0 * etab * beta * epsilon * (1.0 / (rcut * rcut * rcut) - 1.0);
}

double BetaVHS(double r, double rho, double eta)
/*
 * External potential (reduced) inducing constant density profile
 * for hard-sphere fluid (spherical geometry)
 *
 * Parameters:
 *   r    - distance from the center of spherical wall
 *   rho  - radius of the spherical wall
 *   eta  - bulk packing fraction
 *
 * Returns:
 *   BetaVHS - reduced external potential beta * V_ext
 */
{
    if (r < rho + 0.5)
    {
        return M_INF; /* effectively infinite potential */
    }
    return BetaMuHS(eta) + cHS(r, rho, eta) - log(eta * 6.0 / M_PI);
}

double BetaVLJ(double r, double rho, double eta, double beta, double epsilon, double rcut)
/*
 * External potential (reduced) inducing constant density profile
 * for hard-sphere fluid including LJ attractions
 *   - spherical geometry
 *
 * Parameters:
 *   r       - distance from the center
 *   rho     - radius of the spherical wall
 *   eta     - bulk packing fraction
 *   beta    - inverse temperature (1/kT)
 *   epsilon - depth of the LJ potential well
 *   rcut    - cut-off distance of the LJ potential
 *
 * Returns:
 *   BetaVLJ - reduced external potential beta * V_ext
 */
{
    if (r < rho + 0.5)
    {
        return M_INF; /* effectively infinite potential */
    }
    return BetaMuLJ(eta, beta, epsilon, rcut) + cHS(r, rho, eta) - log(eta * 6.0 / M_PI) + catt(r, rho, eta, beta, epsilon, rcut);
}

double BoltzFHS(double r, double rho, double eta)
/*
 * Boltzmann factor exp(-beta V_ext) inducing constant density profile
 * for hard-sphere fluid in spherical geometry
 *
 * Parameters:
 *   r    - distance from the center of spherical wall
 *   rho  - radius of the spherical wall
 *   eta  - bulk packing fraction
 *
 * Returns:
 *   BoltzFHS - Boltzmann factor exp(-beta V_ext)
 */
{
    if (r < rho + 0.5)
    {
        return 0.0; /* effectively infinite potential */
    }
    double val = BetaMuHS(eta) + cHS(r, rho, eta) - log(eta * 6.0 / M_PI);
    return exp(-val);
}

double BoltzFLJ(double r, double rho, double eta, double beta, double epsilon, double rcut)
/*
 * Boltzmann factor exp(-beta V_ext) inducing constant density profile
 * hard-sphere fluid with LJ attractions
 *   - spherical geometry
 *
 * Parameters:
 *   r       - distance from the center
 *   rho     - radius of the spherical wall
 *   eta     - bulk packing fraction
 *   beta    - inverse temperature (1/kT)
 *   epsilon - depth of the LJ potential well
 *   rcut    - cut-off distance of the LJ potential
 *
 * Returns:
 *   BoltzFLJ - Boltzmann factor exp(-beta V_ext)
 */
{
    if (r < rho + 0.5)
    {
        return 0.0; /* effectively infinite potential */
    }
    double val = BetaMuLJ(eta, beta, epsilon, rcut) + cHS(r, rho, eta) - log(eta * 6.0 / M_PI) + catt(r, rho, eta, beta, epsilon, rcut);
    return exp(-val);
}
