#include <math.h>
#include <complex.h>

/*
 * wall_pl
 * External potential inducing constant density profile
 * - Planar geometry
 * - External potential and direct correlation functions
 * - Hard-sphere fluid or hard-sphere fluid with cut LJ attractions (mean-field)
 * - Classical DFT using Rosenfeld's Fundamental Measure Theory (FMT)
 * - For HS contribution, both exact analytical expr. and Padé approximation
 *
 * Author: Jiří Janek
 * Date: January 2026
 *
 */

#define M_INF 1.0e6
#ifndef M_PI
#define M_PI 3.141592653589793238462643383279
#endif

void roots(double etab, double complex *t)
{
    /*
     * Roots of the cubic equation (needed for summation):
     *     2*eta*t^3 - 3*eta*t^2 + 1 = 0
     *
     * Parameters:
     *   etab - bulk packing fraction (input)
     *   t    - output array (length 3) where complex roots will be stored
     */
    double complex rootalpha, rootbeta, rootgamma;
    double complex cetab = etab + 0.0 * I;

    rootalpha = cpow((cetab - 2.0 + 2.0 * csqrt(1.0 - cetab)) * cetab * cetab, 1.0 / 3.0);
    rootbeta = -0.25 * (rootalpha / cetab + cetab / rootalpha) + 0.5;
    rootgamma = 0.25 * sqrt(3.0) * (rootalpha / cetab - cetab / rootalpha);
    t[0] = 0.5 * (rootalpha / cetab + cetab / rootalpha + 1.0);
    t[1] = rootbeta + rootgamma * I;
    t[2] = rootbeta - rootgamma * I;
}

/*
 * Correlation function c^(1)_HS (piecewise defined)
 * - Returns the hard-sphere contribution to the direct correlation function
 * - Piecewise regions:
 *     z > 1.5       : asymptotic bulk value (independent of z)
 *     0.5 <= z <=1.5: analytic form requiring summation over cubic roots
 *     0 <= z < 0.5  : analytic form requiring summation over cubic roots
 *
 * Parameters:
 *   z    - distance from the wall
 *   etab - bulk packing fraction
 *
 * Returns:
 *   cHS   - value of c^(1)_HS at distance z
 */
double cHS(double z, double etab)
{
    double complex t[3];
    double complex cpartial, zc, cetab;
    double crepart, fact;
    int j;

    roots(etab, t);
    zc = z + 0.0 * I;
    cetab = etab + 0.0 * I;

    if (z > 1.5)
    {
        return log(1.0 - etab) - (etab * (14.0 - 13.0 * etab + 5.0 * etab * etab)) / (2.0 * (1.0 - etab) * (1.0 - etab) * (1.0 - etab));
    }
    else if (z >= 0.5)
    {
        cpartial = 0.0 + 0.0 * I;
        for (j = 0; j < 3; j++)
        {
            cpartial += (clog(1.0 - t[j]) - clog(2.0 * zc - 2.0 * t[j] - 1.0)) * (2.0 * z - 1.0) / (96.0 * t[j] * (t[j] - 1.0) * (cetab - 1.0) * (cetab - 1.0)) * ((48.0 * t[j] * t[j] - 16.0 * t[j] - 16.0 * zc - 8.0) * cetab * cetab - (96.0 * t[j] * t[j] + 36.0 * zc * t[j] - 70.0 * t[j] - 22.0 * zc - 23.0) * cetab + 48.0 * t[j] * t[j] - 24.0 * zc + 12.0);
        }
        fact = pow(etab - 1.0, 3.0) * (2.0 + (2.0 * z - 1.0) * (2.0 * z - 1.0) * (z - 2.0) * etab);
        crepart = (64.0 * pow(z, 6.0) - 320.0 * pow(z, 5.0) + 592.0 * pow(z, 4.0) - 544.0 * pow(z, 3.0) + 324.0 * pow(z, 2.0) - 108.0 * z + 13.0) * pow(etab, 4.0) / 8.0 / fact;
        crepart += (128.0 * pow(z, 6.0) - 848.0 * pow(z, 5.0) + 2176.0 * pow(z, 4.0) - 2488.0 * pow(z, 3.0) + 1032.0 * pow(z, 2.0) - 141.0 * z - 36.0) * pow(etab, 3.0) / 16.0 / fact;
        crepart += (128.0 * pow(z, 6.0) - 816.0 * pow(z, 5.0) + 1824.0 * pow(z, 4.0) - 1864.0 * pow(z, 3.0) + 1032.0 * pow(z, 2.0) - 123.0 * z + 90.0) * etab * etab / 16.0 / fact;
        crepart += (4.0 * z * z * z - 12.0 * z * z + 3.0 * z - 5.0) * etab / fact - (z + 0.5) * (log(1.0 - etab));
        crepart += (z - 0.5) * (log(2.0 + (2.0 * z - 1.0) * (2.0 * z - 1.0) * (z - 2.0) * etab) + 2.0 * log(2.0));
        return -creal(cpartial) - crepart;
    }
    else if (z >= 0.0)
    {
        cpartial = 0.0 + 0.0 * I;
        for (j = 0; j < 3; j++)
        {
            cpartial += (clog(-t[j]) - clog(2.0 * zc - 2.0 * t[j] + 1.0)) * (2.0 * z - 1.0) / (96.0 * t[j] * (t[j] - 1.0) * (cetab - 1.0) * (cetab - 1.0)) * ((-48.0 * t[j] * t[j] + 16.0 * t[j] + 16.0 * zc + 8.0) * cetab * cetab + (96.0 * t[j] * t[j] + 36.0 * zc * t[j] - 70.0 * t[j] - 22.0 * zc - 23.0) * cetab - 48.0 * t[j] * t[j] + 24.0 * zc - 12.0);
        }
        fact = 2.0 * etab * (2.0 * z + 1.0) / (pow(etab - 1.0, 2.0) * (2.0 + (2.0 * z + 1.0) * (2.0 * z + 1.0) * (z - 1.0) * etab));
        crepart = (pow(z, 4.0) - 0.5 * pow(z, 3.0) - 0.75 * z * z + 0.125 * z + 0.125) * etab * etab * etab;
        crepart -= (2.5 * pow(z, 4.0) - 1.25 * pow(z, 3.0) - 1.125 * z * z - 0.4375 * z - 0.25) * etab * etab;
        crepart -= (0.75 * pow(z, 4.0) - 2.625 * pow(z, 3.0) - 1.3125 * z * z + 3.65625 * z + 1.03125) * etab - 1.5;
        crepart *= fact;
        crepart -= (z + 0.5) * log(2.0 + (2.0 * z + 1.0) * (2.0 * z + 1.0) * (z - 1.0) * etab) + 2.0 * (z - 1.0) * log(2.0);
        return -creal(cpartial) - crepart;
    }
    else
    {
        return 0.0;
    }

    return 0.0;
}

/*
 * Correlation function c^(1)_HS (Padé approx.)
 * - Returns the hard-sphere contribution to the direct correlation function
 * - Piecewise regions:
 *     z > 1.5       : asymptotic bulk value (independent of z)
 *     0.5 <= z <=1.5: Padé approx.
 *     0 <= z < 0.5  : Padé approx. (imprecise, not needed for Vext)
 *
 * Parameters:
 *   z    - distance from the wall
 *   etab - bulk packing fraction
 *
 * Returns:
 *   cHSPade   - value of c^(1)_HS at distance z
 */
double cHSPade(double z, double etab)
{
    double A0, A1, A2, B1, B2;
    double aux1, aux2, aux3, ln1me;
    double den1, den2, den3;
    double e2, e3, e4, e5, e6;
    double zarg;

    /* Powers of etab */
    e2 = etab * etab;
    e3 = e2 * etab;
    e4 = e2 * e2;
    e5 = e2 * e3;
    e6 = e3 * e3;

    /* Denominators */
    den1 = 2.0 * pow(1.0 - etab, 3.0); // den. for A0, A1 and A2
    den2 = 3.0 * (1.0 - etab) * (etab + 2.0); // den. fro A1, A2 and B1
    den3 = e2 - 2.0 * etab + 4.0; // den. for A2 and B2

    /* Auxiliary quantities */
    aux1 = etab * (5.0 * e2 - 13.0 * etab + 14.0) / den1; // first term of A0, second term in A1

    aux2 = -25.0 * e4 + 4.0 * e3 - 9.0 * e2 - 20.0 * etab - 4.0; // first term in A2, numerator in B2

    aux3 = 107.0 * e6 - 345.0 * e5 + 501.0 * e4 - 253.0 * e3 - 114.0 * e2 + 660.0 * etab - 232.0; // second term in A2

    ln1me = log(1.0 - etab);

    /* Coefficients */
    B1 = (13.0 * e2 + 7.0 * etab - 2.0) / den2;
    A0 = aux1 - ln1me;
    A1 = -B1 * ln1me + aux1 * B1;
    A2 = (aux2 * ln1me + aux3 * etab / den1) / den2 / den3;
    B2 = -aux2 / den2 / den3;

    zarg = z - 1.5;

    /* Piecewise definition */
    if (z > 1.5)
    {
        return log(1.0 - etab) - (etab * (14.0 - 13.0 * etab + 5.0 * etab * etab)) / (2.0 * (1.0 - etab) * (1.0 - etab) * (1.0 - etab)); // bulk limit
    }
    else if (z > 0.5)
    {
        return -(A0 + A1 * zarg + A2 * zarg * zarg) /
               (1.0 + B1 * zarg + B2 * zarg * zarg);
    }
    else if (z >= 0.0)
    {
        return -(A0 + A1 * zarg + A2 * zarg * zarg) /
               (1.0 + B1 * zarg + B2 * zarg * zarg); // irrelevant for Vext
    }
    else
    {
        return 0.0;
    }
}

/*
 * Mean-field Lennard-Jones attractive contribution to the direct correlation function in planar geometry
 *
 * Parameters:
 *   z       - distance from the wall
 *   etab    - bulk packing fraction
 *   beta    - inverse temperature (1/kT)
 *   epsilon - depth of the LJ potential well
 *   rcut    - cut-off distance of the LJ potential
 *
 * Returns:
 *   c^(1)_att(z) - attractive contribution to the direct correlation function
 */
double catt(double z, double etab, double beta, double epsilon, double rcut)
{
    if (z > (0.5 + rcut))
    {
        return -32.0 * etab * beta * epsilon * (1.0 / (rcut * rcut * rcut) - 1.0);
    }
    else if (z > 1.5)
    {
        double term = 1.0 / (6.0 * rcut * rcut * rcut) - 1.0 / 3.0 + (2.0 * z - 1.0) / (16.0 * rcut * rcut * rcut * rcut) - 1.0 / (3.0 * pow(1.0 - 2.0 * z, 3.0));
        return -96.0 * etab * beta * epsilon * term;
    }
    else if (z >= 0.0)
    {
        double term = 1.0 / (3.0 * rcut * rcut * rcut) - 5.0 / 24.0 + (2.0 * z - 1.0) / (8.0 * rcut * rcut * rcut * rcut) - z / 4.0;
        return -48.0 * etab * beta * epsilon * term;
    }
    else
    {
        return 0.0;
    }
    return 0.0;
}

/*
 * Bulk (reduced) chemical potential as a function of density (hard-sphere fluid)
 *
 * Parameters:
 *   etab - bulk packing fraction
 *
 * Returns:
 *   BetaMuHS - reduced chemical potential (beta * mu)
 */
double BetaMuHS(double etab)
{
    return (5.0 * pow(etab, 3.0) - 13.0 * pow(etab, 2.0) + 14.0 * etab) / 2.0 / pow(1.0 - etab, 3.0) - log(1.0 - etab) + log(etab * 6.0 / M_PI);
}

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
double BetaMuLJ(double etab, double beta, double epsilon, double rcut)
{
    return BetaMuHS(etab) + 32.0 * etab * beta * epsilon * (1.0 / (rcut * rcut * rcut) - 1.0);
}

/*
 * External potential (reduced) inducing constant density profile
 * for hard-sphere fluid (planar geometry)
 *
 * Parameters:
 *   z    - distance from the wall
 *   etab - bulk packing fraction
 *
 * Returns:
 *   BetaVHS - reduced external potential beta * V_ext
 */
double BetaVHS(double z, double etab)
{
    if (z < 0.5)
    {
        return M_INF; /* effectively infinite potential */
    }
    return BetaMuHS(etab) + cHS(z, etab) - log(etab * 6.0 / M_PI);
}

/*
 * External potential (reduced) inducing constant density profile
 * for hard-sphere fluid (planar geometry)
 *  - Padé approximation for c^(1)_HS
 *
 * Parameters:
 *   z    - distance from the wall
 *   etab - bulk packing fraction
 *
 * Returns:
 *   BetaVHSPade - reduced external potential beta * V_ext
 */
double BetaVHSPade(double z, double etab)
{
    if (z < 0.5)
    {
        return M_INF; /* effectively infinite potential */
    }
    return BetaMuHS(etab) + cHSPade(z, etab) - log(etab * 6.0 / M_PI);
}

/*
 * External potential (reduced) inducing constant density profile
 * for hard-sphere fluid including LJ attractions
 * planar geometry
 *
 * Parameters:
 *   z       - distance from the wall
 *   etab    - bulk packing fraction
 *   beta    - inverse temperature (1/kT)
 *   epsilon - depth of the LJ potential well
 *   rcut    - cut-off distance of the LJ potential
 *
 * Returns:
 *   BetaVLJ - reduced external potential beta * V_ext
 */
double BetaVLJ(double z, double etab, double beta, double epsilon, double rcut)
{
    if (z < 0.5)
    {
        return M_INF; /* effectively infinite potential */
    }
    return BetaMuLJ(etab, beta, epsilon, rcut) + cHS(z, etab) - log(etab * 6.0 / M_PI) + catt(z, etab, beta, epsilon, rcut);
}

/*
 * External potential (reduced) inducing constant density profile
 * for hard-sphere fluid including LJ attractions
 * planar geometry
 *  - Padé approximation for c^(1)_HS
 *
 * Parameters:
 *   z       - distance from the wall
 *   etab    - bulk packing fraction
 *   beta    - inverse temperature (1/kT)
 *   epsilon - depth of the LJ potential well
 *   rcut    - cut-off distance of the LJ potential
 *
 * Returns:
 *   BetaVLJHS - reduced external potential beta * V_ext
 */
double BetaVLJPade(double z, double etab, double beta, double epsilon, double rcut)
{
    if (z < 0.5)
    {
        return M_INF; /* effectively infinite potential */
    }
    return BetaMuLJ(etab, beta, epsilon, rcut) + cHSPade(z, etab) - log(etab * 6.0 / M_PI) + catt(z, etab, beta, epsilon, rcut);
}

/*
 * Boltzmann factor exp(-beta V_ext) inducing constant density profile
 * for hard-sphere fluid in planar geometry
 *
 * Parameters:
 *   z    - distance from the wall
 *   etab - bulk packing fraction
 *
 * Returns:
 *   BoltzFHS - Boltzmann factor exp(-beta V_ext)
 */
double BoltzFHS(double z, double etab)
{
    if (z < 0.5)
    {
        return 0.0; /* effectively infinite potential */
    }
    double val = BetaMuHS(etab) + cHS(z, etab) - log(etab * 6.0 / M_PI);
    return exp(-val);
}

/*
 * Boltzmann factor exp(-beta V_ext) inducing constant density profile
 * for hard-sphere fluid in planar geometry
 *  - Padé approximation for c^(1)_HS
 *
 * Parameters:
 *   z    - distance from the wall
 *   etab - bulk packing fraction
 *
 * Returns:
 *   BoltzFHSPade - Boltzmann factor exp(-beta V_ext)
 */
double BoltzFHSPade(double z, double etab)
{
    if (z < 0.5)
    {
        return 0.0; /* effectively infinite potential */
    }
    double val = BetaMuHS(etab) + cHSPade(z, etab) - log(etab * 6.0 / M_PI);
    return exp(-val);
}

/*
 * Boltzmann factor exp(-beta V_ext) inducing constant density profile
 * hard-sphere fluid with LJ attractions
 * planar geometry
 *
 * Parameters:
 *   z       - distance from the wall
 *   etab    - bulk packing fraction
 *   beta    - inverse temperature (1/kT)
 *   epsilon - depth of the LJ potential well
 *   rcut    - cut-off distance of the LJ potential
 *
 * Returns:
 *   BoltzFLJ - Boltzmann factor exp(-beta V_ext)
 */
double BoltzFLJ(double z, double etab, double beta, double epsilon, double rcut)
{
    if (z < 0.5)
    {
        return 0.0; /* effectively infinite potential */
    }
    double val = BetaMuLJ(etab, beta, epsilon, rcut) + cHS(z, etab) - log(etab * 6.0 / M_PI) + catt(z, etab, beta, epsilon, rcut);
    return exp(-val);
}

/*
 * Boltzmann factor exp(-beta V_ext) inducing constant density profile
 * hard-sphere fluid with LJ attractions
 * planar geometry
 *  - Padé approximation for c^(1)_HS
 *
 * Parameters:
 *   z       - distance from the wall
 *   etab    - bulk packing fraction
 *   beta    - inverse temperature (1/kT)
 *   epsilon - depth of the LJ potential well
 *   rcut    - cut-off distance of the LJ potential
 *
 * Returns:
 *   BoltzFLJPade - Boltzmann factor exp(-beta V_ext)
 */
double BoltzFLJPade(double z, double etab, double beta, double epsilon, double rcut)
{
    if (z < 0.5)
    {
        return 0.0; /* effectively infinite potential */
    }
    double val = BetaMuLJ(etab, beta, epsilon, rcut) + cHSPade(z, etab) - log(etab * 6.0 / M_PI) + catt(z, etab, beta, epsilon, rcut);
    return exp(-val);
}
