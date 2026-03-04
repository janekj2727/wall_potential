# Wall potential generating constant density profile in classical DFT

This repository provides C, Fortran and Python implementations of the formulae from the article [[1]](#1).
The aim is to provide exact analytical expressions for the potential of the wall that imposes a constant density profile in classical DFT (Density Functional Theory) calculations.
Hard-sphere (HS) fluid and a fluid with HS repulsion and cut LJ attraction are considered.
Besides the planar wall (programs in folder `\planar`), the expressions for the spherical geometry are implemented in `\spherical`.
The original Rosenfeld's FMT (Fundamental Measure Theory) functional [[2]](#2) is assumed.

Provided functions are documented in the code.

## Usage

### C

The code in C is provided as a header file (`wall_pl.h` or `wall_sph.h`) that can be easily included in your code.
The standard `math.h` library is used together with `complex.h` (C99 standard) to deal with complex roots of cubic and quartic polynomials that are needed for calculation of direct correlation function.

In the spherical geometry, roots of quartic polynomials are needed for which we do not provide analytical expressions.
Library `gsl_poly.h` from the GNU Scientific Library [[3]](#3) is used for numerical root finding and must be installed before compiling your code.
The compilation then looks like:
```(bash)
gcc -std=c99 -O2 YOUR_MAIN_PROGRAM.c -lm -lgsl -lgslcblas -o EXECUTABLE_NAME
```

### Fortran

The code in Fortran is provided as a standalone module (defined in `wall_pl.for` or `wall_sph.for` for the planar and spherical geometry, respectively).
Similarly to the C code, module `polyroots_module` from the `polyroots-fortran` GitHub repository by J. Williams [[4]](#4) together with `LAPACK` and `BLAS` libraries is needed for the calculation of the quartic polynomial roots.
The compilation of your main program (`YOUR_PROGRAM_NAME.for`) thus will look like:
```(bash)
gfortran -Wall -O2 -c polyroots_module.f90
gfortran -Wall -Wno-tabs -c -O2 wall_sph.for
gfortran -Wall -Wno-tabs -c -02 YOUR_PROGRAM_NAME.for
gfortran -Wall -Wno-tabs -llapack -lblas -lm -o EXECUTABLE_NAME
```

### Python

Python code defines wall (planar or spherical, `wall_pl.py` or `wall_sph.py`) as a `class`.
Most parameters are set in the constructor. 
The class can easily be imported to your code. 
The only dependency is the `numpy` library. 

## References

<a id="1">[1]</a>
Janek, J.; Pospíšil, M. and Malijevský A.: Phys. Rev. E (*in prep*).

<a id="2">[2]</a>
Rosenfeld, Y: Free-energy model for the inhomogeneous hard-sphere fluid mixture and density-functional theory of freezing, Phys. Rev. Lett. **63**, 980 (1989)

<a id="3">[3]</a>
GNU Scientific Library, [www.gnu.org/software/gsl](https://www.gnu.org/software/gsl/)

<a id="4">[4]</a>
Williams, J.: Polyroots Fortran, [GitHub repository](https://github.com/jacobwilliams/polyroots-fortran/)