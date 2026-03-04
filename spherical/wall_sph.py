import numpy as np
import planar.Python.wall_pl


class wall_sph:
    r"""
    ## Spherical wall class
        - wall inducing constant density profile in spherical geometry
        - using Rosenfeld's Fundamental Measure Theory (FMT) for hard-sphere contribution
        - optionally including Lennard-Jones (LJ) attractive interactions via mean-field approximation
    """
    def __init__(
        self,
        eta: float,
        varrho: float,
        rc: float = 2.5,
        beta: float = 1,
        epsilon: float = 1,
    ):
        r"""
        ## Spherical wall initialization
        - wall inducing constant density profile in spherical geometry
        - using Rosenfeld's Fundamental Measure Theory (FMT) for hard-sphere contribution
        - optionally including Lennard-Jones (LJ) attractive interactions via mean-field approximation
        - setting up parameters for hard-sphere fluid and for LJ interactions
        - setting up polynomials and their roots for direct correlation function calculation
        - setting up planar wall object for Taylor expansion
        
        :param eta: packing fraction of the hard-sphere fluid
        :type eta: float
        :param varrho: radius of the spherical wall
        :type varrho: float
        :param rc: cut-off radius for the LJ potential, defaults to 2.5
        :type rc: float, optional
        :param beta: inverse temperature, defaults to 1
        :type beta: float, optional
        :param epsilon: energy parameter for LJ potential, defaults to 1
        :type epsilon: float, optional
        """
        
        self.R = 0.5  # hard-coded parameter: radius of hard-spheres
        self.varrho = varrho  # spherical wall radius
        self.sigma = 1  # hard-coded parameter: diameter of hard-spheres
        self.eta = eta  # packing fraction
        self.rhob = 6 * eta / np.pi  # bulk number density
        self.beta = beta  # inverse temperature
        self.epsilon = epsilon  # energy parameter for LJ potential
        self.rc = rc  # cut-off radius for the LJ potential
        coeffs4a = (
            3.0 * eta * varrho**2 * (varrho + 1) ** 2,
            -2.0 * (1.0 + eta * varrho * (4.0 * varrho**2 + 6.0 * varrho + 3.0)),
            3.0 * eta * (1.0 + 2.0 * varrho + 2.0 * varrho**2),
            0.0,
            -eta,
        )  # coefficients of polynomial P_4^a
        coeffs4b = (
            2.0 * varrho,
            2.0,
            -3.0 * eta * (2.0 * varrho + 1.0),
            4.0 * varrho * eta,
            eta,
        )  # coefficients of polynomial P_4^b
        polyn4a = np.polynomial.Polynomial(coeffs4a)
        self.roots4a = polyn4a.roots()  # roots of polynomial P_4^a
        polyn4b = np.polynomial.Polynomial(coeffs4b)
        self.roots4b = polyn4b.roots()  # roots of polynomial P_4^b
        self.pl = planar.Python.wall_pl.wall_pl(
            eta, rc, beta, epsilon
        )  # planar wall object for Taylor expansion

    def n3(self, r: float) -> float:
        """
        ## Weighted density n_3 in spherical geometry
        - scalar weighted density

        :param r: radial distance from the center of the spherical wall
        :type r: float
        :return: weighted density n_3 at distance r
        :rtype: float
        """
        if (r >= self.varrho) and (r <= self.varrho + 1.0):
            return (
                self.eta
                / (2 * r)
                * (r - self.varrho) ** 2
                * (
                    3.0 * self.varrho**2
                    - r**2
                    - 2.0 * r * self.varrho
                    + 6.0 * self.varrho
                    + 3.0
                )
            )
        elif r > self.varrho + 1.0:
            return self.eta
        else:
            return 0

    def n2(self, r: float) -> float:
        r"""
        ## Weighted density n_2 in spherical geometry
        - scalar weighted density

        :param r: radial distance from the center of the spherical wall
        :type r: float
        :return: weighted density n_2 at distance r
        :rtype: float
        """
        if (r >= self.varrho) and (r <= self.varrho + 1.0):
            return 3.0 * self.eta / r * (r + self.varrho + 1.0) * (r - self.varrho)
        elif r > self.varrho + 1.0:
            return 6.0 * self.eta
        else:
            return 0

    def n1(self, r: float) -> float:
        r"""
        ## Weighted density n_1 in spherical geometry
        - scalar weighted density

        :param r: radial distance from the center of the spherical wall
        :type r: float
        :return: weighted density n_1 at distance r
        :rtype: float
        """
        if (r >= self.varrho) and (r <= self.varrho + 1.0):
            return (
                3.0
                * self.eta
                / (2.0 * np.pi * r)
                * (r + self.varrho + 1.0)
                * (r - self.varrho)
            )
        elif r > self.varrho + 1.0:
            return self.eta * 3.0 / np.pi
        else:
            return 0

    def n0(self, r: float) -> float:
        r"""
        ## Weighted density n_0 in spherical geometry
        - scalar weighted density

        :param r: radial distance from the center of the spherical wall
        :type r: float
        :return: weighted density n_0 at distance r
        :rtype: float
        """
        if (r >= self.varrho) and (r <= self.varrho + 1.0):
            return (
                3.0
                * self.eta
                / (np.pi * r)
                * (r + self.varrho + 1.0)
                * (r - self.varrho)
            )
        elif r > self.varrho + 1.0:
            return self.eta * 6.0 / np.pi
        else:
            return 0

    def N2(self, r: float) -> float:
        r"""
        ## Weighted density N_2 in spherical geometry
        - vector weighted density (magnitude/radial component)

        :param r: radial distance from the center of the spherical wall
        :type r: float
        :return: weighted density N_2 at distance r
        :rtype: float
        """
        if (r >= self.varrho) and (r <= self.varrho + 1.0):
            return (
                3.0
                * self.eta
                / (2.0 * r**2)
                * (r**2 - (self.varrho + 1.0) ** 2)
                * (r**2 - self.varrho**2)
            )
        else:
            return 0

    def N1(self, r: float) -> float:
        r"""
        ## Weighted density N_1 in spherical geometry
        - vector weighted density (magnitude/radial component)

        :param r: radial distance from the center of the spherical wall
        :type r: float
        :return: weighted density N_1 at distance r
        :rtype: float
        """
        if (r >= self.varrho) and (r <= self.varrho + 1.0):
            return (
                3.0
                * self.eta
                / (4.0 * np.pi * r**2)
                * (r**2 - (self.varrho + 1.0) ** 2)
                * (r**2 - self.varrho**2)
            )
        else:
            return 0

    def dphidn0(self, r: float) -> float:
        r"""
        ## Derivative of Phi with respect to n_0 in spherical geometry
        - Phi: the reduced excess free energy density
        - n_0: scalar weighted density n_0

        :param r: radial distance from the center of the spherical wall
        :type r: float
        :return: derivative of Phi with respect to n_0 at distance r
        :rtype: float
        """
        if r < 0:
            return 0
        return -np.log(1 - self.n3(r))

    def dphidn1(self, r: float) -> float:
        r"""
        ## Derivative of Phi with respect to n_1 in spherical geometry
        - Phi: the reduced excess free energy density
        - n_1: scalar weighted density n_1

        :param r: radial distance from the center of the spherical wall
        :type r: float
        :return: derivative of Phi with respect to n_1 at distance r
        :rtype: float
        """
        if r < 0:
            return 0
        return self.n2(r) / (1 - self.n3(r))

    def dphidn2(self, r: float) -> float:
        r"""
        ## Derivative of Phi with respect to n_2 in spherical geometry
        - Phi: the reduced excess free energy density
        - n_2: scalar weighted density n_2

        :param r: radial distance from the center of the spherical wall
        :type r: float
        :return: derivative of Phi with respect to n_2 at distance r
        :rtype: float
        """
        if r < 0:
            return 0
        x = 1 - self.n3(r)
        return self.n1(r) / x + (self.n2(r) ** 2 - self.N2(r) ** 2) / (8 * np.pi * x**2)

    def dphidn3(self, r: float) -> float:
        r"""
        ## Derivative of Phi with respect to n_3 in spherical geometry
        - Phi: the reduced excess free energy density
        - n_3: scalar weighted density n_3

        :param r: radial distance from the center of the spherical wall
        :type r: float
        :return: derivative of Phi with respect to n_3 at distance r
        :rtype: float
        """
        if r < 0:
            return 0
        x = 1 - self.n3(r)
        return (
            self.n0(r) / x
            + (self.n1(r) * self.n2(r) - self.N1(r) * self.N2(r)) / (x**2)
            + (self.n2(r) ** 3 - 3 * self.n2(r) * self.N2(r) ** 2) / (12 * np.pi * x**3)
        )

    def dphidN1(self, r: float) -> float:
        r"""
        ## Derivative of Phi with respect to N_1 in spherical geometry
        - Phi: the reduced excess free energy density
        - N_1: vector weighted density N_1 (magnitude/radial component)

        :param r: radial distance from the center of the spherical wall
        :type r: float
        :return: derivative of Phi with respect to N_1 at distance r
        :rtype: float
        """
        if r < 0:
            return 0
        return -self.N2(r) / (1 - self.n3(r))

    def dphidN2(self, r: float) -> float:
        r"""
        ## Derivative of Phi with respect to N_2 in spherical geometry
        - Phi: the reduced excess free energy density
        - N_2: vector weighted density N_2 (magnitude/radial component)

        :param r: radial distance from the center of the spherical wall
        :type r: float
        :return: derivative of Phi with respect to N_2 at distance r
        :rtype: float
        """
        if r < 0:
            return 0
        x = 1 - self.n3(r)
        return -self.N1(r) / x - self.n2(r) * self.N2(r) / (4 * np.pi * x**2)

    def __B0(self, t: complex, r: float) -> complex:
        """
        ## Summation part of c0 in spherical geometry
        - auxiliary function, not to be called directly
        - not including the logarithmic part of the summand

        :param t: complex root of the summation polynomial P_4^b
        :type t: complex
        :param r: radial distance from the center of the spherical wall
        :type r: float
        :return: summation part of c0 in spherical geometry
        :rtype: complex
        """
        rho = self.varrho  # abbreviation for the radius of the spherical wall
        eta = self.eta  # abbreviation for the packing fraction eta

        # numerator
        num = (
            (4 * rho * rho + 6 * rho + 3) * eta * t**3
            - 3 * t * t
            - 6 * t * rho
            - 4 * rho * rho
        )
        # denominator
        den = 2 * ((2 * t * t + 6 * (t - 1) * rho - 3) * eta * r * t) + 2 * r

        return num / den

    def c0(self, r: float) -> float:
        r"""
        ## c_0 component of the direct correlation function in spherical geometry
        - uses auxiliary function __B0

        :param r: radial distance from the center of the spherical wall
        :type r: float
        :return: c_0 component of the direct correlation function at distance r
        :rtype: float
        """

        rho = self.varrho
        eta = self.eta

        if r < rho:
            return 0.0

        elif rho <= r < rho + 0.5:
            # Term A
            pref = ((2 * r + 2 * rho + 1) * (2 * r - 2 * rho + 1)) / (16 * r)
            logarg = (
                (4 * r * r + 8 * r * rho - 12 * rho * rho - 20 * rho + 4 * r - 11)
                * eta
                * (2 * r - 2 * rho + 1) ** 2
            ) / (16 * (2 * r + 1)) + 1
            termA = pref * (2 * np.log(logarg) - 3)
            # Term B
            termB = (rho * rho) / (2 * r) * np.log(2 * rho / (2 * r + 1))

            # summation
            sumval = 0.0 + 0.0j
            for t in self.roots4b:
                sumval += self.__B0(t, r) * np.log(
                    (2 * r - 2 * t - 2 * rho + 1) / (-2 * t)
                )

            return termA + termB - sumval.real

        elif rho + 0.5 <= r < rho + 1.5:

            # Term A
            pref = -((2 * r + 2 * rho - 1) * (2 * r - 2 * rho - 1)) / (8 * r)
            logarg = (
                (4 * r * r - 12 * rho * rho + 8 * r * rho - 4 * r - 28 * rho - 11)
                * eta
                * (2 * r - 2 * rho - 1) ** 2
            ) / (16 * (2 * r - 1)) + 1
            termA = pref * np.log(logarg)

            # Term B
            termB = (3 * (2 * r + 2 * rho + 1) * (2 * r - 2 * rho - 3)) / (16 * r)

            # Term C
            termC = (rho * rho) / (2 * r) * np.log((2 * r - 1) / (2 * rho + 2))

            # Term D
            termD = (
                ((2 * r + 2 * rho + 1) * (2 * r - 2 * rho + 1))
                / (8 * r)
                * np.log(1 - eta)
            )

            # Summation
            sumval = 0.0 + 0.0j
            for t in self.roots4b:
                sumval += self.__B0(t, r) * np.log(
                    (2 * r - 2 * t - 2 * rho - 1) / (2 - 2 * t)
                )

            return termA + termB + termC + termD + sumval.real

        elif r >= rho + 1.5:
            return np.log(1 - eta)

        else:
            return 0.0

    def c0T(self, r: float) -> float:
        z = r - self.varrho
        eta = self.eta
        varrho = self.varrho
        T1 = 0.0
        pl = self.pl
        if z >= 1.5:
            return np.log(1 - eta)
        elif z >= 0.5:
            suma = 0.0
            for t in pl.roots3:
                suma += (
                    (7.0 * t**2 * eta - 6.0 * t - 1.0)
                    / (16.0 * t * eta * (t - 1.0))
                    * np.log((2.0 * z - 1.0 - 2.0 * t) / (1.0 - t))
                )
            T1 = (
                45.0 / 32.0
                + (2.0 * z - 1.0) ** 2
                / 8.0
                * np.log(2.0 + (4.0 * z**3 - 12.0 * z**2 + 9.0 * z - 2.0) * eta)
                - (8.0 * z**2 - 8.0 * z - 19.0) / 16.0 * np.log(2.0)
                - (2.0 * z + 1) ** 2 / 8.0 * np.log(1.0 - eta)
                - 9.0 * z**2 / 8.0
                + 3.0 * z / 4.0
                - suma
            )
        elif z >= 0:
            suma = 0
            for t in pl.roots3:
                suma += (
                    (7.0 * t**2 * eta - 6.0 * t - 1)
                    / (16.0 * t * eta * (t - 1))
                    * np.log(-t / (2.0 * z + 1.0 - 2.0 * t))
                )
            T1 = (
                15.0 / 32.0
                - (2.0 * z + 1.0) ** 2
                / 8.0
                * np.log(2.0 + (4.0 * z**3 - 3.0 * z - 1.0) * eta)
                + (8.0 * z**2 + 8.0 * z - 19.0) / 16.0 * np.log(2)
                + 9.0 * z**2 / 8.0
                + 3.0 * z / 2.0
                - suma
            )
        else:
            return 0
        return (1.0 - z / varrho) * pl.c0(z) - T1 / varrho

    def __B1(self, t: complex, r: float) -> complex:
        """
        ## Summation part of c1 in spherical geometry
        - auxiliary function, not to be called directly
        - not including the logarithmic part of the summand

        :param t: complex root of the summation polynomial P_4^b
        :type t: complex
        :param r: radial distance from the center of the spherical wall
        :type r: float
        :return: summation part of c1 in spherical geometry
        :rtype: complex
        """

        rho = self.varrho  # abbreviation for the radius of the spherical wall
        eta = self.eta  # abbreviation for the packing fraction eta

        den = (2 * t * t + 6 * (t - 1) * rho - 3) * eta * r * t + r
        num = (t + 2 * rho + 1) * (t + rho) * t

        return num / den

    def c1(self, r: float) -> float:
        r"""
        ## c_1 component of the direct correlation function in spherical geometry
        - uses auxiliary function __B1

        :param r: radial distance from the center of the spherical wall
        :type r: float
        :return: c_1 component of the direct correlation function at distance r
        :rtype: float
        """

        rho = self.varrho
        eta = self.eta

        if r < rho:
            return 0.0
        elif rho <= r < rho + 0.5:
            sum_val = 0.0 + 0.0j
            for t in self.roots4b:
                logarg = (2 * r - 2 * t - 2 * rho + 1) / (-2 * t)
                sum_val += self.__B1(t, r) * np.log(logarg)

            return -1.5 * eta * np.real(sum_val)

        elif rho + 0.5 <= r < rho + 1.5:
            termA = (3 * (2 * r + 2 * rho + 3) * (2 * r - 2 * rho - 1) * eta) / (
                8 * (eta - 1) * r
            )
            sum_val = 0.0 + 0.0j
            for t in self.roots4b:
                logarg = (2 * r - 2 * t - 2 * rho - 1) / (2 - 2 * t)
                sum_val += self.__B1(t, r) * np.log(logarg)

            return termA + 1.5 * eta * np.real(sum_val)

        elif r >= rho + 1.5:
            return 3 * eta / (eta - 1)

        else:
            return 0.0

    def c1T(self, r: float) -> float:
        z = r - self.varrho
        eta = self.eta
        varrho = self.varrho
        pl = self.pl
        T1 = 0.0
        if z >= 1.5:
            return -3 * eta / (1 - eta)
        elif z >= 0.5:
            suma = 0.0
            for t in pl.roots3:
                suma += (
                    (5.0 * t)
                    / (4.0 * t - 4.0)
                    * np.log((1.0 - t) / (2.0 * z - 2.0 * t - 1.0))
                )
            T1 = (
                15.0 / 4.0 * np.log(2.0)
                + 3
                * eta
                * (
                    (8.0 * z**3 - 8.0 * z**2 - 24.0 * z + 15.0) * (2.0 * z - 1.0) * eta
                    + 4.0 * z**2
                    + 15.0
                )
                * (2.0 * z - 1.0)
                / (8.0 * (2.0 + (2.0 * z - 1.0) ** 2 * (z - 2.0) * eta) * (1.0 - eta))
                + suma
            )
        elif z >= 0.0:
            suma = 0.0
            for t in pl.roots3:
                suma += (
                    5.0 * t / (4.0 * t - 4.0) * np.log((2.0 * z + 1.0 - 2.0 * t) / (-t))
                )
            T1 = (
                -15.0 / 4.0 * np.log(2.0)
                - (3.0 * eta * (2.0 * z - 1.0) * (2.0 * z + 1.0) ** 2)
                / (16.0 + (32.0 * z**3 - 24.0 * z - 8.0) * eta)
                + suma
            )
        else:
            return 0
        return (1.0 - z / varrho) * pl.c1(z) - 0.5 * T1 / varrho

    def __Q2(self, t: complex) -> complex:
        """
        ## Numerator of summand for c2 in spherical geometry
        - auxiliary function, not to be called directly

        :param t: complex root of the summation polynomial P_4^a
        :type t: complex
        :return: numerator of summand for c2 in spherical geometry
        :rtype: complex
        """
        R = self.varrho  # abbreviation for the radius of the spherical wall
        eta = self.eta  # abbreviation for the packing fraction eta

        # eta^3 term
        teta3 = (
            2
            * (
                8 * (4 * R**2 + 6 * R + 3) * t**2
                + (10 * R**2 + 15 * R + 6) * (2 * R + 3) * t
                - 52 * R**4
                - 156 * R**3
                - 177 * R**2
                - 90 * R
                - 18
            )
            * (R - t)
            * R
            * eta**3
        )

        # eta^2 term
        teta2 = (
            8 * (8 * R**3 + 12 * R**2 + 6 * R - 1) * t**3
            + 2 * (48 * R**3 + 72 * R**2 + 36 * R - 3) * t**2
            - 6 * (16 * R**5 + 40 * R**4 + 28 * R**3 - 4 * R**2 - 10 * R - 1) * t
            + (8 * R**3 + 12 * R**2 + 6 * R - 3) * (4 * R**2 + 6 * R + 3) * R
        ) * eta**2

        # eta^1 term
        teta1 = (
            8 * t**3
            + 12 * t**2
            - 6 * (2 * R**2 + 2 * R - 1) * t
            + 3 * (4 * R**2 + 6 * R + 3) * R
            - 1
        ) * eta

        return teta1 + teta2 + teta3 + 1

    def __B2(self, t: complex, r: float) -> complex:
        """
        ## Summation part of c2 in spherical geometry
        - auxiliary function, not to be called directly
        - not including the logarithmic part of the summand

        :param t: complex root of the summation polynomial P_4^a
        :type t: complex
        :param r: radial distance from the center of the spherical wall
        :type r: float
        """
        rho = self.varrho
        eta = self.eta
        termA = 2 * (t + 2 * rho) * (t - rho) - 6 * rho - 3
        termA = termA * eta * (t - rho) + 1
        termB = 2 * (4 * rho**2 + 6 * rho + 3) * eta * rho + 1
        return self.__Q2(t) / (8 * r * termA * termB * (eta - 1))

    def __S2(self, r: float) -> float:
        """
        ## Non-summation part of c2 in spherical geometry
        - auxiliary function, not to be called directly
        - numerator part only
        - for r in [varrho, varrho+1/2]

        :param r: radial distance from the center of the spherical wall
        :type r: float
        """
        R = self.varrho  # abbreviation for the radius of the spherical wall
        eta = self.eta  # abbreviation for the packing fraction eta
        # eta^3 term
        teta3 = (
            8
            * (
                12 * (2 * R**2 + 3 * R + 2) * (R + 1) * r**2
                + 6 * (2 * R**2 + 3 * R + 2) * (4 * R**2 + 8 * R + 5) * r
                - 3 * (12 * R**3 + 40 * R**2 + 49 * R + 22) * (2 * R + 1) * R
                + 6
            )
            * (2 * r - 2 * R + 1) ** 2
            * R
            * eta**3
        )

        # eta^2 term
        teta2 = (
            (
                24 * (8 * R**3 + 12 * R**2 + 8 * R + 3) * r**3
                + 12 * (16 * R**4 + 48 * R**3 + 52 * R**2 + 30 * R + 17) * r**2
                - 6
                * (32 * R**5 + 80 * R**4 + 248 * R**3 + 328 * R**2 + 164 * R - 13)
                * r
                - 3
                * (
                    64 * R**6
                    + 256 * R**5
                    + 352 * R**4
                    + 376 * R**3
                    + 352 * R**2
                    + 178 * R
                    + 1
                )
            )
            * (2 * r - 2 * R + 1)
            * eta**2
        )

        # eta^1 term
        teta1 = (
            3
            * (
                8 * r**3
                + 4 * (2 * R + 3) * r**2
                - (8 * R**2 + 8 * R + 34) * r
                - 8 * R**3
                - 20 * R**2
                - 14 * R
                - 19
            )
            * (2 * r - 2 * R + 1)
            * eta
        )

        return teta3 + teta2 + teta1

    def __T2(self, r: float) -> float:
        """
        ## Non-summation part of c2 in spherical geometry
        - auxiliary function, not to be called directly
        - numerator part only
        - for r in [varrho+1/2, varrho+3/2]

        :param r: radial distance from the center of the spherical wall
        :type r: float
        """

        R = self.varrho  # abbreviation for the radius of the spherical wall
        eta = self.eta  # abbreviation for the packing fraction eta

        # eta^4 term
        teta4 = (
            -(
                384 * (4 * R**2 + 6 * R + 3) * r**5 * (r - 1)
                - 96 * (56 * R**3 + 112 * R**2 + 88 * R + 21) * (2 * R + 1) * r**4
                + 192
                * (64 * R**5 + 208 * R**4 + 288 * R**3 + 192 * R**2 + 48 * R - 3)
                * r**3
                + 24
                * (48 * R**4 + 184 * R**3 + 184 * R**2 + 6 * R - 39)
                * (2 * R + 1) ** 2
                * r**2
                - 24
                * (128 * R**5 + 528 * R**4 + 728 * R**3 + 280 * R**2 - 146 * R - 105)
                * (2 * R + 1) ** 2
                * r
                + 6
                * (48 * R**4 + 136 * R**3 + 80 * R**2 - 46 * R - 45)
                * (2 * R + 3)
                * (2 * R + 1) ** 3
            )
            * R
            * eta**4
        )

        # eta^3 term
        teta3 = (
            -(
                192 * (16 * R**3 + 24 * R**2 + 12 * R + 1) * r**5 * (r - 1)
                - 48
                * (448 * R**5 + 1216 * R**4 + 1472 * R**3 + 932 * R**2 + 280 * R + 7)
                * r**4
                + 96
                * (
                    256 * R**6
                    + 832 * R**5
                    + 1200 * R**4
                    + 992 * R**3
                    + 496 * R**2
                    + 124 * R
                    - 1
                )
                * r**3
                + 12
                * (
                    768 * R**7
                    + 6016 * R**6
                    + 16704 * R**5
                    + 22896 * R**4
                    + 16784 * R**3
                    + 6144 * R**2
                    + 828 * R
                    - 13
                )
                * r**2
                - 12
                * (
                    2048 * R**8
                    + 13568 * R**7
                    + 37888 * R**6
                    + 58048 * R**5
                    + 52880 * R**4
                    + 28176 * R**3
                    + 7688 * R**2
                    + 636 * R
                    - 35
                )
                * r
                + 3
                * (
                    1536 * R**8
                    + 8960 * R**7
                    + 21888 * R**6
                    + 29408 * R**5
                    + 23696 * R**4
                    + 11184 * R**3
                    + 2600 * R**2
                    + 90 * R
                    - 15
                )
                * (2 * R + 3)
            )
            * eta**3
        )

        # eta^2 term
        teta2 = (
            -(
                384 * (r - 1) * r**5
                - 96 * (8 * R**3 + 40 * R**2 + 40 * R + 21) * r**4
                + 192 * (56 * R**3 + 88 * R**2 + 52 * R + 9) * r**3
                + 24
                * (64 * R**5 + 208 * R**4 + 688 * R**3 + 976 * R**2 + 544 * R + 71)
                * r**2
                - 24
                * (448 * R**5 + 1840 * R**4 + 3184 * R**3 + 2640 * R**2 + 968 * R + 65)
                * r
                - 6
                * (
                    64 * R**6
                    + 32 * R**5
                    - 592 * R**4
                    - 1296 * R**3
                    - 1108 * R**2
                    - 398 * R
                    - 15
                )
                * (2 * R + 3)
            )
            * eta**2
        )

        # eta^1 term
        teta1 = (
            96 * (r - 10) * r**3
            - 48 * (4 * R**2 + 4 * R + 15) * r**2
            + 48 * (20 * R**2 + 44 * R + 37) * r
            + 6 * (4 * R**2 - 4 * R - 11) * (2 * R + 3) ** 2
        ) * eta

        return teta4 + teta3 + teta2 + teta1

    def c2(self, r: float) -> float:
        r"""
        ## c_2 component of the direct correlation function in spherical geometry
        - uses auxiliary functions __Q2, __B2, __S2, __T2

        :param r: radial distance from the center of the spherical wall
        :type r: float
        :return: c_2 component of the direct correlation function at distance r
        :rtype: float
        """

        rho = self.varrho
        eta = self.eta
        roots4a = self.roots4a

        if r < rho:
            return 0.0
        elif r < (rho + 0.5):
            term1 = 4 * r**2 + 4 * r + 8 * r * rho - 12 * rho**2 - 20 * rho - 11
            term1 = term1 * eta * (2 * r - 2 * rho + 1) ** 2 + 16 * (2 * r + 1)
            term2 = 2 * (4 * rho**2 + 6 * rho + 3) * eta * rho + 1
            term3 = np.log((2 * r + 1) / 2 / rho) / (8 * r)

            nons = self.__S2(r) / (8 * r * term1 * term2 * (eta - 1)) + term3
            suma = 0.0 + 0.0j
            for t in roots4a:
                suma += self.__B2(t, r) * np.log(
                    (2 * r - 2 * t + 1) / (2 * rho - 2 * t)
                )
            return nons + suma.real
        elif r < (rho + 1.5):
            term1 = 4 * r * r - 4 * r + 8 * r * rho - 12 * rho * rho - 28 * rho - 11
            term2 = 2 * (4 * rho * rho + 6 * rho + 3) * eta * rho + 1
            term3 = (
                16.0
                * r
                * (term1 * eta * (2 * r - 2 * rho - 1) ** 2 + 16.0 * (2 * r - 1))
                * term2
                * (eta - 1) ** 2
            )
            nons = self.__T2(r) / term3 - np.log((2 * r - 1) / (2 * rho + 2)) / (8 * r)
            suma = 0.0 + 0.0j
            for t in roots4a:
                suma += self.__B2(t, r) * np.log(
                    (2 * r - 2 * t - 1) / (2 * rho - 2 * t + 2)
                )
            return nons - suma.real
        else:
            return -3.0 * eta * (2.0 + eta) / (2.0 * (1.0 - eta) ** 2)

    def __Q3(self, t: complex, r: float) -> complex:
        """
        ## Numerator of summand for c3 in spherical geometry
        - auxiliary function, not to be called directly

        :param t: complex root of the summation polynomial P_4^a
        :type t: complex
        :param r: radial distance from the center of the spherical wall
        :type r: float
        :return: complex value of the numerator
        :rtype: complex
        """

        rho = self.varrho
        eta = self.eta

        # first term
        termA = (
            2
            * (
                4 * (4 * rho**2 + 6 * rho + 3) ** 2 * rho**2 * eta**4
                - 4
                * (4 * rho**2 + 6 * rho + 3)
                * (8 * rho**3 + 12 * rho**2 + 6 * rho - 1)
                * rho
                * eta**3
                + (
                    64 * rho**6
                    + 192 * rho**5
                    + 240 * rho**4
                    + 112 * rho**3
                    - 12 * rho**2
                    - 24 * rho
                    + 1
                )
                * eta**2
                + 2 * (8 * rho**3 + 12 * rho**2 + 6 * rho - 1) * eta
                + 1
            )
            * (
                4 * (2 * rho + 1) ** 2 * r**2
                - 96 * (rho + 1) * rho * r
                - 48 * rho**4
                - 96 * rho**3
                - 16 * rho**2
                + 32 * rho
                - 1
            )
            * eta
            * t**3
        )

        # eta^5 term
        teta5 = (
            4
            * (
                4
                * (328 * rho**4 + 948 * rho**3 + 1122 * rho**2 + 639 * rho + 153)
                * r**2
                + 8
                * (
                    88 * rho**5
                    - 52 * rho**4
                    - 590 * rho**3
                    - 879 * rho**2
                    - 546 * rho
                    - 135
                )
                * r
                + (
                    3072 * rho**7
                    + 22560 * rho**6
                    + 61584 * rho**5
                    + 94112 * rho**4
                    + 90600 * rho**3
                    + 56442 * rho**2
                    + 21285 * rho
                    + 3843
                )
            )
            * (rho + 1)
            * rho**3
            * t**2
            - 8
            * (
                12
                * (
                    64 * rho**8
                    + 344 * rho**7
                    + 820 * rho**6
                    + 1122 * rho**5
                    + 979 * rho**4
                    + 576 * rho**3
                    + 234 * rho**2
                    + 63 * rho
                    + 9
                )
                * r**2
                - 8
                * (
                    104 * rho**6
                    + 596 * rho**5
                    + 1198 * rho**4
                    + 1143 * rho**3
                    + 522 * rho**2
                    + 81 * rho
                    - 9
                )
                * (rho + 1)
                * rho
                * r
                + 768 * rho**10
                + 11296 * rho**9
                + 50032 * rho**8
                + 115248 * rho**7
                + 161352 * rho**6
                + 144018 * rho**5
                + 80763 * rho**4
                + 25920 * rho**3
                + 3402 * rho**2
                - 189 * rho
                - 27
            )
            * rho**2
            * t
            + 4
            * (
                4
                * (
                    256 * rho**8
                    + 1480 * rho**7
                    + 3964 * rho**6
                    + 6390 * rho**5
                    + 6849 * rho**4
                    + 5040 * rho**3
                    + 2493 * rho**2
                    + 756 * rho
                    + 108
                )
                * r**2
                - 24
                * (56 * rho**4 + 140 * rho**3 + 106 * rho**2 + 21 * rho - 6)
                * (rho + 1) ** 3
                * rho
                * r
                + 5664 * rho**9
                + 34736 * rho**8
                + 93488 * rho**7
                + 142856 * rho**6
                + 133218 * rho**5
                + 74799 * rho**4
                + 22176 * rho**3
                + 1503 * rho**2
                - 756 * rho
                - 108
            )
            * rho**3
        ) * eta**5

        # eta^4 term
        teta4 = (
            -(
                -(
                    4
                    * (
                        576 * rho**6
                        - 256 * rho**5
                        - 3504 * rho**4
                        - 5280 * rho**3
                        - 3012 * rho**2
                        - 396 * rho
                        + 207
                    )
                    * r**2
                    - 8
                    * (
                        576 * rho**7
                        + 3520 * rho**6
                        + 5360 * rho**5
                        + 1600 * rho**4
                        - 3308 * rho**3
                        - 3052 * rho**2
                        - 611 * rho
                        + 180
                    )
                    * r
                    - (
                        22272 * rho**8
                        + 167424 * rho**7
                        + 462336 * rho**6
                        + 710528 * rho**5
                        + 676320 * rho**4
                        + 405264 * rho**3
                        + 135072 * rho**2
                        + 12564 * rho
                        - 5121
                    )
                )
                * (rho + 1)
                * rho**2
                * t**2
                - 2
                * (
                    12
                    * (
                        704 * rho**9
                        + 3904 * rho**8
                        + 9584 * rho**7
                        + 13424 * rho**6
                        + 11724 * rho**5
                        + 6548 * rho**4
                        + 2289 * rho**3
                        + 437 * rho**2
                        + 12 * rho
                        - 12
                    )
                    * r**2
                    - 8
                    * (
                        576 * rho**8
                        + 4160 * rho**7
                        + 13040 * rho**6
                        + 20992 * rho**5
                        + 18348 * rho**4
                        + 8088 * rho**3
                        + 1173 * rho**2
                        - 168 * rho
                        + 12
                    )
                    * (rho + 1)
                    * rho
                    * r
                    + 8448 * rho**11
                    + 107776 * rho**10
                    + 457216 * rho**9
                    + 1025664 * rho**8
                    + 1388864 * rho**7
                    + 1163552 * rho**6
                    + 562344 * rho**5
                    + 108240 * rho**4
                    - 23427 * rho**3
                    - 11823 * rho**2
                    - 36 * rho
                    + 36
                )
                * rho
                * t
                + (
                    4
                    * (
                        3776 * rho**9
                        + 21056 * rho**8
                        + 53168 * rho**7
                        + 78832 * rho**6
                        + 75172 * rho**5
                        + 46864 * rho**4
                        + 17865 * rho**3
                        + 2859 * rho**2
                        - 585 * rho
                        - 279
                    )
                    * r**2
                    - 24
                    * (
                        576 * rho**6
                        + 2240 * rho**5
                        + 3440 * rho**4
                        + 2368 * rho**3
                        + 596 * rho**2
                        - 52 * rho
                        + 1
                    )
                    * (rho + 1) ** 3
                    * rho
                    * r
                    + 6912 * rho**11
                    + 92928 * rho**10
                    + 419584 * rho**9
                    + 985216 * rho**8
                    + 1379680 * rho**7
                    + 1191248 * rho**6
                    + 603104 * rho**5
                    + 139292 * rho**4
                    - 8973 * rho**3
                    - 8007 * rho**2
                    + 585 * rho
                    + 279
                )
                * rho**2
            )
            * eta**4
        )

        # eta^3 term
        teta3 = (
            (
                4
                * (
                    1536 * rho**6
                    + 4608 * rho**5
                    + 5832 * rho**4
                    + 3132 * rho**3
                    + 306 * rho**2
                    - 279 * rho
                    + 17
                )
                * r**2
                - 8
                * (
                    1920 * rho**6
                    + 5976 * rho**5
                    + 7892 * rho**4
                    + 4546 * rho**3
                    + 647 * rho**2
                    - 324 * rho
                    + 15
                )
                * r
                + (
                    12288 * rho**9
                    + 92160 * rho**8
                    + 258048 * rho**7
                    + 395040 * rho**6
                    + 353520 * rho**5
                    + 177408 * rho**4
                    + 26712 * rho**3
                    - 18942 * rho**2
                    - 9573 * rho
                    + 427
                )
            )
            * (rho + 1)
            * rho
            * t**2
            - (
                (
                    6144 * rho**10
                    + 30720 * rho**9
                    + 66048 * rho**8
                    + 73536 * rho**7
                    + 37728 * rho**6
                    - 3504 * rho**5
                    - 15672 * rho**4
                    - 9000 * rho**3
                    - 2712 * rho**2
                    - 504 * rho
                    + 24
                )
                * r**2
                + 16
                * (
                    1224 * rho**6
                    + 4260 * rho**5
                    + 6258 * rho**4
                    + 4335 * rho**3
                    + 1199 * rho**2
                    - 25 * rho
                    + 1
                )
                * (rho + 1)
                * rho
                * r
                + (
                    6144 * rho**12
                    + 86016 * rho**11
                    + 365568 * rho**10
                    + 795456 * rho**9
                    + 968416 * rho**8
                    + 573184 * rho**7
                    - 45744 * rho**6
                    - 307484 * rho**5
                    - 182498 * rho**4
                    - 30390 * rho**3
                    + 4950 * rho**2
                    + 126 * rho
                    - 6
                )
            )
            * t
            + (
                4
                * (
                    1024 * rho**10
                    + 5632 * rho**9
                    + 14848 * rho**8
                    + 22952 * rho**7
                    + 21428 * rho**6
                    + 10262 * rho**5
                    - 577 * rho**4
                    - 4018 * rho**3
                    - 2233 * rho**2
                    - 378 * rho
                    + 60
                )
                * r**2
                + 24
                * (216 * rho**3 + 388 * rho**2 + 226 * rho - 1)
                * (rho + 1) ** 3
                * rho**2
                * r
                + (
                    18432 * rho**11
                    + 109568 * rho**10
                    + 283808 * rho**9
                    + 400592 * rho**8
                    + 304432 * rho**7
                    + 79480 * rho**6
                    - 52286 * rho**5
                    - 43487 * rho**4
                    - 5918 * rho**3
                    + 2653 * rho**2
                    + 378 * rho
                    - 60
                )
            )
            * rho
        ) * eta**3

        # eta^2 term
        teta2 = (
            (
                4 * (384 * rho**3 + 576 * rho**2 + 288 * rho - 23) * r**2
                - 24 * (160 * rho**3 + 246 * rho**2 + 126 * rho - 9) * r
                + (
                    6144 * rho**6
                    + 27648 * rho**5
                    + 46080 * rho**4
                    + 41088 * rho**3
                    + 19212 * rho**2
                    + 3948 * rho
                    - 853
                )
            )
            * (rho + 1)
            * rho
            * t**2
            + (
                (
                    -1536 * rho**7
                    - 5376 * rho**6
                    - 7296 * rho**5
                    - 4272 * rho**4
                    - 864 * rho**3
                    - 24 * rho**2
                    - 72 * rho
                    + 48
                )
                * r**2
                - 16
                * (288 * rho**3 + 481 * rho**2 + 265 * rho - 1)
                * (rho + 1)
                * rho
                * r
                - (
                    1536 * rho**9
                    + 31488 * rho**8
                    + 110592 * rho**7
                    + 176112 * rho**6
                    + 138992 * rho**5
                    + 38948 * rho**4
                    - 13896 * rho**3
                    - 9038 * rho**2
                    + 262 * rho
                    + 12
                )
            )
            * t
            + (
                (
                    2048 * rho**8
                    + 8192 * rho**7
                    + 15872 * rho**6
                    + 17792 * rho**5
                    + 11804 * rho**4
                    + 3608 * rho**3
                    - 308 * rho**2
                    - 464 * rho
                    + 16
                )
                * r**2
                + 432 * (rho + 1) ** 3 * rho**3 * r
                + (
                    1152 * rho**7
                    + 3904 * rho**6
                    + 4800 * rho**5
                    + 1929 * rho**4
                    - 614 * rho**3
                    - 435 * rho**2
                    + 132 * rho
                    - 4
                )
                * (2 * rho + 1) ** 2
            )
        ) * eta**2

        # eta^1 term
        teta1 = (
            48
            * (2 * r**2 - 5 * r + 20 * rho**3 + 42 * rho**2 + 27 * rho + 9)
            * (rho + 1)
            * rho
            * t**2
            - 2
            * (
                12 * (4 * rho**4 + 8 * rho**3 + 4 * rho**2 + 1) * r**2
                + 144 * (rho + 1) * rho * r
                + (
                    48 * rho**6
                    + 2064 * rho**5
                    + 4932 * rho**4
                    + 4344 * rho**3
                    + 1120 * rho**2
                    - 308 * rho
                    - 3
                )
            )
            * t
            + (
                16
                * (20 * rho**3 + 30 * rho**2 + 15 * rho - 2)
                * (rho**2 + rho + 1)
                * r**2
                + (
                    288 * rho**6
                    + 784 * rho**5
                    + 664 * rho**4
                    + 28 * rho**3
                    - 172 * rho**2
                    - 52 * rho
                    + 8
                )
            )
        ) * eta

        # eta^0 term
        teta0 = 48 * (rho + 1) * rho * t * (t - 4) + 4 * (rho**2 + rho + 1) * (
            2 * r + 1
        ) * (2 * r - 1)

        return termA + teta5 + teta4 + teta3 + teta2 + teta1 + teta0

    def __B3(self, t: complex, r: float) -> complex:
        """
        ## Summation part of c3 in spherical geometry
        - auxiliary function, not to be called directly
        - not including the logarithmic part of the summand

        :param t: complex root of the summation polynomial P_4^a
        :type t: complex
        :param r: radial distance from the center of the spherical wall
        :type r: float
        """

        rho = self.varrho  # abbreviation for the radius of the spherical wall
        eta = self.eta  # abbreviation for the packing fraction eta

        termA = 2 * (t + 2 * rho) * (t - rho) - 6 * rho - 3
        termA = termA * eta * (t - rho) + 1
        termB = (2 * (4 * rho**2 + 6 * rho + 3) * eta * rho + 1) ** 2
        return self.__Q3(t, r) / (
            32 * r * rho * (rho + 1) * termA * termB * (eta - 1) ** 2
        )

    def __S3(self, r: float) -> float:
        """
        ## Non-summation part of c3 in spherical geometry
        - auxiliary function, not to be called directly
        - numerator part only
        - for r in [varrho, varrho+1/2]

        :param r: radial distance from the center of the spherical wall
        :type r: float
        """

        rho = self.varrho  # abbreviation for the radius of the spherical wall
        eta = self.eta  # abbreviation for the packing fraction eta
        prefactor = 2 * r - 2 * rho + 1  # common prefactor

        # eta^6 term
        teta6 = (
            -16
            * (4 * r * (r + 2 * rho + 1) - 12 * rho**2 - 20 * rho - 11)
            * (4 * rho**2 + 6 * rho + 3) ** 2
            * (2 * r - 2 * rho + 1) ** 2
            * (2 * r - 2 * rho - 1)
            * (2 * rho + 1)
            * rho**2
            * eta**6
        )

        # eta^5 term
        teta5 = (
            4
            * (
                32
                * (
                    320 * rho**6
                    + 976 * rho**5
                    + 1200 * rho**4
                    + 636 * rho**3
                    + 72 * rho**2
                    - 54 * rho
                    - 9
                )
                * r**5
                - 16
                * (
                    640 * rho**7
                    + 1376 * rho**6
                    - 1136 * rho**5
                    - 6440 * rho**4
                    - 8236 * rho**3
                    - 4728 * rho**2
                    - 1056 * rho
                    + 27
                )
                * r**4
                - 16
                * (
                    3840 * rho**8
                    + 16832 * rho**7
                    + 32384 * rho**6
                    + 30592 * rho**5
                    + 9168 * rho**4
                    - 8748 * rho**3
                    - 9072 * rho**2
                    - 2718 * rho
                    - 45
                )
                * r**3
                + 8
                * (
                    17920 * rho**9
                    + 78464 * rho**8
                    + 128832 * rho**7
                    + 66112 * rho**6
                    - 76896 * rho**5
                    - 146856 * rho**4
                    - 101124 * rho**3
                    - 33336 * rho**2
                    - 4176 * rho
                    + 63
                )
                * r**2
                - 2
                * (
                    56320 * rho**10
                    + 257792 * rho**9
                    + 385792 * rho**8
                    + 46016 * rho**7
                    - 546624 * rho**6
                    - 686224 * rho**5
                    - 298304 * rho**4
                    + 48980 * rho**3
                    + 91896 * rho**2
                    + 26214 * rho
                    + 81
                )
                * r
                + (
                    30720 * rho**11
                    + 147968 * rho**10
                    + 202496 * rho**9
                    - 70784 * rho**8
                    - 448064 * rho**7
                    - 374048 * rho**6
                    + 47440 * rho**5
                    + 230808 * rho**4
                    + 99516 * rho**3
                    - 12024 * rho**2
                    - 13536 * rho
                    - 99
                )
            )
            * rho
            * eta**5
        )

        # eta^4 term
        teta4 = (
            32
            * (
                384 * rho**7
                + 3456 * rho**6
                + 8832 * rho**5
                + 10992 * rho**4
                + 6864 * rho**3
                + 1812 * rho**2
                - 28 * rho
                - 5
            )
            * r**5
            - 16
            * (
                768 * rho**8
                + 4992 * rho**7
                + 20864 * rho**6
                + 41888 * rho**5
                + 43600 * rho**4
                + 21352 * rho**3
                + 2220 * rho**2
                - 1482 * rho
                + 15
            )
            * r**4
            - 16
            * (
                4608 * rho**9
                + 47616 * rho**8
                + 161664 * rho**7
                + 317248 * rho**6
                + 393792 * rho**5
                + 310176 * rho**4
                + 138336 * rho**3
                + 23792 * rho**2
                - 3220 * rho
                - 25
            )
            * r**3
            + 8
            * (
                21504 * rho**10
                + 207360 * rho**9
                + 862464 * rho**8
                + 1933824 * rho**7
                + 2565824 * rho**6
                + 2048256 * rho**5
                + 922336 * rho**4
                + 177152 * rho**3
                - 10920 * rho**2
                - 5474 * rho
                + 35
            )
            * r**2
            - 2
            * (
                67584 * rho**11
                + 632832 * rho**10
                + 3058688 * rho**9
                + 7361792 * rho**8
                + 9771392 * rho**7
                + 6720960 * rho**6
                + 1059264 * rho**5
                - 1829952 * rho**4
                - 1292528 * rho**3
                - 243324 * rho**2
                + 33516 * rho
                + 45
            )
            * r
            + (
                36864 * rho**12
                + 337920 * rho**11
                + 1855488 * rho**10
                + 4647424 * rho**9
                + 5665792 * rho**8
                + 2344192 * rho**7
                - 1824448 * rho**6
                - 2218048 * rho**5
                - 203296 * rho**4
                + 609032 * rho**3
                + 221516 * rho**2
                - 17638 * rho
                - 55
            )
        ) * eta**4

        # eta^3 term
        teta3 = (
            -2
            * (
                32
                * (
                    192 * rho**6
                    + 576 * rho**5
                    + 672 * rho**4
                    + 36 * rho**3
                    - 444 * rho**2
                    - 281 * rho
                    - 1
                )
                * r**5
                + 16
                * (
                    384 * rho**7
                    - 960 * rho**6
                    - 4800 * rho**5
                    - 6648 * rho**4
                    - 2852 * rho**3
                    + 834 * rho**2
                    + 885 * rho
                    - 67
                )
                * r**4
                - 16
                * (
                    2304 * rho**8
                    + 8448 * rho**7
                    + 19008 * rho**6
                    + 24432 * rho**5
                    + 15792 * rho**4
                    - 2544 * rho**3
                    - 9836 * rho**2
                    - 4817 * rho
                    + 119
                )
                * r**3
                + 8
                * (
                    1536 * rho**9
                    + 39168 * rho**8
                    + 138624 * rho**7
                    + 203232 * rho**6
                    + 112560 * rho**5
                    - 40720 * rho**4
                    - 83096 * rho**3
                    - 34818 * rho**2
                    - 2509 * rho
                    + 243
                )
                * r**2
                + 2
                * (
                    15360 * rho**10
                    - 125952 * rho**9
                    - 577536 * rho**8
                    - 904896 * rho**7
                    - 322816 * rho**6
                    + 722512 * rho**5
                    + 1019792 * rho**4
                    + 459692 * rho**3
                    + 2836 * rho**2
                    - 49481 * rho
                    + 1359
                )
                * r
                - (
                    18432 * rho**11
                    - 70656 * rho**10
                    - 405504 * rho**9
                    - 649344 * rho**8
                    - 275520 * rho**7
                    + 141664 * rho**6
                    - 195376 * rho**5
                    - 712936 * rho**4
                    - 530084 * rho**3
                    - 94034 * rho**2
                    + 35995 * rho
                    - 701
                )
            )
            * eta**3
        )

        # eta^2 term
        teta2 = (
            -(
                96 * (32 * rho**3 + 48 * rho**2 + 22 * rho - 15) * r**5
                + 16
                * (192 * rho**4 - 768 * rho**3 - 1428 * rho**2 - 768 * rho + 149)
                * r**4
                - 16
                * (
                    1152 * rho**5
                    + 2496 * rho**4
                    + 5112 * rho**3
                    + 4692 * rho**2
                    + 1722 * rho
                    - 785
                )
                * r**3
                + 8
                * (
                    3840 * rho**6
                    + 27648 * rho**5
                    + 53616 * rho**4
                    + 39264 * rho**3
                    + 4248 * rho**2
                    - 6312 * rho
                    - 577
                )
                * r**2
                + 2
                * (
                    19968 * rho**7
                    - 123648 * rho**6
                    - 397920 * rho**5
                    - 439056 * rho**4
                    - 100112 * rho**3
                    + 118600 * rho**2
                    + 72082 * rho
                    - 8117
                )
                * r
                - (
                    9216 * rho**8
                    - 61440 * rho**7
                    - 82368 * rho**6
                    + 136704 * rho**5
                    + 360912 * rho**4
                    + 171200 * rho**3
                    - 84340 * rho**2
                    - 82096 * rho
                    + 5587
                )
            )
            * eta**2
        )

        # eta^1 term
        teta1 = (
            -2
            * (
                96 * r**5
                + 48 * (2 * rho - 11) * r**4
                - 48 * (12 * rho**2 + 8 * rho + 33) * r**3
                + 24 * (136 * rho**3 + 372 * rho**2 + 262 * rho - 29) * r**2
                + 2
                * (1776 * rho**4 - 11136 * rho**3 - 16968 * rho**2 - 7824 * rho + 3043)
                * r
                - (
                    288 * rho**5
                    - 3504 * rho**4
                    + 7632 * rho**3
                    + 16968 * rho**2
                    + 9450 * rho
                    - 3055
                )
            )
            * eta
        )

        # eta^0 term
        teta0 = -96 * (2 * r + 2 * rho - 15) * (2 * r + 1)

        return prefactor * (teta6 + teta5 + teta4 + teta3 + teta2 + teta1 + teta0)

    def __T3(self, r: float) -> float:
        """
        ## Non-summation part of c3 in spherical geometry
        - auxiliary function, not to be called directly
        - numerator part only
        - for r in [varrho+1/2, varrho+3/2]

        :param r: radial distance from the center of the spherical wall
        :type r: float
        """

        rho = self.varrho  # abbreviation for the radius of the spherical wall
        eta = self.eta  # abbreviation for the packing fraction eta

        # eta^6 term
        teta6 = (
            16
            * (
                64 * (2 - r) * (4 * rho**2 + 6 * rho + 3) ** 2 * r**7
                + 16
                * (
                    768 * rho**6
                    + 3200 * rho**5
                    + 6256 * rho**4
                    + 7248 * rho**3
                    + 5232 * rho**2
                    + 2196 * rho
                    + 423
                )
                * r**6
                - 32
                * (
                    512 * rho**7
                    + 3008 * rho**6
                    + 7760 * rho**5
                    + 12536 * rho**4
                    + 13924 * rho**3
                    + 10446 * rho**2
                    + 4710 * rho
                    + 981
                )
                * r**5
                - 4
                * (
                    7680 * rho**8
                    + 35840 * rho**7
                    + 83264 * rho**6
                    + 117888 * rho**5
                    + 88080 * rho**4
                    + 1632 * rho**3
                    - 55692 * rho**2
                    - 42768 * rho
                    - 11241
                )
                * r**4
                + 8
                * (
                    12288 * rho**9
                    + 76544 * rho**8
                    + 226688 * rho**7
                    + 434176 * rho**6
                    + 587328 * rho**5
                    + 559456 * rho**4
                    + 356376 * rho**3
                    + 137616 * rho**2
                    + 25848 * rho
                    + 855
                )
                * r**3
                - (
                    102400 * rho**10
                    + 739328 * rho**9
                    + 2481920 * rho**8
                    + 5415680 * rho**7
                    + 8795008 * rho**6
                    + 10944192 * rho**5
                    + 10115488 * rho**4
                    + 6578384 * rho**3
                    + 2805000 * rho**2
                    + 696876 * rho
                    + 75843
                )
                * r**2
                + 2
                * (
                    6144 * rho**9
                    + 43776 * rho**8
                    + 140224 * rho**7
                    + 295328 * rho**6
                    + 494784 * rho**5
                    + 667552 * rho**4
                    + 660492 * rho**3
                    + 432246 * rho**2
                    + 165654 * rho
                    + 28341
                )
                * (2 * rho + 1) ** 2
                * r
                - (
                    1152 * rho**9
                    + 7488 * rho**8
                    + 21376 * rho**7
                    + 42816 * rho**6
                    + 77096 * rho**5
                    + 113852 * rho**4
                    + 116712 * rho**3
                    + 74748 * rho**2
                    + 26793 * rho
                    + 4059
                )
                * (2 * rho + 3)
                * (2 * rho + 1) ** 2
            )
            * rho**2
            * eta**6
        )

        # eta^5 term
        teta5 = (
            2
            * (
                512
                * (2 - r)
                * (4 * rho**2 + 6 * rho + 3)
                * (4 * rho**2 + 2 * rho + 1)
                * (rho + 1)
                * r**7
                + 64
                * (
                    1536 * rho**7
                    + 7232 * rho**6
                    + 12832 * rho**5
                    + 10704 * rho**4
                    + 3552 * rho**3
                    + 72 * rho**2
                    + 204 * rho
                    + 285
                )
                * r**6
                - 64
                * (
                    2048 * rho**8
                    + 13696 * rho**7
                    + 36352 * rho**6
                    + 40384 * rho**5
                    + 7456 * rho**4
                    - 24504 * rho**3
                    - 20208 * rho**2
                    - 3514 * rho
                    + 1317
                )
                * r**5
                - 16
                * (
                    15360 * rho**9
                    + 88320 * rho**8
                    + 167552 * rho**7
                    + 141056 * rho**6
                    + 149152 * rho**5
                    + 336624 * rho**4
                    + 443824 * rho**3
                    + 266444 * rho**2
                    + 53604 * rho
                    - 7491
                )
                * r**4
                + 32
                * (
                    24576 * rho**10
                    + 186368 * rho**9
                    + 548864 * rho**8
                    + 743680 * rho**7
                    + 332288 * rho**6
                    - 312928 * rho**5
                    - 459552 * rho**4
                    - 172392 * rho**3
                    + 21564 * rho**2
                    + 22860 * rho
                    + 603
                )
                * r**3
                - 4
                * (
                    204800 * rho**11
                    + 1811456 * rho**10
                    + 6457856 * rho**9
                    + 10749696 * rho**8
                    + 4965632 * rho**7
                    - 12269760 * rho**6
                    - 25457312 * rho**5
                    - 22228640 * rho**4
                    - 10114944 * rho**3
                    - 2004320 * rho**2
                    + 48548 * rho
                    + 50661
                )
                * r**2
                + 4
                * (
                    98304 * rho**12
                    + 985088 * rho**11
                    + 4075520 * rho**10
                    + 8044544 * rho**9
                    + 4541952 * rho**8
                    - 13309440 * rho**7
                    - 35129088 * rho**6
                    - 40819168 * rho**5
                    - 27170736 * rho**4
                    - 10160248 * rho**3
                    - 1636536 * rho**2
                    + 81714 * rho
                    + 37731
                )
                * r
                - (
                    36864 * rho**12
                    + 356352 * rho**11
                    + 1398784 * rho**10
                    + 2332672 * rho**9
                    - 331264 * rho**8
                    - 8488960 * rho**7
                    - 16720576 * rho**6
                    - 17258112 * rho**5
                    - 10497696 * rho**4
                    - 3555024 * rho**3
                    - 473100 * rho**2
                    + 46152 * rho
                    + 10791
                )
                * (2 * rho + 3)
            )
            * rho
            * eta**5
        )

        # eta^4 term
        teta4 = (
            4
            * (
                64 * (2 - r) * (2 * rho + 1) ** 6 * r**7
                + 16
                * (
                    3072 * rho**8
                    + 14208 * rho**7
                    + 32640 * rho**6
                    + 45600 * rho**5
                    + 40128 * rho**4
                    + 19944 * rho**3
                    + 3936 * rho**2
                    - 416 * rho
                    + 47
                )
                * r**6
                - 32
                * (
                    2048 * rho**9
                    + 13440 * rho**8
                    + 40896 * rho**7
                    + 83232 * rho**6
                    + 114352 * rho**5
                    + 100080 * rho**4
                    + 47040 * rho**3
                    + 6522 * rho**2
                    - 2315 * rho
                    + 109
                )
                * r**5
                - 4
                * (
                    30720 * rho**10
                    + 171520 * rho**9
                    + 489984 * rho**8
                    + 782336 * rho**7
                    + 523264 * rho**6
                    - 278016 * rho**5
                    - 750304 * rho**4
                    - 448248 * rho**3
                    - 34052 * rho**2
                    + 44872 * rho
                    - 1249
                )
                * r**4
                + 8
                * (
                    49152 * rho**11
                    + 362496 * rho**10
                    + 1310208 * rho**9
                    + 3083520 * rho**8
                    + 4941056 * rho**7
                    + 5278912 * rho**6
                    + 3522528 * rho**5
                    + 1257904 * rho**4
                    + 124792 * rho**3
                    - 26752 * rho**2
                    + 6858 * rho
                    + 95
                )
                * r**3
                - (
                    409600 * rho**12
                    + 3520512 * rho**11
                    + 14518272 * rho**10
                    + 39896576 * rho**9
                    + 79227904 * rho**8
                    + 112583424 * rho**7
                    + 109728640 * rho**6
                    + 68306784 * rho**5
                    + 23047216 * rho**4
                    + 1666408 * rho**3
                    - 1144440 * rho**2
                    - 166128 * rho
                    + 8427
                )
                * r**2
                + 2
                * (
                    98304 * rho**13
                    + 956416 * rho**12
                    + 4396032 * rho**11
                    + 13556224 * rho**10
                    + 31287040 * rho**9
                    + 53866368 * rho**8
                    + 66479040 * rho**7
                    + 55691456 * rho**6
                    + 28814144 * rho**5
                    + 6900176 * rho**4
                    - 715544 * rho**3
                    - 686514 * rho**2
                    - 52929 * rho
                    + 3149
                )
                * r
                - (
                    18432 * rho**13
                    + 172032 * rho**12
                    + 751104 * rho**11
                    + 2294016 * rho**10
                    + 5449472 * rho**9
                    + 9633856 * rho**8
                    + 11960512 * rho**7
                    + 9817680 * rho**6
                    + 4782872 * rho**5
                    + 907196 * rho**4
                    - 257442 * rho**3
                    - 137050 * rho**2
                    - 4813 * rho
                    + 451
                )
                * (2 * rho + 3)
            )
            * eta**4
        )

        # eta^3 term
        teta3 = (
            256 * (2 - r) * (16 * rho**3 + 24 * rho**2 + 12 * rho + 1) * r**7
            - 128
            * (
                192 * rho**6
                + 192 * rho**5
                - 480 * rho**4
                - 1596 * rho**3
                - 1638 * rho**2
                - 676 * rho
                + 19
            )
            * r**6
            + 128
            * (
                1344 * rho**6
                + 2976 * rho**5
                + 1008 * rho**4
                - 5784 * rho**3
                - 7834 * rho**2
                - 3437 * rho
                + 198
            )
            * r**5
            + 64
            * (
                2688 * rho**8
                + 6144 * rho**7
                - 8160 * rho**6
                - 43080 * rho**5
                - 54644 * rho**4
                - 16706 * rho**3
                + 15085 * rho**2
                + 10263 * rho
                - 929
            )
            * r**4
            - 32
            * (
                6144 * rho**9
                + 56832 * rho**8
                + 130560 * rho**7
                + 48384 * rho**6
                - 267264 * rho**5
                - 486560 * rho**4
                - 345992 * rho**3
                - 92044 * rho**2
                + 2650 * rho
                - 377
            )
            * r**3
            - 8
            * (
                9216 * rho**10
                - 232448 * rho**9
                - 1256448 * rho**8
                - 2428224 * rho**7
                - 1331808 * rho**6
                + 2314912 * rho**5
                + 4628128 * rho**4
                + 3289308 * rho**3
                + 941094 * rho**2
                + 17232 * rho
                - 8719
            )
            * r**2
            + 8
            * (
                24576 * rho**11
                - 64512 * rho**10
                - 970240 * rho**9
                - 3068160 * rho**8
                - 4269952 * rho**7
                - 1520992 * rho**6
                + 3514576 * rho**5
                + 5456064 * rho**4
                + 3202096 * rho**3
                + 682854 * rho**2
                - 49157 * rho
                - 5952
            )
            * r
            - (
                36864 * rho**11
                - 24576 * rho**10
                - 946176 * rho**9
                - 3121920 * rho**8
                - 4634496 * rho**7
                - 2574464 * rho**6
                + 1765216 * rho**5
                + 3792160 * rho**4
                + 2305160 * rho**3
                + 449576 * rho**2
                - 69262 * rho
                - 2607
            )
            * (2 * rho + 3)
        ) * eta**3

        # eta^2 term
        teta2 = (
            -2
            * (
                128 * (r - 2) * r**7
                + 96 * (32 * rho**3 + 32 * rho**2 - 2 * rho - 35) * r**6
                - 64 * (432 * rho**3 + 582 * rho**2 + 171 * rho - 277) * r**5
                - 8
                * (
                    2688 * rho**5
                    + 3552 * rho**4
                    - 9624 * rho**3
                    - 17988 * rho**2
                    - 7706 * rho
                    + 3411
                )
                * r**4
                + 16
                * (
                    3072 * rho**6
                    + 18816 * rho**5
                    + 29472 * rho**4
                    + 8568 * rho**3
                    - 16428 * rho**2
                    - 11974 * rho
                    - 71
                )
                * r**3
                + 2
                * (
                    4608 * rho**7
                    - 259072 * rho**6
                    - 904608 * rho**5
                    - 1174896 * rho**4
                    - 475504 * rho**3
                    + 231848 * rho**2
                    + 222510 * rho
                    + 11293
                )
                * r**2
                - 4
                * (
                    12288 * rho**8
                    - 48384 * rho**7
                    - 424480 * rho**6
                    - 945360 * rho**5
                    - 923200 * rho**4
                    - 270072 * rho**3
                    + 173214 * rho**2
                    + 121135 * rho
                    - 51
                )
                * r
                + (
                    4608 * rho**8
                    - 5568 * rho**7
                    - 145344 * rho**6
                    - 368496 * rho**5
                    - 401040 * rho**4
                    - 143908 * rho**3
                    + 56548 * rho**2
                    + 51463 * rho
                    - 1467
                )
                * (2 * rho + 3)
            )
            * eta**2
        )

        # eta^1 term
        teta1 = (
            2
            * (
                64 * (29 - 3 * r) * r**5
                + 16 * (84 * rho**2 - 325) * r**4
                - 96 * (80 * rho**3 + 228 * rho**2 + 120 * rho - 77) * r**3
                - 4
                * (144 * rho**4 - 19456 * rho**3 - 33144 * rho**2 - 16416 * rho + 4253)
                * r**2
                + 4
                * (
                    1920 * rho**5
                    - 9456 * rho**4
                    - 48992 * rho**3
                    - 58296 * rho**2
                    - 21104 * rho
                    + 6229
                )
                * r
                - (
                    288 * rho**5
                    + 336 * rho**4
                    - 15216 * rho**3
                    - 22680 * rho**2
                    - 11046 * rho
                    + 2953
                )
                * (2 * rho + 3)
            )
            * eta
        )

        # eta^0 term
        teta0 = -96 * (2 * r + 2 * rho - 15) * (2 * r - 2 * rho - 3) * (2 * r - 1)

        return teta6 + teta5 + teta4 + teta3 + teta2 + teta1 + teta0

    def c3(self, r: float) -> float:
        """
        ## c_3 component of the direct correlation function in spherical geometry
        - uses auxiliary functions __Q3, __B3, __S3, __T3

        :param r: radial distance from the center of the spherical wall
        :type r: float
        """

        rho = self.varrho  # abbreviation for the radius of the spherical wall
        eta = self.eta  # abbreviation for the packing fraction eta
        roots4a = (
            self.roots4a
        )  # abbreviation for the roots of the summation polynomial P_4^a

        if r < rho:
            return 0.0
        elif r < (rho + 0.5):
            term1 = 4 * r**2 + 4 * r + 8 * r * rho - 12 * rho**2 - 20 * rho - 11
            term1 = term1 * eta * (2 * r - 2 * rho + 1) ** 2 + 16 * (2 * r + 1)
            term2 = 2 * (4 * rho**2 + 6 * rho + 3) * eta * rho + 1
            term3 = (
                (2 * r + 1)
                * (2 * r - 1)
                * (rho**2 + rho + 1)
                / (8 * r * rho * (rho + 1))
            )
            term3 = term3 * np.log((2 * r + 1) / 2 / rho)

            nons = self.__S3(r) / (32 * r * term1 * term2**2 * (eta - 1) ** 2) - term3
            suma = 0.0 + 0.0j
            for t in roots4a:
                suma += self.__B3(t, r) * np.log(
                    (2 * r - 2 * t + 1) / (2 * rho - 2 * t)
                )
            return nons + suma.real
        elif r < (rho + 1.5):
            term1 = 4 * r * r - 4 * r + 8 * r * rho - 12 * rho * rho - 28 * rho - 11
            term2 = 2 * (4 * rho**2 + 6 * rho + 3) * eta * rho + 1
            term3 = (
                32.0
                * r
                * (term1 * eta * (2 * r - 2 * rho - 1) ** 2 + 16.0 * (2 * r - 1))
                * term2**2
                * (eta - 1) ** 3
            )
            term4 = (
                (2 * r + 1)
                * (2 * r - 1)
                * (rho**2 + rho + 1)
                / (8 * r * rho * (rho + 1))
            )
            term4 = term4 * np.log((2 * r - 1) / (2 * rho + 2.0))
            nons = self.__T3(r) / term3 + term4
            suma = 0.0 + 0.0j
            for t in roots4a:
                suma += self.__B3(t, r) * np.log(
                    (2 * r - 2 * t - 1) / (2 * rho - 2 * t + 2)
                )
            return nons - suma.real
        else:
            return -eta * (1.0 + eta + eta**2) / ((1.0 - eta) ** 3)

    def __BN1(self, t: complex, r: float) -> complex:
        """
        ## Summation part of cN1 in spherical geometry
        - auxiliary function, not to be called directly
        - not including the logarithmic part of the summand

        :param t: complex root of the summation polynomial P_4^b
        :type t: complex
        :param r: radial distance from the center of the spherical wall
        :type r: float
        """

        rho = self.varrho  # abbreviation for the radius of the spherical wall
        eta = self.eta  # abbreviation for the packing fraction eta

        # denominator common to all terms
        D = (2 * t**2 + 6 * (t - 1) * rho - 3) * eta * t + 1

        term1 = (4 * (4 * rho**2 + 6 * rho + 3) * rho * eta + 4) / D * t**3

        term2 = (
            (
                (8 * rho**2 + 8 * rho + 4) * eta * r**2
                - (8 * rho**4 + 40 * rho**3 + 66 * rho**2 + 46 * rho + 11) * eta
                + 12 * rho
            )
            / D
            * t**2
        )

        term3 = (
            (
                (2 * rho + 1) ** 2 * (2 * rho - 1) * rho * eta
                - 4 * (2 * rho + 1) * rho * eta * r**2
                - 4 * r**2
                + 28 * rho**2
                + 16 * rho
                + 7
            )
            / D
            * t
        )

        term4 = -(4 * r**2 - 20 * rho**2 - 16 * rho - 7) * rho / D

        return term1 + term2 + term3 + term4

    def cN1(self, r: float) -> float:
        """
        ## c_N1 component of the direct correlation function in spherical geometry
        - magnitude/radial part of the vector correlation function

        :param r: radial distance from the center of the spherical wall
        :type r: float
        """

        rho = self.varrho  # abbreviation for the radius of the spherical wall
        roots4b = (
            self.roots4b
        )  # abbreviation for the roots of the summation polynomial P_4^b

        if rho <= r < rho + 0.5:
            termA = -(
                (4 * r**2 - 2 * r - 2 * r * rho - 26 * rho**2 - 25 * rho - 11)
                * (2 * r - 2 * rho + 1)
            ) / (8 * r**2)
            sum_term = 0.0 + 0.0j
            for t in roots4b:
                logarg = (2 * r - 2 * t - 2 * rho + 1) / (-2 * t)
                sum_term += self.__BN1(t, r) * np.log(logarg)
            return termA - (3 / (8 * r**2)) * sum_term
        elif rho + 0.5 <= r < rho + 1.5:
            termB = (
                (4 * r**2 - 2 * r * rho - 26 * rho**2 - 27 * rho - 12)
                * (2 * r - 2 * rho - 3)
            ) / (8 * r**2)
            sum_term = 0.0 + 0.0j
            for t in roots4b:
                logarg = (2 * r - 2 * t - 2 * rho - 1) / (2 - 2 * t)
                sum_term += self.__BN1(t, r) * np.log(logarg)
            return termB + (3 / (8 * r**2)) * sum_term
        else:
            return 0.0

    def __BN2(self, t: complex, r: float) -> complex:
        """
        ## Summation part of cN2 in spherical geometry
        - auxiliary function, not to be called directly
        - not including the logarithmic part of the summand

        :param t: complex root of the summation polynomial P_4^b
        :type t: complex
        :param r: radial distance from the center of the spherical wall
        :type r: float
        """

        rho = self.varrho  # abbreviation for the radius of the spherical wall
        eta = self.eta  # abbreviation for the packing fraction eta

        # common denominators
        denA = (2 * t * t + 6 * (t - 1) * rho - 3) * eta * t + 1
        denB = 4 * denA

        # t^3 term
        tt3 = ((8 * rho**3 + 12 * rho**2 + 6 * rho - 3) * eta + 2) * t**3 / denA

        # t^2 term
        tt2 = (
            (
                4 * (rho**2 + rho + 1) * eta * r * r
                - (4 * rho**4 + 20 * rho**3 + 55 * rho**2 + 54 * rho + 17) * eta
                + 6 * rho
            )
            * t**2
            / denA
        )

        tt1 = (
            -(
                4 * (4 * rho**2 - 2 * rho - 1) * eta * r * r
                - (8 * rho**3 - 24 * rho**2 + 1) * (2 * rho + 1) * eta
                + 8 * r * r
                - 56 * rho**2
                - 32 * rho
                - 46
            )
            * t
            / denB
        )

        tt0 = (
            (
                (4 * r * r - 4 * rho**2 + 1) * eta * (2 * rho + 1)
                - 8 * r * r
                + 40 * rho**2
                + 32 * rho
                + 46
            )
            * rho
        ) / denB

        return tt3 + tt2 + tt1 + tt0

    def cN2(self, r: float) -> float:
        """
        ## c_N2 component of the direct correlation function in spherical geometry
        - magnitude/radial part of the vector correlation function

        :param r: radial distance from the center of the spherical wall
        :type r: float
        """

        rho = self.varrho  # abbreviation for the radius of the spherical wall
        eta = self.eta  # abbreviation for the packing fraction eta
        roots4b = self.roots4b  # abbreviation for the roots of the summ

        if rho <= r < rho + 0.5:
            termA = (
                -3
                * (2 * r + 2 * rho + 3)
                * (2 * r - 2 * rho + 1)
                * (2 * r + 1)
                * eta
                / (
                    (
                        (
                            4 * r * r
                            + 4 * r
                            + 8 * r * rho
                            - 12 * rho * rho
                            - 20 * rho
                            - 11
                        )
                        * eta
                        * (2 * r - 2 * rho + 1) ** 2
                        + 16 * (2 * r + 1)
                    )
                    * r
                )
            )
            termB = (
                -(4 * r * r - 2 * r * (rho + 1) - 26 * rho * rho - 25 * rho - 35)
                * (2 * r - 2 * rho + 1)
                / (8 * r * r)
            )

            # Summation term
            sum_val = 0.0 + 0.0j
            for t in roots4b:
                logarg = (2 * r - 2 * t - 2 * rho + 1) / (-2 * t)
                logterm = np.log(logarg)
                sum_val += self.__BN2(t, r) * logterm
            sum_term = -(3 / (4 * r * r)) * sum_val

            return termA + termB + sum_term

        elif rho + 0.5 <= r < rho + 1.5:
            termA = (
                (4 * r * r - 2 * r * (rho - 1) - 26 * rho * rho - 23 * rho - 29)
                * (2 * r - 2 * rho - 1)
                / (8 * r * r)
            )
            termB = -(
                (12 * r * r * (rho + 2) - 12 * rho**3 - 96 * rho * rho - 105 * rho - 94)
                * eta
                - 12 * r * r
                + 60 * rho * rho
                + 72 * rho
                + 85
            ) / (8 * (eta - 1) * r * r)
            termC = (
                -3
                * (
                    (
                        16 * r**4
                        + 16 * r**3 * rho
                        - 8 * (10 * rho * rho + 11 * rho + 6) * r * r
                        + 12 * (4 * rho * rho + 8 * rho + 3) * rho * r
                        + (6 * rho + 11) * (2 * rho + 1) ** 2
                    )
                    * eta
                    * (2 * r - 2 * rho - 1)
                    + 16 * (2 * r + 1) * (2 * r - 1)
                )
                / (
                    4
                    * (
                        (
                            4 * r * r
                            - 4 * r
                            + 8 * r * rho
                            - 12 * rho * rho
                            - 28 * rho
                            - 11
                        )
                        * eta
                        * (2 * r - 2 * rho - 1) ** 2
                        + 16 * (2 * r - 1)
                    )
                    * r
                    * r
                )
            )

            # Summation term
            sum_val = 0.0 + 0.0j
            for t in roots4b:
                logarg = (2 * r - 2 * t - 2 * rho - 1) / (2 - 2 * t)
                logterm = np.log(logarg)
                sum_val += self.__BN2(t, r) * logterm
            sum_term = (3 / (4 * r * r)) * sum_val

            return termA + termB + termC + sum_term
        else:
            return 0.0

    def chs(self, r: float) -> float:
        """
        ## Hard-sphere contribution to the direct correlation function in spherical geometry
        - sum of c0, c1, c2, c3, cN1 and cN2 components

        :param r: radial distance from the center of the spherical wall
        :type r: float
        :return: value of the hard-sphere contribution to the direct correlation function at distance r
        :rtype: float
        """
        return (
            self.c0(r)
            + self.c1(r)
            + self.c2(r)
            + self.c3(r)
            + self.cN1(r)
            + self.cN2(r)
        )

    def catt(self, r: float) -> float:
        """
        ## Attractive contribution to the direct correlation function in spherical geometry
        - based on the cut Lennard-Jones potential
        - LJ parameters given in the class initialization

        :param r: radial distance from the center of the spherical wall
        :type r: float
        :return: value of the attractive contribution to the direct correlation function at distance r
        :rtype: float
        """

        eta = self.eta  # abbreviation for the packing fraction eta
        varrho = self.varrho  # abbreviation for the radius of the spherical wall
        beta = self.beta  # abbreviation for the inverse temperature beta
        epsilon = self.epsilon  # abbreviation for the LJ energy parameter epsilon
        rc = self.rc  # abbreviation for the LJ cutoff radius rc

        if r >= varrho + 0.5 + rc:
            cat = 32 * (rc**3 - 1) / (rc**3)
        elif r >= varrho + 1.5:
            cat = 6 * (
                8 * (2 * rc**3 - 1) / (3 * rc**3)
                - r / rc**4
                + ((2 * varrho + 1) ** 2 - 8 * rc**2) / (4 * r * rc**4)
                + 4 * (6 * varrho - 2 * r + 3) / (3 * r * (2 * varrho - 2 * r + 1) ** 3)
            )
        elif r >= varrho:
            cat = 6 * (
                8 * (rc**3 - 1) / (3 * rc**3)
                + r * (rc**4 - 1) / rc**4
                + ((2 * varrho + 1) ** 2 - 8 * rc**2) / (4 * r * rc**4)
                - (4 * varrho**2 + 4 * varrho - 7) / (4 * r)
            )
        else:
            cat = 0.0

        return eta * beta * epsilon * cat

    def betamu(self, withLJ: bool = False) -> float:
        """
        ## Reduced bulk chemical potential beta*mu
        - includes attractive contribution if withLJ=True

        :param withLJ: whether to include the attractive contribution, defaults to False
        :type withLJ: bool, optional
        :return: value of the reduced bulk chemical potential
        :rtype: float
        """

        eta = self.eta  # abbreviation for the packing fraction eta
        mu = (
            (5 * eta**3 - 13 * eta**2 + 14 * eta) / 2 / (1 - eta) ** 3
            - np.log(1 - eta)
            + np.log(self.rhob)
        )
        if withLJ:
            return mu - self.catt(self.R + self.rc + 0.5 + self.varrho)
        else:
            return mu

    def Phi(self, r: float) -> float:
        """
        ## Reduced excess free energy density Phi in spherical geometry
        - based on the original Rosenfeld FMT formulation

        :param r: radial distance from the center of the spherical wall
        :type r: float
        :return: value of the reduced excess free energy density at distance r
        :rtype: float
        """
        x = 1 - self.n3(r)
        return (
            -self.n0(r) * np.log(x)
            + (self.n1(r) * self.n2(r) - self.N1(r) * self.N2(r)) / x
            + (self.n3(r) ** 3 - 3 * self.n2(r) * self.N2(r) ** 2) / (24 * np.pi * x**2)
        )

    def betaV(self, r: float, withLJ: bool = False) -> float:
        """
        ## Reduced external potential beta*V in spherical geometry
        - external potential inducing constant density profile
        - assuming attractive contribution if withLJ=True

        :param r: radial distance from the center of the spherical wall
        :type r: float
        :param withLJ: whether to include the attractive contribution, defaults to False
        :type withLJ: bool, optional
        :return: value of the reduced external potential at distance r
        :rtype: float
        """
        # closer than the particle radius: infinite potential
        if r < self.varrho + self.R:
            return np.inf

        # else potential given by imposing constant density profile
        V = self.betamu(withLJ) - np.log(self.rhob) + self.chs(r)
        if withLJ:
            return V + self.catt(r)
        else:
            return V
