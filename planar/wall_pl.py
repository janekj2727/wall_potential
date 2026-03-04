import numpy as np


class wall_pl:
    r"""
    ## Planar wall class
        - wall inducing constant density profile in planar geometry
        - using Rosenfeld's Fundamental Measure Theory (FMT) for hard-sphere contribution
        - optionally including Lennard-Jones (LJ) attractive interactions via mean-field approximation
    """

    def __init__(
        self, eta: float, rc: float = 2.5, beta: float = 1.0, epsilon: float = 1.0
    ):
        r"""
        ## Planar wall initialization
        - wall inducing constant density profile in planar geometry
        - using Rosenfeld's Fundamental Measure Theory (FMT) for hard-sphere contribution
        - optionally including Lennard-Jones (LJ) attractive interactions via mean-field approximation
        - setting up parameters for hard-sphere fluid and for LJ interactions
        - setting up polynomial and its roots for direct correlation function calculation

        :param eta: packing fraction of the hard-sphere fluid
        :type eta: float
        :param rc: cut-off radius for the LJ potential, defaults to 2.5
        :type rc: float, optional
        :param beta: inverse temperature, defaults to 1
        :type beta: float, optional
        :param epsilon: energy parameter for LJ potential, defaults to 1
        :type epsilon: float, optional
        """
        self.R = 0.5
        self.sigma = 1
        self.eta = eta
        self.rhob = 6 * eta / np.pi
        self.beta = beta
        self.epsilon = epsilon
        self.rc = rc
        coeffs3 = (1.0, 0.0, -3.0 * self.eta, 2.0 * self.eta)
        polyn3 = np.polynomial.Polynomial(coeffs3)
        self.roots3 = polyn3.roots()

    def n3(self, z: float) -> float:
        """
        ## Weighted density n_3 in planar geometry
        - scalar weighted density

        :param z: distance from the wall
        :type z: float
        :return: weighted density n_3 at distance z
        :rtype: float
        """
        if (z >= 0) and (z <= 2 * self.R):
            return np.pi * self.rhob * (self.R * z**2 - z**3 / 3)
        elif z > 2 * self.R:
            return 4 / 3 * np.pi * self.R**3 * self.rhob
        else:
            return 0

    def n2(self, z: float) -> float:
        r"""
        ## Weighted density n_2 in planar geometry
        - scalar weighted density

        :param z: distance from the wall
        :type z: float
        :return: weighted density n_2 at distance z
        :rtype: float
        """
        if (z >= 0) and (z <= 2 * self.R):
            return 2 * np.pi * self.R * self.rhob * z
        elif z > 2 * self.R:
            return 4 * np.pi * self.R**2 * self.rhob
        else:
            return 0

    def n1(self, z: float) -> float:
        r"""
        ## Weighted density n_1 in planar geometry
        - scalar weighted density

        :param z: distance from the wall
        :type z: float
        :return: weighted density n_1 at distance z
        :rtype: float
        """
        if (z >= 0) and (z <= 2 * self.R):
            return self.rhob * z / 2
        elif z > 2 * self.R:
            return self.rhob * self.R
        else:
            return 0

    def n0(self, z: float) -> float:
        r"""
        ## Weighted density n_0 in planar geometry
        - scalar weighted density

        :param z: distance from the wall
        :type z: float
        :return: weighted density n_0 at distance z
        :rtype: float
        """
        if (z >= 0) and (z <= 2 * self.R):
            return self.rhob * z / (2 * self.R)
        elif z > 2 * self.R:
            return self.rhob
        else:
            return 0

    def N2(self, z: float) -> float:
        r"""
        ## Weighted density N_2 in planar geometry
        - vector weighted density (component normal to wall)

        :param z: distance from the wall
        :type z: float
        :return: weighted density N_2 at distance z
        :rtype: float
        """
        if (z >= 0) and (z <= 2 * self.R):
            return 2 * np.pi * self.rhob * (z**2 / 2 - z * self.R)
        else:
            return 0

    def N1(self, z: float) -> float:
        r"""
        ## Weighted density N_1 in planar geometry
        - vector weighted density (component normal to wall)

        :param z: distance from the wall
        :type z: float
        :return: weighted density N_1 at distance z
        :rtype: float
        """
        if (z >= 0) and (z <= 2 * self.R):
            return self.rhob / (2 * self.R) * (z**2 / 2 - self.R * z)
        else:
            return 0

    def dphidn0(self, z: float) -> float:
        """
        ## Derivative of Phi with respect to n_0 in planar geometry
        - Phi: the reduced excess free energy density
        - n_0: scalar weighted density n_0

        :param z: distance from the wall
        :type z: float
        :return: derivative of Phi with respect to n_0 at distance z
        :rtype: float
        """
        if z < 0:
            return 0
        return -np.log(1 - self.n3(z))

    def dphidn1(self, z: float) -> float:
        """
        ## Derivative of Phi with respect to n_1 in planar geometry
        - Phi: the reduced excess free energy density
        - n_1: scalar weighted density n_1

        :param z: distance from the wall
        :type z: float
        :return: derivative of Phi with respect to n_1 at distance z
        :rtype: float
        """
        if z < 0:
            return 0
        return self.n2(z) / (1 - self.n3(z))

    def dphidn2(self, z: float) -> float:
        """
        ## Derivative of Phi with respect to n_2 in planar geometry
        - Phi: the reduced excess free energy density
        - n_2: scalar weighted density n_2

        :param z: distance from the wall
        :type z: float
        :return: derivative of Phi with respect to n_2 at distance z
        :rtype: float
        """
        if z < 0:
            return 0
        x = 1 - self.n3(z)
        return self.n1(z) / x + (self.n2(z) ** 2 - self.N2(z) ** 2) / (8 * np.pi * x**2)

    def dphidn3(self, z: float) -> float:
        """
        ## Derivative of Phi with respect to n_3 in planar geometry
        - Phi: the reduced excess free energy density
        - n_3: scalar weighted density n_3

        :param z: distance from the wall
        :type z: float
        :return: derivative of Phi with respect to n_3 at distance z
        :rtype: float
        """
        if z < 0:
            return 0
        x = 1 - self.n3(z)
        return (
            self.n0(z) / x
            + (self.n1(z) * self.n2(z) - self.N1(z) * self.N2(z)) / (x**2)
            + (self.n2(z) ** 3 - 3 * self.n2(z) * self.N2(z) ** 2) / (12 * np.pi * x**3)
        )

    def dphidN1(self, z: float) -> float:
        """
        ## Derivative of Phi with respect to N_1 in planar geometry
        - Phi: the reduced excess free energy density
        - N_1: vector weighted density N_1 (magnitude/normal component)

        :param z: distance from the wall
        :type z: float
        :return: derivative of Phi with respect to N_1 at distance z
        :rtype: float
        """
        if z < 0:
            return 0
        return -self.N2(z) / (1 - self.n3(z))

    def dphidN2(self, z: float) -> float:
        """
        ## Derivative of Phi with respect to N_2 in planar geometry
        - Phi: the reduced excess free energy density
        - N_2: vector weighted density N_2 (magnitude/normal component)

        :param z: distance from the wall
        :type z: float
        :return: derivative of Phi with respect to N_2 at distance z
        :rtype: float
        """
        if z < 0:
            return 0
        x = 1 - self.n3(z)
        return -self.N1(z) / x - self.n2(z) * self.N2(z) / (4 * np.pi * x**2)

    def c0(self, z: float) -> float:
        """
        ## c_0 component of the direct correlation function in planar geometry
        - uses polynomial roots precomputed in __init__

        :param z: distance from the wall
        :type z: float
        :return: c_0 component of the direct correlation function at distance z
        :rtype: float
        """
        if z >= 1.5:
            return np.log(1 - self.eta)
        elif z >= 0.5:
            suma = 0.0
            for t in self.roots3:
                suma += (
                    (np.log(1 - t) - np.log(2 * z - 2 * t - 1))
                    * (self.eta * t**2 - 1)
                    / (2 * t * (t - 1) * self.eta)
                )
            return -(
                suma
                + (z - 0.5) * np.log((z - 2.0) * (2 * z - 1) ** 2 * self.eta + 2)
                - (z + 0.5) * np.log(1 - self.eta)
                - (z - 2) * np.log(2.0)
                + 4.5
                - 3 * z
            )
        elif z >= 0:
            suma = 0.0
            for t in self.roots3:
                suma += (
                    (np.log(2 * z - 2 * t + 1) - np.log(-t))
                    * (self.eta * t**2 - 1)
                    / (2 * t * (t - 1) * self.eta)
                )
            return -(
                suma
                - (z + 0.5) * np.log((z - 1) * (2 * z + 1) ** 2 * self.eta + 2)
                + (z - 1) * np.log(2)
                + 1.5
                + 3 * z
            )
        else:
            return 0

    def c1(self, z: float) -> float:
        """
        ## c_1 component of the direct correlation function in planar geometry

        :param z: distance from the wall
        :type z: float
        :return: c_1 component of the direct correlation function at distance z
        :rtype: float
        """
        if z >= 1.5:
            return -3 * self.eta / (1 - self.eta)
        elif z >= 0.5:
            suma = 0.0
            for t in self.roots3:
                suma += (np.log(1 - t) - np.log(2 * z - 2 * t - 1)) / (2 * (t - 1))
            return -(suma - 3 * (2 * z - 1) * self.eta / (2 * (self.eta - 1)))
        elif z >= 0:
            suma = 0.0
            for t in self.roots3:
                suma += (np.log(2 * z - 2 * t + 1) - np.log(-t)) / (2 * (t - 1))
            return -suma
        else:
            return 0

    def c2(self, z: float) -> float:
        """
        ## c_2 component of the direct correlation function in planar geometry

        :param z: distance from the wall
        :type z: float
        :return: c_2 component of the direct correlation function at distance z
        :rtype: float
        """
        if z >= 1.5:
            return -3 * self.eta * (self.eta + 2) / (2 * (self.eta - 1) ** 2)
        elif z >= 0.5:
            suma = 0.0
            for t in self.roots3:
                suma += ((1 + 2 * t) * (np.log(2 * z - 2 * t - 1) - np.log(1 - t))) / (
                    8 * t * (t - 1) * (self.eta - 1)
                )
            return -(
                suma
                + 3
                * self.eta
                * (2 * z - 1)
                * (
                    (4 * z**2 - 8 * z + 1) * (2 * z - 1) * self.eta**2
                    + (16 * z**3 - 52 * z**2 + 46 * z - 10) * self.eta
                    - 2 * z
                    + 11
                )
                / (
                    8
                    * (self.eta - 1) ** 2
                    * (2 + (z - 2) * (2 * z - 1) ** 2 * self.eta)
                )
            )
        elif z >= 0:
            suma = 0.0
            for t in self.roots3:
                suma += ((1 + 2 * t) * (np.log(2 * z - 2 * t + 1) - np.log(-t))) / (
                    8 * t * (1 - t) * (self.eta - 1)
                )
            return -(
                suma
                - (3 * self.eta * (2 * self.eta * z + self.eta + 1) * (4 * z**2 - 1))
                / (8 * (self.eta - 1) * (2 + (z - 1) * (2 * z + 1) ** 2 * self.eta))
            )
        else:
            return 0

    def c3(self, z: float) -> float:
        """
        ## c_3 component of the direct correlation function in planar geometry

        :param z: distance from the wall
        :type z: float
        :return: c_3 component of the direct correlation function at distance z
        :rtype: float
        """
        if z >= 1.5:
            return -self.eta * (self.eta**2 + self.eta + 1) / (1 - self.eta) ** 3
        elif z >= 0.5:
            suma = 0.0
            aux = (self.eta - 1) ** 3 * (2 + (z - 2) * (2 * z - 1) ** 2 * self.eta)
            for t in self.roots3:
                suma += (
                    (np.log(2 * z - 2 * t - 1) - np.log(1 - t))
                    / (4 * t * (t - 1) * self.eta * (self.eta - 1) ** 2)
                    * (
                        (
                            4 / 3 * z**2
                            + (12 * t**2 - 20 * t + 6) / 3 * z
                            - (8 * t + 1) / 3
                        )
                        * self.eta**3
                        + (
                            (18 * t - 11) / 6 * z**2
                            - (24 * t**2 - 26 * t + 15) / 3 * z
                            + (142 * t + 59) / 24
                        )
                        * self.eta**2
                        + (2 * z**2 + (4 * t**2 - 8 * t) * z - (t + 3)) * self.eta
                        + 2
                    )
                )
            return -(
                suma
                + (
                    32 * z**6
                    - 160 * z**5
                    + 224 * z**4
                    + 40 * z**3
                    - 294 * z**2
                    + 192 * z
                    - 37
                )
                / 4
                / aux
                * self.eta**4
                + (
                    128 * z**6
                    - 848 * z**5
                    + 2656 * z**4
                    - 4840 * z**3
                    + 4872 * z**2
                    - 2433 * z
                    + 558
                )
                / 16
                / aux
                * self.eta**3
                + (
                    128 * z**6
                    - 816 * z**5
                    + 1440 * z**4
                    + 104 * z**3
                    - 2256 * z**2
                    + 1929 * z
                    - 720
                )
                / 16
                / aux
                * self.eta**2
                + (96 * z**4 - 400 * z**3 + 540 * z**2 - 396 * z + 227)
                / 8
                / aux
                * self.eta
                + (6 * z - 9) / aux
                - 3 * z * np.log(2)
            )
        elif z >= 0:
            suma = 0.0
            aux = (
                4
                * (self.eta - 1) ** 2
                * (2 + (z - 2) * (2 * z - 1) ** 2 * self.eta)
                * (2 + (4 * z**3 - 3 * z - 1) * self.eta)
            )
            aux2 = (self.eta - 1) ** 2 * (2 + (z - 1) * (2 * z + 1) ** 2 * self.eta)
            for t in self.roots3:
                suma += (
                    (np.log(-t) - np.log(2 * z - 2 * t + 1))
                    / (t * self.eta * (self.eta - 1) ** 2 * (t - 1))
                    * (
                        0.5
                        + (
                            z**2 / 3
                            + (6 * t**2 - 10 * t + 3) / 6 * z
                            - (8 * t + 1) / 12
                        )
                        * self.eta**3
                        + (
                            (18 * t - 11) / 24 * z**2
                            - (24 * t**2 - 26 * t + 15) / 12 * z
                            + (142 * t + 59) / 96
                        )
                        * self.eta**2
                        + (z**2 / 2 + (t**2 - 2 * t) * z - (t + 3) / 4) * self.eta
                    )
                )
            return -(
                suma
                + (16 * z**5 - 16 * z**3 - 4 * z**2 + 3 * z + 1)
                / 4
                / aux2
                * self.eta**4
                - (80 * z**5 - 96 * z**4 - 128 * z**3 + 28 * z**2 + 75 * z + 23)
                / 8
                / aux2
                * self.eta**3
                - (48 * z**5 + 240 * z**4 + 72 * z**3 - 96 * z**2 - 261 * z - 111)
                / 16
                / aux2
                * self.eta**2
                + (96 * z**4 + 48 * z**3 - 84 * z**2 - 132 * z - 45)
                / 8
                / aux2
                * self.eta
                + (6 * z + 3) / aux2
                + 3 * z * np.log(2)
            )
        else:
            return 0

    def cN1(self, z: float) -> float:
        """
        ## c_N1 component of the direct correlation function in planar geometry
        - component normal to wall

        :param z: distance from the wall
        :type z: float
        :return: c_N1 component of the direct correlation function at distance z
        :rtype: float
        """
        if z >= 1.5:
            return 0.0
        elif z >= 0.5:
            suma = 0.0
            for t in self.roots3:
                suma += (
                    (np.log(1 - t) - np.log(2 * z - 1 - 2 * t))
                    * (1 + self.eta * (2 * z - 1) * t**2 - 2 * self.eta * t * z)
                    / (2 * t * (t - 1) * self.eta)
                )
            return -(suma + 3 * z * (1 + np.log(2)) - 1.5 * (3 + np.log(2)))
        elif z >= 0:
            suma = 0.0
            for t in self.roots3:
                suma += (
                    (np.log(2 * z - 2 * t + 1) - np.log(-t))
                    * (1 + self.eta * (2 * z - 1) * t**2 - 2 * t * self.eta * z)
                    / (2 * self.eta * t * (t - 1))
                )
            return -(suma - 3 * z * (1 + np.log(2)) - 1.5 * (1 - np.log(2)))
        else:
            return 0.0

    def cN2(self, z: float) -> float:
        """
        ## c_N2 component of the direct correlation function in planar geometry
        - component normal to wall

        :param z: distance from the wall
        :type z: float
        :return: c_N2 component of the direct correlation function at distance z
        :rtype: float
        """
        if z >= 1.5:
            return 0.0
        elif z >= 0.5:
            suma = 0.0
            for t in self.roots3:
                suma += (
                    (1 + ((2 * z - 1) * t**2 - 2 * (z + 1) * t + z) * self.eta)
                    * (np.log(1 - t) - np.log(2 * z - 2 * t - 1))
                    / (2 * t * (t - 1) * self.eta)
                )
            return -(
                suma
                + 3
                * (
                    (8 * z**3 - 22 * z**2 + 15 * z - 3) * self.eta**2
                    - (4 * z**3 - 12 * z**2 + 9 * z - 5) * self.eta
                    - 2
                )
                * (2 * z - 3)
                / (2 * (2 + (z - 2) * (2 * z - 1) ** 2 * self.eta) * (self.eta - 1))
                + 3 * (z - 0.5) * np.log(2)
            )
        elif z >= 0:
            suma = 0.0
            for t in self.roots3:
                suma += (
                    (1 + ((2 * z - 1) * t**2 - 2 * (z + 1) * t + z) * self.eta)
                    * (np.log(2 * z - 2 * t + 1) - np.log(-t))
                    / (2 * t * self.eta * (t - 1))
                )
            return -(
                suma
                - 3 * (z - 1 / 2) * np.log(2)
                - (3 * (2 * z + 1) * (2 + (4 * z**3 - 3 * z - 2) * self.eta))
                / (4 + (8 * z**3 - 6 * z - 2) * self.eta)
            )
        else:
            return 0.0

    def chs(self, z: float) -> float:
        """
        ## Hard-sphere contribution to the one-body direct correlation function in planar geometry

        :param z: distance from the wall
        :type z: float
        :return: local hard-sphere contribution at distance z
        :rtype: float
        """
        aux1 = (self.eta - 1) ** 3 * (2 + (2 * z - 1) ** 2 * (z - 2) * self.eta)
        aux2 = (self.eta - 1) ** 2 * (2 + (2 * z + 1) ** 2 * (z - 1) * self.eta)
        if z >= 1.5:
            return -(
                (5 * self.eta**3 - 13 * self.eta**2 + 14 * self.eta)
                / (2 * (1 - self.eta) ** 3)
                - np.log(1 - self.eta)
            )
        elif z >= 0.5:
            suma = 0.0
            for t in self.roots3:
                suma += (
                    (np.log(1 - t) - np.log(2 * z - 2 * t - 1))
                    * (2 * z - 1)
                    / (96 * (t - 1) * t * (self.eta - 1) ** 2)
                    * (
                        (48 * t**2 - 16 * t - 16 * z - 8) * self.eta**2
                        - (96 * t**2 + 36 * t * z - 70 * t - 22 * z - 23) * self.eta
                        + 48 * t**2
                        - 24 * z
                        + 12
                    )
                )
            return -(
                suma.real
                + (
                    64 * z**6
                    - 320 * z**5
                    + 592 * z**4
                    - 544 * z**3
                    + 324 * z**2
                    - 108 * z
                    + 13
                )
                / 8
                / aux1
                * self.eta**4
                + (
                    128 * z**6
                    - 848 * z**5
                    + 2176 * z**4
                    - 2488 * z**3
                    + 1032 * z**2
                    - 141 * z
                    - 36
                )
                / 16
                / aux1
                * self.eta**3
                + (
                    128 * z**6
                    - 816 * z**5
                    + 1824 * z**4
                    - 1864 * z**3
                    + 1032 * z**2
                    - 123 * z
                    + 90
                )
                / 16
                / aux1
                * self.eta**2
                + (4 * z**3 - 12 * z**2 + 3 * z - 5) / aux1 * self.eta
                - (z + 0.5) * np.log(1 - self.eta)
                + (z - 1 / 2)
                * (np.log(2 + (2 * z - 1) ** 2 * (z - 2) * self.eta) + 2 * np.log(2))
            )
        elif z >= 0:
            suma = 0.0
            for t in self.roots3:
                suma += (
                    (np.log(-t) - np.log(2 * z - 2 * t + 1))
                    * (2 * z - 1)
                    / (96 * t * (t - 1) * (self.eta - 1) ** 2)
                    * (
                        (-48 * t**2 + 16 * t + 16 * z + 8) * self.eta**2
                        + (36 * t * z - 22 * z + 96 * t**2 - 70 * t - 23) * self.eta
                        - 48 * t**2
                        + 24 * z
                        - 12
                    )
                )
            return -(
                suma.real
                + 2
                * self.eta
                * (2 * z + 1)
                / aux2
                * (
                    (z**4 - z**3 / 2 - 3 * z**2 / 4 + z / 8 + 1 / 8) * self.eta**3
                    - (5 * z**4 / 2 - 5 * z**3 / 4 - 9 * z**2 / 8 - 7 * z / 16 - 1 / 4)
                    * self.eta**2
                    - (
                        3 * z**4 / 4
                        - 21 * z**3 / 8
                        - 21 * z**2 / 16
                        + 117 * z / 32
                        + 33 / 32
                    )
                    * self.eta
                    + 3 / 2
                )
                - (z + 0.5) * np.log(2 + (2 * z + 1) ** 2 * (z - 1) * self.eta)
                - 2 * (z - 1) * np.log(2)
            )
        else:
            return 0.0
        
    def chsPade(self, z: float) -> float:
        # Powers of eta
        eta = self.eta
        e2 = eta * eta
        e3 = e2 * eta
        e4 = e2 * e2
        e5 = e2 * e3
        e6 = e3 * e3

        # Denominators
        den1 = 2.0 * (1.0 - eta)**3
        den2 = 3.0 * (1.0 - eta) * (eta + 2.0)
        den3 = e2 - 2.0 * eta + 4.0

        # Auxiliary quantities
        aux1 = eta * (5.0 * e2 - 13.0 * eta + 14.0) / den1
        aux2 = (-25.0 * e4 + 4.0 * e3 - 9.0 * e2
                - 20.0 * eta - 4.0)
        aux3 = (107.0 * e6 - 345.0 * e5 + 501.0 * e4
                - 253.0 * e3 - 114.0 * e2
                + 660.0 * eta - 232.0)

        ln1me = np.log(1.0 - eta)

        # Padé coefficients
        B1 = (13.0 * e2 + 7.0 * eta - 2.0) / den2
        A0 = aux1 - ln1me
        A1 = -B1 * ln1me + aux1 * B1
        A2 = (aux2 * ln1me + aux3 * eta / den1) / den2 / den3
        B2 = -aux2 / den2 / den3

        zarg = z - 1.5

        # Piecewise definition
        if z > 1.5:
            return -(
                (5 * self.eta**3 - 13 * self.eta**2 + 14 * self.eta)
                / (2 * (1 - self.eta) ** 3)
                - np.log(1 - self.eta)
            ) # bulk limit - simple expression
        elif z > 0.5:
            return -((A0 + A1 * zarg + A2 * zarg**2)
                / (1.0 + B1 * zarg + B2 * zarg**2))
        elif z >= 0.0:
            return -((A0 + A1 * zarg + A2 * zarg**2)
                / (1.0 + B1 * zarg + B2 * zarg**2)) # irrelevant for Vext
        else:
            return 0.0


    def catt(self, z: float) -> float:
        """
        ## Mean-field Lennard-Jones attractive contribution in planar geometry

        :param z: distance from the wall
        :type z: float
        :return: attractive tail contribution at distance z
        :rtype: float
        """
        if z >= self.R + self.rc:
            cat = 8 * self.R**3 / (3 * self.rc**3) - 1.0 / 3.0
        elif z >= 3 * self.R:
            cat = (
                (4 * self.R**3) / (3 * self.rc**3)
                - 1.0 / 3.0
                - self.R**3 * (self.R - z) / (self.rc**4)
                - self.R**3 / (3 * (self.R - z) ** 3)
            )
        elif z >= 0:
            cat = (
                (4 * self.R**3) / (3 * self.rc**3)
                - 5.0 / 48.0
                - self.R**3 * (self.R - z) / (self.rc**4)
                - z / 16.0 / self.R
            )
        else:
            cat = 0.0
        return -128.0 * np.pi * self.rhob * self.beta * self.epsilon * self.R**3 * cat

    def betamu(self, withLJ: bool = False) -> float:
        """
        ## Bulk reduced chemical potential beta*mu
        - includes attractive contribution if withLJ is True

        :param withLJ: include mean-field LJ tail in chemical potential
        :type withLJ: bool
        :return: reduced chemical potential
        :rtype: float
        """
        mu = (
            (5 * self.eta**3 - 13 * self.eta**2 + 14 * self.eta)
            / 2
            / (1 - self.eta) ** 3
            - np.log(1 - self.eta)
            + np.log(self.rhob)
        )
        if withLJ:
            return mu - self.catt(self.R + self.rc + 0.5)
        else:
            return mu

    def Phi(self, z: float) -> float:
        """
        ## Reduced excess free energy density Phi in planar geometry
        - based on the original Rosenfeld FMT formulation

        :param z: distance from the wall
        :type z: float
        :return: local reduced excess free energy density at distance z
        :rtype: float
        """
        x = 1 - self.n3(z)
        return (
            -self.n0(z) * np.log(x)
            + (self.n1(z) * self.n2(z) - self.N1(z) * self.N2(z)) / x
            + (self.n3(z) ** 3 - 3 * self.n2(z) * self.N2(z) ** 2) / (24 * np.pi * x**2)
        )

    def betaV(self, z: float, withLJ: bool = False) -> float:
        """
        ## External potential beta * V(z) for the planar wall
        - external potential inducing constant density profile
        - assuming attractive contribution if withLJ=True

        :param z: distance from the wall
        :type z: float
        :param withLJ: whether to include LJ attractive tail
        :type withLJ: bool
        :return: external potential (in units of beta) at distance z
        :rtype: float
        """
        V = self.betamu(withLJ) - np.log(self.rhob) + self.chs(z)
        if withLJ:
            return V + self.catt(z)
        else:
            return V
        
    def betaVPade(self, z: float, withLJ: bool = False) -> float:
        """
        ## External potential beta * V(z) for the planar wall
        - external potential inducing constant density profile
        - assuming attractive contribution if withLJ=True
        - approximating c_HS by Padé approximant

        :param z: distance from the wall
        :type z: float
        :param withLJ: whether to include LJ attractive tail
        :type withLJ: bool
        :return: external potential (in units of beta) at distance z
        :rtype: float
        """
        V = self.betamu(withLJ) - np.log(self.rhob) + self.chsPade(z)
        if withLJ:
            return V + self.catt(z)
        else:
            return V
