from __future__ import print_function
from builtins import object
from numpy import pi, sqrt, log
import numpy as np


class TimeDelayLikelihood(object):
    """
    The likelihood of a strong lensing time-delay system as
    modelled in http://arxiv.org/pdf/1306.4732v2.pdf
    and http://arxiv.org/pdf/0910.2773v2.pdf
    and http://arxiv.org/pdf/1607.01790.pdf

    """

    def __init__(self, name, z_d, z_s, lambda_d, mu_d, sigma_d):
        super(TimeDelayLikelihood, self).__init__()
        self.name = name
        self.z_d = z_d
        self.z_s = z_s
        self.lambda_d = lambda_d
        self.mu_d = mu_d
        self.sigma_d = sigma_d

    @classmethod
    def load_catalog(cls, filename, lambdad):
        """
        Create a list of StrongLensLikelihood objects from a catalog
        with rows of form 
        z_d   z_s   Ddt_obs  sigma_d  ...
        """
        cat = np.loadtxt(filename)
        print("#############\t  Reading mock file with", len(cat[:, 0]), "lenses")
        data_cls = []
        for row in cat:
            z_d, z_s, Ddt_obs, sigma_d = row[:4]
            mu_d = np.log(np.abs(Ddt_obs - lambdad))  # lambda_d = 1000. fixed
            drow = cls("mock", z_d, z_s, lambdad, mu_d, sigma_d)
            data_cls.append(drow)
        return data_cls

    def likelihood(self, comovingDistance, omega_k, H0):
        """
        Evaluate the likelihood given a function
        comovingDistance(z) in Mpc, omega_k, and H0 in km/s/Mpc
        """
        # 0910.2773v2 equation 5
        # lambda_d = 1388.8
        # mu_d = 6.4682
        # sigma_d = 0.20560

        D_dt = self.D_deltat(comovingDistance, omega_k, H0)
        d = log(D_dt - self.lambda_d)
        # print D_dt, self.lambda_d
        return -0.5 * (d - self.mu_d)**2 / self.sigma_d**2 - d - log(sqrt(2 * pi) * self.sigma_d)

    def D_deltat(self, comovingDistance, omega_k, H0):
        """
        The nearest thing we have to a "datapoint" for this likelihood.
        A reduced, derived quantity

        Defined under equation 6 of http://arxiv.org/pdf/0910.2773v2.pdf
        """
        z_d = self.z_d
        z_s = self.z_s
        # D_d = angular diameter distance to lens
        # D_s = angular diameter distance to source
        # D_ds = angular diameter distance from lens to source
        c = 299792.4580  # km/s
        D_H = c / H0  # Mpc

        chi_s = comovingDistance(z_s)
        chi_d = comovingDistance(z_d)

        D_s = chi_s / (1 + z_s)
        D_d = chi_d / (1 + z_d)

        f_s = sqrt(1 + omega_k * chi_d**2 / D_H)
        f_d = sqrt(1 + omega_k * chi_s**2 / D_H)
        D_ds = (f_s * chi_s - f_d * chi_d) / (1 + z_s)

        return (1 + z_d) * D_d * D_s / D_ds


class RXJ1131(TimeDelayLikelihood):
    def __init__(self):
        z_d = 0.295
        z_s = 0.654
        lambda_d = 1388.8
        mu_d = 6.4682
        sigma_d = 0.20560
        super(RXJ1131, self).__init__(
            "RXJ1131", z_d, z_s, lambda_d, mu_d, sigma_d)


class B1608(TimeDelayLikelihood):
    def __init__(self):
        z_d = 0.6304
        z_s = 1.394
        lambda_d = 4000.0
        mu_d = 7.0531
        sigma_d = 0.22824
        super(B1608, self).__init__("B1608", z_d, z_s, lambda_d, mu_d, sigma_d)


class HE0435(TimeDelayLikelihood):
    def __init__(self):
        z_d = 0.4546
        z_s = 1.693
        lambda_d = 653.9
        mu_d = 7.5793
        sigma_d = 0.10312
        super(HE0435, self).__init__(
            "HE0435", z_d, z_s, lambda_d, mu_d, sigma_d)


if __name__ == '__main__':
    import astropy.cosmology
    z_d = 0.295
    z_s = 0.658
    lambda_d = 1388.8
    mu_d = 6.4682
    sigma_d = 0.20560
    RXJ = TimeDelayLikelihood("RXJ1131", z_d, z_s, lambda_d, mu_d, sigma_d)

    z_d = 0.6304
    z_s = 1.394
    lambda_d = 4000.0
    mu_d = 7.053
    sigma_d = 0.2282
    B1608 = TimeDelayLikelihood("B1608", z_d, z_s, lambda_d, mu_d, sigma_d)

    z_d = 0.4546
    z_s = 1.693
    lambda_d = 653.9
    mu_d = 7.5793
    sigma_d = 0.10312
    HE0435 = TimeDelayLikelihood("HE0435", z_d, z_s, lambda_d, mu_d, sigma_d)

    H0_values = np.arange(40.0, 100.0, 0.1)
    rxj_likes = np.zeros_like(H0_values)
    b1608_likes = np.zeros_like(H0_values)
    he0435_likes = np.zeros_like(H0_values)

    for i, H0 in enumerate(H0_values):
        omega_k = 0.0
        cosmo = astropy.cosmology.FlatLambdaCDM(H0=H0, Om0=0.3)

        def comovingDistance(z): return (
            cosmo.comoving_distance(z) / astropy.units.megaparsec).value
        rxj_likes[i] = RXJ.likelihood(comovingDistance, omega_k, H0)
        b1608_likes[i] = B1608.likelihood(comovingDistance, omega_k, H0)
        he0435_likes[i] = HE0435.likelihood(comovingDistance, omega_k, H0)

    import pylab
    rxj_likes -= np.nanmax(rxj_likes)
    b1608_likes -= np.nanmax(b1608_likes)
    he0435_likes -= np.nanmax(he0435_likes)

    tdsl_likes = rxj_likes + b1608_likes + he0435_likes
    tdsl_likes -= np.nanmax(tdsl_likes)

    pylab.plot(H0_values, np.exp(rxj_likes))
    pylab.plot(H0_values, np.exp(b1608_likes))
    pylab.plot(H0_values, np.exp(he0435_likes))
    pylab.plot(H0_values, np.exp(tdsl_likes), "black")
    pylab.show()
