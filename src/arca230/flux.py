import numpy as np


class PointSourceFlux:
    """
    Example of single power law point source flux with spectral index gamma and normalisation norm
    """

    def __init__(self, gamma, norm):
        self.gamma = gamma
        self.norm = norm

    def display_parameters(self):
        print(f"Gamma: {self.gamma}")
        print(f"Normalization: {self.norm}")

    def dNdE(self, loge):
        """
        Calculate the number of neutrinos per energy in GeV

        Parameters:
        - logE: log of the true neutrino energy in GeV

        Returns:
        - number of neutrinos per energy in GeV
        """
        return self.norm * np.power(10, -self.gamma * loge)
