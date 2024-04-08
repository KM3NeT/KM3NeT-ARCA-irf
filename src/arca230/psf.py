import pandas as pd
import numpy as np
from scipy.interpolate import interp1d


class PointSpreadFunction:
    """
    Loads the point spread function of the ARCA230 detector for numu selected as track or nue selected as shower. The
    point spread function is stored as the event density dP/dOmega per distance to the source psi. This function is interpolated
    for different true neutrino energy ranges.
    """

    def __init__(self, file_path="../data/psf_numuCC_track.csv"):
        self.file_path = file_path
        self.psf_data = None
        self.load_psf_data()

        self.interpolations = {}
        self.interpolate_psf()

    def load_psf_data(self):
        """
        Loads the data for the energy response
        """
        try:
            self.psf_data = pd.read_csv(self.file_path, delimiter=",")
            print("Point Spread Function data loaded successfully.")
        except FileNotFoundError:
            print(f"Error: File '{self.file_path}' not found.")
        except Exception as e:
            print(f"An error occurred while loading the data: {e}")

    def psf_response(self, logE):
        """
        Filters the point spread function data for the given true neutrino energy
        and returns the distribution of event densitities as a function of log10(psi [degrees])

        Parameters:
        - logE: log of the true neutrino energy in GeV

        Returns:
        - Distribution of event densitities as a function of log10(psi [degrees])
        """
        logE_mask = (self.psf_data["log10(nu_E [GeV]) low"] <= logE) & (self.psf_data["log10(nu_E [GeV]) high"] > logE)

        filtered_rows = self.psf_data[logE_mask].copy()
        return filtered_rows

    def interpolate_psf(self):
        """
        Interpolate the point spread function data and store the results for each unique energy range.
        """
        for _, row in self.psf_data.iterrows():
            logE_low = row["log10(nu_E [GeV]) low"]
            logE_high = row["log10(nu_E [GeV]) high"]

            # Check if the pair is already encountered
            if (logE_low, logE_high) in self.interpolations:
                continue

            filtered_df = self.psf_data[
                (self.psf_data["log10(nu_E [GeV]) low"] == logE_low)
                & (self.psf_data["log10(nu_E [GeV]) high"] == logE_high)
            ]

            self.interpolations[(logE_low, logE_high)] = interp1d(
                x=filtered_df["log10(psi [degrees])"], y=filtered_df["dP/dOmega"], kind="linear"
            )

    def get_interpolation_for_energy(self, logE):
        """
        Get the interpolation function for a specified energy within the available ranges.

        Parameters:
        - logE: the logarithm of the true neutrino energy

        Returns:
        - Interpolation of dP/dOmega versus log(alpha) in degrees
        """

        for (logE_low, logE_high), interp_func in self.interpolations.items():
            if logE_low <= logE < logE_high:
                return interp_func

        raise ValueError(f"No interpolation found for energy {logE}.")

    def integrate_sphere(self, logE, min_angle=1e-4, max_angle=179.99, nsteps=500):
        """
        Numerically integrate the point spread function over the sphere at energy logE

        Parameters:
        - logE: the logarithm of the true neutrino energy in GeV
        - min_angle: in degrees, low edge of the integration, should be larger than 0
        - max_angle: cone size of the integration in degrees
        - nsteps: number of steps in the integration

        Returns:
        - Integration of dP/dOmega over the specified cone
        """
        log_min_angle = 1e-4 if min_angle < 1e-4 else np.log10(min_angle)
        log_max_angle = np.log10(max_angle)

        dla = (log_max_angle - log_min_angle) / nsteps
        c = 0

        for la in np.linspace(log_min_angle, log_max_angle, nsteps + 1):
            P = self.eval(logE, la)
            y = self.d_omega_d_loga(la) * P * dla
            c += y

        return c

    def d_omega_d_loga(self, loga):
        """
        Jacobian d Omega / d log(a), with a angle in degrees

        Parameters:
        - loga: logarithm of the angle with the source in degrees

        Returns:
        - Jacobian d Omega / d log(a)
        """
        a = np.power(10, loga) * np.pi / 180
        if a > np.pi:
            return 0
            # we need to take care about non physical angles > 180 deg
        return np.sin(a) * 2 * np.pi * np.log(10) * a

    def eval(self, logE, loga):
        """
        The central function that returns the interpolated values of dP/dOmega

        Parameters:
        - logE: logarithm of the true neutrino energy in GeV
        - loga: logarithm of the angle with the source in degrees

        Returns:
        - dP/dOmega
        """
        interp_func = self.get_interpolation_for_energy(logE)

        return interp_func(loga)

    def fraction_below_angle(self, logE, angle_max):
        """
        Calculate the fraction of events below a given angle_max in degrees

        Parameters:
        - logE: logarithm of the true neutrino energy in GeV
        - angle_max: cone_size in degrees

        Returns:
        - Fraction of events within angle_max
        """
        c = self.integrate_sphere(logE)
        if c == 0:
            return 0
        else:
            return self.integrate_sphere(logE, 1e-4, angle_max) / c

    def event_table_within_cone(self, event_rate_table, angle_max):
        """
        Copies an event rate table with columns
        'log10(nu_E [GeV]) low', 'log10(nu_E [GeV]) center', 'log10(nu_E [GeV]) high', 'rate[livetime-1]'
        and adds a column 'rate_in_cone[livetime-1]' with the number of events within angle_max

        Parameters:
        - event_rate_table: table with event rate columns as specified
        - angle_max: cone_size in degrees

        Returns:
        - Dataframe with column representing the number of events within the cone
        """

        reconstructed_dataframe = event_rate_table[["log10(nu_E [GeV]) low", "log10(nu_E [GeV]) center", "log10(nu_E [GeV]) high", "rate [livetime^-1]"]].copy()

        reconstructed_dataframe["fraction_in_cone"] = 0.0
        reconstructed_dataframe["rate_in_cone [livetime^-1]"] = 0.0

        for index, row in reconstructed_dataframe.iterrows():
            true_logE_center = row["log10(nu_E [GeV]) center"]
            event_rate = row["rate [livetime^-1]"]
            fraction_below_angle = self.fraction_below_angle(true_logE_center, angle_max)
            reconstructed_dataframe.at[index, "fraction_in_cone"] = fraction_below_angle
            reconstructed_dataframe.at[index, "rate_in_cone [livetime^-1]"] = fraction_below_angle * event_rate

        return reconstructed_dataframe
