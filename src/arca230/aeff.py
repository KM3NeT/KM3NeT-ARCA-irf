import pandas as pd
import numpy as np

from arca230.coordinates import fraction_at_zenith


class EffectiveArea:
    """
    Loads the effective area of the ARCA230 detector for numu selected as track or nue selected as shower
    The effective area is stored as a function of cos(zen) and true neutrino energy
    """

    def __init__(self, file_path="../data/aeff_coszen_numu_track.csv"):
        self.file_path = file_path

        self.effective_area_data = None
        self.load_effective_area_data()

        self.binsize_coszen = 0.05
        self.tolerance = 1e-5  # for comparing floats

    def load_effective_area_data(self):
        """
        Loads the data for the effective area
        """
        try:
            self.effective_area_data = pd.read_csv(self.file_path, delimiter=",")
            print("Effective Area data loaded successfully.")
        except FileNotFoundError:
            print(f"Error: File '{self.file_path}' not found.")
        except Exception as e:
            print(f"An error occurred while loading the data: {e}")

    def effective_area_at_sindec(self, sindec, nsamples=1000):
        """
        Calculate the effective area for a source location. This is obtained by weighting the effective area
        cos(zen) bands with the visibility for that cos(zen)

        Parameters:
        - sindec: Source location
        - nsamples: number of samples in calculation visibility

        Returns:
        - Dataframe with the effective area as a function of true neutrino energy
        """
        if np.abs(sindec) > 1:
            raise ValueError(f"abs(sindec) should be < 1 {sindec}")

        # copy the zenith bands
        zenith_bands = self.effective_area_data[["cos(zen) low", "cos(zen) center", "cos(zen) high"]].copy()
        zenith_bands = zenith_bands.drop_duplicates(keep="first")

        # calculate the visibility for each zenith band
        zenith_bands = fraction_at_zenith(zenith_bands, sindec, nsamples)  # src/coordinates.py

        # make a copy of the effective area for storing the result
        effective_area_source = self.effective_area_data[["log10(nu_E [GeV]) low", "log10(nu_E [GeV]) center", "log10(nu_E [GeV]) high", "aeff [m^2]"]].copy()
        effective_area_source["aeff [m^2]"] = 0
        effective_area_source = effective_area_source.drop_duplicates(keep="first")

        # iterate over each zenith band
        for _, row_zenith in zenith_bands.iterrows():
            # Extract values for calculations
            cos_zen_low = row_zenith["cos(zen) low"]

            # visibility
            weight = row_zenith["weight"]

            effective_area_zenith_band = self.effective_area_data[self.effective_area_data["cos(zen) low"] == cos_zen_low]

            merged_df = pd.merge(
                effective_area_source,
                effective_area_zenith_band,
                on=["log10(nu_E [GeV]) low", "log10(nu_E [GeV]) center", "log10(nu_E [GeV]) high"],
                how="left",
                suffixes=("_df1", "_df2"),
            )

            merged_df["aeff [m^2]"] = merged_df["aeff [m^2]_df1"] + merged_df["aeff [m^2]_df2"] * weight

            effective_area_source = merged_df[
                ["log10(nu_E [GeV]) low", "log10(nu_E [GeV]) center", "log10(nu_E [GeV]) high", "aeff [m^2]"]
            ]

        return effective_area_source[effective_area_source["aeff [m^2]"] > 0]

    def effective_area_zenith_band(self, cos_zen_low, cos_zen_high):
        """
        Calculate the effective area for a zenith band. The instrument response function is already stored as a function
        of cos(zen), but this function is averaged over the full-sky. When selecting a zenith band one needs to average over
        the number of zenith bands.

        Parameters:
        - cos_zen_low: low bound of the cos(zen) band
        - cos_zen_high: high bound of the cos(zen) band

        Returns:
        - Dataframe with the effective area as a function of true neutrino energy for the given zenith band
        """
        rounded_binsize = round(self.binsize_coszen, 5)  # Round binsize to the desired precision

        if (
            abs(round(cos_zen_low / rounded_binsize) - cos_zen_low / rounded_binsize) > self.tolerance
            or abs(round(cos_zen_high / rounded_binsize) - cos_zen_high / rounded_binsize) > self.tolerance
        ):
            raise ValueError(f"cos_zen_low and cos_zen_high need to be factors of the binning of 0.05: {cos_zen_low} {cos_zen_high}")

        if cos_zen_low > cos_zen_high:
            raise ValueError(f"cos_zen_high needs to be higher than cos_zen_low: {cos_zen_low} {cos_zen_high}")

        if abs(cos_zen_low) > 1 or abs(cos_zen_high) > 1:
            raise ValueError(f"cos_zen_low and cos_zen_high need to be within 1 and -1: {cos_zen_low} {cos_zen_high}")

        # how many zenith bands are merged
        number_zenith_bands = round((cos_zen_high - cos_zen_low) / self.binsize_coszen)

        # copy the first zenith band
        aeff_data = self.effective_area_data
        aeff_band = aeff_data[
            (aeff_data["cos(zen) low"] == cos_zen_low)
            & (aeff_data["cos(zen) high"] == round(cos_zen_low + self.binsize_coszen, 2))
        ].copy()
        aeff_band = aeff_band.drop(columns=["cos(zen) low", "cos(zen) center", "cos(zen) high"])

        # iterate over the other zenith bands
        cos_zen_list = np.linspace(cos_zen_low, cos_zen_high, number_zenith_bands + 1)

        for i in range(1, len(cos_zen_list) - 1):
            cos_zen_low_tmp = round(cos_zen_list[i], 2)
            cos_zen_high_tmp = round(cos_zen_list[i] + self.binsize_coszen, 2)

            aeff_band_tmp = aeff_data[(aeff_data["cos(zen) low"] == cos_zen_low_tmp) & (aeff_data["cos(zen) high"] == cos_zen_high_tmp)].copy()
            aeff_band_tmp = aeff_band_tmp.drop(columns=["cos(zen) low", "cos(zen) center", "cos(zen) high"])

            aeff_band_merged = pd.merge(
                aeff_band,
                aeff_band_tmp,
                on=["log10(nu_E [GeV]) low", "log10(nu_E [GeV]) center", "log10(nu_E [GeV]) high"],
                how="left",
                suffixes=("_df1", "_df2"),
            )
            aeff_band_merged["aeff [m^2]"] = aeff_band_merged["aeff [m^2]_df1"] + aeff_band_merged["aeff [m^2]_df2"]
            aeff_band_merged = aeff_band_merged.drop(columns=["aeff [m^2]_df1", "aeff [m^2]_df2"])

            aeff_band = aeff_band_merged.copy()

        aeff_band["aeff [m^2]"] = aeff_band["aeff [m^2]"] / number_zenith_bands

        return aeff_band[aeff_band["aeff [m^2]"] > 0]

    def event_rate(self, flux, sindec, livetime=365.25 * 24 * 60 * 60):
        """
        Calculate the event rate based on a given flux and position in the sky: sin(dec)

        Parameters:
        - flux: PointSourceFlux object (see flux.py)
        - sindec: Source location
        - livetime: in seconds

        Returns:
        - Dataframe with the event rate appended to the selected effective area rows
          corresponding to the source sin(dec)
        """

        effective_area_source = self.effective_area_at_sindec(sindec)

        effective_area_source["energy_bin_width"] = np.power(
            10, effective_area_source["log10(nu_E [GeV]) high"]
        ) - np.power(10, effective_area_source["log10(nu_E [GeV]) low"])

        effective_area_source["rate [livetime^-1]"] = (
            effective_area_source["aeff [m^2]"]
            * flux.dNdE(effective_area_source["log10(nu_E [GeV]) center"])
            * effective_area_source["energy_bin_width"]
            * livetime
        )

        return effective_area_source
