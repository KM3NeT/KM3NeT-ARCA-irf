import pandas as pd


class EnergyResponse:
    """
    Loads the energy response for the ARCA230 detector. The data is stored as a dataframe with
    a reconstructed energy distribution for each true neutrino energy.
    """

    def __init__(self, file_path="../data/energyresponse_numuCC_track.csv"):
        self.file_path = file_path
        self.eresponse_data = None
        self.load_eresponse_data()

    def load_eresponse_data(self):
        """
        Loads the data for the energy response
        """
        try:
            self.eresponse_data = pd.read_csv(self.file_path, delimiter=",")
            print("Energy response data loaded successfully.")
        except FileNotFoundError:
            print(f"Error: File '{self.file_path}' not found.")
        except Exception as e:
            print(f"An error occurred while loading the data: {e}")

    def fraction_between_energy(self, logE, cutoff_low_logerec, cutoff_high_logerec):
        """
        Filters the energy response data for the given true neutrino energy
        and calculates the fraction of events reconstructed within the specified bounds

        Parameters:
        - logE: log of the true neutrino energy in GeV
        - cutoff_low_logerec: lower bound of the log of the reconstructed energy in GeV
        - cutoff_high_logerec: higher bound of the log of the reconstructed energy in GeV

        Returns:
        - Fraction of events within specified reconstructed energy range
        """

        logE_mask = (self.eresponse_data["log10(nu_E [GeV]) low"] <= logE) & (self.eresponse_data["log10(nu_E [GeV]) high"] > logE)

        filtered_rows = self.eresponse_data[logE_mask]

        weight = 0
        norm = 0

        for _, row in filtered_rows.iterrows():
            logE_low = row["log10(reco_E [GeV]) low"]
            logE_high = row["log10(reco_E [GeV]) high"]

            # Calculate the intersection of the ranges
            intersection_low = max(logE_low, cutoff_low_logerec)
            intersection_high = min(logE_high, cutoff_high_logerec)

            # Calculate the fraction in the range
            fraction_in_range = max(0, intersection_high - intersection_low) / (logE_high - logE_low)

            norm += row["dP/dlog10(nu_E [GeV])"]
            weight += row["dP/dlog10(nu_E [GeV])"] * fraction_in_range

        if norm == 0:
            return 0
        return weight / norm

    def energy_response(self, logE):
        """
        Filters the energy response data for the given true neutrino energy
        and returns the normalised distribution of possible reconstructed energies

        Parameters:
        - logE: log of the true neutrino energy in GeV

        Returns:
        - Normalised distribution of possible reconstructed energies
        """
        logE_mask = (self.eresponse_data["log10(nu_E [GeV]) low"] <= logE) & (
            self.eresponse_data["log10(nu_E [GeV]) high"] > logE
        )

        filtered_rows = self.eresponse_data[logE_mask].copy()
        filtered_rows["dP/dlog10(nu_E [GeV])"] = filtered_rows["dP/dlog10(nu_E [GeV])"] / sum(
            filtered_rows["dP/dlog10(nu_E [GeV])"]
        )
        return filtered_rows

    def reconstruct_event_table(self, event_rate_table):
        """
        Copies an event rate table with columns
        'log10(nu_E [GeV]) low', 'log10(nu_E [GeV]) center', 'log10(nu_E [GeV]) high', 'rate[livetime-1]', 'rate_in_cone[livetime-1]'
        and converts the true neutrino energies to the reconstructed energies according to the energy response of the detector

        Parameters:
        - logE: logarithm of the true neutrino energy in GeV
        - angle_max: cone_size in degrees

        Returns:
        - Dataframe with column representing the number of events within the cone
        """

        # copy the input dataframe and move from true to reconstructed variables
        reconstructed_dataframe = event_rate_table[
            [
                "log10(nu_E [GeV]) low",
                "log10(nu_E [GeV]) center",
                "log10(nu_E [GeV]) high",
                "rate [livetime^-1]",
                "rate_in_cone [livetime^-1]",
            ]
        ].copy()

        reconstructed_dataframe.rename(
            columns={
                "log10(nu_E [GeV]) low": "log10(reco_E [GeV]) low",
                "log10(nu_E [GeV]) center": "log10(reco_E [GeV]) center",
                "log10(nu_E [GeV]) high": "log10(reco_E [GeV]) high",
            },
            inplace=True,
        )

        reconstructed_dataframe["rate [livetime^-1]"] = 0
        reconstructed_dataframe["rate [livetime^-1]"] = reconstructed_dataframe["rate [livetime^-1]"].astype(float)

        reconstructed_dataframe["rate_in_cone [livetime^-1]"] = 0
        reconstructed_dataframe["rate_in_cone [livetime^-1]"] = reconstructed_dataframe["rate_in_cone [livetime^-1]"].astype(float)

        # iterate over each true energy
        for _, row_true in event_rate_table.iterrows():
            true_logE_center = row_true["log10(nu_E [GeV]) center"]

            rate = row_true["rate [livetime^-1]"]
            rate_in_cone = row_true["rate_in_cone [livetime^-1]"]

            # iterate over reconstructed energy
            for index_reco, row_reco in reconstructed_dataframe.iterrows():
                reco_E_low = row_reco["log10(reco_E [GeV]) low"]
                reco_E_high = row_reco["log10(reco_E [GeV]) high"]

                # calculates how many events are within the reconstructed energy bin
                reconstructed_rate = rate * self.fraction_between_energy(true_logE_center, reco_E_low, reco_E_high)
                reconstructed_rate_in_cone = rate_in_cone * self.fraction_between_energy(
                    true_logE_center, reco_E_low, reco_E_high
                )

                reconstructed_dataframe.at[index_reco, "rate [livetime^-1]"] += reconstructed_rate
                reconstructed_dataframe.at[index_reco, "rate_in_cone [livetime^-1]"] += reconstructed_rate_in_cone

        return reconstructed_dataframe
