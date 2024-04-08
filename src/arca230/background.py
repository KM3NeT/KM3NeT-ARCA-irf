import pandas as pd
import numpy as np


class BackgroundComponent:
    """
    Loads the expected background dataset of the ARCA230 detector for the track or shower channel
    This function contains the sum of background from atmospheric muons
    and atmospheric neutrinos
    """

    def __init__(self, file_path="../data/bkg_track.csv"):
        self.sindec_binwidth = 0.05
        self.file_path = file_path
        self.background_data = None
        self.load_background_data()

    def load_background_data(self):
        """
        Loads the data for the background
        """
        try:
            self.background_data = pd.read_csv(self.file_path, delimiter=",")
            print("Background data loaded successfully.")
        except FileNotFoundError:
            raise RuntimeError(f"File '{self.file_path}' not found.")
        except Exception as e:
            raise RuntimeError(f"An error occurred while loading the data: {e}")

    def event_rate(self, sindec, angle_max, livetime=365.25 * 24 * 60 * 60):
        """
        Calculate the background event rate as a function of reconstructed energy
        based on a given position in the sky

        Parameters:
        - sindec: Source location
        - angle_max : size of the search cone in degrees
        - livetime: in seconds

        Returns:
        - Dataframe with the event rate per reconstructed energy
          corresponding to the source sin(dec)
        """
        try:
            # Select rows where sindec is within the bounds
            selected_rows = self.background_data[(self.background_data["sin(dec) low"] <= sindec) & (sindec < self.background_data["sin(dec) high"])].copy()

            if not selected_rows.empty:
                fraction_in_cone = (
                    4 * np.pi * np.power(np.sin(angle_max * np.pi / 180 / 2), 2) / (2 * np.pi * self.sindec_binwidth)
                )

                selected_rows["rate [livetime^-1]"] = selected_rows["rate [s^-1]"] * livetime
                selected_rows["rate_in_cone [livetime^-1]"] = selected_rows["rate [livetime^-1]"] * fraction_in_cone

                return selected_rows
            else:
                print(f"No rows found for the given sindec value {sindec}.")
                return None
        except KeyError as e:
            print(f"Error: {e} not found in the loaded data columns.")
        except Exception as e:
            print(f"An error occurred while calculating event rate: {e}")

