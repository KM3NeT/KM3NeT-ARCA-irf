from astropy.coordinates import EarthLocation, AltAz, ICRS
from astropy.time import Time
import astropy.units as u
import numpy as np
import pandas as pd


detector_latitude = 0.633407  # [radians]
detector_longitude = 0.278819  # [radians]
detector_location = EarthLocation.from_geodetic(lat=detector_latitude * u.rad, lon=detector_longitude * u.rad)
observing_time = Time("2020-01-01 02:02:02")  # dummy observing time needed for the functions


def fraction_at_zenith(zenith_dataframe, sindec, nsamples=1000):
    """
    Calculates the relative time spent per day for
    a source at given sindec for the specified zenith bands

    Parameters:
    - zenith_dataframe: Pandas dataframe with zenith bands ('cos(zen) low', 'cos(zen) high')
    - sindec: Source location
    - nsamples: number of samples in right ascension

    Returns:
    - Dataframe with the zenith bands and the fraction of time spent in each zenith band
    """

    right_ascensions = np.array([i * 2 * np.pi / nsamples for i in range(nsamples)])
    declinations = np.array([np.arcsin(sindec) for i in range(nsamples)])

    # Create ICRS coordinates from declination and right ascension
    icrs_coords = ICRS(ra=right_ascensions * u.rad, dec=declinations * u.rad)

    # Transform ICRS coordinates to AltAz coordinates for the observer location and time
    altaz_coords = icrs_coords.transform_to(AltAz(obstime=observing_time, location=detector_location))

    # Get zenith and azimuth from AltAz coordinates
    zeniths = 90.0 * u.deg - altaz_coords.alt
    coszen = np.cos(zeniths.value * np.pi / 180)

    for index, row in zenith_dataframe.iterrows():
        # Update counts for rows where the condition is satisfied
        count = np.sum((coszen >= row["cos(zen) low"]) & (coszen < row["cos(zen) high"]))
        zenith_dataframe.loc[index, "weight"] = count / nsamples

    return zenith_dataframe


def visibility_below_coszen(cos_zen_cut, sindec, nsamples=1000):
    """
    Calculates the visibility - the fraction of time in a day -
    for which the source is below the cos_zen_cut

    Parameters:
    - cos_zen_cut: Cosine of the zenith that defines visible or not
    - sindec: Source location
    - nsamples: number of samples in right ascension

    Returns:
    - Fraction of time spent below cos_zen_cut
    """
    right_ascensions = np.array([i * 2 * np.pi / nsamples for i in range(nsamples)])
    declinations = np.array([np.arcsin(sindec) for i in range(nsamples)])

    # Create ICRS coordinates from declination and right ascension
    icrs_coords = ICRS(ra=right_ascensions * u.rad, dec=declinations * u.rad)

    # Transform ICRS coordinates to AltAz coordinates for the observer location and time
    altaz_coords = icrs_coords.transform_to(AltAz(obstime=observing_time, location=detector_location))

    # Get zenith and azimuth from AltAz coordinates
    zeniths = 90.0 * u.deg - altaz_coords.alt
    coszen = np.cos(zeniths.value * np.pi / 180)

    return np.sum(coszen < cos_zen_cut) / nsamples


if __name__ == "__main__":
    zenith_dataframe = pd.DataFrame(columns=["cos(zen) low", "cos(zen) high", "cos(zen) center"])

    # Define the range for iterations
    num_rows = 40

    # Iterate and add rows to the dataframe
    for i in range(num_rows):
        # Example values, you can replace these with your own logic
        cos_zen_low = -1 + i * 2 / num_rows
        cos_zen_high = -1 + 0.05 + i * 2 / num_rows
        cos_zen_center = (cos_zen_high + cos_zen_low) / 2

        # Add a new row to the dataframe
        zenith_dataframe = zenith_dataframe.append(
            {"cos(zen) low": cos_zen_low, "cos(zen) high": cos_zen_high, "cos(zen) center": cos_zen_center},
            ignore_index=True,
        )

    # Print the resulting dataframe
    fraction_at_zenith(zenith_dataframe, -0.775, 1000)
