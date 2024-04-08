# Analysis

The example analysis in **analysis/analysis.ipynb** calculates the expected number of signal events from a single power law neutrino flux by convolving it with the effective area of the KM3NeT/ARCA230 detector. The obtained event distribution as a function of the true neutrino energy is smeared in energy and direction by the energy response and point spread function. These functions tell us how true neutrino energies and directions are reconstructed. The expected signal events are compared with the number of background events to calculate the confidence level and p-value using Poisson statistics. Furthermore, this repository contains notebooks that plot the instrument response functions.

## Contents

* **analysis.ipynb**: Analysis Jupyter notebook with an example cut-and-count analysis for a neutrino point source.
* **plot_background.ipynb**: Plo the background distribution.
* **plot_effective_area.ipynb**: Example of how to plot the effective area for different zenith bands or declinations.
* **plot_energy_response.ipynb**: Plot the energy response of the detector
* **plot_event_rate.ipynb**: Calculate expected event rates of the detector
* **plot_point_spread_function.ipynb**: Example of how to use the point spread function
* **plot_visibility.ipynb**: Plot the visibility of the KM3NeT/ARCA detector.