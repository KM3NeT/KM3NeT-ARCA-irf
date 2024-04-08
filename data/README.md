[![License](https://img.shields.io/badge/License-CC_BY_4.0-blueviolet.svg)](https://creativecommons.org/licenses/by/4.0/)

# Instrument response functions KM3NeT/ARCA230

Instrument Response Functions (IRFs) for the KM3NeT/ARCA230 detector. All functions are stored as .csv files for $\nu_\mu$ events selected as track or $\nu_e$ events selected as shower. The background is stored for the sum of contributions from all atmospheric neutrino flavours and interactions, and the background contribution from atmospheric muons.

## Effective area

The effective area is stored as a function of $\cos(\theta)$ and $\log_{10}(E_\nu \text{ [GeV]})$  where $\cos(\theta)$ = 1 corresponds to downgoing events. The average was taken over $\nu$ and $\bar{\nu}$ such that it can be multiplied with a combined $\nu$ and $\bar{\nu}$ neutrino flux $\Phi^{\nu+\bar{\nu}}$.

Files:
* **aeff_coszen_allnumuCC_track.csv** 
* **aeff_coszen_allnueCC_shower.csv** 

## Point spread function

The point spread function contains the event density per solid angle $dP/d\Omega$ versus the distance to the source $\log_{10}(\psi$ [degrees]). This function is stored for different energy ranges to account for the energy dependence of the direction reconstruction performance. The point spread function for $\bar{\nu}_\mu$ CC selected as track and $\bar{\nu}_e$ CC selected as shower are shared.

Files:
* **psf_numuCC_track.csv**
* **psf_nueCC_shower.csv**

## Energy response

The energy response is stored as the number of events per reconstructed energy for a given true neutrino energy. The energy response for $\bar{\nu}_\mu$ CC selected as track and $\bar{\nu}_e$ CC selected as shower are shared.

Files: 
* **energyresponse_numuCC_track.csv**
* **energyresponse_nueCC_shower.csv**

## Background

Sum of the background contributions from all atmospheric neutrino flavours and atmospheric muons for a given channel. The background is stored as a function of sin(declination) and reconstructed energy $\log_{10}(E_\nu \text{ [GeV]})$.

Files:
* **bkg_track.csv**
* **bkg_shower.csv**
