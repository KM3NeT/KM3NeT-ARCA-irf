[![License](https://img.shields.io/badge/License-BSD_3--Clause-blueviolet.svg)](https://opensource.org/licenses/BSD-3-Clause)


# KM3NeT/ARCA230 Instrument response functions

This repository contains the instrument response functions of the full KM3NeT/ARCA230 detector. The IRFs are accompanied by a set of scripts to interact with the IRFs and to perform an example cut-and-count analysis to calculate the sensitivity and discovery potential to a neutrino point source.

**N.B.**: The resulting sensitivity and discovery potential is worse than presented in the paper due to:
* The cut-and-count method only looks at the track (or shower) channel instead of combining both,
* This analysis only includes signal from $\nu_\mu$ and $\bar{\nu}_\mu$ CC events selected as track and $\nu_e$ and $\bar{\nu}_e$ CC events selected as shower, instead of all flavours and interactions,
* The paper uses a more sophisticated method than presented here. The paper uses a binned likelihood method and throws pseudo experiments to determine the sensitivity, while in this example we use Poisson statistics for a simple counting experiment.

## Content

* **data/**: Instrument Response Functions (IRFs) for the KM3NeT/ARCA230 detector.
* **analysis/**: Jupyter notebooks with example plots and analysis
* **src/arca230/**:
    * **flux.py**: Class that represents a single power law neutrino point source flux.
    * **aeff.py**: Class that loads the effective area and calculates event rates using a point source flux.
    * **psf.py**: Class that loads the point spread function and calculates probabilities to reconstruct events with a specified search cone size.
    * **energyresponse.py**: Class that loads the energy response and convolves true neutrino energies with the energy response of the detector.
    * **background.py**: Class that calculates expected background rates at different positions in the sky.

## Installation

### Download

The content of the repository is downloaded via `git`:

```sh
git clone git@git.km3net.de:open-data/public-candidates/open-point-source-search.git
```

### Creating the environment

#### Using venv


Create a virtual environment

```sh
python -m venv my_venv
```

Source the virtual environment

```sh
source my_venv/bin/activate
```

### Install

Enter the dowloaded repository
```sh
cd open-point-source-search
```
and install the requirements

```sh
pip install -e .
```

### Running Jupyter notebooks

In order to run the notebooks, you need to install Jupyter, by using `pip install jupyter` or following the instructions at the [Juypter website](https://jupyter.org/install).

### Running the Jupyter kernel

From within your virtual environment, create a Jupyter kernel and launch your notebook:
```sh
python -m ipykernel install --user --name=km3net_ps
jupyter-notebook
```

For `zsh` shell, you need to execute these lines first before installation of the kernel
```zsh
conda install -c conda-forge notebook
conda install -c conda-forge nb_conda_kernels
```

You can then execute the notebooks in your browser following the URL in the stdout.
