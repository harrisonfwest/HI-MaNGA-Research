# HI-MaNGA

This repository houses my thesis research with Professor Karen Masters investigating a 
correlation between galaxy misalignment and HI concentration, and between global gas profile asymmetry and HI concentration, in the SDSS-IV/MaNGA survey.
Galaxy misalignment describes a difference in the position angles of a galaxy's stellar and gas discs, which can range from 0 degrees (well-aligned) to 180 degrees (counter-rotating), 
and HI concentration describes how much non-ionized ('cold') hydrogen gas is present in a galaxy compared to the 
expected (average) amount at the galaxy's given stellar mass. Global gas profile asymmetry refers to deviation from 1 in the ratio
of integrated fluxes on the left and right sides of the galaxy's global gas profile about a given central velocity.
HI mass is particularly easy to detect due to its distinct 21 cm spectral emission line, meaning that we are able to gather 
relatively abundant data for the HI concentrations of galaxies up to significant redshifts.

More information about the MaNGA (Mapping Nearby Galaxies at Apache Point Observatory) survey, part of the Sloan Digital
Sky Survey, can be found on the [MaNGA section of the SDSS website](https://www.sdss4.org/surveys/manga/).

My code can be found in the `S26_Thesis-Code.ipynb` Jupyter Notebook. The galaxy spectra can be found in the `Spectra` folder. The `Data` directory folder contains many iterations of my cross-matched galaxy data,
though the most current version is `HIMaNGADR4Cross_sampgalMA.fits`.

The counter rotating galaxy data is from [Zhou et. al (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.515.5081Z/abstract). Using the [TOPCAT software](https://www.star.bris.ac.uk/~mbt/topcat/), I combined the
misaligned galaxies from the counter-rotating data, the original HI-MaNGA DR4 data, and [Pipe3D](https://www.sdss4.org/dr17/manga/manga-data/manga-pipe3d-value-added-catalog/)
data files, matching their MaNGA IDs.

The survival analysis code used (`survival_analysis.py`) is from [the repository of David Stark](https://github.com/dvstark/survival) 
at the Space Telescope Science Institute (STScI). My repository contains an example of this code's usage: the code
has been copied and slightly modified in the notebook `survival_analysis_example.ipynb` before being implemented into 
the main Jupyter Notebook (`HI-MaNGA_MIsalignmentAndDeficiency.ipynb`). Many thanks to David Stark for providing this 
code, which wraps survival analysis methods from the R Statistical Software into Python, and to my Undergraduate PI 
Karen Masters for referring me to David's repository as well as for her broader research guidance.
