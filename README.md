# Measuring the Production Asymmetry between $D^0$ and  $\bar{D}^0$ mesons in proton-proton collisions at $\sqrt{s} =$ 13 TeV at the LHCb
This project is aimed to calculating the production asymmetry of the $D^0$ meson. This is done in differents regions of the phase space spanned by the transverse momentum ($p_T$) and pseudorapidity ($\eta$) of the $D^0$ meson. The project was done in collaboration with Laxman Seelan.

In this repository there are the necessary tools in order to:
 - Make a selection of the events given a certain criteria,
 - Remove multiple candidates,
 - Perform a global fit of the invariant mass distribution using a binned simultaneous extended maximum likelihood fit,
 - Plot the model of the invariant mass distribution using the global fit parameters,
 - Create a uniform binning across the phase space - ($p_T, \eta$) and individual $p_T$, $\eta$ binning schemes,
 - Plot the distriution of events across the phase-space, overlayed with the ($p_T, \eta$) binning scheme,
 - Perform a local fit of the invariant mass disribution in each of the phase-space bins, by using a simultaneous extended maximum likelihood fit and plotting,
 - Calculate global and local detection asymmetries using the same binning scheme as the raw asymmetry,
 - Calculate a global production asymmetry,
 - Use ```Pythia``` to make a theoretical prediction of the production asymmetry,
 - Plot the asymmetry in the bins of ($p_T, \eta$),
 - Plot the real and simulated production asymmetries against $p_T$ and $\eta$ and compare the results with a pull distribution.

**Warnings:**
 - ```model_fitting.py``` produces a segmentation violation, however this bug does not affect the code and is a memory issue with ```ROOT```,
 - In the event that condor does not recognise ```utils.py``` as a library import, it can be added to a conda envrionment as a pip module by ```pip install -e utils```,
 - ```fit_control_modes.py``` returns an error for versions of ```root``` later than 6.28.0. If the conda environment used uses ```environment.yaml``` in this repository, then this error will be avoided.
 - If the selection and detection asymmetry scipts are to be re-run then all of the file paths need to be changed as you cannot write to a different user's eos,
 - The local detection asymmetry must be run in batches as the eos storage is not large enough to hold all of the ```.root``` files at once.

## How to download
In order to download this package you can use the following commands in your terminal:
```
mkdir ProductionAsymmetry
cd ProductionAsymmetry
git init
git clone git@github.com:sjtaylor-9/LHCb-D0-asymmetry-2.git
git pull origin main
```

## How to use
A miniforge environment containing the libraries from ```ennvironment.yaml``` can be created using ```mamba env create -f environment.yaml -n D0-asymmetry```. ```environment.yaml``` includes all of the necessary ```Python``` libraries needed, apart from ```LHCbStyle```, which can be installed using: ```mamba install -c conda-forge lhcbstyle```. 

The different scripts can be run individually (note that a different set of arguments is required for each), or as a whole using the bash script ```main_processing.sh```.
In order to use ```main_processing.sh``` 5 arguments are required. These are:
- The path where the output should be written,
- The year the data to be used was taken [16, 17, 18],
- The size or amount of data to be used [1-800]. The integers must be in steps of 10,
- Whether a binned fit should be performed (otherwise unbinned fit) [y, Y, n, N],
- The name of the model to be used.

Here is an example of how to call ```main_processing.sh```:
```
bash main_processing.sh 2018 18 100 y Model_1
```
## Credits
A large amount of the scripts uses or is inspired by the code written by Camille Jarvis-Stiggants and Michael England during their MPhys project and Marc Oriol Pérez in his summer internship. The detection asymmetry scripts used by ```HTCondor``` and the ```utils.py``` script were written by Aodhan Burke for his PhD thesis and the ```runpythia.cpp``` script was written by Suzanne Klaver for her study of the $D^\pm$ mesons.


**Authors:** Sam Taylor (samuel.taylor-9@student.manchester.ac.uk) and Laxman Seelan (laxman.seelan@student.manchester.ac.uk)/ **Last modified:** 5th April 2024
