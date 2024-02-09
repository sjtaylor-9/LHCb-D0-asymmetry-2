# Measuring the Production Asymmetry between $D^0$ and  $\bar{D}^0$ mesons due to proton-proton collisions at \sqrt{s} = 13 TeV at the LHCb
This project is aimed to calculating the production asymmetry of the $D^0$ meson. This is done in differents regions of the phase space spanned by the transverse momentum ($p_T$) and pseudorapidity ($\eta$) of the $D^0$ meson. This is the current version (9th February) of the code. The project was done in collaboration with Laxman Seelan.

In this repository there are the necessary tools in order to:
  - (in progress)

The ```environment.yaml``` file inclues all of the necessary python libraries needed, apart from ```LHCbStyle```, which can be installed using: ```mamba install -c conda-forge lhcbstyle```. 

**Warnings:**

While running the code be aware that any change to one of the scripts can lead to a malfunction. In addition, make sure that the directories that will be generated while running the program don't already exist. If they do exist beforehand, this program might not work as intended. ```model_fitting.py``` produces a segmentation violation, however this bug does not affect the code.
In the event that condor does not recognise ```utils.py``` as a library import, it can be added to a conda envrionment as a pip module by ```pip install -e utils```.
```fit_control_modes.py``` returns an error for versions of ```root``` later than 6.28.0. If the conda environment used uses the ```environment.yaml``` file in this repository then this error will be avoided.

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
A miniforge environment can be created with the libraries from ```ennvironment.yaml``` installed using ```mamba env create -f environment.yaml -n D0-asymmetry```.

The different scripts can be run individually (note that a different set of arguments is required for each), or as a whole using the bash script *main.sh*.
In order to use *main.sh* 4 arguments are required. These are:
- The path where the output should be written
- The year the data to be used was taken [16, 17, 18]
- The size or amount of data to be used [small, medium, large, 1, 2, 3, 4, 5, 6, 7, 8]
- Whether a binned fit should be performed (otherwise unbinned fit) [y, Y, n, N]

Here is an example of how to call *main.sh*:
```
bash main.sh example 18 large y
```
This should produce the same output as shown in the folder *example* (still to be implemented).
## Credits
A large amount of the scripts uses or is inspired by the code written by Camille Jarvis-Stiggants and Michael England during their MPhys project and Marc Oriol Pérez in his summer internship.


**Authors:** Sam Taylor (samuel.taylor-9@student.manchester.ac.uk) and Laxman Seelan (laxman.seelan@student.manchester.ac.uk)/ **Last modified:** 9th February 2024
