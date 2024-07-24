# Maxwell Lattice Project

This code package is developed to perform phonon mode analysis, fit spring constants from calculated energy curves, and perform Brownian dynamics (BD) simulation.

Data: 09/2023.

The mode analysis and spring constant parts are developed and tested using Matlab 2021a.  

The BD simulation part is developed and tested using LAMMPS (29sep2021) and Matlab 2021a.

For more information about the project, algorithms, and related publications please refer to the [Chen Group website](https://chenlab.matse.illinois.edu/).  

## Reference
If you find our codes are helpful to your publication, please cite:

Chang Qian, Ethan Stanifer, Zhan Ma, Binbin Luo, Lehan Yao, Chang Liu, Wenxiao Pan,  Xiaoming Mao, Qian Chen, "Nanoscale Imaging of Phonon Modes and Reconfiguration in Topologically-Engineered Nanoparticle Lattices" submitted to _Science_

## Getting started

The mode analysis and spring constant parts of the codes can be used upon downloading. All five examples of particle trajectory are provided for the mode analysis, and 3 energy curves are provided for nearest neighbor (NN) spring and angular (ANG) spring fitting. The example data are included as '.mat' files.
- Mode analysis：
  - Mode analysis/Step1_DriftRotationCorrection_20240709.m: Drift and rotation correction for all data matrices
  - Mode analysis/Step2_1_PMN_rhombic_20240709.m: main code for mode analysis and plotting for rhombic lattice of nanocubes
  - Mode analysis/Step2_2_PMN_hexagonal_20240709.m: mode analysis and plotting for hexagonal lattice of nanorods and nanoprisms
- Spring constants fitting:
  - Spring constants from energy curves/Fit_NNSpring: NN spring fitting
  - Spring constants from energy curves/Fit_ANGSpring: ANG spring fitting
- BD simulation:
  - Please refer to [readme](https://github.com/chenlabUIUC/MaxwellLattice/blob/main/Brownian%20dynamics%20simulation/Readme.md).
