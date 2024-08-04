# Phonon mapping nanoscopy

This code package is developed to perform the phonon mode analysis for nanoparticle superlattices as potential mechanical metamaterials, particularly associated with the manuscript titled, "Nanoscale Phonon Dynamics in Topologically-Engineered, Self-Assembled Nanoparticle Lattices" by C. Qian et al. Particularly, this package includes the following three packages developed by the Qian Chen group at University of Illinois, Xiaoming Mao group at University of Michigan, and Wenxiao Pan group at University of Wisconsin, Madison. 

(1) Mode analysis: to calculate and plot the phonon dispersion relations (the frequency-wavevector relations) of nanoparticle superlattices with the errors of phonon frequency calculated and noted, using the particle centroid positions obtained from liquid-phase TEM as inputs. This code package also fits the calculated phonon dispersion relations to derive the effective nanoscale spring constants in the nanoparticle superlattices. Errors are calculated and output. Here the particle centroid positions can be obtained using the U-Net neural network method presented in L. Yao, Z. Ou, B. Luo, C. Xu, Q. Chen, ACS Central Science 6, 1421 (2020); DOI:10.1021/acscentsci.0c00430. The tracking method is detailed in the Supplementary Information of the manuscript, with codes and examples available at https://github.com/chenlabUIUC/simulatedLPTEM and https://github.com/chenlabUIUC/U-NetExample. The particle centroid positions can also be obtained using the method detailed in the work Z. Ou, Z. Wang, B. Luo, E. Luijten, Q. Chen, Nature Materials 19, 450 (2020); DOI:10.1038/s41563-019-0514-1

(2) Spring constants from energy curves: to fit the spring constants of the nanoparticle superlattices from the colloidal interparticle interaction energies calculated independently from coarse-grained models. Depending on the lattice types, the calculations of both nearest neighboring (NN) springs and the angular (ANG) springs involving beyond-NN are included.

(3) Brownian dynamics simulation: to simulate the self-assembly dynamics of nanoparticles.

Date: 08/2024.

The mode analysis and spring constant parts are developed and tested using Matlab 2021a.  

The BD simulation part is developed and tested using LAMMPS (29sep2021) and Matlab 2021a.

For more information about the project, algorithms, and related publications please refer to the manuscript titled, "Nanoscale Phonon Dynamics in Topologically-Engineered, Self-Assembled Nanoparticle Lattices" by C. Qian (currently under review). 


## Getting started

The mode analysis and spring constant parts of the codes can be used upon downloading. All the data presented in the manuscript by C. Qian, including the particle trajectories of the five lattices, and the interaction energy curves are provided as examples to run and test the codes. These example data are included as '.mat' files. For detailed instructions for each package, please refer to [readme] in each folder.
