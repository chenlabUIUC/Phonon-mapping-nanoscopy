# Mode analysis

This code package is developed for the calculation of phonon dispersion relations of nanoparticle superlattices associated with the manuscript titled, "Nanoscale Phonon Dynamics in Topologically-Engineered, Self-Assembled Nanoparticle Lattices" by C. Qian et al. 

Getting started:

Download the "functions", "Data drift rotation corrected", and "Data raw" folders, and the following three codes. Run the codes directly in Matlab 2021a.

  - Step1_DriftRotationCorrection_20240709.m: Drift and rotation correction for the raw trajectories of nanoparticles. Input: Raw trajectories of nanoparticles in a superlattice in .mat format. Sample inputs are included in the folder of "Data raw". Output: Trajectories of nanoparticles in a superlattice in .mat format, which can be used as inputs for "Step2_1_PMN_rhombic_20240709.m" or "Step2_2_PMN_hexagonal_20240709.m" for mode and spring constant calculation.   
  - Step2_1_PMN_rhombic_20240709.m: main code for mode analysis, plotting, and fitting of nearest neighboring (NN) and angular (ANG) spring constants for rhombic or square lattice. Input: the output from "Step1_DriftRotationCorrection_20240709.m". Sample inputs are included in the folder of "Data drift rotation corrected" on cubes. Output: Plots of phonon dispersion relations with errors in the frequency denoted; fitted NN and ANG springs with errors.
  - Step2_2_PMN_hexagonal_20240709.m: main code for mode analysis, plotting, and fitting of NN spring constant for hexagonal lattice. Input: the output from "Step1_DriftRotationCorrection_20240709.m". Sample inputs included in the folder of "Data drift rotation corrected" on rods and prisms. Output: Plots of phonon dispersion relations with errors in the frequency denoted; fitted NN springs with errors.

Notes: For the Step1 and Step2 codes, inputs are set as the sample data provided in "Data raw" and "Data drift rotation corrected" folders. Comments are given in the codes on how to input your own experimental data following the sample input data formats. The inputs of tracking errors for Step2_1_PMN_rhombic_20240709.m and Step2_2_PMN_hexagonal_20240709.m, can be estimated using the workflow detailed in the Supplementary Information of the manuscript. 
