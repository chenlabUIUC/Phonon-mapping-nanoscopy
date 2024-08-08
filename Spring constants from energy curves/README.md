# Spring constants from energy curves

This code package is developed for the calculation of spring constants based on internanoparticle interactions, associated with the manuscript titled, "Nanoscale Phonon Dynamics in Topologically-Engineered, Self-Assembled Nanoparticle Lattices" by C. Qian et al. 

Getting started:

Download the sample data in 'Data/' and the 4 codes. Run the codes directly in Matlab 2021a.

  - Fit_NNspring.m: Generate fitted constant for the nearest neighboring (NN) springs.
  - Fit_ANGSpring.m: Generate fitted constant for the angular springs considering beyond-NN interactions.
  - FitSpringg_rod_20240806.m: Fit spring constant for rod pairwise interaction.
  - FitSpringg_prism_20240806.m: Fit spring constant for prism column pairwise interaction.

Notes: For the codes, inputs are set as the sample data provided. The colloidal interaction potentials of other nanoparticles can be also used as the inputs in the same data format as our sample data to calculate the spring constants. The choice of fitting range is detailed in the Supplementary Information of the manuscript. 
