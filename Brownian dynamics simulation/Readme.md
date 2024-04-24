# Brownian dynamics simulation
The code package is developed for Brownian dynamic simulation for nanocube assembly in the publication
Date: 09/2023
Test: LAMMPS (29sep2021), Matlab 2021b
For more information about the project, algorithms, and related publications please refer to the [Chen Group website](https://chenlab.matse.illinois.edu/). 

## Getting started
1. Download the code packages.
2. Move the source codes in the folder of ‘Lammps_src/’ to the source code folder of LAMMPS, i.e. Lammps/src. Replace if there are the files with same name.
3. Compile LAMMPS again.
4. Change the line of 5 in ‘in.lammps’ to make sure it refers to the location of folder ‘I02_config’. Change the line of 11 in ‘in.lammps’ to make sure it refers to the location of folder ‘Txt_param_Fsh_NN_I022/’ and ‘Txt_param_Fvdw_Fsh_NN_220227’. 
5. Run ‘in.lammps’ for the simulation.
6. Visualize the results with Matlab by running ‘make_movie.m’.

## Notes
The multilayer perceptron (MLP) model in this work has a layer size of [3, 256, 128, 64, 3], with hyperbolic tangent (tanh) as the activation function. The weights and bias of each layer are stored in ‘Txt_param_Fsh_NN_I022/’ and ‘Txt_param_Fvdw_Fsh_NN_220227’, for van der Waals and electrostatic repulsive interaction, respectively. The models are trained on data generated from coarse-grained (CG) modeling. For more details, please refer to Note S2 and Note S4 in the publication.
