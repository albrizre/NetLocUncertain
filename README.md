R code for the paper *Dealing with location uncertainty for modeling network-constrained lattice data*, by Álvaro Briz-Redón.

- File "1_Simulation_study.R" allows generating the simulated data under the data-generating process described in Section 4 of the paper.
- File "2_Model_uncertain_locations_sim_gamma.R" allows fitting the uncertain segment model described in the paper.
- Folder *Functions* contains the file "sim_pattern_lpp.R" that simulates a point pattern over the road network of Valencia under the data-generating process described in Section 4 of the paper.
- Folder *Data* contains the data files needed for running the scripts: the shapefile of the road network of Valencia (in **linnet** format), and the covariates (X1 and X2) used for the simulation study.
