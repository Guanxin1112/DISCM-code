The project "DISCM-code" is created for the paper "Estimation and Inference for Dynamic Single-Index Varying-Coefficient Models" submitted to Statistica Sinica.

This repository includes codes for two simulation studies in the Manuscript, and one simulation study in the Supplementary materials.

1. The folder "Functions" contains eight R files:

1.1. three_step:   the R code for estimating the varying-coefficient functions and index functions by iteration process described in Section 2.

1.2. myknots.R:  the R code for choosing optimal knots for varying-coefficient function and index function with a specified order of B-spline. 

1.3. betaesti.R:   the R code for parametric estimation and nonparametric function estimation of single index varying-coefficient model.

1.4. betaknot.R:  the R code for selecting optimal knots of index functions for single index varying-coefficient model.

1.5 tune_select.R:  the R code for variable selection procedure in Section 4.1.

1.6. sim1data.R:  the R code for data generation of example 1.

1.7. sim2data.R:  the R code for  data generation of example 2.

1.8. sim3data.R:  the R code for  data generation of example 3.


2.  The folder "Simulation" contains three R files:

2.1.  sim1_main.R:  main procedure for example 1.

2.2.  sim2_main.R:  main procedure for example 2.

2.3.  sim3_main.R:  main procedure for example 3.


