# bayes_cc
This repository contains the code for the paper "A Bayesian framework for case-cohort Cox regression: application to dietary epidemiology", which can be found at http://arxiv.org/abs/2007.12974. Our dataset comes from the EPIC-Norfolk study, as described in Section 4 of the manuscript; we use it to investigate the associations between individual saturated fatty acids and type 2 diabetes. The dataset is publicly available by request. To gain access to the data used in this paper, please send a request and proposal/justification for its use to interact-releases@mrc-epid.cam.ac.uk with reference to the title of our paper. The data is contained in a .csv file; in the code, we use the file name "mydata_up.csv". The results in Section 3 only use simulated data and can therefore be reproduced without the dataset.

# Code

## Abstract

The code provided in the following R script files will reproduce all the numbers and figures in the paper (instructions can be found in the "Instructions" section below):

Section 3 (can be run without data)

* sim.R

Section 4.3 (requires data)

* app.R (requires the provided "prop_list.rds" file)
* data_summ.R
* app_results.R

#### Version of primary software used

R version 3.6.1

#### Libraries and dependencies used by the code

* foreach (1.4.7)
* doParallel (1.0.15)
* CVST (0.2.2)
* MASS (7.3.51.5)
* doRNG (1.8.2)
* survival (3.1.8)
* MCMCpack (1.4.5)
* MBSP (1.0)
* denstrip (1.5.4)
* ggplot2 (3.3.2)
* gridExtra (2.3)
* reshape (0.8.8)

### Instructions

#### Section 3

The file "sim.R" runs the main simulation in Section 3 for a single set of parameter values and prints out the results in Table 1. It runs on 20 cores on a single node for about 20 minutes. To obtain the results for all sets of parameter values, the user could either run the jobs successively, or simultaneously on an array job.

#### Section 4.3

The file "app.R" runs the main analysis for this section. It loads the dataset file and "prop_list.rds" (which contains parameters used in the sampling algorithm). It runs on a single core for about 17 hours and saves the results onto a file called "app_samp.rds".

The file "data_summ.R" creates Figure 1 and the data summary section of Table 2. It requires the dataset file, and it can be run in a short amount of time on a standard desktop computer.

The file "app_results.R" creates Figure 2 and the "Analysis results" section of Table 2. It loads the results from "app.R", and it can be run in a short amount of time on a standard desktop computer.
