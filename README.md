# Counterfactual-Simulation-Eye-tracking
## Project Overview
This study uses eye tracking data taken from Gerstenberg et al. (2017) to investigate how counterfactual simulation is employed during causal judgements. An autogressive logistic mixed effects model with a first-order lag term is used to explore how counterfactual gaze patterns vary across different judgement conditions: causal judgement, counterfactual judgement and outcome judgement. 

## Structure of Repository
This repository contains the following files: 

- script.R

This contains all code for reproducing the study's results. This includes a replication of the Gerstenberg et al. (2017) ANOVA model with post-hoc Tukey tests, an AR(1) logistic mixed effects model, multiple independent samples t-tests and several ggplots.
- sessioninfo.txt:

This documents the R session information at time of analysis, including R version and all package versions used, to facilitate reproducibility.

## Running the Code
### Requirements
The following R packages are required:
- tidyverse
- lme4
- lmerTest
- performance
- arm
- conflicted

Install any missing packages with:
install.packages(c("tidyverse", "lme4", "lmerTest", "performance", "arm", "conflicted"))

### Data
Data is taken from Gerstenberg et al. (2017) and can be downloaded from:
[TobiasGerstenberg (2017)](https://github.com/tobiasgerstenberg/eye_tracking_causality/tree/master). 

Once downloaded, place the RData folder in the same directory as script.R.

### Running
The script can be run from top to bottom in Rstudio. Note that the AR(1) logistic mixed effects model (step 5) may take several minutes to fit due to the size of the dataset.
