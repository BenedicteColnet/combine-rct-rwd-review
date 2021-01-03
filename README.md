# Combine RCT RWD review

This repository contains R notebooks associated with the (submitted) review article: [_Causal inference methods for combining randomized controlled trials and observational studies: a review_](https://arxiv.org/abs/2011.08047).

Four notebooks are available, one for the simulation section and three for the data analysis section.

## Simulation section

- [`simulations-review-paper-samsi-wo-cw.Rmd`](https://gitlab.inria.fr/misscausal/combine-rct-rwd-review/-/blob/master/simulations/simulations-review-paper-samsi-wo-cw.Rmd) reproducing the **simulations and plots** presented in the article.

It needs to be launched with the [`estimators_and_simulations-wo-cw.R`](https://gitlab.inria.fr/misscausal/combine-rct-rwd-review/-/blob/master/simulations/estimators_and_simulations-wo-cw.R) code that contains the functions (IPSW, g-formula, AIPSW, and simulation protocol). 

This notebook does not contain the calibration weighting (CW) method as the function used in the paper is the one from [Lin Dong](https://lynndung.github.io/about/) as an implementations of her [research work](https://arxiv.org/abs/2003.01242). Because of the reviewing process we only publish the rest of the simulations and results.


## Data analysis section

Note that for reproducing the results from this section, an access to the medical data extracted from the [Traumabase registry](http://www.traumabase.eu/en_US) and the [CRASH-3 trial](https://crash3.lshtm.ac.uk/) are necessary and these are not publicly available.

- [`preprocess.Rmd`](https://gitlab.inria.fr/misscausal/combine-rct-rwd-review/-/blob/master/crash3-and-traumabase/preprocess.Rmd) performs the data preprocessing for the joint analysis of CRASH-3 and the Traumabase. It takes as an entry the raw data from each data sets and bind them with proper covariates. The output is the combined data with the raw Traumabase data (with missing values kept). Another similar data frame but with the imputed Traumabase is also produced.

- [`overlap_analysis.Rmd`](https://gitlab.inria.fr/misscausal/combine-rct-rwd-review/-/blob/master/crash3-and-traumabase/overlap_analysis.Rmd) proposes an analysis of the distributional shift between the CRASH-3 and Traumabase patients. The input is the merged table of both the randomized controlled trial and the observational study (corresponding to the output of `preprocess.Rmd`).

- [`combined_analysis_crash3_traumabase.Rmd`](https://gitlab.inria.fr/misscausal/combine-rct-rwd-review/-/blob/master/crash3-and-traumabase/combined_analysis_crash3_traumabase.Rmd) performs the average treatment effect estimation on the preprocessed data for the joint analysis of CRASH-3 and the Traumabase. The input is the merged table of both the randomized controlled trial and the observational study (corresponding to the output of `preprocess.Rmd`).

The key functions to perform the analyses in the latter two notebooks come from the script [`estimators.R`](https://gitlab.inria.fr/misscausal/combine-rct-rwd-review/-/blob/master/crash3-and-traumabase/estimators.R).
