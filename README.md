Bayesian adaptive Lasso for detecting item-trait relationship and DIF in MIRT models

We use a unified framework for detecting item-trait relationship and differential item functioning (DIF) in multidimensional item response theory (MIRT) models. By incorporating DIF effects in MIRT models, these problems can be considered as variable selection for latent/observed variables and their interactions. All codes are written in R for Bayesian estimation.

Uniform-1000SP20.R implements our Bayesian adaptive Lasso procedure for the N=1000, 20% DIF and small DIF condition in simulation study 1.

Nonuniform-1000SP20.R implements our Bayesian adaptive Lasso procedure for the N=1000, 20% DIF and small DIF condition in simulation study 2.

Uniform-example-CFA+DIF.R uses a confirmatory simple-structure MIRT model for DIF detection in the heuristic example. 

Uniform-example-EML1+DIF.R first identifies the item-trait structure by the EML1 method and then uses the structure as confirmatory for DIF detection in the heuristic example. 

Uniform-example-our.R uses our Bayesian adaptive Lasso procedure for simultaneously detecting item-trait relationship and DIF in the heuristic example. 

Uniform-realdata.R implements real data analysis using our Bayesian adaptive Lasso procedure for the uniform DIF condition.
