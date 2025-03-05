# MRegCovD
The minimum regularized covariance determinant estimator
This file contains the data and R code to replicate the results of the manuscript 
"The Regularized Minimum Covariance Determinant Estimator".

The main function is the function mrcd in the script MRCD.R. 

For the applications to the octane and murder rate data, the scripts to use are 
octane.R and murderrate.R.   

The other scripts "helperfunctionsMRCD.R", "optimalh.R", "r6pack.R" and "MRCDreg.R" 
are sourced  when needed by the main scripts "MRCD.R", "octane.R" and "murderrate.R".

We plan to include this code in an R package on CRAN. 