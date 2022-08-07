# The MSTCAR Package
## Introduction
### What is an MSTCAR model?

A Multivariate Spatiotemporal Conditional Auto-Regressive (MSTCAR) model is a Bayesian spatiotemporal model that generates spatially smoothed estimates across discrete regions. It is used to stabilize estimates in regions with low population counts, thereby bringing the inherent spatial structure of your event data to the forefront. It does this by taking information from various aspects of the data:
- The adjacency structure of your data. Neighboring regions consolidate information to stabilize estimates.
- Sociodemographic groups. The MSTCAR model also lends rate information from the sociodemographic groups included in the data.
- Time. The MSTCAR model helps to reveal trends in regions across time that may not be evident from the raw data.

MSTCAR is an R package is used to better understand binomial- or poisson-distributed spatial count data by calculating spatially smoothed estimates from a Gibbs sampler.

### Motivation and Features

It is often very difficult to use an MSTCAR model in practice because it is 1. Very difficult to build without prior knowledge of spatiotemporal Bayesian hierarchical models, and 2. The complexity of the model makes it very computationally intensive. Because of this, There are several key ideas taken into consideration when creating the MSTCAR package:

- Designed to simply import count/population data for discrete regions and receive estimates smoothed by region, sociodemographic groups, and time periods.
- Designed to be easy-to-use and accessible, even for people without experience with Bayesian hierarchical models.
- Designed to be *fast*! Runs approximately 10 times faster than equivalent models written only in R.

## Installation

To install the MSTCAR package, a few dependent packages are needed first. MSTCAR depends on the following:

##### Packages
- [Rcpp](https://cran.r-project.org/package=Rcpp): The most important package for MSTCAR. MSTCAR relies heavily on C++ code for efficiency. Rcpp is a package that allows integration of C++ code into R.
- [RcppArmadillo](https://cran.r-project.org/package=RcppArmadillo): This is an addon for Rcpp that facilitates matrix and array manipulation.
- [RcppDist](https://cran.r-project.org/package=RcppDist): This is an addon for Rcpp that facilitates sampling of multivariate distribution.
- [spdep](https://cran.r-project.org/package=spdep): A package related to spatial dependence. Contains several tools that help articulate the spatial relationships of shapefiles.
- [knitr](https://cran.r-project.org/package=knitr): Facilitates creation of package vignettes.
##### Programs
- [RTools](https://cran.r-project.org/bin/windows/Rtools/): C++ compiler for R. Necessary for Rcpp and its dependencies.

### Installation Instructions
The MSTCAR package can be easily installed with a few lines of R code. From the R console, make sure that all the necessary dependencies are installed by running the following lines:

```sh
# Install dependent packages
install.packages(c("spdep", "Rcpp", "RcppArmadillo", "RcppDist", "knitr"))
# Install RTools
install.packages("installr")
installr::install.Rtools()
```
Then, the package can be easily installed from GitHub:
```sh
install.packages("devtools")
devtools::install_github("dldelara/mstcar")
```
If this is your first time using the MSTCAR package, check the introductory vignette to learn how to use the package:
```sh
browseVignettes("mstcar")
```

### Updates
Updating the MSTCAR package is fairly simple. All that has to be done is remove the package from R, restart R, and re-install from GitHub:
```sh
remove.packages("mstcar")
# Restart R to allow package installation
devtools::install_github("dldelara/mstcar")
```

## Beta
This package is currently in the beta stage. All functions are currently public but will be made private once the package is more stable. In the meantime, each public function and dataset has documentation to help orient you to using the package.

## Thanks!
Thanks for checking out the MSTCAR package. I hope you find it useful in your spatiotemporal modeling! Please feel free to give feedback on bugs, ways to make the package more accessible to use, and features you'd like to see added to the package.
