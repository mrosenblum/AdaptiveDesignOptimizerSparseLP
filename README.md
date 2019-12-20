
<!-- README.md is generated from README.Rmd. Please edit that file -->
AdaptiveDesignOptimizerSparseLP
===============================

<!-- badges: start -->
[![Travis build status](https://travis-ci.com/mrosenblum/AdaptiveDesignOptimizerSparseLP.svg?branch=master)](https://travis-ci.com/mrosenblum/AdaptiveDesignOptimizerSparseLP) [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/mrosenblum/AdaptiveDesignOptimizerSparseLP?branch=master&svg=true)](https://ci.appveyor.com/project/mrosenblum/AdaptiveDesignOptimizerSparseLP) <!-- badges: end -->

The goal of AdaptiveDesignOptimizerSparseLP is to optimize the decision rule and multiple testing procedure for two-stage, two subpopulation, adaptive enrichment trial designs. The methodology used is available here: https://biostats.bepress.com/jhubiostat/paper273/

Installation
------------

The software requires that you have installed a linear program solver, preferably Cplex or Gurobi since they are faster at solving large linear programs (both of which can be obtained free by academics from IBM and Gurobi, respectively), but also possible is to install Matlab or GLPK. You need to ensure that the path to whichever of these you have installed is available to R. You can then select your solver when calling the main function (optimize_design) in this R package.

You can install the released version of AdaptiveDesignOptimizerSparseLP from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("AdaptiveDesignOptimizerSparseLP")
```

Example
-------

This is a basic example which shows you how to solve a common problem:

``` r
library(AdaptiveDesignOptimizerSparseLP)
## basic example code
```
