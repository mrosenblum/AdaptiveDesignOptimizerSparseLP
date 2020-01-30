---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# AdaptiveDesignOptimizerSparseLP

The goal of AdaptiveDesignOptimizerSparseLP is to construct an optimal two-Stage, two subpopulation, adaptive enrichment designs for a given problem. The problem inputs are the desired familywise Type I error rate and power, and set of scenarios (data generating distributions) of interest. The software aims to minimize the expected sample size under the power and Type I error constraints.

## Installation

You can install AdaptiveDesignOptimizerSparseLP using the remotes R package, which can be obtained by typing the following in your R session:

``` r
source("https://install-github.me/r-lib/remotes"
```

and then by entering the following in your R session: 

``` r
remotes::install_github("mrosenblum/AdaptiveDesignOptimizerSparseLP")
```

## Example

This is an example, which involves solving the problem from Example 3.2 as described in Section 5.2 of the manuscript. 


