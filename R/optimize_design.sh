#!/bin/bash

module load conda_R
R CMD BATCH optimize_design.R optimize_design.Rout
