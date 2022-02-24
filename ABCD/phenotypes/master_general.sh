#!/bin/bash

# General scripts to run before any specific phenotypic analyses

module load R/4.0.3
working_dir=/export/home/acarrion/acarrion/projects/resources/datasets/ABCD/data/working_data/phenotypes/
scripts_dir=/export/home/acarrion/acarrion/projects/resources/datasets/ABCD/scripts/phenotypes/

cd ${scripts_dir}
Rscript define_subsets.R --verbose

## next run phenotype specific analyses: e.g. master_reading.sh