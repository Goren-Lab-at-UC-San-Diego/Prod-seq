# Prod-Seq
This repository hosts the source code for the data processing and analysis model of Prod-seq. Please see below for the installation instructions, basic usage, and the documentation of the importable functions.

## Installation
1. Create a conda environment using the prodseq_analysis_environment.yml file (note: this step may take some time for conda to resolve the environment); conda will create an environment named “prodseq”
```
conda env create -f prodseq_analysis_environment.yml
```
**Note:** for users with experience with conda, the packages listed in the .yml file can also be installed manually (the packages are all commonly-used packages)
**Note:** if a different environment name is desired, the user can change the first line of the .yml file to “name: desired_environment_name”

