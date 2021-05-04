# ngs-ml
Next Generation Sequencing and Machine Learning

Using [alignparse](https://jbloomlab.github.io/alignparse/) and [dms_variants](https://jbloomlab.github.io/dms_variants/)

## Workflow
`lt0_mini.ipynb` -> `rib` -> `rib_summary.ipynb`

## Run Notebooks with Conda
1. Download miniconda from https://docs.conda.io/en/latest/miniconda.html
1. Open an Anaconda Prompt and go to the ngs-ml folder (`cd <PATH>`)
1. Type `conda env create -f environment.yml -p ./env` (This may take a while)
1. Actiavte the enviroment with `conda activate ./env` inside the ngs-ml folder
1. Type `jupyter notebook` and open the notebooks!

## Compile *.c files, like rib.c
1. Make sure that `g++` is installed with `g++ --version` in a a shell.
1. In `ngs-ml/code`, type `g++ <FILE.c> -o FILE`.
1. Change permission of the created binary file with `chmod +x <FILE>`.
1. Run the program with `./<FILE>`.
