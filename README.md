# ngs-ml
Next Generation Sequencing and Machine Learning

Using [alignparse](https://jbloomlab.github.io/alignparse/) and [dms_variants](https://jbloomlab.github.io/dms_variants/)

## Workflow
`lt0_mini.ipynb` -> `lut_filter_bc.ipynb` -> `rib` -> `ill_tag1_bins.ipynb`

script | input file(s) | output file(s) | comment
--- | --- | --- | ---
`lt0_mini.ipynb` | .fastq, .gb, .yaml, genSeq | *con_lut.csv*, *var_lut.csv*, *rib_lut.csv* | Generate different look up tables
`lut_filter_bc.ipynb` | *rib_lut.csv* | *p_rib_lut.csv* | Filter lookup table barcodes for only specific region
`lut_stats.ipynb` | *rib_lut.csv*, (*con_lut.csv*, genSeq) | - | Visualize look up table
`rib.c` | .fastq, (*p_rib_lut.csv*) | *Bin#.csv* | *Read Illumina barcodes* from file after cutadapt. Use `rib --help` for more information.
`rib_summary.ipynb` | *Bin#.csv* | - | Summary of the RIB outputs
`rib_collapsed_summary.ipynb` | *Bin#.csv*, c_Bin.csv | - | Summary collapsed RIB outputs
`rib_heatmap.ipynb` | *Bin#.csv* | - | Heatmap of a RIB output
`ill_tag1_bins.ipynb` | *Bin#.csv* | *tag1_Bin#.csv*, *Bins.csv* | Extract usable barcodes from RIB output and merges them into *Bins.csv*
`ill_calc_props.ipynb` | *Bins.csv* | *p_Bins.csv* | Calculate properties for all mutations.
`daoiv` | *p_Binds.csv* | - | _**TODO**_ On [Streamlit](https://share.streamlit.io/aa-schoepfer/daoiv/main/daoiv.py), visualize results
 
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
