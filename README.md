# ngs-ml
Next Generation Sequencing and Machine Learning

Using [minimap2](https://github.com/lh3/minimap2) and [bbmap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/).

## Workflow

# Old workflow:
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

# New workflow
(Pacbio) **Minimap2** -> `ppba` -> **sort** -> `mlut` -| (Illumina) **BBMap** -> `pib` -| (From `mlut` and `pib`) `rib` -> **grep ^1** -> **sort** -> **uniq -c** -> **sed -E 's/^ *//; s/ /\t/'** -|

script | input file(s) | output file(s) | comment
--- | --- | --- | ---
`ppba.c` | *aln.sam* | prepre_lut.tsv | Extract barcode, coressponding mutation(s) and other properties from Minimap2 alignment file.
`mlut.c` | pre_lut.tsv | lut.tsv | Filter and eliminate barcodes to give final lookup table. 
`lut_filter_bc.ipynb` | *lut.tsv* | *p_lut.csv* | Filter lookup table barcodes for only specific region
`lut_stats.ipynb` | *lut.tsv* | - | Visualize look up table
`rib.c` | .fastq, (*p_lut.tsv*) | *Bin#.tsv* | *Read Illumina barcodes* from file after cutadapt. Use `rib --help` for more information.
`rib_summary.ipynb` | *t1sct_Bin#.csv* | - | Summary of the RIB outputs
`rib_heatmap.ipynb` | *t1sct_Bin#.csv* | - | Heatmap of a RIB output
`ill_tag1_bins.ipynb` | *t1sct_Bin#.csv* | *Bins.csv* | Extract usable barcodes from RIB output and merges them into *Bins.csv*
`ill_calc_props.ipynb` | *Bins.csv* | *p_Bins.csv* | _**TODO**_ Calculate properties for all mutations.

## Run Notebooks with Conda
1. Download miniconda from https://docs.conda.io/en/latest/miniconda.html
1. Open an Anaconda Prompt and go to the ngs-ml folder (`cd <PATH>`)
1. Type `conda env create -f environment.yml -p ./env` (This may take a while)
1. Actiavte the enviroment with `conda activate ./env` inside the ngs-ml folder
1. Type `jupyter notebook` and open the notebooks!

## Compile code from source
1. 

## Compile *.c files, like rib.c
1. Make sure that `gcc` is installed with `gcc --version` in a a shell.
1. In `ngs-ml/code`, type `gcc <FILE.c> -o FILE`.
1. Change permission of the created binary file with `chmod +x <FILE>`.
1. Run the program with `./<FILE>`.

Note: For certain files, the `-lm` flag should be added as well, like: `gcc ppba.c -lm -o ppba`
