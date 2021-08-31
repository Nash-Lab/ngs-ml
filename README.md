# ngs-ml
Next Generation Sequencing and Machine Learning

Using [minimap2](https://github.com/lh3/minimap2) and [bbmap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/).

## Requirements
1. [minimap2](https://github.com/lh3/minimap2)
1. [bbmap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/)
1. C compiler such as [gcc](https://gcc.gnu.org/)
1. [Python](https://www.python.org/) ([Miniconda](https://docs.conda.io/en/latest/miniconda.html) is recommended)

## Workflow

# Summary
(Pacbio) **Minimap2** -> `ppba` -> **sort** -> `mlut` -| (Illumina) **BBMap** -> `pib` -| (From `mlut` and `pib`) `rib` -> **grep** ^1 -> **sort** -> **uniq** -c -> **sed** -E 's/^ *//; s/ /\t/' -|

# Found in this repository
script | input file(s) | output file(s) | comment
--- | --- | --- | ---
`ppba.c` | *aln.sam*, *ref.fa* | prepre_lut.tsv | Extract barcode, coressponding mutation(s) and other properties from Minimap2 alignment file.
`mlut.c` | *pre_lut.tsv* | *lut.tsv* | Filter and eliminate barcodes to give final lookup table. 
`lut_filter_bc.ipynb` | *lut.tsv* | *p_lut.csv* | Filter lookup table barcodes for only specific region
`lut_stats.ipynb` | *lut.tsv* | - | Visualize look up table
`pib.c` | *Bin#.sam* | *Bin#.fq* | 
`rib.c` | *Bin#.fq*, (*p_lut.tsv*) | *Bin#.tsv* | *Read Illumina barcodes* from file after cutadapt. Use `rib --help` for more information.
`rib_summary.ipynb` | *t1sct_Bin#.tsv* | - | Summary of the RIB outputs
`rib_heatmap.ipynb` | *t1sct_Bin#.tsv* | - | Heatmap of a RIB output
`ill_tag1_bins.ipynb` | *t1sct_Bin#.tsv* | *Bins.tsv* | Merges all ouputs from RIB into *Bins.tsv*
`ill_calc_props.ipynb` | *Bins.tsv* | *p_Bins.tsv* | Calculate properties for all mutations.

# Typical procedure
```bash
# Pacbio to look up table
minimap2 --cs -ax map-hifi ref.fa pacbio.fq > aln.sam
ppba ref.fa aln.sam > prepre_lut.tsv
sort prepre_lut.tsv > pre_lut.tsv
mlut pre_lut.tsv > m_lut.tsv
sed 's/"/!/g' m_lut.tsv > lut.tsv #Avoid bug in Python pandas

# Illumina to barcode fastq
bbmap.sh in=ill1.fq ref=ref.fa out=Bin1.sam
pib Bin1.sam > Bin1.fq

# Frequency analysis
rib -t lut.tsv -l 102 -q 0 Bin1.fq  > Bin1.tsv
grep ^1 Bin1.tsv > t1_Bin1.tsv
sort t1_Bin1.tsv > t1s_Bin1.tsv
uniq -c t1s_Bin1.tsv > t1sc_Bin1.tsv
sed -E 's/^ *//; s/ /\t/' t1sc_Bin1.tsv > t1sct_Bin1.tsv

# Use ill_tag1_bins.ipynb to get Bins.tsv
# Use ill_calc_props.ipynb to get p_Bins.tsv
```

## Run Notebooks with Conda
1. Download miniconda from https://docs.conda.io/en/latest/miniconda.html
1. Open an Anaconda Prompt and go to the ngs-ml folder (`cd <PATH>`)
1. Type `conda env create -f environment.yml -p ./env` (This may take a while)
1. Actiavte the enviroment with `conda activate ./env` inside the ngs-ml folder
1. Type `jupyter notebook` and open the notebooks!

## Compile *.c files, like rib.c
1. Make sure that `gcc` is installed with `gcc --version` in a a shell.
1. In `ngs-ml/code`, type `gcc <FILE.c> -o FILE`.
1. Change permission of the created binary file with `chmod +x <FILE>`.
1. Run the program with `./<FILE>`.

Note: For certain files, the `-lm` flag should be added as well, like: `gcc ppba.c -lm -o ppba`
