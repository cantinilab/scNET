# Evaluating the reproducibility of single-cell gene regulatory network inference algorithms
We here benchamrk three single-cell network inference algorithms based on their reproducibility, i.e. their ability to infer similar networks once applied to two independent datasets from the same biological condition. 

The benchmarked methods are:
* [PPCOR](https://cran.r-project.org/web/packages/ppcor/index.html)
* [GRNBoost2](https://github.com/aertslab/GRNBoost)
* [GENIE3](https://bioconductor.org/packages/release/bioc/html/GENIE3.html)

The methods are tested in three biological contexes:
1. human retina
2. colorectal cancer (CRC) T-cells
3. different cell types from hematopoiesis

Please note that the notebook executes by default the comparison in human retina. By uncommenting the lines relative to the biolgical context of interest in the cell corresponding to the data loading the user can however change the starting dataset.

## Input data
The input data are available in the `./data/` folder. Each biological context corresponds to a different subfolder.

## Install the software environment

* Install conda from https://docs.conda.io/en/latest/miniconda.html
*

## Run the notebook

* Enter the conda environment: `conda activate XXX`.
* Launch the notebook with `jupyter-notebook`.

##  Cite the work
The preprint describing momix is available in BioRxiv
xxxx
