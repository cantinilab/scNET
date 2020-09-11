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
As detailed in the paper, the data used for this benchmark are the following:
* Retina:[Menon, M. et al.](https://www.nature.com/articles/s41467-019-12780-8) [Lukowski SW. et al.](https://www.embopress.org/doi/10.15252/embj.2018100811)
* CRC T-cells: [Zhang et al.](https://www.nature.com/articles/s41597-019-0131-5) [Li et al.](https://www.nature.com/articles/ng.3818)
* Hematopoiesis:[Hay et al. ](https://www.sciencedirect.com/science/article/pii/S0301472X18308051?via%3Dihub) [Setty et al.](https://www.nature.com/articles/s41587-019-0068-4)

The preprocessed input data are available at XXX. Details on the preprocessing steps are rpovided in the methods of the paper. Download the two independent datasets of the biological context you want to test (retina, CRC T-cells, hematopoiesis) and save them in the `./data/` folder here provided. 

## Install the software environment

* Install conda from https://docs.conda.io/en/latest/miniconda.html
*

## Run the notebook

* Enter the conda environment: `conda activate XXX`.
* Launch the notebook with `jupyter-notebook`.

##  Cite the work
The preprint describing momix is available in BioRxiv
xxxx
