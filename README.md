# Evaluating the reproducibility of single-cell gene regulatory network inference algorithms
We here benchamrk three single-cell network inference algorithms based on their reproducibility, i.e. their ability to infer similar networks once applied to two independent datasets from the same biological condition. 

The benchmarked methods are:
* [PPCOR](https://cran.r-project.org/web/packages/ppcor/index.html)
* [GRNBoost2](https://github.com/aertslab/GRNBoost)
* [GENIE3](https://bioconductor.org/packages/release/bioc/html/GENIE3.html)
* [CLR](https://rdrr.io/bioc/minet/man/clr.html)
* [GeneNet](https://cran.r-project.org/web/packages/GeneNet/index.html)


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

The preprocessed input data are available at https://cloud.biologie.ens.fr/index.php/s/JuJgrIL1jC6yZh4/download. Details on the preprocessing steps are provided in the methods of the paper. 

To access all data:
* Clone or download the scNET repository
* From R terminal or Rstudio, run the following lines

```
setwd('../scNET/')
dataURL= 'https://cloud.biologie.ens.fr/index.php/s/JuJgrIL1jC6yZh4/download'
download.file(dataURL, 'scNET_data.zip')
unzip('scNET_data.zip')
```

* In macOS environment, unzipping the data file from the terminal may be more efficient:

```
cd ~/scNET/
unzip scNET_dat.zip
```


## Install the software environment

* Install conda from https://docs.conda.io/en/latest/miniconda.html
* Create conda environment from yml in scNET repository by entering the following line in a terminal

```
cd scNET
conda env create -f scNET.yml
```

## Run the notebook

* Enter the conda environment: `conda activate scNET`.
* Launch the notebook with `jupyter-notebook`.

## Cite the work
The preprint describing momix is available in BioRxiv
https://www.biorxiv.org/content/10.1101/2020.11.10.375923v1

## Additional Networks
Users can analyze the reproducibility of networks produced by other algorithms using this workflow. 
To do so, save two networks inferred with independant datasets into the scNET Results folder.
Networks must be formatted into 3 columns (colum 1: gene1, column 2: gene2, column 3: interaction weight), in `.tsv` or tab seperated file format
Then, run notebook section `Algorithm reproducibility evaluation`to calculate metrics.


