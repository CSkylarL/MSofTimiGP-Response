# Instructions for Data Download
Before you repeat the analysis, please read the following instructions carefully to create the data folder and download the data. Please note that all the data is for research purposes only.

The `data` folder should contains the following subfolders:
1. `TNBC_scRNA_GSE169246` folder: preprocessed single cell RNA-seq (scRNA-seq) data for TNBC patients, which is downloaded from [GSE169246](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169246) and preprocessed with [Fig1.script1_Preprocess_scRNAseq.R](../Fig1/code/Fig1.script1_Preprocess_scRNAseq.R). We provide the preprocessed data through the zenodo (Version v1.1.1.0-1) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12209783.svg)](https://doi.org/10.5281/zenodo.12209783). Please download and decompress the data to the `data/` folder, rename the downloaded folder to `TNBC_scRNA_GSE169246`. The path to the data should be `data/TNBC_scRNA_GSE169246/*.rda`. 

2. `bulk_transcriptomics` folder: preprocessed bulk transcriptomics data for pan-cancer patients. 

    There are three ways to download the data:

    - We provide the data source for each dataset in `supplementary table 4` of the manuscript, including the original publication, assession number, and the link to the data source. You can download or request the data from the original publication. And preprocess the data according to this [example](https://github.com/CSkylarL/TimiGP/blob/master/inst/extdata/build_TNBCaPD1_RNA_info.R).
    - We provide preprocessed GEO datasets in the zenodo (Version v1.1.1.0-2) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12209773.svg)](https://doi.org/10.5281/zenodo.12209773). These datasets are used in the manuscript and can be used as an example to repeat the analysis. 
    - We deposited all the preprocessed data in the zenodo (Version v1.1.0.0) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12208113.svg)](https://doi.org/10.5281/zenodo.12208113). These datasets have restricted access and can be only used for the replication of the analysis in the manuscript. Please contact the corresponding author for more information. No other usage than repeat Fig.2 is allowed wihout permission.

    After you download or preprocessed the datasets, please save the data or rename the folder as `bulk_transcriptomics`. The path to the data should be `data/bulk_transcriptomics/*.rda`. 

    In addition to this preprocessed data, this repository also contains the intermediate and result files generated with TimiGP-Response. If you want to explore the immune landscape assocaited with the response to immunotherapy, you can use the intermediate files in the [Fig2](../Fig2/result/) folder. Please refer to the [README](../README.md) for more information.

    Note: The data file name is a little different from what displayed in our manuscript. The [config file](../Fig2/code/config.csv) provide the mapping between the data file name ("oldname") and the dataset name in the manuscript ("newname"), and the location to the TimiGP-response result folder.

3. `TNBC_IMC_zenodo.7990870` folder: preprocessed image mass cytometry (IMC) data for TNBC patients, which is downloaded from [zenodo](https://doi.org/10.5281/zenodo.7990870). Please request the data from the original publication. Once you have the data, please rename the downloaded folder as `TNBC_IMC_zenodo.7990870` and save it to the `data` folder. The path to the data should be `data/TNBC_IMC_zenodo.7990870/`.
  
   