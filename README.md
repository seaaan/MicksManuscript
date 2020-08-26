README
================

This document describes how to reproduce the data analysis and figures.

All data go into subfolders of this folder.

## Obtaining the microarray/RNAseq data

Note to reviewers: During review, our data is private. To access our
data, view “data\_access\_for\_reviewers.xlsx” in the same folder as
this file. There is a password for the microarray data. There are
private links for our ELISA, MSD, and ddPCR data. Otherwise all of the
instructions below apply.

First, download four Supplementary files from our data from GEO at
<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142058>. There
are two raw data files and two metadata files:
`GSE142058_GenomeStudio_output_batch_X.txt.gz` and
`GSE142058_Metadata_batch_X.csv.gz`. Place all four files into
`data/raw/`.

Second, download the microarray/RNAseq data from the other studies:

  - From <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80455>,
    download `GSE80455_non_normalized.txt.gz`.
  - From <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86491>,
    download `GSE86491_totalRNA_tissue.counts.txt.gz`.
  - From <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122248>,
    download `GSE122248_counts.matrix.txt.gz`.

Place each of the three above files into subfolders of
`data/raw/other_studies/` named `GSEXXXX/`. For example,
`data/raw/other_studies/GSE80455/GSE80455_non_normalized.txt.gz`.

  - From <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29981>,
    download `GSE29981_RAW.tar`.
  - From <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28044>,
    download `GSE28044_RAW.tar`.
  - From <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6364>,
    download `GSE6364_RAW.tar`.

Untar each of the three above folders and place the resulting files into
folders named `data/raw/other_studies/GSEXXXX/`.

Data from one study (<https://www.ncbi.nlm.nih.gov/pubmed/26750085>) is
not available on GEO and can be obtained by contacting the authors.

## Obtaining the other data

There are several other data types, which are all available in a
figshare collection <https://doi.org/10.6084/m9.figshare.c.4820469>

Download each of the following **raw data** files:

  - `ELISA_Granulysin/ELISA_Granulysin_raw.zip`
  - `ddpcr_Granulysin_YWHAZ/ddpcr_Granulysin_YWHAZ_raw.zip`
  - `MSD/MSD_raw.zip`

Unzip the resulting files. Place all of the unzipped folders into
`data/raw`. For example, you will have:

  - `data/raw/ELISA_Granulysin_raw/ELISA_Granulysin_Plate1.csv`
  - `data/raw/ELISA_Granulysin_raw/ELISA_Granulysin_Plate2.csv`

## Running the analysis

Open `R/source_all.R` and install all of the packages listed there.

Open `make_paper.R`. Assuming that you have placed all of the data files
listed above into into the directories as described, the file paths in
this file will be correct. If you have placed them elsewhere, adjust the
file paths. This assumes that the working directory is set to the top
level of this folder. To reproduce the analysis, source `make_paper.R`.
This will run code to:

  - Load and clean all of the data files
  - Perform statistical analysis
  - Generate all figures
  - Generate all statistics reported in the manuscript
  - Generate the supplementary files

It will take a while to run.

Questions? <smhughes@uw.edu>

## Session info

The analysis for the paper was run using the following versions of R and
the packages:

    ## - Session info -----------------------------------------------------------------------------------
    ##  setting  value                       
    ##  version  R version 4.0.0 (2020-04-24)
    ##  os       Windows 10 x64              
    ##  system   x86_64, mingw32             
    ##  ui       RTerm                       
    ##  language (EN)                        
    ##  collate  English_United States.1252  
    ##  ctype    English_United States.1252  
    ##  tz       America/Los_Angeles         
    ##  date     2020-08-26                  
    ## 
    ## - Packages ---------------------------------------------------------------------------------------
    ##  package              * version    date       lib source        
    ##  abind                  1.4-5      2016-07-21 [1] CRAN (R 4.0.0)
    ##  affy                 * 1.66.0     2020-04-27 [1] Bioconductor  
    ##  affyio                 1.58.0     2020-04-27 [1] Bioconductor  
    ##  annotate               1.66.0     2020-04-27 [1] Bioconductor  
    ##  AnnotationDbi        * 1.50.0     2020-04-27 [1] Bioconductor  
    ##  askpass                1.1        2019-01-13 [1] CRAN (R 4.0.0)
    ##  assertthat             0.2.1      2019-03-21 [1] CRAN (R 4.0.0)
    ##  backports              1.1.6      2020-04-05 [1] CRAN (R 4.0.0)
    ##  base64                 2.0        2016-05-10 [1] CRAN (R 4.0.0)
    ##  beanplot               1.2        2014-09-19 [1] CRAN (R 4.0.0)
    ##  Biobase              * 2.48.0     2020-04-27 [1] Bioconductor  
    ##  BiocFileCache          1.12.0     2020-04-27 [1] Bioconductor  
    ##  BiocGenerics         * 0.34.0     2020-04-27 [1] Bioconductor  
    ##  BiocManager            1.30.10    2019-11-16 [1] CRAN (R 4.0.0)
    ##  BiocParallel           1.22.0     2020-04-27 [1] Bioconductor  
    ##  biomaRt                2.44.0     2020-04-27 [1] Bioconductor  
    ##  Biostrings             2.56.0     2020-04-27 [1] Bioconductor  
    ##  bit                    1.1-15.2   2020-02-10 [1] CRAN (R 4.0.0)
    ##  bit64                  0.9-7      2017-05-08 [1] CRAN (R 4.0.0)
    ##  bitops                 1.0-6      2013-08-17 [1] CRAN (R 4.0.0)
    ##  blob                   1.2.1      2020-01-20 [1] CRAN (R 4.0.0)
    ##  boot                   1.3-24     2019-12-20 [2] CRAN (R 4.0.0)
    ##  broom                  0.5.6      2020-04-20 [1] CRAN (R 4.0.0)
    ##  bumphunter             1.30.0     2020-04-27 [1] Bioconductor  
    ##  car                    3.0-7      2020-03-11 [1] CRAN (R 4.0.0)
    ##  carData                3.0-3      2019-11-16 [1] CRAN (R 4.0.0)
    ##  cellranger             1.1.0      2016-07-27 [1] CRAN (R 4.0.0)
    ##  cli                    2.0.2      2020-02-28 [1] CRAN (R 4.0.0)
    ##  codetools              0.2-16     2018-12-24 [2] CRAN (R 4.0.0)
    ##  colorspace             1.4-1      2019-03-18 [1] CRAN (R 4.0.0)
    ##  conflicted           * 1.0.4      2019-06-21 [1] CRAN (R 4.0.0)
    ##  crayon                 1.3.4      2017-09-16 [1] CRAN (R 4.0.0)
    ##  curl                   4.3        2019-12-02 [1] CRAN (R 4.0.0)
    ##  data.table             1.12.8     2019-12-09 [1] CRAN (R 4.0.0)
    ##  DBI                    1.1.0      2019-12-15 [1] CRAN (R 4.0.0)
    ##  dbplyr                 1.4.3      2020-04-19 [1] CRAN (R 4.0.0)
    ##  DelayedArray           0.14.0     2020-04-27 [1] Bioconductor  
    ##  DelayedMatrixStats     1.10.0     2020-04-27 [1] Bioconductor  
    ##  digest                 0.6.25     2020-02-23 [1] CRAN (R 4.0.0)
    ##  doRNG                  1.8.2      2020-01-27 [1] CRAN (R 4.0.0)
    ##  dplyr                * 0.8.5      2020-03-07 [1] CRAN (R 4.0.0)
    ##  drc                  * 3.0-1      2016-08-30 [1] CRAN (R 4.0.0)
    ##  edgeR                * 3.30.0     2020-04-27 [1] Bioconductor  
    ##  ellipsis               0.3.0      2019-09-20 [1] CRAN (R 4.0.0)
    ##  evaluate               0.14       2019-05-28 [1] CRAN (R 4.0.0)
    ##  fansi                  0.4.1      2020-01-08 [1] CRAN (R 4.0.0)
    ##  forcats              * 0.5.0      2020-03-01 [1] CRAN (R 4.0.0)
    ##  foreach                1.5.0      2020-03-30 [1] CRAN (R 4.0.0)
    ##  foreign                0.8-78     2020-04-13 [2] CRAN (R 4.0.0)
    ##  fs                     1.4.1      2020-04-04 [1] CRAN (R 4.0.0)
    ##  genefilter             1.70.0     2020-04-27 [1] Bioconductor  
    ##  generics               0.0.2      2018-11-29 [1] CRAN (R 4.0.0)
    ##  GenomeInfoDb           1.24.0     2020-04-27 [1] Bioconductor  
    ##  GenomeInfoDbData       1.2.3      2020-05-12 [1] Bioconductor  
    ##  GenomicAlignments      1.24.0     2020-04-27 [1] Bioconductor  
    ##  GenomicFeatures        1.40.0     2020-04-27 [1] Bioconductor  
    ##  GenomicRanges          1.40.0     2020-04-27 [1] Bioconductor  
    ##  GEOquery             * 2.56.0     2020-04-27 [1] Bioconductor  
    ##  ggplot2              * 3.3.0      2020-03-05 [1] CRAN (R 4.0.0)
    ##  ggrepel              * 0.8.2      2020-03-08 [1] CRAN (R 4.0.0)
    ##  glue                   1.4.0      2020-04-03 [1] CRAN (R 4.0.0)
    ##  gtable                 0.3.0      2019-03-25 [1] CRAN (R 4.0.0)
    ##  gtools                 3.8.2      2020-03-31 [1] CRAN (R 4.0.0)
    ##  haven                  2.2.0      2019-11-08 [1] CRAN (R 4.0.0)
    ##  HDF5Array              1.16.0     2020-04-27 [1] Bioconductor  
    ##  hgu133plus2.db       * 3.2.3      2020-05-12 [1] Bioconductor  
    ##  hms                    0.5.3      2020-01-08 [1] CRAN (R 4.0.0)
    ##  htmltools              0.4.0      2019-10-04 [1] CRAN (R 4.0.0)
    ##  httr                   1.4.1      2019-08-05 [1] CRAN (R 4.0.0)
    ##  illuminaHumanv4.db   * 1.26.0     2020-05-12 [1] Bioconductor  
    ##  illuminaio             0.30.0     2020-04-27 [1] Bioconductor  
    ##  IRanges              * 2.22.1     2020-04-28 [1] Bioconductor  
    ##  iterators              1.0.12     2019-07-26 [1] CRAN (R 4.0.0)
    ##  jsonlite               1.6.1      2020-02-02 [1] CRAN (R 4.0.0)
    ##  KernSmooth             2.23-16    2019-10-15 [2] CRAN (R 4.0.0)
    ##  knitr                  1.28       2020-02-06 [1] CRAN (R 4.0.0)
    ##  lattice                0.20-41    2020-04-02 [2] CRAN (R 4.0.0)
    ##  lifecycle              0.2.0      2020-03-06 [1] CRAN (R 4.0.0)
    ##  limma                * 3.44.1     2020-04-28 [1] Bioconductor  
    ##  lme4                 * 1.1-23     2020-04-07 [1] CRAN (R 4.0.0)
    ##  lmerTest             * 3.1-2      2020-04-08 [1] CRAN (R 4.0.0)
    ##  locfit                 1.5-9.4    2020-03-25 [1] CRAN (R 4.0.0)
    ##  lubridate            * 1.7.8      2020-04-06 [1] CRAN (R 4.0.0)
    ##  lumi                 * 2.40.0     2020-04-27 [1] Bioconductor  
    ##  magrittr               1.5        2014-11-22 [1] CRAN (R 4.0.0)
    ##  MASS                 * 7.3-51.5   2019-12-20 [2] CRAN (R 4.0.0)
    ##  Matrix               * 1.2-18     2019-11-27 [2] CRAN (R 4.0.0)
    ##  matrixStats            0.56.0     2020-03-13 [1] CRAN (R 4.0.0)
    ##  mclust                 5.4.6      2020-04-11 [1] CRAN (R 4.0.0)
    ##  memoise                1.1.0      2017-04-21 [1] CRAN (R 4.0.0)
    ##  methylumi              2.34.0     2020-04-27 [1] Bioconductor  
    ##  mgcv                   1.8-31     2019-11-09 [2] CRAN (R 4.0.0)
    ##  minfi                  1.34.0     2020-04-27 [1] Bioconductor  
    ##  minqa                  1.2.4      2014-10-09 [1] CRAN (R 4.0.0)
    ##  modelr                 0.1.7      2020-04-30 [1] CRAN (R 4.0.0)
    ##  multcomp               1.4-13     2020-04-08 [1] CRAN (R 4.0.0)
    ##  multtest               2.44.0     2020-04-27 [1] Bioconductor  
    ##  munsell                0.5.0      2018-06-12 [1] CRAN (R 4.0.0)
    ##  mvtnorm                1.1-0      2020-02-24 [1] CRAN (R 4.0.0)
    ##  nleqslv                3.3.2      2018-05-17 [1] CRAN (R 4.0.0)
    ##  nlme                   3.1-147    2020-04-13 [2] CRAN (R 4.0.0)
    ##  nloptr                 1.2.2.1    2020-03-11 [1] CRAN (R 4.0.0)
    ##  nor1mix                1.3-0      2019-06-13 [1] CRAN (R 4.0.0)
    ##  numDeriv               2016.8-1.1 2019-06-06 [1] CRAN (R 4.0.0)
    ##  openssl                1.4.1      2019-07-18 [1] CRAN (R 4.0.0)
    ##  openxlsx               4.1.5      2020-05-06 [1] CRAN (R 4.0.0)
    ##  org.Hs.eg.db         * 3.11.1     2020-05-12 [1] Bioconductor  
    ##  pander               * 0.6.3      2018-11-06 [1] CRAN (R 4.0.0)
    ##  patchwork            * 1.0.0      2019-12-01 [1] CRAN (R 4.0.0)
    ##  pillar                 1.4.4      2020-05-05 [1] CRAN (R 4.0.0)
    ##  pkgconfig              2.0.3      2019-09-22 [1] CRAN (R 4.0.0)
    ##  plater               * 1.0.2      2020-03-24 [1] CRAN (R 4.0.0)
    ##  plotrix                3.7-8      2020-04-16 [1] CRAN (R 4.0.0)
    ##  plyr                   1.8.6      2020-03-03 [1] CRAN (R 4.0.0)
    ##  preprocessCore         1.50.0     2020-04-27 [1] Bioconductor  
    ##  prettyunits            1.1.1      2020-01-24 [1] CRAN (R 4.0.0)
    ##  progress               1.2.2      2019-05-16 [1] CRAN (R 4.0.0)
    ##  purrr                * 0.3.4      2020-04-17 [1] CRAN (R 4.0.0)
    ##  quadprog               1.5-8      2019-11-20 [1] CRAN (R 4.0.0)
    ##  R6                     2.4.1      2019-11-12 [1] CRAN (R 4.0.0)
    ##  rappdirs               0.3.1      2016-03-28 [1] CRAN (R 4.0.0)
    ##  RColorBrewer           1.1-2      2014-12-07 [1] CRAN (R 4.0.0)
    ##  Rcpp                   1.0.4.6    2020-04-09 [1] CRAN (R 4.0.0)
    ##  RCurl                  1.98-1.2   2020-04-18 [1] CRAN (R 4.0.0)
    ##  readr                * 1.3.1      2018-12-21 [1] CRAN (R 4.0.0)
    ##  readxl                 1.3.1      2019-03-13 [1] CRAN (R 4.0.0)
    ##  reprex                 0.3.0      2019-05-16 [1] CRAN (R 4.0.0)
    ##  reshape                0.8.8      2018-10-23 [1] CRAN (R 4.0.0)
    ##  rhdf5                  2.32.0     2020-04-27 [1] Bioconductor  
    ##  Rhdf5lib               1.10.0     2020-04-27 [1] Bioconductor  
    ##  rio                    0.5.16     2018-11-26 [1] CRAN (R 4.0.0)
    ##  rlang                  0.4.6      2020-05-02 [1] CRAN (R 4.0.0)
    ##  rmarkdown              2.1        2020-01-20 [1] CRAN (R 4.0.0)
    ##  rngtools               1.5        2020-01-23 [1] CRAN (R 4.0.0)
    ##  Rsamtools              2.4.0      2020-04-27 [1] Bioconductor  
    ##  RSQLite                2.2.0      2020-01-07 [1] CRAN (R 4.0.0)
    ##  rstudioapi             0.11       2020-02-07 [1] CRAN (R 4.0.0)
    ##  rtracklayer            1.48.0     2020-04-27 [1] Bioconductor  
    ##  rvest                  0.3.5      2019-11-08 [1] CRAN (R 4.0.0)
    ##  S4Vectors            * 0.26.0     2020-04-27 [1] Bioconductor  
    ##  sandwich               2.5-1      2019-04-06 [1] CRAN (R 4.0.0)
    ##  scales                 1.1.1      2020-05-11 [1] CRAN (R 4.0.0)
    ##  scrime                 1.3.5      2018-12-01 [1] CRAN (R 4.0.0)
    ##  sessioninfo            1.1.1      2018-11-05 [1] CRAN (R 4.0.0)
    ##  siggenes               1.62.0     2020-04-27 [1] Bioconductor  
    ##  statmod                1.4.34     2020-02-17 [1] CRAN (R 4.0.0)
    ##  stringi                1.4.6      2020-02-17 [1] CRAN (R 4.0.0)
    ##  stringr              * 1.4.0      2019-02-10 [1] CRAN (R 4.0.0)
    ##  SummarizedExperiment   1.18.1     2020-04-30 [1] Bioconductor  
    ##  survival               3.1-12     2020-04-10 [2] CRAN (R 4.0.0)
    ##  TH.data                1.0-10     2019-01-21 [1] CRAN (R 4.0.0)
    ##  tibble               * 3.0.1      2020-04-20 [1] CRAN (R 4.0.0)
    ##  tidyr                * 1.0.3      2020-05-07 [1] CRAN (R 4.0.0)
    ##  tidyselect             1.1.0      2020-05-11 [1] CRAN (R 4.0.0)
    ##  tidyverse            * 1.3.0      2019-11-21 [1] CRAN (R 4.0.0)
    ##  vctrs                  0.3.0      2020-05-11 [1] CRAN (R 4.0.0)
    ##  withr                  2.2.0      2020-04-20 [1] CRAN (R 4.0.0)
    ##  writexl              * 1.3        2020-05-05 [1] CRAN (R 4.0.0)
    ##  xfun                   0.13       2020-04-13 [1] CRAN (R 4.0.0)
    ##  XML                    3.99-0.3   2020-01-20 [1] CRAN (R 4.0.0)
    ##  xml2                   1.3.2      2020-04-23 [1] CRAN (R 4.0.0)
    ##  xtable                 1.8-4      2019-04-21 [1] CRAN (R 4.0.0)
    ##  XVector                0.28.0     2020-04-27 [1] Bioconductor  
    ##  yaml                   2.2.1      2020-02-01 [1] CRAN (R 4.0.0)
    ##  zip                    2.0.4      2019-09-01 [1] CRAN (R 4.0.0)
    ##  zlibbioc               1.34.0     2020-04-27 [1] Bioconductor  
    ##  zoo                    1.8-8      2020-05-02 [1] CRAN (R 4.0.0)
    ## 
    ## [1] C:/Users/Sean Hughes/Documents/R/win-library/4.0
    ## [2] C:/Program Files/R/R-4.0.0/library
