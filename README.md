# CTC COAD READ

A [workflowr][] project.

[workflowr]: https://github.com/workflowr/workflowr


## Data pre-processing

Raw data is available at Gene Expression Omnibus (GEO, National Center for Biotechnology Information) under accession codes [GSE347686](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE347686)


##  Data and code availability

To reproduce our analysis, first clone source code from the [GitHub repository](https://github.com/TheAcetoLab/nuesch-ctc-coad-read). This repository is also archived at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXXXXX.svg)](https://doi.org/10.5281/XXXXXXXXXX)

    git clone https://github.com/TheAcetoLab/asawa-cto-pdac

Next, download processed data and output data deposited in [Zenodo](https://zenodo.org/records/22841801) into the cloned project folder and untar the files.

    for file in *.tar.gz; do tar xzvf "${file}" && rm "${file}"; done


## Reproducibility

The results form our analyses are listed below in webpage format. They were generated from R Markdown documents deposited in the [GitHub repository](https://github.com/TheAcetoLab/nuesch-ctc-coad-read). The workflow of the analysis was created using the [workflowr](https://cran.r-project.org/web/packages/workflowr/index.html) R package and can be reproduced in its totality using [workflowr](https://cran.r-project.org/web/packages/workflowr/index.html) [wflow_build](https://jdblischak.github.io/workflowrBeta/reference/wflow_build.html) command after the installation of the proper R packages. Session info, including R and package versions, was automatically included at the end of each analysis file. 

Files containing necessary data and pre-computed results from differential expression or gene-set enrichment analyses were deposited in [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22841801.svg)](https://zenodo.org/records/22841801). In order to generate those files again change the option `eval = FALSE` to `eval = TRUE` in the specific code chunk from the R Markdown file. The expected run time using the data deposited in the zenodo directoy  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21832574.svg)](https://doi.org/10.5281/zenodo.21832574) is approximately 30 minutes in a standard desktop environment.
