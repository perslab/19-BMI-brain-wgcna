# Code to reproduce WGCNA for BMI brain project

## How to run

1. clone the repo
2. download the expression data #TODO
3. install the required R packages
    1. Install Seurat 2.3.3 by following the instructions [here](https://satijalab.org/seurat/install.html) or downloading the source from github and installing from within an R session using `install.packages(pkgs="<path_to_seurat_source>", type=source, repos=NULL)`
    2. Install Bioconductor packages [Biobase](https://www.bioconductor.org/packages/release/bioc/html/Biobase.html) and [STRINGdb](https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html)
    3. Install CRAN packages by installing [renv](https://rstudio.github.io/renv/articles/renv.html) and running either `renv::init()` or `renv::restore(lockfile="./renv.lock")`. Alternatively install the packages manually: [here](https://cran.r-project.org/web/packages/here/index.html), [optparse](https://cran.r-project.org/web/packages/optparse/index.html), [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html), [Matrix](https://cran.r-project.org/web/packages/Matrix/index.html), [parallel](https://www.rdocumentation.org/packages/parallel/versions/3.6.1), [reshape](https://cran.r-project.org/web/packages/reshape/index.html), [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html), [liger](https://cran.r-project.org/web/packages/liger/index.html), [boot](https://cran.r-project.org/web/packages/boot/index.html), [R.utils](https://cran.r-project.org/web/packages/R.utils/index.html), [readr](https://cran.r-project.org/web/packages/readr/index.html).
4. Edit the parameters in `run_wgcna.sh` as desired (run `Rscript ./wgcna.R --help` for information on parameters) and run the analysis using `bash run_wgcna.sh` 
