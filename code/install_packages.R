# default drive for package installation is too small
dir.create('.Rpackages')
.libPaths( c('.Rpackages',.libPaths() ) )

fromCRAN <- c(
   "Rcpp","openssl","CpGassoc", "rmarkdown"
  ,"knitr", "matrixStats","reshape","glmnet"
  ,"statmod","XML"
  ,"pryr", "qqman", "RPMM", "MASS", "sandwich", "lmtest","foreach","doParallel"
  ,"magrittr","purrr","svd","devtools","stringi","data.table"
  )

fromBioC <- c("minfi", "missMethyl", "ENmix","IlluminaHumanMethylation450kanno.ilmn12.hg19",
                      "IlluminaHumanMethylation450kmanifest", "IlluminaHumanMethylationEPICmanifest",
                      "sva", "IlluminaHumanMethylationEPICanno.ilm10b2.hg19","illuminaio", 
                      "DMRcate", "shinyMethyl","bumphunter","wateRmelon","FDb.InfiniumMethylation.hg19","DMRcatedata","ChAMP")

#' install these from CRAN:
toinstallCRAN <- setdiff(fromCRAN, installed.packages()[,1])
if(length(toinstallCRAN >= 1)) {
  install.packages(toinstallCRAN,dependencies=TRUE)
  cat("finished installing new packages from CRAN\n")
} else cat("packages we need from CRAN are already installed\n")


if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.12")

#' install these from BioConductor:
toinstallBioC <- setdiff(fromBioC, installed.packages()[,1])
if(length(toinstallBioC >= 1)) {
  BiocManager::install(toinstallBioC, update = FALSE)
  cat("finished installing new packages from BioConductor\n")
} else cat("packages we need from BioConductor are already installed\n")

devtools::install_github("hhhh5/ewastools@ebc2021")
devtools::install_github("sailalithabollepalli/EpiSmokEr",type = "source") # not on CRAN

#' check that we were successful
if(!all(c(toinstallBioC, toinstallCRAN) %in% installed.packages()[,1])) stop(
  "required packages not installed - please retry script carefully making sure you have already updated R and work through any error messages")

if(!as.numeric(sub("\\.[0-9]$", "", installed.packages()["minfi","Version"])) >= 1.24) stop(
  "you don't have the minfi version needed for this workshop")
