## 20221031
## S. Ning
## argo + sgl, national
## data prep
## input: path_ws, path_data, ili_version, GT_all
## output: state_data
# ready to run

## naming and creating directories
## data input
path_data <- file.path(path_ws, "data")
## ili
ili.folder <- file.path(path_data, ili_version)
population.file <- file.path(path_data,"Population.csv")
## GT
GT_version <- paste0(GT_all, "/")
gt.folder <- file.path(path_data, GT_version)
## output
out.folder <- paste0(ili.folder, GT_version)
if(!exists(out.folder)){
  dir.create(out.folder, showWarnings = FALSE, recursive = TRUE)
}  

## GT format issue after 2020
gt.parser.cur <- argo:::gt.parser.pub.api
## read in data
#setwd("/Users/sn9/Dropbox/Williams/Research/SGL_ARGO/ARGOX_shared2022")
state_data <- load_reg_data(gt.folder=gt.folder,
                            ili.folder=ili.folder,
                            gft.file= file.path(path_data, "GFT.txt"),
                            population.file = population.file,
                            gt.parser = gt.parser.cur)

print("Data input done.")