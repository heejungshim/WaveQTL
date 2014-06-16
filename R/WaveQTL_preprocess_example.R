## `WaveQTL_preprocess_example.R' shows how to preprocess functional data and prepare
## input files for a software WaveQTL. 
## Copyright (C) 2014 Heejung Shim
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.



## read functions for WaveQTL preprocess
source("WaveQTL_preprocess_funcs.R")

## set seed
set.seed(1)

## specify a path to example data which is shown in Figure 2 of Shim and Stephens 2014. 
data.path = "../data/dsQTL/"
output.path = "../test/dsQTL/"


## read functional data
pheno.dat = as.matrix(read.table(paste0(data.path, "chr17.10160989.10162012.pheno.dat")))
dim(pheno.dat)
#[1]   70 1024

## read library read depth
library.read.depth = scan(paste0(data.path, "library.read.depth.dat"))
length(library.read.depth)
#[1] 70

## read Covariates
Covariates = as.matrix(read.table(paste0(data.path, "PC4.dat")))
dim(Covariates)
#[1] 70 4

## preprocess functional phenotype
meanR.thresh = 2
res = WaveQTL_preprocess(Data = pheno.dat, library.read.depth=library.read.depth , Covariates = Covariates, meanR.thresh = meanR.thresh)
str(res)
#List of 2
# $ WCs         : num [1:70, 1:1024] -1.415 0.126 0.347 -1.323 0.674 ...
# $ filtered.WCs: num [1:1024] 1 1 1 1 1 1 1 1 0 0 ...




## save output as files
cat(res$filtered.WCs, file = paste0(output.path, "use.txt"))
write.table(res$WCs, file= paste0(output.path, "WCs.txt"), row.names=FALSE, col.names = FALSE, quote=FALSE)

## produce group information and save it as a file
group.info = generate_Group(dim(res$WCs)[2])
group.info
# [1]   1   2   3   5   9  17  33  65 129 257 513
cat(group.info, file = paste0(output.path, "group.txt"))


## for effect size estimation, we need WCs without QT.
set.seed(1)
res.noQT = WaveQTL_preprocess(Data = pheno.dat, library.read.depth=library.read.depth , Covariates = Covariates, meanR.thresh = meanR.thresh, no.QT = TRUE)

write.table(res.noQT$WCs, file= paste0(output.path, "WCs.no.QT.txt"), row.names=FALSE, col.names = FALSE, quote=FALSE)







