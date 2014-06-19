## `get_effectSizeinDataSpace.R' contains R scripts to obtain estimated effect sizes in data space.
## 
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


## Read DWT matrix 
Wmat_1024 = read.table("../data/DWT/Wmat_1024",as.is = TRUE)
W2mat_1024 = Wmat_1024*Wmat_1024

## We'll look at effect size of 11th SNP in genotype file
sel_geno_IX = 11

## Read posterior mean and varaince of effect sizes in Wavelet space
## Set a path to files
beta_mean_path= "../test/dsQTL/output/test.no.QT.fph.mean.txt"
beta_var_path = "../test/dsQTL/output/test.no.QT.fph.var.txt"

## Read posterior mean in Wavelet space and transform them back to data space 
beta_mean = as.numeric(read.table(beta_mean_path)[sel_geno_IX,2:1025])
beta_dataS = as.vector(-matrix(data=beta_mean, nr = 1, nc = 1024)%*%as.matrix(Wmat_1024))
length(beta_dataS)
## [1] 1024
beta_dataS[1:6]
## [1] 1.352119e-11 1.352119e-11 1.352119e-11 1.352119e-11 1.352119e-11
## [6] 1.352119e-11



## Read posterior variance in Wavelet space, transform them back to data space, and get standard deviation
beta_var = as.numeric(read.table(beta_var_path)[sel_geno_IX,2:1025])
beta_var_dataS = as.vector(matrix(data=beta_var, nr = 1, nc = 1024)%*%as.matrix(W2mat_1024))
beta_sd_dataS = sqrt(beta_var_dataS)
length(beta_sd_dataS)
## [1] 1024
beta_sd_dataS[1:6]
## [1] 3.482833e-11 3.482833e-11 3.482833e-11 3.482833e-11 3.482833e-11
## [6] 3.482833e-11




## Visualize estimated effect size in the data space
ymin_beta = min(beta_dataS - 3*beta_sd_dataS) - abs(min(beta_dataS - 3*beta_sd_dataS))*0.0000000001
ymax_beta = max(beta_dataS + 3*beta_sd_dataS) + abs(max(beta_dataS + 3*beta_sd_dataS))*0.0000000001

beta_l = beta_dataS - 3*beta_sd_dataS
beta_r = beta_dataS + 3*beta_sd_dataS

wh_l = which(beta_l > 0)
wh_r = which(beta_r < 0)
high_wh = sort(unique(union(wh_l, wh_r)))

xval = 1:1024
col_posi = xval[high_wh]

pdf("../test/dsQTL/effectSize.pdf", width = 8, height=3)
par(mar = c(2,4,4,2))
plot(1,1,type="n", xlab = "position", ylab = "Effect size",ylim=c(ymin_beta, ymax_beta),xlim=c(1, 1024),main ="Posterior mean +/-3 posterior standard deviation", axes=FALSE)
axis(2)
if(length(col_posi) > 0){
for(j in 1:length(col_posi)){
	polygon(c(col_posi[j]-0.5, col_posi[j]-0.5, col_posi[j]+0.5, col_posi[j]+0.5), c(ymin_beta-2, ymax_beta+2, ymax_beta+2, ymin_beta-2), col ="pink", border = NA)
}
}

abline(h = 0, col = "red")
points(xval, beta_dataS, col = "blue", type="l")
points(xval, beta_l, col = "skyblue", type="l")
points(xval, beta_r, col = "skyblue", type="l")
box()

dev.off()


