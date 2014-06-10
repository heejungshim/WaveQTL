## `combine_two_strands_with_correction_for_mappability.R' shows how to combine DNase-seq data from two strands while taking mappability into account as we did in Shim and Stephens (2014).
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

#######################################################
## Read DNase-seq data from two strands.
## You probably need to correct a path to the data.
#######################################################
path = "/mnt/lustre/home/shim/WaveQTL/data/dsQTL/DNase.chr17.10160989.10162012.dat"
DNase.dat = read.table(path)
dim(DNase.dat) # 70 by 2048; each column corresponds to each individual; the first (second) 1024 rows contain DNass-seq read count from +(-) strand in each positions; We already masked 5bp surrounding any SNP (i.e., the SNP position and 2bp on either side) to eliminate biases stemming from DNase I seqeucne preference (see the supplementary material of Degner et al 2012 for details). 

#######################################################
## Read mappability.
## You probably need to correct a path to the data.
#######################################################
path = "/mnt/lustre/home/shim/WaveQTL/data/dsQTL/DNase.mappability.chr17.10160989.10162012.dat"
map.dat = read.table(path)
dim(map.dat) # 1 by 2048; the first (second) 1024 rows indicates mappability from +(-) strand in each positions; `1' indicates uniquly mappabile base.

#############################################################
## combine two strands while taking mappability into account
#############################################################
numBPs = dim(DNase.dat)[2]/2
numINDs = dim(DNase.dat)[1]

# take mappability into account
map = rep(0, numBPs*2)
wh = (map.dat[1,] == 1)
map[wh] = 1
dat = matrix(data = 0, nr = numINDs, nc = numBPs*2)
dat[,wh] = as.matrix(DNase.dat[,wh])

# combine two strands
all.dat = dat[,1:numBPs] + dat[,(numBPs+1):(numBPs+numBPs)]
all.map = map[1:numBPs] + map[(numBPs+1):(numBPs+numBPs)]
pheno.dat = matrix(data = 0, nr = numINDs, nc = numBPs)
wh2 = which(all.map > 0)
pheno.dat[,wh2] = t(t(all.dat[,wh2])/all.map[wh2])

# we can check if "pheno.dat" is the same as "/mnt/lustre/home/shim/WaveQTL/data/dsQTL/chr17.10160989.10162012.pheno.dat"
A = as.matrix(read.table("/mnt/lustre/home/shim/WaveQTL/data/dsQTL/chr17.10160989.10162012.pheno.dat"))
sum(pheno.dat != A)
# 0

   
