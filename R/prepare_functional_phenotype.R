## `prepare_functional_phenotype.R' contains R scripts to show 1) reading DNase-seq data from files in hdf5, 2) reading mappability information from file in hdf5, 3) masking 5bp surrounding any SNP (i.e., the SNP position and 2bp on either side) to eliminate biases stemming from DNase I sequence preference, and 4) combining DNase-seq data from two strands while taking mappability into account as we did in Shim and Stephens (2014).
##
##
## Copyright (C) 2014 Heejung Shim and Ester Pantaleo (for the function get.counts.h5)
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

library(rhdf5)

#**************************************************************************
#  function to Get counts from hdf5
#  written by Ester Pantaleo
#**************************************************************************
get.counts.h5 <- function(list_path, chr, locus.start, locus.end, list_name = NULL){
    M <-  NULL
    for (h5file in list_path){
        print(paste0("Loading ", h5file))
        print(paste0("h5read(", h5file, ",", chr, ", index=list(", locus.start, ":", locus.end, ")"))
        M        <- rbind(M, h5read(h5file, chr, index=list(locus.start:locus.end)))
    }
    row.names(M) = list_name
    return(M)
}


#########################################
#### Users need to make changes in paths
#########################################

## Path to directory which contain DNase-seq data as hdf5 format, 
data.path = "/mnt/lustre/data/internal/genome_db/hg18/dnase/"
## Path to mappability information as hdf5 format
path.mapp  = "/mnt/lustre/data/internal/genome_db/hg18/mappability/roger_20bp_mapping_uniqueness.h5"

## Path to the list of individual IDs
inds.IDs = scan("~/WaveQTL/data/Shim_2014_etc/DNaseI.individuals.oneline.txt", what="")
numINDs = length(inds.IDs)



######################################
### specify information on the site
######################################
chrIX = "chr17"
locus.start = 10160989
locus.end = 10162013 - 1 
numBPs = 1024


########################
### 1. Read phenotype data
########################
DNase.hdf5 = matrix(data=NA, nr = numINDs, nc = numBPs*2)
for(i in 1:numINDs){
    
    path.fwd = paste0(data.path, "dnase_", inds.IDs[i], "_fwd.h5")
    DNase.hdf5[i, 1:numBPs] = as.matrix(get.counts.h5(path.fwd, chrIX, locus.start+1, locus.end+1))

    path.rev = paste0(data.path, "dnase_", inds.IDs[i], "_rev.h5")
    DNase.hdf5[i, ((1:numBPs)+numBPs)] = as.matrix(get.counts.h5(path.rev, chrIX, locus.start+1, locus.end+1))
    
}

dim(DNase.hdf5)
# 70 by 2048; each row corresponds to each individual; the first (second) 1024 columns contain DNass-seq read count from +(-) strand in each positions;



###############################
# 2. Read mappability information
###############################
map.hdf5 = matrix(data=NA, nr = 1, nc = numBPs*2)
map.hdf5[1, 1:numBPs] = as.matrix(get.counts.h5(path.mapp, chrIX, locus.start+1, locus.end+1))
map.hdf5[1, ((1:numBPs)+numBPs)] = as.matrix(get.counts.h5(path.mapp, chrIX, locus.start-20+2, locus.end-20+2))
dim(map.hdf5) # 1 by 2048; the first (second) 1024 rows indicates mappability from +(-) strand in each positions; `1' indicates uniquely mappable base.
map.dat = map.hdf5


###############################################################################
##  3. Mask 5bp surrounding any SNP (i.e., the SNP position and 2bp on either side)
##  to eliminate biases stemming from DNase I sequence preference
##  (see the supplementary material of Degner et al 2012 for details). 
###############################################################################

loc_info = rep(NA, numBPs*2)
loc_info[1:numBPs] = locus.start:locus.end
loc_info[(1:numBPs)+numBPs] = locus.start:locus.end

DNase.in = DNase.hdf5
      
## read all SNP information at the site
geno.path = "~/WaveQTL/data/dsQTL/chr17.10160989.10162012.2kb.cis.noMAFfilter.info"

if(file.info(geno.path)$size == 0){

    DNase.out = DNase.in

}else{

    geno = read.table(geno.path, as.is = TRUE)
    
    SNP_posi = as.numeric(geno[,6])
    del_posi = sort(unique(union(union(union(union(SNP_posi-2, SNP_posi -1), SNP_posi), SNP_posi+1), SNP_posi+2)))  
    wh_del = which((loc_info %in% del_posi)==TRUE)
  
    DNase.out = DNase.in
    
    if(length(wh_del) > 0){
        DNase.out[,wh_del] = matrix(data=0, nr = numINDs, nc = length(wh_del))
    }
}


DNase.dat = DNase.out

#############################################################
## 4. Combine two strands while taking mappability into account
#############################################################

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


#out.path = "~/WaveQTL/data/dsQTL/chr17.10160989.10162012.pheno.dat"
#write.table(pheno.dat, file=out.path, row.names=FALSE, col.names = FALSE,quote=FALSE)

