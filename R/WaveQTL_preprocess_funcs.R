## `WaveQTL_preprocess_funcs.R' contains functions to preprocess functional data for WaveQTL
## software. 
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


require("wavethresh")


##' Filter low count WCs.
##'
##'
##' This function filters out some WCs that are computed based on very
##' low counts. This function considers WCs as low count if the total
##' counts used in their computation were less than or equal to `meanR.thresh'
##' per individual on average.
##' @param Data matrix (or a vector when N = 1) with N (# of samples) by T (# of bps in a region); This matrix contains original data before a wavelet transform; Here, T should be a power of 2.  
##' @param meanR.thresh If average reads across individuals <= meanR.thresh,
##' those WCs are filtered out.
##' @return filtered.WCs a vector of length T; t-th element indicates whether t-th WC in output from the function \code{\link{FWT}} is filtered (0) or not (1). 
fiter.WCs <- function(Data, meanR.thresh){

        if(is.vector(Data)){
            dim(Data) <- c(1, length(Data))
        }        
  	numWCs = dim(Data)[2]
  	J = log2(numWCs)

	Mean_R = rep(NA, numWCs)
        Mean_R[1] = mean(apply(Data, 1, sum))
        Mean_R[2] = Mean_R[1]

        posi = 3
        for(ss in 1:(J-1)){
                num_loc = 2^ss
                size_int = numWCs/num_loc
                st = (0:(num_loc-1))*size_int + 1
                en = st + size_int -1

                for(ll in 1:num_loc){
                        Mean_R[posi] = mean(apply(Data[,st[ll]:en[ll]], 1, sum))
                        posi = posi + 1
                }
        }

	filtered.WCs = rep(0, numWCs)
	wh = which(Mean_R > meanR.thresh)

        if(length(wh) > 0){
        	filtered.WCs[wh] = rep(1, length(wh))
        }

	return(filtered.WCs)
}




##' Perform a wavelet transform.
##'
##'
##' This function performs a wavelet transform using a \code{\link{wavethresh}} R package
##' and returns WCs in the order that corresponds to output from the function
##' \code{\link{fiter.WCs}}. For now, the function doesn't allow users to specify the level of wavelet
##' decomposition and uses the maximum level decomposition.
##' @param Data matrix (or a vector when N = 1) with N (# of samples) by T (# of bps in a region);
##' This matrix contains original data to be decomposed; Here, T should be a power of 2.
##' @param filter.number default=1; argument to the function \code{\link{wd}} in the R package
##' \code{\link{wavethresh}}; See their manual for details.
##' @param family default="DaubExPhase"; argument to the function \code{\link{wd}} in
##' the R package \code{\link{wavethresh}}; See their manual for details.
##' @return WCs a matrix with N (# of samples) by T (# of bps in a region); n-th row contains WCs
##' for n-th sample; WCs are ordered from low resolution WC to high resolution WC; For example,
##' with a Haar wavelet transform, the first column contains WC (precisely speaking, scaling
##' coefficient) that corresponds to sum of data in the region. The second column contains WC
##' that contrasts the data in the first half vs second half of the region. The last column
##' contains WC that contrasts the data in the (T-1)-th bp vs T-th bp.
FWT <- function(Data, filter.number=1, family="DaubExPhase"){

        if(is.vector(Data)){
            dim(Data) <- c(1, length(Data))
        }    
  	T = dim(Data)[2]
  	J = log2(T)
  	N = dim(Data)[1]

	dat_D = matrix(data=NA, nr = N, nc = (T - 1))
	dat_C = rep(NA, N)

	dat_W = matrix(data=NA, nr = N, nc = T)

	for(j in 1:N){
		each_WT	= wd(Data[j,], filter.number=filter.number ,family=family) 
		dat_D[j,] = each_WT$D
		dat_C[j] = each_WT$C[length(each_WT$C)]
	}

	dat_W[,1] = dat_C
	dat_W[,2] = -dat_D[,(T -1)]

	st_input = 3
	en_posi = T - 2
	for(k in 1:(J-1)){
		st_posi = en_posi - 2^k + 1
		en_input = st_input + 2^k - 1
		dat_W[,st_input:en_input] = -dat_D[,st_posi:en_posi]
		en_posi = st_posi - 1
		st_input = en_input + 1
	}

	return(list(WCs = dat_W))

}




##' Quantile-transform data to a standard normal distribution. 
##'
##'
##' This function quantile-transforms data to a standard normal distribution. It randomly assign
##' ranks for ties.
##' @param x a vector containing data to be quantile-transformed.
##' @return quantile-transformed data.
QT_randomTie <- function(x) {

	x.rank = rank(x, ties.method="random")
	return(qqnorm(x.rank,plot.it = F)$x)
}


##' Correct for covariates.
##'
##' 
##' This function corrects for covariates.
##' @param x a vector of length N (# of samples) containing data. 
##' @param Covariates a matrix (or a vector if M = 1) with N by M (# of covariates)
##' containing covariates to correct for.
##' @return a vector of length N; covariates corrected data. 
corrected_forCovariates <- function(x, Covariates){
	return(lm(x~Covariates)$residuals)
}



##' Normalize WCs. 
##'
##'
##' This function quantile-transforms WCs to a standard normal distribution.
##' If covariates are provided, it corrects the quantile-transformed WCs for the covariates
##' and quantile-transforms the covariates-corrected WCs to a standard normal distribution.
##' @param WCs a matrix with N (# of samples) by T (# of bps in a region or # of WCs);
##' n-th row contains WCs for n-th sample.
##' @param Covariates default = NULL; a matrix (or a vector if M = 1) with N by M
##' (# of covariates) containing covariates to correct for.
##' @return QT_WCs a matrix with N (# of samples) by T (# of bps in a region or # of WCs);
##' It contains normalized WCs (Quantile-transformed and covariate-corrected WCs).
Normalize.WCs <- function(WCs, Covariates=NULL){

	# QT to a standard normal distribution
	QT_dat = apply(WCs, 2, QT_randomTie)

	# correct for covariates and QT to a standard normal distribution. 
	if(!is.null(Covariates)){
		corrected_QT.dat = apply(QT_dat, 2, corrected_forCovariates, Covariates)
		QT_dat = apply(corrected_QT.dat, 2, QT_randomTie)
	}

	return(list(QT_WCs = QT_dat))

}



##' Generate Group information.  
##'
##'
##' This function generates information on which WCs share hyperparameter \pi
##' in the model described in Shim and Stephens 2014. This information is used
##' as an input in WaveQTL software. As default (group.scale=NULL), the function
##' outputs group information indicating that WCs in the same scale share \pi.
##' To put WCs from multiple scales (only for consecutive scales) in the same group,
##' instruction should be provided in group.scale (see below for details). 
##' @param numWCs a positive number; power of 2; total number of WCs. 
##' @param group.scale default = NULL; a vector of nonnegative numbers; length of
##' the vector is the number of groups; i-th element in the vector indicates
##' the first scale of WCs for i-th group; Suppose group.scale = c(0, 1, 4, 5)
##' and numWCs = 1024. Then, there are four groups. The first group contains WC in
##' the 0-th scale. The second group consists of WCs from the 1st to 3rd scales.
##' WCs in the 4th scale are in the third group and WCs from the 5th
##' to the last (10th) scales are in the fourth group. 
##' @return group a vector of positive numbers; length of the vector is the number
##' of groups; i-th element in the vector indicates the start position of WCs for i-th
##' group from a list of all WCs; Suppose group = c(1, 2, 9, 17) which can be obtained
##' by using group.scale = c(0, 1, 4, 5) and numWCs = 1024 as input. Then, there are
##' four groups of WCs. The first group has one WC that is in the first position from
##' the list of all WCs. The second group has 7 WCs that locate from the second position
##' to 8th position. WCs from 9th to 16th positions are in the third group and WCs
##' from 17th to the end (1024th) positions are in the fourth group. 
generate_Group <- function(numWCs, group.scale=NULL){
    
    J = log(numWCs, 2)
    
    if(is.null(group.scale)){
        group.scale = 0:J
    }

    IX = 2^(group.scale[-1]-1)
    group = c(1, (IX+1))

    return(group)
}



##' Preprocess functional data for a WaveQTL software. 
##'
##'
##' This function preprocesses functional data for a wavelet-based approach
##' implemented in a WaveQTL software. If library.read.depth is provided,
##' the function standardizes the functional data by the library read depth
##' to account for different read depths across individuals. Then, the function
##' decomposes the (standardized) functional data into wavelet coefficients (WCs)
##' using a \code{\link{wavethresh}} R package and normalizes the WCs.
##' If Covariates are provided, the function corrects the WCs for the Covariates
##' during the normalization. See the description of the function
##' \code{\link{Normalize.WCs}} for details of normalization. Users can specify
##' which type of wavelet transform should be applied by using filter.number and
##' family in input arguments. They are arguments to the function \code{\link{wd}}
##' in the R package \code{\link{wavethresh}}; See their manual for details.
##' For now, the function doesn't allow users to specify the level of wavelet
##' decomposition and uses the maximum level decomposition.
##' In addition to wavelet transform, this function filters out some WCs
##' that are computed based on very low counts. The function considers WCs
##' as low count if the total counts used in their computation were less than
##' or equal to `meanR.thresh' per individual on average. 
##' @param Data matrix (or a vector when N = 1) with N (# of samples) by T (# of bps in a region);
##' This matrix contains original functional data to be decomposed;
##' Here, T should be a power of 2.
##' @param library.read.depth default= NULL a vector of length N (# of samples);
##' i-th element contains library read depth for i-th sample.
##' @param Covariates default = NULL; a matrix (or a vector if M = 1)
##' with N by M (# of covariates) containing covariates to correct for.
##' @param meanR.thresh If average reads across individuals <= meanR.thresh,
##' those WCs are filtered out.
##' @param filter.number default=1; argument to the function \code{\link{wd}} in the R package
##' \code{\link{wavethresh}}; See their manual for details.
##' @param family default="DaubExPhase"; argument to the function \code{\link{wd}}
##' in the R package \code{\link{wavethresh}}; See their manual for details.
##' @param no.QT TRUE or FALSE; default=FALSE; if no.QT == FALSE, perform quantile transform during normalization (often for testing); if no.QT == TRUE, just correct WCs for covariates (often for effect size estimation).
##'
##' 
##' @return WCs a matrix with N (# of samples) by T (# of bps in a region); n-th row contains WCs
##' for n-th sample; WCs are ordered from low resolution WC to high resolution WC; For example,
##' with a Haar wavelet transform, the first column contains WC (precisely speaking, scaling
##' coefficient) that corresponds to sum of data in the region. The second column contains WC
##' that contrasts the data in the first half vs second half of the region. The last column
##' contains WC that contrasts the data in the (T-1)-th bp vs T-th bp.
##' @return filtered.WCs a vector of length T; t-th element indicates
##' whether t-th WC in output (WCs) filtered (0) or not (1). 
WaveQTL_preprocess <- function(Data, library.read.depth = NULL, Covariates = NULL, meanR.thresh = 2, no.QT = FALSE, filter.number=1, family="DaubExPhase"){

    
	if(is.vector(Data)){dim(Data)<- c(1,length(Data))} #change Data to matrix
  	if(nrow(Data)==1){Covariates = NULL} #if only one observation, don't correct for covariates

	if(!is.null(Covariates)){
		if(is.vector(Covariates)){dim(Covariates)<- c(1,length(Covariates))} #change C to matrix
	}



  	T = dim(Data)[2]
  	J = log2(T)
  	if(!isTRUE(all.equal(J,trunc(J)))){stop("Error: number of columns of Data must be power of 2")}
  	N = dim(Data)[1]
        
        
	### generate filtered.WCs
	if(!is.null(meanR.thresh)){
		filtered.WCs = fiter.WCs(Data, meanR.thresh)				
	}else{
		filtered.WCs = NULL
	}	

        
        ### corrected for read depth
        if(!is.null(library.read.depth)){
            DataC = Data/library.read.depth
        }else{
            DataC = Data
        }
        
	### Wavelet Transform
        WCs = FWT(DataC, filter.number=filter.number, family=family)$WCs

	
        if(!no.QT){ # Normalize phenotype for testing
            ### Normalize WCs    
            if(N > 1){
		WCs = Normalize.WCs(WCs, Covariates)
            }
            WCs = WCs$QT_WCs
        }else{  # Normalize for effect size estimation 
            ### correct for Covariates 
            if(!is.null(Covariates)){
		WCs = apply(WCs, 2, corrected_forCovariates, Covariates)
            }
        }
        
        return(list(WCs = WCs, filtered.WCs = filtered.WCs))
    } 






