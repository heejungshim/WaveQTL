## `get_Wmat.R' shows how to produce Discrete Wavelet Transform (DWT) matrix, which is required in effect size estimation in data space. 
##
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

## Haar DWT for region of 512 bp
InputD = diag(512)
res = FWT(InputD)
Wmat_512 = t(res$WCs)

## Haar DWT for region of 1024 bp
InputD = diag(1024)
res = FWT(InputD)
Wmat_1024 = t(res$WCs)

## Haar DWT for region of 2048 bp
InputD = diag(2048)
res = FWT(InputD)
Wmat_2048 = t(res$WCs)

## for different DWT for region of 512 bp
InputD = diag(512)
res = FWT(InputD, filter.number = 4, family ="DaubExPhase")
Wmat_512.v2 = t(res$WCs)


