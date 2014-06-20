#! /usr/bin/env python

## `extract.py' is python script to extract genotypes and information for a subset of SNPs which are located in a given region.
##
## Usage: extract.py CHROMOSOME START END
##
## User need to change path to genotype files and informaton files. This script assumes that genotype files are in '/home/shim/wavelets/data/genotype/chrX.YRI.70.mean.genotype.txt' and information files are in '/home/shim/wavelets/data/genotype/chrX.YRI.snpdata.txt'. 
##
## for exampe, "./extract.py 17 10158989 10164012" creates `chr17.10158989.10164012.YRI.mean.genotype.txt' and `chr17.10158989.10164012.YRI.snpdata.txt' in the current working directory.
## 
## Copyright (C) 2014 Shyam Gopalakrishnan
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

import sys

if (len(sys.argv) != 4):
    print "Usage: extract.py CHROMOSOME START END"
    sys.exit(1)

chrom = sys.argv[1]
start = int(sys.argv[2]) 
end = int(sys.argv[3]) 
dirpre = '/home/shim/wavelets/data/genotype/'

ingeno = open(dirpre+'chr'+str(chrom)+'.YRI.70.mean.genotype.txt')
insnp = open(dirpre+'chr'+str(chrom)+'.YRI.snpdata.txt')
outsnp = open('chr'+str(chrom)+'.'+str(start)+'.'+str(end)+'.YRI.snpdata.txt', 'w')
outgeno = open('chr'+str(chrom)+'.'+str(start)+'.'+str(end)+'.YRI.mean.genotype.txt', 'w')

## Remove the first 2 lines of the whole snp file which are just headers :)
line = insnp.readline()
line = insnp.readline()

for line, gine in zip(insnp, ingeno):
    toks = line.strip().split()
    position = int(toks[5])
    if (position > end): break
    if (position > start):
        outsnp.write(line)
        outgeno.write(gine)
insnp.close()
outsnp.close()
ingeno.close()
outgeno.close()
