#! /usr/bin/env python

import sys

if (len(sys.argv) != 4):
    print "Usage: extract.py CHROMOSOME START END"
    sys.exit(1)

chrom = sys.argv[1]
start = int(sys.argv[2]) #- 2000
end = int(sys.argv[3]) #+ 2000
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
