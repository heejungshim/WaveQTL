
This repository contains WaveQTL, a software implementing a wavelet-based approach for genetic association analysis of functional phenotypes (e.g. sequence data arising from high-throughput sequencing assays), described in [Shim and Stephens 2015](https://github.com/heejungshim/WaveQTL/blob/master/doc/paper/Shim_2014.pdf).

We modified source codes in [BIMBAM](http://www.haplotype.org/bimbam.html) to implement a wavelet-based approach.

WaveQTL is a free software, you can redistribute it and/or modify it under
the terms of the GNU General Public License (GPLv3+).

The GNU General Public License does not permit this software to be
redistributed in proprietary programs.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

## WaveQTL

A key factor that distinguishes a wavelet-based approach from typical association analysis is to better exploit high-resolution measurements provided by high-throughput sequencing assays (see [Shim and Stephens 2015](https://github.com/heejungshim/WaveQTL/blob/master/doc/paper/Shim_2014.pdf) for motivations). The implemented wavelet-based methods aim to  

1. test for genetic association of functional phenotype
2. provide explanations for observed associations such as which parts and features of the data are driving the association
3. make applications to large-scale genetic association analyses computationally tractable.

Although wavelet methods were motivated primarily by genetic association studies for high-throughput sequencing data, they could also test for association between functional data and
other covariates, either continuous or discrete. For example, in a genomics context, it could be used to detect differences in gene expression (from RNA-seq data) or TF binding (from ChIP-seq data) measured on two groups (e.g. treatment conditions or cell types). Or it could be used to associate a functional phenotype, such as chromatin accessibility, with a continuous covariate, such as overall expression of a gene. It could also be used for genome-wide association studies of functional phenotypes unrelated to sequencing.
 
## Binary executable file

Binary executable file for Linux is in the `WaveQTL/bin/` directory (complied on 06/16/2014).

## Installation

cd into the `WaveQTL/` directory

    make all

Then, binary executable file will be created in the `WaveQTL/` directory.

## User manual 

User manual is in the `WaveQTL/doc/manual` directory.

## dsQTLs and analysis in Shim and Stephens (2015)

1. The main manuscript and supplementary materials of [Shim and Stephens (2015)](https://github.com/heejungshim/WaveQTL/blob/master/doc/paper/Shim_2014.pdf) are in the `WaveQTL/doc/paper` directory.
2. Information on dsQTLs identified by our analysis in [Shim and Stephens (2015)](https://github.com/heejungshim/WaveQTL/blob/master/doc/paper/Shim_2014.pdf) is in `WaveQTL/Shim_2014/` directory.
3. Information on the top 1% of 1024bp sites with the highest DNase I sensitivity, analyzed in [Shim and Stephens (2015)](https://github.com/heejungshim/WaveQTL/blob/master/doc/paper/Shim_2014.pdf), is in `WaveQTL/Shim_2014/` directory.
4. In the user manual, we describe how to perform an analysis in [Shim and Stephens (2015)](https://github.com/heejungshim/WaveQTL/blob/master/doc/paper/Shim_2014.pdf) for a given site, starting from reading DNase-seq data at the site on 70 individuals from raw data files (downloaded from [http://eqtl.uchicago.edu/dsQTL_data/](http://eqtl.uchicago.edu/dsQTL_data/)) to plotting estimated effect sizes in the data space as shown in the paper. 

## News

See the [`NEWS`](https://github.com/heejungshim/WaveQTL/blob/master/NEWS) file.

## Bug reports

Report bugs as issues on this repository or email the [mailing list](https://groups.google.com/forum/?hl=en#!forum/waveqtlusers).

## How to cite WaveQTL

Heejung Shim and Matthew Stephens (2015). [Wavelet-based genetic association analysis of functional phenotypes arising from high-throughput sequencing assays](https://github.com/heejungshim/WaveQTL/blob/master/doc/paper/Shim_2014.pdf). Ann. Appl. Stat. 9 (2015), no. 2, 665–686. 

## Author

[Heejung Shim](https://github.com/heejungshim) (University of Melbourne)
