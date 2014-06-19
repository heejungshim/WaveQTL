This directory contains information on dsQTLs identified by Shim and Stephens 2014. 


1. `WaveQTL/Shim_2014/entire_dataset_146435` contains information on the top 1% of 1024bp sites with the highest DNase I sensitivity (in total 146,435 sites: 50,000 randomly selected sites in `WaveQTL/Shim_2014/entire_dataset_146435/random_50000` and the remaining sites in `WaveQTL/Shim_2014/entire_dataset_146435/remaining`). Each file contains three columns with chromosome at the first column and start and end positions in the remaining columns. 

2. Information on 3,176 dsQTLs identified by our wavelet-based analysis at an FDR of 10% in Shim and Stephens 2014 are in `WaveQTL/Shim_2014/dsQTL_FDR10`. dsQTLs included in randomly selected 50,000 sites are in `dsQTL.info.from50000.txt` and the remaining dsQTLs are in `dsQTL.info.remaining.txt`. We provide locations (chr, start, end), p-value from the wavelet analysis (wavelets.pvalue), the most strongly associated SNP id (SNP.id), logLR from the SNP (logLR), and indicator whether this dsQTL is newly identified by the wavelet analysis (i.e., not overlapping with the 7,088 100bp windows reported as having dsQTLs in 2kb cis-candidate region from Degner et al. (2012)) (new.dsQTL), additionally p-value from a 100bp window based approach (windows.pvalue) and a 100bp window based approach with 50bp shift (windows.shift.pvalue) for dsQTLs in `dsQTL.info.from50000.txt`.


