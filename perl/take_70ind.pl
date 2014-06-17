## `take_70ind.pl' contains perl script to extract a subset of individuals from mean genotype files. You need to change paths to input and output files.
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


# !/usr/local/bin/perl;
# Input file 
# 1. "index" : indicator to indicate whether ith column in "geno" should be extracted (1) or not 
# 2. "geno" :  full data set from which a subset of columns will be extracted.
# Note : the numbers of lines in "index" file and the number of columns in "geno" file should be the same
#
# Output file 
# 1. "FILE2" : extracted data
# 2. "FILEerror" : error message


$file = '~/WaveQTL/data/Shim_2014_etc/indicator_include.txt';
open(FILE1, $file);
@index = <FILE1>; # 1 : to be used  0 : to be skipped
close(FILE1);

# you need to change this path
$file = '/mnt/lustre/home/jdegner/IMPUTED_YRI_7_JUL_11/FINAL.BIMBAM.INPUT/OUTDIR/output/chr21.YRI.mean.genotype.txt';
open(FILE1, $file);
@geno = <FILE1>;
close(FILE1);

# you need to change these paths
open(FILE2, "> /mnt/lustre/home/shim/wavelets/data/genotype/chr21.YRI.70.mean.genotype.txt");
open(FILEerror, "> /mnt/lustre/home/shim/wavelets/data/genotype/chr21.error.out");


$genoIX = $#geno;
$indexIX = $#index;


print "genoNum: ";
print $genoIX;
print " indexNum: ";
print $indexIX;


foreach $ind  (@geno[0..$genoIX]){

    chop($ind);   
    #@wd = split(/ +/,$ind );
    @wd = split(/\t/,$ind );
    $wdIX = $#wd;


    if($wdIX ne  $indexIX){
	print FILEerror "ERROR! Not matched number! \n";
	print FILEerror $wdIX;
	print FILEerror " ";
	print FILEerror $indexIX;
	print FILEerror "\n";
    }else{
	$i = 0;
	foreach $subwd (@wd[0..$wdIX]){
	    if(@index[$i] == 1){
		print FILE2 $subwd;
		print FILE2 " ";
            }
	    $i++;
        }
	print FILE2 "\n";
    }

}
close(FILE2);
close(FILEerror);


