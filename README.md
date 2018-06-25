# MiMMAl
Using two-component Gaussian mixture modelling to estimate the major allele distribution in genotying data.

## Installation

Using R >= 3.3.2, MiMMAl depends on the following packages:

* ComplexHeatmap (>= 1.12.0)
* circlize (>= 0.3.10)
* cowplot (>= 0.7.0)
* ggplot2 (>= 2.2.1)
* grid (>= 3.3.2)
* mixtools (version 1.0.4)

Once those packages have been installed. MiMMAl can be downloaded and installed in the following way.

`install.packages("MiMMAl", repos = NULL, type="source")`

## Running

To run MiMMAl a tab-separated text file needs to be produces containing four columns containing; chromosome (_chr_), position (_pos_), raw BAF value (_BAF_) and the mean/median value of the segment of which the loci belongs (_BAFseg_), for heterozygous SNPs only.

The minimum requirements for running MiMMAl (`runMiMMAl`) are to include the `samplename`, this will be appended to `.BAFphased.txt`, the output of MiMMAl, that will be produced in the current working directory. You will also have to provide the path and name of the input text file as `inputfile`.

By default MiMMAl will produce plots representing the fits produced, including the results of the initial search of fits using expectation maximisation to search for a range for _sd_, as well as the global and local searches of parameter space including means and _sd_.