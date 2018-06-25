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

By default MiMMAl will produce plots representing the fits produced, including the results of the initial search of fits using expectation maximisation to search for a range for _sd_, as well as the global and local searches of parameter space including means and _sd_ in the current working directory. This can be disabled in the options for MiMMAl.

## Additional parameters

There are some additional parameters that can be set in runMiMMAl as required:

* `min.snps` refers to the minimum number of SNPs required for mixture modelling. Default: 10.
* `sd.width` is the fraction of which the range of _sd_ is set either side of the maxima of kernel density smoothing of the _sd_ of the initial fits of mixture modelling using expectation maximisation. Default: 1/3.
* `preset.sd` if this value is defined, the initial fit will not take place and a range either side of this value as defined by `sd.width` will be used for the global search. Default: NULL.
* `seed` the seed can be set to allow for reproducibility. Default: 1.
* `baf.res` or BAF resolution defines the number of intervals in the BAF values. This is defined as `10^-baf.res`. Therefore `baf.res=2` produces intervals of 0.01 between 0 and 0.5 in the global search for BAF mean. Increasing this increases the number of combinations searched and will increase computational time exponentially. The subsequent local search increases the resolution of the fit further, so MiMMAl will always fit each segment to a higher resolution than this initial global search. Default: 2.
* `use.ks.gate` refers to performing a Kolmogorov-Smirnov test prior to mixture modelling a segment. Default: TRUE.