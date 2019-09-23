# Bayesian Inference of Neuronal Assemblies
This repository contains a C++ implementation of the Bayesian inference method for detecting neuronal assemblies developed in [G. Diana, T. Sainsbury and M. Meyer, bioRxiv 452557](https://doi.org/10.1101/452557)
## Dependencies
* [GNU Scientific Libraries](https://www.gnu.org/software/gsl/)
* [armadillo](http://arma.sourceforge.net/)
## Installation
From command-line run `make`
## How to run
`./gibbsDPA5data <NITER> <BURN_IN> <TRIM> <ASSEMBLIES> <SEED> <BINARY_FILE> <THRESH> <THRESH2> <folder> <continue>`

where

1. `NITER`: number of iterations
1. `BURN_IN`: number of initial MCMC steps excluded
1. `TRIM`: number of MCMC steps between recorded samples
1. `ASSEMBLIES`: initial number of assemblies
1. `SEED`: random seed
1. `BINARY_file`: input file in matrix format [neurons]x[times] where row *i* represents the binary activity of neuron *i*.
1. `THRESH2`: minimum neuronal activity (row sums)
1. `THRESH`: minimum number of synchronously active neurons
1. `folder`: output folder - being created if not already existing 
1. `continue`: uses data from previous run when set equal to 1 otherwise should be set to 0

## Example



