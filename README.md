# Bayesian Inference of Neuronal Assemblies
This repository contains a C++ implementation of the Bayesian inference method for detecting neuronal assemblies developed in [G. Diana, T. Sainsbury and M. Meyer, bioRxiv 452557](https://doi.org/10.1101/452557)
## Dependencies
* [GNU Scientific Libraries](https://www.gnu.org/software/gsl/)
* [armadillo](http://arma.sourceforge.net/)
## Installation
From command-line run `make`
## Run
### Data analysis
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

### Generate testing data
Testing data generated from the assembly model can be simulated by the command

`./generate_data <NCELLS> <TIMES> <ASSEMBLIES> <SEED> <LAMBDA0> <LAMBDA1> <ACTIVITY> <outfile> <lastIsFree>`

where 

1. `NCELLS`: number of cells in the dataset
1. `TIMES`: number of time frames
1. `ASSEMBLIES`: number of assemblies
1. `SEED`: random seed
1. `LAMBDA0`: average asynchrony level
1. `LAMBDA1`: average synchrony level
1. `ACTIVITY`: average activity level
1. `outfile`: result folder containing the binary matrix (`outfile/binary_matrix.dat`) as well as the the original membership for all neurons (`outfile/membership_orig.dat`).
1. `lastIsFree`: optional parameter. `lastIsFree=1` sets equal synchrony and asynchrony for the last assembly, meaning that the last assembly is made by 'free' neurons.


## Example

`./gibbsDPA5data 1000 100 1 

