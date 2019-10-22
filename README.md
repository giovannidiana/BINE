# BINE: Bayesian Inference of Neuronal Ensembles

### Description

This repository contains a C++ implementation of the Bayesian inference method for detecting neuronal ensembles developed in [G. Diana, T. Sainsbury and M. Meyer, bioRxiv 452557](https://doi.org/10.1101/452557) (Algorithm 2). The main code `gibbsDPA5data` requires as input a text file containing the binary matrix of neuronal activity where each row contains the binary trace of a given recorded neuron.  

### Dependencies
* [GNU Scientific Libraries](https://www.gnu.org/software/gsl/)
* [armadillo](http://arma.sourceforge.net/) (version >7.200)

### Installation
- Create local folders `.obj` and `bin`
- run `make`
All binary files will be installed in `bin`.

## Data analysis
The main program for data analysis is `gibbsDPA5data`

### Synopsis
``` 
gibbsDPA5data [OPTIONS]
```

### Options

**-i, --niter**
: number of iterations

**-t, --trim**
: number of MCMC steps between recorded samples

**-u, --burn_in**
: number of initial MCMC steps excluded

**-a, --assemblies**
: initial number of assemblies

**-s, --seed**
: random seed

**-b, --file**
: input file in matrix format [neurons]x[times] where row *i* represents the binary activity of neuron *i*.

**-1, --min_neur**
: minimum number of synchronously active neurons

**-2, --min_act**
: minimum neuronal activity (row sums)

**-f, --folder**
: output folder - being created if not already existing 

**-c, --continue**
: uses data from previous run stored in `folder`

## Generate testing data
Testing data generated from the assembly model can be simulated by the command

`./bin/generate_data <NCELLS> <TIMES> <ASSEMBLIES> <SEED> <LAMBDA0> <LAMBDA1> <ACTIVITY> <outfile> <lastIsFree>`

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


## Test run
Generate testing data of 400 neurons and 1000 time frames organized into 4 assemblies.

 `./bin/generate_data 400 1000 4 1 0.04 0.7 0.3 testing`

where assembly activity, synchrony and asynchrony were set to 30%, 70% and 4% respectively.

To analyze this dataset run the command

`./bin/gibbsDPA5data 1000 100 1 50 1 testing/binary_matrix.dat 0 0 testing 0`


