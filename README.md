# BINE: Bayesian Inference of Neuronal Ensembles

### Description

This repository contains a C++ implementation of the Bayesian inference method for detecting neuronal ensembles developed in [G. Diana, T. Sainsbury and M. Meyer, bioRxiv 452557](https://doi.org/10.1101/452557) (Algorithm 2). The main code `gibbsDPA5data` requires as input a text file containing the binary matrix of neuronal activity where each row contains the binary trace of a given recorded neuron.  

### Dependencies
* [GNU Scientific Libraries](https://www.gnu.org/software/gsl/)
* [armadillo](http://arma.sourceforge.net/) (version >7.200)

### Installation
Run `make` from command line. The makefile will create the folder `bin` for the binary files and `.obj` for objects.

## Data analysis
The main program for data analysis is `gibbsDPA5data`

### Synopsis
``` 
gibbsDPA5data --file=<FILE> --folder=<FOLDER> --niter=<VALUE> --assemblies=<VALUE> [OPTIONS] ...
```

### Required input
**--file [FILE]**

> input file in matrix format [neurons]x[times] where row *i* represents the binary activity of neuron *i*.

**--folder [FOLDER]**

> output folder - being created if not already existing 

**--niter [ITERATIONS]**

> number of iterations of the Markov Chain

**--assemblies [VALUE]**

> initial number of assemblies

### Optional input

**--trim [VALUE]**

> number of MCMC steps between recorded samples. Default 1. 

**--burn_in [VALUE]**

> number of initial MCMC steps excluded

**--seed [VALUE]**

> random seed

**--min_neur [VALUE]**

> minimum number of synchronously active neurons. Default 0.

**--min_act [VALUE]**

> minimum neuronal activity (row sums). Default 0.

**--recorded_assemblies [VALUE]**

> number of assemblies for which properties (activity, synchrony and asynchrony) are written in corresponding output files. Default 100.

**--continue**

> uses data from previous run stored in `folder`

**--verbose**

> show details of the Markov Chain.

## Generate testing data
Testing data generated from the assembly model can be simulated by the command

```
./bin/generate_data <NCELLS> <TIMES> <ASSEMBLIES> <SEED> <LAMBDA0> <LAMBDA1> <ACTIVITY> <outfile> <lastIsFree>
```

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

```
./bin/generate_data 400 1000 4 1 0.04 0.7 0.3 testing
```

where assembly activity, synchrony and asynchrony were set to 30%, 70% and 4% respectively.

To analyze this dataset run the command

```
./bin/gibbsDPA5data --niter=1000 --assemblies=100 --file=testing/binary_matrix.dat --folder=testing
```

## Results
Posterior samples of latent variables and model parameters are stored in dedicated files created in the folder specified in the input. By default, Ensemble properties such as activity, synchrony and asynchrony are written in output files for the first 100 ensembles. This number can be changed by specifying the option `--recorded_assemblies` in the input. Here is a list of files generated:

* **membership_traj.dat**: posterior samples of all neuronal membership by row. 
* **pmu.dat**: posterior samples of activity levels for all ensembles by row.
* **lambda0.dat**: posterior samples of asynchrony levels for all ensembles by row.
* **lambda1.dat**: posterior samples of synchrony levels for all ensembles by row.
* **F.dat**: likelihood for all MCMC samples. This can be used to monitor the stability of the Markov Chain and establish convergence criteria.
* **n.dat**: group sizes of all ensembles by row.
* **omega_traj.dat**: ensemble activity matrix. By default, for each posterior sample, the binary activity of the first 100 ensembles is written as a matrix 100xM where M is the number of synchronous time frames considered. The number of rows can be changed by specifying the option `--recorded_assemblies` in the input.
* **P.dat**: posterior samples of the number of assemblies. 
* **selection.dat**: list of neurons included in the analysis according to the thresholds on activity and number of active neuron per synchronous event.
* **txsweep.dat**: fraction of neurons changing membership across the MCMC.
