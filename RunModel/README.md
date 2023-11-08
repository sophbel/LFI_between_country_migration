### Folder containing necessary files to run the 2 Deme and 4 Deme Models
Included here are input files such as the index files for each country, vcf files for each set of 'neutral' genes, and a folder for the output files for each GPSC. 

## Inputs
The inputs folder exists and contains the vcf for each number of genes and GPSC.

## Outputs
To create the necessary output folders follow the steps below: 

1) create output folder<br />
```mkdir outputs```
2) create folders for each model type<br />
```
cd outputs
mkdir 2Deme
mkdir 4Demes
```
4) create folders for each GPSC<br />
```
cd 2Deme
mkdir GPSC2
```*
4) if testing recapture including simulated data with input migration migration<br />
```cd GPSC2
mkdir truth
mkdir sims
```*

*for each GPSC

