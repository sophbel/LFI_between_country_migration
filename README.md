# Estimating Between Country Migration in Pneumococcal Populations using Likelihood-free Inference 
Repository containing the data for "Estimating Between Country Migration in  Pneumococcal Populations". This repository utilizes the Engine For Likelihood Free Inference and coalescent simulation software msprime.
```
ELFI: https://elfi.readthedocs.io/en/latest/
msprime: https://tskit.dev/msprime/docs/stable/CITATION.html
```

## Setup environments
1) Set up environment to concatenate neutral genes of choice.
```
conda env create -f ./environments/core_aln.yml
```
2) Set up environment to run msprime and and all other dependencies within the migration parameter inference scripts.
```
conda env create -f ./environments/LF_migration.yml
```
# Overview of all steps. For specifics navigate to each folder.  
## Pre-Processing
1) Run panaroo using annotation files to create gene alignments. Panaroo paralellizes very well so multi-threading is a good option. We ran with the sensitive mode. <br>
   Documentation: https://gtonkinhill.github.io/panaroo/#/ <br>
   Manuscript: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02090-4 
3) Select 'neutral' genes
4) activate core_aln environment
```
conda activate core_aln
```
5) Concatenate 'neutral' genes
```
./PreProcessing/GPS_Panaroo/concatenate_core_genome.py
```
7) Remove linked SNP sites
8) Create index file for each deme or country using the outputs from the previous step as described in ./PreProcessing/README.md
```
./PreProcessing/indexalignments_trailHead_GPSCs.R 
```

## Run Asymmetric Models
1) Many paths are hard-coded so adjust accordingly
2) Set flags for python script
3) Create input and output folders for each GPSC/lineage
4) activate LF_migration environment
```
conda activate LF_migration
```
5) Run Models
   Model 1) Asymmetric 2 Deme Model <br>
   Script Usage <br> 
```
   usage: ELFI_2Deme_LDFilts.py [-h] [--gpsc GPSC] [--genes GENES] [--true_data TRUE_DATA] [--country1 COUNTRY1] [--country2 COUNTRY2] [--country3 COUNTRY3]
                             [--country4 COUNTRY4] [--input_dir INPUT_DIR] [--output_dir OUTPUT_DIR] [--evidence EVIDENCE] [--sample SAMPLE] [--bounds BOUNDS]
                             [--suffix SUFFIX] [--sampler SAMPLER] [--prior PRIOR] [--vcf_in VCF_IN]

Running this code to test different GPSC and number of neutral genes impact on migration parameter estimates

optional arguments:
  -h, --help            show this help message and exit
  --gpsc GPSC           Add a GPSC number from options 2,8,5,22,26,10 (default 8)
  --genes GENES         Add number of genes either 81 or 355 (default 81)
  --true_data TRUE_DATA
                        Add either [simulated] or [true] to test with simulated data or get the truth (default simulated)
  --country1 COUNTRY1   Add [sa] for South Africa, [mal] for Malawi, [gam] for Gambia, or [ken] for Kenya (default sa)
  --country2 COUNTRY2   Add [sa] for South Africa, [mal] for Malawi, [gam] for Gambia, or [ken] for Kenya (default mal)
  --country3 COUNTRY3   Add [sa] for South Africa, [mal] for Malawi, [gam] for Gambia, or [ken] for Kenya (do not include in 2Deme)
  --country4 COUNTRY4   Add [sa] for South Africa, [mal] for Malawi, [gam] for Gambia, or [ken] for Kenya (do not include in 2Deme)
  --input_dir INPUT_DIR
                        Path just before GPSCX/ in which the core alignment file is (default local)
  --output_dir OUTPUT_DIR
                        Path just before 2Deme/outputs/ inclusive to write summary stats and plots (default local)
  --evidence EVIDENCE   Number of evidence points for BOLFI (default 200)
  --sample SAMPLE       Number of samples for BOLFI (default 2000)
  --bounds BOUNDS       Upper bound for parameter (default 1)
  --suffix SUFFIX       Suffix for saved files (default "")
  --sampler SAMPLER     Type of sampler. Either [metropolis] or [nuts] (default metropolis)
  --prior PRIOR         Prior distribution. Either [exponential] or [uniform] (default exponential)
  --vcf_in VCF_IN       Full name of VCF file for input

```

Model 2) Symmetric 4 Deme Model <br>
Script Usage <br>
```
usage: ELFI_4Deme_Symmetric.py [-h] [--gpsc GPSC] [--genes GENES] [--true_data TRUE_DATA] [--country1 COUNTRY1] [--country2 COUNTRY2] [--country3 COUNTRY3]
                               [--country4 COUNTRY4] [--input_dir INPUT_DIR] [--output_dir OUTPUT_DIR] [--evidence EVIDENCE] [--sample SAMPLE]
                               [--bounds BOUNDS] [--suffix SUFFIX] [--sampler SAMPLER] [--prior PRIOR] [--epistasis EPISTASIS]

Running this code to test different GPSC and number of neutral genes impact on migration parameter estimates

optional arguments:
  -h, --help            show this help message and exit
  --gpsc GPSC           Add a GPSC number from options 2,8,5,22,26,10 (default 8)
  --genes GENES         Add number of genes either 81 or 355 (default 81)
  --true_data TRUE_DATA
                        Add either [simulated] or [true] to test with simulated data or get the truth (default simulated)
  --country1 COUNTRY1   Add [sa] for South Africa, [mal] for Malawi, [gam] for Gambia, or [ken] for Kenya (default sa)
  --country2 COUNTRY2   Add [sa] for South Africa, [mal] for Malawi, [gam] for Gambia, or [ken] for Kenya (default mal)
  --country3 COUNTRY3   Add [sa] for South Africa, [mal] for Malawi, [gam] for Gambia, or [ken] for Kenya (default ken)
  --country4 COUNTRY4   Add [sa] for South Africa, [mal] for Malawi, [gam] for Gambia, or [ken] for Kenya (default gam)
  --input_dir INPUT_DIR
                        Path just before GPSCX/ in which the core alignment file is (default local)
  --output_dir OUTPUT_DIR
                        Path just before 2Deme/outputs/ inclusive to write summary stats and plots (default local)
  --evidence EVIDENCE   Number of evidence points for BOLFI (default 4000)
  --sample SAMPLE       Number of samples for BOLFI (default 10000)
  --bounds BOUNDS       Upper bound for parameter (default 3)
  --suffix SUFFIX       Suffix string (default "")
  --sampler SAMPLER     Type of sampler. Either [metropolis] or [nuts] (default metropolis)
  --prior PRIOR         Prior distribution. Either [exponential] or [uniform] (default exponential)
  --epistasis EPISTASIS
                        Indicate whether epistatic sites within 1kb with >0.5r2 have been removed

Have fun!
```

## Post-Processing
Navigate to the post processing folder where you can use the output posterior migration parameter files to calculate the relative mean migration, directional migration, and summary tables for each deme set and GPSC.

## Data availability
Metadata is included in the Metadata folder. ENA accession numbers are included in the additional data file "ERRs.csv". 

