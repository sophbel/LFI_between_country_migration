# LFI_between_country_migration
Repository containing the data for "Estimating Between Country Migration in  Pneumococcal Populations"

## Pre-Processing
1) Panaroo
2) Select 'neutral' genes
3) Concatenate 'neutral' genes
4) Remove linked SNP sites
5) Create index file for each deme or country
   
## Run Asymmetric 2-Deme Model
1) Many paths are hard-coded so adjust accordingly
2) Set flags for python script
3) Create input and output folders for each GPSC/lineage
4) Asymmetric 2 Deme Script
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



