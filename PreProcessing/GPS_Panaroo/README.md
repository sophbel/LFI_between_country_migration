## This script concatenates genes from a folder into a single alignment and is extracted from panaroo
core_aln.yml file should be activated into a conda environment for this file.
```
usage: concatenate_core_genome.py [-h] -o OUTPUT_DIR [--core_threshold CORE]
                                  [--filtered-core PREFILTERED]
                                  [--prefix PREFIX]

Concatenate recombination-free core gene alignmets

optional arguments:
  -h, --help            show this help message and exit
  --core_threshold CORE
                        Core-genome sample threshold (default=0.95)
  --filtered-core PREFILTERED
                        Provide a folder in the panaroo output containing a
                        filtered set of core genes
  --prefix PREFIX       Provide outprefix

Input/output:
  -o OUTPUT_DIR, --out_dir OUTPUT_DIR
                        location of the Panaroo output directory
```
