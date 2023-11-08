## Preprocessing so that genes are filtered, concatenated, and indexed.

1) Run Panaroo
2) Extract neutral genes and concatenate
3) Install necessary tools: bcftools and snp-sites (on Sanger systems is as follows)
```
module load bcftools/1.9--h68d8f2e_9
module load snp-sites/2.5.1--hed695b0_0 
```
3) Filter bi-allelic SNP sites
```
###convert to biallelic vcf
snp-sites -v -o PREFIX_core.snp.aln.vcf PREFIX_all_alignment.aln
bcftools view -m2 -M2 -v snps PREFIX_core.snp.aln.vcf > PREFIX_core.snp.aln.biallelic.vcf
```
4) Pull Lane IDs from concatenated genes for indexing by country
```
grep ">" PREFIX_all_alignment.aln > tmp.lanes
### cut the ids so they only include the lane id removing the index, the >, and .velvet
cut -f1 -d';' tmp.lanes | cut -f2 -d'>' | cut -f1 -d'.' > core_PREFIX.laneIDS
```
