## Preprocessing so that genes are filtered, concatenated, and indexed.

1) Run Panaroo
2) Extract neutral genes and concatenate using
   ```
   ./GPS_Panaroo/concatenate_core_genome.py
   ```
4) Install necessary tools: bcftools and snp-sites (on Sanger systems is as follows)
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
5) Remove all sites with >0.5 $r^2$ within 1kb
```
for f in $(cat gpscs);
do
bcftools +prune -l 0.5 -w 1000 ./GPSCs/GPSC${f}/out_gpsc${f}/core_81_gpsc${f}_alignment.snp.aln.biallelic.vcf -o ./RunModel/inputs/GPSC${f}/test.filter_r2_5_1kb.vcf
done
```
where gpscs is:
```
2
5
8
10
22
26
```
