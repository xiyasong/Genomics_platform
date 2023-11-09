## Cohort-level variants findings and ploting
- preparing building index files for all passed VCF files
```
for F in *.vcf.gz ; do tabix -f -p vcf ${F} ; done
```

### For Swedish cohort (n=101)

1.**Per-sample** variants statistics for autosomes is not running because we are not the sequencer

2. **Cohort-level** statistics to get the frequency table
locate in ```/Users/xiyas/pgsc_calc/vcf_swedish_passed_data_101```
- Merge all samples to a whole VCF file
```
bcftools merge *.vcf.gz -Oz -o 101-samples-Merged-add-ref-parameter.vcf.gz --force-samples --missing-to-ref
```
- Split the multi-allelic sites and normalize the VCF file
```
bcftools norm -a -m -any 101-samples-Merged-add-ref-parameter.vcf.gz | bgzip > biallelic-101-samples-Merged-add-ref-parameter.vcf.gz
```
- Annotate the whole VCF file by VEP (UPPMAX) to get the file:```biallelic-101-samples-Merged-add-ref-parameter.vcf.gz_vep_annotated.vcf.gz```
```
vep --cache --dir_cache $VEP_CACHE \
--fork 9 \
--dir_plugins /sw/data/vep/107/Plugins \
-i $1 \
-o ${1}_vep_annotated.vcf.gz \
--compress_output bgzip \
--force_overwrite \
--assembly GRCh38 \
--symbol --vcf --check_existing --variant_class \
--canonical \
--af --af_gnomade --af_gnomadg --max_af \
--custom /proj/snic2020-16-69/nobackup/WGS_SZA/database-file-v2/Whole_out_sorted.vcf.gz,Database,vcf,exact,0,Type,SZAID
```

- After splitting, remove the sites with missing alt allels and get the file:
```biallelic-101-samples-Merged-add-ref-parameter.vcf.gz_vep_annotated.vcf.rm.missing.gz```

```
LC_ALL=C zgrep -E '^#|SNV|insertion|deletion' biallelic-101-samples-Merged-add-ref-parameter.vcf.gz_vep_annotated.vcf.gz > biallelic-101-samples-Merged-add-ref-parameter.vcf.gz_vep_annotated.vcf.rm.missing.gz
```

- Fig.4A get needed column from whole VEP annotated Swedish file
```
bcftools +split-vep biallelic-101-samples-Merged-add-ref-parameter.vcf.gz_vep_annotated.vcf.rm.missing.gz \n
-f '%CHROM %POS %REF %ALT %Existing_varicdation %VARIANT_CLASS \n' -d > small_biallelic-275-samples-Merged-add-ref-parameter.vcf.gz_vep_annotated.vcf.rm.missing.gz

### need to remove the duplicated records.
zless -S small_biallelic-101-samples-Merged-add-ref-parameter.vcf.gz_vep_annotated.vcf.rm.missing.txt.gz | awk '!a[$1$2$3$4]++' > rm_dup_small_biallelic-101-samples-Merged-add-ref-parameter.vcf.gz_vep_annotated.vcf.rm.missing.txt
```
- get Allele frequencies count for Turkish cohort
```
vcftools --gzvcf biallelic-101-samples-Merged-add-ref-parameter.vcf.gz --freq --out freq_biallelic_101
```
- check the amounts and remove missing alt alleles
```
less -S freq_biallelic_101.frq | wc -l
#### A total of ________ lines before removing missing alt alleles
LC_ALL=C grep -Ev '/*' freq_biallelic_101.frq >freq_biallelic_rm_missing_101.frq
#### A total of ________ -1(header line) = ________ distinct alternate alleles
```

2. Fig.4B Consequences drawing, with most severe
```
### now is running this
bcftools +split-vep biallelic-275-samples-Merged-add-ref-parameter.vcf.gz_vep_annotated.vcf.rm.missing.gz -f '%CHROM %POS %REF %ALT %Existing_variation %Consequence %VARIANT_CLASS %IMPACT %MAX_AF %MAX_AF_POPS \n' -s worst > function_biallelic-275-samples-Merged-add-ref-parameter.vcf.gz_vep_annotated.vcf.rm.missing.gz
```
