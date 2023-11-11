## Cohort-level variants findings and ploting
- preparing building index files for all passed VCF files
```
for F in *.vcf.gz ; do tabix -f -p vcf ${F} ; done
```

## For Turkish cohort (n=275)

### **Per-sample** variants statistics for autosomes
Results file in ```/Users/xiyas/pgsc_calc/vcf_turkish-with-merge/autosome_filtered```

- Get variants on autosomes region only
```
for FILE in P001*.vcf.gz;do zless -S $FILE | grep -w '^#\|#CHROM\|chr[1-9]\|chr[1-2][0-9]' | bgzip  > autosome_${FILE};done
for FILE in autosome_P001*.vcf.gz;do bcftools stats $FILE > ${FILE}.txt;done
```
- Get the Ts/Tv ratio and count for autosomes
```
for FILE in *.txt;do less -S $FILE | grep '^TSTV' >> sum_tstv.txt;done
for FILE in *.txt;do less -S $FILE | grep 'number of records:' >> sum_count.txt;done
```

###  **Cohort-level** statistics and get the frequency table

locate in ```/Users/xiyas/pgsc_calc/vcf_turkish_passed_data_275```
- Merge all samples to a whole VCF file
```
bcftools merge *.vcf.gz -Oz -o 275-samples-Merged-add-ref-parameter.vcf.gz --force-samples --missing-to-ref
```
- Split the multi-allelic sites and normalize the VCF file
```
bcftools norm -a -m -any 275-samples-Merged-add-ref-parameter.vcf.gz | bgzip > biallelic-275-samples-Merged-add-ref-parameter.vcf.gz
```
- Annotate the whole VCF file by VEP (UPPMAX) to get the file:```biallelic-275-samples-Merged-add-ref-parameter.vcf.gz_vep_annotated.vcf.gz```
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
```biallelic-275-samples-Merged-add-ref-parameter.vcf.gz_vep_annotated.vcf.rm.missing.gz```

```
LC_ALL=C zgrep -E '^#|SNV|insertion|deletion' biallelic-275-samples-Merged-add-ref-parameter.vcf.gz_vep_annotated.vcf.gz > biallelic-275-samples-Merged-add-ref-parameter.vcf.gz_vep_annotated.vcf.rm.missing.gz
```
- Get a smaller version of annotated VCF, for later figure plotting
```function_biallelic-275-samples-Merged-add-ref-parameter.vcf.gz_vep_annotated.vcf.rm.missing.gz```

```
bcftools +split-vep biallelic-275-samples-Merged-add-ref-parameter.vcf.gz_vep_annotated.vcf.rm.missing.gz -f '%CHROM %POS %REF %ALT %Existing_variation %Consequence %VARIANT_CLASS %IMPACT %MAX_AF %MAX_AF_POPS \n' -s worst > function_biallelic-275-samples-Merged-add-ref-parameter.vcf.gz_vep_annotated.vcf.rm.missing.gz
## check by wc -l == 35725453 lines
```

- Get Allele frequencies count for Turkish cohort
```
vcftools --gzvcf biallelic-275-samples-Merged-add-ref-parameter.vcf.gz --freq --out freq_biallelic
```
- Check the amounts and remove missing alt alleles
```
less -S freq_biallelic.frq | wc -l
#### A total of 41578511 lines before removing missing alt alleles
LC_ALL=C grep -Ev '/*' freq_biallelic.frq >freq_biallelic_rm_missing.frq
#### A total of 35725454 -1(header line) = 35,725,453 distinct alternate alleles
```

## For Swedish cohort (n=101)

### **Per-sample** variants statistics for autosomes 
this is not running because we are not the sequencer

### **Cohort-level** statistics and get the frequency table
locate in ```/Users/xiyas/pgsc_calc/vcf_swedish_passed_data_101```
- Merge all samples to a whole VCF file
```
bcftools merge *.vcf.gz -Oz -o 101-samples-Merged-add-ref-parameter.vcf.gz --force-samples --missing-to-ref
```
- Split the multi-allelic sites and normalize the VCF file
```
bcftools norm -a -m -any 101-samples-Merged-add-ref-parameter.vcf.gz | bgzip > biallelic-101-samples-Merged-add-ref-parameter.vcf.gz
```

- (only for Swedish) remove the old VEP annotation:
```
#rm_old_vep.sh
zgrep -E "^#" $1  > ${1}_header.vcf

echo "get passed variants and also remove previous vep annotation"

zgrep -E "^[^#].*PASS" $1  | perl -pe 's/CSQ=.*?;//' > ${1}.rmvep.PASS.vcf

echo "combine two files"
cat ${1}_header.vcf ${1}.rmvep.PASS.vcf > ${1}.rmvep.merged.PASS.vcf


echo "remove unneeded files"
rm ${1}_header.vcf ${1}.rmvep.PASS.vcf
```

- Annotate the whole VCF file by VEP (UPPMAX) to get the file:```biallelic-101-samples-Merged-add-ref-parameter.vcf.gz.rmvep.merged.PASS.annotated.vcf.gz```

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
LC_ALL=C zgrep -E '^#|SNV|insertion|deletion' biallelic-101-samples-Merged-add-ref-parameter.vcf.gz.rmvep.merged.PASS.annotated.vcf.gz > biallelic-101-samples-Merged-add-ref-parameter.vcf.gz_vep_annotated.vcf.rm.missing.gz
```
- Get a smaller version of annotated VCF, for later figure plotting
```
bcftools +split-vep biallelic-101-samples-Merged-add-ref-parameter.vcf.gz_vep_annotated.vcf.rm.missing.gz -f '%CHROM %POS %REF %ALT %Existing_variation %Consequence %VARIANT_CLASS %IMPACT %MAX_AF %MAX_AF_POPS \n' -s worst > function_biallelic-101-samples-Merged-add-ref-parameter.vcf.gz_vep_annotated.vcf.rm.missing.gz
```

- Get Allele frequencies count for Swedish cohort
```
vcftools --gzvcf biallelic-101-samples-Merged-add-ref-parameter.vcf.gz --freq --out freq_biallelic_101
```
- Check the amounts and remove missing alt alleles
```
less -S freq_biallelic_101.frq | wc -l
#### A total of ________ lines before removing missing alt alleles
grep -Ev '/*' freq_biallelic_101.frq >freq_biallelic_rm_missing_101.frq
#### A total of ________ -1(header line) = ________ distinct alternate alleles
```


