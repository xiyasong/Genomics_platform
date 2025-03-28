## 1. Download the software
```
conda install -c bioconda nextflow
chmod +x nextflow

./nextflow run hello
```
## 2. Running success with the sample test
```
nextflow run pgscatalog/pgsc_calc -profile test,conda
```

## 3. Merge
```
mv *passed.vcf.gz ../passed_data
bcftools merge *.vcf.gz -Oz -o 207-samples-Merged.vcf.gz --force-samples

###Another trying to see whether the scores are different
bcftools merge *.vcf.gz -Oz -o 207-samples-Merged-add-ref-parameter.vcf.gz --force-samples --missing-to-ref
```

###generated:
207-samples-Merged-add-ref-parameter.vcf.gz

## 4. Convert VCF to PLINK data
```
cd /Users/xiyas/software

#1) not adding --missing-to-ref
./plink2 --vcf /Users/xiyas/vcf_turkish/passed_data/207-samples-Merged.vcf.gz \
    --allow-extra-chr \
    --chr 1-22, X, Y, XY \
    -make-pgen --out /Users/xiyas/vcf_turkish/passed_data/207-samples-Merged

#2) adding --missing-to-ref
./plink2 --vcf /Users/xiyas/vcf_turkish/passed_data/207-samples-Merged-add-ref-parameter.vcf.gz \
    --allow-extra-chr \
    --chr 1-22, X, Y, XY \
    -make-pgen --out /Users/xiyas/vcf_turkish/passed_data/207-samples-Merged-add-ref-parameter

./plink2 --vcf /Users/xiyas/vcf-swedish/passed_data/101-samples-Merged.vcf.gz \
    --allow-extra-chr \
		--vcf-half-call m \
    --chr 1-22, X, Y, XY \
    -make-pgen --out /Users/xiyas/vcf-swedish/passed_data/101-samples-Merged
```

### terminal outputs 
```
Start time: Fri Mar 31 12:31:23 2023
16384 MiB RAM detected; reserving 8192 MiB for main workspace.
Using up to 10 threads (change this with --threads).
--vcf: 27257177 variants scanned (501644 skipped).
--vcf: /Users/xiyas/vcf_turkish/passed_data/207-samples-Merged-temporary.pgen +
/Users/xiyas/vcf_turkish/passed_data/207-samples-Merged-temporary.pvar.zst +
/Users/xiyas/vcf_turkish/passed_data/207-samples-Merged-temporary.psam written.
Warning: chrX is present in the input file, but no sex information was
provided; many commands will produce inaccurate results.  You are strongly
encouraged to rerun this import with --psam or --update-sex.  --split-par may
also be appropriate.
207 samples (0 females, 0 males, 207 ambiguous; 207 founders) loaded from
/Users/xiyas/vcf_turkish/passed_data/207-samples-Merged-temporary.psam.
27257177 variants loaded from
/Users/xiyas/vcf_turkish/passed_data/207-samples-Merged-temporary.pvar.zst.
Note: No phenotype data present.
Writing /Users/xiyas/vcf_turkish/passed_data/207-samples-Merged.psam ... done.
Writing /Users/xiyas/vcf_turkish/passed_data/207-samples-Merged.pvar ... done.
Writing /Users/xiyas/vcf_turkish/passed_data/207-samples-Merged.pgen ... done.
End time: Fri Mar 31 12:34:11 2023
(base)
```

## 5. Try the PRS score calculation
```
###T2D diseases risk 
####    error
#nextflow run pgscatalog/pgsc_calc \
#    -profile conda \
##    --input /Users/xiyas/pgsc_calc/samplesheet_test.csv --target_build GRCh38 \
#   --pgs_id PGS000014 \
#	--min_overlap 0.5

###the samplesheet_plink.csv (use plink2 format to calculate)
sampleset,vcf_path,bfile_path,pfile_path,chrom
### as the format as this:
Turkish_cohort,,,/Users/xiyas/vcf_turkish/plink2_converted/207-samples-Merged,

#####Try for a whole afternoon:
## The error happens because of you have created results and 'work' directory last time, so it may conflict;

###If you remove the results directory and work directory: then things worked.
##Every time make a new directory to run inside
nextflow run pgscatalog/pgsc_calc \
    -profile conda \
    --input /Users/xiyas/pgsc_calc/samplesheet_plink.csv \
    --target_build GRCh38 \
    --pgs_id PGS000014
```

<img width="910" alt="image" src="https://github.com/user-attachments/assets/7316e894-3398-4bd8-945f-0943bb878051" />

####Calculate the swedish cohort

###make a bash script for filtering 
###passed_extract.sh

```
### Important: this need to be runned in the /Users/xiyas/vcf-swedish/original_data

for FILE in /Users/xiyas/vcf-swedish/original_data/*filtered.vcf.gz
do
        echo $(basename $FILE)
        zless -S $(basename $FILE) | grep -E "^#|PASS" | bgzip > $(basename $FILE)_passed.vcf.gz
        mv $(basename $FILE)_passed.vcf.gz /Users/xiyas/vcf-swedish/passed_data
done
```
#Then index:
```
for F in *.vcf.gz ; do tabix -f -p vcf ${F}  ; done
```

##Then merge:
```
bcftools merge *.vcf.gz -Oz -o 101-samples-Merged.vcf.gz --force-samples
```


## 6: multiple scores in parallel 
```
cd /Users/xiyas/pgsc_calc/Swedish_4_score
nextflow run pgscatalog/pgsc_calc \
    -profile conda \
    --input /Users/xiyas/pgsc_calc/samplesheet_plink_swedish.csv \
    --target_build GRCh38 \
    --pgs_id PGS000013,PGS000014,PGS000016,PGS002758

cd /Users/xiyas/pgsc_calc/Turkish_4_score
nextflow run pgscatalog/pgsc_calc \
    -profile conda \
    --input /Users/xiyas/pgsc_calc/samplesheet_plink_turkish.csv \
    --target_build GRCh38 \
    --pgs_id PGS000013,PGS000014,PGS000016,PGS002758

```
