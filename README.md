# Genomics_platform

A repository for the preparing manuscript "A platform for genomic data analysis, reporting, and interpretation"

## Overview

------------------------------------------------------------------------

### For Turkish cohort (n=275)

## Per-sample(autosomes) statistics 
Results file in ```/Users/xiyas/pgsc_calc/vcf_turkish-with-merge/autosome_filtered```
**Per-sample** variants statistics for autosomes

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

## Cohort-level variants findings and ploting

**Cohort-level** statistics
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
- Fig.4A Get a smaller version of annotated VCF, for later figure plotting
```function_biallelic-275-samples-Merged-add-ref-parameter.vcf.gz_vep_annotated.vcf.rm.missing.gz```

```
bcftools +split-vep biallelic-275-samples-Merged-add-ref-parameter.vcf.gz_vep_annotated.vcf.rm.missing.gz -f '%CHROM %POS %REF %ALT %Existing_variation %Consequence %VARIANT_CLASS %IMPACT %MAX_AF %MAX_AF_POPS \n' -s worst > function_biallelic-275-samples-Merged-add-ref-parameter.vcf.gz_vep_annotated.vcf.rm.missing.gz
## check by wc -l == 35725453 lines
```

- get Allele frequencies count for Turkish cohort
```
vcftools --gzvcf biallelic-275-samples-Merged-add-ref-parameter.vcf.gz --freq --out freq_biallelic
```
- check the amounts and remove missing alt alleles
```
less -S freq_biallelic.frq | wc -l
#### A total of 41578511 lines before removing missing alt alleles
LC_ALL=C grep -Ev '/*' freq_biallelic.frq >freq_biallelic_rm_missing.frq
#### A total of 35725454 -1(header line) = 35,725,453 distinct alternate alleles
```


### conda (recommended)

To create a conda environment named `btyper3` and install BTyper3 and all of its dependencies:

1. Install [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html), if necessary
2. Create a new environment named `btyper3` by running the following command from your terminal:
   ```console
   conda create -n btyper3
   ```

### pip

1. To run BTyper3, please download and install the following dependencies, if necessary:

   - [Python 3](https://www.python.org/downloads/)
   - [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)


2. [Add BLAST+ to your path](https://unix.stackexchange.com/questions/26047/how-to-correctly-add-a-path-to-path), if necessary (to check if BLAST+ is in your path, try running `makeblastdb -h` and `tblastn -h` from your command line; you should get a help message for each command, with no error messages)

3. Install via `pip` (this will download required Python dependencies as well):
   ```console
   pip install btyper3  
   ```

------------------------------------------------------------------------

* <a href="https://www.tandfonline.com/doi/full/10.1080/10408398.2021.1916735">Review of *B. cereus* group taxonomy/nomenclature</a>

* <a href="https://journals.asm.org/doi/full/10.1128/mBio.00034-20">Standardized nomenclature for the *B. cereus* group</a>

* <a href="https://www.frontiersin.org/articles/10.3389/fmicb.2020.580691/full">Comparison of our standardized nomenclature to other *B. cereus* group typing methods (e.g., MLST, *panC*, ANI-based comparisons to species type strain genomes)</a>

------------------------------------------------------------------------

## Citation

### If you found the BTyper3 tool, its source code, and/or any of its associated databases useful, please cite:

Carroll, Laura M., Martin Wiedmann, Jasna Kovac. 2020. <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7042689/">"Proposal of a Taxonomic Nomenclature for the *Bacillus cereus* Group Which Reconciles Genomic Definitions of Bacterial Species with Clinical and Industrial Phenotypes."</a> *mBio* 11(1): e00034-20; DOI: 10.1128/mBio.00034-20.

Carroll, Laura M., Rachel A. Cheng, Jasna Kovac. 2020. <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7536271/">"No Assembly Required: Using BTyper3 to Assess the Congruency of a Proposed Taxonomic Framework for the *Bacillus cereus* group with Historical Typing Methods."</a> *Frontiers in Microbiology* 11: 580691; DOI: 10.3389/fmicb.2020.580691.

------------------------------------------------------------------------


## Quick Start

For detailed information, check out the <a href="https://github.com/lmc297/BTyper3/wiki">BTyper3 wiki</a>

### Command Structure

```
btyper3 -i [fasta] -o [output directory] [options...]
```

For help, type `btyper3 -h` or `btyper3 --help`

For your current version, type `btyper3 --version`

### Sample Commands

#### Perform all default analyses, using an assembled genome (complete or draft) in (multi-)FASTA format as input (assumes fastANI is in the user's path):

```
btyper3 -i /path/to/genome.fasta -o /path/to/desired/output_directory
```



------------------------------------------------------------------------


Disclaimer: BTyper3 is pretty neat! However, no tool is perfect, and BTyper3 cannot definitively prove whether an isolate is pathogenic or not. As always, interpret your results with caution. We are not responsible for taxonomic misclassifications, misclassifications of an isolate's pathogenic potential or industrial utility, and/or misinterpretations (biological, statistical, or otherwise) of BTyper3 results.
