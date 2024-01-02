#!/bin/bash -l
#SBATCH -A naiss2023-5-247
#SBATCH -p node
#SBATCH -n 9
#SBATCH -t 100:00:00
#SBATCH -J xiya_vep_turkish_cohort_43
#SBATCH --mail-user xiya.song@scilifelab.se
#SBATCH --mail-type=ALL


module load bioinfo-tools
module load vep
###Test 1: 4hours running
###Update: adding CLNREVSTAT, adding fork setting for faster speed, change core mode to node mode when submiting jobs,
### should adding phenotype annotation from HPO? No better separate
####submit this job: sh all_process_uppmax.sh P001_101.hard-filtered.vcf.gz

###SETTING THE INPUT FILES
echo  "The path for storing WGS vcf is /proj/sctatlas/nobackup/xiya-temp-WGS_SZA/vcf_turkish-new-batch"
#echo  "The path for storing WGS vcf is /crex/proj/snic2020-16-69/nobackup/WGS_SZA/ori_200_vcf_data_sza_2023_1/vcf_swedish/dir_rm_old_vep"
echo  "the sample file is $1"

########Process 1: filtered out the not-passed variants in the original vcf files
echo  "start remove the non-pass variants sites in the original vcf files"
#cd /crex/proj/snic2020-16-69/nobackup/WGS_SZA/ori_200_vcf_data_sza_2023_1/vcf_swedish/dir_rm_old_vep
cd /proj/sctatlas/nobackup/xiya-temp-WGS_SZA/vcf_turkish-new-batch
##no need for pass filter
zless $1 | grep -E "^#|PASS" | bgzip > ${1}_PASS.vcf.gz

echo  "succeed of get passed variants"
echo "the filtered file is ${1}_PASS.vcf.gz"
echo "move the filtered file to new filtered_dir"
echo "filtered_dir is /proj/sctatlas/nobackup/xiya-temp-WGS_SZA/vcf_turkish-new-batch/PASS_vcf_sza_turkish"
mv ${1}_PASS.vcf.gz /proj/sctatlas/nobackup/xiya-temp-WGS_SZA/vcf_turkish-new-batch/PASS_vcf_sza_turkish

########Process 2: vep annotation
echo "Start running vep annotation"
echo "vep files stored in /proj/sctatlas/nobackup/xiya-temp-WGS_SZA/New_V2_dir_vep_annotated_vcf_turkish"
#cd /home/szalab/ensembl-vep/
#gnomADe="r2.1.1" gnomADg="v3.1.2"
##Include allele frequency from NHLBI-ESP populations. Must be used with --cache.
##Include allele frequency from ExAC project populations. Must be used with --cache.
## Include allele frequency from Genome Aggregation Database (gnomAD) exome populations.
# Include allele frequency from Genome Aggregation Database (gnomAD) genome populations.
## Report the highest allele frequency observed in any population from 1000 genomes, ESP or gnomAD.
#Add the global allele frequency (AF) from 1000 Genomes Phase 3 data for any known co-located variant to the output
# This is the customized Clinvar file including P, LP and any conflicting interpretations once submitted with P/LP.
#--plugin REVEL,/crex/proj/snic2020-16-69/nobackup/WGS_SZA/database-file/
#Splicing variants db

vep --cache --dir_cache $VEP_CACHE \
--fork 9 \
--dir_plugins /sw/data/vep/107/Plugins \
-i /proj/sctatlas/nobackup/xiya-temp-WGS_SZA/vcf_turkish-new-batch/PASS_vcf_sza_turkish/${1}_PASS.vcf.gz \
-o /proj/sctatlas/nobackup/xiya-temp-WGS_SZA/New_V2_dir_vep_annotated_vcf_turkish/${1}_vep_annotated.vcf \
--force_overwrite \
--assembly GRCh38 \
--symbol --vcf --check_existing --variant_class \
--sift b --polyphen b \
--hgvs \
--fasta /crex/proj/snic2020-16-69/nobackup/WGS_SZA/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
--canonical \
--af --af_gnomade --af_gnomadg --max_af \
--custom /proj/snic2020-16-69/nobackup/WGS_SZA/database-file-v2/clinvar_20221129.P.LP.modify.vcf.gz,ClinVar,vcf,exact,0,ID,CLNSIG,CLNDN,CLNHGVS,CLNSIGINCL,CLNVC,GENEINFO,CLNDISDB,CLNSIGCONF,CLNREVSTAT \
--custom /proj/snic2020-16-69/nobackup/WGS_SZA/database-file-v2/Whole_out_sorted.vcf.gz,Database,vcf,exact,0,Type,SZAID \
--plugin dbscSNV,/crex/proj/snic2020-16-69/nobackup/WGS_SZA/database-file-v2/dbscSNV1.1/dbscSNV1.1_GRCh38.txt.gz \
--plugin REVEL,/crex/proj/snic2020-16-69/nobackup/WGS_SZA/database-file-v2/REVEL/new_tabbed_revel_grch38.tsv.gz \
--plugin dbNSFP,/crex/proj/snic2020-16-69/nobackup/WGS_SZA/database-file-v2/dbNSFP4.3a_new/dbNSFP4.3a_grch38.gz,REVEL_score,REVEL_rankscore,BayesDel_addAF_score,BayesDel_addAF_rankscore,BayesDel_addAF_pred,BayesDel_noAF_score,BayesDel_noAF_rankscore,BayesDel_noAF_pred \
--plugin SpliceAI,snv=/crex/proj/snic2020-16-69/nobackup/WGS_SZA/database-file-v2/spliceAI/spliceai_scores.masked.snv.hg38.vcf.gz,indel=/crex/proj/snic2020-16-69/nobackup/WGS_SZA/database-file/spliceAI/spliceai_scores.masked.indel.hg38.vcf.gz


echo "vep with plugins finished"

########Process 3: python script processing
#echo "The results file will saved in the dir /proj/snic2020-16-69/nobackup/WGS_SZA/python_outputfile"
#python3 GenomeReportingPipeline.py /proj/snic2020-16-69/nobackup/WGS_SZA/vep_annotated_vcf/${1}_vep_annotated.vcf /database-file/GeneDB.txt /database-file/DiseaseDB_with_Description_OMIM_2000.txt /database-file/pheno_OMIM_all.txt /proj/snic2020-16-69/nobackup/WGS_SZA/python_outputfile/${1}_python_outfile.txt

#echo "python processing succeed"

wait
