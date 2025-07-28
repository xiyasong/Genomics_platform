### all traits first annotaed by vep to get the chromosome location and reference allele

###sh file ====================
#module load bioinfo-tools
#module load vep
#vep --cache --dir_cache $VEP_CACHE \
#--fork 9 \
#--dir_plugins /sw/data/vep/107/Plugins \
#-i $1 \
#-o ${1}_vep_annotated.vcf.gz \
#--compress_output bgzip \
#--force_overwrite \
#--assembly GRCh38 \
#--vcf --check_existing --variant_class \
#--af --af_gnomade --af_gnomadg --max_af

# get the file /Users/xiyas/V2_Genome_reporting/V2_database_prepare/vep_traits.txt_vep_annotated.vcf.gz
library(data.table)
trait_vcf <- fread("/Users/xiyas/V2_Genome_reporting/V2_database_prepare/vep_traits.txt_vep_annotated.vcf.gz")
colnames(trait_vcf)[3] <- "variants"

trait_db <- fread("/Users/xiyas/V2_Genome_reporting/database-file-v2/Reports_genome_databases_traits.txt",header = TRUE)
trait_db <- trait_db[1:77,]
Merged_trait_db <- merge(trait_db,trait_vcf,by= "variants",all.x = TRUE)


