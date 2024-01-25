# Process Turkish data -----------------
df_temp_turkish <- read_and_process_files("/Users/xiyas/V2_Genome_reporting/python_output_turkish_275/nodup4_file", "Turkish")
df_temp_turkish <- customize_temp_data(df_temp_turkish)

clinvar_TR <- get_clinvar_variants(df_temp_turkish)
clinvar_TR_unique <- get_unique_SZAID(clinvar_TR)

sillico_pLoFs_TR <- get_other_sillico_pLoFs(df_temp_turkish)
sillico_pLoFs_unique_TR <- get_unique_variants(sillico_pLoFs_TR)

ACMG_TR <- get_ACMG_findings(df_temp_turkish)
#ACMG_TR_unique <- get_unique_SZAID(ACMG_TR)

# clinVar variants 
## before the P_LP filteration
## now sorted_tab_TR = cohort-level unique SZAID clinvar_TR with freq count
sorted_tab_TR <- get_sorted_tab(clinvar_TR, clinvar_TR_unique)

sorted_tab_TR_filter =list(sorted_tab_filter = get_P_LP(sorted_tab_TR$sorted_tab),
                           sorted_tab_old_filter = get_P_LP(sorted_tab_TR$sorted_tab_old))

# sillico variants
sorted_tab_TR_sillico <- get_sorted_tab_sillico(sillico_pLoFs_TR, sillico_pLoFs_unique_TR)

# sillico no need to fileter

# Process Swedish data -----------------
df_temp_swedish <- read_and_process_files("/Users/xiyas/V2_Genome_reporting/python_output_swedish_101/nodup4_file", "Swedish")
df_temp_swedish <- customize_temp_data(df_temp_swedish)

clinvar_SW <- get_clinvar_variants(df_temp_swedish)
clinvar_SW_unique <- get_unique_SZAID(clinvar_SW)

sillico_pLoFs_SW <- get_other_sillico_pLoFs(df_temp_swedish)
sillico_pLoFs_unique_SW <- get_unique_variants(sillico_pLoFs_SW)

ACMG_SW <- get_ACMG_findings(df_temp_swedish)
ACMG_SW_unique <- get_unique_SZAID(ACMG_SW)

# clinVar variants  
## before the P_LP filteration
sorted_tab_SW <- get_sorted_tab(clinvar_SW, clinvar_SW_unique)
## after the P_LP filteration
sorted_tab_SW_filter =list(sorted_tab_filter = get_P_LP(sorted_tab_SW$sorted_tab),
                           sorted_tab_old_filter = get_P_LP(sorted_tab_SW$sorted_tab_old))
# sillico variants
sorted_tab_SW_sillico <- get_sorted_tab_sillico(sillico_pLoFs_SW, sillico_pLoFs_unique_SW)

# Combine the two data frames into one -----------------
combined_df <- rbind(clinvar_TR,clinvar_SW)
combined_sillico_df <- rbind(sillico_pLoFs_TR,sillico_pLoFs_SW)
#write.table(combined_df,file = "combined_df_clinvar.txt",quote = FALSE,col.names = FALSE,sep = "\t")
#combined_df$Population <- factor(combined_df$Population, levels = c("Turkish", "Swedish"))
#clinvar_genes_count_combined <- rbind(clinvar_genes_count_turkish, clinvar_genes_count_swedish)
