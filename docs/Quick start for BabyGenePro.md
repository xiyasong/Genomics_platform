# **Quick start**

## Introduction

BabyGenePro is a freely accessible web portal, aiming to assist the newborn screening processes and support clinical decision-making in neonatal rare disease diagnosis. This web portal streamlines data analysis and generates digital reports for next-generation sequencing data. The platform builds upon the extensively validated Variant Effect Predictor (VEP) framework (McLaren et al., 2016) and incorporates standardized interpretation guidelines from the American College of Medical Genetics and Genomics (ACMG) (Miller et al., 2023). Gene panel curation integrates authoritative sources, including the Recommended Uniform Screening Panel (RUSP) and seven international NBSeq initiatives (Betzler et al., 2024), ensuring comprehensive coverage of clinical practice. The server processes input in variant call format (VCF) files, which contain single nucleotide variants (SNVs), insertions/deletions (INDELs), or structural variants (SVs) using industry-standard annotation tools and databases. Our implementation leverages established clinical genomics resources for variant interpretation, including OMIM (Amberger et al., 2015), ClinVar (Landrum et al., 2014), and GWAS Catalogue (Sollis et al., 2023), while incorporating validated pharmacogenomic databases (Whirl-Carrillo et al., 2021) for paediatric medication guidance. The system generates standardised clinical reports following the ACMG-AMP variant-interpretation guidelines (Richards et al., 2015).

## **Submit jobs**

Users are required to sign up with a personal email address to start analysis for sensitive genome data. BabyGenePro currently only accepts files in the vcf format. Four types of data are available to be analysed in BabyGenePro: 1) single nucleotide polymorphism (SNPs), 2) insertions and deletions (INDELs); 3) copy number variants (CNVs); 4) other structural variants (SVs).

Users can also play with the example file without login. 

![](/Users/xinmengliao/Library/Application%20Support/marktext/images/2025-01-30-16-14-35-image.png)

**References:**

Amberger, J. S., Bocchini, C. A., Schiettecatte, F., Scott, A. F., & Hamosh, A. (2015). [OMIM.org](http://OMIM.org): Online Mendelian Inheritance in Man (OMIM®), an online catalog of human genes and genetic disorders. *Nucleic Acids Res*, *43*(Database issue), D789-798. https://doi.org/10.1093/nar/gku1205

Betzler, I. R., Hempel, M., Mütze, U., Kölker, S., Winkler, E., Dikow, N., Garbade, S. F., Schaaf, C. P., & Brennenstuhl, H. (2024). Comparative analysis of gene and disease selection in genomic newborn screening studies. *J Inherit Metab Dis*, *47*(5), 945-970. https://doi.org/10.1002/jimd.12750

Landrum, M. J., Lee, J. M., Riley, G. R., Jang, W., Rubinstein, W. S., Church, D. M., & Maglott, D. R. (2014). ClinVar: public archive of relationships among sequence variation and human phenotype. *Nucleic Acids Res*, *42*(Database issue), D980-985. [ClinVar: public archive of relationships among sequence variation and human phenotype | Nucleic Acids Research | Oxford Academic](https://doi.org/10.1093/nar/gkt1113)

McLaren, W., Gil, L., Hunt, S. E., Riat, H. S., Ritchie, G. R., Thormann, A., Flicek, P., & Cunningham, F. (2016). The Ensembl Variant Effect Predictor. *Genome Biol*, *17*(1), 122. [The Ensembl Variant Effect Predictor | Genome Biology | Full Text](https://doi.org/10.1186/s13059-016-0974-4)

Miller, D. T., Lee, K., Abul-Husn, N. S., Amendola, L. M., Brothers, K., Chung, W. K., Gollob, M. H., Gordon, A. S., Harrison, S. M., Hershberger, R. E., Klein, T. E., Richards, C. S., Stewart, D. R., & Martin, C. L. (2023). ACMG SF v3.2 list for reporting of secondary findings in clinical exome and genome sequencing: A policy statement of the American College of Medical Genetics and Genomics (ACMG). *Genet Med*, *25*(8), 100866. [Redirecting](https://doi.org/10.1016/j.gim.2023.100866)

Richards, S., Aziz, N., Bale, S., Bick, D., Das, S., Gastier-Foster, J., Grody, W. W., Hegde, M., Lyon, E., Spector, E., Voelkerding, K., & Rehm, H. L. (2015). Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. *Genet Med*, *17*(5), 405-424. [Redirecting](https://doi.org/10.1038/gim.2015.30)

Sollis, E., Mosaku, A., Abid, A., Buniello, A., Cerezo, M., Gil, L., Groza, T., Güneş, O., Hall, P., Hayhurst, J., Ibrahim, A., Ji, Y., John, S., Lewis, E., MacArthur, J. A. L., McMahon, A., Osumi-Sutherland, D., Panoutsopoulou, K., Pendlington, Z.,…Harris, L. W. (2023). The NHGRI-EBI GWAS Catalog: knowledgebase and deposition resource. *Nucleic Acids Res*, *51*(D1), D977-d985. https://doi.org/10.1093/nar/gkac1010

Whirl-Carrillo, M., Huddart, R., Gong, L., Sangkuhl, K., Thorn, C. F., Whaley, R., & Klein, T. E. (2021). An Evidence-Based Framework for Evaluating Pharmacogenomics Knowledge for Personalized Medicine. *Clin Pharmacol Ther*, *110*(3), 563-572. https://doi.org/10.1002/cpt.2350
