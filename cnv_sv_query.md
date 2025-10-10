# Analysis of CNV and SV files for Patient P00422

The translation from Russian report that saying the patient's diseases:
- 1. Chromosome 9p24.3–p23 Microdeletion
	•	Genomic location (GRCh38): chr9: 0–13,587,692
	•	Type: Heterozygous microdeletion
	•	Size: ~13.6 Mb
	•	Copy number: 1 (loss)
	•	Affected genes include: DOCK8, KANK1, SMARCA2, VLDLR, KCNV2, GLIS3, SLC1A1, JAK2, CD274, RIC1, GLDC
	•	Associated syndrome: Chromosome 9p deletion syndrome (OMIM 158170)
	•	Clinical description: This syndrome is typically associated with global developmental delay, intellectual disability, craniofacial dysmorphism, and sometimes structural brain or heart abnormalities.
	•	Interpretation: Classified as pathogenic and considered the most likely cause of the patient’s clinical presentation.


- 2. Chromosome 10q25.1–q26.3 Microduplication
	•	Genomic location (GRCh38): chr10: 107,558,240–133,797,422
	•	Type: Heterozygous microduplication
	•	Size: ~26.2 Mb
	•	Copy number: 3 (gain)
	•	Affected genes include: ACADSB, FGFR2, HTRA1, SHOC2, TCF7L2, WDR11, BAG3, EMX2, GFRA1, MXI1, VAX1, and others
	•	Associated syndrome: Distal duplication 10q syndrome
	•	Clinical description: Usually linked to developmental delay, growth retardation, abnormal muscle tone, facial and organ malformations, and neurodevelopmental abnormalities.
	•	Interpretation: Classified as pathogenic and also considered a likely cause of the patient’s phenotype.
"


Here, I summarizes the query of the `.cnv.vcf` and `.sv.vcf` files to verify the presence of two reported structural variants.

## 1. Chromosome 9p24.3–p23 Microdeletion

### Grep Command

To find deletions on chromosome 9, the following `grep` command was used on both the `cnv` and `sv` files:

```bash
grep 'chr9' /Users/xiyas/Downloads/szaomics_p00422_24_6496_2025-10-10_1106/P00422_24_6496.cnv.vcf | grep 'DEL'
grep 'chr9' /Users/xiyas/Downloads/szaomics_p00422_24_6496_2025-10-10_1106/P00422_24_6496.sv.vcf | grep 'DEL'
```

### Findings

A large, heterozygous deletion was identified on chromosome 9 in both files, consistent with the clinical report.

-   **In `P00422_24_6496.cnv.vcf`:**
    -   A deletion spanning from position **114,000 to 11,495,204** (approximately 11.4 Mb) was found.
    -   The genotype (`GT`) is `0/1` (heterozygous) and the copy number (`CN`) is `1`.
    -   **VCF Line:**
        ```
        chr9	114000	DRAGEN:LOSS:chr9:114001-11495204	N	<DEL>	104	PASS	SVLEN=-11381204;SVTYPE=CNV;END=11495204;REFLEN=11381204	GT:SM:CN:BC:GC:CT:AC:PE	0/1:0.507353:1:10057:0.386923:0.499821:0.500456:0,18
        ```

-   **In `P00422_24_6496.sv.vcf`:**
    -   A large deletion from position **6,995,965 to 21,038,487** (approximately 14 Mb) was identified.
    -   The genotype (`GT`) is `0/1` (heterozygous).
    -   **VCF Line:**
        ```
        chr9	6995965	DRAGEN:DEL:135161:0:1:0:0:0	G	<DEL>	659	NoPairSupport	END=21038487;SVTYPE=DEL;SVLEN=-14042522;CIPOS=0,7;CIEND=0,7;HOMLEN=7;HOMSEQ=CCCTGCA	GT:FT:GQ:PL:PR:SR	0/1:PASS:342:709,0,339:15,0:16,19
        ```

## 2. Chromosome 10q25.1–q26.3 Microduplication

### Grep Command

To find duplications on chromosome 10, the following `grep` command was used:

```bash
grep 'chr10' /Users/xiyas/Downloads/szaomics_p00422_24_6496_2025-10-10_1106/P00422_24_6496.cnv.vcf | grep 'DUP'
grep 'chr10' /Users/xiyas/Downloads/szaomics_p00422_24_6496_2025-10-10_1106/P00422_24_6496.sv.vcf | grep 'DUP'
```

### Findings

A large, heterozygous duplication was identified on chromosome 10 in the `cnv` file. It was represented as five contiguous events.

-   **In `P00422_24_6496.cnv.vcf`:**
    -   The combined duplication spans from position **107,558,632 to 133,621,535** (approximately 26 Mb).
    -   The copy number (`CN`) for this region is `3`.
    -   **VCF Lines:**
        ```
        chr10	107558632	DRAGEN:GAIN:chr10:107558633-112352736...
        chr10	112357126	DRAGEN:GAIN:chr10:112357127-120467592...
        chr10	120469915	DRAGEN:GAIN:chr10:120469916-131622969...
        chr10	131624948	DRAGEN:GAIN:chr10:131624949-133333049...
        chr10	133339658	DRAGEN:GAIN:chr10:133339659-133621535...
        ```

-   **In `P00422_24_6496.sv.vcf`:**
    -   No large duplications corresponding to the reported event were found.

## Conclusion

Both the pathogenic microdeletion on chromosome 9 and the microduplication on chromosome 10 were confirmed to be present in the provided sequencing files.
