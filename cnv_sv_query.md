# Analysis of CNV and SV files for Patient P00422

This document summarizes the analysis of the `.cnv.vcf` and `.sv.vcf` files to verify the presence of two reported structural variants.

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
