# Genomics_platform

A repository for the preparing manuscript "A platform for genomic data analysis, reporting, and interpretation"

## Overview

------------------------------------------------------------------------


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


| Gene     | Variant            | Type       | Zygosity | Individuals | Frequency |
|----------|--------------------|------------|-----------|-------------|-----------|
| TTN      | rs12345, c.123A>T  | ClinVar P  | HET       | 8           | 2.9%      |
| LDLR     | rs67890, c.456C>G  | pLoF       | HOM       | 1           | 0.4%      |
| ABCA4    | rs54321, c.789G>A  | ClinVar LP | HET       | 3           | 1.1%      |
| ...      | ...                | ...        | ...       | ...         | ...       |

------------------------------------------------------------------------


Disclaimer: BTyper3 is pretty neat! However, no tool is perfect, and BTyper3 cannot definitively prove whether an isolate is pathogenic or not. As always, interpret your results with caution. We are not responsible for taxonomic misclassifications, misclassifications of an isolate's pathogenic potential or industrial utility, and/or misinterpretations (biological, statistical, or otherwise) of BTyper3 results.
