# CLIER

This repository contains code and test data for the article [**"Protocol for interpretable and context-specific single-cell informed deconvolution of bulk RNA-seq data,"**](https://www.sciencedirect.com/science/article/pii/S2666166725000760?via%3Dihub) published in STAR Protocols (Malpetti, Mangili et al., 2025).

## Detailed Description of Files

### Code
- **protocol_code.R**: Contains the lines of code included in the manuscript (excluding the processing from FASTQ to TPM).
- **aux_functions.R**: Contains all the R functions necessary for executing the protocol.
- **align_fastq.sh**: Automates the process of downloading FASTQ files and aligning paired-end RNA-Seq data using the STAR aligner.

### Data
- **kidney_atlas_matrix.rds**: Contains the single-cell signatures atlas built in ["A transfer learning framework to elucidate the clinical relevance of altered proximal tubule cell states in kidney disease,"](https://www.sciencedirect.com/science/article/pii/S2589004224004929) published in iScience (Legouis, Rinaldi, Malpetti et al., 2024).
- **kidney_atlas_info.xlsx**: Contains descriptions of the signatures included in the single-cell signatures atlas built in the above study.
- **DKD_tpm.rds**: Contains a processed version (TPM) of the dataset GSE142025, also used in the above study.
- **DKD_clin.rds**: Contains clinical information (fibrosis) regarding the dataset GSE142025.
- **genelength.txt**: Contains genes length (to be used in data processing).

## On the Execution

The code in **protocol_code.R** can be fully executed using the test data provided in this repository. Users who might want to skip the training phase (that takes approximately 9 hours) and test a pre-trained model can find the KCLIER model [here](https://drive.switch.ch/index.php/s/OpvMh1vGRgRmKKf), together with other intermediate files produced during the execution. We share these file separately since, given their large size, they cannot fit on GitHub.
