# Computational investigation of the sequence context of arginine/glycine-rich motifs in the human proteome

This repository contains code and data used in the publication:

**"Computational investigation of the sequence context of arginine/glycine-rich motifs in the human proteome"**  
*Eric Schumbera, Dorothee Dormann, Andreas Walther, Miguel A. Andrade-Navarro*
submitted to BMC Genomics 

## Citation

If you use this repository, please cite our work:
INFORMATION WILL BE ADDED AFTER REVIEW

## Repository structure
```
├── data/                 # Data stored necessary for this project
│ ├── external            # External data
│ ├── processed           # Processed files necessary for the statistical analysis
│ ├── results             # Output files and figures
│ │ ├── subfigures        # Subfigures necessary to recreate the excat figures from the publication
├── 0_preprocessing_of_proteome.ipynb     # Code for processing the raw data and generating the necessary datasets
├── 1_general_analyses.ipynb             # Script for statistical analysis of phys/chem properties of RG motifs and general properties (length, impurity, etc.)
├── 2_IDR_analyses.ipynb                 # Scripts for statistical analysis of IDR-related properties
├── 3_domain_analyses.ipynb              # Scripts for statistical analysis of domain-focused analysis
├── 4_AA_analyses.ipynb                  # Scripts for statistical analysis of detailed amino acid composition
├── LICENSE               # License file
├── environment.yml       # Necessary dependencies for project
└── README.md             # This file
```
## Requirements and Installation

All dependencies for this project are listed in the provided `environment.yml` file.

### Recommended (Using Conda)

Create the environment:

```bash conda env create -f environment.yml ```

### Activate the environment:

```bash conda activate your-environment-name ```

## How to use:

To recreate this analysis, the preprocessing `0_preprocessing_of_proteome.ipynb` is necessary to create/clean/annotate the necessary datasets. The statistical analyses (1-4) can be run in independent order and edited/expanded independently from each other.

## License:

see `LICENSE` file

## Contact

For questions or suggestions, please contact:
Eric Schumbera
e.schumbera@uni-mainz.de