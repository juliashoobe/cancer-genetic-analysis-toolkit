# Cancer Genetic Analysis Toolkit

The **Cancer Genetic Analysis Toolkit** is a Python-based software suite designed to assist cancer researchers analyze genetic data.
This toolkit includes functions to explore and understand the role of mutations and gene expression in cancer. 


## Features
### 1. **Mutation Detection**
- **Description**: Identifies mutations in cancer-related genes by comparing a target DNA sequence with a reference sequence.

### 2. **Gene Expression Analysis**
- **Description**: Analyzes gene expresion data to identify differentially expressed genes in cancerous tissues compared to normal tissues.

### 3. **Mutation-Driven Differential Gene Expression**
- **Description**: Examines whether mutations in specific genes are associated with changes in their expression levels in cancerous versus nromal tissues. 

### 4.**Command-Line Interface (CLI)**
- **Description**: A command-line interface that allows users to run the different analyses and view results direclty from the terminal.


## Software Architecture
- **Main File**: `toolkit.py` contains the main functions for mutation detection, gene expression analysis, and mutation-driven differential gene expression.
- **Tests**: Test files will be separated into a directory `tests/` for unit tests on each feature.
- **CLI**: Uses `argparse` to handle commands and arguments for running different analyses via terminal.
