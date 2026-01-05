[](https://opensource.org/licenses/MIT)

# 3D mRNA Optimization

This repertory narrates a three-dimensional (3D) mRNA optimization algorithm incorporating codon adaptation index (CAI), minimum free energy (MFE), and most notably, codon pair bias (CPB) for enhanced vaccine potency. The algorithm is implemented in the paper:

> **A Three-Dimensional mRNA Optimization Strategy Incorporating Codon Pair Bias Enhances Vaccine Potency**  

## 🧬 Overview

Previously, only the codon optimization strategy focusing on CAI had been intensively adopted in the vaccine sequence design, and merely a few studies included the consideration of secondary structure of RNA (which could be, at least partially represented by MFE). Here, we propose a novel 3D algorithm which takes a step further to introduce **codon pair bias (CPB)** as a third critical dimension, enabling the generation of mRNA sequences with enhanced translational characteristics and improved vaccine immunogenicity. Candidate vaccines designed by this algorithm are likely to trigger more robust immonogeneity and provide more well-rounded protection against infection comparing with prior 1D and 2D optimization. 

## ✨ Key Features

- **Three-dimensional optimization** of CAI, CPB, and MFE parameters
- **Iterative refinement algorithm** for balanced sequence optimization
- **Compatible with standard bioinformatics tools** (RNAfold, EMBOSS CAI)
- **Modular design** for easy integration into existing pipelines
- **Comprehensive validation** against SARS-CoV-2 Spike protein sequences

## 📁 File Description
- src: Source code modules  
- codonpair_optimal.py: Main optimization algorithm
- cpb_calculation.py: CPB score computation
- data: Data files
- CPS_huamn.csv: Codon Pair Score (CPS) for humans 
- replace_codonpair.csv: High-frequency codon pair list
- examples: Usage examples
- SARS2_S_WT.fasta: Example of SARS-CoV-2 prototype spike protein sequence
- requirements.txt: Python dependencies

## 🔧 Modules Description
**1. Codon Pair Optimization** (src/Codonpair_Optimal.py)

Implements the core iterative optimization algorithm that:
- Starts with LinearDesign-generated baseline sequences
- Systematically replaces low-frequency codon pairs
- Maintains optimal CAI and MFE values during optimization

**2. CPB Calculation** (src/Cpb_Calculation.py)

Computes Codon Pair Bias scores based on:
- Human genome codon pair frequency statistics
- Normalized CPB scoring methodology
- Batch processing for multiple sequences

## 📊 External Tool Integration

The algorithm integrates with established bioinformatics tools:
- **MFE Calculation**: RNAfold Web Server
- **CAI Calculation**: EMBOSS CAI Tool
- **Sequence Analysis**: BioPython compatibility

## 🧪 Validation Data

The algorithm has been validated using:
- SARS-CoV-2 Spike protein (prototype and EG.5 variant)
- _In vitro_ protein expression assays
- _In vivo_ immunogenicity studies in mouse models

## 📋 Dependencies
- Python 3.8+
- BioPython >= 1.79
- pandas >= 1.3.0
- NumPy >= 1.21.0

See **requirements.txt** for complete list.

## 🤝 Contributing
We welcome contributions! Please see our **Contributing Guidelines** for details.

## 📜 License
This project is licensed under the MIT License - see the **LICENSE** file for details.

## 📞 Contact
**Gong Cheng** (Corresponding Author) - gongcheng@mail.tsinghua.edu.cn
