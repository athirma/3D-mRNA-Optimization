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

