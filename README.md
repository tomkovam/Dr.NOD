# Dr.NOD
Computational framework for **D**iscovery of **R**egulatory **NO**n-coding **D**rivers

The framework detects non-coding cis-regulatory candidate driver mutations that are predicted to regulate gene expression, using tissue-matched enhancer-gene maps, paired tumour whole-genome sequencing and RNA-sequencing data, and additional genomic and epigenomic data sets. First, we build the background mutagenesis model, then we define regulatory space of every gene expressed in the given tissue, and ultimately we search for genes that are upregulated or downregulated in samples that have mutations in the regulatory space and where the number of mutations in the regulatory space exceeds the expected numbers based on the background mutagenesis model.

The entire framework runs from scriptRunAnalysis.m. It will use the precomputed files from directory save/data. For running from scratch, all the input data files are described in inputParameters.properties.

Requirements: Matlab 2022a, Statistics and Machine Learning Toolbox, Bioinformatics Toolbox
