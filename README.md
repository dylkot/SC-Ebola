# sc-Ebola
Analysis scripts used in the analysis of Seq-Well and CyTOF datasets of Ebola virus infected PBMCs for the manuscript [Single-cell profiling of Ebola virus infection in vivo reveals viral and host transcriptional dynamics](https://www.sciencedirect.com/science/article/pii/S0092867420313088). See https://github.com/dylkot/dropSeqPipe-dak for the pipeline used to map the  single-cell RNA-Seq reads and to generate the expression count matrices.

# Directory structure
- setup.sh          - commands to create the conda environment used for all of the analysis
- Code              - libraries used throughout the analysis
- Analysis          - Jupyter notebooks for running the pipelines and generating the figures in the manuscript
	- InVivo          - Analysis of Seq-Well and CyTOF data from Ebola infected rhesus macaques
	- ExVivo          - Analysis of Seq-Well data from PBMCs extracted from rhesus macaques and infected ex vivo
	- HCA_BoneMarrow  - Comparison of EVD data with healthy human bone marrow and PBMCs
	- Human_EVD_CyTOF - Re-analysis of published CyTOF data of human EVD cases (McElroy et al., 2020)
