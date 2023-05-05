---
title: "cablogrep"
author: "Dániel Gerber"
---

## Cablogrep

### Description

*Cablogrep* is an easy-to-use tool written in R and bash programming languages for inferring mitochondrial haplogroup for horses based on the haplogroup system of [Achilli et al. 2011](https://www.pnas.org/doi/10.1073/pnas.1111637109). The system is updated with haplogroups X and S based on [Der Sarkissian et al. 2015](https://www.cell.com/current-biology/pdfExtended/S0960-9822(15)01003-9) and [Guimaraes et al. 2021](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7494339/), respectively.

### Dependencies

*Cablogrep* is scripted in R 4.2.3, but does not require extra packages, thus it should run without issues on any basic version. The input for the R script a *NUCmer v3.1* output file, thus this program needs to be installed (it is part of the software package *MUMMER*, and can be installed from Ubuntu repository).

### Usage

*Cablogrep* does not need installation, it can be directly run locally with the following options:

    ./cablogrep.sh -p /path/to/fasta_files -o output_name

### Warnings, advices & recommendations

*FASTA* files must have **.fasta** extension

Sample names are going to be the names of the fasta header, not the actual filename, therefore it is recommended to check these

Due to the ambiguity of extant classifications (high number of hotspots, low sample number, etc.), we discarded known hotspots from the reference panel, and offer haplogroup quality and second best guesses as output. Despite these, we observed certain ambiguities in classifications, thus we recommend to check manually haplogroups that has lower quality than 0.6 and/or has non-cladal first and second best hits. For example, we observed that maybe due to previously undescribed hotspots, haplogroup I tend to appear as first best hit for samples that truly belong within haplogroup A'D. In future versions, we are planning to mitigate this issue.

### Output

*Cablogrep* provides two outputs: .HG and .HgCandidates. The former has the sample name and the first and second best hits with haplogroup quality scores, which then can be concatenated to a single table. The latter offers all the possible hits with quality scores, if anyone needed that.

### Cite

Dicső Z., et al. (2023) A genetic study on horsekeeping during the Bronze Age Carpathian Basin
