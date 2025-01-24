This respository corresponds to the manuscript submission for ISMB 2025:

# Versatile clustering of molecular data in diverse scenarios by CliqueSimNet

## Authors
 Piotr Stomma<sup>1</sup>, Sajad Shahbazi<sup>2</sup>, Krzysztof Mnich<sup>2</sup>, Aneta Polewko-Klim<sup>1</sup> and Witold R. Rudnicki<sup>1,2</sup>

 1. Faculty of Computer Science, University of Bialystok
 2. Computational Centre, University of Białystok

# Structure of the repository

This file contains information on input data and implementation of our numerical experiments.
Rendered PDFs of Rmd notebooks contain additional docummentation and some extended results. 

## Data information

### Variable clustering test

`gene_expr_data.rds` - Gene expression profiles of 1394 breast cancer patients, 8673 genes (European Genome-Phenome Archive, accesion ID: EGAS00000000083). References:

- Pereira B et al. The somatic mutation profiles of 2, 433 breast cancers refine their genomic and transcriptomic landscapes. Nature Communications, 2016; 7.1 Preprocessed and prefiltered as in:
- (preprocessing) Polewko-Klim A, Mnich K, Rudnicki W, Robust Data Integration Method for Classification of Biomedical Data. Journal ofMedical Systems, 2021; 45.

### Molecular tumor subtypes problem (Leuk_18)

Preprocessed data is available on zenodo: [https://doi.org/10.1111/j.1365-2141.2008.07261](https://doi.org/10.5281/zenodo.14729079)

`trainLeukemia.RData` - quantile normalized gene expression data (Affymetrix HG-U133 Plus 2.0 array) from the MILE (Microarray Innovations In LEukemia) program<sup>1</sup> (ArrayExpress, ID:GEOD-13159)
Data pre-processing and all analyses were conducted using the open-source statistical software R v3.4.3 (R Core Team, 2017a).
Raw microarray data from stage I of the MILE study included 2,096 .CEL files that contained 16 acute and chronic leukemia subclasses, myelodysplastic syndromes (MDSs), 
and “none of the target classes” control group that included non-leukemia and healthy bone marrow. 

Data preparation included two main steps, namely cleaning and transformation. To perform a quality assessment for Affymetrix GeneChip data, we using the R package affy. 
The homogeneous distribution of signal and intensities and assessment of RNA degradation were checked. The samples were removed with dataset if their value of mean signal intensities was different than:
$\pm3 \sigma$ and if threshold value for the 3. /5 ratios more then 3 for $\beta$-actin, and 1.25 for GAPDH (criterium recommended by Affymetrix).
Next, we calculated the normalized values with robust rma background correction, quantile normalization, and median polish. To achieve a normal or near-normal distribution, 
the log2 transformation of gene expression data was performed . 

After this transformation of the data, the data set (2,077 samples and 54675 probes) was subjected to unsupervised learning. 
References:

1. Kohlmann A, Kipps TJ, Rassenti LZ, Downing JR, Shurtleff SA, Mills KI, Gilkes AF, Hofmann WK, Basso G, Dell'orto MC, Foà R, Chiaretti S, De Vos J, Rauhut S, Papenhausen PR, Hernández JM, Lumbreras E, Yeoh AE, Koay ES, Li R, Liu WM, Williams PM, Wieczorek L, Haferlach T. An international standardization programme towards the application of gene expression profiling in routine leukaemia diagnostics: the Microarray Innovations in LEukemia study prephase. Br J Haematol. 2008 Sep;142(5):802-7. doi: 10.1111/j.1365-2141.2008.07261.x. PMID: 18573112; PMCID: PMC2654477.

# Implementation details

implementation of CliqueSimNet is shared in form of 2 R packages which are included also in this repository for reproduciblity reasons:

- [cliqueClusteR](https://github.com/p100mma/cliqueClusteR)
- [cliquePartitoneR](https://github.com/p100mma/cliquePartitioneR)
- 
Installation of each package and its usage is docummented in respective repositories. This repository contains only instructions on how to reproduce our experimental results and additional information about methods. 

How tp run using docker
  
