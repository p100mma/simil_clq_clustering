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

implementation of CliqueSimNet is shared in form of 2 R packages which are held frozen for reproduciblity reasons:

- [cliqueClusteR](https://github.com/p100mma/cliqueClusteR)
- [cliquePartitoneR](https://github.com/p100mma/cliquePartitioneR)

Installation of each package and its usage is docummented in respective repositories. This repository contains only instructions on how to reproduce our experimental results and additional information about methods. 

# How to run our code

## Downloading the data and code

1. clone the github repository:

```
git clone https://github.com/p100mma/simil_clq_clustering
```

2. Inside the main folder of the repository, download 2 dependencies (CliqueNetSim implementation):

```
cd simil_clq_clustering
git clone https://github.com/p100mma/cliquePartitioneR
git clone https://github.com/p100mma/cliqueClusteR
```

3. For launching tests on  `Leuk_18` dataset, download preprocessed data from [Zenodo](https://doi.org/10.5281/zenodo.14729079) (note the file size requirements: 841 Mb), either manually (and put in the main directory) or by command below (launched inside main folder of respository): 

```
curl --output trainLeukemia.RData https://zenodo.org/records/14729079/files/trainLeukemia.RData
```
## Set un the environment using docker (4.66 GB total)

There are 3 dockerfiles:

- `Dock_R_clique` - R with CliqueNetSim and dependencies
- `Dock_var_clusters` - dependencies for variable clustering
- `Dock_subject_clusters` - dependencies for sample clustering

Contents and package versions contained there are listed at the end of this file.

All of them should be built from the main folder. 
Second and third dockerfile depends on the first.
This will create 3 docker images named `r_cluq, r_var_cl, r_sub_cl`.
```
sudo docker build -f Dock_R_clique -t r_cliq .   # for the package only
sudo docker build -f Dock_var_clusters -t r_var_cl .   # for variable clustering tests
sudo docker build -f Dock_subject_clusters -t r_sub_cl .   # for subject clustering
```

Afterwards, running `docker run` will open up microenvironment in which it is possible to run all the tests.

```
sudo docker run -it -v .:/home/ismb_25 <IMAGE_NAME>
```

Note that for practical reasons (memory limiatations and parallelism of resampling tests), we have run all of the tests
using our computational cluster and [singularity](https://github.com/sylabs/singularity).
One can convert docker images to singularity ones by [docker2singularity](https://github.com/singularityhub/docker2singularity) tool.
We list example commands to convert one of our docker images to singularity below  (creates file in `/tmp/test/` directory):

```
sudo docker run -v /var/run/docker.sock:/var/run/docker.sock \
-v /tmp/test:/output --privileged -t \
 --rm quay.io/singularity/docker2singularity r_sub_cl
```

## Running the computational tests

Scripts starting with `devel_` prefix are meant to be run before the main, computationally heavy tasks.

`v` suffix refers to scripts for variable clustering problem, `s` -- for sample clustering.


### Variable clustering



