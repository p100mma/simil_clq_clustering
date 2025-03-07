This respository corresponds to the manuscript submission for ISMB 2025:

# Versatile clustering of molecular data in diverse scenarios by CliqueSimNet

## Authors
 Piotr Stomma<sup>1</sup>, Sajad Shahbazi<sup>2</sup>, Krzysztof Mnich<sup>2</sup>, Aneta Polewko-Klim<sup>1</sup> and Witold R. Rudnicki<sup>1,2</sup>

 1. Faculty of Computer Science, University of Bialystok
 2. Computational Centre, University of Białystok

# Structure of the repository

This file contains information on input data and implementation of our numerical experiments.

Rendered html of Rmd notebooks in `docs` contain additional docummentation and some extended results. In particular:

1. [CliqueSimNet_overview](https://htmlpreview.github.io/?https://github.com/p100mma/simil_clq_clustering/blob/main/docs/CliqueSimNet_overview.html) contains additional algorithmic details and some motivating examples of properties of cliques found by our methods.
2. [BRCA_initial_selection](https://htmlpreview.github.io/?https://github.com/p100mma/simil_clq_clustering/blob/main/docs/BRCA_initial_selection.html) contains summary of our initial comparison of all different variants on gene coexpression networks.
3. [DBscan_parameters](https://htmlpreview.github.io/?https://github.com/p100mma/simil_clq_clustering/blob/main/docs/DBscan_parameters.html) shows how the `eps` parameter for DB-scan was choosen for all the datasets on which it was tested.
4. [Leuk_18_tables](https://htmlpreview.github.io/?https://github.com/p100mma/simil_clq_clustering/blob/main/docs/Leuk_18_full_ARI_table.html) displays extended table of Adjusted Rand Index of all 14 tested methods under all applicable thresholding scenarios. 

`.rds` file contains one of used datasets (another one is on zenodo repository, see below). Provided `.R` scripts reproduce our results fully from theese files. 

Other subdirectories  are set up to cotain results of scripts that were ran over all combinations/random trials:
 - `B_200`, `B_1394` results from multiple random trials in variable clustering
 - `leuk_UM`, `leuk_DB_scan_kNN_curve` for storing files for initial tuning of DB-scan on `Leuk_18` dataset
 -  `leuk_clusters` results from multiple random trials in sample clustering.

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

Installation of each package and its usage (usage for general purposes) is docummented in respective repositories. This repository contains only instructions on how to reproduce our experimental results and additional information about methods. 

# How to run our code

## Downloading the data and code

1. clone the github repository:

```
git clone https://github.com/p100mma/simil_clq_clustering
```

2. Inside the main folder of the repository, download 2 dependencies (CliqueSimNet implementation):

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

- `Dock_R_clique` - R with CliqueSimNet and dependencies
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

Scripts used with either of the docker images are differentiated by names.

`v` suffix refers to scripts for variable clustering problem, `s` -- for sample clustering.

Other scripts contain internal utility functions.

Running each test (after opening up microenvironment) can be achieved by using `Rscript` (with `<ARGS>` if they are necessary) :

```
Rscript --no-save <SCRIPT-NAME> <ARG>
```

### Variable clustering

First, open up the correct docker image:
```
sudo docker run -it -v.:/home/ismb_25 r_var_cl
```
Afterwards, run each script in the following order (all subsections) to get all of our results:

#### Initial run on whole dataset

- `v1a_similarity_graph.R` -- computes initial correlation network for testing different variants of CliqueSimNet, outputs `BRCA_similarity_graph.rds`.
- `v1b_initial_algo_test.R` -- run with integer argument ranging from 1 to 48 to test combinations listed in [here](https://htmlpreview.github.io/?https://github.com/p100mma/simil_clq_clustering/blob/main/docs/BRCA_initial_selection.html). We ran this part in parallel on computational cluster setting min. memory requirements to 30 GB.
```
Rscript --no-save v1b_initial_algo_test.R 1
Rscript --no-save v1b_initial_algo_test.R 2
...
Rscript --no-save v1b_initial_algo_test.R 48
```
This script outputs `.rds` files named like
```
clusters_t_S=<t_S CHOICE VARIANT>;t_CS=<OBJECTIVE FUNCTION>;mode=<CLIQUE EXPANSION MODE>;join=<CLIQUE RELAXATION TYPE>.rds  
```
each one containing results per each variant. 
- `v1c_initial_summary.R` - aggregates results over all combinations tested in `v1b_`, outputs `devel_method_runtime_df.rds`, `devel_method_similarity.rds`, `devel_label_vecs.rds`.
- `v1d_viable_method_selection.R` - chooses 6 methods we presented based on pairwise similarities between all methods, outputs `vcl_clq_based_variants_to_test.rds` (listing of parameters of variants to test further) and `clq_reference_modules.rds` (corresponding already computed cluster labels in `v1b`. ) . See [here](https://htmlpreview.github.io/?https://github.com/p100mma/simil_clq_clustering/blob/main/docs/BRCA_initial_selection.html)  for visualization of results (contains the same code with textual description).

#### initial WGCNA modules and heatmap comparison

- `v1e_WGCNA.R` - computes clusters for 3 different `minSize` parameters, saves results in `WGCNA_labels_fulldata.rds`.
- `v1f_heatmaps.R <clq_nr> <hmap_size>` - produces gene module heatmap of size `<hmap_size>` visualization for 1 of 6 variants of CliqueSimNet tested in article given by `<clq_nr>` argument (and for WGCNA `minSize=30` ). To get the 2 heatmaps from the article, assuming all previous steps have been completed:
```
Rscript --no-save v1f_heatmaps.R 4 500
```
The script will output 2 pdf files: `4clq_S.pdf`-- heatmap of CliqueSimNet and `1_WGCNAS.pdf` -- of WGCNA.

#### Stability tests by bootstrap

-`v2a_prepare_bootstrap_sets.R` - prepares indexes of samples for 100 repeats bootstrap testing. Outputs 2 lists of indexes corresponding to 2 different bootstrap set sizes: `v_bs_N.rds` (N=1394), `v_bs_200.rds` (N=200).

-`v2b_stability_resample.R <BS_sample_number>` - run batch of tests for all methods for resample number `<BS_sample_number>`. To complete all trials, this would have to be run 100 times with running argument: 
```
Rscript --no-save v2b_stability_resample.R 1
Rscript --no-save v2b_stability_resample.R 2
...
Rscript --no-save v2b_stability_resample.R 100
```
(we ran this in parallel on HPC cluster, setting minimal RAM requirements to 30 GB). Results are saved in subdirectories `B_200`, `B_1394`. Pattern for filename of outputs from CliqueSimNet is the same as in `v1b` script but contains `<BS_sample_number>` prefix, while output for WGCNA is named `<BS_sample_number>_WGCNA_labels.rds`.

-`v2c_aggregate_resamples.R, v2d_aggregate_cores.R` - computes stability based on resamples from previous step, second script limits this to the `cores` of CliqueSimNet ([see here](https://htmlpreview.github.io/?https://github.com/p100mma/simil_clq_clustering/blob/main/docs/CliqueSimNet_overview.html)  for exact definition). Outputs are `v2c_aggregated_stability_results.rds` and `v2d_aggregated_stability_cores.rds` and `v2d_2346_cores.rds` (last file is list of core membership vectors for all of the variants presented in the article computed on original data, for reference).

-`v2e_summarise.R` produces plots:

   - stability based on output of previous 2 scripts -- `v2e_stability.pdf`
 
   - rest of plots in fig. 4 -- `v2e_modularity.pdf`, `v2e_pairwise_sim.pdf`


### Sample clustering

Remember that this test requires downloading Leuk_18 dataset to the main directory:

```
curl --output trainLeukemia.RData https://zenodo.org/records/14729079/files/trainLeukemia.RData
```


First, open up the correct docker image:
```
sudo docker run -it -v.:/home/ismb_25 r_sub_cl
```

Afterwards, run each script in the following order (all subsections) to get all of our results:


#### Initial estimation of parameters of DB-scan on Leuk_18

Running scripts there is not strictly necessary but output provides a rationale for our selection of crucial parameters of DB-scan: `eps` and `minPts`. We provide a shorter visual summary in the [RMD FILE linked here](https://htmlpreview.github.io/?https://github.com/p100mma/simil_clq_clustering/blob/main/docs/DBscan_parameters.html):

- `s1a_initial_UMAP.R <n_genes>` computes and saves initial 3D UMAP embeddings for estimating parameters of DB-scan. Should be run with running argument from 1 to 12 (each integer corresponds to different number of genes used to get the low. dim. space):
```
Rscript --no-save s1a_initial_UMAP.R 1
Rscript --no-save s1a_initial_UMAP.R 2
...
Rscript --no-save s1a_initial_UMAP.R 12
```
Output is stored in files `leuk_UM/n_genes=<n_genes>_UMkD.rds`. 
- `s1b_DB_scan_param_est.R <n_genes>` - calulates kNN distances of each point (k:=3*2-1) for the purpose of picking optimal `eps` by "elbow point" method. This script should be run with running argument exactly like the one above. Results are stored in the files `leuk_DB_scan_kNN_curve/n_genes=<n_genes>_kNN_d_curve.rds`. Additionaly, it produces plots of kNN distances pear each gene number used which are saved as pdf in the same directory.
- `s1c_joint_kNN_plot.R` - uses output of the previous script to prepare a joint plot of all kNN distances (from all gene number cutoffs) and saves in a pdf `leuk_DB_scan_kNN_curve/joint_plot.pdf`. All curves (representing sorted kNN distances in UMAP embedding) are plotted on the same axes, to show that the `eps=0.3` value used is universal for all `n_genes` variants.   

#### 100 trials of clustering by different methods 

- `s2a_clustering.R <n_genes> <trial_number>` - calculates 3D UMAP for gene number `<n_genes>` and trial number `<trial_number>`. This was ran in parallel on our computational cluster (min. memory requirement set to 30 GB). Should be ran with `<n_genes>` running from 1 to 12 and `<trial_number>` running from 1 to 100:
```
Rscript --no-save s2a_clustering.R 1 1
...
Rscript --no-save s2a_clustering.R 1 100
Rscript --no-save s2a_clustering.R 2 1
...
Rscript --no-save s2a_clustering.R 2 100
...
...
Rscript --no-save s2a_clustering.R 12 100
```
Outputs are saved in files `leuk_clusters/n_genes=<n_genes>_trial=<trial_number>_clusters.rds`.
- `s2b_cl_summary.R`- aggregates results from all trials ran from previous script. Outputs a file `leuk_clusters_results_aggregated.rds` containing ARIs computed between method-derived and reference clusters for each method, thresholding strategy and gene number used + numerical values of thresholds.
- `s2c_tables_plots.R` - based on the results of previous scripts above, produces tables of ARI scores on Leuk_18 presented in the publication, along with the plot of comparison of DB-scan with clique based clustering (see [HERE](https://htmlpreview.github.io/?https://github.com/p100mma/simil_clq_clustering/blob/main/docs/Leuk_18_full_ARI_table.html) for full table over all combinations of thresholds and methods). List of outputs:
    - "Raw" tables are contained in files `similarity_based_methods_publication_table.csv` and `metric_based_methods_publication_table.csv`,
    - while their versions formatted for easy inclusion into `.tex` file are contained in files: `s_compARI.csv` and `m_compARI.csv`.
    - Plot is in the generated PDF:  ` ARI_clq_dbscan.pdf`.
    - Additionally: script generates two tables that show per each gene number variant and algorithm combination, which thresholding strategy was the best according to IQR of ARI over 100 trials (NA values indicate no clear winner): `simil_based_best_thr_per_algVn_genes.csv`, `clique_based_best_thr_per_algVn_genes.csv`. 

#### Additional FCPS benchmark tests
- `s3a_dbscan_FCPS.R` - run of DB-scan on chosen benchmark data, `eps` chosen according to method described in depth [here](https://htmlpreview.github.io/?https://github.com/p100mma/simil_clq_clustering/blob/main/docs/DBscan_parameters.html). Outputs `DBscan_FCPS_scores.rds`.
- `s3b_cliques_FCPS.R` - tests of CliqueSimNet on FCPS benchmark data. Outputs `cliquesOnFCPS.rds`.
- `s3c_fcps_summary.R` - aggregate results from above tests of clique based clustering and DB-scan into tables and generate visualization plots of the point cloud data. Generates `fcps_results.pdf` and `FCPS_pointClouds.pdf`.
