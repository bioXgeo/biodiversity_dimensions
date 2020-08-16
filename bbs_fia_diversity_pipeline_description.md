# Bird and tree diversity pipelines used in Read et al. 2020 GEB paper

QDR, 16 Aug 2020

## High-level description

To calculate alpha, beta, and gamma diversity across the three flavors of diversity (taxonomic = TD, phylogenetic = PD, and functional = FD), you need to do the following steps. The subdirectories in the [nasabio github repo](https://github.com/qdread/nasabio) correspond to these steps.

1. Process raw phylogenetic and trait data to cleaned form, including imputation. *scripts in `trait_phylo_data_processing` folder*
2. Create phylogenetic and functional distance matrices, and site-by-species matrices, that are used for the diversity calculations. *scripts in `prep_diversity_files` folder*
3. Actually calculate the alpha, beta, and gamma diversity for points, pairs, and radii respectively, across TD, PD, and FD, then get the average value for each point and radius. *scripts in `run_compile_diversity` folder*

Following are the specific scripts used to execute these three steps for trees and birds.

## Tree (FIA) diversity pipeline

`nasabio/trait_phylo_data_processing/allfiaspp_master_processing.r`: This reads the raw tree data I got from Andy Finley, as well as the true unfuzzed coordinates (which you have to have locally), some raw data with species codes and lookup tables, the raw phylogeny, and raw trait data. It does the imputations of missing values using a multiple imputation procedure, and then saves the cleaned and processed trait and phylogeny data to a folder.

`nasabio/prep_diversity_files/fiabeta_prep_wholeusa.r`: This calculates the FD and PD distance matrices, cleans them up and gets rid of duplicates, and sums up the basal areas of trees in each plot to use as abundances for the site-by-species matrices for diversity calculations. It also adds a QC step of getting rid of plantation plots. All that stuff is saved as an R object so that it can be loaded in the diversity calculation step.

`nasabio/run_compile_diversity/fiawholeusa/<several scripts>.r`: These scripts will do the alpha, beta, and gamma diversity for all the plots. First, you calculate alpha diversity for each point, beta diversity for each pairwise combination within the highest radius we looked at, and gamma diversity for each radius. Then you average the diversities by radius to get alpha and beta (gamma is already done by radius). All of these scripts are set up to run in parallel, but if you are only doing a small area that will likely not be needed.

## Bird (BBS) diversity pipeline

`nasabio/trait_phylo_data_processing/bbs_consensustree.r`: This loads raw phylogenetic trees and gets a single consensus tree for use in PD calculations. (Even if you are only doing FD, the phylogeny is still needed because the traits are imputed with phylogenetic information).

`nasabio/trait_phylo_data_processing/bbs_impute.r`: Here we load trait data that has already been cleaned up using scripts from aquaxterra, as well as the consensus tree produced in the previous step. Use that tree to impute missing values. Save the imputed traits to a folder.

`nasabio/prep_diversity_files/bbsbeta_byroute_prep.r`: This script prepares FD and PD distance matrices for birds, as well as the route-by-species abundance matrices used for diversity calculation. (There is also a `bbs_bystop_prep.r` which makes additional abundance matrices for stop x species).

`nasabio/run_compile_diversity/bbs1year`, `bbseachyear`, and `bbswithinoute` folders: These folders contain scripts to calculate bbs alpha, beta, and gamma diversity three different ways: with all years pooled into one by route, with each year done separately by route, and with all years pooled but among stops within routes, respectively. The procedure is the same as for trees described above. The `bbs1year` folder is the most up to date, because it is the one used for the GEB paper.
