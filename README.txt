This is a github repository for my Master's final thesis, created on 12/12/2024.

For now, the goal of this repository is for me to get familiarized with the use of github, and to provide access for my TFM tutor to my scripts. 

At the moment, the repository contains 3 folders, one for each part of the analysis performed so far. The names of the folders are provisional.

1- GetCDS_ProtSeq_DisDom: Contains the workflow (Snakefile) to compute the disorder rate and functional domains of exons of the proteomes of differents species, as well as the associated config file (yaml) and an example LOG file that the completion of this workflow generates. Finished script.

2- InfoXspecies_toplot: Snakefile and config file used to extract and combine information from the results of the previous analysis in order to generate tables that facilitate the generation of plots (in R, manually). This is a provisional file, intended for personal use. 

3- Compare_TSrAS: Snakefile, config file and bin directory with an R script. Workflow under development. Used to define different exon sets (TS, PanAS, High_PSI, Low_PSI), obtain their disorder rate and associated domains, their corresponding genes and proteins... per species and to combine the results in summary files (information from all species). This workflow is currently in used but far from being completed. At the moment, it is intended for personal use, but will be part of the final memory of the TFM. 


The final version of the repository will include extra files, such as a specific README for each folder with details of the workflow and mat&met information, a yaml file for the conda environment beign used, json files specifying the resources needed for each workflow and additional scripts.

   
