# SV_Impute
project to benchmark imputation tools for the impotation of structural variants in human genomes

Relevant links: \
[Project Proposal to NCIG Board](<https://docs.google.com/document/d/1kotjFJI86qZz3PTGbJsK2hUQc9vtO0KzouX426ZW2uI/edit?usp=sharing>) \
[Project Outline that will become the manuscript](<https://docs.google.com/document/d/13Ug62GrY3RKd6AoOQfrVfWusoaQicKBtcgUPxR89NbQ/edit?usp=sharing>)

`lab_book.md` is essentially my "lab book" of working on the project.

Locations of files on NCI cluster:

Current working directory of the project within the te53 project directory: \
`/g/data/te53/sv_imputation`
In this directory you will find `dir_map.md` that explains all the files and directories.

`main.nf` nextflow script for the pipeline. 

`nextflow.config` Configuration file to run the pipeline on NCI cluster, defines resource limiations for each step
