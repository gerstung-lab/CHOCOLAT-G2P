| Script Name                   | Input                                            | Output                           | Description                                                              |
|-------------------------------|--------------------------------------------------|----------------------------------|--------------------------------------------------------------------------|
| 00_compose_data.py            | - Raw Visium files (#01)                         | - 11 AnnData objects (#03)       | - Composes raw Visium files and nodule histopathological annotations    |
|                               | - Nodule histopathological annotations (#02)     |                                  |                                                                          |
| 01_nodule_aggregates.py      | - AnnData objects (#03)                          | - 1 AnnData object (#04)         | - Computes nodule-level normalized expression value aggregates          |
|                               |                                                  |                                  |   - Phenotype and TME binarization                                        |
|-------------------------------|--------------------------------------------------|----------------------------------|--------------------------------------------------------------------------|
| 02_gene_coexp_proproc.py     | - AnnData objects (#03)                          | - 1 AnnData object (#05)         | - Concatenates normalized and gene-subset expression matrices            |
|-------------------------------|--------------------------------------------------|----------------------------------|--------------------------------------------------------------------------|
| 02_gene_coexp_viz.Rmd        | - AnnData object (#05)                           | - Plot (EDF 13)                   | - Generates gene-coexpression network visualization using R Markdown    |
|-------------------------------|--------------------------------------------------|----------------------------------|--------------------------------------------------------------------------|
| 03_heatmaps_phenotypes.Rmd   | - AnnData object (#03)                           | - Plot (Fig. 3, EDF 11, EDF 16)  | - Produces global heatmaps for ROIs and slides using R Markdown          |
|                               |                                                  |                                  |   - Separated heatmaps per slide                                         |
|-------------------------------|--------------------------------------------------|----------------------------------|--------------------------------------------------------------------------|
| 04_nodule_phenotyping_factors.py | - AnnData objects (#03)                       | - 1 AnnData object (#06)         | - Computes factor and gene loadings per spot for all slides             |
|                                  |                                                  |                                  |   - Generates factor to cell type signature map                          |
|-------------------------------|--------------------------------------------------|----------------------------------|--------------------------------------------------------------------------|
| 05_phenotype_spatial_maps.py | - AnnData object containing expression aggregates and binarization (#04) | - Plot (Fig. 3, EDF 14) | - Generates spatial maps of aggregates and single gene expressions      |
|                               | - AnnData object containing factors (#06)       |                                  |                                                                          |
|-------------------------------|--------------------------------------------------|----------------------------------|--------------------------------------------------------------------------|
| utils_functions.py           |                                                  |                                  | - Auxiliary functions used for R scripts                                 |
|-------------------------------|--------------------------------------------------|----------------------------------|--------------------------------------------------------------------------|
| utils_dependencies.py        |                                                  |                                  | - Packages used for R scripts                                             |
