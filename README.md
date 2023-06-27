# PCR data accurately predict infectious virus: a characterization of SARS-CoV-2 in non-human primates

Celine Snedden(1), James Lloyd-Smith(1, 2)

1. Department of Ecology and Evolutionary Biology, University of California, Los Angeles, Los Angeles, CA, USA
2. Department of Computational Medicine, University of California, Los Angeles, Los Angeles, CA, USA

The preprint is available at: https://doi.org/10.1101/2023.06.23.546114


## Repository information

This repository provides the data and code to reproduce the analyses presented in the associated paper. Additional analyses mentioned but not shown in display elements of that paper are also included. 

Instructions to recreate all analyses are summarized below, including explanations of adapting the hurdle model framework to your own data. 

## Article abstract

Researchers and clinicians often rely on molecular assays like PCR to identify and monitor viral infections instead of the resource-prohibitive gold standard of viral culture. However, it remains unclear when (if ever) PCR measurements of viral load are reliable indicators of replicating or infectious virus. Here, we compare total RNA, subgenomic RNA, and viral culture results from 24 studies of SARS-CoV-2 in non-human primates using bespoke statistical models. On out-of-sample data, our best models predict subgenomic RNA from total RNA with 91% accuracy, and they predict culture positivity with 85% accuracy. Total RNA and subgenomic RNA showed equivalent performance as predictors of culture positivity. Multiple cofactors, including exposure conditions and host traits, influence culture predictions for total RNA quantities spanning twelve orders of magnitude. Our model framework can be adapted to compare any assays, in any host species, and for any virus, to support laboratory analyses, medical decisions, and public health guidelines. 

## Citation information

If you use any of the code or data associated with this manuscript, please cite our work. (_Citation information to come_). 

## File structure and naming conventions

All content is organized into the following folders:

- `code`: contains all code, including processing data, model definitions, model selection, final model fitting, and figure and table generation. 
  - primary analyses are ordered by increasing numbers (e.g., `01-prep-clean-data.R`, `02-PCR-model-selection.R`) representing the sequence in which they were run for the investigations presented in the paper.
  - unless otherwise specified, each file should run without needing to run any other files prior
  - files that generate certain figures or tables are preceded with `fig` or `tbl` and the corresponding number from the paper (if presented).

- `data`: contains data as `.csv` files, including raw data (`raw_data.csv`) and cleaned data (`clean-data.csv`). 
  - additional data files (`pred-sg-data.csv`, `pred-culture-data.csv`) are based on the cleaned data but contain additional columns with predictions generated by the final (best) sgRNA and culture models for each relevant sample.

- `outputs`: contains all output files, including figures and tables for the main text, supplementary text, additional visualizations, and model fits.
  - just as for code files, each output is marked as a figure (`fig`) or a table (`table`) preceding the corresponding number. The code generating the figure can be found with the corresponding name in the code file. 
  - `outputs/fits`: contains model fits for the best and simplest models, including those generated for informative and non-informative priors and any other relevant analyses.

Any file marked with a preceding `EA-` indicates **E**xtra **A**nalyses, which are not presented formally as display elements in the article but are mentioned in passing. This convention holds for extra code, figures, and table files. 

## Further database information

`tbl-database-summary.csv` includes additional information relevant for each article, including whether data was digitized, etc.

## Installing software dependencies

To run our analyses, you will need to install certain software, described below. 

- We use the statistical programming language R, for which you can find installation instructions at: https://cran.r-project.org/doc/manuals/r-release/R-admin.html. 

- We use both CmdStanR and RStan for Bayesian model compilation. We recommend that you: 

  1. review and follow the installation instructions at https://mc-stan.org/cmdstanr/articles/cmdstanr.html and https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
  
  2. run the example models to confirm performance is as expected before proceeding. 

- We also use various other R packages. Each code file contains prompts to download and attach any required external packages. 
