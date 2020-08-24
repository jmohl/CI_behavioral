## Modeling of behavioral causal inference in an audio-visual localization task
This project is contains all code and data used for the manuscript "Monkeys and Humans Implement Causal Inference to Simultaneously Localize Auditory and Visual Stimuli" (Mohl, Pearson, & Groh 2019, bioRxiv doi:10.1101/823385, in press at Journal of Neurophysiology: https://doi.org/10.1152/jn.00046.2020)

### Description of contents
src: contains all source code relevant to this project, including code for preparing data, fitting models, and generating plots
- Running analysis: model fitting and analysis is initiated from master_script.m
- Plotting: all plots used in the manuscript are included in figure_script.m
- lautils: borrows utility functions for numerical integration from https://github.com/lacerbi/lautils-mat

data: contains data used for this paper, including all human and monkey behavioral files in as tables in .mat format

results: not included in repo, but all figures and intermediate data files will be generated in this folder when the code is run.

doc: not included in repo, published manuscript contains relevant documentation https://doi.org/10.1152/jn.00046.2020 (paywalled, if this is a problem for you please email me and I can send you a copy)

### Known issues
- Code was written for MATLAB version 2018b and will not work on version 2018a or earlier. This is likely do to the change of the 'sum()' function in that version
  - workaround: change all instances of sum(...,'all') to sum(sum(...)) or sum(sum(sum(...))) depending on the number of dimensions needed.
