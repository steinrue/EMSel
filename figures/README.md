# Recreating main text figures

This directory contains scripts to recreate the data-derived figures in the main text and Supplementary Material of Fine and Steinrücken (2024), or generate similar plots for a dataset of your choosing. Step-by-step instructions to generate each plot can be found in the README of each subfolder.

To run the figures scripts, additional packages are needed - they can be installed by running the command
```
pip install "emsel[plots] @ git+https://github.com/steinrue/EMSel"
```

## simulation subdirectory

The subdirectory [simulation/](simulation/) contains scripts to generate the boxplots, strip plots, AUC plots, Q-Q plots, and confusion matrices found in Sections 3 and 4.2 of the manuscript, as well as Figures S.1-S.11 of the Supplementary Material.

Detailed descriptions on how to run the scripts can be found in the subdirectory.

## gb_dataset subdirectory

The subdirectory [gb_dataset/](gb_dataset/) contains scripts to generate all plots found in Section 4.3 of the manuscript, as well as Figures S.12-S.32 and Tables S.1-S.3 of the Supplementary Material. These plots _do_ require a VCF as input. 

Detailed descriptions on how to run the scripts can be found in the subdirectory.
