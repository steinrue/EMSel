# Recreating manuscript figures

This directory contains scripts to recreate the data-derived figures in the main text and Supplementary Material of Fine and Steinrücken (2024), or generate similar plots for a dataset of your choosing. Step-by-step instructions to generate each plot can be found in the README of each subfolder.

To run the figures scripts, additional packages are needed - they can be installed by running the command
```
pip install "emsel[plots] @ git+https://github.com/steinrue/EMSel"
```

Additionally, note that the code to recreate the figures is set up to use a SLURM-based cluster to parallelize computation. If this is not available, commands to run each analysis manually are provided, though note that this will be significantly slower.

## simulation subdirectory

The subdirectory [simulation/](simulation/) contains scripts to generate the boxplots, strip plots, AUC plots, Q-Q plots, and confusion matrices found in Section 3.1 and 3.2.2 of the manuscript, as well as Figures S.1-S.11 and S.22A of the Supplementary Material.

Detailed descriptions on how to run the scripts can be found in the subdirectory.

## gb_dataset subdirectory

The subdirectory [gb_dataset/](gb_dataset/) contains scripts to generate all plots found in Section 3.2.3-3.2.5 of the manuscript, as well as Figures S.10-S.38 and Tables S.1-S.3 of the Supplementary Material. These plots _do_ require a VCF as input. 

Detailed descriptions on how to run the scripts can be found in the subdirectory.

## horse_dataset subdirectory

The subdirectory [horse_dataset/](horse_dataset/) contains scripts to generate Figure 15 and Table 2 of the manuscript and carry out all related analyses.
