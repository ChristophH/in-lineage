Data and code associated with the paper 

*Developmental diversification of cortical inhibitory interneurons*

by Christian Mayer#, Christoph Hafemeister#, Rachel C. Bandler#, Robert Machold, Renata Batista Brito, Xavier Jaglin, Kathryn Allaway, Andrew Butler, Gord Fishell\* and Rahul Satija\*

\# Equal contribution  
\* Corresponding authors

Preprint on [bioRxiv](https://doi.org/10.1101/105312)

### How to run

There are individual scripts for the different parts of the analysis. Run them in this order:

1. Rscript R/maturation_trajectory.R
2. Rscript R/mitotic_cells.R
3. Rscript R/branch_analysis.R CGE
4. Rscript R/branch_analysis.R LGE
5. Rscript R/branch_analysis.R MGE

### R/maturation_trajectory.R

this script performs the following steps:  

1. read in drop-seq digital expression and normalize using regularized NB regression  
2. cluster all cells and remove contaminating cell populations  
3. fit a maturation trajectory through the remaining cells  
4. identify maturation score cutoff to separate mitotic from post-mitotic cells  
5. visualize results  
6. create smooth expression (as function of maturation score) for visualization later on  

these files are created in the results directory:  

* all\_samples\_normalized\_expression.Rds  
* all\_samples\_maturation\_trajectory\_meta\_data.Rds  
* all\_samples\_maturation\_trajectory.pdf  
* all\_samples\_smooth\_expression.Rds  

running times on an Intel Xeon Processor E5-2697 v3 @ 2.6 to 3.6 GHz
(using 6 cores for some of the steps)  

* steps 1-5: ca. 75 minutes
* step 6: ca. 95 minutes

### R/mitotic_cells.R

this script performs the following steps:

1. run differential expression test of CGE vs MGE in early mitotic cells
2. identify genes that are associated with maturation in all eminences and all mitotic cells

note that you need to run maturation_trajectory.R first

these files are created in results directory:

* all\_samples\_differential\_expression\_early\_mitotic\_TFs\_CGE\_vs\_MGE.csv
* all\_samples\_temporal\_mitotic\_genes.pdf

running times on an Intel Xeon Processor E5-2697 v3 @ 2.6 to 3.6 GHz
(using 6 cores for some of the steps)

* ca. 6 minutes

### R/branch_analysis.R

this script performs the following steps:

1. load maturation trajectory results and keep only post-mitotic cells from
   one specific eminence (CGE by default, or first command line argument if present)
2. re-normalize data and perform dimensionality reduction
3. use bootstrapped minimum spanning trees to create new cell-to-cell distances
4. use consensus tree to identify branches
5. run differential expression tests between terminal branches to identify marker genes

these files are created in results directory (in the case of CGE):

* all\_samples\_CGE\_branch\_analysis.pdf
* all\_samples\_CGE\_branch\_analysis\_meta\_data.Rds
* all\_samples\_CGE\_branch\_analysis\_top\_de\_genes.pdf
* all\_samples\_CGE\_branch\_analysis\_top\_marker\_genes.csv

running times on an Intel Xeon Processor E5-2697 v3 @ 2.6 to 3.6 GHz
(using 6 cores for some of the steps)

* CGE: ca. 10 minutes
* LGE: ca. 9 minutes
* MGE: ca. 6 minutes
