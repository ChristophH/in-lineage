Data and code associated with the paper 

*Developmental diversification of cortical inhibitory interneurons*

by Christian Mayer#, Christoph Hafemeister#, Rachel C. Bandler#, Robert Machold, Renata Batista Brito, Xavier Jaglin, Kathryn Allaway, Andrew Butler, Gord Fishell\* and Rahul Satija\*

\# Equal contribution  
\* Corresponding authors

[*Nature*, Advanced Online Publication 05 March 2018, doi:10.1038/nature25999](https://dx.doi.org/10.1038/nature25999)

[Free read-only version of the full published article through Springer Nature SharedIt](http://rdcu.be/JA5l)

[Preprint on bioRxiv](https://www.biorxiv.org/content/early/2017/09/13/105312)  

### How to run

There are individual scripts for the different parts of the analysis. Run them in this order:

1. Rscript R/maturation_trajectory.R
2. Rscript R/mitotic_cells.R
3. Rscript R/branch_analysis.R CGE
4. Rscript R/branch_analysis.R LGE
5. Rscript R/branch_analysis.R MGE
6. Rscript R/map_10x_E14_to_branches.R
7. Rscript R/variance_explained.R

Code related to the integrated analysis linking heterogenity in our data to heterogeneity in adult cells (Figure 4) is provided by Christian Mayer [in his repository](https://github.com/mayer-lab/Mayer-et-al-2018_IntegratedAnalysis).

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
* all\_samples\_CGE\_branch\_analysis\_de\_genes.Rds
* all\_samples\_CGE\_branch\_analysis\_expr\_branch\_avg.Rds


running times on an Intel Xeon Processor E5-2697 v3 @ 2.6 to 3.6 GHz
(using 6 cores for some of the steps)

* CGE: ca. 10 minutes
* LGE: ca. 9 minutes
* MGE: ca. 6 minutes

### R/map_10x_E14_to_branches.R

this script performs the following steps:

1. load 10x data and subset to E14.5 Lhx6neg (CGE)  
2. re-normalize data and run maturation trajectory analysis to isolate post-mitotic cells  
3. cluster cells and remove Lhx6 positive clusters (contamination)  
4. also load the 10x E14.5 Lhx6pos (MGE) cells  
5. map both set of cells to the branches identified in E13.5 dropseq data  

these files are created in results directory:

* CGE\_E14.5\_Lhx6neg\_maturation\_trajectory.pdf
* CGE\_E14.5\_Lhx6neg\_postmitotic\_clusters\_lhx6\_detection\_rate.pdf
* CGE\_E14.5\_Lhx6neg\_mapped.Rds
* MGE\_E14.5\_Lhx6pos\_mapped.Rds
* CGE\_E14.5\_Lhx6neg\_mapped\_single\_cell\_heatmap.pdf
* MGE\_E14.5\_Lhx6pos\_mapped\_single\_cell\_heatmap.pdf

running times on an Intel Xeon Processor E5-2697 v3 @ 2.6 to 3.6 GHz
(using 6 cores for some of the steps)

*  ca. 33 minutes

### R/variance_explained.R

this script performs the following steps:

1. load mitotic cells from CGE and MGE dropseq experiments
2. quantify variance explained by individual factors
3. load E14.5 10x data and subset to postmitotic cells
4. quantify variance explained by individual factors

these files are created in results directory:

* variance\_explained\_dropseq\_mitotic\_CGE\_MGE.pdf
* variance\_explained\_10x\_E14.5\_CGE\_MGE.pdf

running times on an Intel Xeon Processor E5-2697 v3 @ 2.6 to 3.6 GHz
(using 6 cores for some of the steps)

* ca. 40 minutes
