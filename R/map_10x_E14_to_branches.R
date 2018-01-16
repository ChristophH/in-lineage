# this script performs the following steps:
# 1) load 10x data and subset to E14.5 Lhx6neg (CGE)
# 2) re-normalize data and run maturation trajectory analysis to isolate post-mitotic cells
# 3) cluster cells and remove Lhx6 positive clusters (contamination)
# 4) also load the 10x E14.5 Lhx6pos (MGE) cells
# 5) map both set of cells to the branches identified in E13.5 dropseq data

# these files are created in results directory:
# - CGE_E14.5_Lhx6neg_maturation_trajectory.pdf
# - CGE_E14.5_Lhx6neg_postmitotic_clusters_lhx6_detection_rate.pdf
# - CGE_E14.5_Lhx6neg_mapped.Rds
# - MGE_E14.5_Lhx6pos_mapped.Rds
# - CGE_E14.5_Lhx6neg_mapped_single_cell_heatmap.pdf
# - MGE_E14.5_Lhx6pos_mapped_single_cell_heatmap.pdf

# running times on an Intel Xeon Processor E5-2697 v3 @ 2.6 to 3.6 GHz
# (using 6 cores for some of the steps)
# ca. 33 minutes

# support functions are defined in library
source('R/lib.R')

set.seed(42)
options(mc.cores = 6)

# set basename for loading and saving results
#input.bn <- 'all_samples'
#result.bn <- 'all_samples'

# load digital expression matrix and meta data
cm <- readRDS('data/10x_digitial_expression.Rds')
md <- readRDS('data/10x_meta_data.Rds')

# since the CGE_E14.5_Lhx6neg sample contains a mix of mitotic and post-mitotic cells
# we need to do a maturation trajectory analysis (as for the dropseq data)
# in order to isolate the post-mitotic cells
md.cge <- subset(md, sample.name == 'CGE_E14.5_Lhx6neg')
cm.cge <- cm[, rownames(md.cge)]

# get cell-cycle score
cc <- get.cc.score(cm.cge, seed=42)
md.cge$cc <- cc$score
md.cge$cc.phase <- cc$phase

# normalize data 
genes <- rownames(cm.cge)[apply(cm.cge > 0, 1, mean) >= 0.005 & apply(cm.cge > 0, 1, sum) >= 3]
cat('Normalizing', length(genes), 'genes that are present in at least 0.5% of the cells AND in at least 3 cells\n')

expr <- norm.nb.reg(cm.cge[genes, ], md.cge[, c('reads', 'mols.per.gene', 'cc')], pr.th = 30)
# save the normalized expression data (this could be a rather large file)
#saveRDS(expr, file = sprintf('results/%s_normalized_expression.Rds', result.bn))

# keep only protein coding genes
pcg <- read.table('annotation/Mus_musculus.GRCm38.84.protein_coding_genes.txt', stringsAsFactors=FALSE)$V1
expr <- expr[rownames(expr) %in% pcg, ]
cm.cge <- cm.cge[rownames(cm.cge) %in% pcg, ]

# cluster cells and remove contaminating populations
cat('Doing initial clustering\n')
cl <- cluster.the.data.simple(cm.cge, expr, 9, seed=42)
md.cge$init.cluster <- cl$clustering
# detection rate per cluster for some marker genes
goi <- c('Igfbp7', 'Col4a1', 'Neurod2', 'Neurod6')
det.rates <- apply(cm.cge[goi, ] > 0, 1, function(x) aggregate(x, by=list(cl$clustering), FUN=mean)$x)
contam.clusters <- sort(unique(cl$clustering))[apply(det.rates > 1/3, 1, any)]
use.cells <- !(cl$clustering %in% contam.clusters)
cat('Of the', ncol(cm.cge), 'cells', sum(!use.cells), 'are determined to be part of a contaminating cell population.\n')
cm.cge <- cm.cge[, use.cells]
expr <- expr[, use.cells]
md.cge <- md.cge[use.cells, ]

# fit maturation trajectory
mat.traj <- maturation.trajectory(cm.cge, md.cge, expr, pricu.f=1/3)
md.cge <- mat.traj$md
pc.line <- mat.traj$pc.line
mt.th <- mat.traj$mt.th

# visualize result
pdf(sprintf('results/%s_maturation_trajectory.pdf', 'CGE_E14.5_Lhx6neg'), width = 7, height = 5)

g <- ggplot(md.cge, aes(DMC1, DMC2)) + geom_point(aes(color=maturation.score.smooth), size=1, shape=16) + 
  scale_color_gradientn(colours=my.cols.RYG, name='Maturation score') +
  stat_density2d(n=111, na.rm=TRUE, color='black', size=0.33, alpha=0.5) +
  geom_line(data=pc.line, color='deeppink', size=0.77) +
  theme_grey(base_size=12) + labs(x='DMC1', y='DMC2')
plot(g)

g <- ggplot(md.cge, aes(DMC1, DMC2)) + geom_point(aes(color=rank(maturation.score.smooth)), size=1, shape=16) + 
  scale_color_gradientn(colours=my.cols.RYG, name='Maturation score rank') +
  stat_density2d(n=111, na.rm=TRUE, color='black', size=0.33, alpha=0.5) +
  geom_line(data=pc.line, color='deeppink', size=0.77) +
  theme_grey(base_size=12) + labs(x='DMC1', y='DMC2')
plot(g)

g <- ggplot(md.cge, aes(maturation.score.smooth, cc.phase.fit)) + geom_point(aes(color=postmitotic), size=2) +
  geom_hline(yintercept =  mean(md.cge$in.cc.phase)/2) + geom_vline(xintercept = mt.th) +
  theme_grey(base_size=12) + labs(x='Maturation score', y='Fraction of cells in G2/M or S phase')
plot(g)

g <- ggplot(md.cge, aes(DMC1, DMC2)) + geom_point(aes(color=postmitotic), size=1, shape=16) +
  theme_grey(base_size=12) + labs(x='DMC1', y='DMC2')
plot(g)
dev.off()

# keep only the post-mitotic cells
md.cge <- subset(md.cge, postmitotic)
cm.cge <- cm.cge[, rownames(md.cge)]
expr <- expr[, rownames(md.cge)]

# cluster and get Lhx6 detection rate
cl <- cluster.the.data.simple(cm.cge, expr, 9, seed=42)
det.rate <- data.frame(aggregate(cm.cge['Lhx6', ] > 0, by=list(cl$clustering), FUN=mean))
colnames(det.rate) <- c('cluster', 'Lhx6.detection.rate')
det.rate$cluster <- factor(det.rate$cluster, levels=det.rate$cluster, ordered=TRUE)
g <- ggplot(det.rate, aes(cluster, Lhx6.detection.rate)) + geom_bar(stat='identity') 
ggsave('results/CGE_E14.5_Lhx6neg_postmitotic_clusters_lhx6_detection_rate.pdf', width=8, height=6, units='cm')

# discard Lhx6+ cells and those that are in Lhx6+ clusters
keep <- cm.cge['Lhx6', ] <= 0 & cl$clustering %in% as.numeric(as.character(subset(det.rate, Lhx6.detection.rate <= 0.2)$cluster))
cm.cge <- cm.cge[, keep]
md.cge <- md.cge[keep, ]
expr.cge <- expr[, keep]



# get the Lhx6 positive (MGE) data set
md.mge <- subset(md, sample.name == 'MGE_E14.5_Lhx6pos')
cm.mge <- cm[, rownames(md.mge)]

# get cell-cycle score
cc <- get.cc.score(cm.mge, seed=42)
md.mge$cc <- cc$score
md.mge$cc.phase <- cc$phase

# normalize data 
genes <- rownames(cm.mge)[apply(cm.mge > 0, 1, mean) >= 0.005 & apply(cm.mge > 0, 1, sum) >= 3]
cat('Normalizing', length(genes), 'genes that are present in at least 0.5% of the cells AND in at least 3 cells\n')
expr.mge <- norm.nb.reg(cm.mge[genes, ], md.mge[, c('reads', 'mols.per.gene', 'cc')], pr.th = 30)


# map the cells
de.genes.cge <- readRDS(file = 'results/all_samples_CGE_branch_analysis_de_genes.Rds')
de.genes.mge <- readRDS(file = 'results/all_samples_MGE_branch_analysis_de_genes.Rds')
de.genes.both <- union(de.genes.cge, de.genes.mge)

# get branch average for CGE dropseq data (only genes DE between branches)
expr.avg <- readRDS(file = 'results/all_samples_CGE_branch_analysis_expr_branch_avg.Rds')
de.genes <- Reduce(intersect, list(de.genes.both, rownames(expr.avg), rownames(expr.cge)))
# map the CGE cells
map.cge <- transfer.label.cor(expr.avg[de.genes, ], colnames(expr.avg), expr.cge[de.genes, ])
table(map.cge$label, map.cge$cor.p.adjust < 0.1)

# get branch average for MGE dropseq data (only genes DE between branches)
expr.avg <- readRDS(file = 'results/all_samples_MGE_branch_analysis_expr_branch_avg.Rds')
de.genes <- Reduce(intersect, list(de.genes.both, rownames(expr.avg), rownames(expr.mge)))
# map the MGE cells
map.mge <- transfer.label.cor(expr.avg[de.genes, ], colnames(expr.avg), expr.mge[de.genes, ])
table(map.mge$label, map.mge$cor.p.adjust < 0.1)

# save results
saveRDS(map.cge, 'results/CGE_E14.5_Lhx6neg_mapped.Rds')
saveRDS(map.mge, 'results/MGE_E14.5_Lhx6pos_mapped.Rds')

# plot results
map.cge$label[map.cge$cor.p.adjust >= 0.1] <- NA
hm <- sc.cluster.heatmap(expr.cge, map.cge$label, 15, 'results/CGE_E14.5_Lhx6neg_mapped_single_cell_heatmap.pdf')
map.mge$label[map.mge$cor.p.adjust >= 0.1] <- NA
hm <- sc.cluster.heatmap(expr.mge, map.mge$label, 15, 'results/MGE_E14.5_Lhx6pos_mapped_single_cell_heatmap.pdf')



