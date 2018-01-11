# this script performs the following steps:
# 1) read in drop-seq digital expression and normalize using regularized NB regression
# 2) cluster all cells and remove contaminating cell populations
# 3) fit a maturation trajectory through the remaining cells
# 4) identify maturation score cutoff to separate mitotic from post-mitotic cells
# 5) visualize results
# 6) create smooth expression (as function of maturation score) for visualization later on

# these files are created in results directory:
# - all_samples_normalized_expression.Rds
# - all_samples_maturation_trajectory_meta_data.Rds
# - all_samples_maturation_trajectory.pdf
# - all_samples_smooth_expression.Rds

# running times on an Intel Xeon Processor E5-2697 v3 @ 2.6 to 3.6 GHz
# (using 6 processors for some of the steps)
# steps 1-5: ca. 75 minutes
# step 6: ca. 95 minutes

# support functions are defined in library
source('R/lib.R')

set.seed(42)
options(mc.cores = 6)


# load digital expression matrix and meta data
cm <- readRDS('data/dropseq_digitial_expression.Rds')
md <- readRDS('data/dropseq_meta_data.Rds')

# set basename for result files
result.bn <- 'all_samples'

# get cell-cycle score
cc <- get.cc.score(cm, seed=42)
md$cc <- cc$score
md$cc.phase <- cc$phase

# normalize data 
genes <- rownames(cm)[apply(cm > 0, 1, mean) >= 0.005 & apply(cm > 0, 1, sum) >= 3]
cat('Normalizing', length(genes), 'genes that are present in at least 0.5% of the cells AND in at least 3 cells\n')

md$mols.per.gene <- md$mols / md$genes
expr <- norm.nb.reg(cm[genes, ], md[, c('reads', 'mols.per.gene', 'cc')], pr.th = 30)
# save the normalized expression data (this could be a rather large file)
saveRDS(expr, file = sprintf('results/%s_normalized_expression.Rds', result.bn))

# keep only protein coding genes
pcg <- read.table('annotation/Mus_musculus.GRCm38.84.protein_coding_genes.txt', stringsAsFactors=FALSE)$V1
expr <- expr[rownames(expr) %in% pcg, ]
cm <- cm[rownames(cm) %in% pcg, ]

# cluster cells and remove contaminating populations
cat('Doing initial clustering\n')
cl <- cluster.the.data.simple(cm, expr, 9, seed=42)
md$init.cluster <- cl$clustering
# detection rate per cluster for some marker genes
goi <- c('Igfbp7', 'Col4a1', 'Neurod2', 'Neurod6')
det.rates <- apply(cm[goi, ] > 0, 1, function(x) aggregate(x, by=list(cl$clustering), FUN=mean)$x)
contam.clusters <- sort(unique(cl$clustering))[apply(det.rates > 1/3, 1, any)]
use.cells <- !(cl$clustering %in% contam.clusters)
cat('Of the', ncol(cm), 'cells', sum(!use.cells), 'are determined to be part of a contaminating cell population.\n')
cm <- cm[, use.cells]
expr <- expr[, use.cells]
md <- md[use.cells, ]


# fit maturation trajectory
cat('Fitting maturation trajectory\n')
genes <- apply(cm[rownames(expr), ] > 0, 1, mean) >= 0.02 & apply(cm[rownames(expr), ] > 0, 1, sum) >= 3
cat('Using', length(genes), 'genes in diffusion map\n')
rd <- dim.red(expr[genes, ], max.dim=50, ev.red.th=0.04)
# for a consisten look use Nes expression to orient each axis
for (i in 1:ncol(rd)) {
  if (cor(expr['Nes', ], rd[, i]) > 0) {
    rd[, i] <- -rd[, i]
  }
}

md <- cbind(md, rd)
# save the updated meta data
saveRDS(md, file = sprintf('results/%s_maturation_trajectory_meta_data.Rds', result.bn))

# fit principal curve through the data
pricu <- principal.curve(rd, smoother='lowess', trace=TRUE, f=1/3, stretch=333)
pc.line <- as.data.frame(pricu$s[order(pricu$lambda), ])
md$maturation.score <- pricu$lambda/max(pricu$lambda)

# orient maturation score using Nes expression
if (cor(md$maturation.score, expr['Nes', ]) > 0) {
  md$maturation.score <- -(md$maturation.score - max(md$maturation.score))
}

# use 1% of neighbor cells to smooth maturation score
md$maturation.score.smooth <- nn.smooth(md$maturation.score, rd[, 1:2], round(ncol(expr)*0.01, 0))


# pick maturation score cutoff to separate mitotic from post-mitotic cells
in.phase <- md$cc.phase != 0
fit <- loess(as.numeric(in.phase) ~ md$maturation.score.smooth, span=0.5, degree=2)
md$cc.phase.fit <- fit$fitted
# pick MT threshold based on drop in cc.phase cells
# ignore edges of MT because of potential outliers
mt.th <- max(subset(md, cc.phase.fit > mean(in.phase)/2 & maturation.score.smooth >= 0.2 & maturation.score.smooth <= 0.8)$maturation.score.smooth)

md$postmitotic <- md$maturation.score.smooth > mt.th


# save the meta data including the maturation trajectory results
saveRDS(md, file = sprintf('results/%s_maturation_trajectory_meta_data.Rds', result.bn))


# visualize result
pdf(sprintf('results/%s_maturation_trajectory.pdf', result.bn), width = 7, height = 5)

g <- ggplot(md, aes(DMC1, DMC2)) + geom_point(aes(color=maturation.score.smooth), size=1, shape=16) + 
  scale_color_gradientn(colours=my.cols.RYG, name='Maturation score') +
  stat_density2d(n=111, na.rm=TRUE, color='black', size=0.33, alpha=0.5) +
  geom_line(data=pc.line, color='deeppink', size=0.77) +
  theme_grey(base_size=12) + labs(x='DMC1', y='DMC2')
plot(g)

g <- ggplot(md, aes(DMC1, DMC2)) + geom_point(aes(color=rank(maturation.score.smooth)), size=1, shape=16) + 
  scale_color_gradientn(colours=my.cols.RYG, name='Maturation score rank') +
  stat_density2d(n=111, na.rm=TRUE, color='black', size=0.33, alpha=0.5) +
  geom_line(data=pc.line, color='deeppink', size=0.77) +
  theme_grey(base_size=12) + labs(x='DMC1', y='DMC2')
plot(g)

g <- ggplot(md, aes(maturation.score.smooth, cc.phase.fit)) + geom_point(aes(color=postmitotic), size=2) +
  geom_hline(yintercept =  mean(in.phase)/2) + geom_vline(xintercept = mt.th) +
  theme_grey(base_size=12) + labs(x='Maturation score', y='Fraction of cells in G2/M or S phase')
plot(g)

g <- ggplot(md, aes(DMC1, DMC2)) + geom_point(aes(color=postmitotic), size=1, shape=16) +
  theme_grey(base_size=12) + labs(x='DMC1', y='DMC2')
plot(g)
dev.off()


# create smooth expression (as function of maturation score) for visualization later on
x.pred <- seq(min(md$maturation.score.smooth), max(md$maturation.score.smooth), length.out=100)
c.cge <- which(md$eminence == 'CGE')
c.lge <- which(md$eminence == 'LGE')
c.mge <- which(md$eminence == 'MGE')
expr.cge.fit <- smooth.expr(expr[, c.cge], md$maturation.score.smooth[c.cge], x.pred)
expr.lge.fit <- smooth.expr(expr[, c.lge], md$maturation.score.smooth[c.lge], x.pred)
expr.mge.fit <- smooth.expr(expr[, c.mge], md$maturation.score.smooth[c.mge], x.pred)
expr.all.fit <- smooth.expr(expr, md$maturation.score.smooth, x.pred)
fit.lst <- list(all=expr.all.fit, CGE=expr.cge.fit, LGE=expr.lge.fit, MGE=expr.mge.fit)
# save the smooth expression data
saveRDS(fit.lst, file = sprintf('results/%s_smooth_expression.Rds', result.bn))

