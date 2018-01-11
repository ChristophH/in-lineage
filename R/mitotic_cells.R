# this script performs the following steps:
# 1) run differential expression test of CGE vs MGE in early mitotic cells
# 2) identify genes that are associated with maturation in all eminences and all mitotic cells

# note that you need to run maturation_trajectory.R first

# these files are created in results directory:
# - all_samples_differential_expression_early_mitotic_TFs_CGE_vs_MGE.csv
# - all_samples_temporal_mitotic_genes.pdf

# running times on an Intel Xeon Processor E5-2697 v3 @ 2.6 to 3.6 GHz
# (using 6 processors for some of the steps)
# ca. 6 minutes


# support functions are defined in library
source('R/lib.R')

set.seed(42)
options(mc.cores = 6)

# set basename for loading results of the maturation trajectory analysis
# results of this script will use the same prefix
result.bn <- 'all_samples'


# load maturation trajectory results and keep only earliest cells, e.g. those with maturation score < 0.3
md <- readRDS(file = sprintf('results/%s_maturation_trajectory_meta_data.Rds', result.bn))
cm <- readRDS('data/dropseq_digitial_expression.Rds')

md <- subset(md, maturation.score.smooth < 0.3)
cm <- cm[, rownames(md)]

# run differential expression test of one eminence vs another
de.out <- de.nb.reg.helper(cm, md, c('reads', 'mols.per.gene', 'cc'), 'eminence', 'CGE', 'MGE')

mouse.tfs <- sort(unique(read.table('annotation/Mus_musculus_TFs_mart_export.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')$Associated.Gene.Name))
de.out$is.TF <- de.out$gene %in% mouse.tfs

# filter to only show significantly different TFs
sig.p.th <- 1e-4
sig.fc.th <- 1
df <- subset(de.out, (freq1 > 0.1 | freq2 > 0.1) & abs(log.fc) > sig.fc.th & is.TF & adj.pval < sig.p.th)
write.csv(df[, c('gene', 'log.fc', 'freq1', 'freq2', 'adj.pval')], row.names = FALSE,
          file = sprintf('results/%s_differential_expression_early_mitotic_TFs_CGE_vs_MGE.csv', result.bn))


# find genes that have high mutual information with maturation score
# keep only mitotic cells and protein coding genes
pcg <- read.table('annotation/Mus_musculus.GRCm38.84.protein_coding_genes.txt', stringsAsFactors=FALSE)$V1
md <- readRDS(file = sprintf('results/%s_maturation_trajectory_meta_data.Rds', result.bn))
expr <- readRDS(file = sprintf('results/%s_normalized_expression.Rds', result.bn))
md <- subset(md, !postmitotic)
expr <- expr[rownames(expr) %in% pcg, rownames(md)]

# get mutual information between gene expression and maturation score
mi.out <- mi(t(expr), md$maturation.score.smooth, nbins=13)
# shuffle maturation score to get distribution of random background mutual information
set.seed(42)
mi.out.r <- Reduce(cbind, lapply(1:13, function(i) mi(t(expr), sample(md$maturation.score.smooth), nbins=13)))
# z-score mutual information with respect to random background
mi.z <- (mi.out - apply(mi.out.r, 1, mean)) / apply(mi.out.r, 1, sd)

# load the list with smooth expression (expression as a function of maturation score)
fit.lst <- readRDS(file = sprintf('results/%s_smooth_expression.Rds', result.bn))
x.pred <- as.numeric(colnames(fit.lst[[1]]))
# use only those maturation score bins with at least 13 cells in each eminence
mss <- md$maturation.score.smooth
use.pred <- sapply(x.pred, function(x) sum(mss < x & md$eminence == 'CGE')) >= 13 & 
  sapply(x.pred, function(x) sum(mss < x & md$eminence == 'LGE')) >= 13 &
  sapply(x.pred, function(x) sum(mss < x & md$eminence == 'MGE')) >= 13 & x.pred < max(mss)

# select the genes with high mutual information with maturation score
# and high conserved pattern across eminences
dyn.genes <- rownames(mi.z)[mi.z > 20]
cvec.cge <- diag(cor(t(fit.lst[['all']][dyn.genes, use.pred]), t(fit.lst[['CGE']][dyn.genes, use.pred])))
cvec.lge <- diag(cor(t(fit.lst[['all']][dyn.genes, use.pred]), t(fit.lst[['LGE']][dyn.genes, use.pred])))
cvec.mge <- diag(cor(t(fit.lst[['all']][dyn.genes, use.pred]), t(fit.lst[['MGE']][dyn.genes, use.pred])))
plot.genes <- dyn.genes[cvec.cge > 0.90 & cvec.lge > 0.90 & cvec.mge > 0.90]

# order the genes by when they reach 80% of their max first
max.pt <- apply(fit.lst[['all']][plot.genes, use.pred], 1, function(x) which(x >= quantile(x, probs = 0.8))[1])
# show heatmap of selected genes
mat <- t(scale(t(fit.lst[['all']][plot.genes[order(max.pt)], use.pred])))
mat[mat > 2] <- 2
mat[mat < -2] <- -2
rownames(mat) <- rep('', nrow(mat))
colnames(mat) <- round(as.numeric(colnames(mat)), 2)
pdf(sprintf('results/%s_temporal_mitotic_genes.pdf', result.bn), width = 7, height = 10)
hm.out <- heatmap.2(mat, trace='none', Colv=NA, Rowv=NA, dendrogram = 'none',
                    scale='row', symbreaks=TRUE, cexCol = 0.5,
                    col=colorRampPalette(c("#fff7fb", "#fff7fb", "#d0d1e6", "#74a9cf", "#0570b0", "#023858", "#023858"))(31))
dev.off()

