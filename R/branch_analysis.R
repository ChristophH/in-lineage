# this script performs the following steps:
# 1) load maturation trajectory results and keep only post-mitotic cells from
#    one specific eminence (CGE by default, or first command line argument if present)
# 2) re-normalize data and perform dimensionality reduction
# 3) use bootstrapped minimum spanning trees to create new cell-to-cell distances
# 4) use consensus tree to identify branches
# 5) run differential expression tests between terminal branches to identify marker genes

# these files are created in results directory (in the case of CGE):
# - all_samples_CGE_branch_analysis.pdf
# - all_samples_CGE_branch_analysis_meta_data.Rds
# - all_samples_CGE_branch_analysis_top_de_genes.pdf
# - all_samples_CGE_branch_analysis_top_marker_genes.csv
# - all_samples_CGE_branch_analysis_de_genes.Rds
# - all_samples_CGE_branch_analysis_expr_branch_avg.Rds

# running times on an Intel Xeon Processor E5-2697 v3 @ 2.6 to 3.6 GHz
# (using 6 cores for some of the steps)
# CGE: ca. 10 minutes
# LGE: ca. 9 minutes
# MGE: ca. 6 minutes


# support functions are defined in library
source('R/lib.R')

set.seed(42)
options(mc.cores = 6)

# set basename for loading and saving results
input.bn <- 'all_samples'
# set eminence to do branch analysis on
emin <- 'CGE'
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  emin <- args[1]
}
if (!(emin %in% c('CGE', 'LGE', 'MGE'))) {
  stop('First command line argument must be either CGE, LGE, or MGE')
}
result.bn <- paste(input.bn, emin, sep = '_')

# load maturation trajectory results and keep only post-mitotic cells from one eminence
md <- readRDS(file = sprintf('results/%s_maturation_trajectory_meta_data.Rds', input.bn))
cm <- readRDS('data/dropseq_digitial_expression.Rds')

md <- subset(md, postmitotic & eminence == emin)
cm <- cm[, rownames(md)]
cat('There are', ncol(cm), 'post-mitotic cells\n')

# keep only protein coding genes
pcg <- read.table('annotation/Mus_musculus.GRCm38.84.protein_coding_genes.txt', stringsAsFactors=FALSE)$V1
cm <- cm[rownames(cm) %in% pcg, ]

# re-normalize data 
# since we are now only looking at the post-mitotic subset of the data
genes <- rownames(cm)[apply(cm>0, 1, sum) >= 5]
cat('Normalizing', length(genes), 'genes that have been detected in at least 5 post-mitotic cells\n')

if (length(unique(md$sample.name)) > 1) {
  expr <- norm.nb.reg(cm[genes, ], md[, c('reads', 'mols.per.gene', 'cc', 'sample.name')], min.theta=0.01, pr.th=30, bins=64)
} else {
  expr <- norm.nb.reg(cm[genes, ], md[, c('reads', 'mols.per.gene', 'cc')], min.theta=0.01, pr.th=30, bins=64)
}

# use only variable genes 
vg <- genes[scale(sqrt(apply(expr^2, 1, sum)))[, 1] > 1]

# run diffusion map and return significant dimensions
dm <- dim.red(expr[vg, ], max.dim=50, ev.red.th=0.04, do.scale.result = TRUE)

# generate set of bootstrapped minimum spanning trees
# and average distance and adjacency matrices
set.seed(42)
mst.avg.res <- mst.avg(as.matrix(dist(dm)), 50, floor(0.66*ncol(expr)))

d <- mst.avg.res$dmats.avg  # new cell-to-cell distances

# pick as root the cell that maximizes correlation of maturation score and distance relative to all other cells
root <- which.max(cor(d, md$maturation.score.smooth, method='spearman')[, 1])
# get each cell's distance from root 
dist.from.root <- d[root, ]

# calculate 2D MDS layout of cell-to-cell distances
mds <- igraph::layout.mds(igraph::graph.full(nrow(d)), dist=d)

#consensus.mst <- igraph::mst(igraph::graph_from_adjacency_matrix(d, mode='undirected', weighted=TRUE, add.colnames=NA))
consensus.mst <- igraph::mst(igraph::graph_from_adjacency_matrix(as.matrix(dist(mds)), mode='undirected', weighted=TRUE, add.colnames=NA))

# the consensus MST is undirected - use root to convert it to a directed acyclic graph
dag <- tree.to.dag(consensus.mst, root)
branch <- branch.id(dag, w.th=0.08)

# identify terminal branches (the rest is considered part of the root)
dag.amat <- igraph::get.adjacency(dag, type='both', attr=NULL, sparse=FALSE)
branch.label <- rep(0, length(branch))
child.br <- apply(dag.amat, 1, function(x) branch[which(x==1)])
term.br <- c()
for (b in sort(unique(branch))) {
  if (length(unique(unlist(child.br[branch==b]))) == 1) {
    #term.br <- c(term.br, b)
    branch.label[branch == b] <- max(branch.label) + 1
  }
}
branch.label[branch.label == 0] <- NA

# for consistent branch naming, use expression pf some marker genes
# to order the branch IDs
new.branch.label <- factor(branch.label)
goi <- c('Tcf4', 'Ebf1')
goi.expr.avg <- sapply(goi, function(gene) aggregate(expr[gene, ], by=list(branch.label), FUN=mean)$x)
tmp <- data.frame(apply(goi.expr.avg, 2, function(x) -(x==max(x))))
levels(new.branch.label) <- do.call(order, tmp)
branch.label <- as.numeric(as.character(new.branch.label))

# append results to meta data
md$root <- (1:nrow(md)) == root
md$dfr <- dist.from.root
md$MDS1 <- mds[, 1]
md$MDS2 <- mds[, 2]
md$branch <- branch.label

# save the meta data including the branch analysis results
saveRDS(md, file = sprintf('results/%s_branch_analysis_meta_data.Rds', result.bn))

# visualize result
pdf(sprintf('results/%s_branch_analysis.pdf', result.bn), width = 7, height = 5, useDingbats = FALSE)
g <- ggplot(md, aes(MDS1, MDS2)) + geom_point(aes(color=rank(dfr)), size=2, shape=16) + 
  scale_color_gradientn(colours=my.cols.RYG, name='Distance from\nroot (ranked)') + 
  geom_density_2d(color='black', size=0.5, alpha=0.4) +
  theme_bw(base_size=12) + labs(x='MDS 1', y='MDS 2') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks=element_blank(), axis.text=element_blank())
plot(g)
g <- ggplot(md, aes(MDS1, MDS2)) + geom_point(aes(color=factor(branch)), size=2, shape=16) + 
  scale_color_discrete(name='Branch') + 
  theme_bw(base_size=12) + labs(x='MDS 1', y='MDS 2') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks=element_blank(), axis.text=element_blank())
plot(g)
g <- ggplot(md, aes(MDS1, MDS2)) + geom_point(aes(color=dfr), size=2, shape=16) + 
  scale_color_gradientn(colours=my.cols.RYG, name='Distance from\nroot') + 
  geom_density_2d(color='black', size=0.5, alpha=0.4) +
  theme_bw(base_size=12) + labs(x='MDS 1', y='MDS 2') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks=element_blank(), axis.text=element_blank())
plot(g)
g <- ggplot(md, aes(MDS1, MDS2)) + geom_point(aes(color=rank(maturation.score.smooth)), size=2, shape=16) + 
  scale_color_gradientn(colours=my.cols.RYG, name='Maturation score\n(ranked)') + 
  geom_density_2d(color='black', size=0.5, alpha=0.4) +
  theme_bw(base_size=12) + labs(x='MDS 1', y='MDS 2') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks=element_blank(), axis.text=element_blank())
plot(g)
g <- ggplot(md, aes(MDS1, MDS2)) + geom_point(aes(color=maturation.score.smooth), size=2, shape=16) + 
  scale_color_gradientn(colours=my.cols.RYG, name='Maturation\nscore') + 
  geom_density_2d(color='black', size=0.5, alpha=0.4) +
  theme_bw(base_size=12) + labs(x='MDS 1', y='MDS 2') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks=element_blank(), axis.text=element_blank())
plot(g)
dev.off()


# run differential expression test for each branch vs the oters
# keep track of genes higher in branch
branches <- sort(unique(md$branch))
de.genes <- list()
for (br in branches) {
  de.out <- de.nb.reg.helper(cm, md, c('reads', 'mols.per.gene', 'cc', 'sample.name'), 'branch', br, setdiff(branches, br))
  de.genes[[br]] <- subset(de.out, adj.pval < 0.001 & log.fc < -1 & freq1 >= 0.1)$gene
}

# for the top marker genes, average expression per branch and show on heatmap
de.genes.short <- unique(unlist(lapply(de.genes, function(x) x[1:10])))
tmp <- t(apply(expr[de.genes.short, ], 1, function(x) aggregate(x, by=list(branch=md$branch), FUN=mean)$x))
colnames(tmp) <- paste('Branch', 1:ncol(tmp))
tmp <- t(scale(t(tmp)))
tmp[tmp > 1.5] <- 1.5
tmp[tmp < -1.5] <- -1.5
pdf(sprintf('results/%s_branch_analysis_top_de_genes.pdf', result.bn), width = 5, height = 7, useDingbats = FALSE, pointsize = 12)
hm <- gplots::heatmap.2(tmp, scale = 'none', trace='none', col=colorRampPalette(c("#998ec3", "black", "#f1a340"))(31), 
                distfun = dist.p.cor, hclustfun = my.hclust.av,  Colv=NA, Rowv = NA, dendrogram = 'none',
                cexCol = 12/12, cexRow = 12/12, density.info = 'none',
                key.title = NA, key.xlab = 'Row-scaled avg expression')
dev.off()

branch.markers <- data.frame(gene=rownames(tmp), branch=apply(tmp, 1, which.max))
write.csv(branch.markers, row.names = FALSE,
          file = sprintf('results/%s_branch_analysis_top_marker_genes.csv', result.bn))

# save the following for later
# full set of differential expression results
de.genes.union <- unique(unlist(de.genes))
saveRDS(de.genes.union, file = sprintf('results/%s_branch_analysis_de_genes.Rds', result.bn))
# branch-averaged expression
expr.branch.avg <- t(apply(expr, 1, function(x) aggregate(x, by=list(branch=md$branch), FUN=mean)$x))
colnames(expr.branch.avg) <- paste('Branch', sort(unique(md$branch)))
saveRDS(expr.branch.avg, file = sprintf('results/%s_branch_analysis_expr_branch_avg.Rds', result.bn))
