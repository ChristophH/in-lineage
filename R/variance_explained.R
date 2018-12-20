# this script performs the following steps:
# 1) load mitotic cells from CGE and MGE dropseq experiments
# 2) quantify variance explained by individual factors
# 3) load E14.5 10x data and subset to postmitotic cells
# 4) quantify variance explained by individual factors

# these files are created in results directory:
# - variance_explained_dropseq_mitotic_CGE_MGE.pdf
# - variance_explained_10x_E14.5_CGE_MGE.pdf

# running times on an Intel Xeon Processor E5-2697 v3 @ 2.6 to 3.6 GHz
# (using 6 cores for some of the steps)
# ca. 40 minutes

# support functions are defined in library
source('R/lib.R')

set.seed(42)
options(mc.cores = 6)

# set basename for loading and saving results
input.bn <- 'all_samples'

# load maturation trajectory results and keep only mitotic cells from CGE and MGE
md <- readRDS(file = sprintf('results/%s_maturation_trajectory_meta_data.Rds', input.bn))
cm <- readRDS('data/dropseq_digitial_expression.Rds')

md <- subset(md, !postmitotic & eminence %in% c('CGE', 'MGE'))
cm <- cm[, rownames(md)]
cat('There are', ncol(cm), 'mitotic CGE and MGE cells\n')

# specify which factors to regress out
norm.fac <- c('reads', 'mols.per.gene')
# quantify the variance explained for the following factors
foi <- c('maturation.score.smooth', 'cc', 'eminence', 'reads', 'mols')
varex.out <- varex(cm, md, norm.fac, foi)

# how much variance does the first PC explain
pc1.varex <- varex.out$x.pca.sdev[1]^2/sum(varex.out$x.pca.sdev^2)

tmp <- sort(unlist(varex.out$fraction.explained), decreasing = TRUE)
df2 <- data.frame(fac=factor(names(tmp), ordered=TRUE, levels=names(tmp)), var.expl=tmp)
g2 <- ggplot(df2) + geom_bar(aes(x=fac, y=var.expl/pc1.varex), stat='identity') +
  ylab('Rel. variance explained') + xlab('Factor') +  theme_bw(base_size = 10) +
  theme(legend.position="none") + 
  theme(panel.grid.major = element_line(colour = "grey90", size = 0.3), panel.grid.minor = element_blank()) +
  theme(panel.grid.major.x = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
ggsave(g2, filename = 'results/variance_explained_dropseq_mitotic_CGE_MGE.pdf', width = 4, height = 8 , units = 'cm')


# do the same analysis for the E14.5 10x cells

# get the branch assignments from the mapping analysis
map.cge <- readRDS('results/CGE_E13.5_Lhx6neg_mapped.Rds')
map.mge <- readRDS('results/MGE_E13.5_Lhx6pos_mapped.Rds')
map.both <- rbind(map.cge, map.mge)

cm <- cbind(readRDS('data/10x_digitial_expression_E13.Rds'),
            readRDS('data/10x_digitial_expression_E18.Rds'),
            readRDS('data/10x_digitial_expression_P10.Rds'))
md <- rbind(readRDS('data/10x_meta_data_E13.Rds'),
            readRDS('data/10x_meta_data_E18.Rds'),
            readRDS('data/10x_meta_data_P10.Rds'))

sel <- intersect(rownames(md), rownames(map.both))
md <- md[sel, ]
cm <- cm[, rownames(md)]
cc <- get.cc.score(cm, seed=42)
md$cc <- cc$score
md$branch <- factor(map.both[sel, 'label'])

norm.fac <- c('reads', 'mols.per.gene')
# quantify the variance explained for the following factors
foi <- c('branch', 'cc', 'eminence', 'reads', 'mols')
varex.out <- varex(cm, md, norm.fac, foi)

# how much variance does the first PC explain
pc1.varex <- varex.out$x.pca.sdev[1]^2/sum(varex.out$x.pca.sdev^2)

tmp <- sort(unlist(varex.out$fraction.explained), decreasing = TRUE)
df2 <- data.frame(fac=factor(names(tmp), ordered=TRUE, levels=names(tmp)), var.expl=tmp)
g2 <- ggplot(df2) + geom_bar(aes(x=fac, y=var.expl/pc1.varex), stat='identity') +
  ylab('Rel. variance explained') + xlab('Factor') +  theme_bw(base_size = 10) +
  theme(legend.position="none") + 
  theme(panel.grid.major = element_line(colour = "grey90", size = 0.3), panel.grid.minor = element_blank()) +
  theme(panel.grid.major.x = element_blank(), axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
ggsave(g2, filename = 'results/variance_explained_10x_E13.5_CGE_MGE.pdf', width = 4, height = 8 , units = 'cm')

