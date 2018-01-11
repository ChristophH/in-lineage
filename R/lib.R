
library('Matrix')
library('parallel')
library('MASS')
library('diffusionMap')
library('FNN')
library('igraph')
library('princurve')
library('ggplot2')
library('inline')
library('gplots')


###############################################################################
# for cell cycle score

get.bg.lists <- function(goi, N, expr.bin) {
  res <- list()
  goi.bin.tab <- table(expr.bin[goi])
  for (i in 1:N) {
    res[[i]] <- unlist(lapply(names(goi.bin.tab), function(b) {
      sel <- which(expr.bin == as.numeric(b) & !(names(expr.bin) %in% goi))
      sample(names(expr.bin)[sel], goi.bin.tab[b])
    }))
  }
  return(res)
}

enr.score <- function(expr, goi, bg.lst) {
  goi.mean <- apply(expr[goi, ], 2, mean)
  bg.mean <- sapply(1:length(bg.lst), function(i) apply(expr[bg.lst[[i]], ], 2, mean))
  return((goi.mean - apply(bg.mean, 1, mean)) / apply(bg.mean, 1, sd))
}


get.cc.score <- function(cm, N=100, seed=42) {
  set.seed(seed)
  cat('get.cc.score, ')
  cat('number of random background gene sets set to', N, '\n')
  
  min.cells <- 5
  
  cells.mols <- apply(cm, 2, sum)
  gene.cells <- apply(cm>0, 1, sum)
  cm <- cm[gene.cells >= min.cells, ]
  
  gene.mean <- apply(cm, 1, mean)
  
  breaks <- unique(quantile(log10(gene.mean), probs = seq(0,1, length.out = 50)))
  gene.bin <- cut(log10(gene.mean), breaks = breaks, labels = FALSE)
  names(gene.bin) <- rownames(cm)
  gene.bin[is.na(gene.bin)] <- 0
  
  regev.s.genes <- read.table(file='./annotation/s_genes.txt', header=FALSE, stringsAsFactors=FALSE)$V1
  regev.g2m.genes <- read.table(file='./annotation/g2m_genes.txt', header=FALSE, stringsAsFactors=FALSE)$V1
  
  goi.lst <- list('S'=rownames(cm)[!is.na(match(toupper(rownames(cm)), regev.s.genes))],
                  'G2M'=rownames(cm)[!is.na(match(toupper(rownames(cm)), regev.g2m.genes))])
  
  n <- min(40, min(sapply(goi.lst, length)))
  goi.lst <- lapply(goi.lst, function(x) x[order(gene.mean[x], decreasing = TRUE)[1:n]])
  
  bg.lst <- list('S'=get.bg.lists(goi.lst[['S']], N, gene.bin),
                 'G2M'=get.bg.lists(goi.lst[['G2M']], N, gene.bin))
  
  all.genes <- sort(unique(c(unlist(goi.lst, use.names=FALSE), unlist(bg.lst, use.names=FALSE))))
  
  expr <- log10(cm[all.genes, ]+1)
  
  s.score <- enr.score(expr, goi.lst[['S']], bg.lst[['S']])
  g2m.score <- enr.score(expr, goi.lst[['G2M']], bg.lst[['G2M']])
  
  phase <- as.numeric(g2m.score > 2 & s.score <= 2)
  phase[g2m.score <= 2 & s.score > 2] <- -1
  
  return(data.frame(score=s.score-g2m.score, s.score, g2m.score, phase))
}


###############################################################################
# for normalization

theta.reg <- function(cm, regressors, min.theta=0.01, bins=64) {
  b.id <- (1:nrow(cm)) %% max(1, bins, na.rm=TRUE) + 1
  cat(sprintf('get regularized theta estimate for %d genes and %d cells\n', nrow(cm), ncol(cm)))
  cat(sprintf('processing %d bins with ca %d genes in each\n', bins, round(nrow(cm)/bins, 0)))
  theta.estimate <- rep(NA, nrow(cm))
  for (bin in sort(unique(b.id))) {
    sel.g <- which(b.id == bin)
    bin.theta.estimate <- unlist(mclapply(sel.g, function(i) {
      as.numeric(theta.ml(cm[i, ], glm(cm[i, ] ~ ., data = regressors, family=poisson)$fitted))
    }), use.names = FALSE)
    theta.estimate[sel.g] <- bin.theta.estimate
    cat(sprintf('%d ', bin))
  }
  cat('done\n')
  raw.mean <- apply(cm, 1, mean)
  log.raw.mean <- log10(raw.mean)
  var.estimate <- raw.mean + raw.mean^2/theta.estimate
  
  fit <- loess(log10(var.estimate) ~ log.raw.mean, span=0.33)
  theta.fit <- raw.mean^2 / (10^fit$fitted - raw.mean)
  
  to.fix <- theta.fit <= min.theta | is.infinite(theta.fit)
  if (any(to.fix)) {
    cat('Fitted theta below', min.theta, 'for', sum(to.fix), 'genes, setting them to', min.theta, '\n')
    theta.fit[to.fix] <- min.theta
  }
  names(theta.fit) <- rownames(cm)
  return(theta.fit)
}

nb.residuals.glm <- function(y, regression.mat, fitted.theta, gene) {
  fit <- 0
  try(fit <- glm(y ~ ., data = regression.mat, family=negative.binomial(theta=fitted.theta)), silent=TRUE)
  if (class(fit)[1] == 'numeric') {
    message(sprintf('glm and family=negative.binomial(theta=%f) failed for gene %s; falling back to scale(log10(y+1))', 
                    fitted.theta, gene))
    return(scale(log10(y+1))[, 1])
  }
  return(residuals(fit, type='pearson'))
}

norm.nb.reg <- function(cm, regressors, min.theta=0.01, bins=64, theta.fit=NA, pr.th=NA, save.theta.fit=c()) {
  cat('Normalizing data using regularized NB regression\n')
  cat('explanatory variables:', colnames(regressors), '\n')
  if (any(is.na(theta.fit))) {
    theta.fit <- theta.reg(cm, regressors, min.theta, bins)
    if (is.character(save.theta.fit)) {
      save(theta.fit, file=save.theta.fit)
    }
  }
  
  b.id <- (1:nrow(cm)) %% max(1, bins, na.rm=TRUE) + 1
  cat('Running NB regression\n')
  res <- matrix(NA, nrow(cm), ncol(cm), dimnames=dimnames(cm))
  for (bin in sort(unique(b.id))) {
    sel.g <- rownames(cm)[b.id == bin]
    expr.lst <- mclapply(sel.g, function(gene) nb.residuals.glm(cm[gene, ], regressors, theta.fit[gene], gene), mc.preschedule = TRUE)
    res[sel.g, ] <- do.call(rbind, expr.lst)
    cat(sprintf('%d ', bin))
  }
  cat('done\n')
  if (!any(is.na(pr.th))) {
    res[res > pr.th] <- pr.th
    res[res < -pr.th] <- -pr.th
  }
  attr(res, 'theta.fit') <- theta.fit
  return(res)
}

###############################################################################
# for clustering

dim.red <- function(expr, max.dim, ev.red.th, plot.title=NA, do.scale.result=FALSE) {
  cat('Dimensionality reduction via diffusion maps using', nrow(expr), 'genes and', ncol(expr), 'cells\n')
  if (sum(is.na(expr)) > 0) {
    dmat <- 1 - cor(expr, use = 'pairwise.complete.obs')
  } else {
    dmat <- 1 - cor(expr)
  }
  
  max.dim <- min(max.dim, nrow(dmat)/2)
  dmap <- diffuse(dmat, neigen=max.dim, maxdim=max.dim)
  ev <- dmap$eigenvals
  
  ev.red <- ev/sum(ev)
  evdim <- rev(which(ev.red > ev.red.th))[1]
  
  if (is.character(plot.title)) {
    plot(ev, ylim=c(0, max(ev)), main = plot.title)
    abline(v=evdim + 0.5, col='blue')
  }
  
  evdim <- max(2, evdim, na.rm=TRUE)
  cat('Using', evdim, 'significant DM coordinates\n')
  
  colnames(dmap$X) <- paste0('DMC', 1:ncol(dmap$X))
  res <- dmap$X[, 1:evdim]
  if (do.scale.result) {
    res <- scale(dmap$X[, 1:evdim])
  } 
  return(res)
}

# jaccard similarity
# rows in 'mat' are cells
jacc.sim <- function(mat, k) {
  # generate a sparse nearest neighbor matrix
  nn.indices <- get.knn(mat, k)$nn.index
  j <- as.numeric(t(nn.indices))
  i <- ((1:length(j))-1) %/% k + 1
  nn.mat <- sparseMatrix(i=i, j=j, x=1)
  rm(nn.indices, i, j)
  # turn nn matrix into SNN matrix and then into Jaccard similarity
  snn <- nn.mat %*% t(nn.mat)
  snn.summary <- summary(snn)
  snn <- sparseMatrix(i=snn.summary$i, j=snn.summary$j, x=snn.summary$x/(2*k-snn.summary$x))
  rm(snn.summary)
  return(snn)
}


cluster.the.data.simple <- function(cm, expr, k, sel.g=NA, min.mean=0.001, 
                                    min.cells=3, z.th=1, ev.red.th=0.02, seed=NULL, 
                                    max.dim=50) {
  if (all(is.na(sel.g))) {
    # no genes specified, use most variable genes
    goi <- rownames(expr)[apply(cm[rownames(expr), ]>0, 1, sum) >= min.cells & apply(cm[rownames(expr), ], 1, mean) >= min.mean]
    sspr <- apply(expr[goi, ]^2, 1, sum)
    sel.g <- goi[scale(sqrt(sspr)) > z.th]
  }
  cat(sprintf('Selected %d variable genes\n', length(sel.g)))
  sel.g <- intersect(sel.g, rownames(expr))
  cat(sprintf('%d of these are in expression matrix.\n', length(sel.g)))
  
  if (is.numeric(seed)) {
    set.seed(seed)
  }
  
  dm <- dim.red(expr[sel.g, ], max.dim, ev.red.th, do.scale.result = TRUE)
  
  sim.mat <- jacc.sim(dm, k)
  
  gr <- graph_from_adjacency_matrix(sim.mat, mode='undirected', weighted=TRUE, diag=FALSE)
  cl <- as.numeric(membership(cluster_louvain(gr)))
  
  results <- list()
  results$dm <- dm
  results$clustering <- cl
  results$sel.g <- sel.g
  results$sim.mat <- sim.mat
  results$gr <- gr
  cat('Clustering table\n')
  print(table(results$clustering))
  return(results)
}

######################################################################
# for maturation trajectory

# fit maturation trajectory
maturation.trajectory <- function(cm, md, expr) {
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
  
  md <- md[, !grepl('^DMC', colnames(md))]
  md <- cbind(md, rd)
  
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
  md$in.cc.phase <- md$cc.phase != 0
  fit <- loess(as.numeric(md$in.cc.phase) ~ md$maturation.score.smooth, span=0.5, degree=2)
  md$cc.phase.fit <- fit$fitted
  # pick MT threshold based on drop in cc.phase cells
  # ignore edges of MT because of potential outliers
  mt.th <- max(subset(md, cc.phase.fit > mean(md$in.cc.phase)/2 & maturation.score.smooth >= 0.2 & maturation.score.smooth <= 0.8)$maturation.score.smooth)
  
  md$postmitotic <- md$maturation.score.smooth > mt.th
  return(list(md=md, pricu=pricu, mt.th=mt.th))
}


# for smoothing maturation score

nn.smooth <- function(y, coords, k) {
  knn.out <- FNN::get.knn(coords, k)
  w <- 1 / (knn.out$nn.dist+.Machine$double.eps)
  w <- w / apply(w, 1, sum)
  v <- apply(knn.out$nn.index, 2, function(i) y[i])
  return(apply(v*w, 1, sum))
}

# maturation score colors
my.cols.RYG <- colorRampPalette(c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee08b",
                                  "#ffffbf", "#d9ef8b", "#a6d96a", "#66bd63", "#1a9850", "#006837"))(11)


# for differential expression testing

de.nb.reg <- function(y, theta, md, com.fac, grp.fac) {
  fit2 <- 0
  fit4 <- 0
  try(fit2 <- glm(y ~ ., data=md[, com.fac, drop=FALSE], family=negative.binomial(theta=theta)), silent=TRUE)
  try(fit4 <- glm(y ~ ., data=md[, c(com.fac, grp.fac)], family=negative.binomial(theta=theta)), silent=TRUE)
  tab <- as.matrix(table(y > 0, md[, grp.fac]))
  freqs <- tab['TRUE', ] / apply(tab, 2, sum)
  if (class(fit2)[1] == 'numeric' | class(fit4)[1] == 'numeric') {
    message('One of the glm.nb calls failed')
    return(c(rep(NA, 5), freqs))
  }
  pval <- anova(fit2, fit4, test='Chisq')$'Pr(>Chi)'[2]
  grp.coef <- coef(fit4)[grp.fac]
  log.fc <- log2(exp(grp.coef)) #log.fc <- log2(1/exp(coef(fit4)[foi]))
  #print(coef(fit4))
  #print(grp.fac)
  #browser()
  return(c(fit2$deviance, fit4$deviance, pval, grp.coef, log.fc, freqs))
}
de.nb.reg.helper <- function(cm, md, com.fac, grp.fac, grp.val1, grp.val2, min.frac=0.02, bins=32) {
  sel.1 <- which(md[, grp.fac] %in% grp.val1)
  sel.2 <- which(md[, grp.fac] %in% grp.val2)
  goi <- (apply(cm[, sel.1] > 0, 1, mean) >= min.frac | apply(cm[, sel.2] > 0, 1, mean) >= min.frac)
  cm <- cm[goi, c(sel.1, sel.2)]
  tmp.md <- md[c(sel.1, sel.2), ]
  #tmp.md$group <- factor(as.character(tmp.md[, grp.fac] %in% grp.val1), ordered=TRUE, levels=c('TRUE', 'FALSE'))
  tmp.md$de.nb.reg.group <- c(rep(0, length(sel.1)), rep(1, length(sel.2)))
  cat(sprintf('Running differential expression test on %d genes and %d + %d cells.\n', nrow(cm), length(sel.1), length(sel.2)))
  theta.fit <- theta.reg(cm, tmp.md[, com.fac, drop=FALSE], bins=bins)
  tmp <- mclapply(1:nrow(cm), function(i) de.nb.reg(cm[i, ], theta.fit[i], tmp.md, com.fac, 'de.nb.reg.group'))
  res <- as.data.frame(matrix(unlist(tmp), ncol=length(tmp[[1]]), byrow=TRUE))
  dimnames(res) <- list(rownames(cm), c('dev1', 'dev2', 'pval', 'coef', 'log.fc', 'freq1', 'freq2'))
  res$adj.pval <- p.adjust(res$pval, method='fdr')
  res <- res[order(res$pval, -abs(res$log.fc)), ]
  res$gene <- rownames(res)
  return(res)
}


###############################################################################
# for binned mutual information

nmi.cont <- function(a, b, nbins=10) {
  a <- discretize(a, nbins)
  b <- discretize(b, nbins)
  n <- length(a)
  cooc <- as.matrix(table(a, b))
  cooc[cooc == 0] <- NA
  arf <- apply(cooc / n, 1, sum, na.rm=TRUE)
  brf <- apply(cooc / n, 2, sum, na.rm=TRUE)
  tab2 <- arf %*% t(brf)
  iab <- sum((cooc / n) * log((cooc / n) / tab2), na.rm=TRUE)
  ha <- -sum(arf * log(arf))
  hb <- -sum(brf * log(brf))
  return(iab/max(ha, hb))
}

discretize <- function(X, nbins) {
  N <- length(X)
  X.min <- min(X)
  X.max <- max(X)
  tiny <- max(.Machine$double.eps * (X.max - X.min), .Machine$double.eps)
  X.disc <- floor((X - X.min) / (X.max - X.min + tiny) * nbins)
  return(as.integer(X.disc))
}


# compute MI of each colum in x to each column in y
# you'll want the larger matrix to be x
mi <- function(x, y, nbins=10) {
  cpu.n <- as.numeric(options('mc.cores'))
  x <- as.matrix(x)
  y <- as.matrix(y)
  # discretize the columns
  x <- apply(x, 2, discretize, nbins)
  y <- apply(y, 2, discretize, nbins)
  
  m <- ncol(x)
  n <- ncol(y)
  ret <- matrix(0, m, n, dimnames=list(colnames(x), colnames(y)))
  
  if (cpu.n == 1 | cpu.n > m) {  
    s <- nrow(x)
    return(matrix(mi.in.c4(s, nbins, x, m, y, n, ret)$mi, m, dimnames=list(colnames(x), colnames(y))))
  }
  
  col.list <- split(1:m, cut(1:m, cpu.n, labels=F))
  ret.list <- mclapply(col.list, mi.thread, x, y, nbins)
  for (i in 1:length(ret.list)) {
    ret[col.list[[i]], ] <- ret.list[[i]]
  }
  return(ret)
}


mi.thread <- function(cols, x, y, nbins) {
  x <- x[, cols, drop = FALSE]
  m <- ncol(x)
  n <- ncol(y)
  s <- nrow(x)
  ret <- matrix(0, m, n)
  return(matrix(mi.in.c4(s, nbins, x, m, y, n, ret)$mi, m))
}

mi.in.c4.sig <- signature(n = 'integer', k = 'integer', x = 'integer', 
                          ncolx = 'integer', y = 'integer', 
                          ncoly = 'integer', mi = 'numeric')
mi.in.c4.code <- '
  int i, j, cx, cy;
  double *pxy, *px, *py;

  py = (double *)calloc(*k * *ncoly, sizeof(double));
  for (cy = 0; cy < *ncoly; ++cy) {
    for (i = 0; i < *n; ++i) {
      ++py[cy * *k + y[cy * *n + i]];
    }
  }

  for (cx = 0; cx < *ncolx; ++cx) {
    for (cy = 0; cy < *ncoly; ++cy) {
      pxy = (double *)calloc(*k * *k, sizeof(double));
      px = (double *)calloc(*k, sizeof(double));

      for (i = 0; i < *n; ++i) {
        ++pxy[(y[cy * *n + i]) * (*k) + (x[cx * *n + i])];
        ++px[x[cx * *n + i]];
      }
      for (i = 0; i < *k; ++i) {
        for (j = 0; j < *k; ++j) {
          if (pxy[j * (*k) + i] > 0) {
            mi[cy * *ncolx + cx] += pxy[j * (*k) + i] / *n * 
                                    log(pxy[j * (*k) + i] / *n / 
                                       (px[i] / *n * py[cy * *k + j] / *n));
          }
        }
      }
      free(px);
      free(pxy);
    }
  }
  free(py);'
mi.in.c4 <- cfunction(mi.in.c4.sig, mi.in.c4.code, convention=".C")


###############################################################################
# for smooth expression as function of maturation score
smooth.expr <- function(mat, x, x.pred) {
  n <- length(x.pred)
  nd <- data.frame(x=x.pred)
  tmp <- mclapply(1:nrow(mat), function(i) {
    fit <- loess(mat[i, ] ~ x, span=0.5, degree=2)
    predict(fit, newdata=nd, se=FALSE)
  })
  expr.fit <- t(matrix(unlist(tmp), n))
  rownames(expr.fit) <- rownames(mat)
  colnames(expr.fit) <- nd$x
  return(expr.fit)
}


###############################################################################
# for branch analysis

mst.avg <- function(dmap.dmat, n, ss, sp.weights=NULL) {
  cat('Bootstrapped minimum spanning trees will be based on random bootstraps of size', ss, 'cells\n')
  cat('Will sample until every pair of cells has been sampled at least', n, 'times\n')
  cat('Current bootstrap: ')
  dmat.sum <- dmap.dmat * 0
  adj.mat.sum <- dmap.dmat * 0
  count.mat <- dmap.dmat * 0
  sel.list <- list()
  #g.list <- list()
  
  i <- 0
  while (sum(count.mat < n) > 0) {
    i <- i + 1
    cat(sprintf('%d ', i))
    
    p <- exp(-(count.mat+1))
    vals <- sample.int(length(count.mat), ss, replace=TRUE, prob=p)
    sel <- sort(unique(c(vals %/% nrow(count.mat) + 1, vals %% nrow(count.mat)))[1:ss])
    
    sel.list[[i]] <- sel
    
    g <- igraph::mst(igraph::graph_from_adjacency_matrix(dmap.dmat[sel, sel], mode='undirected', weighted=TRUE, add.colnames=NA))
    #g.list[[i]] <- g
    
    sp <- igraph::shortest.paths(g, weights=sp.weights)
    
    dmat.sum[sel, sel] <- dmat.sum[sel, sel] + sp
    adj.mat.sum[sel, sel] <- adj.mat.sum[sel, sel] + igraph::get.adjacency(g, sparse = FALSE)
    
    count.mat[sel, sel] <- count.mat[sel, sel] + 1
    #cat('count mat zeros', sum(count.mat == 0), '\n')
  }
  cat('\n')
  # average all distances
  dmats.avg <- dmat.sum / count.mat
  adj.mat.avg <- adj.mat.sum / count.mat
  
  return(list(dmats.avg=dmats.avg, adj.mat.avg=adj.mat.avg, sel.list=sel.list))#, g.list=g.list))
}


tree.to.dag <- function(tree, root, normalize.weights=TRUE) {
  tmp <- igraph::unfold.tree(igraph::as.directed(tree), mode='out', roots=root)  # unfortunately this removes the weight attribute
  dag <- igraph::delete.vertices(tmp$tree, which(duplicated(tmp$vertex_index)))
  di.we <- igraph::get.adjacency(tree, type='both', attr='weight', sparse=FALSE)
  di.we[igraph::get.adjacency(dag, sparse=FALSE) == 0] <- 0
  if (normalize.weights) {
    di.we <- di.we / sum(di.we)
  }
  dag <- igraph::graph_from_adjacency_matrix(di.we, mode='directed', weighted=TRUE)
  return(dag)
}

# traverse directed acyclic graph and detect branching events
branch.id <- function(dag, w.th=0.2, weighted=FALSE, max.splits=Inf) {
  if (weighted) {
    dag.amat <- igraph::get.adjacency(dag, type='both', attr='weight', sparse=FALSE)
  } else {
    dag.amat <- igraph::get.adjacency(dag, type='both', attr=NULL, sparse=FALSE)
  }
  dag.amat <- dag.amat / sum(dag.amat)
  in.weight <- apply(dag.amat, 2, sum, na.rm=TRUE)
  dag.sp <- igraph::shortest.paths(dag, weights=NULL, mode='out')
  dag.sp[is.infinite(dag.sp)] <- NA
  branch <- rep(1, nrow(dag.amat))
  split.pot <- c()
  split.cnt <- 0
  cum.dist <- 0
  to.visit <- which(apply(dag.amat == 0, 2, all))
  while (length(to.visit) > 0 & split.cnt < max.splits) {
    kids <- which(dag.amat[to.visit[1], ] > 0)
    #cat(to.visit[1], '\t', to.visit.br[1], '\t', kids)
    weights <- c()
    for (k in kids) {
      #weights <- c(weights, sum(dag.sp[k, ] + dag.sp[to.visit[1], k], na.rm=TRUE) + dag.sp[to.visit[1], k])
      weights <- c(weights, sum(in.weight[c(k, which(dag.sp[k, ] > 0))]))
    }
    if (length(weights) > 1) {
      split.pot <- c(split.pot, sort(weights, decreasing=TRUE)[2])
    }
    if (sum(weights > w.th) > 1) {
      split.cnt <- split.cnt + 1
      for (k in kids[weights > w.th]) {
        branch[c(k, which(dag.sp[k, ] > 0))] <- max(branch) + 1
      }
    }
    #cat('\t', weights, '\t', new.branch, '\n')
    to.visit <- c(to.visit[-1], kids)
  }
  print(head(sort(split.pot, decreasing=TRUE)))
  return(branch)
}


###############################################################################
# for variance explained analysis

frac.var.expl <- function(fac.vec, x, x.pca) {
  fac.vec <- fac.vec - mean(fac.vec)
  fac.vec <- fac.vec / sqrt(sum(fac.vec^2))
  proj <- t(t(x - x.pca$center)) %*% fac.vec
  return(var(proj)/sum(x.pca$sdev^2))
}

varex <- function(cm, md, norm.fac, foi, x=c()) {
  genes <- rownames(cm)[apply(cm > 0, 1, sum) >= 5]
  cat('Normalizing', length(genes), 'genes that are present in at least 5 cells\n')
  expr <- norm.nb.reg(cm[genes, ], md[, norm.fac, drop=FALSE], pr.th = 30)
  vg <- rownames(expr)[scale(sqrt(apply(expr^2, 1, sum))) > 1]
  x <- expr[vg, ]
  x.pca <- prcomp(x, center = TRUE, scale. = FALSE)
  fraction.explained <- list()
  for (fac in foi) {
    fac.vec <- md[, fac]
    if (!is.numeric(fac.vec)) {
      fac.vec <- as.factor(fac.vec)
      fac.mat <- as.matrix(model.matrix(~ 1 + fac.vec))
      tmp.pca <- prcomp(fac.mat, center = TRUE, scale. = FALSE)
      fraction.explained[[fac]] <- sum(apply(tmp.pca$x, 2, function(fac.vec) frac.var.expl(fac.vec, x, x.pca)), na.rm = TRUE)
    } else {
      fraction.explained[[fac]] <- frac.var.expl(fac.vec, x, x.pca)
    }
  }
  return(list(x.pca.sdev=x.pca$sdev, fraction.explained=fraction.explained))
}
