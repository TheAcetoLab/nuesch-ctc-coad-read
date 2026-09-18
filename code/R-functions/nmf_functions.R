
#' Obtained from https://github.com/yanailab/PanCancer. See also
#'  [Barkley et al. (2022)](https://www.nature.com/articles/s41588-022-01141-9)
#' Extracts a list of modules from a NMF object.  To define nonoverlapping gene modules,
#' a gene ranking algorithm was implemented92. Beginning with the matrix of the contribution
#' of genes (rows) to the factors (columns), two ranking matrices were constructed: list 1 ranking
#' the gene contributions to each factor, and list 2 ranking for each gene the factors to which it
#' contributes. For each factor, genes were added in the order of their contribution (list 1), until a
#' gene was reached that contributed more to another factor, that is, its rank across factors (list 2)
#' was not 1. Factors that yielded fewer than five genes were removed.
#' @param dat A NMFfit object as a result of running NMF::nmf
#' @param gmin minimum number of genes in a factor. Factors with less genes than gmin will be removed
#'
#' @examples
#' res <- nmf(mat, rank = 10, nrun = 1, seed = 'ica', method = 'nsNMF')
#' modules <-  NMFToModules(res, gmin = 5)
#'
NMFToModules = function(
    res,
    gmin = 5
){

  scores = basis(res)
  coefs = coefficients(res)

  # Remove if fewer than gmin genes
  ranks_x = t(apply(-t(t(scores) / apply(scores, 2, mean)), 1, rank))
  ranks_y = apply(-t(t(scores) / apply(scores, 2, mean)), 2, rank)
  for (i in 1:ncol(scores)){
    ranks_y[ranks_x[,i] > 1,i] = Inf
  }
  modules = apply(ranks_y, 2, function(m){
    a = sort(m[is.finite(m)])
    a = a[a == 1:length(a)]
    names(a)
  })
  l = sapply(modules, length)
  keep = (l >= gmin)
  scores = scores[, keep]
  coefs = coefs[keep, ]

  # Find modules
  ranks_x = t(apply(-t(t(scores) / apply(scores, 2, mean)), 1, rank))
  ranks_y = apply(-t(t(scores) / apply(scores, 2, mean)), 2, rank)
  for (i in 1:ncol(scores)){
    ranks_y[ranks_x[,i] > 1,i] = Inf
  }
  modules = apply(ranks_y, 2, function(m){
    a = sort(m[is.finite(m)])
    a = a[a == 1:length(a)]
    names(a)
  })

  names(modules) = sapply(modules, '[', 1)
  names(modules) = paste('m', names(modules), sep = '_')
  names(modules) = gsub('-','_',names(modules))

  return(modules)
}
