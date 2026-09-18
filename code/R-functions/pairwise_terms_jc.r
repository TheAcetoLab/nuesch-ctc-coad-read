#' Function to calculate pairwise jaccard similarity index between terms in a
#' gene set collection list. Returns a similarity matrix
#' @param gset_collection GeneSetCollection object as obtained from
#' GSEABase::getGmt or a list of gene sets
#' @examples
#' res <- nmf(mat, rank = 10, nrun = 1, seed = 'ica', method = 'nsNMF')
#' modules <-  NMFToModules(res, gmin = 5)
#'
pairwise_terms_jc <- function(gset_collection){
  terms_comb <- combn(names(gset_collection), 2, simplify = FALSE)
  jc_df <- foreach(i=terms_comb, .combine = rbind) %do% {
    if(class(gset_collection)=="GeneSetCollection"){
      jc_i <- jaccard(
        gset_collection[[i[1]]]@geneIds,
        gset_collection[[i[2]]]@geneIds
      )
    } else if (class(gset_collection)=='list') {
      jc_i <- jaccard(
        gset_collection[[i[1]]],
        gset_collection[[i[2]]]
      )
    }

    data.frame(
      a = i,
      b = rev(i),
      jc = c(jc_i, jc_i)
    )
  }
  res <- jc_df %>%
    pivot_wider(names_from = b, values_from = jc) %>%
    column_to_rownames("a")
  res <- res[, rownames(res)]
  res[is.na(res)] <- 1
  return(as.matrix(res))
}


jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}
