
##' @rdname eigencorval
##' @importFrom PCAtools::eigencorplot
##' Extract correlation data from PCAtools::eigencorplot.
##' @params Check PCAtools::eigencorplot for description
##' @return list with the correlation matrix and corresponding p-values
##'
eigencorval <- function (
    pcaobj,
    components = getComponents(pcaobj, seq_len(10)),
    metavars,
    scale = TRUE, main = "",
    corFUN = "pearson",
    corUSE = "pairwise.complete.obs",
    corMultipleTestCorrection = "none"
    )
{
  data <- pcaobj$rotated
  metadata <- pcaobj$metadata
  for (i in seq_len(length(components))) {
    if (!is.numeric(data[, components[i]])) {
      warning(components[i], " is not numeric - please check the source data",
              " as everything will be converted to a matrix")
    }
  }
  for (i in seq_len(length(metavars))) {
    if (!is.numeric(metadata[, metavars[i]])) {
      warning(metavars[i], " is not numeric - please check the source data",
              " as non-numeric variables will be coerced to numeric")
    }
  }
  xvals <- data.matrix(data[, which(colnames(data) %in% components),
                            drop = FALSE])
  yvals <- metadata[, which(colnames(metadata) %in% metavars),
                    drop = FALSE]
  chararcter_columns = unlist(lapply(yvals, is.numeric))
  chararcter_columns = !chararcter_columns
  chararcter_columns = names(which(chararcter_columns))
  for (c in chararcter_columns) {
    yvals[, eval(quote(c))] = as.numeric(as.factor(yvals[,
                                                         eval(quote(c))]))
  }
  yvals <- data.matrix(yvals)
  corvals <- cor(xvals, yvals, use = corUSE, method = corFUN)

  N <- ncol(xvals) * ncol(yvals)
  pvals <- data.frame(pval = numeric(N), i = numeric(N), j = numeric(N))
  k <- 0
  for (i in seq_len(ncol(xvals))) {
    for (j in seq_len(ncol(yvals))) {
      k <- k + 1
      pvals[k, "pval"] <- cor.test(xvals[, i], yvals[,
                                                     j], use = corUSE, method = corFUN)$p.value
      pvals[k, "i"] <- colnames(xvals)[i]
      pvals[k, "j"] <- colnames(yvals)[j]
    }
  }
  if (corMultipleTestCorrection != "none") {
    pvals$pval <- p.adjust(pvals$pval, method = corMultipleTestCorrection)
  }
  pvals <- reshape2::dcast(pvals, i ~ j, value.var = "pval")
  rownames(pvals) <- pvals$i
  pvals$i <- NULL
  pvals <- pvals[match(rownames(corvals), rownames(pvals)),
  ]
  pvals <- pvals[colnames(corvals)]

  return(
    list(
      corvals = corvals,
      pvals = pvals
    )
  )
}
