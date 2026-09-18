#' Function to generate a volcano plot from a set of p-values and fold changes
#' @description Taken from StackOverflow.
#'
#' @param x
#' @param y
#' @param rasterize Rasterize plot using ggrastr
#'
#'
#' @author Francesc Castro-Giner, \email{fcastroginer@proton.me}
#'
#' @export
#'
#' @examples
#'

plot.volcano <- function(x, y, label = NULL, annot = NULL, y.lab = expression(paste("-lo", g[10],"(Adjusted P value)")), x.lab = expression(paste("lo", g[2],"(Fold change)")),
                         pval.thr = 0.05, top.labels = 10, selected.labels = NULL, point.size = 2, point.alpha = 0.3,
                         point.stroke = 0.01, point.shape = 21, point.colour = 'white', show.legend = FALSE,
                         label.size = 2, fill.palette = c("gray60", "red"), hline.type = 2, hline.size = 1, hline.color =  'grey60',
                         fill.by = "is.significant", fill.title = NULL, symmetric = TRUE,
                         interactive = FALSE, rasterize = FALSE, rastr.dpi = 600){
  require(dplyr)
  require(ggplot2)
  require(ggrepel)

  # Select ggplot2 point function
  if(rasterize){
    require(ggrastr)
    gg_geom_point <- ggrastr::geom_point_rast
    options(ggrastr.default.dpi = rastr.dpi)
  } else {
    gg_geom_point <- ggplot2::geom_point
  }

  # Prepare dataframe
  dat <- cbind(x = x, y = y) %>%
    data.frame %>%
    dplyr::mutate(
      is.significant = y <= pval.thr,
      y = -log10(y)
    )

  # Add annotation
  if(!is.null(annot))
    dat <- cbind(dat, annot = annot)

  # Add labels
  if(!is.null(label)){
    dat <- cbind(dat, label = label) %>%
      dplyr::arrange(desc(y))
    dat$label.plot <- FALSE
    if(top.labels > 0)
      dat$label.plot[1:top.labels] <- TRUE
    if(top.labels == 0 & !is.null(selected.labels))
      dat$label.plot <- dat$label %in% selected.labels
  }

  # Arrange by annotation, in a way that some group of points are brought to the top of the point clouds
  if(!is.null(annot))
    dat <- dat %>%
    dplyr::arrange(annot)

  # Check for errors
  if(!is.null(fill.by)){
    if(!fill.by %in% names(dat))
      stop("Error: Argument fill.by not found in the data")
  }

  # Set X limits
  if(symmetric) {
    use_xlim <- xlim(-max(abs(dat$x), na.rm = TRUE), max(abs(dat$x), na.rm = TRUE))
  } else {
    use_xlim <- xlim(min(dat$x, na.rm = TRUE), max(dat$x, na.rm = TRUE))
  }

  # Volcano plot
  if(!is.null(fill.by) & !interactive){
    res <- ggplot(dat, aes(.data[['x']], .data[['y']], fill = .data[[fill.by]])) +
      gg_geom_point(size = point.size, stroke = point.stroke, shape = point.shape,
                 colour = point.colour, alpha = point.alpha) +
      use_xlim
  } else if(interactive & !is.null(label)) {
    res <- ggplot(dat, aes(.data[['x']], .data[['y']], fill = .data[[fill.by]], label = 'label')) +
      gg_geom_point(size = point.size, stroke = point.stroke, shape = point.shape,
                 colour = point.colour, alpha = point.alpha) +
      use_xlim
  }else{
    res <- ggplot(dat, aes(.data[['x']], .data[['y']])) +
      gg_geom_point(size = point.size, stroke = point.stroke, shape = point.shape,
                 colour = point.colour, alpha = point.alpha, fill = fill.palette[1])+
      use_xlim
  }

  res <- res +
    geom_hline(yintercept = -log10(pval.thr), linetype = hline.type, linewidth = hline.size, color = hline.color) +
    labs(
      fill = fill.title,
      x = x.lab,
      y = y.lab
    )
  if(!is.null(fill.by) & !is.null(fill.palette))
    res <- res + scale_fill_manual(values=fill.palette)
  if(!show.legend)
    res <- res + theme(legend.position="none")

  # Add labels to points
  if(!is.null(label) & top.labels > 0 & !interactive){
    res <- res +
      geom_text_repel(
        data = subset(dat, is.significant & label.plot),
        aes(label = label),
        size = label.size,
        min.segment.length = 0,
        segment.size  = 0.2,
        box.padding = 0.5, max.overlaps = Inf
      )

  }
  if(!is.null(label) & top.labels == 0 & !is.null(selected.labels) & !interactive){
    res <- res +
      geom_text_repel(
        data = subset(dat, label.plot),
        aes(label = label),
        size = label.size,
        min.segment.length = 0,
        segment.size  = 0.2,
        box.padding = 0.5, max.overlaps = Inf
      )

  }
  return(res)
}



#' Function to generate a MA plot
#' @description Taken from StackOverflow.
#'
#' @param x
#' @param y
#' @param rasterize Rasterize plot using ggrastr
#'
#'
#' @author Francesc Castro-Giner, \email{fcastroginer@proton.me}
#'
#' @export
#'
#' @examples
#'
plot.ma <- function(x, y, pval = NULL, annot = NULL, label = NULL,
                    y.lab = expression(paste("lo", g[2],"(Fold change)")), x.lab = 'Average of log normalized CPM',
                    pval.thr = 0.05, top.labels = 10, selected.labels = NULL, point.size = 2, point.alpha = 0.3, point.stroke = 0.01,
                    point.shape = 21, point.colour = 'white', show.legend = FALSE, label.size = 2,
                    fill.palette = c("gray60", "red"), hline.type = 2, hline.size = 1,
                    fill.by = "is.significant", fill.title = NULL, symmetric = TRUE, rasterize = FALSE){
  require(dplyr)
  require(ggplot2)
  require(ggrepel)

  # Select ggplot2 point function
  if(rasterize){
    require(ggrastr)
    gg_geom_point <- ggrastr::geom_point_rast
  } else {
    gg_geom_point <- ggplot2::geom_point
  }

  # Prepare dataframe
  dat <- cbind(x = x, y = y) %>%
    data.frame %>%
    dplyr::mutate(
      is.significant = FALSE
    )

  # Add pvalue
  if(!is.null(pval))
    dat <- cbind(dat, pval = pval)  %>%
    dplyr::mutate(
      is.significant = pval <= pval.thr
    )

  # Add labels
  if(!is.null(label))
    dat <- cbind(dat, label = label)

  # Add annotations
  if(!is.null(annot))
    dat <- cbind(dat, annot = annot)

  # Arrange data
  dat <- dat %>% dplyr::arrange(desc(abs(y)))
  if(!is.null(pval))
    dat <- dat %>% dplyr::arrange(pval)


  # Top-labels to plot after arrangement
  if(!is.null(label)) {
    dat$label.plot <- FALSE
    if(top.labels > 0)
      dat$label.plot[1:top.labels] <- TRUE
    if(top.labels == 0 & !is.null(selected.labels))
      dat$label.plot <- dat$label %in% selected.labels
  }

  # Arrange by annotation, in a way that some group of points are brought to the top of the point clouds
  if(!is.null(annot))
    dat <- dat %>%
    dplyr::arrange(annot)

  # Check for errors
  if(!is.null(fill.by)){
    if(!fill.by %in% names(dat))
      stop("Error: Argument fill.by not found in the data")
  }

  # Set Y limits
  if(symmetric) {
    use_ylim <- ylim(-max(abs(dat$y), na.rm = TRUE), max(abs(dat$y), na.rm = TRUE))
  } else {
    use_ylim <- ylim(min(dat$y, na.rm = TRUE), max(dat$y, na.rm = TRUE))
  }


  # MA plot
  if(!is.null(fill.by)){
    res <- ggplot(dat, aes(.data[['x']], .data[['y']], fill = .data[[fill.by]])) +
      gg_geom_point(size = point.size, stroke = point.stroke, shape = point.shape,
                 colour = point.colour, alpha = point.alpha) +
      use_ylim
  } else {
    res <- ggplot(dat, aes(.data[['x']], .data[['y']])) +
      gg_geom_point(size = point.size, stroke = point.stroke, shape = point.shape,
                 colour = point.colour, alpha = point.alpha, fill = fill.palette[1]) +
      use_ylim
  }

  res <- res +
    geom_hline(yintercept = 0, linetype = hline.type, linewidth = hline.size) +
    labs(
      fill = fill.title,
      x = x.lab,
      y = y.lab
    )
  if(!is.null(fill.palette))
    res <- res + scale_fill_manual(values=fill.palette)
  if(!show.legend)
    res <- res + theme(legend.position="none")

  # Add labels to points
  if(!is.null(label) & top.labels > 0){
    res <- res +
      geom_text_repel(
        data = subset(dat, label.plot),
        aes(label = label),
        size = label.size,
        min.segment.length = 0,
        segment.size  = 0.2,
        box.padding = 0.5, max.overlaps = Inf
      )

  }
  if(!is.null(label) & top.labels == 0 & !is.null(selected.labels)){
    res <- res +
      geom_text_repel(
        data = subset(dat, label.plot),
        aes(label = label),
        size = label.size,
        min.segment.length = 0,
        segment.size  = 0.2,
        box.padding = 0.5, max.overlaps = Inf
      )

  }

  return(res)
}



p.density.plot <- function(pval, annot = NULL, pval.thr = 0.05, line.type = 1, line.size = 1, line.color = 'black',
                           vline.type = 2, vline.size = 1, vline.color = 'grey50', y.lab = 'Density',
                           x.lab = '-log10(Adjusted P value)', color.by = NULL, color.palette = c("black", "red"),
                           color.title = NULL, show.legend = FALSE){
  # Prepare results data
  dat <- cbind(pval = pval) %>%
    data.frame %>%
    dplyr::mutate(
      log10pval = -log10(pval),
      is.significant = pval <= pval.thr
    )

  # Add annotations
  if(!is.null(annot))
    dat <- cbind(dat, annot = annot)

  # Check for errors
  if(!is.null(color.by)){
    if(!color.by %in% names(dat))
      stop("Error: Argument color.by not found in the data")
  }

  # Generate plot
  if(!is.null(color.by)){
    res <- ggplot(dat, aes(.data[['log10pval']], colour = color.by)) +
      geom_density(linewidth = line.size, linetype = line.type)
    if(!is.null(color.palette))
      res <- res + scale_color_manual(values=color.palette)
  } else {
    res <- ggplot(dat, aes(.data[['log10pval']])) +
      geom_density(linewidth = line.size, linetype = line.type, color = line.color)
  }
  if(!show.legend)
    res <- res + theme(legend.position="none")

  res <- res +
    geom_vline(xintercept = -log10(pval.thr), linetype = vline.type, linewidth = vline.size, color =vline.color) +
    labs(
      color = color.title,
      x = x.lab,
      y = y.lab
    )
  return(res)
}

fc.density.plot <- function(fc, annot = NULL, pval = NULL, pval.thr = 0.05, line.type = 1, line.size = 1, line.color = 'black',
                            vline.type = 2, vline.size = 1, vline.color = 'grey50', y.lab = 'Density',
                            x.lab = 'log2(Fold Change)', color.palette = c("black", "red"), show.legend = FALSE,
                            color.by = 'is.significant', color.title = NULL, symmetric = TRUE){
  # Prepare results data
  dat <- cbind(fc = fc) %>%
    data.frame

  # Add annotations
  if(!is.null(annot))
    dat <- cbind(dat, annot = annot)

  # Add pvalue
  if(!is.null(pval))
    dat <- cbind(dat, pval = pval)  %>%
      dplyr::mutate(
        is.significant = pval <= pval.thr
      )

  # Set X limits
  if(symmetric) {
    use_xlim <- xlim(-max(abs(fc)), max(abs(fc)))
  } else {
    use_xlim <- xlim(min(fc), max(fc))
  }

  # generate plot
  if(!is.null(color.by)){
    res <- ggplot(dat, aes(.data[['fc']], colour = color.by)) +
      geom_density(linewidth = line.size, linetype = line.type) +
      use_xlim
    if(!is.null(color.palette))
      res <- res + scale_color_manual(values=color.palette)
  } else {
    res <- ggplot(dat, aes(.data[['fc']])) +
      geom_density(linewidth = line.size, linetype = line.type, color = line.color) +
      use_xlim
  }

  if(!show.legend)
    res <- res + theme(legend.position="none")

  res <- res  +
    geom_vline(xintercept = 0, linetype = vline.type, linewidth = vline.size, color = vline.color) +
    labs(
      color = color.title,
      x = x.lab,
      y = y.lab
    )




  return(res)
}





