# https://gist.github.com/MichaelJW/b4a3dd999a47b137d12f42a8f7562b11
# https://stackoverflow.com/questions/15365829/dynamic-height-and-width-for-knitr-plots
#' Function to create different heights and widths in RMarkdown Code Chunks which are Created automatically
#'
#' @description Taken from StackOverflow.
#'
#' @param g A plot.
#' @param fig.height The Figure height as in an RMarkdown document.
#' @param fig.width The Figure width as in an RMarkdown document.
#' @param dpi Dots per inch for figure.
#' @param echo As in an RMarkdown document.
#' @param message As in an RMarkdown document.
#' @param warning As in an RMarkdown document.
#'
#' @author Andreas Scharmueller, \email{andschar@@proton.me}
#'
#' @export
#'
#' @examples
#' subchunkify(plot(iris, Sepal.Length ~ Sepal.Width))
#'
#'
subchunkify = function(g,
                       fig.height = 7,
                       fig.width = 5,
                       dpi = 72,
                       echo = FALSE,
                       message = FALSE,
                       warning = FALSE,
                       suffix = NULL,
                       prefix = NULL) {
  g_deparsed = paste0(deparse(
    function() {g}
  ), collapse = '')

  if(is.null(prefix)) {
    prefix <- 'sub_chunk_'
  }
  if(is.null(suffix)) {
    suffix <- floor(runif(1) * 1e5)
  }
  sub_chunk = paste0("```{r ",prefix,"-", suffix,
                     ", fig.height=", fig.height,
                     ", fig.width=", fig.width,
                     ", dpi=", dpi,
                     ", echo=", echo,
                     ", warning=", warning,
                     ", message=", message,
                     " }",
                     "\n(",
                     g_deparsed
                     , ")()",
                     "\n```
                     ")

  cat(knitr::knit(text = knitr::knit_expand(text = sub_chunk), quiet = TRUE))
}
