#' Helper function to multiplot signal tracks
#' @param fn character vector of bigWig file names
#' @param p a pgParams object
#' @param y initial y-position (numeric)
#' @param labs character vector of labels
#' @param cols character vector of length(fn) or
#'  single string for the color of labels & tracks.
#' @returns a multi-plotted signal track normalized
#'  to max of y-range.
#'
#' @nord
customMultiPlot <- \(fn, p, y, labs, cols) {
    library(plotgardener)
    ypos <- pageLayoutRow(
        y = y,
        height = 0.25,
        space = 0.1,
        n = length(fn)
    )
    sigData <- lapply(fn, \(x) {
        readBigwig(file = x, params = p)
    })
    sigRange <- c(0, ceiling(max(unlist(
        lapply(sigData, "[[", "score")
    ))))

    purrr::pmap(
        .l = list(sigData, ypos, labs),
        .f = \(d, y, l) {
            plotSignal(
                params = p,
                data = d,
                range = sigRange,
                y = y,
                height = 0.25,
                linecolor = cols,
                fill = cols
            )
            plotText(
                label = paste0(l, " [0-", sigRange[2], "]"),
                fontsize = 6,
                fontcolor = cols,
                x = p$x,
                y = y,
                just = c("left", "top")
            )
        }
    )
}
