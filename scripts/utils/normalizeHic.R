#' Need a function to visually normalize the Hi-C
#' data that might be sequenced to different depths.
#' 
#' @param files A named (optional) character vector
#'  specifying the path to each .hic file.
#' @inheritParams readHic
#' 
#' @importFrom plotgardener readHic
#' @returns A list of data.frames with Hi-C data that
#'  has been normalized to the local difference in
#'  diagonal strength.
#' 
normalizeHic <- function(
  files,
  chrom,
  chromstart = NULL,
  chromend = NULL,
  altchrom = NULL,
  altchromstart = NULL,
  altchromend = NULL,
  assembly = "hg38",
  resolution = "auto",
  res_scale = "BP",
  zrange = NULL,
  norm = "KR",
  matrix = "observed",
  params = NULL,
  quiet = FALSE
) {
  
  ## Read in Hi-C data
  args <- as.list(match.call()[-1])[-1]
  hicData <- lapply(files, \(f) {
    args2 <- append(args, list(file=f, zrange=NULL, quiet=TRUE))
    do.call(readHic, args2)
  })
  
  ## Get sum of diagonal values
  diagSum <- vapply(
    X=hicData, 
    FUN=\(x) {sum(x[x[,1] == x[,2], 3])},
    FUN.VALUE=numeric(1L)
  )
  
  ## Calculate scaling factors
  scaleFactors <- vapply(
    X=diagSum,
    FUN=\(x) {mean(diagSum) / x},
    FUN.VALUE=numeric(1L)
  )
  
  ## Adjust counts with scaling factor
  ans <- mapply(
    FUN=\(x,y) {x[,3] <- x[,3]*y; x},
    x=hicData,
    y=scaleFactors,
    SIMPLIFY=FALSE
  )

  ans
}

