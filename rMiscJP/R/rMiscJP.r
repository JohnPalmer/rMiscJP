#' convert a matrix to a named vector in which the names are taken from the matrix row and column names
matrix2vector = function(matrix) {unlist(lapply(1:ncol(x), function(i) lapply(1:nrow(x), function(j) {
  result = x[j,i]
  names(result) = paste(rownames(x)[j], colnames(x)[i], sep="_")
  return(result)
})))
}

#' create a timestamp of the current time in format that can be used in filename 
timestampNow <- function() paste(unlist(strsplit(unlist(strsplit(as.character(Sys.time()), " ")), ":")), collapse="-")
