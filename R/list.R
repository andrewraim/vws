#' Access member variable from all elements of a list
#'
#' @param x A list.
#' @param name Name of variable to access from each element of the list.
#' @param unlist If \code{TRU}, unlist the result before returning.
#'
#' @examples
#' x = list()
#' x[[1]] = list(a = 1:10, b = 1:5)
#' x[[2]] = list(a = 1:5, b = 0, c = -1)
#' Access(x, "a")
#' Access(x, "a", unlist = FALSE)
#' Access(x, "b")
#' Access(x, "c")
#'
#' @export
Access = function(x, name, unlist = TRUE) {
	out = Map(function(z) { z[`name`] }, x)
	if (unlist) { return(unlist(out)) } else { return(out) }
}

