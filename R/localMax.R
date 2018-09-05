#' Get the local maximizers of local diagnostic function
#'
#' This function finds the local maximizers of local diagnostic function |D(x,h )| and returns the value of position index for all the local maximizers.
#' @param y the list of local diagnostic values within a small neighborhood [x-h, x+h]
#' @param span the bandwidth to find local Maximizer of local diagnostic function
#' @return numeric vector of position index for all the local maximizers of local diagnostic function
#' @details Local maximizer is defined as follows. For any number x, the interval (x-h, x+h) is called the h-neighborhood of x. And, x is an h-local maximizer of function f(.) if reaches the maximum at x within the h-neighborhood at x. In other words, f(x)>=f(x') for all x' in (x-h, x+h).
#' @export
localMax <-function(y, span = 5){
    if (length(y) < span * 2 + 1)  return(NULL)
    n  = length(y)
    index = NULL
    for (i in (span + 1) : (n - span) ) {
         if ( y[i] == max(y[(i - span) : (i + span)]) ) index = c(index, i)
    }
    return (index)
}
