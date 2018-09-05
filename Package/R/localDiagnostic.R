#' Calculate the value for local diagnostic function
#'
#' This function calculates local diagnostic function D(x,h) at each point x which depends only on observations in a small neighborhood [x-h,x+h].
#' @param y the numeric vector of the intensities of markers
#' @param h the bandwidth for the screening procedure
#' @details
#' Local diagnostic function reflects of position x being or neighborhooding a change-point. A reasonable
#' local diagnostic is\deqn{D(x)=\frac{\Sigma_{k=1}^{h} y_{x+k}-\Sigma_{k=1}^{h} y_{x+1-k}}{h}}
#' which is the difference between averages of h data points on the left side and right side of x.
#' Suppose the errors \eqn{\varepsilon_{i}=0} which means \eqn{Y=\mu} is a piecewise constant vector and D(x) is piecewise linear function.
#' Based on this function, we proposed a recursive formula\deqn{D(x+1)_{h}=D(x)_{h}+\frac{Y_{x-h+1}+Y_{x+h+1}-2Y_{x+1}}{h}}
#' @return This function generates a numeric vector of local diagnostic function values D(x,h) at each point x
#' @useDynLib modSaRa2
#' @export
localDiagnostic <-function(y,h) {
    yy = c(rep(0, h-1), y, rep(0, h))            # add some zeros on the head and tail
    n = length(y)
    z = rep(0, n)
    ans<-.C("diagnosticValue", yy=as.double(yy), h=as.integer(h), n=as.integer(n), z=as.double(z))

    return (ans$z)
}
