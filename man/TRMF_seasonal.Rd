\name{TRMF_seasonal}
\alias{TRMF_seasonal}


\title{Add seasonal regularization model to a TRMF object}
\description{Creates a regularization scheme that favors seasonally varying solutions and adds it to a TRMF object. In matrix optimization form, it adds the following term to the TRMF cost function:  \eqn{R(x) = lambdaD^2||w(DX_s)||^2 + lambdaA^2||X_s||^2}  where \eqn{X_s} is sub-set of the Xm matrix controlled by this model and D is a (with a lag of freq) finite difference matrix.}

\usage{
TRMF_seasonal(obj,numTS = 1,freq = 12,sumFirst=FALSE,lambdaD=1,lambdaA=0.0001,weight=1)
}
\arguments{
  \item{obj}{A TRMF object}
  \item{numTS}{number of latent time series in this model}
  \item{lambdaD}{regularization parameter for temporal constraint matrix}
  \item{lambdaA}{regularization parameter to apply simple L2 regularization to this time series model}
  \item{weight}{optional vector of weights to weight constraints, i.e. R(x) = lambdaD^2*||w*(D\%*\%X)||^2}
  \item{freq}{The frequency of the seasonal time series model. Minimize the differences of lag = \code{freq}}
  \item{sumFirst}{minimize the sum of first \code{freq} elements in time series}
}
\details{
 \code{TRMF_seasonal(freq=N)} fits a lag N random walk. For monthly data, use freq=12, for quarterly data, freq=4. If \code{sumFirst = TRUE}, the sum of the first \code{freq} elements in the latent time series are also minimized.This can be used to help force the seasonal component to vary around a zero mean.
}
\value{
Returns an updated object of class TRMF.
}
\references{
Yu, Hsiang-Fu, Nikhil Rao, and Inderjit S. Dhillon. "High-dimensional time series prediction with missing values." arXiv preprint arXiv:1509.08333 (2015).

}
\author{
Chad Hammerquist
}

\seealso{
\code{\link{create_TRMF}}, \code{\link{TRMF_coefficients}},\code{\link{TRMF_simple}}, \code{\link{TRMF_trend}}}

\examples{
# create test data
tm = 3*poly(x = (-20:20)/10,degree=3)
sm = diffinv(rnorm(29,0,.1),lag=12,xi=(-5:6)/6)
xm = cbind(sm,tm)
fm = matrix(runif(40),4,10)
Am = xm\%*\%fm+rnorm(410,0,.1)

# create model
obj = create_TRMF(Am)
obj = TRMF_coefficients(obj,reg_type ="interval")
obj = TRMF_trend(obj,numTS=3,order=2)
obj = TRMF_seasonal(obj,numTS=1,freq=12,lambdaD=5)
out = train(obj)
plot(out)

}
