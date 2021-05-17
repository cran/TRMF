\name{TRMF_es}
\alias{TRMF_es}

\title{Add exponential smoothing regularization model to a TRMF object.}
\description{Creates a regularization scheme that favors exponentially smoothed solutions and adds it to a TRMF object. In matrix optimization form, it adds the following term to the TRMF cost function:  \eqn{R(x) = lambdaD^2||w(DX_s)||^2 + lambdaA^2||X_s||^2}  where \eqn{X_s} is sub-set of the Xm matrix controlled by this model and D is a matrix with weights from exponential smoothing.}

\usage{
TRMF_es(obj,numTS = 1,alpha=1,es_type=c("single","double"),
        lambdaD=1,lambdaA=0.0001,weight=1)
}

\arguments{
  \item{obj}{A TRMF object}
  \item{numTS}{number of latent time series in this model}
  \item{lambdaD}{regularization parameter for temporal constraint matrix}
  \item{lambdaA}{regularization parameter to apply simple L2 regularization to this time series model}
  \item{weight}{optional vector of weights to weight constraint, i.e. R(x) = lambdaD^2*||w*(D\%*\%X)||^2}
  \item{es_type}{type of exponential smoothing. \code{"double"} does Brown's double exponential smoothing.}
  \item{alpha}{exponential smoothing parameter, constrained to be in the interval [0,1]}
}
\details{
This creates a non-sparse constraint matrix which could slow training down for longer time series.
}
\value{
Returns an updated object of class TRMF.
}
\references{
Yu, Hsiang-Fu, Nikhil Rao, and Inderjit S. Dhillon. "High-dimensional time series prediction with missing values." arXiv preprint arXiv:1509.08333 (2015).

https://en.wikipedia.org/wiki/Exponential_smoothing
}
\author{
Chad Hammerquist
}

\seealso{
\code{\link{create_TRMF}}, \code{\link{TRMF_columns}}, \code{\link{TRMF_trend}},\code{\link{TRMF_seasonal}}}

\examples{

# create test data
xm = cbind(cumsum(rnorm(20)),cumsum(rnorm(20)))
fm = matrix(runif(40),2,10)
Am = xm\%*\%fm+rnorm(200,0,.2)

# create model
obj = create_TRMF(Am)
obj = TRMF_es(obj,numTS=2,alpha=0.5)
out = train(obj)
plot(out)


}
