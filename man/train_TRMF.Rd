%\name{train_TRMF}
\name{train.TRMF}
\alias{train.TRMF}

\title{Train a TRMF model}

\description{
This function is the "engine"" of the TRMF package. It takes a previously created TRMF object and fits it to the data using an alternating least squares algorithm.
}

\usage{\method{train}{TRMF}(x, numit = 10, ...)}

\arguments{
  \item{x}{A TRMF object to be fit.}
  \item{numit}{Number of alternating least squares iterations}
  \item{...}{ignored}
}
\details{If a coefficient model is not present in \code{object}, it adds a L2 regularization model. If no time series models have been added to \code{object}, it adds a simple model using \code{\link{TRMF_simple}}.
}

\value{
\code{train} returns a fitted object of \code{\link{class}} "\code{TRMF}" that contains the data, all added models, matrix factorization and fitted model. The matrix factors Xm, Fm 
are stored in \code{object$Factors$Xm} and \code{object$Factors$Fm} respectively. Use \code{\link{fitted}} to get fitted model, use \code{\link{resid}} to get residuals, use \code{\link{coef}} to get coefficients (Fm matrix) and \code{\link{components}} to get \code{Xm} or \code{Fm}.
} 
\references{
Yu, Hsiang-Fu, Nikhil Rao, and Inderjit S. Dhillon. "High-dimensional time series prediction with missing values." arXiv preprint arXiv:1509.08333 (2015).
}

\author{
Chad Hammerquist
}



\seealso{
\code{\link{create_TRMF}}, \code{\link{TRMF_coefficients}}, \code{\link{TRMF_trend}}}

\examples{
# create test data
xm = poly(x = (-10:10)/10,degree=4)
fm = matrix(rnorm(40),4,10)
Am = xm\%*\%fm+rnorm(210,0,.2)

# create model
obj = create_TRMF(Am)
out = train(obj)
plot(out)

}

