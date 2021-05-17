\name{components.TRMF}
\alias{components.TRMF}


\title{
 Access TRMF factors
}
\description{
This function returns the factors (Xm, Fm) from a trained TRMF object
}
\usage{
\method{components}{TRMF}(object, XorF = c("Xm","Fm"), ...)
}

\arguments{
  \item{object}{trained TRMF object}
  \item{XorF}{ which factor to return}
  \item{...}{ignored}
}
\details{
Returns the matrix factors. Could also use \code{object$Factors$Xm, object$Factors$Fm} . If matrix normalization was used in \code{\link{create_TRMF}}, \code{Xm\%*\%Fm} could look much different than the input data matrix. 
}
\value{
A matrix.
}
\author{
Chad Hammerquist
}



\seealso{
\code{\link{create_TRMF}}, \code{\link{TRMF_columns}}, \code{\link{TRMF_trend}}}

\examples{
# create test data
xm = poly(x = (-10:10)/10,degree=4)
fm = matrix(rnorm(40),4,10)
Am = xm\%*\%fm+rnorm(210,0,.2)

# create model
obj = create_TRMF(Am)
out = train(obj)
plot(out)
components(out,"Xm")

}


