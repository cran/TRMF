\name{TRMF_columns}
\alias{TRMF_columns}
\alias{TRMF_coefficients}
%- also NEED an '\alias' for EACH other topic documented here.
\title{
Add a column regularization model to TRMF object
}
\description{
Adds a regularization model to TRMF object created by \code{create_TRMF()} to constrain the fitting process of the coefficient matrix. 

\code{TRMF_coefficient} is a (soon to be deprecated) alias for \code{TRMF_columns}.
}
\usage{
TRMF_columns(obj,
    reg_type = c("l2", "nnls", "constrain", "interval", "none"), lambda = 0.0001)
TRMF_coefficients(obj,
    reg_type = c("l2", "nnls", "constrain", "interval", "none"), lambda = 0.0001)
}

\arguments{
  \item{obj}{TRMF object created by \code{create_TRMF()}}
  \item{reg_type}{regularization type to apply when fitting TRMF model. \code{l2} regularizes by simple sum of squares, \code{nnls}
  forces coefficients to be non-negative. \code{constrain} constrains coefficients to be non-negative and to sum to 1. \code{interval} constrains coefficients
  to the interval [0-1]
}
  \item{lambda}{L2 regularization parameter used for all regularization types. If NULL, uses lambda set in create_TRMF().
}
}
\details{
This function doesn't do any computations, it just sets up regularization parameters for the coefficient matrix. This function should only be called once on a TRMF object. If called twice, it will overwrite previous model with a warning.
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
\code{\link{train.TRMF}}, \code{\link{create_TRMF}}, \code{\link{TRMF_trend}}
}
\examples{

# create test data
xm = poly(x = (-10:10)/10,degree=4)
fm = matrix(abs(rnorm(40)),4,10)
Am = xm\%*\%fm+rnorm(210,0,.2)

# create model
obj = create_TRMF(Am)
obj = TRMF_columns(obj,reg_type ="nnls")
out = train(obj)
plot(out)


}
