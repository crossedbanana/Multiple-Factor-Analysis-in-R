% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mfa.R
\name{CVD}
\alias{CVD}
\title{Contributions of a Variable to a Dimension}
\usage{
CVD(object)
}
\arguments{
\item{object}{an object of class \code{"mfa"}}
}
\value{
an object of class \code{"list"}
}
\description{
Outputs values that help interpret how the variables contribute to a dimension.
}
\examples{

simulate_data <- matrix(rnorm(200), ncol = 10, byrow = TRUE)
test <- mfa_const(data = simulate_data, sets = list(1:3, 4:5, 6:10), ncomps = 2)
CVD(test)

}

