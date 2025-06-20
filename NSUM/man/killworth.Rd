\name{killworth}
\alias{killworth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Killworth Estimates
}
\description{
This function calculates the Killworth estimates for unknown subpopulation sizes based on NSUM data.
}
\usage{
killworth(dat, known, N)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
a matrix of non-negagtive integers, the \code{(i,k)}-th entry represents the number of people that the \code{i}-th individual knows from the \code{k}-th subpopulation.
}
  \item{known}{
a vector of positive numbers, the sizes of known subpopulations. All additional columns of \code{dat} are treated as unknown.
}
  \item{N}{
a positive number, the (known) total population size.
}
}
\details{
The function \code{killworth} allows for the estimation of subpopulation sizes from Killworth's network scale-up model. These estimates can be used to compare with the MCMC results in this package. For reasonable starting values for the MCMC function \code{nsum.mcmc}, see the function \code{killworth.start}.
}
\value{
A vector of positive numbers with length equal to the number of unknown subpopulations, the Killworth estimates of the subpopulation sizes.
}
\author{
Rachael Maltiel and Aaron J. Baraff

Maintainer: Aaron J. Baraff <ajbaraff at uw.edu>
}
\references{
Killworth, P., Johnsen, E., McCarty, C., Shelley, G., and Bernard, H. (1998a), "A Social Network Approach to Estimating Seroprevalence in the United States," Social Networks, 20, 23-50.

Killworth, P., McCarty, C., Bernard, H., Shelley, G., and Johnsen, E. (1998b), "Estimation of Seroprevalence, Rape, and Homelessness in the United States using a Social Network Approach," Evaluation Review, 22, 289-308.
}
\seealso{
\code{\link[NSUM:NSUM-package]{killworth.start}}
}
\examples{
## load data
data(McCarty)

## simulate from model with barrier effects
sim.bar <- with(McCarty, nsum.simulate(100, known, unknown, N, model="barrier",
                                       mu, sigma, rho))

## estimate unknown population sizes
dat.bar <- sim.bar$y
NK.killworth <- with(McCarty, killworth(dat.bar, known, N))
}