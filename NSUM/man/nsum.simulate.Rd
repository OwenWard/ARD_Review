\name{nsum.simulate}
\alias{nsum.simulate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate NSUM Data
}
\description{
This function simulates data from one of the four NSUM models.
}
\usage{
nsum.simulate(n, known, unknown, N, model = "degree", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
a non-negative integer, the number respondents in the sample.
}
  \item{known}{
a vector of positive numbers, the sizes of known subpopulations.
}
  \item{unknown}{
a vector of positive numbers, the sizes of unknown subpopulations.
}
  \item{N}{
a positive number, the (known) total population size.
}
  \item{model}{
a character string, the model to be simulated from. This must be one of \code{"degree"}, \code{"barrier"}, \code{"transmission"}, or \code{"combined"}, with default \code{"degree"}.
}
  \item{\dots}{
additional arguments to be passed to methods, such as starting values, prior parameters, and tuning parameters. Many methods will accept the following arguments:
    \describe{
    \item{\code{mu}}{
a real number, the location parameter for the log-normal distribution of network degrees, with default 5.
    }
    \item{\code{sigma}}{
a positive number, the scale parameter for the log-normal distribution of network degrees, with default 1.
    }
    \item{\code{rho}}{
a vector of numbers between 0 and 1 with length equal to the total number of subpopulations, known and unknown, the dispersion parameters for the barrier effects, with defaults 0.1.
    }
    \item{\code{tauK}}{
a vector of numbers between 0 and 1 with length equal to the total number of unknown subpopulations, the multipliers for the transmission biases, with defaults 1.
    }
  }
}
}
\details{
The function \code{nsum.simulate} allows for the simulation of data from a random degree model based upon the Network Scale Up Method (NSUM). Options allow for the inclusion of barrier and transmission effects, both separately and combined, resulting in four models altogether. Each call to the function results in the simulation of a single realization of data.
}
\value{
A list with two components:
\item{y}{
a matrix of non-negagtive integers, the \code{(i,k)}-th entry represents the number of people that the \code{i}-th individual knows from the \code{k}-th subpopulation with the columns representing known subpopulations coming before the columns representing unknown subpopulations.
}
\item{d}{
a vector of positive numbers, the network degrees of the individuals. Only the integer parts were used for simulation.
}
}
\author{
Rachael Maltiel and Aaron J. Baraff

Maintainer: Aaron J. Baraff <ajbaraff at uw.edu>
}
\references{
Maltiel, R., Raftery, A. E., McCormick, T. H., and Baraff, A. J., "Estimating Population Size Using the Network Scale Up Method." CSSS Working Paper 129. Retrieved from \url{https://www.csss.washington.edu/Papers/2013/wp129.pdf}
}
\seealso{
\code{\link[NSUM:NSUM-package]{nsum.mcmc}}
}
\examples{
## load data
data(McCarty)

## simulate from model with barrier effects
sim.bar <- with(McCarty, nsum.simulate(100, known, unknown, N, model="barrier", 
                                       mu, sigma, rho))

## simulate from model with both barrier effects and transmission biases
sim.comb <- with(McCarty, nsum.simulate(100, known, unknown, N, model="combined",
                                        mu, sigma, rho, tauK))

## extract data for use in MCMC
dat.bar <- sim.bar$y

## view degree distribution
hist(sim.bar$d)
}