\name{nsum.mcmc}
\alias{nsum.mcmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run MCMC for NSUM Parameters
}
\description{
This function produces an MCMC sample from the posterior distributions of the subpopulation size parameters from an NSUM model.
}
\usage{
nsum.mcmc(dat, known, N, indices.k = (length(known)+1):(dim(dat)[2]), 
          iterations = 1000, burnin = 100, size = iterations, 
          model = "degree", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
a matrix of non-negagtive integers, the \code{(i,k)}-th entry represents the number of people that the \code{i}-th individual knows from the \code{k}-th subpopulation with the columns representing known subpopulations coming before the columns representing unknown subpopulations.
}
  \item{known}{
a vector of positive numbers, the sizes of known subpopulations.
}
  \item{N}{
a positive number, the (known) total population size.
}
  \item{indices.k}{
a vector of positive integers, the indices of the columns of \code{dat} representing the unknown subpopulations of interest, with defaults of all unknown subpopulations in \code{dat}.
}
  \item{iterations}{
a positive integer, the total number of MCMC iterations after burn-in, with default 1000.
}
  \item{burnin}{
a non-negative integer, the number of burn-in MCMC iterations, with default 100.
}
  \item{size}{
a positive integer, the number of MCMC iterations kept after thinning, with default equal to \code{iterations}. 
}
  \item{model}{
a character string, the model to be simulated from. This must be one of \code{"degree"}, \code{"barrier"}, \code{"transmission"}, or \code{"combined"}, with default \code{"degree"}.
}
  \item{\dots}{
additional arguments to be passed to methods, such as starting values, prior parameters, and tuning parameters. Many methods will accept the following arguments:
    \describe{
    \item{\code{NK.start}}{
a vector of positive numbers with length equal to the total number of unknown subpopulations, the starting values for the sizes of the unknown subpopulations, with defaults based on the Killworth estimates.
    }
    \item{\code{d.start}}{
a vector of positive numbers with length equal to the number of individuals, the starting values for the network degrees, with defaults based on the Killworth estimates.
    }
    \item{\code{mu.start}}{
a real number, the starting value for the location parameter for the log-normal distribution of network degrees, with default based on the Killworth estimates.
    }
    \item{\code{sigma.start}}{
a positive number, the starting value for the scale parameter for the log-normal distribution of network degrees, with default based on the Killworth estimates.
    }
    \item{\code{rho.start}}{
a vector of numbers between 0 and 1 with length equal to the total number of subpopulations, known and unknown, the starting values for the dispersion parameters for the barrier effects, with defaults 0.1.
    }
    \item{\code{tauK.start}}{
a vector of numbers between 0 and 1 with length equal to the total number of unknown subpopulations, the starting values for the multipliers for the transmission biases, with defaults 0.5.
    }
    \item{\code{q.start}}{
a matrix of numbers between 0 and 1, the \code{(i,k)}-th entry is the starting value for the binomial probability of the number of people that the \code{i}-th individual knows from the \code{k}-th subpopulation, with defaults of simple proportions based on the known subpopulation sizes and the Killworth estimates for unknown population sizes.
    }
    \item{\code{mu.prior}}{
a vector of two real numbers, the parameters of the uniform prior for the location parameter of the log-normal distribution of network degrees, with default \code{c(3,8)}.
    }
    \item{\code{sigma.prior}}{
a vector of two positive numbers, the parameters of the uniform prior for the scale parameter of the log-normal distribution of network degrees, with default \code{c(1/4,2)}.
    }
    \item{\code{rho.prior}}{
a vector of two numbers between 0 and 1, the parameters of the uniform prior for the dispersion parameters for the barrier effects, with default \code{c(0,1)}.
    }
    \item{\code{tauK.prior}}{
a matrix of numbers between 0 and 1 with two columns and rows equal to the total number of unknown subpopulations, the parameters of the beta priors for the multipliers for the transmission biases, with defaults \code{c(1,1)}. 
    }
    \item{\code{NK.tuning}}{
a vector of positive numbers with length equal to the total number of unknown subpopulations, the standard deviations of the normal MCMC transitions for the sizes of the unknown subpopulations, with defaults of 0.25 times the starting values.
    }
    \item{\code{d.tuning}}{
a vector of positive numbers with length equal to the number of individuals, the standard deviation of the normal MCMC transitions for the network degrees, with defaults of 0.25 times the starting values.
    }
    \item{\code{rho.tuning}}{
a vector of numbers between 0 and 1 with length equal to the total number of subpopulations, known and unknown, the standard deviations of the normal MCMC transitions for the dispersion parameters for the barrier effects, with defaults of 0.25 times the starting values.
    }
    \item{\code{tauK.tuning}}{
a vector of numbers between 0 and 1 with length equal to the total number of unknown subpopulations, the standard deviations of the normal MCMC transitions for the multipliers for the transmission biases, with defaults of 0.25 times the starting values.
    }
    \item{\code{q.tuning}}{
a matrix of numbers between 0 and 1, the \code{(i,k)}-th entry is the standard deviation of the normal MCMC transitions for the binomial probability of the number of people that the \code{i}-th individual knows from the \code{k}-th subpopulation, with defaults of 0.25 times the starting values.
    }
  }
}
}
\details{
The function \code{nsum.mcmc} allows for the estimation of the various parameters from a random degree model based upon the Network Scale Up Method (NSUM) by producing Markov chain Monte Carlo (MCMC) samples from their posterior distributions. Options allow for the inclusion of barrier and transmission effects, both separately and combined, resulting in four models altogether. A large number of iterations may be required for accurate inference due to slow mixing, so the resulting chain can be thinned using the \code{size} argument. It should be noted that subpopulation size estimation in the presence of transmission bias can be greatly improved when the priors for the multipliers \code{tauK} are highly informative.
}
\value{
A list with up to nine components:
\item{NK.values}{
a matrix of positive numbers with a row for each unknown subpopulation, the thinned MCMC chains representing the posterior distributions of the sizes of the unknown subpopulations.
}
\item{d.values}{
a matrix of positive numbers with a row for each individual, the thinned MCMC chains representing the posterior distributions of the network degrees.
}
\item{mu.values}{
a vector of real numbers, the thinned MCMC chain representing the posterior distribution of the location parameter of the log-normal distribution of network degrees.
}
\item{sigma.values}{
a vector of positive numbers, the thinned MCMC chain representing the posterior distribution of the scale parameter of the log-normal distribution of network degrees.
}
\item{rho.values}{
a matrix of numbers between 0 and 1 with a row for each subpopulation, known and unknown, the thinned MCMC chains representing the posterior distributions of the dispersion parameters for the barrier effects.
}
\item{tauK.values}{
a matrix of numbers between 0 and 1 with a row for each unknown subpopulation, the thinned MCMC chains representing the posterior distributions of the multipliers for the transmission biases.
}
\item{q.values}{
a three-dimensional array of numbers between 0 and 1 with a row for each pairing of individual and subpopulation, the thinned MCMC chains representing the binomial probabilities of the number of people that the individual knows from the subpopulation.
}
\item{NK.values}{
a matrix of positive numbers with a row for each unknown subpopulation, the thinned MCMC chains representing the posterior distributions of the sizes of the unknown subpopulations.
}
\item{\code{iterations}}{
a positive integer, the total number of MCMC iterations after burn-in.
}
\item{\code{burnin}}{
a non-negative integer, the number of burn-in MCMC iterations.
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
\code{\link[NSUM:NSUM-package]{killworth.start}}
}
\examples{
## load data
data(McCarty)

## simulate from model with barrier effects
sim.bar <- with(McCarty, nsum.simulate(100, known, unknown, N, model="barrier", 
                                       mu, sigma, rho))

## estimate unknown population size
dat.bar <- sim.bar$y
mcmc <- with(McCarty, nsum.mcmc(dat.bar, known, N, model="barrier", iterations=100, 
                                burnin=50))

## view posterior distribution of subpopulation sizes for the first subpopulation
hist(mcmc$NK.values[1,])

## view posterior distribution of barrier effect parameters for the first subpopulation
hist(mcmc$rho.values[1,])
}