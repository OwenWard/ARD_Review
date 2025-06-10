\name{NSUM-package}
\alias{NSUM-package}
\alias{NSUM}
\docType{package}
\title{
NSUM: Network Scale Up Method
}
\description{
A Bayesian framework for subpopulation size estimation using the Network Scale Up Method (NSUM). Size estimates are based on a random degree model and include options to adjust for barrier and transmission effects.
}
\details{
\tabular{ll}{
Package: \tab NSUM\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-12-17\cr
License: \tab GPL-2 | GPL-3\cr
}
The main estimation function is \code{\link[NSUM:NSUM-package]{nsum.mcmc}}. It produces a Markov chain Monte Carlo (MCMC) sample from the posterior distributions of the subpopulation size parameters from a random degree model based upon the Network Scale Up Method (NSUM). Options allow for the inclusion of barrier and transmission effects, both separately and combined, resulting in four models altogether. Also included are functions to simulate data from any of these four models (\code{\link[NSUM:NSUM-package]{nsum.simulate}}) and to estimate reasonable starting values for the MCMC sampler (\code{\link[NSUM:NSUM-package]{killworth.start}}). Two data sets have been provided for testing purposes (\code{\link[NSUM:NSUM-package]{McCarty}} and \code{\link[NSUM:NSUM-package]{Curitiba}}).
}
\author{
Rachael Maltiel and Aaron J. Baraff

Maintainer: Aaron J. Baraff <ajbaraff at uw.edu>
}
\references{
Killworth, P., Johnsen, E., McCarty, C., Shelley, G., and Bernard, H. (1998a), "A Social Network Approach to Estimating Seroprevalence in the United States," Social Networks, 20, 23-50.

Killworth, P., McCarty, C., Bernard, H., Shelley, G., and Johnsen, E. (1998b), "Estimation of Seroprevalence, Rape, and Homelessness in the United States using a Social Network Approach," Evaluation Review, 22, 289-308.

Maltiel, R., Raftery, A. E., McCormick, T. H., and Baraff, A. J., "Estimating Population Size Using the Network Scale Up Method." CSSS Working Paper 129. Retrieved from \url{https://www.csss.washington.edu/Papers/2013/wp129.pdf}

McCarty, C., Killworth, P. D., Bernard, H. R., Johnsen, E. C., and Shelley, G. A. (2001), "Comparing Two Methods for Estimating Network Size," Human Organization, 60, 28-39.

Salganik, M., Fazito, D., Bertoni, N., Abdo, A., Mello, M., and Bastos, F. (2011a), "Assessing Network Scale-up Estimates for Groups Most at Risk of HIV/AIDS: Evidence From a Multiple-Method Study of Heavy Drug Users in Curitiba, Brazil," American Journal of Epidemiology, 174, 1190-1196.
}
\keyword{ package }
\seealso{
\code{\link[NSUM:NSUM-package]{killworth.start}}, \code{\link[NSUM:NSUM-package]{nsum.mcmc}}, \code{\link[NSUM:NSUM-package]{nsum.simulate}}
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

## view posterior distribution
hist(mcmc$NK.values[1,])
}