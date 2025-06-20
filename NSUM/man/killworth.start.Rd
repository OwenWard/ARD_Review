\name{killworth.start}
\alias{killworth.start}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Killworth Starting Values for MCMC
}
\description{
This function uses the Killworth estimates to calculate reasonable starting values for the MCMC estimation.
}
\usage{
killworth.start(dat, known, N)
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
The function \code{killworth.start} allows for the estimation reasonable starting values for many of the parameters in the MCMC function \code{nsum.mcmc} based on Killworth's network scale-up model. These are the default starting values where applicable. For simple subpopulation size estimation using Killworth's model, see the function \code{killworth}.
}
\value{
A list with four components:
\item{NK.start}{
a vector of positive numbers with length equal to the total number of unknown subpopulations, the starting values for the sizes of the unknown subpopulations\.
}
\item{d.start}{
a vector of positive numbers with length equal to the number of individuals, the starting values for the network degrees.
}
\item{mu.start}{
a real number, the starting value for the location parameter for the log-normal distribution of network degrees.
}
\item{sigma.start}{
a positive number, the starting value for the scale parameter for the log-normal distribution of network degrees.
}
}
\author{
Rachael Maltiel and Aaron J. Baraff

Maintainer: Aaron J. Baraff <ajbaraff at uw.edu>
}
\references{
Killworth, P., Johnsen, E., McCarty, C., Shelley, G., and Bernard, H. (1998a), "A Social Network Approach to Estimating Seroprevalence in the United States," Social Networks, 20, 23-50.

Killworth, P., McCarty, C., Bernard, H., Shelley, G., and Johnsen, E. (1998b), "Estimation of Seroprevalence, Rape, and Homelessness in the United States using a Social Network Approach," Evaluation Review, 22, 289-308.

Maltiel, R., Raftery, A. E., McCormick, T. H., and Baraff, A. J., "Estimating Population Size Using the Network Scale Up Method." CSSS Working Paper 129. Retrieved from \url{https://www.csss.washington.edu/Papers/2013/wp129.pdf}
}
\seealso{
\code{\link[NSUM:NSUM-package]{killworth.start}}, \code{\link[NSUM:NSUM-package]{nsum.mcmc}}
}
\examples{
## load data
data(McCarty)

## simulate from model with barrier effects
sim.bar <- with(McCarty, nsum.simulate(100, known, unknown, N, model="barrier",
                                       mu, sigma, rho))

## estimate Killworth starting values
dat.bar <- sim.bar$y
start <- with(McCarty, killworth.start(dat.bar, known, N))

## estimate unknown population size from MCMC
mcmc <- with(McCarty, nsum.mcmc(dat.bar, known, N, model="barrier", iterations=100,
                                burnin=50, NK.start=start$NK.start, d.start=start$d.start, 
                                mu.start=start$mu.start, sigma.start=start$sigma.start))
}