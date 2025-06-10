\name{McCarty}
\alias{McCarty}
\docType{data}
\title{
McCarty Dataset
}
\description{
This dataset contains the subpopulation sizes and parameters used for simulations involving the McCarty data.
}
\usage{data("McCarty")}
\format{
  A list with the following 7 variables.
  \describe{
    \item{known}{
a vector of positive numbers, the sizes of known subpopulations.
    }
    \item{unknown}{
a vector of positive numbers, the sizes of unknown subpopulations.
    }
    \item{N}{
a positive number, the (known) total population size.
    }
    \item{mu}{
a real number, the location parameter for the log-normal distribution of network degrees, with default 5.
    }
    \item{sigma}{
a positive number, the scale parameter for the log-normal distribution of network degrees, with default 1.
    }
    \item{rho}{
a vector of numbers between 0 and 1 with length equal to the total number of subpopulations, known and unknown, the dispersion parameters for the barrier effects, with defaults 0.1.
    }
    \item{tauK}{
a vector of numbers between 0 and 1 with length equal to the total number of unknown subpopulations, the multipliers for the transmission biases, with defaults 1.
    }
}
}
\details{
The McCarty data set was obtained through random digit dialing within the United
States. It contains responses from 1,375 adults from two surveys: survey 1 with 801 responses
conducted in January 1998 and survey 2 with 574 responses conducted in January
1999.
}
\source{
Killworth, P., Johnsen, E., McCarty, C., Shelley, G., and Bernard, H. (1998a), "A Social Network Approach to Estimating Seroprevalence in the United States," Social Networks, 20, 23-50.

Killworth, P., McCarty, C., Bernard, H., Shelley, G., and Johnsen, E. (1998b), "Estimation of Seroprevalence, Rape, and Homelessness in the United States using a Social Network Approach," Evaluation Review, 22, 289-308.
}
\examples{
## load data
data(McCarty)

## simulate from model with barrier effects
sim.bar <- with(McCarty, nsum.simulate(100, known, unknown, N, model="barrier", 
                                       mu, sigma, rho))
}
\keyword{datasets}
