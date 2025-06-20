\name{Curitiba}
\alias{Curitiba}
\docType{data}
\title{
Curitiba Dataset
}
\description{
This dataset contains the subpopulation sizes and parameters used for simulations involving the Curitiba data.
}
\usage{data("Curitiba")}
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
The Curitiba dataset consists of 500 adult residents of Curitiba, Brazil and was collected
through a household-based random sample in 2010.
}
\source{
Salganik, M., Fazito, D., Bertoni, N., Abdo, A., Mello, M., and Bastos, F. (2011a), "Assessing Network Scale-up Estimates for Groups Most at Risk of HIV/AIDS: Evidence From a Multiple-Method Study of Heavy Drug Users in Curitiba, Brazil," American Journal of Epidemiology, 174, 1190-1196.
}
\examples{
## load data
data(Curitiba)

## simulate from model with transmission bias
sim.trans <- with(Curitiba, nsum.simulate(100, known, unknown, N, model="transmission", 
                                          mu, sigma, tauK))
}
\keyword{datasets}
