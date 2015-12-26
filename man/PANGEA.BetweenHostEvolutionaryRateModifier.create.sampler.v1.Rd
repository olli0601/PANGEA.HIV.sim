\name{PANGEA.BetweenHostEvolutionaryRateModifier.create.sampler.v1}
\alias{PANGEA.BetweenHostEvolutionaryRateModifier.create.sampler.v1}
\title{Create sampler of the Between Host Evolutionary Rate Multipliers}
\usage{
  PANGEA.BetweenHostEvolutionaryRateModifier.create.sampler.v1(bwerm.mu = 1.5,
    bwerm.sigma = 0.12)
}
\arguments{
  \item{bwerm.mu}{mean of the log normal density}

  \item{bwerm.sigma}{standard deviation of the log normal
  density}
}
\value{
  R function
}
\description{
  Returns a function to draw the between host evolutionary
  rate multipliers. Currently modelled with a log normal
  density.
}

