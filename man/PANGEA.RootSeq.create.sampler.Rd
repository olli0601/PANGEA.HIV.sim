\name{PANGEA.RootSeq.create.sampler}
\alias{PANGEA.RootSeq.create.sampler}
\title{Create starting sequence sampler}
\usage{
  PANGEA.RootSeq.create.sampler(root.ctime.grace = 0.5,
    sample.grace = 3)
}
\arguments{
  \item{root.ctime.grace}{Sample a starting sequence with
  time that matches a query times +- this grace}

  \item{sample.grace}{Internal parameter to make sure the
  requested number of samples is obtained. Internally
  oversample by this multiplier to the sample size, and
  then check if sequences are unique.}
}
\value{
  list of the sampler \code{rANCSEQ} and its arguments
  \code{rANCSEQ.args}
}
\description{
  Returns a function and function arguments to draw
  ancestral sequences.
}

