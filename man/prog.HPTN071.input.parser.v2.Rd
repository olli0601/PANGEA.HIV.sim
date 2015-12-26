\name{prog.HPTN071.input.parser.v2}
\alias{prog.HPTN071.input.parser.v2}
\title{HPTN071 parser (version 2, includes simulation of imports)}
\usage{
  prog.HPTN071.input.parser.v2()
}
\value{
  NULL. Saves simulations to file.
}
\description{
  Reads files from the epi simulator in directory
  \code{indir} and writes csv files in directory
  \code{outdir} for the virus tree simulator. The program
  samples sequences according to an exponentially
  increasing sampling fraction in the same way as
  \code{prog.HPTN071.input.parser.v1}. In addition,
  transmissions are broken and treated as imported from
  outside the simulated population. The infected of a
  broken transmission chain is considered a new index case
  of a transmission chain within the simulated population.
  All input arguments are specified via the \code{argv}
  string, see the Examples.
}
\examples{
##--------------------------------------------------------------------------------------------------------
##	simulate sequence sampling for epi simulation
##	sequences are sampled assuming an exponentially growing sequence sampling rate so that
##	1\% is sampled in 1980
##	25\% is sampled in 2020
##	In addition, 10\% of transmissions are broken and treated as imported from outside the simulated population.
##	The infected of a broken transmission chain is considered a new index case of a transmission chain within the 
##	simulated population
##--------------------------------------------------------------------------------------------------------
\\dontrun{
indir			<- system.file(package="rPANGEAHIVsim", "misc")
#	re-name the following:
tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140914'
dir.create(tmpdir, showWarnings=FALSE)
tmpdir.HPTN071	<- paste(tmpdir,'/TrChains',sep='')
dir.create(tmpdir.HPTN071, showWarnings=FALSE)
#	simulation input files from the epi-simulator
infile.ind		<- '140716_RUN001_IND.csv'
infile.trm		<- '140716_RUN001_TRM.csv'
#	input arguments for the pipeline
pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42, s.PREV.min=0.01, s.PREV.max=0.20, epi.dt=1/48, epi.import=0.1)
infile.args		<- paste(tmpdir,'/',substr(infile.ind, 1, nchar(infile.ind)-7), 'PipeArgs.R',sep='')
save(pipeline.args, file=infile.args)
#	get input into 'argv'. this is needed because the input parser is usually called from the command line, and 'argv' mimics the way input is provided when the parser is called from the command line
cmd				<- cmd.HPTN071.input.parser.v2(indir, infile.trm, infile.ind, infile.args, tmpdir.HPTN071,  infile.trm, infile.ind)				 
argv			<<- unlist(strsplit(cmd,' '))
#	run the parser
prog.HPTN071.input.parser.v2()
}
}

