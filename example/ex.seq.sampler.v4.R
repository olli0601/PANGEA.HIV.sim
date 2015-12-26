##--------------------------------------------------------------------------------------------------------
##	simulate sequence sampling for epi simulation
##	sequences are sampled assuming an exponentially growing sequence sampling rate so that
##	1% is sampled in 1980
##	25% is sampled in 2020
##	In addition, 10% of transmissions are broken and treated as imported from outside the simulated population.
##	The infected of a broken transmission chain is considered a new index case of a transmission chain within the 
##	simulated population
##--------------------------------------------------------------------------------------------------------
\dontrun{
indir			<- '/Users/Oliver/git/HPTN071sim/source/rPANGEAHIVsim/inst/misc'
#	re-name the following:
tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp150201-scA'
dir.create(tmpdir, showWarnings=FALSE)
tmpdir.HPTN071	<- paste(tmpdir,'/TrChains',sep='')
dir.create(tmpdir.HPTN071, showWarnings=FALSE)
#	simulation input files from the epi-simulator
infile.ind		<- '150129_HPTN071_scA_IND.csv'
infile.trm		<- '150129_HPTN071_scA_TRM.csv'
#	input arguments for the pipeline
pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42, s.MODEL='Prop2DiagB4I',
			s.PREV.max=0.08, s.INTERVENTION.start=2015, s.INTERVENTION.mul= 2, s.ARCHIVAL.n=50,
			epi.model='HPTN071', epi.dt=1/48, epi.import=0.05,
			v.N0tau=1, v.r=2.851904, v.T50=-2,
			wher.mu=log(0.00447743)-0.5^2/2-0.5^2/2, wher.sigma=0.5, bwerm.mu=log(0.002239075)-0.13^2/2, bwerm.sigma=0.13, er.gamma=4,
			dbg.GTRparam=0, dbg.rER=0, index.starttime.mode='fix1955', startseq.mode='many', seqtime.mode='AtDiag')								
infile.args		<- paste(tmpdir,'/',substr(infile.ind, 1, nchar(infile.ind)-7), 'PipeArgs.R',sep='')
save(pipeline.args, file=infile.args)
#	get input into 'argv'. this is needed because the input parser is usually called from the command line, and 'argv' mimics the way input is provided when the parser is called from the command line
cmd				<- cmd.HPTN071.input.parser.v4(indir, infile.trm, infile.ind, infile.args, tmpdir.HPTN071,  infile.trm, infile.ind)				 
argv			<<- unlist(strsplit(cmd,' '))
#	run the parser
prog.HPTN071.input.parser.v4()
}
