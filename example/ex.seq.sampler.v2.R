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
indir			<- system.file(package="rPANGEAHIVsim", "misc")
#	re-name the following:
tmpdir			<- '/Users/Oliver/git/HPTN071sim/HscABase140920'
dir.create(tmpdir, showWarnings=FALSE)
tmpdir.HPTN071	<- paste(tmpdir,'/TrChains',sep='')
dir.create(tmpdir.HPTN071, showWarnings=FALSE)
#	simulation input files from the epi-simulator
infile.ind		<- '180914_HPTN071_scA_rep1_IND.csv'
infile.trm		<- '180914_HPTN071_scA_rep1_TRM.csv'
#	input arguments for the pipeline
pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
												s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.1, 
												epi.model='HPTN071', epi.dt=1/48, epi.import=0,
												v.N0tau=1, v.r=2.851904, v.T50=-2,
												wher.mu=log(0.00447743)-0.3^2/2, wher.sigma=0.3,
												bwerm.mu=log(0.002239075)-0.13^2/2, bwerm.sigma=0.13)	
infile.args		<- paste(tmpdir,'/',substr(infile.ind, 1, nchar(infile.ind)-7), 'PipeArgs.R',sep='')
save(pipeline.args, file=infile.args)
#	get input into 'argv'. this is needed because the input parser is usually called from the command line, and 'argv' mimics the way input is provided when the parser is called from the command line
cmd				<- cmd.HPTN071.input.parser.v2(indir, infile.trm, infile.ind, infile.args, tmpdir.HPTN071,  infile.trm, infile.ind)				 
argv			<<- unlist(strsplit(cmd,' '))
#	run the parser
prog.HPTN071.input.parser.v2()
}
