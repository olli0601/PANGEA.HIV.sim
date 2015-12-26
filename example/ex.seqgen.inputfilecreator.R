##--------------------------------------------------------------------------------------------------------
##	1)
##	simulate sequence sampling for epi simulation
##	sequences are sampled assuming an exponentially growing sequence sampling rate so that
##	1% is sampled in 1980
##	25% is sampled in 2020
##	In addition, 10% of transmissions are broken and treated as imported from outside the simulated population.
##	The infected of a broken transmission chain is considered a new index case of a transmission chain within the 
##	simulated population
##	2)
##	Call virus tree simulator with input args as below
##	3)
##	Call SeqGen input file creator
##--------------------------------------------------------------------------------------------------------
\dontrun{
indir			<- system.file(package="rPANGEAHIVsim", "misc")
indir			<- ifelse(indir=='','/Users/Oliver/git/HPTN071sim/raw_trchain',indir)
#	re-name the following:
tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140910'
dir.create(tmpdir, showWarnings=FALSE)
tmpdir.HPTN071	<- paste(tmpdir,'/HPTN071parser',sep='')
dir.create(tmpdir.HPTN071, showWarnings=FALSE)
#	simulation input files from the epi-simulator
infile.ind		<- '140716_RUN001_IND.csv'
infile.trm		<- '140716_RUN001_TRM.csv'
#	input arguments for the pipeline
pipeline.args	<- rPANGEAHIVsim.pipeline.args()
infile.args		<- paste(tmpdir,'/',substr(infile.ind, 1, nchar(infile.ind)-7), 'PipeArgs.R',sep='')
save(pipeline.args, file=infile.args)
#	get input into 'argv'. this is needed because the input parser is usually called from the command line, and 'argv' mimics the way input is provided when the parser is called from the command line
cmd				<- cmd.HPTN071.input.parser.v2(indir, infile.trm, infile.ind, infile.args, tmpdir.HPTN071,  infile.trm, infile.ind)				 
argv			<<- unlist(strsplit(cmd,' '))
#
#	step: sequence sampler
#
prog.HPTN071.input.parser.v2()
#	
#	step: virus tree sampler
#
tmpdir.VTS		<- paste(tmpdir,'/VirusTreeSimulator',sep='')
dir.create(tmpdir.VTS, showWarnings=FALSE)
#	output file prefix for virus tree sampler
outfile			<- substr(infile.ind, 1, nchar(infile.ind)-7)
#	logistic growth within host coalescent model
prog.args		<- paste('-demoModel Logistic -N0 ',pipeline.args['v.N0tau',][,v] ,' -growthRate ', pipeline.args['v.r',][,v],' -t50 ',pipeline.args['v.T50',][,v], sep='')
cmd				<- cmd.VirusTreeSimulator(tmpdir.HPTN071, infile.trm, infile.ind, tmpdir.VTS, outfile, prog.args=prog.args)	
#	TODO run this command
cat(cmd)
#	
#	step: Seq-Gen input file creator
#
tmpdir.SG		<- paste(tmpdir,'/SeqGen',sep='')
dir.create(tmpdir.SG, showWarnings=FALSE)
infile.epi		<- paste( substr(infile.ind, 1, nchar(infile.ind)-7),'SAVE.R', sep='' )
infile.vts		<- substr(infile.ind, 1, nchar(infile.ind)-7)
cmd				<- cmd.SeqGen.createInputFiles(tmpdir.HPTN071, infile.epi, tmpdir.VTS, infile.vts, infile.args, tmpdir.SG)
argv			<<- unlist(strsplit(cmd,' '))
prog.PANGEA.SeqGen.createInputFile()
}
