##--------------------------------------------------------------------------------------------------------
##	example pipeline to simulate sequences for a given DSPS epi simulation  
##--------------------------------------------------------------------------------------------------------
\dontrun{
indir			<- system.file(package="rPANGEAHIVsim", "misc")
#	re-name the following:
tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140912'
dir.create(tmpdir, showWarnings=FALSE)
#	simulation input files from the epi-simulator
infile.trm		<- '140911_DSPS_RUN002_TRM.csv'
infile.ind		<- NULL
#	input arguments for the pipeline
pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42, s.PREV.min=0.01, s.PREV.max=0.70, epi.model='DSPS', epi.dt=1/48, epi.import=0.1 )	
infile.args		<- paste(tmpdir,'/',substr(infile.trm, 1, nchar(infile.trm)-7), 'PipeArgs.R',sep='')
save(pipeline.args, file=infile.args)
#	
#	call simulation pipeline
#	this generates a UNIX batch file if no HPC system is detected, or
#	this generates and runs a qsub file if an HPC system is detected 
#
file			<- rPANGEAHIVsim.pipeline(indir, infile.ind, infile.trm, tmpdir, pipeline.args=pipeline.args)
cat(file)
}
