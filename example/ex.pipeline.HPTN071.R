##--------------------------------------------------------------------------------------------------------
##	example pipeline to simulate sequences for a given HPTN071 epi simulation  
##--------------------------------------------------------------------------------------------------------
\dontrun{
#	re-name the following:
outdir			<- '/Users/Oliver/duke/2015_various'
#	input arguments for the pipeline
pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1985, yr.end=2020, seed=42, s.MODEL='Fixed2Prop', report.prop.recent=1.0,
			s.PREV.max.n=1600, s.INTERVENTION.prop=0.5, s.INTERVENTION.start=2015, s.INTERVENTION.mul=NA, s.ARCHIVAL.n=50,
			epi.model='HPTN071', epi.acute='high', epi.intervention='fast', epi.dt=1/48, epi.import=0.05, root.edge.fixed=0,
			v.N0tau=1, v.r=2.851904, v.T50=-2,
			wher.mu=log(0.00447743)-0.5^2/2, wher.sigma=0.5, bwerm.mu=log(0.002239075)-0.3^2/2, bwerm.sigma=0.3, er.gamma=4,
			dbg.GTRparam=0, dbg.rER=0, index.starttime.mode='fix1970', startseq.mode='one', seqtime.mode='AtDiag')									
#	
#	call simulation pipeline
#	this generates a UNIX batch file if no HPC system is detected, or
#	this generates and runs a qsub file if an HPC system is detected 
#
file			<- rPANGEAHIVsim.pipeline(outdir, pipeline.args=pipeline.args)
cat(file)
}
