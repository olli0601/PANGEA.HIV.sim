######################################################################################
#	Program to add gaps into sequences  	
#	olli originally written 23-06-2015
######################################################################################
prog.PANGEA.AddGaps.run.v1<- function()
{	
	indir.simu		<- '/Users/Oliver/git/HPTN071sim/treec150623/nogaps'
	indir.gap		<-	'~/git/HPTN071sim/treec150623/PANGEAreal'
	infile.simu		<- '150227_HPTN071_TRAIN1_SIMULATED'
	infile.gap		<- '150623_GlobalAlignment_cov1.fasta'
	outdir			<- '/Users/Oliver/git/HPTN071sim/treec150623/withgaps'
	outfile.cov		<- regmatches(infile.gap,regexpr('cov[0-9]+',basename(infile.gap)))	
	gap.country		<- 'ZA'
	gap.symbol		<- '?'
	gap.seed		<- 42
	outfile			<- paste(infile.simu, '_', gap.country, outfile.cov, '.fa', sep='')
	verbose			<- 1
	#
	#	read args
	#
	if(exists("argv"))
	{
		#	args input
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									indir.simu= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir.simu<- tmp[1]	
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									indir.gap= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir.gap<- tmp[1]	
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.gap= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.gap<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									infile.simu= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.simu<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,8),
									outfile= return(substr(arg,10,nchar(arg))),NA)	}))
		if(length(tmp)>0) outfile<- tmp[1]
		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									gap.country= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) gap.country<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									gap.symbol= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) gap.symbol<- tmp[1]					
	}
	if(verbose)
	{
		cat('\ninput args\n',paste(indir.simu, indir.gap, infile.simu, infile.gap, outdir, outfile.cov, outfile, gap.country, gap.symbol, gap.seed	, sep='\n'))
	}
	sgp		<- PANGEA.add.gaps(indir.simu, indir.gap, infile.simu, infile.gap, gap.country, gap.symbol, gap.seed, verbose=1)
	write.dna(sgp, file=paste(outdir, outfile, sep='/'), format='fasta', colsep='', nbcol=-1)	
}

##--------------------------------------------------------------------------------------------------------
##	olli originally written 10-09-2014
##--------------------------------------------------------------------------------------------------------
#' @title Arguments for rPANGEAHIV simulation pipeline
#' @description Construct table with all input arguments to the rPANGEAHIV simulation pipeline.
#' The demographic within host coealescent model is N0tau*(1+exp(-r*T50))/(1+exp(-r*(T50-t))). Default parameters are set
#' so that the curve asymptotes at Ne*tau=3e5 and reaches half its asymptotic value 1 year post infection. Branch lengths are 
#' multiplied by within host evolutionary rates, and within host branch lengths ending in a transmission event are multiplied
#' in addition with a between-host evolutionary rate multiplier. 
#' @param yr.start				Start year of the epi simulation (default 1980)
#' @param yr.end				First year after the epi simulation (default 2020)
#' @param seed					Random number seed
#' @param s.MODEL				Sampling model to use
#' @param s.INC.recent			Proportion of incident cases sampled (not used)
#' @param s.INC.recent.len		Number of last year in which an exact proportion of incident cases is sampled (not used) 
#' @param s.PREV.min			Proportion of infected cases sampled at start of the simulation
#' @param s.PREV.max			Proportion of infected cases sampled at the end of the simulation
#' @param s.PREV.max.n			Number of infected cases sampled (usually NA, only used for Prop2Untreated)
#' @param s.INTERVENTION.prop	Proportion of infected cases that are sampled from after intervention start
#' @param s.INTERVENTION.start	Year in which the community intervention starts
#' @param s.INTERVENTION.mul	Multiplier of number of sequences sampled per year after start of the intervention
#' @param s.ARCHIVAL.n			Total number of sequences sampled before diagnosis
#' @param epi.model				The epi model used to create the epi simulation (default HPTN071, alternatively DSPS)
#' @param epi.dt				Time increment of the epi simulation (default 1/48)
#' @param epi.import			Proportion of imported cases of all cases (default 0.1)
#' @param root.edge.fixed		Boolean; fix evolutionary rate of root edges to mean rate if true 
#' @param v.N0tau				Parameter of BEAST::LogisticGrowthN0::N0 (default 3.58e4) 
#' @param v.r					Parameter of BEAST::LogisticGrowthN0::r (default 2)
#' @param v.T50					Parameter of BEAST::LogisticGrowthN0::T50
#' @param wher.mu				Mean within host evolutionary rate of log normal density
#' @param wher.sigma			Standard deviation in within host evolutionary rate of log normal density
#' @param bwerm.mu				Mean between host evolutionary rate multiplier of log normal density
#' @param bwerm.sigma			Standard deviation in between host evolutionary rate multiplier of log normal density
#' @param sp.prop.of.sexactive	Proportion of population sampled in seroprevalence survey
#' @param report.prop.recent	Proportion of individuals for whom recency of infection should be reported
#' @param dbg.GTRparam			debug flag
#' @param dbg.rER				debug flag
#' @param startseq.mode			Number of different starting sequences either 'many' or 'one'.
#' @param index.starttime.mode	distribution to sample times for starting sequence: Normal(1960,7) or Unif(1959.75, 1960.25)
#' @return data.table
#' @example example/ex.pipeline.HPTN071.R
#' @export
rPANGEAHIVsim.pipeline.args<- function(	yr.start=1980, yr.end=2020, seed=42,										
										s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.25, s.PREV.base=exp(1), s.INTERVENTION.start=2015, 
										s.INTERVENTION.mul= 2, s.ARCHIVAL.n=50, s.MODEL='Prop2DiagB4I', s.PREV.max.n=NA, s.INTERVENTION.prop=NA,
										epi.model='HPTN071', epi.acute='high', epi.intervention='fast', epi.dt=1/48, epi.import=0.1, root.edge.fixed=0,
										v.N0tau=3.58e4, v.r=2, v.T50=-1,
										wher.mu=0.005, wher.sigma=0.8,
										bwerm.mu=1.5, bwerm.sigma=0.12, er.gamma=4, 
										sp.prop.of.sexactive= 0.05, 
										report.prop.recent=0.2,
										dbg.GTRparam=0, dbg.rER=0, 
										index.starttime.mode='normal', startseq.mode='many', seqtime.mode='Gamma9')									
{
	#explore within host Neff*tau model to choose pars
	if(0)
	{
		#visualize dependence of size after 4 years on growth rate
		r		<- seq(1,16,0.1)
		plot(r, ( 1+exp(r) ) / ( 1+exp(-3*r) ), type='l')
		#estimate growth rates for desired sizes	T50=-1
		f		<- function(r, b){	abs(( 1+exp(r) ) / ( 1+exp(-3*r) )-b)	}		
		optimize(f=f, interval=c(2, 15), b=3e5 )	#12.61152
		optimize(f=f, interval=c(2, 15), b=3e4 )	#10.30891
		optimize(f=f, interval=c(2, 15), b=3e3 )	#8.006034
		optimize(f=f, interval=c(2, 15), b=3e2 )	#5.70044
		#estimate growth rates for desired sizes	T50=-2
		f		<- function(r, b){	abs(( 1+exp(2*r) ) / ( 1+exp(-2*r) )-b)	}		
		optimize(f=f, interval=c(2, 15), b=3e5 )	#6.305779
		optimize(f=f, interval=c(2, 15), b=3e4 )	#5.154461
		optimize(f=f, interval=c(2, 15), b=3e3 )	#4.003191
		optimize(f=f, interval=c(2, 15), b=3e2 )	#2.851904		
		f		<- function(r, b){	abs(( 1+exp(2*r) ) / ( 1+exp(-r*(-2+10)) )-b)	}
		optimize(f=f, interval=c(2, 15), b=3e2 )	#2.850242
		optimize(f=f, interval=c(2, 15), b=5e4 )	#5.409859
		optimize(f=f, interval=c(2, 15), b=2e5 )	#6.103021
		optimize(f=f, interval=c(2, 15), b=150 )	#2.501955
		optimize(f=f, interval=c(0.1, 15), b=10 )	#1.098685
		
		Net			<- function(t, N0tau, r, T50){  N0tau*(1+exp(-r*T50))/(1+exp(-r*(T50-t)))	}		
		x			<- seq(-10,0,0.001)
		#tmp			<- data.table(x=x, y5=Net(x, 1, 12.61152, -1), y4=Net(x, 1, 10.30891, -1), y3=Net(x, 1, 8.006034, -1), y2=Net(x, 1, 5.70044, -1))
		tmp			<- data.table(x=x, y5=Net(x, 1, 6.305779, -2), y4=Net(x, 1, 5.154461, -2), y3=Net(x, 1, 4.003191, -2), y2=Net(x, 1, 2.850242, -2), y1=Net(x, 1, 1.098685, -2))
		tmp			<- melt(tmp, id.var='x')
		ggplot(tmp, aes(x=x, y=value, group=variable, colour=variable)) + geom_line() + scale_y_log10(breaks=c(3e2,3e3,3e4,3e5)) + scale_x_continuous(breaks=seq(-20,10,1))
	}
	#	plot increasing sampling fraction to choose pars
	if(0)
	{
		s.PREV.base	<- exp(1)
		yr.start	<- 1980
		yr.end		<- 2020
		s.PREV.min	<- 0.01
		df.sample	<- as.data.table(expand.grid(YR=seq.int(yr.start, yr.end), s.PREV.max=c(0.11, 0.15, 0.185)))		
		#	exponential rate of increasing s.TOTAL (total sampling rate) per year
		set(df.sample, NULL, 'r', df.sample[, log( s.PREV.max/s.PREV.min, base=s.PREV.base ) / diff(range(YR)) ] )
		set(df.sample, NULL, 's.fraction', df.sample[, s.PREV.base^( r*(YR-min(YR)) ) * s.PREV.min ] )
		set(df.sample, NULL, 's.PREV.max', df.sample[, factor(s.PREV.max, levels=c('0.11','0.15','0.185'), labels=c('A','B','C'))]) 			
		ggplot(df.sample, aes(x=YR, y=s.fraction, colour=s.PREV.max)) + geom_line() + scale_y_continuous(breaks=seq(0, 0.16, 0.02)) + labs(colour='scenario', x='year', y='sampling fraction')		
	}
	stopifnot(epi.model=='HPTN071', epi.acute%in%c('low','high'), epi.intervention%in%c('none','slow','fast'), epi.dt==1/48)
	pipeline.args	<- data.table(	stat= 	c('yr.start','yr.end','s.MODEL','s.INC.recent','s.INC.recent.len', 's.PREV.min', 's.PREV.max', 's.PREV.max.n', 's.INTERVENTION.prop', 's.INTERVENTION.start', 's.INTERVENTION.mul', 's.ARCHIVAL.n', 's.seed', 'index.starttime.mode', 'startseq.mode', 'seqtime.mode','root.edge.fixed','epi.model', 'epi.acute', 'epi.intervention', 'epi.dt', 'epi.import','v.N0tau','v.r','v.T50','wher.mu','wher.sigma','bwerm.mu','bwerm.sigma','er.gamma','sp.prop.of.sexactive','report.prop.recent','dbg.GTRparam','dbg.rER'), 
									v	=	c(yr.start, yr.end, s.MODEL, s.INC.recent, s.INC.recent.len, s.PREV.min, s.PREV.max, s.PREV.max.n, s.INTERVENTION.prop, s.INTERVENTION.start, s.INTERVENTION.mul, s.ARCHIVAL.n, seed, index.starttime.mode, startseq.mode, seqtime.mode, root.edge.fixed, epi.model, epi.acute, epi.intervention, epi.dt, epi.import, v.N0tau, v.r, v.T50, wher.mu, wher.sigma, bwerm.mu, bwerm.sigma, er.gamma, sp.prop.of.sexactive, report.prop.recent, dbg.GTRparam, dbg.rER) )
	setkey(pipeline.args, stat)	
	pipeline.args
}
##--------------------------------------------------------------------------------------------------------
##	olli originally written 06-08-2015
##--------------------------------------------------------------------------------------------------------
pipeline.various<- function()
{
	if(1)	#align sequences in fasta file with Clustalo
	{
		cmd			<- cmd.various()
		cmd			<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeelab', hpc.walltime=771, hpc.mem="5000mb")
		cat(cmd)		
		outdir		<- paste(HOME,"tmp",sep='/')
		outfile		<- paste("vrs",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
		cmd.hpccaller(outdir, outfile, cmd)
		quit("no")		
	}	
	if(0)
	{
		mfile		<- paste(DATA,'model_150816a.R',sep='/')
		indir.st	<- paste(DATA,'contigs_150408_wref_cutstat',sep='/')
		indir.al	<- paste(DATA,'contigs_150408_wref',sep='/')
		outdir		<- paste(DATA,'contigs_150408_model150816a',sep='/')
		trainfile	<- paste(DATA,'contigs_150408_trainingset_subsets.R',sep='/')
		batch.n		<- 200
		for(batch.id in seq.int(1,14))
		{			
			cmd			<- cmd.haircut.call(indir.st, indir.al, outdir, mfile, trainfile=trainfile, batch.n=batch.n, batch.id=batch.id, prog=PR.HAIRCUT.CALL )
			cmd			<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeelab', hpc.walltime=21, hpc.mem="5000mb")
			cat(cmd)		
			outdir		<- paste(HOME,"tmp",sep='/')
			outfile		<- paste("cntcall",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			cmd.hpccaller(outdir, outfile, cmd)				
		}	
		quit("no")
	}
}
##--------------------------------------------------------------------------------------------------------
##	olli originally written 07-11-2015
##--------------------------------------------------------------------------------------------------------
prog.treecomparison.metrics<- function()
{
	file	<- '/work/or105/Gates_2014/tree_comparison/submitted_151101.rda'
	treedist.quartets.add(file=file, with.save=1)
	#treedist.billera.add(file=file, with.save=1)
	quit("no")
}
##--------------------------------------------------------------------------------------------------------
##	olli originally written 08-09-2014
##--------------------------------------------------------------------------------------------------------
#' @title rPANGEAHIV simulation pipeline
#' @description Reads two files \code{infile.ind} and \code{infile.trm} from the epi simulator in directory \code{indir} and produces a UNIX batch
#' file that contains all the simulation steps in directory \code{outdir}. 
#' @param indir				Input directory
#' @param infile.ind		Input file with individual metavariables
#' @param infile.trm		Input file with transmission events
#' @param pipeline.args		Input arguments for the simulation
#' @return file name of qsub or UNIX batch file
#' @example example/ex.pipeline.HPTN071.R
#' @example example/ex.pipeline.DSPS.R
#' @export
rPANGEAHIVsim.pipeline<- function(outdir, pipeline.args=rPANGEAHIVsim.pipeline.args() )
{
	verbose			<- 1
	#
	#	pipeline start
	#	
	##	sample sequences and draw imports 
	cmd				<- "#######################################################
#######################################################
#######################################################
#
# start: run rPANGEAHIVsim.pipeline
#
#######################################################
#######################################################
#######################################################"		
	outdir.EpiSim	<- paste(outdir,'/EpiSim',sep='')
	##	run HPTN071 IBM
	cmd				<- paste(cmd, '\nmkdir -p ', outdir.EpiSim,sep='')
	cmd				<- paste(cmd, '\ncd ', outdir.EpiSim,sep='')
	cmd				<- paste(cmd, cmd.HPTN071.simulator(outdir.EpiSim, seed=pipeline.args['s.seed',][, as.numeric(v)], opt.acute=pipeline.args['epi.acute',][, v], opt.intervention=pipeline.args['epi.intervention',][, v]), sep='\n')
	cmd				<- paste(cmd, 'cd ..', sep='')
	outdir.TrChain	<- paste(outdir,'/TrChains',sep='')
	cmd				<- paste(cmd, '\nmkdir -p ', outdir.TrChain,sep='')
	stopifnot(pipeline.args['epi.model'][,v]=='HPTN071')
	##	get infile names
	infile.ind		<- basename(substring(regmatches(cmd, regexpr('mv phylogenetic_individualdata[^\n]+', cmd)), 33))
	infile.trm		<- basename(substring(regmatches(cmd, regexpr('mv phylogenetic_transmission[^\n]+', cmd)), 31))
	##	save pipeline pars
	infile.args		<- paste(outdir,'/',substr(infile.trm, 1, nchar(infile.trm)-7), 'PipeArgs.R',sep='')
	save(pipeline.args, file=infile.args)
	##	run input parser
	cmd				<- paste(cmd, cmd.HPTN071.input.parser.v4(outdir.EpiSim, infile.trm, infile.ind, infile.args, outdir.TrChain,  infile.trm, infile.ind), sep='\n')			
	##	run virus tree simulator
	outdir.VTS		<- paste(outdir,'/VirusTreeSimulator',sep='')
	cmd				<- paste(cmd, 'mkdir -p ', outdir.VTS,sep='')
	outfile			<- substr(infile.ind, 1, nchar(infile.ind)-7)
	prog.args		<- paste('-seed ',pipeline.args['s.seed',][, v],' -demoModel Logistic -N0 ',pipeline.args['v.N0tau',][,v] ,' -growthRate ', pipeline.args['v.r',][,v],' -t50 ',pipeline.args['v.T50',][,v], sep='')	
	cmd				<- paste(cmd, cmd.VirusTreeSimulator(outdir.TrChain, infile.trm, infile.ind, outdir.VTS, outfile, prog.args=prog.args), sep='\n')	
	##	create seq gen input files 
	outdir.SG		<- paste(outdir,'/SeqGen',sep='')
	cmd				<- paste(cmd, 'mkdir -p ', outdir.SG,sep='')
	infile.epi		<- paste( substr(infile.ind, 1, nchar(infile.ind)-7),'SAVE.R', sep='' )
	infile.vts		<- substr(infile.ind, 1, nchar(infile.ind)-7)
	cmd				<- paste(cmd, cmd.SeqGen.createInputFiles(outdir.TrChain, infile.epi, outdir.VTS, infile.vts, infile.args, outdir.SG), sep='\n')
	##	run SeqGen	
	outfile			<- substr(infile.ind, 1, nchar(infile.ind)-7)
	cmd				<- paste(cmd, cmd.SeqGen.run(outdir.TrChain, infile.epi, outdir.SG, outfile, infile.args, outdir), sep='')
	##	clean up
	cmd				<- paste(cmd,'rm -rf ',outdir.EpiSim, ' ',outdir.TrChain,' ', outdir.VTS,' ', outdir.SG,'\n',sep='')
	cmd				<- paste(cmd,"#######################################################
#######################################################
#######################################################
#
# end: run rPANGEAHIVsim.pipeline
#
#######################################################
#######################################################
#######################################################\n",sep='')
	if(verbose)
		cat(cmd)
	outfile			<- paste("pngea",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')	
	cmd.hpccaller(outdir, outfile, cmd)	
}