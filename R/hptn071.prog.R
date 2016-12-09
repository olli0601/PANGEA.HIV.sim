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
#' @title Input arguments to simulate sequences and viral phylogenetic trees\cr under the Regional Transmission and Intervention Model
#' @description Set all input arguments to the simulation. Please see for model details \url{https://github.com/olli0601/PANGEA.HIV.sim}
#' @seealso \code{\link{sim.regional}}
#' @import data.table
#' @param yr.start				Start year of the epi simulation (default: 1985)
#' @param yr.end				First year after the epi simulation (default: 2020; optional 2018)
#' @param seed					Random number seed (default: 42; optional integer)
#' @param s.MODEL				Sampling model to use (default: 'Fixed2Prop'). See \url{https://github.com/olli0601/PANGEA.HIV.sim#sequence-sampling-model}.
#' @param s.INC.recent			Deprecated
#' @param s.INC.recent.len		Deprecated 
#' @param s.PREV.min			Deprecated
#' @param s.PREV.max			Deprecated
#' @param s.PREV.max.n			Number of infected cases sampled (default: 1600; optional 3200). See \url{https://github.com/olli0601/PANGEA.HIV.sim#sequence-sampling-model}.
#' @param s.INTERVENTION.prop	Proportion of infected cases that are sampled from after intervention start (default: 0.5; optional 0.85). See \url{https://github.com/olli0601/PANGEA.HIV.sim#sequence-sampling-model}.
#' @param s.INTERVENTION.start	Year in which the community intervention starts (default: 2015). See \url{https://github.com/olli0601/PANGEA.HIV.sim#sequence-sampling-model}.
#' @param s.INTERVENTION.mul	Deprecated
#' @param s.ARCHIVAL.n			Total number of sequences sampled at random from infected individuals before 2000 (default: 50). See \url{https://github.com/olli0601/PANGEA.HIV.sim#sequence-sampling-model}.
#' @param seqtime.mode			Time at which sequences are sampled (default: 'AtDiag'). See \url{https://github.com/olli0601/PANGEA.HIV.sim#sequence-sampling-model}.
#' @param epi.model				The epi model used to create the epi simulation (default: 'HPTN071')
#' @param epi.acute				The proportion of transmissions in 2015 from individuals in their first three months of infection (default: 'high'/40\%, optional: 'low'/10\%)
#' @param epi.intervention		The scale of the intervention package (default:'fast', optional: 'slow', 'none')				 
#' @param epi.dt				Time increment of the epi simulation (default: 1/48, this is fixed to the underlying epidemic simulation)
#' @param epi.import			Proportion of viral introductions amongst new cases per year (default: 0.05, optional: 0.2). See \url{https://github.com/olli0601/PANGEA.HIV.sim#viral-introductions}.
#' @param root.edge.fixed		Fix evolutionary rate of root edges of new lineages in the population to mean evolutionary rate (default: 1, optional: 0). See \url{https://github.com/olli0601/PANGEA.HIV.sim#overall-evolutionary-rates-of-transmission-and-non-transmission-lineages}. 
#' @param v.N0tau				Within host effective population size at time of infection (default 1). See \url{https://github.com/olli0601/PANGEA.HIV.sim#dated-viral-phylogenies}. 
#' @param v.r					Parameter logistic growth model (default 2.851904). See \url{https://github.com/olli0601/PANGEA.HIV.sim#dated-viral-phylogenies}.
#' @param v.T50					Parameter logistic growth model (default -2). See \url{https://github.com/olli0601/PANGEA.HIV.sim#dated-viral-phylogenies}.
#' @param wher.mu				Mean of within host evolutionary rate (default: log(0.00447743)-0.5^2/2). See \url{https://github.com/olli0601/PANGEA.HIV.sim#overall-evolutionary-rates-of-transmission-and-non-transmission-lineages}.
#' @param wher.sigma			Standard deviation of within host evolutionary (default: 0.5). See \url{https://github.com/olli0601/PANGEA.HIV.sim#overall-evolutionary-rates-of-transmission-and-non-transmission-lineages}.
#' @param bwerm.mu				Mean of between host evolutionary rate (default: log(0.002239075)-0.3^2/2). See \url{https://github.com/olli0601/PANGEA.HIV.sim#overall-evolutionary-rates-of-transmission-and-non-transmission-lineages}.
#' @param bwerm.sigma			Standard deviation of between host evolutionary rate (default: 0.3). See \url{https://github.com/olli0601/PANGEA.HIV.sim#overall-evolutionary-rates-of-transmission-and-non-transmission-lineages}.
#' @param startseq.mode			Number of different starting sequences (default: 'one', optional 'many'). See \url{https://github.com/olli0601/PANGEA.HIV.sim#seed-sequences}.
#' @param index.starttime.mode	Sampling distribution of sampling times of starting sequences (default: 'fix1970', optional 'fix19XX' where XX between 45-79, 'normal' meaning Normal(1955,7)). See \url{https://github.com/olli0601/PANGEA.HIV.sim#seed-sequences}.
#' @param er.gamma				Site heterogeneity parameter (default: 4, optional 0)
#' @param er.gtr				Evolution model identifier. Default: GTR_VARIABLE. Alternatively: GTR_POL_CP1.
#' @param sp.prop.of.sexactive	Proportion of population sampled in seroprevalence survey (default: 0.05)
#' @param report.prop.recent	Proportion of individuals in seroprevalence survey for whom recency of infection should be reported (default: 1)
#' @param dbg.GTRparam			Use mean GTR rates instead of samples from empirical prior (default: 0, optional 1)
#' @param dbg.rER				Use mean evolutionary rate instead of samples from empirical prior (default: 0, optional 1)
#' @return data.table with columns 'stat' (name of input argument) and 'v' (value of input argument).
#' @export
sim.regional.args<- function(			yr.start=1985, yr.end=2020, seed=42,										
										s.INC.recent=NA, s.INC.recent.len=NA, s.PREV.min=NA, s.PREV.max=NA, s.PREV.base=NA, s.INTERVENTION.start=2015, 
										s.INTERVENTION.mul= NA, s.ARCHIVAL.n=50, s.MODEL='Fixed2Prop', s.PREV.max.n=1600, s.INTERVENTION.prop=0.5,
										epi.model='HPTN071', epi.acute='high', epi.intervention='fast', epi.dt=1/48, epi.import=0.05, root.edge.fixed=0,
										v.N0tau=1, v.r=2.851904, v.T50=-2,
										wher.mu=log(0.00447743)-0.5^2/2, wher.sigma=0.5,
										bwerm.mu=log(0.002239075)-0.3^2/2, bwerm.sigma=0.3, er.gamma=4, er.gtr='GTR_VARIABLE', 
										sp.prop.of.sexactive= 0.05, 
										report.prop.recent=1.0,
										dbg.GTRparam=0, dbg.rER=0, 
										index.starttime.mode='fix1970', startseq.mode='one', seqtime.mode='AtDiag')									
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
		
		ggplot(subset(tmp, variable=='y2'), aes(x=-x, y=value, group=variable)) + geom_line() + scale_y_continuous(breaks=c(1e2, 2e2, 3e2, 4e2), limits=c(0,4e2)) + scale_x_continuous(breaks=seq(-20,20,1)) + labs(x='\ntime since infection\n(years)',y='effective population size\n')
		ggsave(file='~/git/PANGEA.HIV.sim/man/fig_EffPopSize.png', width=8, height=4)
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
	pipeline.args	<- data.table(	stat= 	c('yr.start','yr.end','s.MODEL','s.INC.recent','s.INC.recent.len', 's.PREV.min', 's.PREV.max', 's.PREV.max.n', 's.INTERVENTION.prop', 's.INTERVENTION.start', 's.INTERVENTION.mul', 's.ARCHIVAL.n', 's.seed', 'index.starttime.mode', 'startseq.mode', 'seqtime.mode','root.edge.fixed','epi.model', 'epi.acute', 'epi.intervention', 'epi.dt', 'epi.import','v.N0tau','v.r','v.T50','wher.mu','wher.sigma','bwerm.mu','bwerm.sigma','er.gamma','er.gtr','sp.prop.of.sexactive','report.prop.recent','dbg.GTRparam','dbg.rER'), 
									v	=	c(yr.start, yr.end, s.MODEL, s.INC.recent, s.INC.recent.len, s.PREV.min, s.PREV.max, s.PREV.max.n, s.INTERVENTION.prop, s.INTERVENTION.start, s.INTERVENTION.mul, s.ARCHIVAL.n, seed, index.starttime.mode, startseq.mode, seqtime.mode, root.edge.fixed, epi.model, epi.acute, epi.intervention, epi.dt, epi.import, v.N0tau, v.r, v.T50, wher.mu, wher.sigma, bwerm.mu, bwerm.sigma, er.gamma, er.gtr, sp.prop.of.sexactive, report.prop.recent, dbg.GTRparam, dbg.rER) )
	setkey(pipeline.args, stat)	
	pipeline.args
}
##--------------------------------------------------------------------------------------------------------
##	olli originally written 06-08-2015
##--------------------------------------------------------------------------------------------------------
pipeline.various<- function()
{
	if(0)	#submit various
	{
		cmd			<- cmd.various()
		cmd			<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q=NA, hpc.walltime=71, hpc.mem="63000mb")
		#cmd			<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeelab', hpc.walltime=71, hpc.mem="64000mb")
		cat(cmd)		
		outdir		<- paste(HOME,"tmp",sep='/')
		outfile		<- paste("vrs",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
		cmd.hpccaller(outdir, outfile, cmd)
		quit("no")		
	}	
	#	split into individual distances and run MVR on HPC because MDS needs a ton of RAM
	if(0)
	{
		require(big.phylo)
		
		wdir	<- '/work/or105/Gates_2014/tree_comparison/mvr'
		#wdir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/tree_mvr'
		infile	<- '150701_Regional_TRAIN4_SIMULATED_tps.rda'
		infile	<- '150701_Regional_TRAIN2_SIMULATED_tps.rda'
		load(file.path(wdir, infile))
		
		name.gd.col					<- 'GD_MEAN'
		complete.distance.matrix	<- 0
		#loop.method					<- c('MVR','BioNJ')
		loop.method					<- 'BioNJ'
		#loop.rep					<- tp[, unique(REP)]
		loop.rep					<- setdiff( tp[, unique(REP)], c(1) )
		loop.gene					<- setdiff(tp[, unique(GENE)],'env')
		#loop.gene					<- "gag+pol+env"
		for(gene in loop.gene)
			for(rep in loop.rep)
				for(method in loop.method)
				{
					tps		<- subset(tp, GENE==gene & REP==rep)
					tmp		<- file.path(wdir, gsub('\\.rda',paste('_GENE_',gene,'_REP_',rep,'_M_',method,'_C_',complete.distance.matrix,'.rda',sep=''), infile))
					save(tps, file=tmp)
					cmd		<- cmd.mvr(tmp, outfile=gsub('\\.rda','',tmp), method=method, complete.distance.matrix=complete.distance.matrix, name.gd.col=name.gd.col)
					cmd	<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeelab', hpc.walltime=8, hpc.mem="16000mb")
					#cmd		<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.walltime=8, hpc.mem="16000mb")
					cat(cmd)		
					cmd.hpccaller(paste(HOME,"tmp",sep='/'), paste("mvr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.'), cmd)
					#quit('no')
				}
		quit('no')
	}
	if(0)	#run LSD
	{
		require(ape)
		require(data.table)
		require(big.phylo)
		#wdir	<- '~/duke/tmp'		 
		#file	<- file.path('~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation','/','submitted_160713_05QD.rda')		
		wdir	<- '/work/or105/Gates_2014/tree_comparison/lsd'
		file	<- '/work/or105/Gates_2014/tree_comparison/submitted_161123.rda'
		load(file)				
		
		setkey(submitted.info, TEAM, SC, GENE, RUNGAPS)

		ds		<- subset(submitted.info, MODEL=='R')
		#ds		<- subset(submitted.info, MODEL=='R' & IDX%in%c(76,79,106,143,145,146,176,232,331,356,358,359,361))
		#ds		<- subset(submitted.info, TEAM=='IQTree' & SC=='150701_REGIONAL_TRAIN1' & OTHER=='N' & GENE%in%c('GAG','GAG+POL+ENV'))[, list(IDX=IDX[1]), by='GENE']
		#ds		<- rbind(ds, subset(submitted.info, TEAM=='IQTree' & SC=='150701_REGIONAL_TRAIN2' & OTHER=='N' & GENE%in%c('GAG','GAG+POL+ENV'))[, list(IDX=IDX[1:5]), by='GENE'])
		#ds		<- rbind(ds, subset(submitted.info, TEAM=='IQTree' & SC=='150701_REGIONAL_TRAIN4' & OTHER=='N' & GENE%in%c('GAG','GAG+POL+ENV'))[, list(IDX=IDX[1:5]), by='GENE'])		
		#ds		<- rbind(ds, subset(submitted.info, TEAM=='RUNGAPS_ExaML' & SC=='150701_REGIONAL_TRAIN2' & OTHER=='N' & GENE%in%c('P17','GAG','GAG+POL+ENV'))[, list(IDX=IDX[seq.int(1,length(IDX),3)]), by='GENE'])		
		#ds		<- merge(ds, submitted.info, by=c('IDX','GENE'))	
		ds[, {	
					#ph				<- strs_rtt[[IDX]]
					ph				<- strs[[IDX]]
					tmp				<- data.table(TAXA=ph$tip.label, SEQ_T= sapply(strsplit(ph$tip.label, '|', fixed=1),'[[', 4))
					tmp				<- tmp[,  list(STR=paste(TAXA,' ',SEQ_T,sep='')), by='TAXA'][, paste(STR, collapse='\n')]
					tmp				<- paste(Ntip(ph),'\n',tmp,sep='')
					#infile.tree		<- file.path(wdir, paste(gsub('\\.treefile','',basename(FILE)),'_IDX_',IDX,'_GAPS_',GAPS,'_RUNGAPS_',RUNGAPS,'.newick',sep=''))				
					#infile.dates	<- file.path(wdir, paste(gsub('\\.treefile','',basename(FILE)),'_IDX_',IDX,'_GAPS_',GAPS,'_RUNGAPS_',RUNGAPS,'.txt',sep=''))
					#outfile			<- file.path(wdir, paste(gsub('\\.treefile','',basename(FILE)),'_IDX_',IDX,'_GAPS_',GAPS,'_RUNGAPS_',RUNGAPS,'_LSD',sep=''))					
					infile.tree		<- file.path(wdir, paste('IDX_',IDX,'.newick',sep=''))				
					infile.dates	<- file.path(wdir, paste('IDX_',IDX,'.txt',sep=''))
					outfile			<- file.path(wdir, paste('IDX_',IDX,'_LSD',sep=''))
					cat(tmp, file=infile.dates)
					write.tree(ph, file=infile.tree)
					ali.nrow		<- 6800
					if(GENE=='GAG')
						ali.nrow	<- 1440	
					if(GENE=='GAG+PARTIALPOL')
						ali.nrow	<- 3000
					if(GENE=='P17')
						ali.nrow	<- 400	
					if(GENE=='POL')
						ali.nrow	<- 2850						
					cmd				<- cmd.lsd(infile.tree, infile.dates, ali.nrow, outfile=outfile, pr.args='-v 2 -c -b 10 -r as')
					#cmd				<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeelab', hpc.walltime=40, hpc.mem="11800mb")
					#cmd				<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeph', hpc.walltime=40, hpc.mem="3600mb")
					cmd				<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q=NA, hpc.walltime=3, hpc.mem="1800mb")
					cat(cmd)		
					cmd.hpccaller(wdir, paste("lsd",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.'), cmd)					
					#quit('no')
				}, by='IDX']		
		quit('no')
	}
	if(1)	#	run ExamML with partition for tree comparison
	{		
		require(big.phylo)		
		#indir.wgaps	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/running_gaps_simulations3'
		#indir.wgaps	<- '/work/or105/Gates_2014/tree_comparison/rungaps3'
		indir.wgaps	<- '/work/or105/Gates_2014/tree_comparison/rungaps4'
		#indir.wgaps	<- '/work/or105/Gates_2014/tree_comparison/rungaps'
		infiles		<- data.table(FILE=list.files(indir.wgaps, pattern='\\.fasta$|\\.fa$'))
		infiles[, PARTITION:= gsub('\\.fasta|\\.fa','_gene.txt',FILE)]
		infiles	<- subset(infiles, grepl('TRAIN64',FILE))
		outdir		<- indir.wgaps	 			
		infiles[, {	
					args.starttree.type	<- 'parsimony'
					#args.starttree.type	<- 'random'
					args.parser			<- "-m DNA "
					if(file.exists(file.path(indir.wgaps, PARTITION)))
						args.parser		<- paste("-m DNA -q",PARTITION)					
					cmd		<- cmd.examl.single(indir.wgaps, FILE, outdir=outdir, args.parser=args.parser, args.starttree.type=args.starttree.type, args.examl="-m GAMMA -f d -D", verbose=1)
					#cmd		<- cmd.hpcwrapper(cmd, hpc.walltime=571, hpc.q="pqeelab", hpc.mem="23850mb", hpc.nproc=4)
					cmd		<- cmd.hpcwrapper(cmd, hpc.walltime=71, hpc.q=NA, hpc.mem="15850mb", hpc.nproc=24)
					#cmd		<- cmd.hpcwrapper(cmd, hpc.walltime=571, hpc.q="pqeelab", hpc.mem="5850mb", hpc.nproc=1)
					signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
					outfile	<- paste("exr4",signat,sep='.')
					cat(cmd)
					cmd.hpccaller(outdir, outfile, cmd)
					Sys.sleep(1)
					#stop()
				}, by='FILE']	
		quit('no')
	}
	if(0)	#	run ExamML with partition for tree comparison in replicate of 10
	{		
		require(big.phylo)		
		indir.wgaps	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/simulations_simpleGTR'
		indir.wgaps	<- '/work/or105/Gates_2014/tree_comparison/simpleGTR'
		infiles		<- data.table(FILE=list.files(indir.wgaps, pattern='\\.fasta$'))
		infiles[, PARTITION:= gsub('\\.fasta','_gene.txt',FILE)]
		outdir		<- indir.wgaps		
		infiles		<- subset(infiles, grepl('FULL',FILE))
		infiles[, {					
					for(i in 1:10)
					{
						args.parser		<- "-m DNA "
						if(file.exists(file.path(indir.wgaps, PARTITION)))
							args.parser	<- paste("-m DNA -q",PARTITION)
						#args.examl		<- "-m GAMMA -f d -D"
						args.examl		<- "-m GAMMA"
						cmd				<- cmd.examl.single(indir.wgaps, FILE, outdir=outdir, outfile=gsub('\\.fasta',paste('_REP',i,sep=''),FILE), args.parser=args.parser, args.examl=args.examl, verbose=1)
						cmd				<- cmd.hpcwrapper(cmd, hpc.walltime=410, hpc.q="pqeelab", hpc.mem="5900mb", hpc.nproc=1)
						signat			<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
						outfile			<- paste("ex",signat,sep='.')
						cat(cmd)
						hivc.cmd.hpccaller(outdir, outfile, cmd)
						Sys.sleep(1)
						#stop()	
					}					
				}, by='FILE']					
	}
	if(0)	#	run ExamML with partition for prev unsuccessful runs
	{		
		require(big.phylo)		
		outdir		<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/simulations'
		outdir		<- '/work/or105/Gates_2014/tree_comparison/simpleGTR'
		indir.wgaps	<- outdir
		infiles		<- data.table(OF=c("161121_GTRFIXED2_FULL_SIMULATED_REP5", "161121_GTRFIXED3_FULL_SIMULATED_REP1", "161121_GTRFIXED3_FULL_SIMULATED_REP2", "161121_GTRFIXED3_FULL_SIMULATED_REP5"))
		infiles[, IF:=paste(gsub('_REP[0-9]+','',OF),'.fasta',sep='')]
		infiles[, PARTITION:= gsub('\\.fasta','_gene.txt',IF)]
		infiles[, {					
					args.parser		<- "-m DNA "
					if(file.exists(file.path(indir.wgaps, PARTITION)))
						args.parser	<- paste("-m DNA -q",PARTITION) 
					cmd				<- cmd.examl.single(indir.wgaps, IF, outdir=outdir, outfile=OF, args.parser=args.parser, args.examl="-m GAMMA -f d -D", verbose=1)
					cmd				<- cmd.hpcwrapper(cmd, hpc.walltime=241, hpc.q="pqeelab", hpc.mem="5900mb", hpc.nproc=1)
					signat			<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
					outfile			<- paste("ex",signat,sep='.')
					cat(cmd)
					hivc.cmd.hpccaller(outdir, outfile, cmd)
					Sys.sleep(1)
					#stop()						
				}, by='OF']
	}
	if(0)	#	run ExamML with partition for prev unsuccessful runs
	{		
		require(big.phylo)		
		outdir		<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/simulations'
		outdir		<- '/work/or105/Gates_2014/tree_comparison/rungaps2'
		indir.wgaps	<- outdir
		infiles		<- data.table(OF=c("150701_Regional_TRAIN21750_FULL_SIMULATED","150701_Regional_TRAIN21780_FULL_SIMULATED",
						"150701_Regional_TRAIN22350_FULL_SIMULATED","150701_Regional_TRAIN23070_FULL_SIMULATED",
						"150701_Regional_TRAIN23090_FULL_SIMULATED","150701_Regional_TRAIN23150_FULL_SIMULATED",
						"150701_Regional_TRAIN23750_FULL_SIMULATED","150701_Regional_TRAIN24860_FULL_SIMULATED"))
		infiles[, IF:=paste(gsub('_REP[0-9]+','',OF),'.fa',sep='')]
		infiles[, PARTITION:= gsub('\\.fa','_gene.txt',IF)]
		infiles[, {					
					args.parser		<- "-m DNA "
					if(file.exists(file.path(indir.wgaps, PARTITION)))
						args.parser	<- paste("-m DNA -q",PARTITION) 
					cmd				<- cmd.examl.single(indir.wgaps, IF, outdir=outdir, outfile=OF, args.parser=args.parser, args.examl="-m GAMMA -f d -D", verbose=1)
					cmd				<- cmd.hpcwrapper(cmd, hpc.walltime=41, hpc.q="pqeelab", hpc.mem="5900mb", hpc.nproc=1)
					signat			<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
					outfile			<- paste("ex",signat,sep='.')
					cat(cmd)
					hivc.cmd.hpccaller(outdir, outfile, cmd)
					Sys.sleep(1)
					#stop()						
				}, by='IF']
	}
	if(0)	#calculate genetic distances in alignment + get bootstrap variance
	{
		#batch.i		<- 1
		indir		<- file.path(DATA, 'gds')	
		infile.fa	<- '150701_Regional_TRAIN2_SIMULATED.fa'
		infile.ge	<- '150701_Regional_TRAIN2_SIMULATED_gene.txt'
		outdir		<- indir
		for(batch.i in 1:400)
		{
			outfile	<- paste(gsub('.fa','',infile.fa),'_GDS_BATCH',batch.i,'.rda',sep='')
			cmd		<- cmd.gendist(indir, infile.fa, infile.ge, outdir, outfile, batch.i)		
			cmd		<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeelab', hpc.walltime=58, hpc.mem="5000mb")
			cat(cmd)		
			cmd.hpccaller(paste(HOME,"tmp",sep='/'), paste("vrs",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.'), cmd)				
		}			
		quit("no")
	}	
	if(0)	#submit ExaML runs
	{
		require(ape)
		require(big.phylo)	
		
		indir		<- '~/git/HPTN071sim/treecomparison/withgaps_160729'		
		indir		<- file.path(DATA, 'gpsi')
		infile		<- data.table(FILE=list.files(indir, pattern='\\.R'))
		infile[, SC:= gsub('\\..*','',FILE)]
		tmp			<- data.table(PARTITION=list.files(indir, pattern='\\.txt'))
		tmp[, SC:= gsub('_gene.txt','',PARTITION)]
		infile		<- merge(infile, tmp, by='SC', all.x=1)
		#	
		invisible(infile[, {
					#SC<- '150701_Regional_TRAIN202_FULL_SIMULATED'; PARTITION<- '150701_Regional_TRAIN202_FULL_SIMULATED_gene.txt'
					args.parser			<- paste("-m DNA")
					args.examl			<- "-f d -D -m GAMMA"
					args.starttree.type	<- 'parsimony'
					if(!is.na(PARTITION))						
						args.parser			<- paste("-m DNA -q",PARTITION)
					if(!grepl('FULL', SC))
						args.starttree.type	<- 'random'
					
					cmd		<- cmd.examl.single(indir, SC, outdir=indir, args.starttree.type=args.starttree.type, args.parser=args.parser, args.examl=args.examl, bs.seed=42)
					cmd		<- cmd.hpcwrapper(cmd, hpc.walltime=129, hpc.q="pqeelab", hpc.mem="5800mb", hpc.nproc=1)
					signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
					cat(cmd)
					cmd.hpccaller(indir, paste("exa",signat,sep='.'), cmd)
					Sys.sleep(1)
				}, by='SC'])		
	}
	if(0)	#submit ExaML of partial sequences 
	{
		require(ape)
		require(big.phylo)	
		
		infile		<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/simulations/150701_Regional_TRAIN1_SIMULATED.fa'
		infile		<- '/work/or105/Gates_2014/tree_comparison/partiallen/150701_Regional_TRAIN1_SIMULATED.fa'
		indir		<- dirname(infile)
		infile		<- as.data.table(expand.grid(PARTIAL_LEN=round(seq(0.05, .99, 0.01)*6807), IF=infile, stringsAsFactors=FALSE))
		infile[, OF:= paste(gsub('TRAIN1.*', '',gsub('150701','161125',basename(IF))),paste('TRAIN1_PL',PARTIAL_LEN,'_SIMULATED',sep=''), sep='')]
		infile		<- subset(infile, PARTIAL_LEN<=2900)
		invisible(infile[, {																					
							seq					<- read.dna(IF, format='fa')
							ifp					<- paste(OF, '.fasta', sep='')
							write.dna(seq[,1:PARTIAL_LEN], file=file.path(dirname(IF), ifp), format='fa', , colsep='', nbcol=-1)
													
							args.examl			<- "-f d -D -m GAMMA"
							#args.starttree.type	<- 'parsimony'
							args.starttree.type	<- 'random'
							
							partition			<- paste(OF,'_gene.txt',sep='')
							if(PARTIAL_LEN>=4500)
								cat(paste('DNA, gag = 1-1440\nDNA, pol = 1441-4284\nDNA, env = 4285-',PARTIAL_LEN,'\n',sep=''), file=file.path(indir, partition))
							if(PARTIAL_LEN>=1700 & PARTIAL_LEN<4500)
								cat(paste('DNA, gag = 1-1440\nDNA, pol = 1441-',PARTIAL_LEN,'\n',sep=''), file=file.path(indir, partition))
							if(PARTIAL_LEN<1700)
								partition		<- NA
							
							args.parser			<- paste("-m DNA")
							if(!is.na(partition))						
								args.parser		<- paste("-m DNA -q",partition)
							
							cmd		<- cmd.examl.single(	dirname(IF), 
															ifp, 
															outdir=dirname(IF), 
															outfile= OF, 
															args.starttree.type=args.starttree.type, 
															args.parser=args.parser, 
															args.examl=args.examl, 
															bs.seed=42)
							cmd		<- cmd.hpcwrapper(cmd, hpc.walltime=129, hpc.q="pqeelab", hpc.mem="5800mb", hpc.nproc=1)
							signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
							cat(cmd)
							cmd.hpccaller(indir, paste("exa",signat,sep='.'), cmd)
							Sys.sleep(1)
						}, by='OF'])				
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
prog.treecomparison<- function()
{
	if(1)
	{
		#file	<- '/work/or105/Gates_2014/tree_comparison/submitted_151101.rda'
		#file	<- '/work/or105/Gates_2014/tree_comparison/submitted_160627.rda'
		#file	<- '/work/or105/Gates_2014/tree_comparison/submitted_160713.rda'
		file	<- '/work/or105/Gates_2014/tree_comparison/submitted_161123.rda'
		#treedist.quartets.add(file=file, with.save=1)
		treecomparison.submissions.160627.stuffoncluster(file)
		#treedist.billera.add(file=file, with.save=1)		
	}
	if(0)
	{
		indir	<- wdir	<- file.path(DATA,'gds')
		treecomparison.bootstrap.sd.vs.coverage(indir, wdir)
	}	
	quit("no")
}
##--------------------------------------------------------------------------------------------------------
##	olli originally written 08-09-2014
##--------------------------------------------------------------------------------------------------------
#' @title Simulate sequences and viral phylogenetic trees\cr under the Regional Transmission and Intervention Model
#' @description The Regional Transmission and Intervention Model captures individual-level HIV transmission dynamics in a regional population of ~80,000 invididuals,  
#' that is broadly similar to a site (cluster) of the HPTN071/PopART HIV prevention trial in South Africa (Hayes et al., 2014).
#' \cr
#' \cr 
#' As of December 2015, this individual-level model is unpublished. 
#' The model builds on work of the HIV modelling consortium (Eaton et al., 2012) as well as
#' an earlier compartmental model that was used to inform the design of the HPTN071/PopART trial (Cori et al., 2014). 
#' An extended version of this individual-level model will be used to help evaluate results from the HPTN071/PopART trial.
#' \cr
#' \cr 
#' Please see for model details \url{https://github.com/olli0601/PANGEA.HIV.sim}
#' @import data.table 
#' @param outdir			Output directory. Must have write access. Directory name must not contain whitespace, brackets, etc
#' @param pipeline.args		Input arguments for the simulation, in form of a \code{data.table}, see \code{\link{sim.regional.args}}.
#' @return File name of qsub or UNIX batch file.  
#' @seealso \code{\link{sim.regional.args}}
#' @example example/ex.pipeline.HPTN071.R
#' @export
#' @references  Hayes R, Ayles H, Beyers N, Sabapathy K, Floyd S, et al. (2014) HPTN 071 (PopART): rationale and design of a cluster-randomised trial of the population impact of an HIV combination prevention intervention including universal testing and treatment - a study protocol for a cluster randomised trial. Trials 15: 57.
#' @references	Eaton JW, Johnson LF, Salomon JA, Barnighausen T, Bendavid E, et al. (2012) HIV Treatment as Prevention: Systematic Comparison of Mathematical Models of the Potential Impact of Antiretroviral Therapy on HIV Incidence in South Africa. PLoS Med 9: e1001245.
#' @references	Cori A, Ayles H, Beyers N, Schaap A, Floyd S, et al. (2014) HPTN 071 (PopART): a cluster-randomized trial of the population impact of an HIV combination prevention intervention including universal testing and treatment: mathematical model. PLoS One 9: e84511.
sim.regional<- function(outdir, pipeline.args=sim.regional.args() )
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
# start: run sim.regional
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
# end: run sim.regional
#
#######################################################
#######################################################
#######################################################\n",sep='')
	if(verbose)
		cat(cmd)
	outfile			<- paste("pngea",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')	
	cmd.hpccaller(outdir, outfile, cmd)	
}