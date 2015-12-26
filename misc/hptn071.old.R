##--------------------------------------------------------------------------------------------------------
#	return ancestral sequence sampler	
#	olli originally written 10-09-2014
##--------------------------------------------------------------------------------------------------------
PANGEA.RootSeq.create.sampler.v2<- function(root.ctime.grace= 0.5, sample.grace= 3)
{	
	#tree.id.labelsep		<- '|'
	#tree.id.labelidx.ctime	<- 4
	file			<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfgBwhRc-_140902_n390_AncSeq.R')
	cat(paste('\nLoading starting sequences from file', file))
	load(file)		#expect "anc.seq.gag"  "anc.seq.pol"  "anc.seq.env"  "anc.seq.info"
	setkey(anc.seq.info, CALENDAR_TIME)
	rANCSEQ.args	<- list(	root.ctime.grace=root.ctime.grace, sample.grace=sample.grace, anc.seq.info=anc.seq.info, anc.seq.gag=anc.seq.gag, anc.seq.pol=anc.seq.pol, anc.seq.env=anc.seq.env)	
	
	rANCSEQ<- function(root.ctime, rANCSEQ.args)
	{		
		tmp		<- lapply(seq_along(root.ctime), function(i)
				{
					tmp	<- subset(rANCSEQ.args$anc.seq.info, CALENDAR_TIME>root.ctime[i]-rANCSEQ.args$root.ctime.grace &  CALENDAR_TIME<=root.ctime[i]+rANCSEQ.args$root.ctime.grace)
					if(nrow(tmp)<rANCSEQ.args$sample.grace*100)
					{
						warning( paste('\nFor root',i,': number of samples is n=',nrow(tmp),'. safe pool size is n=',rANCSEQ.args$sample.grace*100) )
					}
					#print( paste('\n',nrow(tmp),'\t',rANCSEQ.args$sample.grace*100) )
					data.table( LABEL= tmp[, sample( LABEL, rANCSEQ.args$sample.grace ) ], CALENDAR_TIME=root.ctime[i], DRAW=i )
				})
		tmp		<- do.call('rbind',tmp)
		#	get unique seqs
		setkey(tmp, LABEL)
		tmp				<- unique(tmp)	
		anc.seq.draw	<- do.call( 'cbind', list( rANCSEQ.args$anc.seq.gag[tmp[, LABEL], ], rANCSEQ.args$anc.seq.pol[tmp[, LABEL], ], rANCSEQ.args$anc.seq.env[tmp[, LABEL], ] ) ) 
		anc.seq.draw	<- seq.unique(anc.seq.draw)
		#	check if at least one seq for each draw
		tmp				<- merge( data.table(LABEL=rownames(anc.seq.draw)), tmp, by='LABEL' )	
		stopifnot( !length(setdiff( seq_along(root.ctime), tmp[, unique(DRAW)] )) )
		#	take first seq for each draw	
		tmp				<- tmp[, list(LABEL= LABEL[1]), by='DRAW']
		setkey(tmp, DRAW)
		anc.seq.draw	<- anc.seq.draw[ tmp[, LABEL], ]
		#	set new calendar time for sequences
		#set(tmp, NULL, 'LABEL_NEW', tmp[, as.numeric( sapply( strsplit(LABEL,tree.id.labelsep,fixed=TRUE), '[[', tree.id.labelidx.ctime) ) ])
		#set(tmp, NULL, 'LABEL_NEW', tmp[, LABEL_NEW+sample.shift])	
		#tmp			<- tmp[,	{
		#							z							<- strsplit(LABEL,tree.id.labelsep,fixed=TRUE)[[1]]
		#							z[tree.id.labelidx.ctime]	<- LABEL_NEW
		#							list(LABEL_NEW=paste(z, collapse=tree.id.labelsep,sep=''))
		#						}, by='LABEL']
		#setkey(tmp, LABEL)
		#rownames(anc.seq.draw)		<- tmp[ rownames(anc.seq.draw), ][, LABEL_NEW]		
		anc.seq.draw
	}
	list(rANCSEQ=rANCSEQ, rANCSEQ.args=rANCSEQ.args)
}
##--------------------------------------------------------------------------------------------------------
#	return ancestral sequence sampler	
#	olli originally written 22-08-2014
##--------------------------------------------------------------------------------------------------------
PANGEA.RootSeq.create.sampler.v1<- function(root.ctime.grace= 0.5, sample.grace= 3, sample.shift= 40)
{	
	#tree.id.labelsep		<- '|'
	#tree.id.labelidx.ctime	<- 4
	file			<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfgBwhRc-_140811_n390_AncSeq.R')
	cat(paste('\nLoading starting sequences from file', file))
	load(file)		#expect "anc.seq.gag"  "anc.seq.pol"  "anc.seq.env"  "anc.seq.info"
	setkey(anc.seq.info, CALENDAR_TIME)
	rANCSEQ.args<- list(	root.ctime.grace=root.ctime.grace, sample.grace=sample.grace, sample.shift=sample.shift, 
			anc.seq.info=anc.seq.info, anc.seq.gag=anc.seq.gag, anc.seq.pol=anc.seq.pol, anc.seq.env=anc.seq.env)	
	
	rANCSEQ<- function(root.ctime, rANCSEQ.args)
	{				
		tmp		<- lapply(seq_along(root.ctime), function(i)
				{
					tmp	<- subset(rANCSEQ.args$anc.seq.info, CALENDAR_TIME+rANCSEQ.args$sample.shift>root.ctime[i]-rANCSEQ.args$root.ctime.grace &  CALENDAR_TIME+rANCSEQ.args$sample.shift<=root.ctime[i]+rANCSEQ.args$root.ctime.grace)
					if(nrow(tmp)<=rANCSEQ.args$sample.grace*100)
					{
						print(c(nrow, rANCSEQ.args$sample.grace*100))
						stop()
					}
					data.table( LABEL= tmp[, sample( LABEL, rANCSEQ.args$sample.grace ) ], CALENDAR_TIME=root.ctime[i], DRAW=i )
				})
		tmp		<- do.call('rbind',tmp)
		#	get unique seqs
		setkey(tmp, LABEL)
		tmp				<- unique(tmp)	
		anc.seq.draw	<- do.call( 'cbind', list( rANCSEQ.args$anc.seq.gag[tmp[, LABEL], ], rANCSEQ.args$anc.seq.pol[tmp[, LABEL], ], rANCSEQ.args$anc.seq.env[tmp[, LABEL], ] ) ) 
		anc.seq.draw	<- seq.unique(anc.seq.draw)
		#	check if at least one seq for each draw
		tmp				<- merge( data.table(LABEL=rownames(anc.seq.draw)), tmp, by='LABEL' )	
		stopifnot( !length(setdiff( seq_along(root.ctime), tmp[, unique(DRAW)] )) )
		#	take first seq for each draw	
		tmp				<- tmp[, list(LABEL= LABEL[1]), by='DRAW']
		setkey(tmp, DRAW)
		anc.seq.draw	<- anc.seq.draw[ tmp[, LABEL], ]
		#	set new calendar time for sequences
		#set(tmp, NULL, 'LABEL_NEW', tmp[, as.numeric( sapply( strsplit(LABEL,tree.id.labelsep,fixed=TRUE), '[[', tree.id.labelidx.ctime) ) ])
		#set(tmp, NULL, 'LABEL_NEW', tmp[, LABEL_NEW+sample.shift])	
		#tmp			<- tmp[,	{
		#							z							<- strsplit(LABEL,tree.id.labelsep,fixed=TRUE)[[1]]
		#							z[tree.id.labelidx.ctime]	<- LABEL_NEW
		#							list(LABEL_NEW=paste(z, collapse=tree.id.labelsep,sep=''))
		#						}, by='LABEL']
		#setkey(tmp, LABEL)
		#rownames(anc.seq.draw)		<- tmp[ rownames(anc.seq.draw), ][, LABEL_NEW]		
		anc.seq.draw
	}
	list(rANCSEQ=rANCSEQ, rANCSEQ.args=rANCSEQ.args)
}
##--------------------------------------------------------------------------------------------------------
#	return ancestral sequence sampler	
#	olli originally written 11-09-2014
##--------------------------------------------------------------------------------------------------------
PANGEA.Seqsampler.v1<- function(df.ind, df.trm, pipeline.args, outfile.ind, outfile.trm, with.plot=1)
{	
	#	TODO is IDTR integer?
	# 	compute prevalence and incidence by year
	#	sampling model
	epi.adult	<- 13
	suppressWarnings( df.trm[, YR:= df.trm[, floor(TIME_TR)]] )
	df.epi		<- df.trm[, list(INC=length(IDREC), INC_ACUTE=length(which(TR_ACUTE=='Yes')),IMPORT=length(which(IDTR<0))), by='YR']
	tmp			<- df.epi[, 	{
				sexactive	<- which( floor(df.ind[['DOB']]+epi.adult)<=YR  &  ceiling(df.ind[['DOD']])>YR )
				infected	<- which( floor(df.ind[['DOB']])<=YR  &  floor(df.ind[['DOD']])>YR  &  floor(df.ind[['TIME_TR']])<=YR )
				infdead		<- which( floor(df.ind[['DOD']])==YR  &  floor(df.ind[['TIME_TR']])<=YR )
				list(POP=length(sexactive), PREV=length(infected), PREVDIED=length(infdead))				
			},by='YR']
	df.epi		<- merge( tmp, df.epi, by='YR' )	
	set(df.epi, NULL, 'PREVp', df.epi[, PREV/(POP-PREV)])	
	set(df.epi, NULL, 'INCp', df.epi[, INC/(POP-PREV)])
	set(df.epi, NULL, 'IMPORTp', df.epi[, IMPORT/INC])
	set(df.epi, NULL, 'ACUTEp', df.epi[, INC_ACUTE/INC])
	# 	SAMPLING PROBABILITIES and TOTALS PER YEAR
	#
	#	Can we detect a 25% or 50% reduction in HIV incidence in the most recent 2 or 3 years 
	#	with 1%, 5%, 10% of all recent incident cases sampled?
	#
	#	suppose exponentially increasing sampling over time
	#	the number of incident cases sampled is the total sampled in that year * the proportion of incident cases out of all non-sampled cases to date
	#	TODO this needs to be changed to fix the proportion of sequences sampled from incident
	s.PREV.base	<- exp(1)
	df.sample	<- subset( df.epi, YR>= pipeline.args['yr.start',][, as.numeric(v)] & YR<pipeline.args['yr.end',][, as.numeric(v)] )
	#	exponential rate of increasing s.TOTAL (total sampling rate) per year
	tmp			<- log( pipeline.args['s.PREV.max',][, as.numeric(v)]/pipeline.args['s.PREV.min',][, as.numeric(v)], base=s.PREV.base ) / df.sample[, diff(range(YR))]
	tmp			<- df.sample[, s.PREV.base^( tmp*(YR-min(YR)) ) * pipeline.args['s.PREV.min',][, as.numeric(v)] ]	
	#tmp			<- log( 1+pipeline.args['s.PREV.max',][, as.numeric(v)]-pipeline.args['s.PREV.min',][, as.numeric(v)], base=s.PREV.base ) / df.sample[, diff(range(YR))]
	#tmp			<- df.sample[, s.PREV.base^( tmp*(YR-min(YR)) ) - 1 + pipeline.args['s.PREV.min',][, as.numeric(v)] ]
	set(df.sample, NULL, 's.CUMTOTAL', tmp)		
	set(df.sample, NULL, 's.n.CUMTOTAL', df.sample[, round(PREV*s.CUMTOTAL)])
	set(df.sample, NULL, 's.n.TOTAL', c(df.sample[1, s.n.CUMTOTAL], df.sample[, diff(s.n.CUMTOTAL)]))	
	set(df.sample, NULL, 's.n.INC', df.sample[, round(INC/(PREV-s.n.CUMTOTAL) * s.n.TOTAL)])
	set(df.sample, NULL, 's.n.notINC', df.sample[, round(s.n.TOTAL-s.n.INC)])
	stopifnot(df.sample[, all(s.n.TOTAL>0)])
	cat(paste('\n total number of sequences sampled=', df.sample[, sum( s.n.TOTAL )]))
	cat(paste('\n prop of sequences sampled among HIV+=', df.sample[, sum( s.n.TOTAL )] / df.sample[, rev(PREV)[1]]))		
	cat(paste('\n total number of incident sequences to sample=', df.sample[, sum( s.n.INC )]))
	cat(paste('\n total number of non-incident sequences to sample=', df.sample[, sum( s.n.notINC )]))	
	#	SAMPLE INFECTED INDIVIDUALS BASED ON NUMBERS PER YEAR
	#
	#	sample incident cases by year
	df.inds	<- copy(df.ind)
	tmp		<- subset(df.trm, YR>= pipeline.args['yr.start',][, as.numeric(v)])
	tmp		<- tmp[, {
				z	<- df.sample[['s.n.INC']][ which(df.sample[['YR']]==YR) ]
				z	<- sample(seq_along(IDREC), z)
				list( 	IDPOP=IDREC[z], TIME_TR=TIME_TR[z], 
						TIME_SEQ=TIME_TR[z]+rexp(length(z), rate=1/(3*30))/365, 
						INCIDENT_SEQ=rep('Y',length(z) ) )
			}, by='YR']
	df.inds	<- merge(df.inds, subset(tmp, select=c(IDPOP, TIME_SEQ, INCIDENT_SEQ)), by='IDPOP', all.x=1)
	cat(paste('\n total number of incident sequences sampled=', df.inds[, length(which(!is.na(TIME_SEQ)))] ))	
	#	sample non-incident cases by year
	for(yr in df.sample[, YR])		#TODO took out [-1] because there are s.n.notINC for DSPS in 1980
	{
		#	of all infected and not incident and not yet sampled, sample
		cat(paste('\nadd non-incident samples in year',yr))
		tmp		<- subset(df.inds, is.na(TIME_SEQ) & floor(DOB)<=yr & ceiling(DOD)>yr & floor(TIME_TR)<yr)
		cat(paste('\navailable non-sampled non-incident cases in year=',nrow(tmp)))
		tmp		<- subset(tmp, T1_SEQ-0.5<=yr & T1_SEQ+0.5>yr & T1_SEQ-0.5-TIME_TR>1)
		cat(paste('\navailable with T1_SEQ +- 6mo, n=',nrow(tmp)))
		tmp2	<- df.sample[['s.n.notINC']][ which(df.sample[['YR']]==yr) ]
		stopifnot(tmp2<=nrow(tmp))
		tmp2	<- sample(seq_len(nrow(tmp)), tmp2)
		#	set variables in df.inds
		tmp		<- data.table(IDPOP= tmp[tmp2, IDPOP], TIME_SEQ=runif(length(tmp2), min=yr, max=yr+1), INCIDENT_SEQ=rep('N',length(tmp2) ))
		cat(paste('\nsampled non-incident cases in year=',nrow(tmp)))
		tmp2	<- sapply(tmp[,IDPOP], function(x) df.inds[,which(IDPOP==x)])
		set(df.inds, tmp2, 'TIME_SEQ', tmp[,TIME_SEQ])
		set(df.inds, tmp2, 'INCIDENT_SEQ', tmp[,INCIDENT_SEQ])		
	}
	cat(paste('\n total number of HIV+ in df.inds=', nrow(subset(df.inds, !is.na(TIME_TR)))))
	cat(paste('\n total number of sampled HIV+ in df.inds=', nrow(subset(df.inds, !is.na(TIME_TR) & !is.na(TIME_SEQ)))))
	#
	#	check that allocation OK
	#	
	set(df.inds, NULL, 'TIME_SEQYR', df.inds[, floor(TIME_SEQ)])
	tmp	<- subset(df.inds, !is.na(TIME_SEQ))[, list(s.n.TOTAL=length(IDPOP)), by='TIME_SEQYR']
	setkey(tmp, TIME_SEQYR)
	set(tmp,NULL,'s.n.CUMTOTAL',tmp[, cumsum(s.n.TOTAL)])
	stopifnot(  tmp[,tail(s.n.CUMTOTAL,1)]==df.sample[, tail(s.n.CUMTOTAL,1)] ) 
	#	set sampling in df.trm
	tmp		<- subset( df.inds, !is.na(TIME_SEQ), select=c(IDPOP, TIME_SEQ) )
	setnames(tmp, c('IDPOP','TIME_SEQ'), c('IDREC','SAMPLED_REC'))
	df.trms	<- merge(df.trm, tmp, by='IDREC', all.x=TRUE)
	setnames(tmp, c('IDREC','SAMPLED_REC'), c('IDTR','SAMPLED_TR'))
	df.trms	<- merge(df.trms, tmp, by='IDTR', all.x=TRUE)	
	#
	#	TRANSMISSION NETWORKS
	#
	require(igraph)
	#	cluster with index case
	tmp			<- subset(df.trms, select=c(IDTR, IDREC))			
	tmp			<- graph.data.frame(tmp, directed=TRUE, vertices=NULL)
	tmp			<- data.table(IDPOP=as.integer(V(tmp)$name), CLU=clusters(tmp, mode="weak")$membership)
	tmp2		<- tmp[, list(CLU_SIZE=length(IDPOP)), by='CLU']
	setkey(tmp2, CLU_SIZE)
	tmp2[, IDCLU:=rev(seq_len(nrow(tmp2)))]
	tmp			<- subset( merge(tmp, tmp2, by='CLU'), select=c(IDPOP, IDCLU) )
	setnames(tmp, 'IDPOP', 'IDREC')
	df.trms		<- merge( df.trms, tmp, by='IDREC', all.x=TRUE )
	stopifnot( nrow(subset(df.trms, is.na(IDCLU)))==0 )
	cat(paste('\nFound transmission clusters, n=', df.trms[, length(unique(IDCLU))]))
	#	add IDCLU to df.inds
	tmp			<- subset( df.trms, select=c(IDREC, IDTR, IDCLU) )
	tmp			<- subset( melt(tmp, id.var='IDCLU', value.name='IDPOP'), select=c(IDPOP, IDCLU))
	setkey(tmp, IDPOP, IDCLU)
	tmp			<- unique(tmp)
	df.inds		<- merge( df.inds, tmp, by='IDPOP', all.x=TRUE )
	#
	#	PLOTS
	#
	if(with.plot)
	{
		require(ggplot2)
		require(reshape2)
		#	plot numbers sampled, prevalent, incident
		set(df.sample, NULL, 'POP', df.sample[, as.real(POP)])
		set(df.sample, NULL, 'PREV', df.sample[, as.real(PREV)])
		set(df.sample, NULL, 'INC', df.sample[, as.real(INC)])
		tmp	<- data.table(	stat=c('POP','PREV','INC','ACUTEp','IMPORTp','s.n.TOTAL','s.n.INC','s.n.notINC'), 
				stat.long=c('population size','HIV infected', 'HIV incident', 'Proportion incident\nfrom acute in study', 'Proportion incident\nimported', 'Total\nsequenced', 'Total\nincident\nsequenced', 'Total\nnon-incident\nsequenced'))
		tmp	<- merge(	melt(df.sample, id.vars='YR', measure.vars=c('POP','PREV','INC','ACUTEp','IMPORTp','s.n.TOTAL','s.n.INC','s.n.notINC'), variable.name='stat', value.name='v'),
				tmp, by='stat' )
		ggplot(tmp, aes(x=YR, y=v, group=stat.long)) + geom_point() +
				scale_x_continuous(name='year', breaks=seq(1980,pipeline.args['yr.end',][, as.numeric(v)],2)) + scale_y_continuous(name='total')	+
				facet_grid(stat.long ~ ., scales='free_y', margins=FALSE)
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_Totals.pdf',sep='')
		cat(paste('\nPlotting to file',file))
		ggsave(file=file, w=16, h=10)		 
		#	plot distribution between transmission time and sequencing time
		tmp	<- subset(df.inds, !is.na(TIME_SEQ))
		set(tmp, NULL, 'TIME_TO_SEQ', tmp[, TIME_SEQ-TIME_TR])
		ggplot(tmp, aes(x=TIME_TO_SEQ)) + geom_histogram(binwidth=1) + 
				scale_x_continuous(name='time from transmission to sequence sampling\n(years)', breaks=seq(0,100,2))
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_Time2Seq.pdf',sep='')
		ggsave(file=file, w=8, h=8)
		#	plot distribution of #recipients per infector
		tmp	<- df.trms[, list(N= length(IDREC)), by='IDTR']
		ggplot(tmp, aes(x=N)) + geom_histogram(binwidth=1, aes(y= ..density..)) +
				scale_x_continuous(name='recipients per source case\n(number)', breaks=seq(1,100,1)+0.5, label=seq(1,100,1)) +
				scale_y_continuous(breaks=seq(0,1,0.05))
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_RecSource.pdf',sep='')
		ggsave(file=file, w=8, h=8)
		#	plot time to death for infected
		tmp	<- subset(df.inds, !is.na(TIME_TR) & IDPOP>0 & DOD<max(DOD, na.rm=1), select=c(TIME_TR, DOD))
		ggplot(tmp, aes(x=DOD-TIME_TR)) + geom_histogram(binwidth=1, aes(y= ..density..)) +
				scale_x_continuous(name='time to death for HIV+\n(years)', breaks=seq(1,100,2))
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_T2DeathForInf.pdf',sep='')
		ggsave(file=file, w=8, h=8)		
		#	plot transmission network
		file		<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_TrNetworks.pdf',sep='')
		cat(paste('\nPlotting to file',file))
		pdf(file=file, w=20, h=20)
		dummy	<- sapply( df.inds[, sort(na.omit(unique(IDCLU)))], function(clu)
				{
					cat(paste('\nprocess cluster no',clu))
					tmp					<- subset(df.inds, IDCLU==clu & IDPOP>=0, select=c(IDPOP, GENDER, TIME_SEQ))
					tmp[, IS_SEQ:= tmp[, factor(!is.na(TIME_SEQ), label=c('N','Y'), levels=c(FALSE, TRUE))]]
					clu.igr				<- graph.data.frame(subset(df.trms, IDCLU==clu & IDTR>=0, select=c(IDTR, IDREC)), directed=TRUE, vertices=subset(tmp, select=c(IDPOP, GENDER, IS_SEQ)))
					V(clu.igr)$color	<- ifelse( get.vertex.attribute(clu.igr, 'IS_SEQ')=='Y', 'green', 'grey90' )
					V(clu.igr)$shape	<- ifelse( get.vertex.attribute(clu.igr, 'GENDER')=='M', 'circle', 'square' )
					
					par(mai=c(0,0,1,0))
					plot(clu.igr, main=paste('IDCLU=',clu,sep=''), vertex.size=2, vertex.label.cex=0.25, edge.arrow.size=0.5, layout=layout.fruchterman.reingold(clu.igr, niter=1e3) )
					legend('bottomright', fill=c('green','grey90'), legend=c('sequence sampled','sequence not sampled'), bty='n')
					legend('bottomleft', legend=c('square= Female','circle= Male'), bty='n')				
				})
		dev.off()
		#ggplot(df.trms, aes(x=IDTR, y=TIME_TR)) + geom_point()
		#ggplot(df.trms, aes(x=IDTR, y=IDCLU)) + geom_point()
	}
	#
	#	SAVE SAMPLED RECIPIENTS AND TRANSMISSIONS TO SAMPLED RECIPIENTS
	#
	#	save for us
	file		<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'SAVE.R',sep='')
	save(file=file, df.epi, df.trms, df.inds, df.sample)
	#	save for virus tree simulator
	#	exclude columns that are not needed	
	df.inds	<- subset(df.inds, !is.na(TIME_TR))
	if('RISK'%in%colnames(df.inds))
		df.inds[, RISK:=NULL]
	if('INCIDENT_SEQ'%in%colnames(df.inds))
		df.inds[, INCIDENT_SEQ:=NULL]
	if('TIME_SEQYR'%in%colnames(df.inds))
		df.inds[, TIME_SEQYR:=NULL]	
	if('TR_ACUTE'%in%colnames(df.trms))
		df.trms[, TR_ACUTE:=NULL]
	if('YR'%in%colnames(df.trms))
		df.trms[, YR:=NULL]	
	#	add columns that the virus tree simulator needs
	if(!'IDTR_TIME_INFECTED'%in%colnames(df.trms))
	{
		tmp		<- subset( df.inds, !is.na(TIME_TR), c(IDPOP, TIME_TR) )
		setnames(tmp, c('IDPOP','TIME_TR'), c('IDTR','IDTR_TIME_INFECTED') )
		df.trms	<- merge(df.trms, tmp, by='IDTR', all.x=TRUE)		
	}
	cat(paste('\nwrite to file',outfile.ind))
	write.csv(file=outfile.ind, df.inds)
	cat(paste('\nwrite to file',outfile.trm))
	write.csv(file=outfile.trm, df.trms)
}
##--------------------------------------------------------------------------------------------------------
#	return ancestral sequence sampler	
#	olli originally written 24-10-2014
##--------------------------------------------------------------------------------------------------------
PANGEA.Seqsampler.v3<- function(df.ind, df.trm, pipeline.args, outfile.ind, outfile.trm, with.plot=1)
{	
	#	TODO is IDTR integer?
	# 	compute prevalence and incidence by year
	#	sampling model
	epi.adult	<- 13
	suppressWarnings( df.trm[, YR:= df.trm[, floor(TIME_TR)]] )
	df.epi		<- df.trm[, list(INC=length(IDREC), INC_ACUTE=length(which(TR_ACUTE=='Yes')),IMPORT=length(which(IDTR<0))), by='YR']
	tmp			<- df.epi[, 	{
				sexactive	<- which( floor(df.ind[['DOB']]+epi.adult)<=YR  &  ceiling(df.ind[['DOD']])>YR )
				infected	<- which( floor(df.ind[['DOB']])<=YR  &  floor(df.ind[['DOD']])>YR  &  floor(df.ind[['TIME_TR']])<=YR )
				infdead		<- which( floor(df.ind[['DOD']])==YR  &  floor(df.ind[['TIME_TR']])<=YR )
				list(POP=length(sexactive), PREV=length(infected), PREVDIED=length(infdead))				
			},by='YR']
	df.epi		<- merge( tmp, df.epi, by='YR' )	
	set(df.epi, NULL, 'PREVp', df.epi[, PREV/(POP-PREV)])	
	set(df.epi, NULL, 'INCp', df.epi[, INC/(POP-PREV)])
	set(df.epi, NULL, 'IMPORTp', df.epi[, IMPORT/INC])
	set(df.epi, NULL, 'ACUTEp', df.epi[, INC_ACUTE/INC])
	set(df.epi, NULL, 'GROWTHr', c(NA_real_, df.epi[, diff(log(PREV))]))	
	# 	SAMPLING PROBABILITIES and TOTALS PER YEAR
	#
	#	Can we detect a 25% or 50% reduction in HIV incidence in the most recent 2 or 3 years 
	#	with 1%, 5%, 10% of all recent incident cases sampled?
	#
	#	suppose exponentially increasing sampling over time
	#	the number of incident cases sampled is the total sampled in that year * the proportion of incident cases out of all non-sampled cases to date
	#	TODO this needs to be changed to fix the proportion of sequences sampled from incident
	s.PREV.base	<- exp(1)
	df.sample	<- subset( df.epi, YR>= pipeline.args['yr.start',][, as.numeric(v)] & YR<pipeline.args['yr.end',][, as.numeric(v)] )
	#	exponential rate of increasing s.TOTAL (total sampling rate) per year
	tmp			<- log( pipeline.args['s.PREV.max',][, as.numeric(v)]/pipeline.args['s.PREV.min',][, as.numeric(v)], base=s.PREV.base ) / df.sample[, diff(range(YR))]
	tmp			<- df.sample[, s.PREV.base^( tmp*(YR-min(YR)) ) * pipeline.args['s.PREV.min',][, as.numeric(v)] ]	
	#tmp			<- log( 1+pipeline.args['s.PREV.max',][, as.numeric(v)]-pipeline.args['s.PREV.min',][, as.numeric(v)], base=s.PREV.base ) / df.sample[, diff(range(YR))]
	#tmp			<- df.sample[, s.PREV.base^( tmp*(YR-min(YR)) ) - 1 + pipeline.args['s.PREV.min',][, as.numeric(v)] ]
	set(df.sample, NULL, 's.CUMTOTAL', tmp)		
	set(df.sample, NULL, 's.n.CUMTOTAL', df.sample[, round(PREV*s.CUMTOTAL)])
	set(df.sample, NULL, 's.n.TOTAL', c(df.sample[1, s.n.CUMTOTAL], df.sample[, diff(s.n.CUMTOTAL)]))	
	set(df.sample, NULL, 's.n.INC', df.sample[, round(INC/(PREV-s.n.CUMTOTAL) * s.n.TOTAL)])
	set(df.sample, NULL, 's.n.notINC', df.sample[, round(s.n.TOTAL-s.n.INC)])
	stopifnot(df.sample[, all(s.n.TOTAL>=0)])
	cat(paste('\n total number of sequences sampled=', df.sample[, sum( s.n.TOTAL )]))
	cat(paste('\n prop of sequences sampled among HIV+=', df.sample[, sum( s.n.TOTAL )] / df.sample[, rev(PREV)[1]]))		
	cat(paste('\n total number of incident sequences to sample=', df.sample[, sum( s.n.INC )]))
	cat(paste('\n total number of non-incident sequences to sample=', df.sample[, sum( s.n.notINC )]))	
	#	SAMPLE INFECTED INDIVIDUALS BASED ON NUMBERS PER YEAR
	#
	#	sample incident cases by year
	df.inds	<- copy(df.ind)
	tmp		<- subset(df.trm, YR>= pipeline.args['yr.start',][, as.numeric(v)])
	tmp		<- tmp[, {
				z	<- df.sample[['s.n.INC']][ which(df.sample[['YR']]==YR) ]
				z	<- sample(seq_along(IDREC), z)
				list( 	IDPOP=IDREC[z], TIME_TR=TIME_TR[z], 
						TIME_SEQ=TIME_TR[z]+rexp(length(z), rate=1/(3*30))/365, 
						INCIDENT_SEQ=rep('Y',length(z) ) )
			}, by='YR']
	df.inds	<- merge(df.inds, subset(tmp, select=c(IDPOP, TIME_SEQ, INCIDENT_SEQ)), by='IDPOP', all.x=1)
	tmp		<- df.inds[, which(TIME_SEQ>DOD)]
	set(df.inds, tmp, 'TIME_SEQ', df.inds[tmp, TIME_TR+(DOD-TIME_TR)/2])
	cat(paste('\n total number of incident sequences sampled=', df.inds[, length(which(!is.na(TIME_SEQ)))] ))	
	#	sample non-incident cases by year
	for(yr in df.sample[, YR])		#TODO took out [-1] because there are s.n.notINC for DSPS in 1980
	{
		#	of all infected and not incident and not yet sampled, sample
		cat(paste('\nadd non-incident samples in year',yr))
		tmp		<- subset(df.inds, is.na(TIME_SEQ) & floor(DOB)<=yr & ceiling(DOD)>yr+1 & yr-TIME_TR>1)
		cat(paste('\navailable non-sampled non-incident cases in year=',nrow(tmp)))
		tmp		<- subset(tmp, yr-TIME_TR<12)
		cat(paste('\navailable with TIME_TR at most 12 years earlier, n=',nrow(tmp)))		
		tmp		<- subset(tmp, is.na(T1_SEQ) | (T1_SEQ-0.5<=yr & T1_SEQ+0.5>yr & T1_SEQ-0.5-TIME_TR>1))
		cat(paste('\navailable with T1_SEQ +- 6mo, n=',nrow(tmp)))
		tmp2	<- df.sample[['s.n.notINC']][ which(df.sample[['YR']]==yr) ]
		stopifnot(tmp2<=nrow(tmp))
		tmp2	<- sample(seq_len(nrow(tmp)), tmp2)
		#	set variables in df.inds		
		tmp		<- data.table(IDPOP= tmp[tmp2, IDPOP], TIME_SEQ=runif(length(tmp2), min=yr, max=yr+1), INCIDENT_SEQ=rep('N',length(tmp2) ))
		cat(paste('\nsampled non-incident cases in year=',nrow(tmp)))
		tmp2	<- sapply(tmp[,IDPOP], function(x) df.inds[,which(IDPOP==x)])
		set(df.inds, tmp2, 'TIME_SEQ', tmp[,TIME_SEQ])
		set(df.inds, tmp2, 'INCIDENT_SEQ', tmp[,INCIDENT_SEQ])		
	}
	stopifnot( df.inds[, !any(TIME_SEQ>DOD, na.rm=TRUE)] )
	stopifnot( df.inds[, !any(TIME_SEQ<TIME_TR, na.rm=TRUE)] )	
	cat(paste('\n total number of HIV+ in df.inds=', nrow(subset(df.inds, !is.na(TIME_TR)))))
	cat(paste('\n total number of sampled HIV+ in df.inds=', nrow(subset(df.inds, !is.na(TIME_TR) & !is.na(TIME_SEQ)))))	
	#
	#	interpolate CD4 count at time of sampling
	#
	tmp		<- subset(df.inds, !is.na(TIME_SEQ), select=c(IDPOP, TIME_TR, INCIDENT_SEQ, TIME_SEQ, T1_CD4_500, T1_CD4_350, T1_CD4_200, DOD))
	set(tmp, NULL, 'TIME_TR', tmp[, TIME_TR-0.05])
	set(tmp, NULL, 'DOD', tmp[, DOD+0.05])
	setnames(tmp, c('TIME_TR','DOD'), c('T1_CD4_1000','T1_CD4_100'))
	tmp		<- melt(tmp, id.vars=c('IDPOP','TIME_SEQ','INCIDENT_SEQ'), value.name='T1_CD4', variable.name='CD4')
	set(tmp, NULL, 'CD4', tmp[, as.numeric(sapply(strsplit(as.character(CD4),'_'),'[[',3)) ])
	tmp2	<- tmp[, which(CD4==1000)]
	set(tmp, tmp2, 'CD4', runif(length(tmp2), min=700, max=1000))
	#tmp		<- tmp[, list(TIME_SEQ=TIME_SEQ[1], CD4_SEQ=exp(approx(T1_CD4, log(CD4), xout=TIME_SEQ[1], rule=1)$y)), by='IDPOP']
	tmp		<- tmp[, list(TIME_SEQ=TIME_SEQ[1], CD4_SEQ=round(approx(T1_CD4, CD4, xout=TIME_SEQ[1], rule=1)$y)), by='IDPOP']
	stopifnot( tmp[, !any(is.na(CD4_SEQ))] )
	df.inds	<- merge(df.inds, tmp, all.x=TRUE, by=c('IDPOP','TIME_SEQ'))
	#
	#	check that allocation OK
	#	
	set(df.inds, NULL, 'TIME_SEQYR', df.inds[, floor(TIME_SEQ)])
	tmp	<- subset(df.inds, !is.na(TIME_SEQ))[, list(s.n.TOTAL=length(IDPOP)), by='TIME_SEQYR']
	setkey(tmp, TIME_SEQYR)
	set(tmp,NULL,'s.n.CUMTOTAL',tmp[, cumsum(s.n.TOTAL)])
	stopifnot(  tmp[,tail(s.n.CUMTOTAL,1)]==df.sample[, tail(s.n.CUMTOTAL,1)] ) 
	#	set sampling in df.trm
	tmp		<- subset( df.inds, !is.na(TIME_SEQ), select=c(IDPOP, TIME_SEQ) )
	setnames(tmp, c('IDPOP','TIME_SEQ'), c('IDREC','SAMPLED_REC'))
	df.trms	<- merge(df.trm, tmp, by='IDREC', all.x=TRUE)
	setnames(tmp, c('IDREC','SAMPLED_REC'), c('IDTR','SAMPLED_TR'))
	df.trms	<- merge(df.trms, tmp, by='IDTR', all.x=TRUE)	
	#
	#	TRANSMISSION NETWORKS
	#
	require(igraph)
	#	cluster with index case
	tmp			<- subset(df.trms, select=c(IDTR, IDREC))			
	tmp			<- graph.data.frame(tmp, directed=TRUE, vertices=NULL)
	tmp			<- data.table(IDPOP=as.integer(V(tmp)$name), CLU=clusters(tmp, mode="weak")$membership)
	tmp2		<- tmp[, list(CLU_SIZE=length(IDPOP)), by='CLU']
	setkey(tmp2, CLU_SIZE)
	tmp2[, IDCLU:=rev(seq_len(nrow(tmp2)))]
	tmp			<- subset( merge(tmp, tmp2, by='CLU'), select=c(IDPOP, IDCLU) )
	setnames(tmp, 'IDPOP', 'IDREC')
	df.trms		<- merge( df.trms, tmp, by='IDREC', all.x=TRUE )
	stopifnot( nrow(subset(df.trms, is.na(IDCLU)))==0 )
	cat(paste('\nFound transmission clusters, n=', df.trms[, length(unique(IDCLU))]))
	#	add IDCLU to df.inds
	tmp			<- subset( df.trms, select=c(IDREC, IDTR, IDCLU) )
	tmp			<- subset( melt(tmp, id.var='IDCLU', value.name='IDPOP'), select=c(IDPOP, IDCLU))
	setkey(tmp, IDPOP, IDCLU)
	tmp			<- unique(tmp)
	df.inds		<- merge( df.inds, tmp, by='IDPOP', all.x=TRUE )
	#
	#	PLOTS
	#
	if(with.plot)
	{
		require(ggplot2)
		require(reshape2)
		#	plot numbers sampled, prevalent, incident
		set(df.sample, NULL, 'POP', df.sample[, as.numeric(POP)])
		set(df.sample, NULL, 'PREV', df.sample[, as.numeric(PREV)])
		set(df.sample, NULL, 'INC', df.sample[, as.numeric(INC)])
		tmp	<- data.table(	stat=c('POP','PREV','INC','ACUTEp','IMPORTp','s.n.TOTAL','s.n.INC','s.n.notINC','GROWTHr'), 
				stat.long=c('population size','HIV infected', 'HIV incident', 'Proportion incident\nfrom acute in study', 'Proportion incident\nimported', 'Total\nsequenced', 'Total\nincident\nsequenced', 'Total\nnon-incident\nsequenced','Annual growth rate'))
		tmp	<- merge(	melt(df.sample, id.vars='YR', measure.vars=c('POP','PREV','INC','ACUTEp','IMPORTp','s.n.TOTAL','s.n.INC','s.n.notINC','GROWTHr'), variable.name='stat', value.name='v'),
				tmp, by='stat' )
		ggplot(tmp, aes(x=YR, y=v, group=stat.long)) + geom_point() +
				scale_x_continuous(name='year', breaks=seq(1980,pipeline.args['yr.end',][, as.numeric(v)],2)) + scale_y_continuous(name='total')	+
				facet_grid(stat.long ~ ., scales='free_y', margins=FALSE)
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_Totals.pdf',sep='')
		cat(paste('\nPlotting to file',file))
		ggsave(file=file, w=16, h=16)	
		#	plot CD4 counts at time of sampling
		ggplot(subset(df.inds, !is.na(TIME_SEQ)), aes(x=TIME_SEQ, y=CD4_SEQ, colour=INCIDENT_SEQ)) + geom_point() + geom_smooth() + scale_y_continuous(breaks=seq(100,1000,100))		
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_CD4.pdf',sep='')
		cat(paste('\nPlotting to file',file))
		ggsave(file=file, w=8, h=8)					
		ggplot(subset(df.inds, !is.na(TIME_SEQ)), aes(x=floor(TIME_SEQ), fill=cut(CD4_SEQ, breaks=c(0,200,350,500,700,2000)))) + geom_bar(binwidth=1, position='fill') + 
				labs(fill='CD4 category', y='percent', x='year sequenced') + scale_y_continuous(breaks=seq(0,1,0.1)) + scale_fill_brewer(palette='Set1')		
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_CD4PROP.pdf',sep='')
		cat(paste('\nPlotting to file',file))
		ggsave(file=file, w=8, h=8)			
		#	plot distribution between transmission time and sequencing time
		tmp	<- subset(df.inds, !is.na(TIME_SEQ))
		set(tmp, NULL, 'TIME_TO_SEQ', tmp[, TIME_SEQ-TIME_TR])
		ggplot(tmp, aes(x=TIME_TO_SEQ)) + geom_histogram(binwidth=1) + 
				scale_x_continuous(name='time from transmission to sequence sampling\n(years)', breaks=seq(0,100,2))
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_Time2Seq.pdf',sep='')
		ggsave(file=file, w=8, h=8)
		#	plot distribution of #recipients per infector
		tmp	<- df.trms[, list(N= length(IDREC)), by='IDTR']
		ggplot(tmp, aes(x=N)) + geom_histogram(binwidth=1, aes(y= ..density..)) +
				scale_x_continuous(name='recipients per source case\n(number)', breaks=seq(1,100,1)+0.5, label=seq(1,100,1)) +
				scale_y_continuous(breaks=seq(0,1,0.05))
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_RecSource.pdf',sep='')
		ggsave(file=file, w=8, h=8)
		#	plot time to death for infected
		tmp	<- subset(df.inds, !is.na(TIME_TR) & IDPOP>0 & DOD<max(DOD, na.rm=1), select=c(TIME_TR, DOD))
		ggplot(tmp, aes(x=DOD-TIME_TR)) + geom_histogram(binwidth=1, aes(y= ..density..)) +
				scale_x_continuous(name='time to death for HIV+\n(years)', breaks=seq(1,100,2))
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_T2DeathForInf.pdf',sep='')
		ggsave(file=file, w=8, h=8)		
		#	plot transmission network
		file		<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_TrNetworks.pdf',sep='')
		cat(paste('\nPlotting to file',file))
		pdf(file=file, w=20, h=20)
		dummy	<- sapply( df.inds[, sort(na.omit(unique(IDCLU)))], function(clu)
				{
					cat(paste('\nprocess cluster no',clu))
					tmp					<- subset(df.inds, IDCLU==clu & IDPOP>=0, select=c(IDPOP, GENDER, TIME_SEQ))
					tmp[, IS_SEQ:= tmp[, factor(!is.na(TIME_SEQ), label=c('N','Y'), levels=c(FALSE, TRUE))]]
					clu.igr				<- graph.data.frame(subset(df.trms, IDCLU==clu & IDTR>=0, select=c(IDTR, IDREC)), directed=TRUE, vertices=subset(tmp, select=c(IDPOP, GENDER, IS_SEQ)))
					V(clu.igr)$color	<- ifelse( get.vertex.attribute(clu.igr, 'IS_SEQ')=='Y', 'green', 'grey90' )
					V(clu.igr)$shape	<- ifelse( get.vertex.attribute(clu.igr, 'GENDER')=='M', 'circle', 'square' )
					
					par(mai=c(0,0,1,0))
					plot(clu.igr, main=paste('IDCLU=',clu,sep=''), vertex.size=2, vertex.label.cex=0.25, edge.arrow.size=0.5, layout=layout.fruchterman.reingold(clu.igr, niter=1e3) )
					legend('bottomright', fill=c('green','grey90'), legend=c('sequence sampled','sequence not sampled'), bty='n')
					legend('bottomleft', legend=c('square= Female','circle= Male'), bty='n')				
				})
		dev.off()
		#ggplot(df.trms, aes(x=IDTR, y=TIME_TR)) + geom_point()
		#ggplot(df.trms, aes(x=IDTR, y=IDCLU)) + geom_point()
	}
	#
	#	SAVE SAMPLED RECIPIENTS AND TRANSMISSIONS TO SAMPLED RECIPIENTS
	#
	#	save for us
	file		<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'SAVE.R',sep='')
	save(file=file, df.epi, df.trms, df.inds, df.sample)
	#	save for virus tree simulator
	#	exclude columns that are not needed	
	df.inds	<- subset(df.inds, !is.na(TIME_TR))
	if('RISK'%in%colnames(df.inds))
		df.inds[, RISK:=NULL]
	if('INCIDENT_SEQ'%in%colnames(df.inds))
		df.inds[, INCIDENT_SEQ:=NULL]
	if('TIME_SEQYR'%in%colnames(df.inds))
		df.inds[, TIME_SEQYR:=NULL]	
	if('TR_ACUTE'%in%colnames(df.trms))
		df.trms[, TR_ACUTE:=NULL]
	if('YR'%in%colnames(df.trms))
		df.trms[, YR:=NULL]	
	#	add columns that the virus tree simulator needs
	if(!'IDTR_TIME_INFECTED'%in%colnames(df.trms))
	{
		tmp		<- subset( df.inds, !is.na(TIME_TR), c(IDPOP, TIME_TR) )
		setnames(tmp, c('IDPOP','TIME_TR'), c('IDTR','IDTR_TIME_INFECTED') )
		df.trms	<- merge(df.trms, tmp, by='IDTR', all.x=TRUE)		
	}
	cat(paste('\nwrite to file',outfile.ind))
	write.csv(file=outfile.ind, df.inds)
	cat(paste('\nwrite to file',outfile.trm))
	write.csv(file=outfile.trm, df.trms)
}
##--------------------------------------------------------------------------------------------------------
#	simulate guide to sequence sampling times. if not NA, then in every year, individuals in +-6mo to the guide are sampled	
#	olli originally written 24-10-2014
##--------------------------------------------------------------------------------------------------------
PANGEA.Seqsampler.SimulateGuideToSamplingTimes<- function(df.ind, seqtime.mode)
{
	stopifnot(seqtime.mode%in%c('Gamma3','Gamma9','Unif12'))
	cat(paste('\nUsing seqtime.mode=', seqtime.mode ))
	if(seqtime.mode=='Gamma3')
	{
		library(distr) 
		tmp		<- Gammad(shape=3, scale=2)
		tmp		<- Truncate(tmp, lower=0, upper=8)
		df.ind[, T1_SEQ:= df.ind[, r(tmp)(nrow(df.ind)) + TIME_TR]]			
	}
	if(seqtime.mode=='Gamma9')
	{
		df.ind[, T1_SEQ:= df.ind[, rgamma(nrow(df.ind),shape=9,scale=0.25 ) + TIME_TR]]	
	}
	if(seqtime.mode=='Unif12')
	{
		df.ind[, T1_SEQ:= runif(nrow(df.ind), 0, 12) + df.ind[,TIME_TR] ]
	}
	df.ind
}
##--------------------------------------------------------------------------------------------------------
#	return distribution of GTR parameters	
#	olli originally written 10-09-2014
##--------------------------------------------------------------------------------------------------------
PANGEA.GTR.params.v2<- function()
{		
	file		<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfgBwhRc-_140902_n390_BEASTlog.R')	
	cat(paste('\nreading GTR parameters from file',file))
	load(file)	# expect log.df
	#	exclude odd BEAST runs
	log.df		<- subset(log.df, !(GENE=='GAG' & FILE=='pool1'))
	log.df		<- subset(log.df, !(GENE=='POL' & FILE=='pool2'))
	#	all ENV are a bit odd ... 
	#	exclude cols
	log.df[, ucld.mean:=NULL]
	log.df[, ucld.stdev:=NULL]
	log.df[, coefficientOfVariation:=NULL]
	log.df[, treeModel.rootHeight:=NULL]
	#	set mean meanRate and put all variation into the mu's
	tmp		<- log.df[, mean(meanRate)]
	set(log.df, NULL, 'mu', log.df[, mu * meanRate / tmp])
	set(log.df, NULL, 'meanRate', tmp)
	log.df
}
##--------------------------------------------------------------------------------------------------------
#	return distribution of GTR parameters	
#	olli originally written 09-09-2014
##--------------------------------------------------------------------------------------------------------
PANGEA.GTR.params.v1<- function()
{	
	file		<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfgBwhRc-_140811_n390_BEASTlog.R')		
	cat(paste('\nreading GTR parameters from file',file))
	load(file)	# expect log.df
	log.df[, state:=NULL]
	log.df[, ucldmean:=NULL]
	log.df[, ucldstdev:=NULL]
	log.df[, treeLikelihood:=NULL]
	log.df[, FILE:=NULL]
	log.df
}
##--------------------------------------------------------------------------------------------------------
#	return distribution of GTR parameters	
#	olli originally written 14-09-2014
##--------------------------------------------------------------------------------------------------------
PANGEA.GTR.params.v3<- function()
{		
	file		<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfgBwhRc-_140907_n390_BEASTlog.R')	
	cat(paste('\nreading GTR parameters from file',file))
	load(file)	# expect log.df
	#	exclude odd BEAST runs
	log.df		<- subset(log.df, !(GENE=='ENV' & FILE=='pool3'))
	#	exclude cols
	log.df[, ucld.mean:=NULL]
	log.df[, ucld.stdev:=NULL]
	log.df[, coefficientOfVariation:=NULL]
	log.df[, treeModel.rootHeight:=NULL]
	#	get mean rate. need to be a bit careful here: rate is per site, so length of gene does not matter
	#	but we may have different numbers of samples for each gene
	tmp		<- log.df[, list(meanRate=mean(meanRate)), by='GENE'][, mean(meanRate)]	
	#	set mean meanRate and put all variation into the mu's
	set(log.df, NULL, 'mu', log.df[, mu * meanRate / tmp])
	set(log.df, NULL, 'meanRate', tmp)
	log.df
}
##--------------------------------------------------------------------------------------------------------
#	simulate imports during the epidemic	
#	olli originally written 13-09-2014
##--------------------------------------------------------------------------------------------------------
PANGEA.ImportSimulator.SimulateStartingTimeOfIndexCase<- function(df.ind, df.trm, index.starttime.mode='normal')
{
	stopifnot(index.starttime.mode%in%c('normal','fix','fix45','shift'))
	#	add transmission time for index case -- this is 40 years back in time so we can sample a starting sequence 
	#	and then generate a long branch to the transmission chain in the population. No hack. :-)
	tmp			<- subset( df.trm, IDTR<0, select=IDTR )	
	if( index.starttime.mode == 'normal' )
	{
		cat(paste('\nUsing index.starttime.mode rnorm(n, 1955, 7)'))
		tmp2		<- rnorm(2*nrow(tmp), 1955, 7)
		tmp2		<- tmp2[ tmp2>1946 & tmp2<1980]
		stopifnot( nrow(tmp)<=length(tmp2) )
		length(tmp2)<- nrow(tmp)			
		set(tmp, NULL, 'IDTR_TIME_INFECTED.new', tmp2 )
	}
	if( index.starttime.mode == 'fix' || index.starttime.mode == 'shift')
	{
		cat(paste('\nUsing index.starttime.mode runif( n, 1954.75, 1955.25 )'))
		tmp2		<- runif( nrow(tmp), 1954.75, 1955.25 )
		set(tmp, NULL, 'IDTR_TIME_INFECTED.new', tmp2 )
	}
	if( index.starttime.mode == 'fix45')
	{
		cat(paste('\nUsing index.starttime.mode runif( n, 1945, 1945.5 )'))
		tmp2		<- runif( nrow(tmp), 1945, 1945.5 )
		set(tmp, NULL, 'IDTR_TIME_INFECTED.new', tmp2 )
	}
	#
	df.trm		<- merge( df.trm, tmp, by='IDTR', all.x=TRUE )
	tmp2		<- df.trm[, which(!is.na(IDTR_TIME_INFECTED.new) & IDTR_TIME_INFECTED<IDTR_TIME_INFECTED.new)]
	set(df.trm, tmp2, 'IDTR_TIME_INFECTED.new', df.trm[tmp2, IDTR_TIME_INFECTED])
	#	
	tmp2		<- df.trm[, which(!is.na(IDTR_TIME_INFECTED.new))]
	set(df.trm, tmp2, 'IDTR_TIME_INFECTED', df.trm[tmp2, IDTR_TIME_INFECTED.new])
	df.trm[, IDTR_TIME_INFECTED.new:=NULL]
	#	
	stopifnot( nrow(subset(df.trm, TIME_TR<IDTR_TIME_INFECTED))==0 )
	stopifnot( nrow(subset(df.trm, is.na(TIME_TR)))==0 )
	stopifnot( nrow(subset(df.trm, is.na(IDTR_TIME_INFECTED)))==0 )
	#
	tmp			<- subset( df.trm, IDTR<0, select=c(IDTR, IDTR_TIME_INFECTED) )
	setnames(tmp, c('IDTR','IDTR_TIME_INFECTED'), c('IDPOP','IDTR_TIME_INFECTED.new'))
	df.ind		<- merge(df.ind, tmp, by='IDPOP', all.x=TRUE)
	tmp2		<- df.ind[, which(!is.na(IDTR_TIME_INFECTED.new))]
	set(df.ind, tmp2, 'TIME_TR', df.ind[tmp2, IDTR_TIME_INFECTED.new])
	df.ind[, IDTR_TIME_INFECTED.new:=NULL]
	list(df.ind=df.ind, df.trm=df.trm)
}
##--------------------------------------------------------------------------------------------------------
#	return ancestral sequence sampler	
#	olli originally written 14-09-2014
##--------------------------------------------------------------------------------------------------------
PANGEA.RootSeq.create.sampler.v3<- function(root.ctime.grace= 0.5, sample.grace= 3)
{	
	file			<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfgBwhRc-_140907_n390_AncSeq.R')
	cat(paste('\nLoading starting sequences from file', file))
	load(file)		#expect "anc.seq.gag"  "anc.seq.pol"  "anc.seq.env"  "anc.seq.info"
	setkey(anc.seq.info, CALENDAR_TIME)
	rANCSEQ.args	<- list(	root.ctime.grace=root.ctime.grace, sample.grace=sample.grace, anc.seq.info=anc.seq.info, anc.seq.gag=anc.seq.gag, anc.seq.pol=anc.seq.pol, anc.seq.env=anc.seq.env)	
	
	rANCSEQ<- function(root.ctime, rANCSEQ.args)
	{		
		tmp		<- lapply(seq_along(root.ctime), function(i)
				{
					tmp	<- subset(rANCSEQ.args$anc.seq.info, CALENDAR_TIME>root.ctime[i]-rANCSEQ.args$root.ctime.grace &  CALENDAR_TIME<=root.ctime[i]+rANCSEQ.args$root.ctime.grace)
					if(nrow(tmp)<rANCSEQ.args$sample.grace*100)
					{
						warning( paste('\nFor root',i,': number of samples is n=',nrow(tmp),'. safe pool size is n=',rANCSEQ.args$sample.grace*100) )
					}					
					data.table( LABEL= tmp[, sample( LABEL, rANCSEQ.args$sample.grace ) ], CALENDAR_TIME=root.ctime[i], DRAW=i )
				})
		tmp		<- do.call('rbind',tmp)
		#	get unique seqs
		setkey(tmp, LABEL)
		tmp				<- unique(tmp)	
		anc.seq.draw	<- do.call( 'cbind', list( rANCSEQ.args$anc.seq.gag[tmp[, LABEL], ], rANCSEQ.args$anc.seq.pol[tmp[, LABEL], ], rANCSEQ.args$anc.seq.env[tmp[, LABEL], ] ) ) 
		anc.seq.draw	<- seq.unique(anc.seq.draw)
		#	check if at least one seq for each draw
		tmp				<- merge( data.table(LABEL=rownames(anc.seq.draw)), tmp, by='LABEL' )	
		stopifnot( !length(setdiff( seq_along(root.ctime), tmp[, unique(DRAW)] )) )
		#	take first seq for each draw	
		tmp				<- tmp[, list(LABEL= LABEL[1]), by='DRAW']
		setkey(tmp, DRAW)
		anc.seq.draw	<- anc.seq.draw[ tmp[, LABEL], ]
		anc.seq.draw
	}
	list(rANCSEQ=rANCSEQ, rANCSEQ.args=rANCSEQ.args)
}
##--------------------------------------------------------------------------------------------------------
##	command line generator for 'prog.HPTN071.input.parser.v1'
##	olli originally written 19-08-2014
##--------------------------------------------------------------------------------------------------------
#' @title Command line generator for \code{HPTN071.input.parser.v1}
#' @return command line string
#' @example example/ex.seq.sampler.v1.R
cmd.HPTN071.input.parser.v1<- function(indir, infile.trm, infile.ind, infile.args, outdir, outfile.trm, outfile.ind, prog=PR.HPTN071.INPUT.PARSER1 )	
{
	cmd<- "#######################################################
# start: run HPTN071.input.parser.v1
			#######################################################"
	cmd		<- paste(cmd, paste("\necho \'run ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd, paste(prog,' -indir=', indir,' -infile.trm=',infile.trm,' -infile.ind=',infile.ind,' -infile.args=',infile.args,' -outdir=',outdir,' -outfile.ind=',outfile.ind,' -outfile.trm=',outfile.trm,' \n', sep=''))
	cmd		<- paste(cmd,paste("echo \'end ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd,"#######################################################
# end: run HPTN071.input.parser.v1
					#######################################################\n",sep='')
	cmd
}
##--------------------------------------------------------------------------------------------------------
##	command line generator for 'prog.HPTN071.input.parser.v2'
##	olli originally written 08-09-2014
##--------------------------------------------------------------------------------------------------------
#' @title Command line generator for \code{HPTN071.input.parser.v2}
#' @return command line string
#' @example example/ex.seq.sampler.v2.R
cmd.HPTN071.input.parser.v2<- function(indir, infile.trm, infile.ind, infile.args, outdir, outfile.trm, outfile.ind, prog=PR.HPTN071.INPUT.PARSER2 )	
{
	cmd<- "#######################################################
# start: run HPTN071.input.parser.v2
			#######################################################"
	cmd		<- paste(cmd, paste("\necho \'run ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd, paste(prog,' -indir=', indir,' -infile.trm=',infile.trm,' -infile.ind=',infile.ind,' -infile.args=',infile.args,' -outdir=',outdir,' -outfile.ind=',outfile.ind,' -outfile.trm=',outfile.trm,' \n', sep=''))
	cmd		<- paste(cmd,paste("echo \'end ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd,"#######################################################
# end: run HPTN071.input.parser.v2
					#######################################################\n",sep='')
	cmd
}
##--------------------------------------------------------------------------------------------------------
##	command line generator for 'prog.HPTN071.input.parser.v3'
##	olli originally written 08-09-2014
##--------------------------------------------------------------------------------------------------------
#' @title Command line generator for \code{HPTN071.input.parser.v2}
#' @return command line string
#' @example example/ex.seq.sampler.v2.R
cmd.HPTN071.input.parser.v3<- function(indir, infile.trm, infile.ind, infile.args, outdir, outfile.trm, outfile.ind, prog=PR.HPTN071.INPUT.PARSER3 )	
{
	cmd<- "#######################################################
# start: run HPTN071.input.parser.v3
			#######################################################"
	cmd		<- paste(cmd, paste("\necho \'run ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd, paste(prog,' -indir=', indir,' -infile.trm=',infile.trm,' -infile.ind=',infile.ind,' -infile.args=',infile.args,' -outdir=',outdir,' -outfile.ind=',outfile.ind,' -outfile.trm=',outfile.trm,' \n', sep=''))
	cmd		<- paste(cmd,paste("echo \'end ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd,"#######################################################
# end: run HPTN071.input.parser.v3
					#######################################################\n",sep='')
	cmd
}
##--------------------------------------------------------------------------------------------------------
##	command line generator for 'prog.HPTN071.input.parser.v2'
##	olli originally written 08-09-2014
##--------------------------------------------------------------------------------------------------------
#' @title Command line generator for \code{DSPS.input.parser.v2}
#' @return command line string
#' @example example/ex.seq.sampler.DSPS.v2.R
cmd.DSPS.input.parser.v2<- function(indir, infile.trm, infile.args, outdir, outfile.trm, outfile.ind, prog=PR.DSPS.INPUT.PARSER2 )	
{
	cmd<- "#######################################################
# start: run DSPS.input.parser.v2
			#######################################################"
	cmd		<- paste(cmd, paste("\necho \'run ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd, paste(prog,' -indir=', indir,' -infile.trm=',infile.trm,' -infile.args=',infile.args,' -outdir=',outdir,' -outfile.ind=',outfile.ind,' -outfile.trm=',outfile.trm,' \n', sep=''))
	cmd		<- paste(cmd,paste("echo \'end ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd,"#######################################################
# end: run DSPS.input.parser.v2
					#######################################################\n",sep='')
	cmd
}
##--------------------------------------------------------------------------------------------------------
##	olli originally written 08-09-2014
##--------------------------------------------------------------------------------------------------------
#' @title DSPS parser (version 2, includes simulation of imports)
#' @description Reads files from the DSPS epi simulator in directory \code{indir} and writes csv files
#' in directory \code{outdir} for the virus tree simulator. The program samples sequences according to
#' an exponentially increasing sampling fraction in the same way as \code{prog.HPTN071.input.parser.v1}.
#' In addition, transmissions are broken and treated as imported from outside the simulated population.
#' The infected of a broken transmission chain is considered a new index case of a transmission chain within the 
#' simulated population. All input arguments are specified via the \code{argv} 
#' string, see the Examples.
#' @return NULL. Saves simulations to file.
#' @example example/ex.seq.sampler.DSPS.v2.R
prog.DSPS.input.parser.v2<- function()
{
	require(data.table)
	verbose			<- 1
	with.plot		<- 1	
	pipeline.args	<- NULL
	indir			<- system.file(package="rPANGEAHIVsim", "misc")	
	outdir			<- '/Users/Oliver/git/HPTN071sim/tmp140911'
	infile.trm		<- '140911_DSPS_RUN001_TRM.csv'
	infile.args		<- NA
	outfile.ind		<- '140911_DSPS_RUN001_IND.csv'
	outfile.trm		<- '140911_DSPS_RUN001_TRM.csv'
	#
	if(exists("argv"))
	{
		#	args input
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.trm= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.trm<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									infile.args= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.args<- tmp[1]
		#	args output
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									outfile.ind= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) outfile.ind<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									outfile.trm= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) outfile.trm<- tmp[1]		
	}
	if(verbose)
	{
		cat('\ninput args\n',paste(indir, infile.trm, outdir, outfile.ind, outfile.trm, sep='\n'))
	}	
	if(!is.na(infile.args))
	{
		load(infile.args)	#expect 'pipeline.args'
	}
	if(is.null(pipeline.args))
	{
		cat('\nCould not find pipeline.args, generating default')
		pipeline.args	<- rPANGEAHIVsim.pipeline.args()
	}
	stopifnot( all( c('yr.start', 'yr.end', 's.seed', 's.PREV.min', 's.PREV.max', 'epi.dt', 'epi.import')%in%pipeline.args[, stat] ) )
	#
	#infile.ind	<- paste(indir, '/', infile.ind, sep='')
	infile.trm	<- paste(indir, '/', infile.trm, sep='')
	outfile.ind	<- paste(outdir, '/', outfile.ind, sep='')
	outfile.trm	<- paste(outdir, '/', outfile.trm, sep='')
	
	#	set seed
	set.seed( pipeline.args['s.seed',][, as.numeric(v)] )
	#
	#	prepare transmission data.table
	#
	df.trm	<- as.data.table(read.csv(infile.trm, stringsAsFactors=FALSE, sep=',', dec='.'))
	#	for now ignore SAMPLING: use sampling model consistent with HPTN071 sampling model
	df.trm	<- subset(df.trm, EventType%in%c('INFECTION','DEATH','BIRTH'))
	#	rescale time 
	set(df.trm, NULL, 'ActionTime', df.trm[, ActionTime+1950])
	#	since - is confusing with -1 change to _
	set(df.trm, NULL, 'FromDeme.FromHost', df.trm[, gsub('-','_',FromDeme.FromHost)])
	set(df.trm, NULL, 'ToDeme.ToHost', df.trm[, gsub('-','_',ToDeme.ToHost)])
	tmp		<- df.trm[, which(is.na(FromDeme.FromHost))]
	cat(paste('\nFound TRM entries with NA FromDeme.FromHost, n=', length(tmp)))
	set(df.trm, tmp, 'FromDeme.FromHost', paste('HouseNA_',-seq_along(tmp),'_NA',sep=''))
	#	get ID of transmitters
	set(df.trm, NULL, 'IDTR', as.numeric(df.trm[, sapply( strsplit(FromDeme.FromHost,'_'), '[[', 2 )]))
	#	get ID of recipients
	set(df.trm, NULL, 'IDREC', as.numeric(df.trm[, sapply( strsplit(ToDeme.ToHost,'_'), '[[', 2 )]))
	#	
	#	prepare metavariable data.table
	#
	df.ind	<- subset(df.trm, select=c(IDREC, ToDeme.ToHost, EventType, ActionTime))
	setnames(df.ind, c('IDREC','ToDeme.ToHost'), c('IDPOP','INFO'))
	tmp		<- subset(df.trm, select=c(IDTR, FromDeme.FromHost))
	setnames(tmp, c('IDTR','FromDeme.FromHost'), c('IDPOP','INFO'))
	df.ind	<- rbind(df.ind, tmp, useNames=TRUE, fill=TRUE)
	df.ind	<- df.ind[-nrow(df.ind), ]
	df.ind[, V1:=NULL]
	#	removing unrealistic IDPOP
	cat(paste('\nFound IDPOP>1e6, remove, n=',nrow(subset(df.ind, IDPOP>1e6)),'2Emma: these are birth events'))
	df.ind	<- subset(df.ind, IDPOP<1e6)
	#	get House ID
	set(df.ind, NULL, 'HOUSE', df.ind[, sapply( strsplit(INFO,'_'), '[[', 1 )])
	suppressWarnings( set(df.ind, NULL, 'HOUSE', df.ind[, as.integer(substr(HOUSE,6,nchar(HOUSE)))]) )
	#	get GENDER
	set(df.ind, NULL, 'GENDER', df.ind[, sapply( strsplit(INFO,'_'), '[[', 3 )])
	set(df.ind, NULL, 'GENDER', df.ind[, factor(GENDER, levels=c('FEMALE','MALE'), labels=c('F','M'))])	
	#	add DOD	date of death if any
	tmp		<- subset(df.ind, EventType=='DEATH', select=c(IDPOP, ActionTime))
	setkey(tmp, IDPOP)
	df.ind	<- merge( subset(df.ind, is.na(EventType) | EventType!='DEATH'), unique(tmp), by='IDPOP', all.x=TRUE)	
	setnames(df.ind, c('ActionTime.x','ActionTime.y'), c('ActionTime','DOD'))
	#	add DOB	date of birth if any
	tmp		<- subset(df.ind, EventType=='BIRTH', select=c(IDPOP, ActionTime))
	setkey(tmp, IDPOP)
	df.ind	<- merge( subset(df.ind, is.na(EventType) | EventType!='BIRTH'), unique(tmp), by='IDPOP', all.x=TRUE)
	setnames(df.ind, c('ActionTime.x','ActionTime.y'), c('ActionTime','DOB'))
	#	add time of infection if any
	tmp		<- subset(df.ind, EventType=='INFECTION', select=c(IDPOP, ActionTime))
	setkey(tmp, IDPOP)
	#tmp		<- subset(df.trm, EventType=='INFECTION', select=c(IDREC, ActionTime) )
	#setnames(tmp, c('IDREC','ActionTime'), c('IDPOP', 'TIME_TR'))	
	#merge( subset(df.ind, is.na(EventType) | EventType!='INFECTION'), unique(tmp), by='IDPOP', all.x=TRUE, all.y=TRUE)
	df.ind	<- merge( subset(df.ind, is.na(EventType) | EventType!='INFECTION'), unique(tmp), by='IDPOP', all.x=TRUE, all.y=TRUE)
	setnames(df.ind, c('ActionTime.x','ActionTime.y'), c('ActionTime','TIME_TR'))	
	#	nothing else to process
	stopifnot( df.ind[, all(is.na(EventType))] )
	df.ind[, EventType:=NULL]
	df.ind[, ActionTime:=NULL]	
	#
	#	clean up df.ind
	#
	setkey(df.ind, IDPOP)	
	df.ind	<- unique(df.ind)
	stopifnot( nrow(subset(df.ind, DOD<=DOB))==0 )
	stopifnot( nrow(subset(df.ind, TIME_TR<=DOB))==0 )
	set(df.ind, df.ind[, which(TIME_TR<=DOB)], 'DOB', NA_real_)
	stopifnot( nrow(subset(df.ind, TIME_TR>=DOD))==0 )
	df.ind[, INFO:=NULL]
	set( df.ind, NULL, 'IDPOP', df.ind[, as.integer(IDPOP)] )
	df.ind	<- subset( df.ind, !is.na(DOB) | !is.na(DOD) | !is.na(TIME_TR) | IDPOP<0 )
	tmp		<- df.ind[, which(IDPOP>0 & is.na(DOD))]
	cat(paste('\nFound individuals alive at simulation end, n=', length(tmp)))	
	set(df.ind, tmp, 'DOD', df.ind[, ceiling(max(max(DOD, na.rm=TRUE), max(TIME_TR, na.rm=TRUE))+1)] )
	tmp		<- df.ind[, which(IDPOP>0 & is.na(DOB))]
	cat(paste('\nFound individuals with no birth date, n=', length(tmp)))
	set(df.ind, tmp, 'DOB', df.ind[, floor(min(TIME_TR, na.rm=TRUE))-1] )
	cat(paste('\nFound individuals with a valid record, n=', nrow(df.ind)))
	cat(paste('\nFound individuals with an infection event, n=', nrow(subset(df.ind,!is.na(TIME_TR)))))
	cat(paste('\nFound index cases, n=', nrow(subset(df.ind,IDPOP<0))))
	#
	#	clean up df.trm
	#
	df.trm	<- subset(df.trm, EventType=='INFECTION')
	setnames(df.trm, c("ActionTime"), c('TIME_TR'))		
	df.trm[, EventType:=NULL]
	set( df.trm, NULL, 'IDTR', df.trm[, as.integer(IDTR)] )
	set( df.trm, NULL, 'IDREC', df.trm[, as.integer(IDREC)] )
	#	check that transmission happen at unique times;	this seems to be the case in the DSPS model
	stopifnot( nrow(df.trm)==df.trm[, length(unique(TIME_TR))] )	
	#	add time of infection of transmitter to df.trm	
	tmp		<- subset(df.ind, select=c(IDPOP, TIME_TR))
	setnames(tmp, c('IDPOP','TIME_TR'), c('IDTR','IDTR_TIME_INFECTED') )
	setkey(tmp, IDTR)
	df.trm	<- merge(df.trm, unique(tmp), by='IDTR', all.x=TRUE)	
	#	
	df.trm[, FromDeme.FromHost:=NULL]
	df.trm[, ToDeme.ToHost:=NULL]
	cat(paste('\nFound transmissions, n=', nrow(df.trm)))
	cat(paste('\nTotal transmitters, n=', df.trm[, length(unique(IDTR))]))
	#
	#	reduce to time frame of interest
	#
	tmp		<- subset( df.trm, TIME_TR>=as.numeric( pipeline.args['yr.end',][, as.numeric(v)] ) )[, IDREC]
	df.trm	<- subset( df.trm, TIME_TR<as.numeric( pipeline.args['yr.end',][, as.numeric(v)] ) )
	df.ind	<- subset(df.ind, !IDPOP%in%tmp)
	df.ind	<- subset(df.ind, is.na(DOB) | DOB<pipeline.args['yr.end',][, as.numeric(v)] )
	cat(paste('\nFound individuals born before',pipeline.args['yr.end',][, as.numeric(v)],', n=', nrow(df.ind)))
	cat(paste('\nFound transmissions before',pipeline.args['yr.end',][, as.numeric(v)],', n=', nrow(df.trm)))
	cat(paste('\nTotal transmitters, n=', df.trm[, length(unique(IDTR))]))		
	stopifnot( length(setdiff( df.trm[, IDTR], df.ind[, IDPOP] ))==0 )
	stopifnot( length(setdiff( df.trm[, IDREC], df.ind[, IDPOP] ))==0 )
	#	optional: multiple baseline index cases
	if(0)
	{
		#	consider only transmissions from between yr.start and yr.end
		df.trm	<- subset( df.trm, TIME_TR>=as.numeric( pipeline.args['yr.start',][, as.numeric(v)] ) & TIME_TR<as.numeric( pipeline.args['yr.end',][, as.numeric(v)] ) )
		#
		#	all transmitters infected before yr.start are re-declared as from 'outside the study'. This is coded with negative POPIDs.
		#
		setkey(df.trm, TIME_TR)
		tmp		<- unique( subset( df.trm, IDTR_TIME_INFECTED<pipeline.args['yr.start',][, as.numeric(v)], IDTR) )	
		tmp[, IDTR.new:=-rev(seq_len(nrow(tmp)))]
		df.trm	<- merge(df.trm, tmp, by='IDTR', all.x=TRUE)
		setnames(tmp, 'IDTR', 'IDPOP')
		df.ind	<- merge(df.ind, tmp, by='IDPOP', all.x=TRUE)
		tmp		<- df.ind[, which(!is.na(IDTR.new))]
		set(df.ind, tmp, 'IDPOP', df.ind[tmp, IDTR.new])
		tmp		<- df.trm[, which(!is.na(IDTR.new))]
		set(df.trm, tmp, 'IDTR', df.trm[tmp, IDTR.new]) 	
		df.trm[, IDTR.new:=NULL]
		df.ind[, IDTR.new:=NULL]
		#	delete all POPID that don t appear as IDTR or IDREC in df.trm
		tmp		<- melt(subset(df.trm,select=c(IDTR,IDREC)), measure.vars=c('IDTR','IDREC'), value.name='IDPOP')
		tmp		<- unique(subset(tmp, select=IDPOP))
		df.ind	<- merge(df.ind, tmp, by='IDPOP')
		#
		cat(paste('\nFound transmissions between',pipeline.args['yr.start',][, as.numeric(v)],'-',pipeline.args['yr.end',][, as.numeric(v)],', n=', nrow(df.trm)))
		cat(paste('\nTotal transmitters, n=', df.trm[, length(unique(IDTR))]))	
		cat(paste('\nTotal transmitters before',pipeline.args['yr.start',][, as.numeric(v)],', these are treated as index cases before study start, n=', subset(df.trm, IDTR<0)[, length(unique(IDTR))] ))		
	}
	#	simulate a fraction of transmissions to be imports
	tmp		<- PANGEA.ImportSimulator.SimulateIndexCase(df.ind, df.trm, epi.import= pipeline.args['epi.import',][,as.numeric(v)])
	df.trm	<- tmp$df.trm
	df.ind	<- tmp$df.ind
	#
	#	set infection times for index case
	#
	#	add IDTR_TIME_INFECTED for -1
	stopifnot( subset(df.trm, is.na(IDTR_TIME_INFECTED))[, unique(IDTR)]==-1 )
	tmp		<- df.trm[, which(is.na(IDTR_TIME_INFECTED))]	
	set( df.trm, tmp, 'IDTR_TIME_INFECTED', df.trm[tmp, runif(1, TIME_TR-5, TIME_TR)] )	
	tmp		<- PANGEA.ImportSimulator.SimulateStartingTimeOfIndexCase(df.ind, df.trm)
	df.trm	<- tmp$df.trm
	df.ind	<- tmp$df.ind	
	#
	#	sample sequences 
	#
	PANGEA.Seqsampler(df.ind, df.trm, pipeline.args, outfile.ind, outfile.trm, with.plot=with.plot)
	#
	return(1)
}
##--------------------------------------------------------------------------------------------------------
##	olli originally written 26-07-2014
##--------------------------------------------------------------------------------------------------------
#' @title HPTN071 parser (version 1)
#' @description Reads files from the epi simulator in directory \code{indir} and writes csv files
#' in directory \code{outdir} for the virus tree simulator. The program samples sequences according to
#' an exponentially increasing sampling fraction. All input arguments are specified via the \code{argv} 
#' string, see the Examples.
#' @return NULL. Saves simulations to file.
#' @example example/ex.seq.sampler.v1.R
prog.HPTN071.input.parser.v1<- function()	
{
	require(data.table)
	verbose			<- 1
	with.plot		<- 1	
	pipeline.args	<- NULL
	indir			<- '/Users/Oliver/git/HPTN071sim/raw_trchain'
	infile.ind		<- '140716_RUN001_IND.csv'
	infile.trm		<- '140716_RUN001_TRM.txt'
	infile.args		<- NA
	outdir			<- '/Users/Oliver/git/HPTN071sim/sim_trchain'
	outfile.ind		<- '140716_RUN001_IND.csv'
	outfile.trm		<- '140716_RUN001_TRM.csv'
	#	
	if(exists("argv"))
	{
		#	args input
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.ind= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.ind<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.trm= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.trm<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									infile.args= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.args<- tmp[1]		
		#	args output
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									outfile.ind= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) outfile.ind<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									outfile.trm= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) outfile.trm<- tmp[1]		
	}
	if(verbose)
	{
		cat('\ninput args\n',paste(indir, infile.ind, infile.trm, outdir, outfile.ind, outfile.trm, sep='\n'))
	}	
	if(!is.na(infile.args))
	{
		load(infile.args)	#expect 'pipeline.args'
	}
	if(is.null(pipeline.args))
	{
		cat('\nCould not find pipeline.args, generating default')
		pipeline.args	<- rPANGEAHIVsim.pipeline.args()
	}
	stopifnot( all( c('yr.start', 'yr.end', 's.seed', 's.PREV.min', 's.PREV.max', 'epi.dt', 'epi.import')%in%pipeline.args[, stat] ) )
	#
	infile.ind	<- paste(indir, '/', infile.ind, sep='')
	infile.trm	<- paste(indir, '/', infile.trm, sep='')
	outfile.ind	<- paste(outdir, '/', outfile.ind, sep='')
	outfile.trm	<- paste(outdir, '/', outfile.trm, sep='')	
	#	set seed
	set.seed( pipeline.args['s.seed',][, as.numeric(v)] )
	#
	df.trm	<- as.data.table(read.csv(infile.trm, stringsAsFactors=FALSE, sep=' ', dec='.'))
	setnames(df.trm, c("IdInfector","IdInfected","TimeOfInfection","IsInfectorAcute"), c('IDTR','IDREC','TIME_TR','TR_ACUTE'))		
	#	transmissions happen either at baseline, or at unique times.
	#	the epi simulation allocates transmissions in 1/48 of a year, so draw a uniform number if there are more transmission per TIME_TR
	df.trm	<- df.trm[, {
				z<- TIME_TR
				if(TIME_TR>pipeline.args['yr.start',][, as.numeric(v)] & length(IDTR)>1)
					z<- sort(runif(length(IDTR), z, min(pipeline.args['yr.end',][, as.numeric(v)], z+pipeline.args['epi.dt',][, as.numeric(v)])))
				list(IDTR=IDTR, IDREC=IDREC, TIME_TR.new=z, TR_ACUTE=TR_ACUTE, l=length(IDTR))
			}, by='TIME_TR']
	df.trm[, TIME_TR:=NULL]
	setnames(df.trm, 'TIME_TR.new', 'TIME_TR')
	set(df.trm, NULL, 'YR', df.trm[, floor(TIME_TR)])
	#	check that all transmission times except baseline are unique
	tmp		<- subset(df.trm, TIME_TR>pipeline.args['yr.start',][, as.numeric(v)])
	stopifnot( nrow(tmp)==tmp[,length(unique(TIME_TR))] )
	
	
	df.ind	<- as.data.table(read.csv(infile.ind, stringsAsFactors=FALSE))		
	setnames(df.ind, c("Id","Gender","DoB","DateOfDeath","RiskGroup","Circumcised"), c('IDPOP','GENDER','DOB','DOD','RISK','CIRCM'))
	set(df.ind, df.ind[, which(CIRCM=='')], 'CIRCM', NA_character_)
	set(df.ind, NULL, 'CIRCM', df.ind[, factor(CIRCM)])
	set(df.ind, NULL, 'GENDER', df.ind[, factor(GENDER)])
	set(df.ind, NULL, 'RISK', df.ind[, factor(RISK)])	
	set(df.ind, df.ind[, which(DOD==-1)], 'DOD', pipeline.args['yr.end',][, as.numeric(v)]+1.)		
	tmp			<- subset(df.trm, select=c(IDREC, TIME_TR))
	setnames(tmp, 'IDREC','IDPOP')
	df.ind		<- merge(df.ind, tmp, by='IDPOP', all.x=TRUE)
	
	
	# compute prevalence and incidence by year	
	df.epi		<- df.trm[, list(INC=length(IDREC)), by='YR']
	tmp			<- df.epi[, 	{
				alive		<- which( floor(df.ind[['DOB']])<=YR  &  ceiling(df.ind[['DOD']])>YR )
				infected	<- which( floor(df.ind[['DOB']])<=YR  &  ceiling(df.ind[['DOD']])>YR  &  floor(df.ind[['TIME_TR']])<=YR )
				list(POP=length(alive), PREV=length(infected))				
			},by='YR']
	df.epi		<- merge( tmp, df.epi, by='YR' )	
	set(df.epi, NULL, 'PREVp', df.epi[, PREV/POP])	
	set(df.epi, NULL, 'INCp', df.epi[, INC/POP])			
	# 	SAMPLING PROBABILITIES and TOTALS PER YEAR
	#
	#	Can we detect a 25% or 50% reduction in HIV incidence in the most recent 2 or 3 years 
	#	with 1%, 5%, 10% of all recent incident cases sampled?
	#
	#	suppose exponentially increasing sampling over time
	#	the number of incident cases sampled is the total sampled in that year * the proportion of incident cases out of all non-sampled cases to date
	#	TODO this needs to be changed to fix the proportion of sequences sampled from incident
	df.sample	<- subset( df.epi, YR>= pipeline.args['yr.start',][, as.numeric(v)] & YR<pipeline.args['yr.end',][, as.numeric(v)] )
	#	exponential rate of increasing s.TOTAL (total sampling rate) per year
	tmp			<- log( 1+pipeline.args['s.PREV.max',][, as.numeric(v)]-pipeline.args['s.PREV.min',][, as.numeric(v)] ) / df.sample[, diff(range(YR))]
	tmp			<- df.sample[, exp( tmp*(YR-min(YR)) ) - 1 + pipeline.args['s.PREV.min',][, as.numeric(v)] ]
	set(df.sample, NULL, 's.CUMTOTAL', tmp)		
	set(df.sample, NULL, 's.n.CUMTOTAL', df.sample[, round(PREV*s.CUMTOTAL)])
	set(df.sample, NULL, 's.n.TOTAL', c(df.sample[1, s.n.CUMTOTAL], df.sample[, diff(s.n.CUMTOTAL)]))	
	set(df.sample, NULL, 's.n.INC', df.sample[, round(INC/(PREV-s.n.CUMTOTAL) * s.n.TOTAL)])
	set(df.sample, NULL, 's.n.notINC', df.sample[, round(s.n.TOTAL-s.n.INC)])	
	cat(paste('\n total number of sequences sampled=', df.sample[, sum( s.n.TOTAL )]))
	cat(paste('\n prop of sequences sampled among HIV+=', df.sample[, sum( s.n.TOTAL )] / df.sample[, rev(PREV)[1]]))		
	cat(paste('\n total number of incident sequences sampled=', df.sample[, sum( s.n.INC )]))
	cat(paste('\n total number of non-incident sequences sampled=', df.sample[, sum( s.n.notINC )]))	
	#	SAMPLE INFECTED INDIVIDUALS BASED ON NUMBERS PER YEAR
	#
	#	sample incident cases by year
	df.inds	<- copy(df.ind)
	tmp		<- df.trm[, {
				tmp<- df.sample[['s.n.INC']][ which(df.sample[['YR']]==YR) ]
				tmp<- sample(seq_along(IDREC), tmp)
				list( 	IDPOP=IDREC[tmp], TIME_TR=TIME_TR[tmp], 
						TIME_SEQ=TIME_TR[tmp]+rexp(length(tmp), rate=1/(3*30))/365, 
						INCIDENT_SEQ=rep('Y',length(tmp) ) )
			}, by='YR']
	df.inds	<- merge(df.inds, subset(tmp, select=c(IDPOP, TIME_SEQ, INCIDENT_SEQ)), by='IDPOP', all.x=1)			
	#	sample non-incident cases by year
	for(yr in df.sample[, YR][-1])
	{
		#	of all infected and not incident and not yet sampled, sample
		cat(paste('\nadd non-incident samples in year',yr))
		tmp		<- subset(df.inds, is.na(TIME_SEQ) & floor(DOB)<=yr & ceiling(DOD)>yr & floor(TIME_TR)<yr)
		cat(paste('\navailable non-sampled non-incident cases in year=',nrow(tmp)))
		tmp2	<- df.sample[['s.n.notINC']][ which(df.sample[['YR']]==yr) ]
		tmp2	<- sample(seq_len(nrow(tmp)), tmp2)
		#	set variables in df.inds
		tmp		<- data.table(IDPOP= tmp[tmp2, IDPOP], TIME_SEQ=runif(length(tmp2), min=yr, max=yr+1), INCIDENT_SEQ=rep('N',length(tmp2) ))
		cat(paste('\nsampled non-incident cases in year=',nrow(tmp)))
		tmp2	<- sapply(tmp[,IDPOP], function(x) df.inds[,which(IDPOP==x)])
		set(df.inds, tmp2, 'TIME_SEQ', tmp[,TIME_SEQ])
		set(df.inds, tmp2, 'INCIDENT_SEQ', tmp[,INCIDENT_SEQ])		
	}
	cat(paste('\n total number of HIV+ in df.inds=', nrow(subset(df.inds, !is.na(TIME_TR)))))
	cat(paste('\n total number of sampled HIV+ in df.inds=', nrow(subset(df.inds, !is.na(TIME_TR) & !is.na(TIME_SEQ)))))
	#
	#	check that allocation OK
	#
	set(df.inds, NULL, 'TIME_SEQYR', df.inds[, floor(TIME_SEQ)])
	tmp	<- subset(df.inds, !is.na(TIME_SEQ))[, list(s.n.TOTAL=length(IDPOP)), by='TIME_SEQYR']
	setkey(tmp, TIME_SEQYR)
	set(tmp,NULL,'s.n.CUMTOTAL',tmp[, cumsum(s.n.TOTAL)])
	stopifnot(  tmp[,tail(s.n.CUMTOTAL,1)]==df.sample[, tail(s.n.CUMTOTAL,1)] ) 
	#	set sampling in df.trm
	tmp		<- subset( df.inds, !is.na(TIME_SEQ), select=c(IDPOP, TIME_SEQ) )
	setnames(tmp, c('IDPOP','TIME_SEQ'), c('IDREC','SAMPLED_REC'))
	df.trms	<- merge(df.trm, tmp, by='IDREC', all.x=TRUE)
	setnames(tmp, c('IDREC','SAMPLED_REC'), c('IDTR','SAMPLED_TR'))
	df.trms	<- merge(df.trms, tmp, by='IDTR', all.x=TRUE)	
	#
	#	TRANSMISSION NETWORKS
	#
	require(igraph)
	#	need a unique index number for every cluster
	setkey(df.trms, TIME_TR)
	tmp			<- df.trms[, which(IDTR<0)]
	set(df.trms, tmp, 'IDTR', rev(-seq_along(tmp)))
	#	cluster with index case
	tmp			<- subset(df.trms, select=c(IDTR, IDREC))			
	tmp			<- graph.data.frame(tmp, directed=TRUE, vertices=NULL)
	tmp			<- data.table(IDPOP=as.integer(V(tmp)$name), CLU=clusters(tmp, mode="weak")$membership)
	tmp2		<- tmp[, list(CLU_SIZE=length(IDPOP)), by='CLU']
	setkey(tmp2, CLU_SIZE)
	tmp2[, IDCLU:=rev(seq_len(nrow(tmp2)))]
	tmp			<- subset( merge(tmp, tmp2, by='CLU'), select=c(IDPOP, IDCLU) )
	df.inds		<- merge( df.inds, tmp, by='IDPOP', all.x=TRUE )
	setnames(tmp, 'IDPOP', 'IDREC')
	df.trms		<- merge( df.trms, tmp, by='IDREC', all.x=TRUE )
	stopifnot( nrow(subset(df.trms, is.na(IDCLU)))==0 )
	cat(paste('\nFound transmission clusters, n=', df.trms[, length(unique(IDCLU))]))
	#
	#	PLOTS
	#
	if(with.plot)
	{
		require(ggplot2)
		require(reshape2)
		#	plot numbers sampled, prevalent, incident
		set(df.sample, NULL, 's.n.INC', df.sample[, as.integer(s.n.INC)])
		set(df.sample, NULL, 's.n.notINC', df.sample[, as.integer(s.n.notINC)])
		set(df.sample, NULL, 's.n.TOTAL', df.sample[, as.integer(s.n.TOTAL)])
		tmp	<- data.table(	stat=c('POP','PREV','INC','s.n.TOTAL','s.n.INC','s.n.notINC'), 
				stat.long=c('population size','HIV infected', 'HIV incident', 'Total\nsequenced', 'Total\nincident\nsequenced', 'Total\nnon-incident\nsequenced'))
		tmp	<- merge(	melt(df.sample, id.vars='YR', measure.vars=c('POP','PREV','INC','s.n.TOTAL','s.n.INC','s.n.notINC'), variable.name='stat', value.name='v'),
				tmp, by='stat' )
		ggplot(tmp, aes(x=YR, y=v, group=stat.long)) + geom_point() +
				scale_x_continuous(name='year', breaks=seq(1980,pipeline.args['yr.end',][, as.numeric(v)],2)) + scale_y_continuous(name='total')	+
				facet_grid(stat.long ~ ., scales='free_y', margins=FALSE)
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_Totals.pdf',sep='')
		ggsave(file=file, w=16, h=8)
		#	plot distribution between transmission time and sequencing time
		tmp	<- subset(df.inds, !is.na(TIME_SEQ))
		set(tmp, NULL, 'TIME_TO_SEQ', tmp[, TIME_SEQ-TIME_TR])
		ggplot(tmp, aes(x=TIME_TO_SEQ)) + geom_histogram(binwidth=1) + 
				scale_x_continuous(name='time from transmission to sequence sampling\n(years)', breaks=seq(0,100,2))
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_Time2Seq.pdf',sep='')
		ggsave(file=file, w=8, h=8)
		#	plot transmission network
		file		<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_TrNetworks.pdf',sep='')
		pdf(file=file, w=20, h=20)
		dummy	<- sapply( df.inds[, sort(na.omit(unique(IDCLU)))], function(clu)
				{
					cat(paste('\nprocess cluster no',clu))
					tmp					<- subset(df.inds, IDCLU==clu, select=c(IDPOP, GENDER, TIME_SEQ))
					tmp[, IS_SEQ:= tmp[, factor(!is.na(TIME_SEQ), label=c('N','Y'), levels=c(FALSE, TRUE))]]
					clu.igr				<- graph.data.frame(subset(df.trms, IDCLU==clu & IDTR>=0, select=c(IDTR, IDREC)), directed=TRUE, vertices=subset(tmp, select=c(IDPOP, GENDER, IS_SEQ)))
					V(clu.igr)$color	<- ifelse( get.vertex.attribute(clu.igr, 'IS_SEQ')=='Y', 'green', 'grey90' )
					V(clu.igr)$shape	<- ifelse( get.vertex.attribute(clu.igr, 'GENDER')=='M', 'circle', 'square' )
					
					par(mai=c(0,0,1,0))
					plot(clu.igr, main=paste('IDCLU=',clu,sep=''), vertex.size=2, vertex.label.cex=0.25, edge.arrow.size=0.5, layout=layout.fruchterman.reingold(clu.igr, niter=1e3) )
					legend('bottomright', fill=c('green','grey90'), legend=c('sequence sampled','sequence not sampled'), bty='n')
					legend('bottomleft', legend=c('square= Female','circle= Male'), bty='n')				
				})
		dev.off()	
	}	
	#
	#	SAVE SAMPLED RECIPIENTS AND TRANSMISSIONS TO SAMPLED RECIPIENTS
	#
	#	save for us
	file		<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'SAVE.R',sep='')
	save(file=file, df.epi, df.trms, df.inds, df.sample)
	#	save for virus tree simulator
	#	exclude columns that are not needed	
	df.inds	<- subset(df.inds, !is.na(TIME_TR))
	df.inds[, RISK:=NULL]
	df.inds[, INCIDENT_SEQ:=NULL]
	df.inds[, TIME_SEQYR:=NULL]	
	df.trms[, TR_ACUTE:=NULL]
	df.trms[, YR:=NULL]	
	cat(paste('\nwrite to file',outfile.ind))
	write.csv(file=outfile.ind, df.inds)
	cat(paste('\nwrite to file',outfile.trm))
	write.csv(file=outfile.trm, df.trms)
}

##--------------------------------------------------------------------------------------------------------
##	olli originally written 08-09-2014
##--------------------------------------------------------------------------------------------------------
#' @title HPTN071 parser (version 2, includes simulation of imports)
#' @description Reads files from the epi simulator in directory \code{indir} and writes csv files
#' in directory \code{outdir} for the virus tree simulator. The program samples sequences according to
#' an exponentially increasing sampling fraction in the same way as \code{prog.HPTN071.input.parser.v1}.
#' In addition, transmissions are broken and treated as imported from outside the simulated population.
#' The infected of a broken transmission chain is considered a new index case of a transmission chain within the 
#' simulated population. All input arguments are specified via the \code{argv} 
#' string, see the Examples.
#' @return NULL. Saves simulations to file.
#' @example example/ex.seq.sampler.v2.R
prog.HPTN071.input.parser.v2<- function()	
{
	require(data.table)
	verbose			<- 1
	with.plot		<- 1	
	pipeline.args	<- NULL
	indir			<- system.file(package="rPANGEAHIVsim", "misc")
	indir			<- ifelse(indir=='','/Users/Oliver/git/HPTN071sim/raw_trchain',indir)
	outdir			<- '/Users/Oliver/git/HPTN071sim/tmp140908'
	infile.ind		<- '140716_RUN001_IND.csv'
	infile.trm		<- '140716_RUN001_TRM.csv'
	infile.args		<- NA
	outfile.ind		<- '140716_RUN001_IND.csv'
	outfile.trm		<- '140716_RUN001_TRM.csv'
	#
	if(exists("argv"))
	{
		#	args input
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.ind= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.ind<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.trm= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.trm<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									infile.args= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.args<- tmp[1]
		#	args output
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									outfile.ind= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) outfile.ind<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									outfile.trm= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) outfile.trm<- tmp[1]		
	}
	if(verbose)
	{
		cat('\ninput args\n',paste(indir, infile.ind, infile.trm, outdir, outfile.ind, outfile.trm, sep='\n'))
	}	
	if(!is.na(infile.args))
	{
		load(infile.args)	#expect 'pipeline.args'
	}
	if(is.null(pipeline.args))
	{
		cat('\nCould not find pipeline.args, generating default')
		pipeline.args	<- rPANGEAHIVsim.pipeline.args()
	}
	stopifnot( all( c('yr.start', 'yr.end', 's.seed', 's.PREV.min', 's.PREV.max', 'epi.dt', 'epi.import')%in%pipeline.args[, stat] ) )
	#
	infile.ind	<- paste(indir, '/', infile.ind, sep='')
	infile.trm	<- paste(indir, '/', infile.trm, sep='')
	outfile.ind	<- paste(outdir, '/', outfile.ind, sep='')
	outfile.trm	<- paste(outdir, '/', outfile.trm, sep='')
	
	#	set seed
	set.seed( pipeline.args['s.seed',][, as.numeric(v)] )
	#
	#	prepare transmissions
	#	
	df.trm	<- as.data.table(read.csv(infile.trm, stringsAsFactors=FALSE, sep='', dec='.'))
	setnames(df.trm, c("IdInfector","IdInfected","TimeOfInfection","IsInfectorAcute"), c('IDTR','IDREC','TIME_TR','TR_ACUTE'))
	stopifnot( df.trm[, !any(is.na(TR_ACUTE))] )
	set(df.trm, df.trm[, which(TR_ACUTE<0)], 'TR_ACUTE', NA_integer_)
	set(df.trm, NULL, 'TR_ACUTE', df.trm[, factor(TR_ACUTE, levels=c(0,1), labels=c('No','Yes'))])
	#	transmissions happen either at baseline, or at unique times.
	#	the epi simulation allocates transmissions in 1/48 of a year, so draw a uniform number if there are more transmission per TIME_TR
	tmp		<- df.trm[, range(TIME_TR)]
	df.trm	<- df.trm[, {
				z<- TIME_TR
				if(TIME_TR>tmp[1] & length(IDTR)>1)
					z<- sort(runif(length(IDTR), z, min(ceiling(tmp[2]), z+pipeline.args['epi.dt',][, as.numeric(v)])))
				list(IDTR=IDTR, IDREC=IDREC, TIME_TR.new=z, TR_ACUTE=TR_ACUTE, l=length(IDTR))
			}, by='TIME_TR']
	df.trm[, TIME_TR:=NULL]
	setnames(df.trm, 'TIME_TR.new', 'TIME_TR')	
	#	set baseline cases as negative ID
	tmp		<- df.trm[, which(IDTR=='-1')]
	cat(paste('\nFound index cases, n=', length(tmp)))
	set(df.trm, tmp, 'IDTR', rev(-seq_along(tmp)))
	#	check that all transmission times except baseline are unique
	tmp		<- subset(df.trm, TIME_TR>min(TIME_TR))
	stopifnot( nrow(tmp)==tmp[,length(unique(TIME_TR))] )	
	cat(paste('\nFound transmissions, n=', nrow(df.trm)))
	cat(paste('\nTotal transmitters, n=', df.trm[, length(unique(IDTR))]))
	cat(paste('\nTotal index cases, n=', df.trm[, length(which(unique(IDTR)<0))]))
	#
	#	prepare patient metavariables
	#
	df.ind	<- as.data.table(read.csv(infile.ind, stringsAsFactors=FALSE))		
	setnames(df.ind, c("Id","Gender","DoB","DateOfDeath","RiskGroup","T1_CD4350","Circumcised"), c('IDPOP','GENDER','DOB','DOD','RISK','T1_CD4350','CIRCM'))
	set(df.ind, df.ind[, which(CIRCM=='')], 'CIRCM', NA_character_)
	set(df.ind, NULL, 'CIRCM', df.ind[, factor(CIRCM)])
	set(df.ind, NULL, 'GENDER', df.ind[, factor(GENDER)])
	set(df.ind, NULL, 'RISK', df.ind[, factor(RISK)])	
	set(df.ind, df.ind[, which(DOD==-1)], 'DOD', ceiling(max(df.trm$TIME_TR))+1.)
	set(df.ind, df.ind[, which(T1_CD4350==-10)],'T1_CD4350',NA_real_)
	set(df.ind, df.ind[, which(T1_CD4350%in%c(-8,-9,-11,-12))],'T1_CD4350',Inf)	
	tmp			<- subset(df.trm, select=c(IDREC, TIME_TR))
	stopifnot( df.ind[, !any(is.na(DOB))] )
	setnames(tmp, 'IDREC','IDPOP')
	df.ind		<- merge(df.ind, tmp, by='IDPOP', all.x=TRUE)
	cat(paste('\nFound individuals with a valid record, n=', nrow(df.ind)))
	cat(paste('\nFound individuals with an infection event, n=', nrow(subset(df.ind,!is.na(TIME_TR)))))
	#	add time of infection of transmitter to df.trm	
	tmp		<- subset(df.ind, select=c(IDPOP, TIME_TR))
	setnames(tmp, c('IDPOP','TIME_TR'), c('IDTR','IDTR_TIME_INFECTED') )
	setkey(tmp, IDTR)
	df.trm	<- merge(df.trm, unique(tmp), by='IDTR', all.x=TRUE)		
	#	simulate time individual ready for sequencing
	#	ignore T1_CD4350 for now
	df.ind[, T1_CD4350:=NULL]
	df.ind[, T1_SEQ:= df.ind[, rgamma(nrow(df.ind),shape=9,scale=0.25 ) + TIME_TR]]
	#
	#	reduce to time frame of interest
	#
	#tmp		<- subset( df.trm, TIME_TR>=as.numeric( pipeline.args['yr.end',][, as.numeric(v)] ) )[, IDREC]
	df.trm	<- subset( df.trm, TIME_TR<as.numeric( pipeline.args['yr.end',][, as.numeric(v)] ) )
	#df.ind	<- subset(df.ind, !IDPOP%in%tmp)
	df.ind	<- subset(df.ind, is.na(DOB) | DOB<pipeline.args['yr.end',][, as.numeric(v)] )
	df.ind	<- subset(df.ind, is.na(DOD) | DOD >= floor(min(df.trm$TIME_TR)) )
	cat(paste('\nFound individuals born before',pipeline.args['yr.end',][, as.numeric(v)],', n=', nrow(df.ind)))
	cat(paste('\nFound transmissions before',pipeline.args['yr.end',][, as.numeric(v)],', n=', nrow(df.trm)))
	cat(paste('\nTotal transmitters, n=', df.trm[, length(unique(IDTR))]))		
	stopifnot( length(setdiff( subset(df.trm, IDTR>0)[, IDTR], df.ind[, IDPOP] ))==0 )
	stopifnot( length(setdiff( df.trm[, IDREC], df.ind[, IDPOP] ))==0 )
	#	simulate a fraction of transmissions to be imports
	tmp		<- PANGEA.ImportSimulator.SimulateIndexCase(df.ind, df.trm, epi.import= pipeline.args['epi.import',][,as.numeric(v)])
	df.trm	<- tmp$df.trm
	df.ind	<- tmp$df.ind
	#
	#	set infection times for index case
	#
	#	add IDTR_TIME_INFECTED for baseline cases
	tmp		<- df.trm[, which(is.na(IDTR_TIME_INFECTED))]	
	set( df.trm, tmp, 'IDTR_TIME_INFECTED', df.trm[tmp, runif(length(tmp), TIME_TR-5, TIME_TR)] )	
	tmp		<- PANGEA.ImportSimulator.SimulateStartingTimeOfIndexCase(df.ind, df.trm, index.starttime.mode= pipeline.args['index.starttime.mode',][,v])
	df.trm	<- tmp$df.trm
	df.ind	<- tmp$df.ind	
	#
	#	sample sequences 
	#
	PANGEA.Seqsampler(df.ind, df.trm, pipeline.args, outfile.ind, outfile.trm, with.plot=with.plot)
	#
	return(1)
}
##--------------------------------------------------------------------------------------------------------
##	olli originally written 23-10-2014
##--------------------------------------------------------------------------------------------------------
#' @title HPTN071 parser (version 3, uses CD4 counts)
#' @description Reads files from the epi simulator in directory \code{indir} and writes csv files
#' in directory \code{outdir} for the virus tree simulator. The program samples sequences according to
#' an exponentially increasing sampling fraction in the same way as \code{prog.HPTN071.input.parser.v1}.
#' In addition, transmissions are broken and treated as imported from outside the simulated population.
#' The infected of a broken transmission chain is considered a new index case of a transmission chain within the 
#' simulated population. All input arguments are specified via the \code{argv} 
#' string, see the Examples.
#' @return NULL. Saves simulations to file.
#' @example example/ex.seq.sampler.v2.R
prog.HPTN071.input.parser.v3<- function()	
{
	require(data.table)
	verbose			<- 1
	with.plot		<- 1	
	pipeline.args	<- NULL
	indir			<- system.file(package="rPANGEAHIVsim", "misc")
	indir			<- ifelse(indir=='','/Users/Oliver/git/HPTN071sim/raw_trchain',indir)
	outdir			<- '/Users/Oliver/git/HPTN071sim/tmp140908'
	infile.ind		<- '140716_RUN001_IND.csv'
	infile.trm		<- '140716_RUN001_TRM.csv'
	infile.args		<- NA
	outfile.ind		<- '140716_RUN001_IND.csv'
	outfile.trm		<- '140716_RUN001_TRM.csv'
	#
	if(exists("argv"))
	{
		#	args input
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.ind= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.ind<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.trm= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.trm<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									infile.args= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.args<- tmp[1]
		#	args output
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									outfile.ind= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) outfile.ind<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									outfile.trm= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) outfile.trm<- tmp[1]		
	}
	if(verbose)
	{
		cat('\ninput args\n',paste(indir, infile.ind, infile.trm, outdir, outfile.ind, outfile.trm, sep='\n'))
	}	
	if(!is.na(infile.args))
	{
		load(infile.args)	#expect 'pipeline.args'
	}
	if(is.null(pipeline.args))
	{
		cat('\nCould not find pipeline.args, generating default')
		pipeline.args	<- rPANGEAHIVsim.pipeline.args()
	}
	stopifnot( all( c('yr.start', 'yr.end', 's.seed', 's.PREV.min', 's.PREV.max', 'epi.dt', 'epi.import')%in%pipeline.args[, stat] ) )
	#
	infile.ind	<- paste(indir, '/', infile.ind, sep='')
	infile.trm	<- paste(indir, '/', infile.trm, sep='')
	outfile.ind	<- paste(outdir, '/', outfile.ind, sep='')
	outfile.trm	<- paste(outdir, '/', outfile.trm, sep='')
	
	#	set seed
	set.seed( pipeline.args['s.seed',][, as.numeric(v)] )
	#
	#	prepare transmissions
	#	
	df.trm	<- as.data.table(read.csv(infile.trm, stringsAsFactors=FALSE, sep='', dec='.'))
	setnames(df.trm, c("IdInfector","IdInfected","TimeOfInfection","IsInfectorAcute"), c('IDTR','IDREC','TIME_TR','TR_ACUTE'))
	stopifnot( df.trm[, !any(is.na(TR_ACUTE))] )
	set(df.trm, df.trm[, which(TR_ACUTE<0)], 'TR_ACUTE', NA_integer_)
	set(df.trm, NULL, 'TR_ACUTE', df.trm[, factor(TR_ACUTE, levels=c(0,1), labels=c('No','Yes'))])
	#	transmissions happen either at baseline, or at unique times.
	#	the epi simulation allocates transmissions in 1/48 of a year, so draw a uniform number if there are more transmission per TIME_TR
	tmp		<- df.trm[, range(TIME_TR)]
	df.trm	<- df.trm[, {
				z<- TIME_TR
				if(length(IDTR)>1)
					z<- sort(runif(length(IDTR), z, min(ceiling(tmp[2]), z+pipeline.args['epi.dt',][, as.numeric(v)])))
				list(IDTR=IDTR, IDREC=IDREC, TIME_TR.new=z, TR_ACUTE=TR_ACUTE, l=length(IDTR))
			}, by='TIME_TR']
	df.trm[, TIME_TR:=NULL]
	setnames(df.trm, 'TIME_TR.new', 'TIME_TR')	
	#	stop if not all transmission times are unique, except at the end which does not matter
	stopifnot( subset(df.trm, TIME_TR<max(TIME_TR))[, length(unique(TIME_TR))==length(TIME_TR)] )
	#	set baseline cases as negative ID
	tmp		<- df.trm[, which(IDTR=='-1')]
	cat(paste('\nFound index cases, n=', length(tmp)))
	set(df.trm, tmp, 'IDTR', rev(-seq_along(tmp)))
	cat(paste('\nFound transmissions, n=', nrow(df.trm)))
	cat(paste('\nTotal transmitters, n=', df.trm[, length(unique(IDTR))]))
	cat(paste('\nTotal index cases, n=', df.trm[, length(which(unique(IDTR)<0))]))
	#
	#	prepare patient metavariables
	#
	df.ind	<- as.data.table(read.csv(infile.ind, stringsAsFactors=FALSE))		
	setnames(df.ind, c(	"Id","Gender","DoB","DateOfDeath","RiskGroup","Circumcised"), c('IDPOP','GENDER','DOB','DOD','RISK','CIRCM'))
	setnames(df.ind, c(	"t.cd4below350","t.cd4.500","t.cd4.350","t.cd4.200","t.aidsdeath","SPVL.category","t.sc"), c('T1_CD4_l350','T1_CD4_500','T1_CD4_350','T1_CD4_200',"DOAD","SPVL","DOSC"))		
	set(df.ind, df.ind[, which(CIRCM=='')], 'CIRCM', NA_character_)
	set(df.ind, NULL, 'CIRCM', df.ind[, factor(CIRCM)])
	set(df.ind, NULL, 'GENDER', df.ind[, factor(GENDER)])
	set(df.ind, NULL, 'RISK', df.ind[, factor(RISK)])	
	set(df.ind, df.ind[, which(DOD==-1)], 'DOD', ceiling(max(df.trm$TIME_TR))+1.)
	set(df.ind, df.ind[, which(DOAD==-10)], 'DOAD', ceiling(max(df.trm$TIME_TR))+1.)
	set(df.ind, df.ind[, which(DOSC==-10)], 'DOSC', NA_real_)
	set(df.ind, df.ind[, which(SPVL==-1)], 'SPVL', NA_integer_)
	set(df.ind, NULL, 'SPVL', df.ind[, factor(SPVL, levels=c(0, 1, 2, 3), labels=c('le40','le45','le50','g50'))])	
	set(df.ind, df.ind[, which(T1_CD4_500%in%c(-1,-10))],'T1_CD4_500',NA_real_)
	stopifnot( df.ind[, !any(T1_CD4_500<0, na.rm=TRUE)])
	set(df.ind, df.ind[, which(T1_CD4_350==-10)],'T1_CD4_350',NA_real_)
	stopifnot( df.ind[, !any(T1_CD4_350<0, na.rm=TRUE)] )
	set(df.ind, df.ind[, which(T1_CD4_200==-10)],'T1_CD4_200',NA_real_)
	stopifnot( df.ind[, !any(T1_CD4_200<0, na.rm=TRUE)] )
	set(df.ind, df.ind[, which(T1_CD4_l350==-10)],'T1_CD4_l350',NA_real_)
	set(df.ind, df.ind[, which(T1_CD4_l350%in%c(-8,-9,-11,-12))],'T1_CD4_l350',Inf)	
	stopifnot( df.ind[, !any(T1_CD4_l350<0, na.rm=TRUE)] )
	stopifnot( df.ind[, !any(DOD>DOAD+0.125)] )
	set(df.ind, NULL, c('T1_CD4_l350','DOAD'), NULL)
	#	add transmission time
	tmp			<- subset(df.trm, select=c(IDREC, TIME_TR))
	stopifnot( df.ind[, !any(is.na(DOB))] )
	setnames(tmp, 'IDREC','IDPOP')
	df.ind		<- merge(df.ind, tmp, by='IDPOP', all.x=TRUE)
	cat(paste('\nFound individuals with a valid record, n=', nrow(df.ind)))
	cat(paste('\nFound individuals with an infection event, n=', nrow(subset(df.ind,!is.na(TIME_TR)))))
	#	reset DOSC if needed, because TIME_TR got randomized by a small bit above
	tmp			<- df.ind[, which(DOSC<TIME_TR)]
	set(df.ind, tmp, 'DOSC', df.ind[tmp, TIME_TR+(TIME_TR-DOSC)])
	tmp			<- df.ind[, which(T1_CD4_500<TIME_TR)]
	set(df.ind, tmp, 'T1_CD4_500', df.ind[tmp, TIME_TR+(TIME_TR-T1_CD4_500)])
	tmp			<- df.ind[, which(T1_CD4_350<TIME_TR)]
	set(df.ind, tmp, 'T1_CD4_350', df.ind[tmp, TIME_TR+(TIME_TR-T1_CD4_350)])
	tmp			<- df.ind[, which(T1_CD4_200<TIME_TR)]
	set(df.ind, tmp, 'T1_CD4_200', df.ind[tmp, TIME_TR+(TIME_TR-T1_CD4_200)])
	#	reset DOD if needed, because TIME_TR got randomized by a small bit above
	tmp			<- df.ind[, which(DOD<TIME_TR)]
	set(df.ind, tmp, 'DOD', df.ind[tmp, TIME_TR+(TIME_TR-DOD)])
	#	set T1_CD4 if starting with lower count than 500
	tmp			<- df.ind[, which(!is.na(TIME_TR) & !is.na(T1_CD4_200) & is.na(T1_CD4_350))]
	set(df.ind, tmp, 'T1_CD4_350', df.ind[tmp, TIME_TR])
	set(df.ind, tmp, 'T1_CD4_500', df.ind[tmp, TIME_TR-0.01])
	tmp			<- df.ind[, which(!is.na(TIME_TR) & !is.na(T1_CD4_350) & is.na(T1_CD4_500))]
	set(df.ind, tmp, 'T1_CD4_500', df.ind[tmp, TIME_TR])
	#	set T1_CD4 if ending with lower count than 200
	tmp			<- df.ind[, which(!is.na(TIME_TR) & !is.na(T1_CD4_350) & is.na(T1_CD4_200))]
	set(df.ind, tmp, 'T1_CD4_200', df.ind[tmp, DOD])
	tmp			<- df.ind[, which(!is.na(TIME_TR) & !is.na(T1_CD4_500) & is.na(T1_CD4_350))]
	set(df.ind, tmp, 'T1_CD4_350', df.ind[tmp, DOD])
	set(df.ind, tmp, 'T1_CD4_200', df.ind[tmp, DOD+0.01])
	tmp			<- df.ind[, which(!is.na(TIME_TR) & is.na(T1_CD4_500))]
	set(df.ind, tmp, 'T1_CD4_500', df.ind[tmp, DOD])
	set(df.ind, tmp, 'T1_CD4_350', df.ind[tmp, DOD+0.01])
	set(df.ind, tmp, 'T1_CD4_200', df.ind[tmp, DOD+0.02])
	#	T1_CD4 is only after acute phase, and before we interpolate from 1500. 
	#	add time of infection of transmitter to df.trm	
	tmp		<- subset(df.ind, select=c(IDPOP, TIME_TR))
	setnames(tmp, c('IDPOP','TIME_TR'), c('IDTR','IDTR_TIME_INFECTED') )
	setkey(tmp, IDTR)
	df.trm	<- merge(df.trm, unique(tmp), by='IDTR', all.x=TRUE)
	stopifnot( df.trm[, !any(TIME_TR<=IDTR_TIME_INFECTED, na.rm=TRUE)] )
	#	simulate time individual ready for sequencing
	df.ind	<- PANGEA.Seqsampler.SimulateGuideToSamplingTimes(df.ind, seqtime.mode= pipeline.args['seqtime.mode',][,v])
	#
	#
	#	reduce to time frame of interest
	#
	#tmp		<- subset( df.trm, TIME_TR>=as.numeric( pipeline.args['yr.end',][, as.numeric(v)] ) )[, IDREC]
	df.trm	<- subset( df.trm, TIME_TR<as.numeric( pipeline.args['yr.end',][, as.numeric(v)] ) )
	#df.ind	<- subset(df.ind, !IDPOP%in%tmp)
	df.ind	<- subset(df.ind, is.na(DOB) | DOB<pipeline.args['yr.end',][, as.numeric(v)] )
	df.ind	<- subset(df.ind, is.na(DOD) | DOD >= floor(min(df.trm$TIME_TR)) )
	cat(paste('\nFound individuals born before',pipeline.args['yr.end',][, as.numeric(v)],', n=', nrow(df.ind)))
	cat(paste('\nFound transmissions before',pipeline.args['yr.end',][, as.numeric(v)],', n=', nrow(df.trm)))
	cat(paste('\nTotal transmitters, n=', df.trm[, length(unique(IDTR))]))		
	stopifnot( length(setdiff( subset(df.trm, IDTR>0)[, IDTR], df.ind[, IDPOP] ))==0 )
	stopifnot( length(setdiff( df.trm[, IDREC], df.ind[, IDPOP] ))==0 )
	#	simulate a fraction of transmissions to be imports
	tmp		<- PANGEA.ImportSimulator.SimulateIndexCase(df.ind, df.trm, epi.import= pipeline.args['epi.import',][,as.numeric(v)])
	df.trm	<- tmp$df.trm
	df.ind	<- tmp$df.ind
	#
	#	set infection times for index case
	#
	#	add IDTR_TIME_INFECTED for baseline cases
	tmp		<- df.trm[, which(is.na(IDTR_TIME_INFECTED))]	
	set( df.trm, tmp, 'IDTR_TIME_INFECTED', df.trm[tmp, runif(length(tmp), TIME_TR-5, TIME_TR)] )	
	tmp		<- PANGEA.ImportSimulator.SimulateStartingTimeOfIndexCase(df.ind, df.trm, index.starttime.mode= pipeline.args['index.starttime.mode',][,v])
	df.trm	<- tmp$df.trm
	df.ind	<- tmp$df.ind	
	#
	#	sample sequences 
	#
	PANGEA.Seqsampler.v3(df.ind, df.trm, pipeline.args, outfile.ind, outfile.trm, with.plot=with.plot)
	#
	return(1)
}
##--------------------------------------------------------------------------------------------------------
##	program to generate files for Seq Gen from output of Matt s VirusTreeSimulator
##	olli originally written 17-08-2014
##--------------------------------------------------------------------------------------------------------
prog.PANGEA.SeqGen.createInputFile.v2<- function()
{
	stop()
	verbose			<- 1
	with.plot		<- 1
	label.sep		<- '|'	
	#
	#	read I/O
	#
	indir.epi		<- '/Users/Oliver/git/HPTN071sim/tmp140914/TrChains'
	infile.epi		<- '140716_RUN001_SAVE.R'	
	indir.vts		<- '/Users/Oliver/git/HPTN071sim/tmp140914/VirusTreeSimulator'
	infile.prefix	<- '140716_RUN001_'	
	infile.args		<- '/Users/Oliver/git/HPTN071sim/tmp140914/140716_RUN001_PipeArgs.R'
	outdir.sg		<- '/Users/Oliver/git/HPTN071sim/tmp140914/SeqGen'	
	if(exists("argv"))
	{
		#	args input
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									indir.epi= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir.epi<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.epi= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.epi<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									indir.vts= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir.vts<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.vts= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.prefix<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									infile.args= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.args<- tmp[1]		
		#	args output
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir.sg<- tmp[1]		
		
	}
	if(verbose)
	{
		cat('\ninput args\n',paste(indir.epi, infile.epi, indir.vts, infile.prefix, outdir.sg, sep='\n'))
	}	
	if(!is.na(infile.args))
	{
		load(infile.args)	#expect 'pipeline.args'
	}
	if(is.null(pipeline.args))
	{
		cat('\nCould not find pipeline.args, generating default')
		pipeline.args	<- rPANGEAHIVsim.pipeline.args()
	}	
	stopifnot( all( c('s.seed','wher.mu','wher.sigma','bwerm.mu')%in%pipeline.args[, stat] ) )
	#
	#	setup samplers
	#
	cat(paste('\ncreate sampler of evolutionary rates'))
	#	create sampler of within host evolutionary rates
	rER.pol			<- PANGEA.WithinHostEvolutionaryRate.create.sampler.v1(wher.mu=pipeline.args['wher.mu',][, as.numeric(v)], wher.sigma=pipeline.args['wher.sigma',][, as.numeric(v)])
	#	create sampler of between host evolutionary rate dampener
	tmp				<- PANGEA.TransmissionEdgeEvolutionaryRate.create.sampler(er.shift=pipeline.args['bwerm.mu',][, as.numeric(v)])
	rERbw			<- tmp$rERbw
	rERbw.args		<- tmp$rERbw.args
	#	create sampler of ancestral sequences
	cat(paste('\ncreate sampler of ancestral sequences'))
	tmp				<- PANGEA.RootSeq.create.sampler(root.ctime.grace= 0.5, sample.grace= 3)
	rANCSEQ			<- tmp$rANCSEQ
	rANCSEQ.args	<- tmp$rANCSEQ.args 	
	#	read GTR parameters
	log.df			<- PANGEA.GTR.params( )
	if( pipeline.args['dbg.rER',][, as.numeric(v)] )
	{
		cat(paste('\nSetting mus to mean per gene and codon_pos'))
		tmp		<- log.df[, list(mu= mean(mu)), by=c('GENE','CODON_POS')]
		#tmp[, ER:=mu*log.df$meanRate[1]]
		log.df	<- merge( subset(log.df, select=which(colnames(log.df)!='mu')), tmp, by=c('GENE','CODON_POS'))		
	}	
	#
	#
	#
	file		<- paste(indir.epi, '/', infile.epi, sep='')
	load(file)	#expect "df.epi"    "df.trms"   "df.inds"   "df.sample"
	#	
	infiles		<- list.files(indir.vts)
	tmp			<- paste('^',infile.prefix,'.*nex$',sep='')
	infiles		<- infiles[ grepl(tmp, infiles)  ]	
	if(!length(infiles))	stop('cannot find files matching criteria')
	#
	set.seed( pipeline.args['s.seed',][, as.numeric(v)] )
	#
	#	read from VirusTreeSimulator and convert branch lengths in time to branch lengths in subst/site
	#
	df.ph			<- vector('list', length(infiles))
	df.nodestat		<- vector('list', length(infiles))	
	for(i in seq_along(infiles))
	{				
		#i<- 4
		infile			<- infiles[i]
		cat(paste('\nprocess file',i,infile))
		file			<- paste(indir.vts, '/', infile, sep='')
		#	read brl, units from annotated nexus file. attention: () may not contain two nodes
		tmp				<- hivc.beast2out.read.nexus.and.stats(file, method.node.stat='any.node')
		ph				<- tmp$tree
		node.stat		<- tmp$node.stat
		node.stat		<- subset(node.stat, STAT=='Unit')
		set(node.stat, NULL, 'VALUE', node.stat[, gsub('\"','',VALUE)])
		node.stat[, IDPOP:= as.integer(node.stat[,substr(VALUE, 4, nchar(VALUE))])]
		node.stat		<- merge(subset(df.inds, select=c(IDPOP, GENDER, DOB, TIME_SEQ, IDCLU)), subset(node.stat, select=c(IDPOP, NODE_ID)), by='IDPOP')	
		#
		#	create collapsed Newick tree with expected substitutions / site for each branch 
		#
		#	draw evolutionary rates for every individual in the transmission chain
		node.stat		<- merge(node.stat, data.table( IDPOP=node.stat[, unique(IDPOP)], ER= rER.pol(node.stat[, length(unique(IDPOP))]) ), by='IDPOP')
		#	smaller ER for transmission edge
		node.stat[, BWM:= rERbw(nrow(node.stat), rERbw.args)]
		set(node.stat, NULL, 'BWM', node.stat[, ER/BWM])	#re-set to previous notation
		#	set BWM to 1 for all non-transmission edges
		tmp				<- ph$edge
		colnames(tmp)	<- c('NODE_ID','NODE_ID_TO')
		tmp				<- merge( subset(node.stat, select=c(NODE_ID, IDPOP)), data.table(tmp, EDGE_ID=seq_len(nrow(tmp))), by='NODE_ID') 
		setnames(tmp, c('NODE_ID','IDPOP','NODE_ID_TO'), c('NODE_ID_FROM','IDPOP_FROM','NODE_ID'))
		tmp				<- merge( subset(node.stat, select=c(NODE_ID, IDPOP)), tmp, by='NODE_ID') 
		tmp				<- subset( tmp, IDPOP!=IDPOP_FROM)[, NODE_ID_FROM]	#ER slows down in transmitter leading to infection. want to slow down edge leading to NODE_FROM.
		set( node.stat, node.stat[, which(!NODE_ID%in%tmp)], 'BWM', 1. )
		#	no ER possible for root node - there s no edge leading to it
		set(node.stat, node.stat[, which(NODE_ID==Ntip(ph)+1)], c('ER','BWM'), NA_real_)		
		#	set root edge evolutionary rate to overall mean between-host rate
		#	get NODE_ID of edge from root
		tmp				<- ph$edge[match(Ntip(ph)+1, ph$edge[1, ]), 2]
		tmp				<- node.stat[, which(NODE_ID==tmp)]		
		set(node.stat, tmp, 'ER', log.df[1,meanRate] )		
		set(node.stat, tmp, 'BWM', 1. )		# no need to further slow down root edge
		#	check calendar time of root in simulated phylogeny for consistency
		tmp				<- seq.collapse.singles(ph)
		tmp2			<- regmatches(tmp$tip.label[1], regexpr('ID_[0-9]+',tmp$tip.label[1]))
		tmp2			<- as.numeric(substr(tmp2, 4, nchar(tmp2)))
		tmp2			<- subset(node.stat, IDPOP==tmp2)[1, TIME_SEQ]
		root.ctime		<- ifelse(Nnode(tmp), tmp2 - (node.depth.edgelength(tmp)[1] + tmp$root.edge), tmp2-tmp$root.edge)		
		tmp				<- subset(node.stat, IDPOP<0)[, unique(IDPOP)]
		stopifnot(length(tmp)==1)
		stopifnot(subset(df.trms, IDTR==tmp)[, round(IDTR_TIME_INFECTED, d=1)]==round(root.ctime, d=1))
		#	set expected numbers of substitutions per branch within individual IDPOP
		setkey(node.stat, NODE_ID)
		ph$edge.length	<- ph$edge.length * node.stat[ ph$edge[, 2], ][, ER / BWM]		 
		stopifnot(all(!is.na(ph$edge.length)))		
		#	once expected number of substitutions / site are simulated, can collapse singleton nodes
		ph				<- seq.collapse.singles(ph)	
		#	set tip label so that IDPOP can be checked for consistency	
		if(pipeline.args['epi.model',][,v]=='HPTN071')
			node.stat[, LABEL:= node.stat[, paste('IDPOP_',IDPOP,label.sep,GENDER,label.sep,'DOB_',round(DOB,d=3),label.sep,round(TIME_SEQ,d=3),sep='')]]
		if(pipeline.args['epi.model',][,v]=='DSPS')
			node.stat[, LABEL:= node.stat[, paste('IDPOP_',IDPOP,label.sep,GENDER,label.sep,'DOB_',NA,label.sep,round(TIME_SEQ,d=3),sep='')]]		
		setkey(node.stat, NODE_ID)
		ph$tip.label	<- node.stat[seq_len(Ntip(ph)), ][, LABEL]
		#
		df.nodestat[[i]]<- node.stat
		tmp				<- ifelse(Nnode(ph), write.tree(ph, digits = 10), paste( '(',ph$tip.label,':',ph$root.edge,',NOEXIST_NA|NA|DOB_NA|',root.ctime,':0):0;', sep='')	)
		df.ph[[i]]		<- data.table(ROOT_CALENDAR_TIME= root.ctime, IDCLU=node.stat[, unique(IDCLU)], NEWICK=tmp)
		#readline()
	}
	df.ph		<- do.call('rbind',df.ph)
	df.nodestat	<- do.call('rbind',df.nodestat)
	#	check that we have exactly one root edge with overall mean between host rate per cluster
	stopifnot( df.nodestat[, length(unique(IDCLU))]==nrow(subset(df.nodestat, ER==log.df$meanRate[1])) )
	
	if(with.plot)
	{
		#ggplot(subset(df.nodestat, ER!=log.df$meanRate[1]) , aes(x=ER/BWM)) + geom_histogram(binwidth=0.001)
		#ggplot(subset(df.nodestat, ER!=log.df$meanRate[1] & BWM!=1.) , aes(x=ER/BWM)) + geom_histogram(binwidth=0.0001)
		#ggplot(subset(df.nodestat, ER!=log.df$meanRate[1]), aes(x=ER, y=BWM)) + geom_point()	
		#	plot used within-host ERs
		ggplot(subset(df.nodestat, ER!=log.df$meanRate[1]), aes(x=ER/BWM)) + geom_histogram(binwidth=0.001)	+ labs(x='simulated within-host evolutionary rate') +
				scale_x_continuous(breaks= seq(0, 0.02, 0.002))
		file	<- paste(indir.epi, '/', substr(infile.epi,1,nchar(infile.epi)-6),'INFO_sg_ER.pdf', sep='')
		ggsave(file, w=6, h=6)
		#	plot used between host modifiers
		ggplot(subset(df.nodestat, ER!=log.df$meanRate[1] & BWM==1) , aes(x=ER/BWM)) + geom_histogram(binwidth=0.001) + labs(x='simulated within-host rate evolutionary rate\nwithout transmission edges') +
				scale_x_continuous(breaks= seq(0, 0.02, 0.002))
		file	<- paste(indir.epi, '/', substr(infile.epi,1,nchar(infile.epi)-6),'INFO_sg_BWM.pdf', sep='')
		ggsave(file, w=6, h=6)
		#	plot used ERs along transmission edges
		ggplot(subset(df.nodestat, ER!=log.df$meanRate[1] & BWM!=1) , aes(x=ER/BWM)) + geom_histogram(binwidth=0.0001) + labs(x='simulated evolutionary rates along transmission edges') +
				scale_x_continuous(breaks= seq(0, 0.02, 0.0005))
		file	<- paste(indir.epi, '/', substr(infile.epi,1,nchar(infile.epi)-6),'INFO_sg_BWER.pdf', sep='')
		ggsave(file, w=6, h=6)		
	}	
	#
	#	draw ancestral sequences and add to df.ph
	#
	root.ctime		<- df.ph[, ROOT_CALENDAR_TIME]
	ancseq			<- rANCSEQ(root.ctime, rANCSEQ.args)
	ancseq			<- data.table(ANCSEQ= apply(as.character(ancseq),1,function(x) paste(x, collapse='')) )		
	df.ph			<- cbind(df.ph, ancseq)
	#
	#	create SEQ-GEN input file
	#	
	partition.len	<- c( ncol(rANCSEQ.args$anc.seq.gag), ncol(rANCSEQ.args$anc.seq.pol), ncol(rANCSEQ.args$anc.seq.env) )
	#partition.er	<- c( 2.5, 4, 5 )
	#	split ancestral sequence into GENE and CODON_POS
	df.ph[, ANCSEQ.GAG:= substr(ANCSEQ, 1, partition.len[1])]
	df.ph[, ANCSEQ.POL:= substr(ANCSEQ, partition.len[1]+1, partition.len[1]+partition.len[2])]
	df.ph[, ANCSEQ.ENV:= substr(ANCSEQ, partition.len[1]+partition.len[2]+1, partition.len[1]+partition.len[2]+partition.len[3])]	
	stopifnot( all( strsplit( df.ph[1, ANCSEQ], '')[[1]] == strsplit( paste( df.ph[1, ANCSEQ.GAG], df.ph[1, ANCSEQ.POL], df.ph[1, ANCSEQ.ENV], sep=''), '' )[[1]] ) )
	df.ph[, ANCSEQ:=NULL]
	df.ph	<- melt(df.ph, id.var=c('ROOT_CALENDAR_TIME','IDCLU','NEWICK'), value.name='ANCSEQ', variable.name='GENE', variable.factor=FALSE)
	set(df.ph, NULL, 'GENE', df.ph[, substr(GENE, 8, nchar(GENE))])	
	df.ph	<- df.ph[,  {
				tmp	<- strsplit(ANCSEQ, '')
				list( 	ANCSEQ.CP1= sapply(tmp, function(x) paste(x[seq.int(1,nchar(ANCSEQ[1]),3)],collapse='')  ),
						ANCSEQ.CP2= sapply(tmp, function(x) paste(x[seq.int(2,nchar(ANCSEQ[1]),3)],collapse='') ),
						ANCSEQ.CP3= sapply(tmp, function(x) paste(x[seq.int(3,nchar(ANCSEQ[1]),3)],collapse='')  ), ROOT_CALENDAR_TIME=ROOT_CALENDAR_TIME, IDCLU=IDCLU, NEWICK=NEWICK	)				
			}, by='GENE']
	df.ph	<- melt(df.ph, id.var=c('ROOT_CALENDAR_TIME','IDCLU','NEWICK','GENE'), value.name='ANCSEQ', variable.name='CODON_POS', variable.factor=FALSE)
	set(df.ph, NULL, 'CODON_POS', df.ph[, substr(CODON_POS, 8, nchar(CODON_POS))])
	#
	#	save to file all we need to call SeqGen
	#
	file	<- paste(outdir.sg,'/',infile.prefix, 'seqgen.R',sep='')
	cat(paste('\nsave to file=',file))
	df.seqgen	<- df.ph
	save(df.seqgen, log.df, df.nodestat, file=file)
}
##--------------------------------------------------------------------------------------------------------
##	program to generate files for Seq Gen from output of Matt s VirusTreeSimulator
##	olli originally written 26-08-2014
##--------------------------------------------------------------------------------------------------------
prog.PANGEA.SeqGen.createInputFile.v1<- function()
{
	stop()
	verbose			<- 1
	with.plot		<- 1
	label.sep		<- '|'	
	#
	#	read I/O
	#
	indir.epi		<- '/Users/Oliver/git/HPTN071sim/tmp140914/TrChains'
	infile.epi		<- '140716_RUN001_SAVE.R'	
	indir.vts		<- '/Users/Oliver/git/HPTN071sim/tmp140914/VirusTreeSimulator'
	infile.prefix	<- '140716_RUN001_'	
	infile.args		<- '/Users/Oliver/git/HPTN071sim/tmp140914/140716_RUN001_PipeArgs.R'
	outdir.sg		<- '/Users/Oliver/git/HPTN071sim/tmp140914/SeqGen'	
	if(exists("argv"))
	{
		#	args input
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									indir.epi= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir.epi<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.epi= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.epi<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									indir.vts= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir.vts<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.vts= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.prefix<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									infile.args= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.args<- tmp[1]		
		#	args output
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir.sg<- tmp[1]		
		
	}
	if(verbose)
	{
		cat('\ninput args\n',paste(indir.epi, infile.epi, indir.vts, infile.prefix, outdir.sg, sep='\n'))
	}	
	if(!is.na(infile.args))
	{
		load(infile.args)	#expect 'pipeline.args'
	}
	if(is.null(pipeline.args))
	{
		cat('\nCould not find pipeline.args, generating default')
		pipeline.args	<- rPANGEAHIVsim.pipeline.args()
	}	
	stopifnot( all( c('s.seed','wher.mu','wher.sigma','bwerm.mu','bwerm.sigma')%in%pipeline.args[, stat] ) )
	#
	#	setup samplers
	#
	cat(paste('\ncreate sampler of evolutionary rates'))
	#	create sampler of within host evolutionary rates
	rER.pol			<- PANGEA.WithinHostEvolutionaryRate.create.sampler.v1(wher.mu=pipeline.args['wher.mu',][, as.numeric(v)], wher.sigma=pipeline.args['wher.sigma',][, as.numeric(v)])
	#	create sampler of between host evolutionary rate dampener
	rER.bwm			<- PANGEA.BetweenHostEvolutionaryRateModifier.create.sampler.v1(bwerm.mu=pipeline.args['bwerm.mu',][, as.numeric(v)], bwerm.sigma=pipeline.args['bwerm.sigma',][, as.numeric(v)])
	#	create sampler of ancestral sequences
	cat(paste('\ncreate sampler of ancestral sequences'))
	tmp				<- PANGEA.RootSeq.create.sampler(root.ctime.grace= 0.5, sample.grace= 3)
	rANCSEQ			<- tmp$rANCSEQ
	rANCSEQ.args	<- tmp$rANCSEQ.args 	
	#	read GTR parameters
	log.df			<- PANGEA.GTR.params()	
	#
	#
	#
	file		<- paste(indir.epi, '/', infile.epi, sep='')
	load(file)	#expect "df.epi"    "df.trms"   "df.inds"   "df.sample"
	#	
	infiles		<- list.files(indir.vts)
	tmp			<- paste('^',infile.prefix,'.*nex$',sep='')
	infiles		<- infiles[ grepl(tmp, infiles)  ]	
	if(!length(infiles))	stop('cannot find files matching criteria')
	#
	set.seed( pipeline.args['s.seed',][, as.numeric(v)] )
	#
	#	read from VirusTreeSimulator and convert branch lengths in time to branch lengths in subst/site
	#
	df.ph			<- vector('list', length(infiles))
	df.nodestat		<- vector('list', length(infiles))	
	for(i in seq_along(infiles))
	{				
		#i<- 11
		infile			<- infiles[i]
		cat(paste('\nprocess file',i,infile))
		file			<- paste(indir.vts, '/', infile, sep='')
		#	read brl, units from annotated nexus file. attention: () may not contain two nodes
		tmp				<- hivc.beast2out.read.nexus.and.stats(file, method.node.stat='any.node')
		ph				<- tmp$tree
		node.stat		<- tmp$node.stat
		node.stat		<- subset(node.stat, STAT=='Unit')
		set(node.stat, NULL, 'VALUE', node.stat[, gsub('\"','',VALUE)])
		node.stat[, IDPOP:= as.integer(node.stat[,substr(VALUE, 4, nchar(VALUE))])]
		node.stat		<- merge(subset(df.inds, select=c(IDPOP, GENDER, DOB, TIME_SEQ, IDCLU)), subset(node.stat, select=c(IDPOP, NODE_ID)), by='IDPOP')	
		#
		#	create collapsed Newick tree with expected substitutions / site for each branch 
		#
		#	draw evolutionary rates for every individual in the transmission chain
		node.stat		<- merge(node.stat, data.table( IDPOP=node.stat[, unique(IDPOP)], ER= rER.pol(node.stat[, length(unique(IDPOP))]) ), by='IDPOP')
		#	smaller ER for transmission edge
		node.stat[, BWM:= rER.bwm(nrow(node.stat))]
		#	set BWM to 1 for all non-transmission edges
		tmp				<- ph$edge
		colnames(tmp)	<- c('NODE_ID','NODE_ID_TO')
		tmp				<- merge( subset(node.stat, select=c(NODE_ID, IDPOP)), data.table(tmp, EDGE_ID=seq_len(nrow(tmp))), by='NODE_ID') 
		setnames(tmp, c('NODE_ID','IDPOP','NODE_ID_TO'), c('NODE_ID_FROM','IDPOP_FROM','NODE_ID'))
		tmp				<- merge( subset(node.stat, select=c(NODE_ID, IDPOP)), tmp, by='NODE_ID') 
		tmp				<- subset( tmp, IDPOP!=IDPOP_FROM)[, NODE_ID_FROM]	#ER slows down in transmitter leading to infection. want to slow down edge leading to NODE_FROM.
		set( node.stat, node.stat[, which(!NODE_ID%in%tmp)], 'BWM', 1. )
		#	no ER possible for root node - there s no edge leading to it
		set(node.stat, node.stat[, which(NODE_ID==Ntip(ph)+1)], c('ER','BWM'), NA_real_)		
		#	set root edge evolutionary rate to overall mean between-host rate
		#	get NODE_ID of edge from root
		tmp				<- ph$edge[match(Ntip(ph)+1, ph$edge[1, ]), 2]
		tmp				<- node.stat[, which(NODE_ID==tmp)]		
		set(node.stat, tmp, 'ER', log.df[1,meanRate] )		
		stopifnot( node.stat[tmp, BWM]>1)
		set(node.stat, tmp, 'BWM', 1. )		# no need to further slow down root edge
		#	check calendar time of root in simulated phylogeny for consistency
		tmp				<- seq.collapse.singles(ph)
		tmp2			<- regmatches(tmp$tip.label[1], regexpr('ID_[0-9]+',tmp$tip.label[1]))
		tmp2			<- as.numeric(substr(tmp2, 4, nchar(tmp2)))
		tmp2			<- subset(node.stat, IDPOP==tmp2)[1, TIME_SEQ]
		root.ctime		<- ifelse(Nnode(tmp), tmp2 - (node.depth.edgelength(tmp)[1] + tmp$root.edge), tmp2-tmp$root.edge)		
		tmp				<- subset(node.stat, IDPOP<0)[, unique(IDPOP)]
		stopifnot(length(tmp)==1)
		stopifnot(subset(df.trms, IDTR==tmp)[, round(IDTR_TIME_INFECTED, d=1)]==round(root.ctime, d=1))
		#	set expected numbers of substitutions per branch within individual IDPOP
		setkey(node.stat, NODE_ID)
		ph$edge.length	<- ph$edge.length * node.stat[ ph$edge[, 2], ][, ER / BWM]
		stopifnot(all(!is.na(ph$edge.length)))
		#	once expected number of substitutions / site are simulated, can collapse singleton nodes
		ph				<- seq.collapse.singles(ph)	
		#	set tip label so that IDPOP can be checked for consistency	
		if(pipeline.args['epi.model',][,v]=='HPTN071')
			node.stat[, LABEL:= node.stat[, paste('IDPOP_',IDPOP,label.sep,GENDER,label.sep,'DOB_',round(DOB,d=3),label.sep,round(TIME_SEQ,d=3),sep='')]]
		if(pipeline.args['epi.model',][,v]=='DSPS')
			node.stat[, LABEL:= node.stat[, paste('IDPOP_',IDPOP,label.sep,GENDER,label.sep,'DOB_',NA,label.sep,round(TIME_SEQ,d=3),sep='')]]		
		setkey(node.stat, NODE_ID)
		ph$tip.label	<- node.stat[seq_len(Ntip(ph)), ][, LABEL]
		#
		df.nodestat[[i]]<- node.stat
		tmp				<- ifelse(Nnode(ph), write.tree(ph, digits = 10), paste( '(',ph$tip.label,':',ph$root.edge,',NOEXIST_NA|NA|DOB_NA|',root.ctime,':0):0;', sep='')	)
		df.ph[[i]]		<- data.table(ROOT_CALENDAR_TIME= root.ctime, IDCLU=node.stat[, unique(IDCLU)], NEWICK=tmp)
		#readline()
	}
	df.ph		<- do.call('rbind',df.ph)
	df.nodestat	<- do.call('rbind',df.nodestat)
	#	check that we have exactly one root edge with overall mean between host rate per cluster
	stopifnot( df.nodestat[, length(unique(IDCLU))]==nrow(subset(df.nodestat, ER==log.df$meanRate[1])) )
	
	if(with.plot)
	{
		#ggplot(subset(df.nodestat, ER!=log.df$meanRate[1]) , aes(x=ER/BWM)) + geom_histogram(binwidth=0.001)
		#ggplot(subset(df.nodestat, ER!=log.df$meanRate[1]), aes(x=ER, y=BWM)) + geom_point()	
		#	plot used within-host ERs
		ggplot(subset(df.nodestat, ER!=log.df$meanRate[1]), aes(x=ER)) + geom_histogram(binwidth=0.001)	+ labs(x='simulated within-host evolutionary rate') +
				scale_x_continuous(breaks= seq(0, 0.02, 0.002))
		file	<- paste(indir.epi, '/', substr(infile.epi,1,nchar(infile.epi)-6),'INFO_sg_ER.pdf', sep='')
		ggsave(file, w=6, h=6)
		#	plot used between host modifiers
		ggplot(subset(df.nodestat, BWM!=1) , aes(x=BWM)) + geom_histogram(binwidth=0.05) + labs(x='simulated between-host rate multiplier') +
				scale_x_continuous(breaks= seq(0.8, 2.4, 0.2))
		file	<- paste(indir.epi, '/', substr(infile.epi,1,nchar(infile.epi)-6),'INFO_sg_BWM.pdf', sep='')
		ggsave(file, w=6, h=6)
		#	plot used ERs along transmission edges
		ggplot(subset(df.nodestat, BWM!=1) , aes(x=ER/BWM)) + geom_histogram(binwidth=0.001) + labs(x='simulated evolutionary rates along transmission edges') +
				scale_x_continuous(breaks= seq(0, 0.02, 0.002))
		file	<- paste(indir.epi, '/', substr(infile.epi,1,nchar(infile.epi)-6),'INFO_sg_BWER.pdf', sep='')
		ggsave(file, w=6, h=6)		
	}	
	#
	#	draw ancestral sequences and add to df.ph
	#
	root.ctime		<- df.ph[, ROOT_CALENDAR_TIME]
	ancseq			<- rANCSEQ(root.ctime, rANCSEQ.args)
	ancseq			<- data.table(ANCSEQ= apply(as.character(ancseq),1,function(x) paste(x, collapse='')) )		
	df.ph			<- cbind(df.ph, ancseq)
	#
	#	create SEQ-GEN input file
	#	
	partition.len	<- c( ncol(rANCSEQ.args$anc.seq.gag), ncol(rANCSEQ.args$anc.seq.pol), ncol(rANCSEQ.args$anc.seq.env) )
	#partition.er	<- c( 2.5, 4, 5 )
	#	split ancestral sequence into GENE and CODON_POS
	df.ph[, ANCSEQ.GAG:= substr(ANCSEQ, 1, partition.len[1])]
	df.ph[, ANCSEQ.POL:= substr(ANCSEQ, partition.len[1]+1, partition.len[1]+partition.len[2])]
	df.ph[, ANCSEQ.ENV:= substr(ANCSEQ, partition.len[1]+partition.len[2]+1, partition.len[1]+partition.len[2]+partition.len[3])]	
	stopifnot( all( strsplit( df.ph[1, ANCSEQ], '')[[1]] == strsplit( paste( df.ph[1, ANCSEQ.GAG], df.ph[1, ANCSEQ.POL], df.ph[1, ANCSEQ.ENV], sep=''), '' )[[1]] ) )
	df.ph[, ANCSEQ:=NULL]
	df.ph	<- melt(df.ph, id.var=c('ROOT_CALENDAR_TIME','IDCLU','NEWICK'), value.name='ANCSEQ', variable.name='GENE', variable.factor=FALSE)
	set(df.ph, NULL, 'GENE', df.ph[, substr(GENE, 8, nchar(GENE))])	
	df.ph	<- df.ph[,  {
				tmp	<- strsplit(ANCSEQ, '')
				list( 	ANCSEQ.CP1= sapply(tmp, function(x) paste(x[seq.int(1,nchar(ANCSEQ[1]),3)],collapse='')  ),
						ANCSEQ.CP2= sapply(tmp, function(x) paste(x[seq.int(2,nchar(ANCSEQ[1]),3)],collapse='') ),
						ANCSEQ.CP3= sapply(tmp, function(x) paste(x[seq.int(3,nchar(ANCSEQ[1]),3)],collapse='')  ), ROOT_CALENDAR_TIME=ROOT_CALENDAR_TIME, IDCLU=IDCLU, NEWICK=NEWICK	)				
			}, by='GENE']
	df.ph	<- melt(df.ph, id.var=c('ROOT_CALENDAR_TIME','IDCLU','NEWICK','GENE'), value.name='ANCSEQ', variable.name='CODON_POS', variable.factor=FALSE)
	set(df.ph, NULL, 'CODON_POS', df.ph[, substr(CODON_POS, 8, nchar(CODON_POS))])
	#
	#	save to file all we need to call SeqGen
	#
	file	<- paste(outdir.sg,'/',infile.prefix, 'seqgen.R',sep='')
	cat(paste('\nsave to file=',file))
	df.seqgen	<- df.ph
	save(df.seqgen, log.df, df.nodestat, file=file)
}
######################################################################################
#	Program to simulate sequences with Seq-Gen-1.3.2 via phyclust 	
#	olli originally written 15-09-2014
######################################################################################
prog.PANGEA.SeqGen.run.WINDOWScompatible<- function()
{	
	indir.epi			<- '/Users/Oliver/git/HPTN071sim/tmp140914/TrChains'
	infile.epi			<- '140716_RUN001_SAVE.R'		
	indir.sg			<- '/Users/Oliver/git/HPTN071sim/tmp140908/SeqGen'
	infile.prefix		<- '140716_RUN001_'
	infile.args			<- NA
	outdir				<- '/Users/Oliver/git/HPTN071sim/tmp140908'
	with.plot			<- 1	
	verbose				<- 1
	label.idx.codonpos	<- 1
	label.idx.gene		<- 2
	label.idx.clu		<- 3
	treelabel.idx.idpop	<- 1
	treelabel.idx.sep	<- '|'	
	#
	if(exists("argv"))
	{
		#	args input
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									indir.epi= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir.epi<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.epi= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.epi<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									indir.sg= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir.sg<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									infile.sg= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.prefix<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									infile.args= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.args<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]				
	}
	if(verbose)
	{
		cat('\ninput args\n',paste(indir.sg, infile.prefix, sep='\n'))
	}
	if(!is.na(infile.args))
	{
		load(infile.args)	#expect 'pipeline.args'
	}
	if(is.null(pipeline.args))
	{
		cat('\nCould not find pipeline.args, generating default')
		pipeline.args	<- rPANGEAHIVsim.pipeline.args()
	}	
	stopifnot( all( c('s.seed')%in%pipeline.args[, stat] ) )	
	set.seed(pipeline.args['s.seed',][, as.numeric(v)])
	#
	file		<- paste(indir.epi, '/', infile.epi, sep='')
	load(file)	#expect "df.epi"    "df.trms"   "df.inds"   "df.sample"
	#
	file		<- paste(indir.sg,'/',infile.prefix, 'seqgen.R',sep='')
	cat(paste('\nLoad seqgen R input file, file=',file))
	load(file)	#expect df.seqgen, gtr.central, log.df, df.nodestat
	#
	log.df[, IDX:= seq_len(nrow(log.df))]
	log.df[, FILE:=NULL]	 
	#	sample GTR parameters from posterior
	tmp			<- df.seqgen[, {
				tmp	<- log.df[['IDX']][ which( log.df[['GENE']]==GENE & log.df[['CODON_POS']]==CODON_POS ) ]
				stopifnot(length(tmp)>1)
				list( IDCLU=IDCLU, IDX=sample(tmp, length(IDCLU), replace=FALSE) )	
			}, by=c('GENE','CODON_POS')]
	tmp			<- merge(tmp, log.df, by=c('GENE', 'CODON_POS', 'IDX'))
	df.seqgen	<- merge(df.seqgen, tmp, by=c('GENE','CODON_POS','IDCLU'))
	#
	if(with.plot)
	{
		tmp			<- subset(df.nodestat, ER!=log.df$meanRate[1], select=c(IDPOP, ER, BWM, IDCLU))
		tmp			<- merge(tmp, subset(df.seqgen, select=c(GENE, CODON_POS, IDCLU, mu)), by='IDCLU', allow.cartesian=TRUE)
		set(tmp, NULL, 'ER', tmp[, ER*mu])		
		
		ggplot(tmp, aes(x=CODON_POS, y=ER, colour=CODON_POS, group=CODON_POS)) + geom_boxplot() +				
				facet_grid(.~GENE, scales='free_y') +
				scale_colour_discrete(guide=FALSE) +
				scale_y_continuous(breaks= seq(0, 0.05, 0.002)) + labs(linetype='Gene', y='simulated within-host evolutionary rate', x='codon position')
		file		<- paste(indir.sg,'/',infile.prefix, 'ER_by_gene.pdf',sep='')
		ggsave(file=file, w=6, h=6)						
	}			
	#
	#	run SeqGen-1.3.2 (no seed) 
	#
	df.ph.out	<- df.seqgen[,	{									
				file	<- paste(indir.sg,'/',infile.prefix, IDCLU, '_', GENE, '_', CODON_POS,'.phy',sep='')
				cat(paste('\nwrite to file',file))
				opts	<- paste('-n1 -k1 -or -mGTR -a',alpha,' -g4 -i0 -s',mu,' -f',a,',',c,',',g,',',t,' -r',ac,',',ag,',',at,',',cg,',',1,',',gt,sep='')
				input 	<- c(paste(" 1 ", nchar(ANCSEQ), sep=''), paste('ANCSEQ\t', toupper(ANCSEQ), sep = ''), 1, NEWICK)
				z		<- seqgen(opts, input=input, temp.file=file)				
			}, by=c('GENE','CODON_POS','IDCLU')]
	#
	#	process SeqGen runs
	#
	#	collect simulated sequences
	infiles		<- list.files(indir.sg)
	infiles		<- infiles[ grepl('*phy$', infiles)  ]	
	if(!length(infiles))	stop('cannot find files matching criteria')		
	infile.df	<- data.table(FILE=infiles)
	tmp			<- infile.df[, strsplit(FILE, '_') ]
	infile.df[, CODON_POS:= sapply(tmp, function(x) rev(x)[label.idx.codonpos])]
	infile.df[, GENE:= sapply(tmp, function(x) rev(x)[label.idx.gene])]
	infile.df[, IDCLU:= sapply(tmp, function(x) rev(x)[label.idx.clu])]
	set(infile.df, NULL, 'CODON_POS', infile.df[, substr(CODON_POS,1,3)])
	cat(paste('\nFound sequences for clusters, nclu=', infile.df[, length(unique(IDCLU))]))
	#
	#	read simulated sequences
	#
	df.seq		<- infile.df[,	{
				cat(paste('\nread seq in file',FILE))									
				file	<- paste(indir.sg,'/',FILE,sep='')
				tmp		<- as.character(read.dna(file, format='sequential'))
				tmp		<- tmp[!grepl('NOEXIST',rownames(tmp)), , drop=FALSE]
				list( SEQ=apply(tmp,1,function(x) paste(x, collapse='')), LABEL=rownames(tmp) )				
			}, by='FILE']
	df.seq		<- merge(df.seq, infile.df, by='FILE')
	#
	#	reconstruct genes from codon positions
	#
	df.seq[, STAT:=paste(GENE,CODON_POS,sep='.')]		
	df.seq		<- dcast.data.table(df.seq, IDCLU + LABEL ~ STAT, value.var="SEQ")
	#	check that seq of correct size
	stopifnot( df.seq[, length(unique(LABEL))]==nrow(df.seq) )
	stopifnot( df.seq[, all( nchar(GAG.CP1)==nchar(GAG.CP2) & nchar(GAG.CP1)==nchar(GAG.CP3) )] )
	stopifnot( df.seq[, all( nchar(POL.CP1)==nchar(POL.CP2) & nchar(POL.CP1)==nchar(POL.CP3) )] )
	stopifnot( df.seq[, all( nchar(ENV.CP1)==nchar(ENV.CP2) & nchar(ENV.CP1)==nchar(ENV.CP3) )] )
	#
	df.seq		<- df.seq[, {
				tmp		<- do.call('rbind',sapply(list(ENV.CP1,ENV.CP2,ENV.CP3), strsplit, ''))
				env		<- paste(as.vector(tmp), collapse='')
				tmp		<- do.call('rbind',sapply(list(GAG.CP1,GAG.CP2,GAG.CP3), strsplit, ''))
				gag		<- paste(as.vector(tmp), collapse='')
				tmp		<- do.call('rbind',sapply(list(POL.CP1,POL.CP2,POL.CP3), strsplit, ''))
				pol		<- paste(as.vector(tmp), collapse='')
				list(GAG=gag, POL=pol, ENV=env, IDCLU=IDCLU)
			}, by=c('LABEL')]
	#	check that we have indeed seq for all sampled individuals
	df.seq		<- subset( df.seq, !grepl('NOEXIST',LABEL) )	
	tmp			<- df.seq[, sapply( strsplit(LABEL, treelabel.idx.sep, fixed=TRUE), '[[', treelabel.idx.idpop )]
	df.seq[, IDPOP:=as.integer(substr(tmp,7,nchar(tmp)))]	
	stopifnot( setequal( subset( df.inds, !is.na(TIME_SEQ) )[, IDPOP], df.seq[,IDPOP]) )
	#	merge simulated data
	simu.df		<- merge(subset(df.seq, select=c(IDPOP, GAG, POL, ENV)), subset( df.inds, !is.na(TIME_SEQ) ), by='IDPOP', all.x=TRUE)
	#
	#	save simulated data -- internal
	#
	file			<- paste(outdir, '/', infile.prefix, 'SIMULATED_INTERNAL.R', sep='')
	cat(paste('\nwrite to file', file))
	save(df.epi, df.trms, df.inds, df.sample, df.seq, file=file)
	#
	#	save simulated data -- to be shared
	#	
	if(pipeline.args['epi.model'][,v]=='HPTN071')
		tmp			<- subset( df.inds, !is.na(TIME_SEQ), select=c(IDPOP, GENDER, CIRCM, DOB, DOD, TIME_SEQ ) )
	if(pipeline.args['epi.model'][,v]=='DSPS')
		tmp			<- subset( df.inds, !is.na(TIME_SEQ), select=c(IDPOP, GENDER, TIME_SEQ ) )
	file			<- paste(outdir, '/', infile.prefix, 'SIMULATED_metadata.csv', sep='')
	cat(paste('\nwrite to file', file))
	write.csv(tmp, file)		
	tmp				<- tolower(do.call('rbind',strsplit(df.seq[, GAG],'')))
	rownames(tmp)	<- df.seq[, LABEL]
	tmp				<- as.DNAbin(tmp)
	file			<- paste(outdir, '/', infile.prefix, 'SIMULATED_gag.fa', sep='')
	write.dna(tmp, file, format = "fasta")	
	tmp				<- tolower(do.call('rbind',strsplit(df.seq[, POL],'')))
	rownames(tmp)	<- df.seq[, LABEL]
	tmp				<- as.DNAbin(tmp)
	file			<- paste(outdir, '/', infile.prefix, 'SIMULATED_pol.fa', sep='')
	write.dna(tmp, file, format = "fasta")	
	tmp				<- tolower(do.call('rbind',strsplit(df.seq[, ENV],'')))
	rownames(tmp)	<- df.seq[, LABEL]
	tmp				<- as.DNAbin(tmp)
	file			<- paste(outdir, '/', infile.prefix, 'SIMULATED_env.fa', sep='')
	write.dna(tmp, file, format = "fasta")
	#
	#	zip simulated files
	#
	tmp				<- c( paste(outdir, '/', infile.prefix, 'SIMULATED_metadata.csv', sep=''), paste(outdir, '/', infile.prefix, 'SIMULATED_env.fa', sep=''), paste(outdir, '/', infile.prefix, 'SIMULATED_pol.fa', sep=''), paste(outdir, '/', infile.prefix, 'SIMULATED_gag.fa', sep='') )
	zip( paste(outdir, '/', infile.prefix, 'SIMULATED.zip', sep=''), tmp, flags = "-FSr9XTj")
	dummy			<- file.remove(tmp)
	#
	#	zip internal files
	#
	tmp				<- list.files(outdir, pattern='*pdf$', recursive=TRUE, full.names=TRUE) 
	file.copy(tmp, outdir, overwrite=TRUE)	 
	tmp				<- c( paste(outdir, '/', infile.prefix, 'SIMULATED_INTERNAL.R', sep=''), list.files(outdir, pattern='*pdf$', recursive=FALSE, full.names=TRUE) ) 	
	zip( paste(outdir, '/', infile.prefix, 'INTERNAL.zip', sep=''), tmp, flags = "-FSr9XTj")
	dummy			<- file.remove(tmp)
	
	return(1)
}
######################################################################################
#	Program to simulate sequences with Seq-Gen-1.3.3 	
#	olli originally written 09-09-2014
#	used up to version prog.HPTN071.input.parser.v3
######################################################################################
#' @title Program to simulate gene sequences
#' @description \code{prog.PANGEA.SeqGen.run} reads file \code{infile.sg} in directory \code{indir.sg} that was
#' created with the \code{SeqGen} input file creator. The simulated partial sequences are collected, coerced back
#' into Gag, Pol, Env genes, and written in fasta format to directory \code{outdir}. Patient Metavariables are 
#' stored in the same directory, and zip files are created.
#' @return NULL. Saves zip files with simulations.
#' @example example/ex.seqgen.run.R
prog.PANGEA.SeqGen.run<- function()
{	
	indir.epi			<- '/Users/Oliver/git/HPTN071sim/tmp140914/TrChains'
	infile.epi			<- '140716_RUN001_SAVE.R'		
	indir.sg			<- '/Users/Oliver/git/HPTN071sim/tmp140908/SeqGen'
	infile.prefix		<- '140716_RUN001_'
	infile.args			<- NA
	outdir				<- '/Users/Oliver/git/HPTN071sim/tmp140908'
	with.plot			<- 1	
	verbose				<- 1
	label.idx.codonpos	<- 1
	label.idx.gene		<- 2
	label.idx.clu		<- 3
	treelabel.idx.idpop	<- 1
	treelabel.idx.sep	<- '|'	
	#
	if(exists("argv"))
	{
		#	args input
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									indir.epi= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir.epi<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.epi= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.epi<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									indir.sg= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir.sg<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									infile.sg= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.prefix<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									infile.args= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.args<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]				
	}
	if(verbose)
	{
		cat('\ninput args\n',paste(indir.sg, infile.prefix, sep='\n'))
	}
	if(!is.na(infile.args))
	{
		load(infile.args)	#expect 'pipeline.args'
	}
	if(is.null(pipeline.args))
	{
		cat('\nCould not find pipeline.args, generating default')
		pipeline.args	<- rPANGEAHIVsim.pipeline.args()
	}	
	#
	file		<- paste(indir.epi, '/', infile.epi, sep='')
	load(file)	#expect "df.epi"    "df.trms"   "df.inds"   "df.sample"
	#
	file			<- paste(indir.sg,'/',infile.prefix, 'seqgen.R',sep='')
	cat(paste('\nLoad seqgen R input file, file=',file))
	load(file)	#expect df.seqgen, gtr.central, log.df, df.nodestat
	#
	if( pipeline.args['index.starttime.mode',][,v]=='shift' )		
		root.edge.rate	<- 1e-6
	if( pipeline.args['index.starttime.mode',][,v]!='shift' )
		root.edge.rate	<- log.df[1,meanRate]		
	stopifnot( all( c('s.seed')%in%pipeline.args[, stat] ) )	
	set.seed(pipeline.args['s.seed',][, as.numeric(v)])	
	#
	#	create SeqGen input files
	#
	df.ph.out	<- df.seqgen[,	{
				#print( table( strsplit(ANCSEQ, ''), useNA='if') )
				file	<- paste(indir.sg,'/',infile.prefix, IDCLU, '_', GENE, '_', CODON_POS,'.seqgen',sep='')
				cat(paste('\nwrite to file',file))
				txt		<- paste('1\t',nchar(ANCSEQ),'\n',sep='')
				txt		<- paste(txt, 'ANCSEQ\t',toupper(ANCSEQ),'\n',sep='')					
				txt		<- paste(txt, '1\n',sep='')
				txt		<- paste(txt, NEWICK, '\n',sep='')
				cat(txt, file=file)
				list(FILE= paste(infile.prefix, IDCLU, '_', GENE, '_', CODON_POS,'.seqgen',sep='') )
				# ./seq-gen -mHKY -t3.0 -f0.3,0.2,0.2,0.3 -n1 -k1 -on < /Users/Oliver/git/HPTN071sim/tmp/SeqGen/140716_RUN001_50.seqgen > example.nex
			}, by=c('GENE','CODON_POS','IDCLU')]
	df.ph.out	<- df.ph.out[, 	{
				tmp	<- log.df[['IDX']][ which( log.df[['GENE']]==GENE & log.df[['CODON_POS']]==CODON_POS ) ]
				stopifnot(length(tmp)>1)
				list( FILE=FILE, IDCLU=IDCLU, IDX=sample(tmp, length(FILE), replace=FALSE) )
			}, by=c('GENE','CODON_POS')]
	#	sample GTR parameters from posterior
	df.ph.out	<- merge(df.ph.out, log.df, by=c('GENE', 'CODON_POS', 'IDX'))
	#
	if(with.plot)
	{
		tmp			<- subset(df.nodestat, ER!=root.edge.rate, select=c(IDPOP, ER, BWM, IDCLU))
		tmp			<- merge(tmp, subset(df.ph.out, select=c(GENE, CODON_POS, IDCLU, mu)), by='IDCLU', allow.cartesian=TRUE)
		set(tmp, NULL, 'ER', tmp[, ER*mu])		
		
		ggplot(tmp, aes(x=CODON_POS, y=ER, colour=CODON_POS, group=CODON_POS)) + geom_boxplot() +				
				facet_grid(.~GENE, scales='free_y') +
				scale_colour_discrete(guide=FALSE) +
				scale_y_continuous(breaks= seq(0, 0.05, 0.002)) + labs(linetype='Gene', y='simulated within-host evolutionary rate', x='codon position')
		file		<- paste(indir.sg,'/',infile.prefix, 'ER_by_gene.pdf',sep='')
		ggsave(file=file, w=6, h=6)						
	}	
	#
	#	call SeqGen command line
	#
	cat(paste('\nUsing Gamma rate variation, gamma=',pipeline.args['er.gamma',][, as.numeric(v)]))
	tmp		<- df.ph.out[, {	
				cat(paste('\nProcess', IDCLU, GENE, CODON_POS))				
				cmd	<- cmd.SeqGen(indir.sg, FILE, indir.sg, gsub('seqgen','phy',FILE), prog=PR.SEQGEN, prog.args=paste('-n',1,' -k1 -or -z',pipeline.args['s.seed',][, as.numeric(v)],sep=''), 
						alpha=alpha, gamma=pipeline.args['er.gamma',][, as.numeric(v)], invariable=0, scale=mu, freq.A=a, freq.C=c, freq.G=g, freq.T=t,
						rate.AC=ac, rate.AG=ag, rate.AT=at, rate.CG=cg, rate.CT=1, rate.GT=gt)
				system(cmd)				
				list(CMD=cmd)							
			}, by='FILE']
	#
	#	process SeqGen runs
	#
	#	collect simulated sequences
	infiles		<- list.files(indir.sg)
	infiles		<- infiles[ grepl('*phy$', infiles)  ]	
	if(!length(infiles))	stop('cannot find files matching criteria')		
	infile.df	<- data.table(FILE=infiles)
	tmp			<- infile.df[, strsplit(FILE, '_') ]
	infile.df[, CODON_POS:= sapply(tmp, function(x) rev(x)[label.idx.codonpos])]
	infile.df[, GENE:= sapply(tmp, function(x) rev(x)[label.idx.gene])]
	infile.df[, IDCLU:= sapply(tmp, function(x) rev(x)[label.idx.clu])]
	set(infile.df, NULL, 'CODON_POS', infile.df[, substr(CODON_POS,1,3)])
	cat(paste('\nFound sequences for clusters, nclu=', infile.df[, length(unique(IDCLU))]))
	#
	#	read simulated sequences
	#
	df.seq		<- infile.df[,	{
				cat(paste('\nread seq in file',FILE))
				file	<- paste(indir.sg,'/',FILE,sep='')
				tmp		<- as.character(read.dna(file, format='sequential'))
				list( SEQ=apply(tmp,1,function(x) paste(x, collapse='')), LABEL=rownames(tmp) )				
			}, by='FILE']
	df.seq		<- merge(df.seq, infile.df, by='FILE')
	#
	#	reconstruct genes from codon positions
	#
	df.seq[, STAT:=paste(GENE,CODON_POS,sep='.')]		
	df.seq		<- dcast.data.table(df.seq, IDCLU + LABEL ~ STAT, value.var="SEQ")
	#	check that seq of correct size
	stopifnot( df.seq[, all( nchar(GAG.CP1)==nchar(GAG.CP2) & nchar(GAG.CP1)==nchar(GAG.CP3) )] )
	stopifnot( df.seq[, all( nchar(POL.CP1)==nchar(POL.CP2) & nchar(POL.CP1)==nchar(POL.CP3) )] )
	stopifnot( df.seq[, all( nchar(ENV.CP1)==nchar(ENV.CP2) & nchar(ENV.CP1)==nchar(ENV.CP3) )] )
	#
	df.seq		<- df.seq[, {
				tmp		<- do.call('rbind',sapply(list(ENV.CP1,ENV.CP2,ENV.CP3), strsplit, ''))
				env		<- paste(as.vector(tmp), collapse='')
				tmp		<- do.call('rbind',sapply(list(GAG.CP1,GAG.CP2,GAG.CP3), strsplit, ''))
				gag		<- paste(as.vector(tmp), collapse='')
				tmp		<- do.call('rbind',sapply(list(POL.CP1,POL.CP2,POL.CP3), strsplit, ''))
				pol		<- paste(as.vector(tmp), collapse='')
				list(GAG=gag, POL=pol, ENV=env, IDCLU=IDCLU)
			}, by=c('LABEL')]
	#	check that we have indeed seq for all sampled individuals
	df.seq		<- subset( df.seq, !grepl('NOEXIST',LABEL) )	
	tmp			<- df.seq[, sapply( strsplit(LABEL, treelabel.idx.sep, fixed=TRUE), '[[', treelabel.idx.idpop )]
	df.seq[, IDPOP:=as.integer(substr(tmp,7,nchar(tmp)))]	
	stopifnot( setequal( subset( df.inds, !is.na(TIME_SEQ) )[, IDPOP], df.seq[,IDPOP]) )
	#
	#	save simulated data -- internal
	#
	file			<- paste(outdir, '/', infile.prefix, 'SIMULATED_INTERNAL.R', sep='')
	cat(paste('\nwrite to file', file))
	save(df.epi, df.trms, df.inds, df.sample, df.seq, file=file)
	#
	#	save simulated data -- to be shared
	#	
	if(pipeline.args['epi.model'][,v]=='HPTN071')
	{
		tmp			<- subset( df.inds, !is.na(TIME_SEQ), select=c(IDPOP, GENDER, CIRCM, DOB, DOD, TIME_SEQ, CD4_SEQ, INCIDENT_SEQ ) )
		setnames(tmp, 'INCIDENT_SEQ', 'INCIDENT_WITHIN1YEAR_SEQ')
	}		
	if(pipeline.args['epi.model'][,v]=='DSPS')
		tmp			<- subset( df.inds, !is.na(TIME_SEQ), select=c(IDPOP, GENDER, TIME_SEQ ) )
	file			<- paste(outdir, '/', infile.prefix, 'SIMULATED_metadata.csv', sep='')
	cat(paste('\nwrite to file', file))
	write.csv(tmp, file)		
	tmp				<- tolower(do.call('rbind',strsplit(df.seq[, GAG],'')))
	rownames(tmp)	<- df.seq[, LABEL]
	df.seq.gag		<- as.DNAbin(tmp)
	file			<- paste(outdir, '/', infile.prefix, 'SIMULATED_gag.fa', sep='')
	write.dna(df.seq.gag, file, format = "fasta")	
	tmp				<- tolower(do.call('rbind',strsplit(df.seq[, POL],'')))
	rownames(tmp)	<- df.seq[, LABEL]
	df.seq.pol		<- as.DNAbin(tmp)
	file			<- paste(outdir, '/', infile.prefix, 'SIMULATED_pol.fa', sep='')
	write.dna(df.seq.pol, file, format = "fasta")	
	tmp				<- tolower(do.call('rbind',strsplit(df.seq[, ENV],'')))
	rownames(tmp)	<- df.seq[, LABEL]
	df.seq.env		<- as.DNAbin(tmp)
	file			<- paste(outdir, '/', infile.prefix, 'SIMULATED_env.fa', sep='')
	write.dna(df.seq.env, file, format = "fasta")
	if(with.plot)
	{
		#
		#	create and plot NJ tree on conc seq
		#			
		#	load outgroup sequences
		file			<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfg_HXB2outgroup.R')
		cat(paste('\nLoading outgroup seq from file', file))
		load(file)		#expect "outgroup.seq.gag" "outgroup.seq.pol" "outgroup.seq.env"
		#	concatenate sequences
		tmp				<- tolower(do.call('rbind',strsplit(df.seq[, paste(GAG,POL,ENV,sep='')],'')))
		rownames(tmp)	<- df.seq[, paste(IDCLU,treelabel.idx.sep,LABEL,sep='')]	
		seq				<- as.DNAbin(tmp)
		tmp				<- cbind(outgroup.seq.gag[,1:ncol(df.seq.gag)], outgroup.seq.pol, outgroup.seq.env)
		seq				<- rbind(seq,tmp)	
		seq.ph			<- nj(dist.dna(seq, model='raw'))		
		tmp				<- which(seq.ph$tip.label=="HXB2")
		seq.ph			<- reroot(seq.ph, tmp, seq.ph$edge.length[which(seq.ph$edge[,2]==tmp)])
		seq.ph			<- ladderize(seq.ph)
		#	plot
		file			<- paste(indir.sg, '/', infile.prefix, 'INFO_NJconc.pdf', sep='')				
		pdf(file=file, w=10, h=80)
		plot(seq.ph, show.tip=TRUE, cex=0.5)
		dev.off()		
	}
	#
	#	zip simulated files
	#
	tmp				<- c( paste(outdir, '/', infile.prefix, 'SIMULATED_metadata.csv', sep=''), paste(outdir, '/', infile.prefix, 'SIMULATED_env.fa', sep=''), paste(outdir, '/', infile.prefix, 'SIMULATED_pol.fa', sep=''), paste(outdir, '/', infile.prefix, 'SIMULATED_gag.fa', sep='') )
	zip( paste(outdir, '/', infile.prefix, 'SIMULATED.zip', sep=''), tmp, flags = "-FSr9XTj")
	dummy			<- file.remove(tmp)
	#
	#	zip internal files
	#
	tmp				<- list.files(outdir, pattern='*pdf$', recursive=TRUE, full.names=TRUE) 
	file.copy(tmp, outdir, overwrite=TRUE)	 
	tmp				<- c( paste(outdir, '/', infile.prefix, 'SIMULATED_INTERNAL.R', sep=''), list.files(outdir, pattern='*pdf$', recursive=FALSE, full.names=TRUE) ) 	
	zip( paste(outdir, '/', infile.prefix, 'INTERNAL.zip', sep=''), tmp, flags = "-FSr9XTj")
	dummy			<- file.remove(tmp)
	
	return(1)
}