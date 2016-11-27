##--------------------------------------------------------------------------------------------------------
##	olli 23.10.15
##--------------------------------------------------------------------------------------------------------
#
#	compute path differences on complete trees
#
treedist.pathdifference.wrapper<- function(df, ttrs, s, use.brl=FALSE, use.weight=FALSE)
{
	#tmp			<- subset(submitted.info, IDX==463)[1,]
	#IDX<- 463
	#IDX_T<-7
	tmp				<- df[, {
				cat('\nAt IDX', IDX)
				tidx		<- ifelse(use.brl, SUB_IDX_T, TIME_IDX_T)
				#IDX<- 241; TIME_IDX_T<- 6; IDX<- 1; TIME_IDX_T<- 1 
				stree		<- s[[IDX]]
				otree		<- ttrs[[tidx]]								
				z			<- setdiff(otree$tip.label, stree$tip.label)
				stopifnot( length(z)==abs(diff(c(Ntip(otree), Ntip(stree)))) )
				if(length(z))
					otree	<- drop.tip(otree, z)				
				z			<- path.dist(otree, stree, use.weight=use.weight)
				list(PD=z, NPD=z/choose(Ntip(otree),2), NPDSQ=z/sqrt(choose(Ntip(otree),2)), TAXA_NJ=Ntip(otree))
			}, by='IDX']
	tmp
}
#--------------------------------------------------------------------------------------------------------
#	MSE between true distances and patristic distances (units time) in reconstructed tree
#--------------------------------------------------------------------------------------------------------
treedist.MSE.wrapper<- function(df, s, tbrl, tinfo, use.brl=TRUE)
{		
	ans	<- df[, {
				#IDX<- 724; TIME_IDX_T<- 13; SUB_IDX_T<- 2
				#IDX<- 241; TIME_IDX_T<- 6; SUB_IDX_T<- 4
				cat('\nLSD distances IDX at', IDX)
				tidx	<- ifelse(use.brl, SUB_IDX_T, TIME_IDX_T)				
				stree	<- s[[IDX]]				
				#	mean squared error and mean absolute error of all pairwise distances
				tmp3	<- cophenetic.phylo(stree)
				#tmp3	<- distTips(stree, seq_len(Ntip(stree)), method='patristic', useC=TRUE)				
				#tmp3	<- as.matrix(tmp3)
				tmp3[upper.tri(tmp3, diag=TRUE)]	<- NA_real_
				tmp3	<- as.data.table(melt(tmp3))								
				setnames(tmp3, c('Var1','Var2','value'),c('TAXA1','TAXA2','PD_SIM'))
				tmp3	<- subset(tmp3, !is.na(PD_SIM))
				tmp2	<- subset(tbrl, IDX_T==tidx)
				set(tmp2,NULL,'TAXA1',tmp2[, as.character(TAXA1)])
				set(tmp2,NULL,'TAXA2',tmp2[, as.character(TAXA2)])
				set(tmp3,NULL,'TAXA1',tmp3[, as.character(TAXA1)])
				set(tmp3,NULL,'TAXA2',tmp3[, as.character(TAXA2)])
				tmp2	<- merge(tmp3, tmp2, by=c('TAXA1','TAXA2'))
				stopifnot(nrow(tmp2)==Ntip(stree)*(Ntip(stree)-1)/2)
				mse		<- tmp2[, mean((PD-PD_SIM)*(PD-PD_SIM))]
				mae		<- tmp2[, mean(abs(PD-PD_SIM))]
				#	mean squared error and mean absolute error of pairwise distances of sampled transmission pairs
				set(tmp2,NULL,'TAXA1',tmp2[, as.integer(gsub('IDPOP_','',gsub('\\|.*','',as.character(TAXA1))))])
				set(tmp2,NULL,'TAXA2',tmp2[, as.integer(gsub('IDPOP_','',gsub('\\|.*','',as.character(TAXA2))))])
				setnames(tmp2,c('TAXA1','TAXA2'),c('IDPOP','IDTR'))				
				tmp3	<- unique(subset(tinfo, IDX_T==tidx & IDTR_SAMPLED=='Y', select=c(IDPOP, IDTR)))
				tmp		<- copy(tmp3)
				setnames(tmp, c('IDPOP','IDTR'), c('IDTR','IDPOP'))
				tmp3	<- rbind(tmp3, tmp, use.names=TRUE)
				set(tmp3,NULL,'IDPOP',tmp3[, as.integer(gsub('IDPOP_','',IDPOP))])
				set(tmp3,NULL,'IDTR',tmp3[, as.integer(gsub('IDPOP_','',IDTR))])					
				tmp2	<- merge(tmp2,tmp3,by=c('IDPOP','IDTR'))
				mse.tp	<- tmp2[, mean((PD-PD_SIM)*(PD-PD_SIM))]
				mae.tp	<- tmp2[, mean(abs(PD-PD_SIM))]
				list(MSE=mse, MAE=mae, MSE_TP=mse.tp, MAE_TP=mae.tp, TAXA_NJ=Ntip(stree), EDGE_NJ=nrow(tmp2))				
			}, by='IDX']
	ans
}
#--------------------------------------------------------------------------------------------------------
#	MSE between true distances and patristic distances (units time) for each cluster in reconstructed tree
#--------------------------------------------------------------------------------------------------------
treedist.MSE.clusters.wrapper<- function(df, s, tbrl, tinfo, use.brl=TRUE)
{
	#
	setkey(tinfo, IDX_T)	
	ans		<- df[, {
				cat('\nAt IDX', IDX)
				#	IDX<- 724; TIME_IDX_T<- 13; SUB_IDX_T<- 2
				tidx		<- ifelse(use.brl, SUB_IDX_T, TIME_IDX_T)
				stree		<- s[[IDX]]								
				#	get all clusters of this true tree (with IDX_T) that are of size>=3 (use "tinfo" for that)
				z			<- subset(tinfo, CLU_N>3 & IDX_T==tidx)	
				setkey(z, TAXA)
				z			<- unique(z)				
				z			<- merge(z, data.table(TAXA=stree$tip.label, IN_STREE=1), by='TAXA', all.x=1)
				#	calculate the size of these clusters in the simulated tree
				z			<- merge(z, z[, list(CLU_NS= length(which(IN_STREE==1))), by='IDCLU'], by='IDCLU')
				#	get all clusters of size >= 3 in both the simulated and true tree
				z			<- subset(z, CLU_NS>3)
				#	precompute what can be precomputed before next loop
				tmp2	<- subset(tbrl, IDX_T==tidx)
				set(tmp2,NULL,'TAXA1',tmp2[, as.character(TAXA1)])
				set(tmp2,NULL,'TAXA2',tmp2[, as.character(TAXA2)])
				set(tmp2,NULL,'IDPOP',tmp2[, as.integer(gsub('IDPOP_','',gsub('\\|.*','',as.character(TAXA1))))])
				set(tmp2,NULL,'IDTR',tmp2[, as.integer(gsub('IDPOP_','',gsub('\\|.*','',as.character(TAXA2))))])
				#
				tmp4		<- unique(subset(tinfo, IDX_T==tidx & IDTR_SAMPLED=='Y', select=c(IDPOP, IDTR)))
				tmp			<- copy(tmp4)
				setnames(tmp, c('IDPOP','IDTR'), c('IDTR','IDPOP'))
				tmp4		<- rbind(tmp4, tmp, use.names=TRUE)	
				set(tmp4,NULL,'IDPOP',tmp4[, as.integer(gsub('IDPOP_','',IDPOP))])
				set(tmp4,NULL,'IDTR',tmp4[, as.integer(gsub('IDPOP_','',IDTR))])									
				#	if there any such clusters, calculate the quartet distance
				if(nrow(z))
				{
					#IDCLU	<- 3; TAXA	<- subset(z, IDCLU==3)[, TAXA]
					ans		<- z[, {								
								sclu	<- drop.tip(stree, setdiff(stree$tip.label,TAXA))								
								#	mean squared error of all pairwise distances
								tmp3	<- cophenetic.phylo(sclu)
								tmp3[upper.tri(tmp3, diag=TRUE)]	<- NA_real_
								tmp3	<- as.data.table(melt(tmp3))								
								setnames(tmp3, c('Var1','Var2','value'),c('TAXA1','TAXA2','PD_SIM'))
								tmp3	<- subset(tmp3, !is.na(PD_SIM))
								set(tmp3,NULL,'TAXA1',tmp3[, as.character(TAXA1)])
								set(tmp3,NULL,'TAXA2',tmp3[, as.character(TAXA2)])								
								#	this merges to the intersection of taxa in sclu and the corresponding observed clu
								tmp3	<- merge(tmp3, tmp2, by=c('TAXA1','TAXA2'))								
								mse		<- tmp3[, mean((PD-PD_SIM)*(PD-PD_SIM))]
								mae		<- tmp3[, mean(abs(PD-PD_SIM))]
								#	mean squared error and mean absolute error of pairwise distances of sampled transmission pairs																				
								tmp3	<- merge(tmp3,tmp4,by=c('IDPOP','IDTR'))
								mse.tp	<- tmp3[, mean((PD-PD_SIM)*(PD-PD_SIM))]
								mae.tp	<- tmp3[, mean(abs(PD-PD_SIM))]
								list(MSE=mse, MAE=mae, MSE_TP=mse.tp, MAE_TP=mae.tp, TAXA_NC=Ntip(sclu), EDGE_NC=nrow(tmp2))								
							}, by='IDCLU']	
				}
				if(!nrow(z))
					ans		<- data.table(MSE=NA_real_, MAE=NA_real_, MSE_TP=NA_real_, MAE_TP=NA_real_, TAXA_NC=NA_integer_, EDGE_NC=NA_integer_)
				ans			
			}, by='IDX']	
	ans	
}
#--------------------------------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------------------------------
treedist.pathdifference.clusters.wrapper<- function(df, ttrs, s, tinfo, use.brl=TRUE, use.weight=FALSE)
{
	#
	setkey(tinfo, IDX_T)
	tmp		<- df[, {
				cat('\nAt IDX', IDX)
				#	IDX<- 1; SUB_IDX_T<- 1
				stree		<- s[[IDX]]
				tidx		<- ifelse(use.brl, SUB_IDX_T, TIME_IDX_T)
				otree		<- ttrs[[tidx]]				
				#	get all clusters of this true tree (with IDX_T) that are of size>=3 (use "tinfo" for that)
				z			<- subset(tinfo, CLU_N>3 & IDX_T==tidx)	
				setkey(z, TAXA)
				z			<- unique(z)				
				z			<- merge(z, data.table(TAXA=stree$tip.label, IN_STREE=1), by='TAXA', all.x=1)
				#	calculate the size of these clusters in the simulated tree
				z			<- merge(z, z[, list(CLU_NS= length(which(IN_STREE==1))), by='IDCLU'], by='IDCLU')
				#	get all clusters of size >= 3 in both the simulated and true tree
				z			<- subset(z, CLU_NS>3)
				#	if there any such clusters, calculate the quartet distance
				if(nrow(z))
				{
					#IDCLU	<- 6; TAXA	<- subset(z, IDCLU==6)[, TAXA]
					ans		<- z[, {								
								sclu	<- drop.tip(stree, setdiff(stree$tip.label,TAXA))
								oclu	<- drop.tip(otree, union( setdiff(otree$tip.label, stree$tip.label), setdiff(otree$tip.label,TAXA)))
								z		<- path.dist(oclu, sclu, use.weight=use.weight)
								list(PD=z, NPD=z/choose(Ntip(oclu),2), NPDSQ=z/sqrt(choose(Ntip(oclu),2)), TAXA_NC=Ntip(oclu))								
							}, by='IDCLU']	
				}
				if(!nrow(z))
					ans		<- data.table(PD=NA_real_, NPD=NA_real_, NPDSQ=NA_real_, TAXA_NC=NA_integer_)
				ans			
			}, by='IDX']	
	tmp	
}
##--------------------------------------------------------------------------------------------------------
##	olli 01.08.16	compute quartet differences on complete trees
##--------------------------------------------------------------------------------------------------------
treedist.quartetdifference.wrapper<- function(submitted.info, ttrs, strs_rtt)
{
	tmp				<- submitted.info[, {
				cat('\nAt IDX', IDX)
				#IDX<- 241; TIME_IDX_T<- 6
				stree		<- unroot(strs_rtt[[IDX]])
				otree		<- unroot(ttrs[[TIME_IDX_T]])								
				#print(stree)
				#print(otree)
				z			<- setdiff(otree$tip.label, stree$tip.label)
				stopifnot( length(z)==abs(diff(c(Ntip(otree), Ntip(stree)))) )
				if(length(z))
					otree	<- unroot(drop.tip(otree, z))				
				z			<- quartets.distance.cmd(otree, stree)
				list(TAXA_NJ=z['TAXA_NJ'], NQD=z['NQD'])
			}, by='IDX']
	tmp
}
#--------------------------------------------------------------------------------------------------------
#	olli 01.08.16	compute quartet differences on sub trees
#--------------------------------------------------------------------------------------------------------
treedist.quartetdifference.clusters.wrapper<- function(submitted.info, ttrs, strs_rtt, tinfo)
{
	setkey(tinfo, IDX_T)
	tmp		<- subset(submitted.info, MODEL=='R')[, {
				cat('\nAt IDX', IDX)
				#	IDX<- 1; SUB_IDX_T<- 1
				stree		<- strs_rtt[[IDX]]
				otree		<- ttrs[[SUB_IDX_T]]				
				#	get all clusters of this true tree (with IDX_T) that are of size>=3 (use "tinfo" for that)
				z			<- subset(tinfo, CLU_N>3 & IDX_T==SUB_IDX_T)	
				setkey(z, TAXA)
				z			<- unique(z)				
				z			<- merge(z, data.table(TAXA=stree$tip.label, IN_STREE=1), by='TAXA', all.x=1)
				#	calculate the size of these clusters in the simulated tree
				z			<- merge(z, z[, list(CLU_NS= length(which(IN_STREE==1))), by='IDCLU'], by='IDCLU')
				#	get all clusters of size >= 3 in both the simulated and true tree
				z			<- subset(z, CLU_NS>3)
				#	if there any such clusters, calculate the quartet distance
				if(nrow(z))
				{
					#IDCLU	<- 6; TAXA	<- subset(z, IDCLU==6)[, TAXA]
					ans		<- z[, {								
								sclu	<- unroot(drop.tip(stree, setdiff(stree$tip.label,TAXA)))
								oclu	<- unroot(drop.tip(otree, union( setdiff(otree$tip.label, stree$tip.label), setdiff(otree$tip.label,TAXA))))
								z		<- quartets.distance.cmd(oclu, sclu)
								list(NQDC=z['NQD'], TAXA_NC=Ntip(oclu))
							}, by='IDCLU']	
				}
				if(!nrow(z))
					ans		<- data.table(IDCLU=NA_integer_, NQDC=NA_real_, TAXA_NC=NA_integer_)
				ans			
			}, by='IDX']	
	tmp			
}
#--------------------------------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------------------------------
project.PANGEA.visualize.call.patterns.align150623<- function()
{
	outdir			<- '/Users/Oliver/Dropbox (Infectious Disease)/2016_PANGEA_treecomp/figures'
	min.coverage	<- 600
	min.depth		<- 10
	#with.gaps		<- 1
	#outfile			<- '150623_PANGEAGlobal2681_C5_wgaps.pdf'
	with.gaps		<- 0
	outfile			<- '150623_PANGEAGlobal2681_C10.pdf'
	
	
	infile			<- '/Users/Oliver/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/PANGEA_150623/150623_PANGEAGlobal2681_C10.fa'
	so				<- read.dna(infile, format='fasta')
	rownames(so)	<- gsub('-','_',rownames(so)) 
	#	drop LTR from sequences
	#	need to load latest alignment, get LTR start, and translate into one common sequence
	infile			<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/alignments_160110/PANGEA_HIV_n5003_Imperial_v160110_GlobalAlignment.rda'
	load(infile)	#loads sqi, sq
	tmp				<- sq[which(grepl('HXB2',rownames(sq))),]
	pattern			<- 'a-*t-*g-*g-*g-*t-*g-*c-*g-*a-*g-*a-*g-*c-*g-*t-*c-*a'
	idx.st.in.sq	<- regexpr(pattern, paste(as.character( tmp ),collapse=''))
	tx.common		<- intersect( rownames(so) , rownames(sq) )[1]
	pattern			<- paste(as.character( sq[tx.common, idx.st.in.sq:(idx.st.in.sq+10)]),collapse='')	
	idx.st.in.so	<- regexpr(pattern, paste(as.character(so[tx.common,]), collapse=''))	
	so				<- so[, seq.int(idx.st.in.so, ncol(so))]
	#	drop sites outside the HIV reference compendium
	#tx.common		<- tail( intersect( rownames(so) , rownames(sq) ), 1 )
	tx.common		<- 'PG14_ZA100095_S01002'
	tmp				<- gsub('-|?','', paste( as.character(sq[tx.common,]), collapse='') )	
	pattern			<- substring(tmp, nchar(tmp)-20, nchar(tmp) )
	pattern			<- paste(strsplit(pattern, '')[[1]], collapse='-*')
	tmp				<- so[tx.common, ]
	idx.end.in.so	<- regexpr(pattern, paste(as.character( tmp ),collapse=''))
	idx.end.in.so	<- as.integer(idx.end.in.so+attr(idx.end.in.so,'match.length')-1L)
	so				<- so[, seq.int(1, idx.end.in.so)]
	# 	drop gap only columns
	if(!with.gaps)
	{
		tmp			<- apply( as.character(so), 2, function(x) !all(x%in%c('?','-','n')) ) 
		so			<- so[, tmp]		
	}
	#	convert into chunks
	ch				<- lapply(seq_len(nrow(so)), function(i)
			{
				z	<- gregexpr('1+', paste(as.numeric( !as.character( so[i,] )%in%c('-','?','n') ), collapse='') )[[1]]
				data.table(PANGEA_ID= rownames(so)[i], POS=as.integer(z), DEPTH=min.depth, REP=attr(z,"match.length"))
			})
	ch				<- do.call('rbind',ch)	
	set(ch, NULL, 'SITE', ch[, regmatches(PANGEA_ID, regexpr('_[A-Z]+',PANGEA_ID))])
	set(ch, NULL, 'SITE', ch[, substring(SITE,2)])
	ch				<- merge(ch, ch[, list(COV=sum(REP)), by='PANGEA_ID'], by='PANGEA_ID')
	#	select min.coverage, select min.depth
	ch			<- subset(ch, !is.na(SITE) & COV>=min.coverage & DEPTH>=min.depth)
	#	define chunks
	ch[, POS_NEXT:= POS+REP]	
	ch		<- ch[, list(SITE=SITE, POS=POS, DEPTH=DEPTH, REP=REP, CHUNK=cumsum(as.numeric(c(TRUE, POS[-1]!=POS_NEXT[-length(POS_NEXT)])))), by='PANGEA_ID']
	ch		<- ch[, list(SITE=SITE[1], POS_CH=min(POS), REP_CH=sum(REP), DEPTH_CH= sum(DEPTH*REP)/sum(REP) ), by=c('PANGEA_ID','CHUNK')]
	ch[, DEPTH_MIN:=min.depth]
	set(ch, NULL, 'SITE', ch[, factor(SITE, levels=c('BW', 'ZA', 'UG'), labels=c('Botswana', 'South Africa', 'Uganda'))])
	ch				<- merge(ch, ch[, list(COV=sum(REP_CH)), by='PANGEA_ID'], by='PANGEA_ID')
	ch[, COVP:= COV/ncol(so)]	
	#	
	setkey(ch, PANGEA_ID)	
	dcast.data.table(unique(ch)[, list(P=seq(0,1,0.1), Q=quantile(COV, p=seq(0,1,0.1))), by='SITE'], P~SITE, value.var='Q')
	dcast.data.table(unique(ch)[, list(Q=seq(1e3,8e3,1e3), P=ecdf(COV)(seq(1e3,8e3,1e3))), by='SITE'], Q~SITE, value.var='P')
	#	proportion of non-gaps in alignment
	#unique(ch)[, sum(COV)] / (nrow(unique(ch))*ncol(so))
	nrow(unique(ch))
	#	C5: 	2367
	#	C10: 	2126
	unique(ch)[, mean(COVP)]
	#	C5: 	0.6201531
	#	C10:	0.6510485
	unique(ch)[, mean(COV)]
	#	C5: 	5472.231
	#	C10:	5744.852
	
	#	plot chunks
	require(viridis)
	ch		<- merge(ch, ch[, list(POS_CHF=min(POS_CH)), by='PANGEA_ID'], by='PANGEA_ID')	
	setkey(ch, SITE, PANGEA_ID)	
	tmp		<- unique(ch)
	setkey(tmp, POS_CHF)	
	if(min.depth==5)
		tmp[, PLOT:=ceiling(seq_len(nrow(tmp))/600)]
	if(min.depth==10)
		tmp[, PLOT:=ceiling(seq_len(nrow(tmp))/535)]	
	ch		<- merge(ch, subset(tmp, select=c(PANGEA_ID, PLOT)), by='PANGEA_ID')
	set(ch, NULL, 'PLOT', ch[, factor(PLOT, levels=c(4,3,2,1), labels=c(4,3,2,1))])
	setkey(ch, POS_CH, SITE)
	set(ch, NULL, 'PANGEA_ID', ch[, factor(PANGEA_ID, levels=unique(PANGEA_ID), labels=unique(PANGEA_ID))])
	ggplot(ch, aes(y=PANGEA_ID, yend=PANGEA_ID, x=POS_CH, xend=POS_CH+REP_CH-1L, colour=SITE)) +
			scale_x_continuous(expand=c(0,0), breaks=seq(0,10e3,1e3), minor_breaks=seq(0,10e3,100)) +
			scale_colour_manual(values=c('Botswana'="#1B0C42FF", 'South Africa'="#CF4446FF", 'Uganda'="#781C6DFF")) +			
			geom_segment() + theme_bw() + 
			facet_wrap(~PLOT, scales='free_y', ncol=4) +
			labs(x='\nalignment position', y='PANGEA-HIV sequences\n', colour='sampling\nlocation') + 			
			theme(	axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank(), legend.position='bottom',
					strip.text= element_blank(), strip.background=element_blank()) +
			guides(colour=guide_legend(override.aes=list(size=5)))	
	ggsave(file=file.path(outdir,outfile), w=15, h=10, limitsize = FALSE) 
}
#--------------------------------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------------------------------
project.PANGEA.visualize.call.patterns.align160110<- function()
{
	min.coverage	<- 600
	min.depth		<- 10
	
	infile			<- '/Users/Oliver/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/readlengths/bam_stats_150218.rda'
	load(infile)
	setnames(bam.cov, c('FILE_ID','COV'), c('SANGER_ID','DEPTH'))	
	wdir			<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/explaingaps'
	wfile			<- 'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.rda'
	load(file.path(wdir,wfile))	#loads sqi, sq, dm 
	outdir			<- '~/Dropbox (Infectious Disease)/2016_PANGEA_treecomp/figures'	
	#
	#	convert sqp into chunks
	#
	ch				<- lapply(seq_len(nrow(sqp)), function(i)
			{
				z	<- gregexpr('1+', paste(as.numeric( !as.character( sqp[i,] )%in%c('-','?','n') ), collapse='') )[[1]]
				data.table(PANGEA_ID= rownames(sqp)[i], POS=as.integer(z), DEPTH=min.depth, REP=attr(z,"match.length"))
			})
	ch				<- do.call('rbind',ch)	
	#	define SITE
	tmp				<- subset(dm, select=c(TAXA, COHORT))
	setnames(tmp, c('TAXA','COHORT'),c('PANGEA_ID','SITE'))
	ch				<- merge(ch, tmp, by='PANGEA_ID')
	#	get coverage		
	ch				<- merge(ch, ch[, list(COV=sum(REP)), by='PANGEA_ID'], by='PANGEA_ID')
	#	select min.coverage, select min.depth
	ch			<- subset(ch, COV>=min.coverage & DEPTH>=min.depth)
	#	define chunks
	ch[, POS_NEXT:= POS+REP]	
	ch		<- ch[, list(SITE=SITE, POS=POS, DEPTH=DEPTH, REP=REP, CHUNK=cumsum(as.numeric(c(TRUE, POS[-1]!=POS_NEXT[-length(POS_NEXT)])))), by='PANGEA_ID']
	ch		<- ch[, list(SITE=SITE[1], POS_CH=min(POS), REP_CH=sum(REP), DEPTH_CH= sum(DEPTH*REP)/sum(REP) ), by=c('PANGEA_ID','CHUNK')]
	ch[, DEPTH_MIN:=min.depth]
	set(ch, NULL, 'SITE', ch[, factor(SITE, levels=c("AC_Resistance","BW-Mochudi","UG-MRC","RCCS"), labels=c('Africa Centre\n(Resistance Cohort)', 'Botswana\n(Mochudi)', 'Uganda-\nMRC','Rakai Community\nCohort Study'))])	
	ch				<- merge(ch, ch[, list(COV=sum(REP_CH)), by='PANGEA_ID'], by='PANGEA_ID')
	ch[, COVP:= COV/ncol(sq)]	
	ch		<- merge(ch, ch[, list(POS_CHF=min(POS_CH)), by='PANGEA_ID'], by='PANGEA_ID')
	#		
	#	plot chunks
	#
	require(viridis)		
	setkey(ch, SITE, PANGEA_ID)	
	tmp		<- unique(ch)
	setkey(tmp, POS_CHF)	
	tmp[, PLOT:=ceiling(seq_len(nrow(tmp))/ceiling(nrow(tmp)/4))]	
	ch		<- merge(ch, subset(tmp, select=c(PANGEA_ID, PLOT)), by='PANGEA_ID')
	set(ch, NULL, 'PLOT', ch[, factor(PLOT, levels=c(4,3,2,1), labels=c(4,3,2,1))])
	setkey(ch, POS_CH, SITE)
	set(ch, NULL, 'PANGEA_ID', ch[, factor(PANGEA_ID, levels=unique(PANGEA_ID), labels=unique(PANGEA_ID))])
	
	dpani	<- subset(dpan, !is.na(START))[, list(START=START[1], END=START[1]+max(IDX)-1L), by='PR']
	dgeni	<- subset(dgene, GENE%in%c('GAG','POL','ENV'))
	ggplot(ch) +
			scale_x_continuous(expand=c(0,0), breaks=seq(0,10e3,1e3), minor_breaks=seq(0,10e3,100)) +
			scale_colour_manual(values=c('Botswana\n(Mochudi)'="#33638DFF", 'Africa Centre\n(Resistance Cohort)'="#CF4446FF", 'Rakai Community\nCohort Study'="#A8327DFF", 'Uganda-\nMRC'="#29AF7FFF")) +			
			geom_segment(aes(y=PANGEA_ID, yend=PANGEA_ID, x=POS_CH, xend=POS_CH+REP_CH-1L, colour=SITE)) +
			geom_rect(data=dpani, aes(xmin=START, xmax=END, ymin=-Inf, ymax=Inf), fill="black", size=5) +
			geom_text(data=dpani, aes(x=START, y=seq_len(length(START))*10+10, label=PR), colour="black", hjust=-.2, size=2) +
			geom_vline(data=dgeni, aes(xintercept=START)) +
			geom_vline(data=dgeni, aes(xintercept=END)) +
			geom_text(data=dgeni, aes(x=START, y=seq_len(length(START))*10+10, label=GENE), colour="blue", hjust=-.2, size=2) +
			geom_text(data=dgeni, aes(x=END, y=seq_len(length(END))*10+10, label=GENE), colour="blue", hjust=-.2, size=2) +
			theme_bw() + 
			facet_wrap(~PLOT, scales='free_y', ncol=4) +
			labs(x='\nalignment position', y='PANGEA-HIV sequences\n', colour='sampling\nlocation') + 			
			theme(	axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank(), legend.position='bottom',
					strip.text= element_blank(), strip.background=element_blank()) +
			guides(colour=guide_legend(override.aes=list(size=5)))	
	ggsave(file=file.path(outdir,'PANGEA_HIV_n5003_Imperial_v160110_selectedchunks_160905.pdf'), w=15, h=15, limitsize=FALSE)
	#
	#	for comparison, ATHENA data set
	#
	infile		<- "/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013/data/tmp/ATHENA_2013_03_NoDRAll+LANL_Sequences_Thu_Aug_01_17:05:23_2013.fasta"
	sq			<- read.dna(infile, format='fa')
	sq			<- sq[1:1293,]
	sqi			<- data.table(FASTA= rownames(sq))
	sqi[, FRGN:= grepl('^TN', FASTA) | grepl('PROT+P51', FASTA)]
	sq			<- sq[ subset(sqi, !FRGN)[, FASTA], ]
	#	remove all gap cols
	sq			<- as.character(sq)
	tmp			<- apply(sq, 2, function(x) !all(x%in%c('-','n','?')))
	sq			<- sq[,tmp]
	#	count coverage
	sqi			<- subset(sqi, !FRGN)[, list(COV=length(which( !sq[FASTA,]%in%c('-','n','?') ))), by='FASTA']
	sqi[, COVP:= COV/ncol(sq)]
	sqi[, mean(COVP)]
	#	0.96
}
##--------------------------------------------------------------------------------------------------------
##	olli 01.11.15
##	from https://www.cs.upc.edu/~valiente/comput-biol/
##--------------------------------------------------------------------------------------------------------
quartet <- function (t,i,j,k,l) 
{
	unroot(drop.tip(t,setdiff(t$tip.label,c(i,j,k,l))))
}
quartets.distance <- function (t1,t2) {
	L <- sort(t1$tip.label)
	d <- 0
	for (i in L[1:(length(L)-3)]) {
		for (j in L[(match(i,L)+1):(length(L)-2)]) {
			for (k in L[(match(j,L)+1):(length(L)-1)]) {
				for (l in L[(match(k,L)+1):length(L)]) {
					q1 <- quartet(t1,i,j,k,l)
					q2 <- quartet(t2,i,j,k,l)
					if (dist.topo(q1,q2) != 0) { d <- d + 2 }
				}
			}
		}
	}
	d
}
quartets.distance.cmd	<- function(otree, stree, PROG.QDIST='/apps/qdist/2.0/bin/qdist')
{
	file1	<- paste( getwd(), '/tree1', sep='')
	file2	<- paste( getwd(), '/tree2', sep='')
	write.tree(otree,file=file1)
	write.tree(stree,file=file2)
	cmd<- paste(PROG.QDIST,file1,file2)
	tmp		<- system(cmd, intern=TRUE)
	c('TAXA_NJ'=Ntip(otree), 'NQD'=as.numeric(tail(unlist(strsplit(tmp[2],'\t')),1)))
}
treedist.get.tree.100bs<- function(phr, phs, tmpdir)
{
	require(ggtree)
	require(phangorn)	
	#	write to tmp directory
	write.tree(phr, file=file.path(tmpdir,'phr.newick'))
	invisible(sapply(seq_along(phs), function(i){		write.tree(phs[[i]], file=file.path(tmpdir, paste('phs',i,'.newick', sep='')))		}))
	cmd				<- paste('for i in $(seq 1 ',length(phs),'); do cat ',file.path(tmpdir,'phs$i.newick'),' >> ',file.path(tmpdir,'phc.newick'),'; done', sep='')
	system(cmd)
	cmd				<- paste('cd ',tmpdir,'\nraxmlHPC-AVX -f b -t phr.newick -z phc.newick -m GTRCAT -s alg -n TEST', sep='')
	system(cmd)
	#	prepare RAXML tree for plotting
	raxml 	<- read.raxml(file.path(tmpdir,'RAxML_bipartitionsBranchLabels.TEST'))
	pho		<- raxml@phylo
	tmp		<- as.data.table(raxml@bootstrap)
	z		<- subset(tmp, bootstrap==100)[, {
				#node	<- 1597
				desc	<- Descendants(pho, node, type='all')
				if( any(desc%in%subset(tmp, bootstrap<100)[, node]) )
					desc<- NA_integer_
				list(desc=desc)
			}, by='node']
	z		<- subset(z, !is.na(desc))
	z[, bootstrap:=100]
	z[, node:=NULL]
	setnames(z, 'desc', 'node')
	z		<- rbind(z, data.table(node= setdiff(seq.int(Nnode(pho, internal.only=FALSE)), z[, node]), bootstrap=0))
	raxml@bootstrap	<- as.data.frame(z)
	#	delete stuff
	invisible(file.remove(list.files(tmpdir, full.names=TRUE)))
	raxml
}
treedist.quartets.add<- function(submitted.info=NULL, ttrs=NULL, strs=NULL, file=NULL, with.save=0)
{
	#	file<- '/work/or105/Gates_2014/tree_comparison/submitted_151101.rda'
	require(ape)
	require(data.table)	
	stopifnot( !is.null(submitted.info) || !is.null(file))	
	if(is.null(submitted.info))
	{
		load(file)
		with.save	<- 1	
	}		
	#tmp			<- subset(submitted.info, IDX==463)[1,]
	#IDX<- 463
	#IDX_T<-7
	tmp				<- submitted.info[, {
				cat('\nAt IDX', IDX)
				#IDX<- 241; TIME_IDX_T<- 6
				stree		<- unroot(strs_rtt[[IDX]])
				otree		<- unroot(ttrs[[TIME_IDX_T]])								
				#print(stree)
				#print(otree)
				z			<- setdiff(otree$tip.label, stree$tip.label)
				stopifnot( length(z)==abs(diff(c(Ntip(otree), Ntip(stree)))) )
				if(length(z))
					otree	<- unroot(drop.tip(otree, z))				
				z			<- quartets.distance.cmd(otree, stree)
				list(TAXA_NJ=z['TAXA_NJ'], NQD=z['NQD'])
			}, by='IDX']
	submitted.info	<- merge(submitted.info, tmp, by='IDX')
	#
	setkey(tinfo, IDX_T)
	tmp		<- subset(submitted.info, MODEL=='R')[, {
				cat('\nAt IDX', IDX)
				#	IDX<- 1; SUB_IDX_T<- 1
				stree		<- strs_rtt[[IDX]]
				otree		<- ttrs[[SUB_IDX_T]]				
				#	get all clusters of this true tree (with IDX_T) that are of size>=3 (use "tinfo" for that)
				z			<- subset(tinfo, CLU_N>3 & IDX_T==SUB_IDX_T)	
				setkey(z, TAXA)
				z			<- unique(z)				
				z			<- merge(z, data.table(TAXA=stree$tip.label, IN_STREE=1), by='TAXA', all.x=1)
				#	calculate the size of these clusters in the simulated tree
				z			<- merge(z, z[, list(CLU_NS= length(which(IN_STREE==1))), by='IDCLU'], by='IDCLU')
				#	get all clusters of size >= 3 in both the simulated and true tree
				z			<- subset(z, CLU_NS>3)
				#	if there any such clusters, calculate the quartet distance
				if(nrow(z))
				{
					#IDCLU	<- 6; TAXA	<- subset(z, IDCLU==6)[, TAXA]
					ans		<- z[, {								
								sclu	<- unroot(drop.tip(stree, setdiff(stree$tip.label,TAXA)))
								oclu	<- unroot(drop.tip(otree, union( setdiff(otree$tip.label, stree$tip.label), setdiff(otree$tip.label,TAXA))))
								z		<- quartets.distance.cmd(oclu, sclu)
								list(NQDC=z['NQD'])
							}, by='IDCLU']	
				}
				if(!nrow(z))
					ans		<- data.table(IDCLU=NA_integer_, NQDC=NA_real_)
				ans			
			}, by='IDX']	
	sclu.info	<- merge(sclu.info, tmp, by=c('IDX','IDCLU'))		
	
	if(with.save)
		save(strs, strs_rtt, ttrs, tinfo, tfiles, submitted.info, sclu.info, lba, file=gsub('\\.rda','_QD\\.rda',file))
}
##--------------------------------------------------------------------------------------------------------
##	olli 13.07.16
##--------------------------------------------------------------------------------------------------------
treedist.closest.ind.obs<- function(tinfo, gd.thresh=Inf, rtn.pairs=FALSE)
{
	tmp				<- subset(tinfo, BRL_T=='subst')[, {
				#z<- subset(tinfo, BRL_T=='subst' & IDX_T==1); IDPOP<- z$IDPOP; IDTR<- z$IDTR; IDREC<- z$IDREC; IDPOP_CL<- z$IDPOP_CL; IDPOP_CL_GD<- z$IDPOP_CL_GD
				gds	<- data.table(IDPOP=as.integer(gsub('IDPOP_','',IDPOP)), IDTR= as.integer(IDTR), IDREC=IDREC, IDPOP_CL=as.integer(IDPOP_CL), IDPOP_CL_GD=IDPOP_CL_GD)
				gds[, CL_PH_PAIR:= gds[, as.character(as.numeric(IDPOP<IDPOP_CL))]]
				tmp		<- gds[, which(CL_PH_PAIR=='1')]
				set(gds, tmp, 'CL_PH_PAIR', gds[tmp, paste(IDPOP,IDPOP_CL,sep=',')])
				tmp		<- gds[, which(CL_PH_PAIR=='0')]
				set(gds, tmp, 'CL_PH_PAIR', gds[tmp, paste(IDPOP_CL,IDPOP,sep=',')])
				setkey(gds, CL_PH_PAIR)
				ans		<- subset(gds, IDPOP_CL_GD<=gd.thresh)				
				#	define tr->POP pairs
				ans[, TR_PAIR:= ans[, as.character(as.numeric(IDPOP<IDTR))]]
				tmp		<- ans[, which(TR_PAIR=='1')]
				set(ans, tmp, 'TR_PAIR', ans[tmp, paste(IDPOP,IDTR,sep=',')])
				tmp		<- ans[, which(TR_PAIR=='0')]
				set(ans, tmp, 'TR_PAIR', ans[tmp, paste(IDTR,IDPOP,sep=',')])
				#	define POP->rec pairs
				ans[, REC_PAIR:= ans[, as.character(as.numeric(IDPOP<IDREC))]]
				tmp		<- ans[, which(REC_PAIR=='1')]
				set(ans, tmp, 'REC_PAIR', ans[tmp, paste(IDPOP,IDREC,sep=',')])
				tmp		<- ans[, which(REC_PAIR=='0')]
				set(ans, tmp, 'REC_PAIR', ans[tmp, paste(IDREC,IDPOP,sep=',')])
				#	determine matches
				ans[, IN:= CL_PH_PAIR==TR_PAIR]
				tmp		<- ans[, which(!is.na(IDREC) & !IN)]
				set(ans, tmp, 'IN', ans[tmp, CL_PH_PAIR==REC_PAIR])
				ans		<- ans[, list(TRUE_PAIR=any(IN)), by='CL_PH_PAIR']
				if(!rtn.pairs)
				{
					#	calculate proportion of phylog closest pairs that are true transmission pairs					
					#ans		<- list(TR_REC_perc_T= ans[, mean(as.numeric(TRUE_PAIR))]  )
					ans		<- list(TPAIR_PHCL_T= ans[, length(which(TRUE_PAIR))], NTPAIR_PHCL_T= ans[, length(which(!TRUE_PAIR))])
				}
				ans
			}, by='IDX_T']
	setnames(tmp, 'IDX_T','SUB_IDX_T')
	tmp
}
##--------------------------------------------------------------------------------------------------------
##	olli 13.07.16
##--------------------------------------------------------------------------------------------------------
treedist.closest.ind.reconstructed<- function(submitted.info, tinfo, strs, gd.thresh)
{
	sucl			<- subset(submitted.info, MODEL=='R')[, {
				print(IDX)
				#IDX<- 557; SUB_IDX_T<-2; SC<- '150701_REGIONAL_TRAIN2'
				ph			<- strs[[IDX]]
				model.reg	<- grepl('REGIONAL',SC)
				gds			<- treedist.closest.ind(ph, model.reg)
				gds			<- subset(gds, IDPOP_CL_GD<=gd.thresh)
				ans			<- list(TPAIR_PHCL=NA_integer_, NTPAIR_PHCL=NA_integer_)
				if(nrow(gds)>0)
				{
					tmp			<- subset(tinfo, IDX_T==SUB_IDX_T, c(IDPOP, IDTR, IDREC))
					ans			<- merge(gds, tmp, by='IDPOP')
					set(ans, NULL, 'IDPOP', ans[, as.integer(gsub('IDPOP_','',IDPOP))])
					set(ans, NULL, 'IDPOP_CL', ans[, as.integer(gsub('IDPOP_','',IDPOP_CL))])
					set(ans, NULL, 'IDTR', ans[, as.integer(IDTR)])
					set(ans, NULL, 'IDREC', ans[, as.integer(IDREC)])
					#	get phylogenetically closest pairs			
					ans[, CL_PH_PAIR:= ans[, as.character(as.numeric(IDPOP<IDPOP_CL))]]
					tmp		<- ans[, which(CL_PH_PAIR=='1')]
					set(ans, tmp, 'CL_PH_PAIR', ans[tmp, paste(IDPOP,IDPOP_CL,sep=',')])
					tmp		<- ans[, which(CL_PH_PAIR=='0')]
					set(ans, tmp, 'CL_PH_PAIR', ans[tmp, paste(IDPOP_CL,IDPOP,sep=',')])
					setkey(ans, CL_PH_PAIR)
					#	define tr->POP pairs
					ans[, TR_PAIR:= ans[, as.character(as.numeric(IDPOP<IDTR))]]
					tmp		<- ans[, which(TR_PAIR=='1')]
					set(ans, tmp, 'TR_PAIR', ans[tmp, paste(IDPOP,IDTR,sep=',')])
					tmp		<- ans[, which(TR_PAIR=='0')]
					set(ans, tmp, 'TR_PAIR', ans[tmp, paste(IDTR,IDPOP,sep=',')])
					#	define POP->rec pairs
					ans[, REC_PAIR:= ans[, as.character(as.numeric(IDPOP<IDREC))]]
					tmp		<- ans[, which(REC_PAIR=='1')]
					set(ans, tmp, 'REC_PAIR', ans[tmp, paste(IDPOP,IDREC,sep=',')])
					tmp		<- ans[, which(REC_PAIR=='0')]
					set(ans, tmp, 'REC_PAIR', ans[tmp, paste(IDREC,IDPOP,sep=',')])
					#	determine matches
					ans[, IN:= CL_PH_PAIR==TR_PAIR]
					tmp		<- ans[, which(!is.na(IDREC) & !IN)]
					set(ans, tmp, 'IN', ans[tmp, CL_PH_PAIR==REC_PAIR])
					#	calculate proportion of phylog closest pairs that are true transmission pairs
					ans		<- ans[, list(TRUE_PAIR=any(IN)), by='CL_PH_PAIR']
					ans		<- list(TPAIR_PHCL= ans[, length(which(TRUE_PAIR))], NTPAIR_PHCL= ans[, length(which(!TRUE_PAIR))])					 
				}
				ans				
			}, by=c('IDX')]
	sucl
}
##--------------------------------------------------------------------------------------------------------
##	olli 13.07.16
##--------------------------------------------------------------------------------------------------------
treedist.closest.ind.reconstructed.oftruepairs<- function(submitted.info, tinfo.pairs, strs)
{
	tmp			<- subset(submitted.info, MODEL=='R')[, {
				print(IDX)
				#IDX<- 557; SUB_IDX_T<-2; SC<- '150701_REGIONAL_TRAIN2'
				ph			<- strs[[IDX]]
				model.reg	<- grepl('REGIONAL',SC)
				gds			<- treedist.closest.ind(ph, model.reg)
				gds			<- subset(gds, IDPOP_CL_GD<=Inf)
				#	get correct phylogenetically closest pairs
				z			<- SUB_IDX_T
				z			<- subset(tinfo.pairs, SUB_IDX_T==z)		
				z			<- subset(z, TRUE_PAIR_Inf)
				#	get phylogenetically closest pairs in simulation
				ans			<- copy(gds)
				set(ans, NULL, 'IDPOP', ans[, as.integer(gsub('IDPOP_','',IDPOP))])
				set(ans, NULL, 'IDPOP_CL', ans[, as.integer(gsub('IDPOP_','',IDPOP_CL))])								
				ans[, CL_PH_PAIR:= ans[, as.character(as.numeric(IDPOP<IDPOP_CL))]]
				tmp		<- ans[, which(CL_PH_PAIR=='1')]
				set(ans, tmp, 'CL_PH_PAIR', ans[tmp, paste(IDPOP,IDPOP_CL,sep=',')])
				tmp		<- ans[, which(CL_PH_PAIR=='0')]
				set(ans, tmp, 'CL_PH_PAIR', ans[tmp, paste(IDPOP_CL,IDPOP,sep=',')])
				setkey(ans, CL_PH_PAIR)					
				#	calculate which of the correct phylogenetically closest pairs are also in the simulation
				ans		<- merge(z, unique(ans), by='CL_PH_PAIR', all.x=1)				
				list(TR_PAIR_rec= ans[, mean(!is.na(IDPOP))]  )									
			}, by=c('IDX')]
	tmp
}
##--------------------------------------------------------------------------------------------------------
##	olli 30.04.16
##--------------------------------------------------------------------------------------------------------
treedist.closest.ind<- function(ph, model.reg)
{
	tmp			<- cophenetic.phylo(ph)
	diag(tmp)	<- Inf
	ans			<- data.table(IDPOP=rownames(tmp), IDPOP_CL=colnames(tmp)[apply(tmp, 1, which.min)])
	ans			<- merge(ans, ans[,  list(IDPOP_CL_GD=tmp[IDPOP, IDPOP_CL]), by='IDPOP'], by='IDPOP')
	if( !model.reg )
	{
		set(ans, NULL, 'IDPOP', ans[, toupper(gsub('-FEMALE|-MALE','',sapply(strsplit(IDPOP,'_',fixed=1),'[[',1)))] ) 
		set(ans, NULL, 'IDPOP_CL', ans[, toupper(gsub('-FEMALE|-MALE','',sapply(strsplit(IDPOP_CL,'_',fixed=1),'[[',1)))] )		
	}
	if( model.reg )
	{
		set(ans, NULL, 'IDPOP', ans[, sapply(strsplit(IDPOP,'|',fixed=1),'[[',1)] ) 
		set(ans, NULL, 'IDPOP_CL', ans[, sapply(strsplit(IDPOP_CL,'|',fixed=1),'[[',1)] )		
	}
	ans
}
##--------------------------------------------------------------------------------------------------------
##	olli 03.12.15
##--------------------------------------------------------------------------------------------------------
treedist.robinsonfould.wrapper<- function(submitted.info, ttrs, strs, check.binary.sim=TRUE)
{	
	setkey(submitted.info, IDX)
	#tmp				<- subset(submitted.info, IDX==321)[1,]
	#IDX<- 321;	TIME_IDX_T<-1
	tmp				<- submitted.info[, {
				cat('\nAt IDX', IDX)
				stree		<- unroot(strs[[IDX]])
				otree		<- unroot(multi2di(ttrs[[TIME_IDX_T]], random=FALSE))				
				if(check.binary.sim && !is.binary.tree(stree))
				{
					cat('\nFound non-binary tree at IDX',IDX)
					stree	<- multi2di(stree, random=FALSE)
				}
				#print(stree)
				#print(otree)
				z			<- setdiff(otree$tip.label, stree$tip.label)				
				stopifnot( length(z)==abs(diff(c(Ntip(otree), Ntip(stree)))) )
				if(length(z))
					otree	<- unroot(drop.tip(otree, z))				
				#https://groups.google.com/forum/#!topic/raxml/JgvxgknTeqw
				#normalize with 2n-6		
				rf			<- RF.dist(otree, stree, check.labels=TRUE)
				list(RF=rf, NRF=rf/(2*Ntip(otree)-6), TAXA_NJ=Ntip(otree))
			}, by='IDX']
	tmp
}
##--------------------------------------------------------------------------------------------------------
##	olli 03.12.15
##--------------------------------------------------------------------------------------------------------
treedist.robinsonfouldclusters.wrapper<- function(submitted.info, ttrs, strs, tinfo)
{
	#tmp		<- subset(submitted.info, MODEL=='Model: Regional')[1,]	
	#IDX<- 1; TIME_IDX_T<-12		
	setkey(tinfo, IDX_T)
	tmp		<- subset(submitted.info, MODEL=='R')[, {
				#IDX<- 208; TIME_IDX_T<- 16
				cat('\nAt IDX', IDX)
				stree		<- strs[[IDX]]
				otree		<- ttrs[[TIME_IDX_T]]				
				z			<- TIME_IDX_T
				z			<- subset(tinfo, CLU_N>3 & IDX_T==z)
				setkey(z, TAXA)
				z			<- unique(z)
				z			<- merge(z, data.table(TAXA=stree$tip.label, IN_STREE=1), by='TAXA', all.x=1)
				z			<- merge(z, z[, list(CLU_NS= length(which(IN_STREE==1))), by='IDCLU'], by='IDCLU')
				z			<- subset(z, CLU_NS>3)
				if(nrow(z))
				{
					#IDCLU	<- 6
					#TAXA	<- subset(z, IDCLU==6)[, TAXA]
					ans		<- z[, {								
								sclu	<- unroot(drop.tip(stree, setdiff(stree$tip.label,TAXA)))
								oclu	<- unroot(drop.tip(otree, union( setdiff(otree$tip.label, stree$tip.label), setdiff(otree$tip.label,TAXA))))
								rf		<- RF.dist(oclu, sclu, check.labels=TRUE)
								list(TAXA_NC=Ntip(oclu), RFC=rf, NRFC=rf/(2*Ntip(oclu)-6))
							}, by='IDCLU']	
				}
				if(!nrow(z))
					ans		<- data.table(IDCLU=NA_integer_, TAXA_NC=NA_integer_, RFC=NA_integer_, NRFC=NA_real_)
				ans			
			}, by='IDX']
	tmp
}
##--------------------------------------------------------------------------------------------------------
##	olli 07.11.15
##--------------------------------------------------------------------------------------------------------
treedist.billera.add<- function(submitted.info=NULL, ttrs=NULL, strs=NULL, file=NULL, with.save=0)
{
	#	file<- '/work/or105/Gates_2014/tree_comparison/submitted_151101.rda'
	require(ape)
	require(data.table)
	require(distory)
	stopifnot( !is.null(submitted.info) || !is.null(file))	
	if(is.null(submitted.info))
	{
		load(file)
		with.save	<- 1	
	}		
	tmp				<- submitted.info[, {
									tmp	<- lapply(IDX, function(i) strs[[i]]$tip.label)
									list(POSTHOC_ROOT=Reduce(intersect, tmp)[1])				
								}, by='IDX_T']
	submitted.info	<- merge(submitted.info, tmp, by='IDX_T')
	#tmp			<- subset(submitted.info, IDX==65)[1,]
	#IDX<- 65;	IDX_T<-3; POSTHOC_ROOT<-'IDPOP_101537|F|DOB_2000.58|2019.48'
	#POSTHOC_ROOT<-'HOUSE3326-7343-FEMALE_SAMPLED_30.4797372334259'
	tmp				<- submitted.info[, {
				cat('\nFT: At IDX', IDX)
				stree		<- strs[[IDX]]
				otree		<- multi2di(ttrs[[IDX_T]], random=FALSE)				
				if(!is.binary.tree(stree))
				{
					cat('\nFound non-binary tree at IDX',IDX)
					stree	<- multi2di(stree, random=FALSE)
				}					
				#print(stree)
				#print(otree)
				z			<- setdiff(otree$tip.label, stree$tip.label)
				stopifnot( length(z)==abs(diff(c(Ntip(otree), Ntip(stree)))) )				
				if(length(z))
					otree	<- drop.tip(otree, z)				
				otree		<- root(otree, outgroup=POSTHOC_ROOT, resolve.root=TRUE)
				stree		<- root(stree, outgroup=POSTHOC_ROOT, resolve.root=TRUE)
				tmp			<- data.table(TAXA=otree$tip.label, TAXA_NEW=seq_len(Ntip(otree)))
				otree$tip.label	<- tmp[, TAXA_NEW]
				setkey(tmp, TAXA)				
				stree$tip.label	<- tmp[stree$tip.label, ][, TAXA_NEW]				
				tmp			<- dist.multiPhylo( list(otree, stree) )[1]				
				list(BILL=tmp)
			}, by='IDX']
	submitted.info	<- merge(submitted.info, tmp, by='IDX')
	#
	if(with.save)
		save(strs, ttrs, tinfo, submitted.info, sclu.info, file=gsub('\\.rda','_BL\\.rda',file))
	#
	setkey(tinfo, IDX_T)
	#	IDX_T<- IDX<- 1
	tmp		<- subset(submitted.info, MODEL=='R')[, {
				cat('\nCT: At IDX', IDX)
				stree		<- strs[[IDX]]
				otree		<- ttrs[[IDX_T]]
				z			<- IDX_T
				z			<- subset(tinfo, CLU_N>3 & IDX_T==z)
				z			<- merge(z, data.table(TAXA=stree$tip.label, IN_STREE=1), by='TAXA', all.x=1)
				z			<- merge(z, z[, list(CLU_NS= length(which(IN_STREE==1))), by='IDCLU'], by='IDCLU')
				z			<- subset(z, CLU_NS>3)
				if(nrow(z))
				{
					#TAXA	<- subset(z, IDCLU==22)[, TAXA]
					ans		<- z[, {	
								print(IDCLU)
								sclu			<- drop.tip(stree, setdiff(stree$tip.label,TAXA), rooted=TRUE)
								oclu			<- drop.tip(otree, union( setdiff(otree$tip.label, stree$tip.label), setdiff(otree$tip.label,TAXA)), rooted=TRUE)								
								tmp				<- data.table(TAXA=oclu$tip.label, TAXA_NEW=seq_len(Ntip(oclu)))
								oclu$tip.label	<- tmp[, TAXA_NEW]
								setkey(tmp, TAXA)				
								sclu$tip.label	<- tmp[sclu$tip.label, ][, TAXA_NEW]
								if(!is.rooted(sclu) | !is.rooted(oclu))
								{
									sclu		<- root(sclu, outgroup='1',resolve.root=1)
									oclu		<- root(oclu, outgroup='1',resolve.root=1)
								}
								tmp				<- dist.multiPhylo( list(oclu, sclu) )[1]
								#print(tmp)								
								list(BILL=tmp)								
							}, by='IDCLU']	
				}
				if(!nrow(z))
					ans		<- data.table(IDCLU=NA_integer_, BILL=NA_real_)
				ans			
			}, by='IDX']	
	sclu.info	<- merge(submitted.info, tmp, by='IDX')		
	
	if(with.save)
		save(strs, ttrs, tinfo, submitted.info, sclu.info, file=gsub('\\.rda','_BL\\.rda',file))
}
##--------------------------------------------------------------------------------------------------------
##	olli 19.11.15
##--------------------------------------------------------------------------------------------------------
treecomparison.submissions.151119<- function()	
{
	require(data.table)
	require(ape)
	require(phangorn)
	#
	#	get true trees
	#
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim_internal/freeze_July15'
	tfiles	<- list.files(indir, pattern='newick$', full.names=TRUE)
	tfiles	<- data.table( FILE_T= tfiles[ grepl('SUBSTTREE', tfiles) | grepl('Vill_99', tfiles) | grepl('Vill.*DATEDTREE', tfiles) ] )
	tfiles[, SC:= toupper(gsub('_SUBSTTREE|_DATEDTREE','',gsub('.newick','',basename(FILE_T))))]
	tmp		<- rbind( subset(tfiles, SC=='VILL_99_APR15'), subset(tfiles, SC=='VILL_99_APR15'), subset(tfiles, SC=='VILL_99_APR15') )
	set(tmp, NULL, 'SC', c('150701_VILL_SCENARIO-C','150701_VILL_SCENARIO-D','150701_VILL_SCENARIO-E'))
	tfiles	<- rbind(tfiles, tmp)
	tmp		<- list.files(indir, pattern='newick$', full.names=TRUE)
	tmp		<- data.table( FILE_T= tmp[ grepl('Reg.*DATEDTREE', tmp) ] )
	tmp[, SC:= toupper(gsub('_SUBSTTREE|_DATEDTREE','',gsub('.newick','',basename(FILE_T))))]
	tfiles	<- rbind(tfiles, tmp)
	tfiles[, BRL_T:= 'time']	
	set(tfiles, tfiles[, which(grepl('REG',SC) & grepl('SUBST',FILE_T))], 'BRL_T', 'subst')	
	ttrs	<- lapply(tfiles[, FILE_T], function(x)	read.tree(file=x) )
	names(ttrs)	<- tfiles[, SC]	
	for(z in c('VILL_99_APR15','150701_VILL_SCENARIO-C','150701_VILL_SCENARIO-D','150701_VILL_SCENARIO-E'))
		ttrs[[z]]	<- root(ttrs[[z]], node=Ntip(ttrs[[z]])+2, resolve.root=1)	
	tfiles[, IDX_T:=seq_along(ttrs)]
	tfiles[, TAXAN_T:= sapply(ttrs, Ntip)]
	#	info on true trees
	tinfo	<- merge(tfiles, do.call('rbind',lapply(seq_along(ttrs), function(i) data.table(TAXA=ttrs[[i]]$tip.label, IDX_T=i))), by='IDX_T')	
	tinfo[, IDPOP:=NA_character_]
	tmp		<- tinfo[, which(grepl('REGIONAL',SC))]
	set(tinfo, tmp, 'IDPOP', tinfo[tmp,regmatches(TAXA, regexpr('IDPOP_[0-9]+',TAXA))])
	tmp		<- tinfo[, which(!grepl('REGIONAL',SC))]
	set(tinfo, tmp, 'IDPOP', tinfo[tmp, regmatches(TAXA, regexpr('HOUSE[0-9]+-[0-9]+|House[0-9]+-[0-9]+',TAXA))])		
	stopifnot(subset(tinfo, grepl('VILL',SC))[, length(which(substring(TAXA,1,10)!=substring(IDPOP,1,10)))]==0)	
	stopifnot( tinfo[, length(which(is.na(IDPOP)))==0] )	
	set(tinfo, NULL, 'IDPOP', tinfo[,toupper(IDPOP)])
	set(tinfo, NULL, 'TAXA', tinfo[,toupper(TAXA)])
	#	read cluster membership from DATEDCLUTREES	
	tmp		<- list.files(indir, pattern='DATEDCLUTREES', full.names=TRUE)
	tmp		<- data.table( 	FILE_CLU_T= tmp, 
							SC= toupper(gsub('_DATEDCLUTREES','',gsub('.newick','',basename(tmp)))),
							BRL_T= 'time') 
	tfiles	<- merge(tfiles, tmp, by=c('SC','BRL_T'), all=1)	
	tmp		<- subset(tfiles, !is.na(FILE_CLU_T))[, {
				z		<- read.tree(FILE_CLU_T)
				do.call('rbind',lapply(seq_along(z), function(i) data.table(IDCLU=i, TAXA=z[[i]]$tip.label)))				
			}, by=c('SC','BRL_T')]	
	tinfo	<- merge(tinfo, tmp, by=c('SC','BRL_T','TAXA'), all=1)
	tmp		<- subset(tinfo, !is.na(IDCLU))[, list(CLU_N= length(IDPOP)), by=c('SC','BRL_T','IDCLU')]
	tinfo	<- merge(tinfo, tmp, by=c('SC','BRL_T','IDCLU'), all=1)
	#	read sequences and determine %gappiness
	tmp		<- list.files(indir, pattern='fa$|fasta$', full.names=TRUE)
	tmp		<- data.table( FILE_SEQ_T= tmp, SC= toupper(gsub('_SIMULATED','',gsub('.fa','',basename(tmp)))))
	z		<- subset(tmp, SC=='VILL_99_APR15')
	set(z, NULL, 'SC', '150701_VILL_SCENARIO-C')
	tmp		<- rbind( tmp, z )	
	tfiles	<- merge(tfiles, tmp, by='SC', all=1)
	tmp		<- subset(tfiles, !is.na(FILE_SEQ_T))[, {
				z		<- read.dna(FILE_SEQ_T, format='fasta')	
				ans		<- sapply(seq_len(nrow(z)), function(i) base.freq(z[i,], all=1))
				ans		<- apply(ans[c('n','-','?'),], 2, sum)
				list(TAXA=rownames(z), GPS=ans)				
			}, by=c('SC','BRL_T')]
	tinfo	<- merge(tinfo, tmp, by=c('SC','BRL_T','TAXA'), all.x=1)
	#
	#	get submitted trees
	#	
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/IQTree/IQTree201507'
	infiles	<- list.files(indir, pattern='treefile$', recursive=1, full.names=1)
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/IQTree/IQTree201510'
	infiles	<- c(infiles, list.files(indir, pattern='treefile$', recursive=1, full.names=1))	
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/PhyML'
	infiles	<- c(infiles, list.files(indir, pattern='*tree*', recursive=1, full.names=1))
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/RAxML'
	infiles	<- c(infiles, list.files(indir, pattern='*RAxML_bestTree*', recursive=1, full.names=1))
	infiles	<- c(infiles, list.files(indir, pattern="best_tree.newick", recursive=1, full.names=1))
	infiles	<- data.table(FILE=infiles)
	strs	<- lapply(infiles[, FILE], function(x)
			{
				cat(x)
				read.tree(file=x)	
			})
	names(strs)	<- infiles[, FILE]
	
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/MetaPIGA'
	tmp		<-  list.files(indir, pattern='*result*', recursive=1, full.names=1)
	tmp		<- data.table(FILE=tmp)	
	tmp.trees			<- lapply(tmp[, FILE], function(x)
			{
				cat(x)
				read.nexus(file=x)	
			})
	sapply(tmp.trees, length)
	MetaPIGA.trees			<- c(lapply(tmp.trees, '[[', 1), lapply(tmp.trees, '[[', 2), lapply(tmp.trees, '[[', 3), lapply(tmp.trees, '[[', 4))
	names(MetaPIGA.trees)	<- c(sapply(tmp.trees, function(x) paste(names(x)[1],'_use',sep='')), sapply(tmp.trees, function(x) names(x)[2]), sapply(tmp.trees, function(x) names(x)[3]), sapply(tmp.trees, function(x) names(x)[4]))	
	names(MetaPIGA.trees)	<- gsub("'",'',names(MetaPIGA.trees), fixed=1)	
	strs					<- c(strs, MetaPIGA.trees)	
	submitted.info			<- data.table(FILE=names(strs))
	#
	#
	#	
	submitted.info[, IDX:=seq_along(strs)]	
	submitted.info[, TEAM:=NA_character_]
	set(submitted.info, submitted.info[, which(grepl('RAXML|RAxML',FILE))], 'TEAM', 'RAXML')
	set(submitted.info, submitted.info[, which(grepl('IQTree',FILE))], 'TEAM', 'IQTree')
	set(submitted.info, submitted.info[, which(grepl('MetaPIGA|Consensus pruning|Best individual of population',FILE))], 'TEAM', 'MetaPIGA')
	set(submitted.info, submitted.info[, which(grepl('PhyML',FILE))], 'TEAM', 'PhyML')	
	stopifnot( submitted.info[, length(which(is.na(TEAM)))==0] )
	#
	#	scenario
	#	
	submitted.info[, SC:=NA_character_]
	tmp		<- submitted.info[, which(grepl('150701_Regional_TRAIN[0-9]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, regmatches(FILE, regexpr('150701_Regional_TRAIN[0-9]',FILE))])
	tmp		<- submitted.info[, which(grepl('150701_Vill_SCENARIO-[A-Z]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, regmatches(FILE, regexpr('150701_Vill_SCENARIO-[A-Z]',FILE))])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('TRAIN[0-9]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, paste('150701_Regional_',regmatches(FILE, regexpr('TRAIN[0-9]',FILE)),sep='')])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('scenario[A-Z]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, paste('150701_Vill_',regmatches(FILE, regexpr('scenario[A-Z]',FILE)),sep='')])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('150701_regional_train[0-9]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, regmatches(FILE, regexpr('150701_regional_train[0-9]',FILE))])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('150701_vill_scenario-[A-Z]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, regmatches(FILE, regexpr('150701_vill_scenario-[A-Z]',FILE))])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('Vill_99_Apr15', FILE))]
	set(submitted.info, tmp, 'SC', 'Vill_99_Apr15')	
	set(submitted.info, NULL, 'SC', submitted.info[, toupper(SC)])
	tmp		<- submitted.info[, which(grepl('150701_VILL_SCENARIO[A-Z]', SC))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, gsub('150701_VILL_SCENARIO','150701_VILL_SCENARIO-',SC)])
	stopifnot( submitted.info[, length(which(is.na(SC)))==0] )
	#
	#	set covariates of scenarios
	#
	tmp		<- data.table(	SC=		c("150701_REGIONAL_TRAIN1","150701_REGIONAL_TRAIN2","150701_REGIONAL_TRAIN3","150701_REGIONAL_TRAIN4" ,"150701_REGIONAL_TRAIN5", "150701_VILL_SCENARIO-A", "150701_VILL_SCENARIO-B", "VILL_99_APR15","150701_VILL_SCENARIO-C", "150701_VILL_SCENARIO-D", "150701_VILL_SCENARIO-E"),
			MODEL=	c('R','R','R','R','R','V','V','V','V','V','V'),
			SEQCOV= c(0.16, 0.16, 0.16, 0.16, 0.16, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6),
			ACUTE=	c('low', 'low', 'high', 'low', 'high', 'high', 'high', 'high', 'high', 'high', 'high'),
			GAPS=	c('none', 'low', 'low', 'high', 'high', 'low', 'high', 'none', 'none', 'low', 'high'), 
			ART=	c('none', 'none', 'none', 'none', 'none', 'none', 'none', 'fast', 'fast', 'fast', 'fast'),
			EXT= 	c('5pc', '5pc', '5pc', '5pc', '5pc', '~0pc', '~0pc', '~0pc', '~0pc', '~0pc', '~0pc')
	)
	submitted.info	<- merge(submitted.info, tmp, by='SC')
	#
	#	best tree for each scenario
	#
	submitted.info[, BEST:='N']
	set(submitted.info, submitted.info[, which(grepl('RAxML', FILE) & grepl('best_tree', FILE))], 'BEST', 'Y')	
	#	copied from ListOfBestTrees_IQTree150818.txt
	#	there are several best trees for some scenarios
	tmp	<- c( '150701_Vill_SCENARIO-A_IQTree150814_partition_12_3_07',	
			'150701_Vill_SCENARIO-A_IQTree150814_partition_12_3_04.',	
			'150701_Vill_SCENARIO-B_IQTree150814_partition_12_3_03.',	
			'Vill_99_Apr15_IQTree150814_partition_123.',			
			'150701_Vill_SCENARIO-D_IQTree150814_partition_12_3.',	
			'150701_Vill_SCENARIO-E_IQTree150814_partition_12_3.',
			'150701_Vill_SCENARIO-A_IQTree150814_pol_partition_12_3.',	
			'150701_Vill_SCENARIO-B_IQTree150814_pol_partition_12_3_05.',
			'Vill_99_Apr15_IQTree150814_pol_partition_12_3_09.',		
			'Vill_99_Apr15_IQTree150814_pol_partition_12_3_10.',		
			'150701_Vill_SCENARIO-D_IQTree150814_pol_partition_12_3_05.',
			'150701_Vill_SCENARIO-D_IQTree150814_pol_partition_12_3_06.',
			'150701_Vill_SCENARIO-D_IQTree150814_pol_partition_12_3_09.',
			'150701_Vill_SCENARIO-E_IQTree150814_pol_partition_12_3_06.',
			'150701_Regional_TRAIN1_IQTree150818_partition_123_03.',	
			'150701_Regional_TRAIN1_IQTree150818_pol_partition_123_05.')
	tmp	<- sapply(tmp, function(x) submitted.info[, which((grepl('IQTree150814/', FILE, fixed=1) | grepl('IQTree150818/', FILE, fixed=1)) & grepl(x, FILE, fixed=1))] )
	set(submitted.info, tmp, 'BEST', 'Y')
	tmp	<- c('150701_Regional_TRAIN2_IQTree151019_partition_123_10', 			
			'150701_Regional_TRAIN3_IQTree151019_partition_123_03',	
			'150701_Regional_TRAIN4_IQTree151019_partition_123_10',	
			'150701_Regional_TRAIN5_IQTree151019_partition_123_01',
			'150701_Regional_TRAIN2_IQTree151019_pol_partition_123_08',
			'150701_Regional_TRAIN3_IQTree151019_pol_partition_123_08',	
			'150701_Regional_TRAIN4_IQTree151019_pol_partition_123_05',
			'150701_Regional_TRAIN5_IQTree151019_pol_partition_123_10')
	tmp	<- sapply(tmp, function(x) submitted.info[, which((grepl('IQTree151019', FILE, fixed=1)) & grepl(x, FILE, fixed=1))] )
	set(submitted.info, tmp, 'BEST', 'Y')
	#	PhyML no replicates: all files best
	set(submitted.info, submitted.info[, which(TEAM=='PhyML')], 'BEST', 'Y')		
	#
	#	set OTHER (ie old or some preliminary/unknown tree)
	#
	submitted.info[, OTHER:='N']
	#	MetaPIGA tree to be used is first in nexus list (which was tagged with best above)
	set(submitted.info, submitted.info[, which(TEAM=='MetaPIGA' & !grepl('use', FILE))], 'OTHER', 'Y')
	#	IQTree did several uploads, use only most recent in main analysis
	set(submitted.info, submitted.info[, which(grepl('150701_Regional_TRAIN1_IQTree150814', FILE))], 'OTHER', 'Y')
	set(submitted.info, submitted.info[, which(TEAM=='IQTree' & MODEL=='R' & !grepl('TRAIN1', SC) & grepl('201507/',FILE,fixed=1))], 'OTHER', 'Y')
	#
	#	set which gene used to construct tree (either pol or concatenated gag+pol+env)
	#
	submitted.info[, GENE:=NA_character_]
	set(submitted.info, submitted.info[, which(TEAM=='RAXML' & grepl('full', FILE))], 'GENE', 'GAG+POL+ENV')
	set(submitted.info, submitted.info[, which(TEAM=='RAXML' & grepl('pol', FILE))], 'GENE', 'POL')
	stopifnot(nrow(subset(submitted.info, TEAM=='RAXML' & is.na(GENE)))==0)
	set(submitted.info, submitted.info[, which(TEAM=='PhyML')], 'GENE', 'POL')
	set(submitted.info, submitted.info[, which(TEAM=='MetaPIGA')], 'GENE', 'GAG+POL+ENV')	
	set(submitted.info, submitted.info[, which(TEAM=='IQTree' & grepl('[0-9]_partition', FILE))], 'GENE', 'GAG+POL+ENV')
	set(submitted.info, submitted.info[, which(TEAM=='IQTree' & grepl('[0-9]_pol_partition', FILE))], 'GENE', 'POL')
	stopifnot(nrow(subset(submitted.info, TEAM=='IQTree' & is.na(GENE)))==0)
	#
	#	number taxa in tree
	#
	setkey(submitted.info, IDX)
	submitted.info[, TAXAN:= sapply(strs, Ntip)]
	#
	#	are trees rooted?
	#
	setkey(submitted.info, IDX)
	submitted.info[, ROOTED:=factor(sapply(strs, is.rooted),levels=c(TRUE,FALSE),labels=c('Y','N'))]
	#
	#	add BRL_UNITS
	#
	submitted.info[, BRL:='subst']
	#
	#	add index of true tree
	#
	require(phangorn)
	tmp				<- subset(tfiles, select=c('SC','IDX_T','BRL_T'))
	tmp				<- dcast.data.table(tmp, SC~BRL_T, value.var='IDX_T')
	setnames(tmp, c('subst','time'), c("SUB_IDX_T","TIME_IDX_T"))
	submitted.info	<- merge(submitted.info, tmp, by='SC')
	submitted.info	<- merge(submitted.info, unique(subset(tfiles, select=c('SC','TAXAN_T'))), by='SC')	
	stopifnot(nrow(subset(submitted.info, TAXAN>TAXAN_T))==0)
	#
	#	fix taxa names that teams have changed
	#
	tmp		<- subset(submitted.info, TEAM=='IQTree' & MODEL=='R')[, IDX]
	for(i in tmp)
	{
		strs[[i]]$tip.label	<- sapply(strsplit(strs[[i]]$tip.label,'_'), function(x)	paste(x[1],'_',x[2],'|',x[3],'|',x[4],'_',x[5],'|',x[6],sep='')	)
	}
	for(i in seq_along(strs))
	{
		strs[[i]]$tip.label	<- toupper(strs[[i]]$tip.label)
	}
	for(i in seq_along(ttrs))
	{
		ttrs[[i]]$tip.label	<- toupper(ttrs[[i]]$tip.label)
	}	
	###
	tmp		<- subset(submitted.info, TEAM=='PhyML' & MODEL=='R')[, IDX]
	for(i in tmp)
	{
		
		z	<- data.table(IDX=seq_along(strs[[i]]$tip.label), IDPOP=regmatches(strs[[i]]$tip.label, regexpr('IDPOP_[0-9]+',strs[[i]]$tip.label)), SC=subset(submitted.info, IDX==i)[,SC])
		z	<- merge(subset(tinfo, BRL_T=='time', select=c(IDPOP,SC,TAXA)), z, by=c('IDPOP','SC'))
		setkey(z, IDX)
		stopifnot(nrow(z)==Ntip(strs[[i]]))
		strs[[i]]$tip.label	<- z[, TAXA]
	}
	tmp		<- subset(submitted.info, TEAM=='PhyML' & MODEL=='V')[, IDX]
	for(i in tmp)
	{
		
		z	<- data.table(IDX=seq_along(strs[[i]]$tip.label), IDPOP=regmatches(strs[[i]]$tip.label, regexpr('HOUSE[0-9]+-[0-9]+|House[0-9]+-[0-9]+',strs[[i]]$tip.label)), SC=subset(submitted.info, IDX==i)[,SC])
		z	<- merge(subset(tinfo, BRL_T=='time', select=c(IDPOP,SC,TAXA)), z, by=c('IDPOP','SC'))
		stopifnot(nrow(z)==length(strs[[i]]$tip.label))
		setkey(z, IDX)
		strs[[i]]$tip.label	<- z[, TAXA]
	}
	#
	#	compute Robinson Fould of complete tree
	#
	tmp				<- treedist.robinsonfould.wrapper(submitted.info, ttrs, strs)
	submitted.info	<- merge(submitted.info, tmp, by='IDX')
	#	compute Robinson Fould of clusters, then take sum
	tmp				<- treedist.robinsonfouldclusters.wrapper(submitted.info, ttrs, strs, tinfo)
	sclu.info		<- merge(subset(submitted.info, select=c("IDX","SC","FILE","TEAM","MODEL","SEQCOV","ACUTE","GAPS","ART","EXT","BEST","OTHER","GENE","TAXAN","ROOTED","BRL","SUB_IDX_T","TIME_IDX_T","TAXAN_T")), tmp, by='IDX')
	#
	strs.new			<- strs
	ttrs.new			<- ttrs
	tinfo.new			<- copy(tinfo)
	submitted.info.new	<- copy(submitted.info)
	sclu.info.new		<- copy(sclu.info)
	outdir		<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation'
	z			<- load(paste(outdir, 'submitted_151119_S.rda', sep='/'))	
	stopifnot( nrow(subset(merge(subset(submitted.info, select=c('FILE','IDX')), subset(submitted.info.new, select=c('FILE','IDX')), by='FILE'), IDX.x!=IDX.y))==0 )	
	submitted.info		<- merge(submitted.info.new, subset(submitted.info, select=c('IDX','NQD','lm_intercept','lm_slope','lm_rsq')), by='IDX')
	strs				<- strs.new
	ttrs				<- ttrs.new
	sclu.info			<- merge(sclu.info.new, subset(sclu.info, select=c('IDX','IDCLU','NQDC')), by=c('IDX','IDCLU'))
	#
	outfile	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation/submitted_151119_SRFQD.rda'
	save(strs, strs_lsd_brl, strs_lsd_date, ttrs, tinfo, submitted.info, sclu.info, file=outfile)
}
##--------------------------------------------------------------------------------------------------------
##	olli 25.07.16
##--------------------------------------------------------------------------------------------------------
treecomparison.explaingaps.bams.160817<- function()
{
	require(scales)
	require(ggplot2)
	#
	#	collect information on batches
	#
	wdir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/explaingaps'
	load(file.path(wdir, 'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.rda'))
	load(file.path(wdir, 'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_logistic.rda'))
	bami	<- subset(dpand, PR=='1R' & POS=='PR_1' & !is.na(PANGEA_ID) & !is.na(SANGER_ID), select=c(PANGEA_ID, STUDY_ID, SANGER_ID))
	#	extract location of samples in batch
	bami[, BATCH:= sapply(strsplit(SANGER_ID,'_'),'[[',1)]
	bami[, LOC:= substr(gsub('PG[0-9]+-','',PANGEA_ID),1,2)]
	bami	<- bami[, list(LOC= paste(unique(LOC),collapse='-')), by='BATCH']
	#	extract significant batches
	tmp		<- subset(m2f.1.or, u95<0.95)
	tmp[, PR:='2F']
	tmp2	<- subset(m2r.1.or, u95<0.95)
	tmp2[, PR:='2R']
	tmp		<- rbind(tmp, tmp2)
	tmp		<- subset(tmp, grepl('BATCH', COEF))
	set(tmp, NULL, 'BATCH', tmp[,gsub('BATCH','',COEF)])
	#	
	setkey(tmp, PR, OR)
	ggplot(tmp, aes(x=BATCH, y=OR, ymin=l95, ymax=u95)) + geom_point() + geom_errorbar() + facet_grid(~PR) + coord_flip()
	ggsave(file=file.path(wdir, 'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_batch_oddsratio.pdf'), w=7, h=10)
	#	batches are typically significant in both 2F and 2R	
	tmp2	<- dcast.data.table(tmp, BATCH~PR, value.var='OR')
	setnames(tmp2, c('2F','2R'), c('OR_2F','OR_2R'))
	bami	<- merge(bami, tmp2, by='BATCH',all.x=1)
	tmp2	<- dcast.data.table(tmp, BATCH~PR, value.var='u95')
	setnames(tmp2, c('2F','2R'), c('ORu95_2F','ORu95_2R'))
	bami	<- merge(bami, tmp2, by='BATCH',all.x=1)
	tmp2	<- dcast.data.table(tmp, BATCH~PR, value.var='l95')
	setnames(tmp2, c('2F','2R'), c('ORl95_2F','ORl95_2R'))
	bami	<- merge(bami, tmp2, by='BATCH',all.x=1)	
	bami	<- melt(bami, id.vars=c('BATCH','LOC'))	
	set(bami, bami[, which(is.na(value))],'value', 1)
	bami	<- dcast.data.table(bami, BATCH+LOC~variable, value.var='value')
	#
	#	get length of short reads
	#	
	infile	<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/readlengths/bam_stats_150218.rda'
	load(infile)	
	bam.len[, BATCH:= sapply(strsplit(FILE,'_'),'[[',1)]
	tmp		<- subset(bam.len,QU>=40 & QU<320)[, list(CDF.bamlen= mean(CDF)), by=c('BATCH','QU')]
	tmp		<- merge(tmp, bami, by='BATCH')
		
	ggplot(tmp, aes(y=CDFm, x=QU, colour=LOC, group=BATCH)) +
			geom_line() +
			#geom_boxplot() + 
			scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
			theme_bw() + labs(y='average frequency in run\n(cumulated)\n', x='\nlength of quality-trimmed short reads\n(nt)', fill='sequence run') +
			theme(legend.position='bottom') 
	ggsave(file=file.path(wdir, 'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_bamlen_by_run.pdf'), w=14, h=7)
	
	ggplot(tmp, aes(y=CDF.bamlen, x=QU, colour=OR_2F, group=BATCH)) +
			geom_line() +
			scale_x_continuous(breaks=seq(0,5e2,50), expand=c(0,0)) +
			scale_colour_gradientn(colours=c('red','grey50')) +
			scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
			theme_bw() + labs(y='average frequency of read length of individuals in batch\n(cumulated)\n', x='\nlength of quality-trimmed short reads\n(nt)', colour='odds ratio for batch effect\nin 2F primer analysis\n(0 is strong batch effect)      ') +
			theme(legend.position='bottom') 
	ggsave(file=file.path(wdir, 'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_bamlen_by_OR2F.pdf'), w=14, h=7)
	
	ggplot(tmp, aes(y=CDF.bamlen, x=QU, colour=OR_2R, group=BATCH)) +
			geom_line() +
			scale_x_continuous(breaks=seq(0,5e2,50), expand=c(0,0)) + 
			scale_colour_gradientn(colours=c('red','grey50')) +
			scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
			theme_bw() + 
			labs(y='average frequency of read length of individuals in batch\n(cumulated)\n', x='\nlength of quality-trimmed short reads\n(nt)', colour='odds ratio for batch effect\nin 2R primer analysis\n(0 is strong batch effect)      ') +
			theme(legend.position='bottom') 
	ggsave(file=file.path(wdir, 'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_bamlen_by_OR2R.pdf'), w=14, h=7)
	#
	#	OK so those significant on the 2R primer are also a bit shorter than others
	#
	#
	bam.cov[, BATCH:= sapply(strsplit(FILE_ID,'_'),'[[',1)]
	bam.cov	<- merge(bam.cov, bami, by='BATCH')
	setkey(bam.cov, BATCH, FILE_ID, POS)
	bam.covs<- lapply(bam.cov[, unique(BATCH)], function(batch)
			{
				cat('\nprocess run', batch)
				tmp		<- subset(bam.cov, BATCH==batch)
				tmp		<- tmp[, 	{
							z	<- rep(COV,REP)
							list(COV=z, POS=seq_along(z), REF=REF[1])
						}, by=c('FILE_ID','BATCH')]
				tmp		<- tmp[, list(QU=paste('QU',100*c(0.25, 0.5, 0.75),sep=''), COV=quantile(COV, p=c(0.25, 0.5, 0.75))), by=c('BATCH','REF','POS')]
				tmp		<- dcast.data.table(tmp, REF+BATCH+POS~QU, value.var='COV')
				tmp
			})
	bam.covs<- do.call('rbind', bam.covs)
	stopifnot( bam.covs[, length(unique(REF))]==1 )	#if not the same reference, then the positions do not share the same coordinate system
	bam.covs<- merge(bam.covs, bami, by='BATCH')
	
	ggplot(bam.covs, aes(y=QU50, x=POS, colour=OR_2R, group=BATCH)) +
			geom_line() +			 
			scale_colour_gradientn(colours=c('red','grey50')) +
			scale_x_continuous(expand=c(0,0), breaks=seq(0,2e4,1e3)) +
			scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
			theme_bw() + labs(y='average number of reads per individual in batch\n', x='\nnt position in reference', colour='odds ratio for batch effect\nin 2R primer analysis\n(0 is strong batch effect)      ') +
			theme(legend.position='bottom') 
	ggsave(file=file.path(wdir, 'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_bamcov_by_OR2R.pdf'), w=14, h=7)
	ggplot(bam.covs, aes(y=QU50, x=POS, colour=OR_2F, group=BATCH)) +
			geom_line() +			 
			scale_colour_gradientn(colours=c('red','grey50')) +
			scale_x_continuous(expand=c(0,0), breaks=seq(0,2e4,1e3)) +
			scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
			theme_bw() + labs(y='average number of reads per individual in batch\n', x='\nnt position in reference', colour='odds ratio for batch effect\nin 2F primer analysis\n(0 is strong batch effect)      ') +
			theme(legend.position='bottom') 
	ggsave(file=file.path(wdir, 'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_bamcov_by_OR2F.pdf'), w=14, h=7)
	
	
	tmp		<- subset(bam.covs, POS>=1500 & POS<=4500)
	tmp[, OR_2F:= factor(OR_2F<1, levels=c(TRUE,FALSE), labels=c('Yes','No'))]
	tmp[, OR_2R:= factor(OR_2R<1, levels=c(TRUE,FALSE), labels=c('Yes','No'))]
	tmp[, POSc:= as.numeric(as.character(cut(POS, breaks=seq(1500,4500,250), labels=seq(1501,4500,250))))]
	tmp		<- subset(tmp, !is.na(POSc))
	set(tmp, NULL, 'POSc', tmp[, paste(POSc,'-',POSc+249,sep='')])
	ggplot(tmp, aes(x=POSc, y=QU50, fill=OR_2F)) + geom_boxplot(outlier.shape=NA) +
			scale_fill_manual(values=c('Yes'='red','No'='grey50')) +
			scale_x_discrete(expand=c(0,0)) +
			coord_trans(y='log10p1') +
			scale_y_continuous(breaks=c(1,10,100, 1e3,1e4, 1e5)) +
			theme_bw() + 
			labs(y='average number of reads per individual in batch\n', x='\nnt position in reference', fill='significant odds ratio for batch effect\nin 2F primer analysis      ') +
			theme(legend.position='bottom')
	ggsave(file=file.path(wdir, 'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_bamcov_by_OR2F_boxpol.pdf'), w=12, h=6)
	ggplot(tmp, aes(x=POSc, y=QU50, fill=OR_2R)) + geom_boxplot(outlier.shape=NA) +
			scale_fill_manual(values=c('Yes'='red','No'='grey50')) +
			scale_x_discrete(expand=c(0,0)) +
			coord_trans(y='log10p1') +
			scale_y_continuous(breaks=c(1,10,100, 1e3,1e4, 1e5)) +
			theme_bw() + 
			labs(y='average number of reads per individual in batch\n', x='\nnt position in reference', fill='significant odds ratio for batch effect\nin 2R primer analysis     ') +
			theme(legend.position='bottom')
	ggsave(file=file.path(wdir, 'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_bamcov_by_OR2R_boxpol.pdf'), w=12, h=6)
	
	tmp		<- copy(bam.covs)
	tmp[, PLOT_ROW:= factor(POS<=4500, levels=c(TRUE,FALSE),labels=c('1st half of genome','2nd half of genome'))]
	tmp[, OR_2F:= factor(OR_2F<1, levels=c(TRUE,FALSE), labels=c('Yes','No'))]
	tmp[, OR_2R:= factor(OR_2R<1, levels=c(TRUE,FALSE), labels=c('Yes','No'))]
	tmp[, POSc:= as.numeric(as.character(cut(POS, breaks=seq(0,9250,250), labels=seq(1,9250,250))))]
	tmp		<- subset(tmp, !is.na(POSc))
	setkey(tmp, POSc)
	tmp2	<- tmp[, unique(POSc)]
	set(tmp, NULL, 'POSc', tmp[, factor(POSc, levels=tmp2, labels=paste(tmp2,'-',tmp2+249,sep=''))])
	ggplot(tmp, aes(x=POSc, y=log10(QU50+.1), fill=OR_2R)) + geom_boxplot(outlier.shape=NA) +
			scale_fill_manual(values=c('Yes'='red','No'='grey50')) +
			scale_x_discrete(expand=c(0,0)) +			
			scale_y_continuous(breaks=c(0,1,2,3,4,5), labels=c(1,10,100, 1e3,1e4, 1e5)) +
			theme_bw() + 
			facet_wrap(~PLOT_ROW, ncol=1, scale='free_x') +
			labs(y='average number of reads per individual in batch\n', x='\nnt position in reference', fill='significant odds ratio for batch effect   ') +
			theme(legend.position='bottom')
	ggsave(file=file.path(wdir, 'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_bamcov_by_OR2R_boxall.pdf'), w=16, h=12)
	
}
##--------------------------------------------------------------------------------------------------------
##	olli 25.07.16
##--------------------------------------------------------------------------------------------------------
treecomparison.samplecharacteristics.160817<- function()
{
	require(ape)
	require(scales)
	require(ggplot2)
	require(data.table)
	require(Hmisc)
	
	wdir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/explaingaps'
	wfile	<- 'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.rda'
	load(file.path(wdir,wfile))
	#
	#	table on sequence characteristics by cohort
	#
	pngi	<- dpand[, list(NT_DIFF_SUM=sum(NT_DIFF)), by=c('TAXA','PR')]
	set(pngi, NULL, 'PR', pngi[, paste('NT_DIFF_',PR,sep='')])	
	set(pngi, pngi[, which(NT_DIFF_SUM>=1)], 'NT_DIFF_SUM', 1L)
	set(pngi, pngi[, which(is.na(NT_DIFF_SUM))], 'NT_DIFF_SUM', 2L)
	set(pngi, NULL, 'PR_NTDIFFc', pngi[, as.character(factor(NT_DIFF_SUM, levels=c(2,0,1), labels=c('unassembled','no mutation','at least one mutation')))])	
	pngi	<- dcast.data.table(pngi, TAXA~PR, value.var='PR_NTDIFFc')	
	pngi	<- merge(pngi, dm, by='TAXA')
	#	sex by cohort
	set(pngi, pngi[, which(is.na(SEX))], 'SEX', 'missing')
	pngs	<- pngi[, {
							z<- table(SEX)
							list(STAT='SEX', N= as.numeric(z), LABEL= attr(z, 'dimnames')[[1]])
						}, by='COHORT']
	pngs	<- merge(pngs, as.data.table(expand.grid(COHORT=pngs[, unique(COHORT)], STAT='SEX', LABEL= pngs[, unique(LABEL)], stringsAsFactors=FALSE)), by=c('COHORT','STAT','LABEL'), all=1)
	set(pngs, NULL, 'LABEL', pngs[, factor(LABEL, levels=c("F", "M", "missing"))])
	#	if only sampling year known, set to midpoint
	tmp		<- pngi[, which(SAMPLEDATE==floor(SAMPLEDATE))]
	set(pngi, tmp, 'SAMPLEDATE', pngi[tmp, SAMPLEDATE+.5])
	#	sampling years
	pngi[, SAMPLEYR:= as.character(floor(SAMPLEDATE))]
	set(pngi, pngi[, which(is.na(SAMPLEYR))], 'SAMPLEYR', 'missing')
	tmp		<- pngi[, {
							z<- table(SAMPLEYR)
							list(STAT='SAMPLEYR', N= as.numeric(z), LABEL= attr(z, 'dimnames')[[1]])
						}, by='COHORT']
	tmp		<- merge(tmp, as.data.table(expand.grid(COHORT=tmp[, unique(COHORT)], STAT='SAMPLEYR', LABEL= tmp[, unique(LABEL)], stringsAsFactors=FALSE)), by=c('COHORT','STAT','LABEL'), all=1)
	set(tmp, NULL, 'LABEL', tmp[, factor(LABEL, levels=c("2010","2011", "2012", "2013","2014","missing"))])
	pngs	<- rbind(pngs, tmp)		
	#	recent viral load around sampling time
	pngi[, DUMMY:= 'missing']
	tmp		<- pngi[, which(!is.na(RECENTVLDATE) & !is.na(SAMPLEDATE) & abs(RECENTVLDATE-SAMPLEDATE)<1)]
	set(pngi, tmp, 'DUMMY', pngi[tmp, cut(as.numeric(RECENTVL), right=FALSE, breaks=c(-Inf, 1e4, 2e4, 4e4, 1e5, Inf), labels=c('<10,000','10,000-19,999','20,000-39,999','40,000-99,999','>=100,000'))])
	tmp		<- pngi[, {
				z<- table(DUMMY)
				list(STAT='RECENTVL', N= as.numeric(z), LABEL= attr(z, 'dimnames')[[1]])
			}, by='COHORT']
	tmp		<- merge(tmp, as.data.table(expand.grid(COHORT=tmp[, unique(COHORT)], STAT='RECENTVL', LABEL= tmp[, unique(LABEL)], stringsAsFactors=FALSE)), by=c('COHORT','STAT','LABEL'), all=1)
	set(tmp, NULL, 'LABEL', tmp[, factor(LABEL, levels=c('<10,000','10,000-19,999','20,000-39,999','40,000-99,999','>=100,000',"missing"))])
	pngs	<- rbind(pngs, tmp)		
	#	ART
	pngi[, DUMMY:= NULL]
	pngi[, DUMMY:= 'N']
	tmp		<- pngi[, which(ARTSTART==floor(ARTSTART))]
	set(pngi, tmp, 'ARTSTART', pngi[tmp, ARTSTART+.5])
	tmp		<- pngi[, which(ARTSTART<=SAMPLEDATE | PREVARTUSE=='Y' | COHORT=='BW-Mochudi' & CURRENTLYONART=='Y' | (everSelfReportArt==1 & FirstSelfReportArt<SAMPLEDATE))]
	set(pngi, tmp, 'DUMMY', 'Y')
	tmp		<- pngi[, which( (COHORT=='RCCS' & (is.na(everSelfReportArt))) | 	#RCCS does not code CURRENTLYONART, and missing ARTSTART may indicate ART not yet started   
							 (COHORT=='BW-Mochudi' & is.na(CURRENTLYONART)) |						#BW codes already on ART 
							 (COHORT=='AC_Resistance' & is.na(ARTSTART)) |							#in the resistance cohort, there must be an ART start date
							 (COHORT=='UG-MRC' & is.na(CURRENTLYONART)))]
	set(pngi, tmp, 'DUMMY', 'missing')
	tmp		<- pngi[, {
				z<- table(DUMMY)
				list(STAT='PREVARTATSAMPLING', N= as.numeric(z), LABEL= attr(z, 'dimnames')[[1]])
			}, by='COHORT']
	tmp		<- merge(tmp, as.data.table(expand.grid(COHORT=tmp[, unique(COHORT)], STAT='PREVARTATSAMPLING', LABEL= tmp[, unique(LABEL)], stringsAsFactors=FALSE)), by=c('COHORT','STAT','LABEL'), all=1)
	set(tmp, NULL, 'LABEL', tmp[, factor(LABEL, levels=c('Y','N',"missing"))])
	pngs	<- rbind(pngs, tmp)	
	#	Age
	tmp		<- pngi[, which(is.na(AGE) & !is.na(DOB) & !is.na(SAMPLEDATE))]
	set(pngi, tmp, 'AGE', pngi[tmp, SAMPLEDATE-DOB])
	pngi[, DUMMY:= NULL]
	pngi[, DUMMY:= 'missing']
	tmp		<- pngi[, which(!is.na(AGE))]
	set(pngi, tmp, 'DUMMY', pngi[tmp, cut(AGE, breaks=c(-Inf, 25, 30, 35, 40, Inf), labels=c('<25','<30','<35','<40','>=40'))])
	tmp		<- pngi[, {
				z<- table(DUMMY)
				list(STAT='AGE', N= as.numeric(z), LABEL= attr(z, 'dimnames')[[1]])
			}, by='COHORT']	
	tmp		<- merge(tmp, as.data.table(expand.grid(COHORT=tmp[, unique(COHORT)], STAT='AGE', LABEL= tmp[, unique(LABEL)], stringsAsFactors=FALSE)), by=c('COHORT','STAT','LABEL'), all=1)
	set(tmp, NULL, 'LABEL', tmp[, factor(LABEL, levels=c('<25','<30','<35','<40','>=40',"missing"))])
	pngs	<- rbind(pngs, tmp)	
	#	Subtype
	tmp		<- pngi[, {
			z<- table(COMET_CONS)
			list(STAT='SUBTYPE', N= as.numeric(z), LABEL= attr(z, 'dimnames')[[1]])
		}, by='COHORT']	
	tmp		<- merge(tmp, as.data.table(expand.grid(COHORT=tmp[, unique(COHORT)], STAT='SUBTYPE', LABEL= tmp[, unique(LABEL)], stringsAsFactors=FALSE)), by=c('COHORT','STAT','LABEL'), all=1)
	pngs	<- rbind(pngs, tmp)
	#	Primer assembled
	tmp		<- melt(subset(pngi, select=c('COHORT',colnames(pngi)[grepl('NT_DIFF',colnames(pngi))])), id.vars='COHORT')
	tmp		<- tmp[, {
				z<- table(value)
				list(STAT='PRIMER', N= as.numeric(z), LABEL= attr(z, 'dimnames')[[1]])	
			}, by=c('COHORT','variable')]
	set(tmp, NULL, 'STAT', tmp[, paste(STAT,gsub('NT_DIFF','',variable),sep='')])
	tmp[, variable:=NULL]	
	pngs	<- rbind(pngs, tmp)
	#
	#	add proportions
	#
	set(pngs, pngs[, which(is.na(N))], 'N', 0)
	tmp		<- subset(pngs, STAT=='SEX')[, list(STAT='TOTAL', LABEL='TOTAL', N=sum(N)), by='COHORT']	
	pngs	<- merge(pngs, dcast.data.table(tmp, COHORT~STAT, value.var='N'), by='COHORT')
	pngs[, P:= 100*round(N/TOTAL, d=2)]
	#
	#	make table
	#
	pngst	<- dcast.data.table(pngs, STAT+LABEL~COHORT, value.var='P')
	tmp		<- subset(pngs, STAT=='SEX')[, list(STAT='TOTAL', LABEL='TOTAL', N=sum(N)), by='COHORT']
	tmp		<- dcast.data.table(tmp, STAT~COHORT, value.var='N')
	tmp[, LABEL:='']
	pngst	<- rbind(subset(tmp, select=c('STAT', 'LABEL', 'AC_Resistance', 'BW-Mochudi', 'RCCS', 'UG-MRC')), pngst,use.names=TRUE,fill=TRUE)
	#
	save(pngi, pngs, pngst, file=file.path(wdir, gsub('\\.rda','_samplecharacteristics.rda', wfile)))
	write.csv(pngst, row.names=FALSE, file=file.path(wdir, gsub('\\.rda','_samplecharacteristics.csv', wfile)))
	#
	#	primers, assembled
	#
	tmp		<- subset(pngs, grepl('PRIMER',STAT) & LABEL!='unassembled')[, list(LABEL=LABEL, PCA= round(N/sum(N),d=3)), by=c('STAT','COHORT')] 
	tmp		<- dcast.data.table(tmp, STAT+LABEL~COHORT, value.var='PCA')
	
}
##--------------------------------------------------------------------------------------------------------
##	olli 25.07.16
##--------------------------------------------------------------------------------------------------------
treecomparison.explaingaps.mutationspectrum.160725<- function()
{
	require(ape)
	require(scales)
	require(data.table)
	require(Hmisc)
	
	wdir			<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/explaingaps'
	wfile			<- 'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.rda'
	load(file.path(wdir, wfile))	
	#
	#	compare mutations of taxa that amplified to those that did not
	#	
	z		<- merge(dgd, dm, by='TAXA')	
	dsmut	<- subset(z, GENE%in%c('1R-3F-firsthalf','1R-3F-secondhalf','2F-1R','2R-4F','4F-3R-firsthalf') & UNASS<0.1 , select=c(TAXA, GENE, COHORT))
	dsmut[, TYPE:= 'Coverage_>90pc']
	tmp		<- subset(z, GENE%in%c('1R-3F-firsthalf','1R-3F-secondhalf','2F-1R','2R-4F','4F-3R-firsthalf') & UNASS>0.95 , select=c(TAXA, GENE, COHORT))
	tmp[, TYPE:= 'Coverage_<5pc']
	dsmut	<- rbind(dsmut, tmp)
	tmp		<- subset(z, GENE%in%c('1R-3F-firsthalf','1R-3F-secondhalf','2F-1R','2R-4F','4F-3R-firsthalf') & UNASS>0.95 & COMET_Region1=='A1', select=c(TAXA, GENE, COHORT))
	tmp[, TYPE:= 'Coverage_<5pc_A1']
	dsmut	<- rbind(dsmut, tmp)
	tmp		<- subset(z, GENE%in%c('1R-3F-firsthalf','1R-3F-secondhalf','2F-1R','2R-4F','4F-3R-firsthalf') & UNASS>0.95 & COMET_Region1=='D', select=c(TAXA, GENE, COHORT))
	tmp[, TYPE:= 'Coverage_<5pc_D']
	dsmut	<- rbind(dsmut, tmp)
	tmp		<- subset(z, GENE%in%c('1R-3F-firsthalf','1R-3F-secondhalf','2F-1R','2R-4F','4F-3R-firsthalf') & UNASS>0.95 & COMET_Region1=='C', select=c(TAXA, GENE, COHORT))
	tmp[, TYPE:= 'Coverage_<5pc_C']
	dsmut	<- rbind(dsmut, tmp)
	#	merge with primer mutations
	dsmut	<- merge(dsmut, subset(dpand, PR%in%c('1R','2F','2R','3F','4F')), by='TAXA',allow.cartesian=TRUE)
	#	
	dsmut	<- subset(dsmut, 	(GENE=='1R-3F-firsthalf' & PR=='2F') | 
								(GENE=='1R-3F-secondhalf' & PR=='2R') |
								(GENE=='2F-1R' & PR=='1R') |
								(GENE=='2R-4F' & PR=='3F') |
								(GENE=='4F-3R-firsthalf' & PR=='4F'))
	set(dsmut, NULL, 'POS', dsmut[,as.integer(gsub('PR_','',POS))])
	set(dsmut, NULL, 'PR', dsmut[, paste('PR_',PR, sep='')])
	#	plot	
	tmp		<- dsmut[,	{
							z	<- round(as.numeric(binconf(length(which(NT_DIFF==1)), length(which(!is.na(NT_DIFF))))), d=3)
							list(EST=c('central','l95','u95'), VAL= z)				
						}, by=c('GENE','TYPE','COHORT','PR','POS')]
	tmp		<- dcast.data.table(tmp, PR+POS+TYPE+COHORT~EST, value.var='VAL')
	tmp		<- merge(tmp, as.data.table(expand.grid(PR=tmp[, unique(PR)], POS=tmp[, unique(POS)], TYPE=tmp[, unique(TYPE)], COHORT=tmp[, unique(COHORT)], stringsAsFactors=FALSE)), by=c('PR','POS','TYPE','COHORT'), all=1, allow.cartesian=TRUE)
	set(tmp, tmp[, which(is.na(central) | central<=0.001)],'l95',NA_real_)
	set(tmp, tmp[, which(is.na(central) | central<=0.001)],'u95',NA_real_)
	set(tmp, tmp[, which(is.na(central) | central<=0.001)],'central',0)
	#set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, levels=c("PANGEA_All","Rakai_All","Rakai_Coverage_>90pc","Rakai_Coverage_<5pc","Rakai_VL_<2e4_Coverage_>90pc", "Rakai_VL_>4e4_Coverage_<5pc","Rakai_D_Coverage_<5pc", "Rakai_A1_Coverage_<5pc", "Rakai_C_Coverage_<5pc"))])	
	#ggplot(subset(tmp, TYPE%in%c('PANGEA_All',"Rakai_Coverage_>90pc","Rakai_Coverage_<5pc","Rakai_D_Coverage_<5pc", "Rakai_A1_Coverage_<5pc", "Rakai_C_Coverage_<5pc")), aes(x=POS, fill=TYPE)) +
	setkey(tmp, COHORT, TYPE)
	set(tmp, NULL, 'COHORT', tmp[, factor(COHORT, levels=c("AC_Resistance","BW-Mochudi","RCCS","UG-MRC"), labels=c("South Africa\nResistance Cohort","Mochudi Prevention\nProject","Rakai Community\nCohort Study","Uganda\nMRC"))])
	ggplot(subset(tmp, TYPE%in%c('Coverage_>90pc','Coverage_<5pc') & COHORT!="South Africa\nResistance Cohort"), aes(x=POS, fill=TYPE)) + 
			geom_bar(aes(y=central), stat='identity', width=0.7, position=position_dodge(0.8)) + 
			facet_grid(PR+COHORT~., scales='free_y') + 
			geom_errorbar(aes(ymin= l95, ymax=u95), position=position_dodge(0.8), width=0.3, size=0.4) +
			theme_bw() + theme(legend.position='bottom') +
			scale_x_continuous(breaks=tmp[, seq_len(max(POS))]) +
			scale_y_continuous(labels=percent) +
			scale_fill_brewer(palette='Set1') +
			scale_alpha_manual(values=c('Coverage_<5pc'=1, 'Coverage_>90pc'=0.7)) +
			#coord_cartesian(ylim=c(0,.2)) +
			labs(x='\nNucleotide position in primer (always forward sense)', y='proportion of assembled sequences with mutation from primer\n', fill='group of\nsequences') 
	ggsave(file=file.path(wdir,gsub('.rda','_gapsprimers_1R2F2R3F4F_eval1.pdf',wfile)), h=25, w=10, limitsize = FALSE)
	
	ggplot(subset(tmp, PR%in%c('PR_2F','PR_2R') & TYPE%in%c('Coverage_>90pc','Coverage_<5pc_A1','Coverage_<5pc_C','Coverage_<5pc_D') & COHORT!="Africa Centre\nResistance Cohort"), aes(x=POS, fill=TYPE)) + 
			geom_bar(aes(y=central), stat='identity', width=0.7, position=position_dodge(0.8)) + 
			facet_grid(PR+COHORT~., scales='free_y') + 
			geom_errorbar(aes(ymin= l95, ymax=u95), position=position_dodge(0.8), width=0.3, size=0.4) +
			theme_bw() + theme(legend.position='bottom') +
			scale_x_continuous(breaks=tmp[, seq_len(max(POS))]) +
			scale_y_continuous(labels=percent) +
			scale_fill_brewer(palette='Set1') +			
			#coord_cartesian(ylim=c(0,.2)) +
			labs(x='\nNucleotide position in primer (always forward sense)', y='proportion of assembled sequences with mutation from primer\n', fill='group of\nsequences') 
	ggsave(file=file.path(wdir,gsub('.rda','_gapsprimers_2F2R_eval2.pdf',wfile)), h=10, w=10, limitsize = FALSE)
	#
	#
	#	old stuff
	#
	#
	
	dsmut	<- subset(dpand,	PR=='1R' & POS=='PR_1' &	
					UNASS_HALF_INDIR_P<0.1 &
					RECENTVL<2e4 & abs(SAMPLEDATE-RECENTVLDATE)<.5, select=c(TAXA, PR, UNASS_HALF_INDIR_P, UNASS_TO_NEXTPRIMER_P))
	dsmut[, TYPE:= 'Rakai_VL_<2e4_Coverage_>90pc']
	tmp		<- subset(dpand,	PR=='3F' & POS=='PR_1' &	
					UNASS_HALF_INDIR_P<0.1 &
					RECENTVL<2e4 & abs(SAMPLEDATE-RECENTVLDATE)<.5, select=c(TAXA, PR, UNASS_HALF_INDIR_P, UNASS_TO_NEXTPRIMER_P))
	tmp[, TYPE:= 'Rakai_VL_<2e4_Coverage_>90pc']		
	dsmut	<- rbind(dsmut, tmp)
	#	high viral load and no coverage
	tmp		<- subset(dpand,	PR=='1R' & 	POS=='PR_1' &
					UNASS_HALF_INDIR_P>0.95 &
					(is.na(ARTSTART) | SAMPLEDATE<ARTSTART) &
					(!everSelfReportArt | everSelfReportArt & SAMPLEDATE<FirstSelfReportArt) &
					RECENTVL>4e4 & abs(SAMPLEDATE-RECENTVLDATE)<.5, select=c(TAXA, PR, UNASS_HALF_INDIR_P, UNASS_TO_NEXTPRIMER_P))
	tmp[, TYPE:= 'Rakai_VL_>4e4_Coverage_<5pc']		
	dsmut	<- rbind(dsmut, tmp)
	tmp		<- subset(dpand,	PR=='3F' & 	POS=='PR_1' &
					UNASS_HALF_INDIR_P>0.95 &
					(is.na(ARTSTART) | SAMPLEDATE<ARTSTART) &
					(!everSelfReportArt | everSelfReportArt & SAMPLEDATE<FirstSelfReportArt) &
					RECENTVL>4e4 & abs(SAMPLEDATE-RECENTVLDATE)<.5, select=c(TAXA, PR, UNASS_HALF_INDIR_P, UNASS_TO_NEXTPRIMER_P))
	tmp[, TYPE:= 'Rakai_VL_>4e4_Coverage_<5pc']	
	dsmut	<- rbind(dsmut, tmp)
	#	coverage
	tmp		<- subset(dpand,	PR=='1R' & POS=='PR_1' & !is.na(STUDY_ID) &
					UNASS_HALF_INDIR_P<0.1, select=c(TAXA, PR, UNASS_HALF_INDIR_P, UNASS_TO_NEXTPRIMER_P))
	tmp[, TYPE:= 'Rakai_Coverage_>90pc']			
	dsmut	<- rbind(dsmut, tmp)
	tmp		<- subset(dpand,	PR=='1R' & 	POS=='PR_1' &
					UNASS_HALF_INDIR_P>0.95, select=c(TAXA, PR, UNASS_HALF_INDIR_P, UNASS_TO_NEXTPRIMER_P))
	tmp[, TYPE:= 'Rakai_Coverage_<5pc']			
	dsmut	<- rbind(dsmut, tmp)
	tmp		<- subset(dpand,	PR=='3F' & POS=='PR_1' & !is.na(STUDY_ID) &
					UNASS_HALF_INDIR_P<0.1, select=c(TAXA, PR, UNASS_HALF_INDIR_P, UNASS_TO_NEXTPRIMER_P))
	tmp[, TYPE:= 'Rakai_Coverage_>90pc']			
	dsmut	<- rbind(dsmut, tmp)
	tmp		<- subset(dpand,	PR=='3F' & 	POS=='PR_1' &
					UNASS_HALF_INDIR_P>0.95, select=c(TAXA, PR, UNASS_HALF_INDIR_P, UNASS_TO_NEXTPRIMER_P))
	tmp[, TYPE:= 'Rakai_Coverage_<5pc']			
	dsmut	<- rbind(dsmut, tmp)	
	#	population
	tmp		<- subset(dpand,	PR=='1R' & 	POS=='PR_1' & !is.na(PANGEA_ID), select=c(TAXA, PR, UNASS_HALF_INDIR_P, UNASS_TO_NEXTPRIMER_P))
	tmp[, TYPE:= 'PANGEA_All']
	dsmut	<- rbind(dsmut, tmp)
	tmp		<- subset(dpand,	PR=='1R' & 	POS=='PR_1' & !is.na(STUDY_ID), select=c(TAXA, PR, UNASS_HALF_INDIR_P, UNASS_TO_NEXTPRIMER_P))
	tmp[, TYPE:= 'Rakai_All']
	dsmut	<- rbind(dsmut, tmp)
	tmp		<- subset(dpand,	PR=='3F' & 	POS=='PR_1' & !is.na(PANGEA_ID), select=c(TAXA, PR, UNASS_HALF_INDIR_P, UNASS_TO_NEXTPRIMER_P))
	tmp[, TYPE:= 'PANGEA_All']
	dsmut	<- rbind(dsmut, tmp)
	tmp		<- subset(dpand,	PR=='3F' & 	POS=='PR_1' & !is.na(STUDY_ID), select=c(TAXA, PR, UNASS_HALF_INDIR_P, UNASS_TO_NEXTPRIMER_P))
	tmp[, TYPE:= 'Rakai_All']
	dsmut	<- rbind(dsmut, tmp)	
	#	subtypes
	tmp		<- subset(dpand,	PR=='1R' & POS=='PR_1' & COMET_Region1=='A1' &
					UNASS_HALF_INDIR_P>0.95, select=c(TAXA, PR, UNASS_HALF_INDIR_P, UNASS_TO_NEXTPRIMER_P))
	tmp[, TYPE:= 'Rakai_A1_Coverage_<5pc']			
	dsmut	<- rbind(dsmut, tmp)
	tmp		<- subset(dpand,	PR=='1R' & POS=='PR_1' & COMET_Region1=='C' &
					UNASS_HALF_INDIR_P>0.95, select=c(TAXA, PR, UNASS_HALF_INDIR_P, UNASS_TO_NEXTPRIMER_P))
	tmp[, TYPE:= 'Rakai_C_Coverage_<5pc']			
	dsmut	<- rbind(dsmut, tmp)
	tmp		<- subset(dpand,	PR=='1R' & POS=='PR_1' & COMET_Region1=='D' &
					UNASS_HALF_INDIR_P>0.95, select=c(TAXA, PR, UNASS_HALF_INDIR_P, UNASS_TO_NEXTPRIMER_P))
	tmp[, TYPE:= 'Rakai_D_Coverage_<5pc']			
	dsmut	<- rbind(dsmut, tmp)
	tmp		<- subset(dpand,	PR=='3F' & POS=='PR_1' & COMET_Region1=='A1' &
					UNASS_HALF_INDIR_P>0.95, select=c(TAXA, PR, UNASS_HALF_INDIR_P, UNASS_TO_NEXTPRIMER_P))
	tmp[, TYPE:= 'Rakai_A1_Coverage_<5pc']			
	dsmut	<- rbind(dsmut, tmp)
	tmp		<- subset(dpand,	PR=='3F' & POS=='PR_1' & COMET_Region1=='C' &
					UNASS_HALF_INDIR_P>0.95, select=c(TAXA, PR, UNASS_HALF_INDIR_P, UNASS_TO_NEXTPRIMER_P))
	tmp[, TYPE:= 'Rakai_C_Coverage_<5pc']			
	dsmut	<- rbind(dsmut, tmp)
	tmp		<- subset(dpand,	PR=='3F' & POS=='PR_1' & COMET_Region1=='D' &
					UNASS_HALF_INDIR_P>0.95, select=c(TAXA, PR, UNASS_HALF_INDIR_P, UNASS_TO_NEXTPRIMER_P))
	tmp[, TYPE:= 'Rakai_D_Coverage_<5pc']			
	dsmut	<- rbind(dsmut, tmp)
	
	tmp		<- merge(subset(dpand, PR=='2F', select=c(TAXA, PR, POS, NT_DIFF, LOC, COMM_NUM, HH_NUM, SEX, AGE, COMET_Region1)), subset(dsmut, PR=='1R', c(TAXA, UNASS_HALF_INDIR_P, TYPE)), by=c('TAXA'), allow.cartesian=TRUE)	
	dsmut	<- rbind(tmp, merge(subset(dpand, PR=='2R', select=c(TAXA, PR, POS, NT_DIFF, LOC, COMM_NUM, HH_NUM, SEX, AGE, COMET_Region1)), subset(dsmut, PR=='3F', c(TAXA, UNASS_HALF_INDIR_P, TYPE)), by=c('TAXA'), allow.cartesian=TRUE))
	
	set(dsmut, NULL, 'POS', dsmut[,gsub('PR_','',POS)])
	set(dsmut, NULL, 'PR', dsmut[, paste('PR_',PR, sep='')])
	dsmut	<- subset(dsmut, !is.na(POS))
	#dcast.data.table(dsmut, TYPE+TAXA+PR1R_UNASS_TO_NEXTPRIMER_P+LOC+COMM_NUM+HH_NUM+SEX+AGE+COMET_Region1  ~  PR+POS, value.var='NT_DIFF')
	
	tmp		<- dsmut[,{
				z	<- round(as.numeric(binconf(length(which(NT_DIFF==1)), length(which(!is.na(NT_DIFF))))), d=3)
				list(EST=c('central','l95','u95'), VAL= z)				
			}, by=c('TYPE','PR','POS')]
	tmp		<- dcast.data.table(tmp, PR+POS+TYPE~EST, value.var='VAL')
	set(tmp, NULL, 'POS', tmp[, as.integer(POS)])
	set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, levels=c("PANGEA_All","Rakai_All","Rakai_Coverage_>90pc","Rakai_Coverage_<5pc","Rakai_VL_<2e4_Coverage_>90pc", "Rakai_VL_>4e4_Coverage_<5pc","Rakai_D_Coverage_<5pc", "Rakai_A1_Coverage_<5pc", "Rakai_C_Coverage_<5pc"))])
	
	#ggplot(subset(tmp, TYPE%in%c('PANGEA_All',"Rakai_Coverage_>90pc","Rakai_Coverage_<5pc","Rakai_D_Coverage_<5pc", "Rakai_A1_Coverage_<5pc", "Rakai_C_Coverage_<5pc")), aes(x=POS, fill=TYPE)) +
	ggplot(subset(tmp, TYPE%in%c('PANGEA_All',"Rakai_Coverage_>90pc","Rakai_Coverage_<5pc")), aes(x=POS, fill=TYPE)) + 
			geom_bar(aes(y=central), stat='identity', width=0.7, position=position_dodge(0.8)) + 
			facet_grid(PR~.) + 
			geom_linerange(aes(ymin= l95, ymax=u95), position=position_dodge(0.8)) +
			theme_bw() + theme(legend.position='bottom') +
			scale_x_continuous(breaks=tmp[, seq_len(max(POS))]) +
			scale_y_continuous(labels=percent, expand=c(0,0)) +
			coord_cartesian(ylim=c(0,.2)) +
			labs(x='\nNucleotide position in primer\n(always forward sense)', y='PANGEA sequences with mutation from primer\n', fill='selected sequences') 
	ggsave(file=file.path(wdir,gsub('.rda','_gapsprimers_2F2R_eval1.pdf',wfile)), h=8, w=10, limitsize = FALSE)
	
	
	ggplot(subset(tmp, TYPE%in%c('PANGEA_All',"Rakai_Coverage_<5pc", "Rakai_VL_>4e4_Coverage_<5pc")), aes(x=POS, fill=TYPE)) + 
			geom_bar(aes(y=central), stat='identity', width=0.7, position=position_dodge(0.8)) + 
			facet_grid(PR~.) + 
			geom_linerange(aes(ymin= l95, ymax=u95), position=position_dodge(0.8)) +
			theme_bw() + theme(legend.position='bottom') +
			scale_x_continuous(breaks=tmp[, seq_len(max(POS))]) +
			scale_y_continuous(labels=percent, expand=c(0,0)) +
			coord_cartesian(ylim=c(0,.5)) +
			labs(x='\nNucleotide position in primer\n(always forward sense)', y='PANGEA sequences with mutation from primer\n', fill='selected sequences') 
	ggsave(file=file.path(wdir,gsub('.rda','_gapsprimers_2F2R_eval2.pdf',wfile)), h=10, w=12, limitsize = FALSE)
	#
	#	plot selected data sets just to make sure I selected correctly
	#	
	tmp		<- unique(subset(dsmut, TYPE%in%c('Rakai_VL_high_Coverage_none','Rakai_VL_low_Coverage_high'), select=c(TAXA, TYPE)))
	chr		<- merge(ch, tmp, by='TAXA', all.x=1)			
	setkey(chr, TAXA)	
	tmp		<- unique(chr)
	setkey(tmp, TYPE, COVP, TAXA)	
	tmp		<- tmp[, list(TAXA=TAXA, PLOT=TYPE, PLOT_ID=seq_along(TAXA)), by='TYPE']
	chr		<- merge(chr, subset(tmp, select=c(TAXA, PLOT_ID)), by='TAXA')	
	ggplot(chr) +
			geom_segment(aes(y=PLOT_ID, yend=PLOT_ID, x=POS_CH, xend=POS_CH+REP_CH-1L, colour=TYPE)) +  
			geom_rect(data=dpani, aes(xmin=START, xmax=END, ymin=-Inf, ymax=Inf), fill="black") +			
			facet_wrap(~TYPE, scales='free_y', ncol=6) +
			scale_x_continuous(expand=c(0,0), breaks=dpani$START, labels=dpani$PR) +
			scale_y_continuous(expand=c(0,0)) +
			scale_colour_brewer(palette='Dark2') +						
			labs(x='\nalignment position', y='Rakai PANGEA-HIV sequences\n', colour='region') +
			theme_bw() +
			theme(	legend.position='bottom', strip.text= element_blank(), strip.background=element_blank()) +
			guides(colour=guide_legend(override.aes=list(size=5)))
}
##--------------------------------------------------------------------------------------------------------
##	olli 25.07.16
##--------------------------------------------------------------------------------------------------------
treecomparison.explaingaps.regressions.160804<- function()
{
	require(ape)	
	require(data.table)	
	
	wdir			<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/explaingaps'
	wfile			<- 'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.rda'
	load(file.path(wdir, wfile))
	#	find reference batch for three gene regions '1R-3F-firsthalf' '1R-3F-secondhalf' '1R-3F'
	if(0)
	{
		tmp		<- subset(dm, select=c(TAXA, SANGER_ID))
		tmp[, BATCH:= regmatches(SANGER_ID,regexpr('^[0-9]+', SANGER_ID))]
		tmp		<- merge(dgd, tmp, by='TAXA')
		dbatch	<- subset(tmp, GENE%in%c('1R-3F','1R-3F-firsthalf','1R-3F-secondhalf'))[, list(ASS_AVG= mean(1-UNASS)), by=c('GENE','BATCH')]
		
		dbatch	<- dbatch[, {
					z<- sort(abs(ASS_AVG-0.66),index.return=TRUE)$ix[1:10]				
					list(REFBATCH= BATCH[z], REF_ASS_AVG=ASS_AVG[z])
				}, by='GENE']
		dbatch	<- dcast.data.table(dbatch, REFBATCH~GENE, value.var='REF_ASS_AVG')
		#	14683 is among the 10 best for all gene regions -- take this one (from BW)
		
	}		
	#	calculate total mutational distance in all primers
	dr		<- dpand[, list(NT_DIFF_SUM=sum(NT_DIFF)), by=c('TAXA','PR')]
	#dr[, length(TAXA), by='PR']
	set(dr, NULL, 'PR', dr[, paste('NT_DIFF_',PR,sep='')])
	#	merge with gene regions
	tmp	<- subset(dgd, grepl('half',GENE) | (!grepl('GAG|POL|ENV',GENE) & !grepl('half',GENE)))
	set(tmp, NULL, c('START','END','ACTG'), NULL)	
	dr		<- merge(dr, tmp, by='TAXA',allow.cartesian=TRUE)
	#	select only meaningful comparisons between primer and gene region
	dr		<- subset(dr, 	(PR%in%c("NT_DIFF_2R","NT_DIFF_2F")&GENE=='1R-3F')	| 
							(PR=="NT_DIFF_1R"&GENE=='2F-1R') |
							(PR=="NT_DIFF_3F"&GENE=='2R-4F')) #| 
							#(PR=="NT_DIFF_3R"&GENE=='4F-3R-secondhalf'))
	#
	#	prepare regression factors
	#	
	tmp		<- subset(dm, select=c(PANGEA_ID, TAXA, STUDY_ID, EXTRACT_ID, SANGER_ID, SAMPLEDATE, ARTSTART, selfReportArt, everSelfReportArt, FirstSelfReportArt, CURRENTLYONART, RECENTVL, RECENTVLDATE, COMET_1F1R, COMET_3F4F, COMET_4F3R, COMET_CONS, COMET_1F1R_N, COMET_3F4F_N, COMET_4F3R_N, COMET_CONS_N, COHORT, LOC, COMM_NUM))
	set(tmp, NULL, 'RECENTVL', tmp[, as.numeric(RECENTVL)])
	dr		<- merge(dr, tmp, by='TAXA')
	#	batches
	dr[, BATCH:=NA_character_]
	tmp		<- dr[, which(!is.na(SANGER_ID))]
	set(dr, tmp, 'BATCH', dr[tmp, regmatches(SANGER_ID,regexpr('^[0-9]+', SANGER_ID))])
	dr[, BATCH2:= as.integer(BATCH)]
	#	add locations to BATCH id
	tmp		<- dr[, list(BATCH_LABEL=paste(BATCH, ' (', paste(unique(COHORT),collapse='+'), ')',sep='')), by='BATCH']
	dr		<- merge(dr, tmp, by='BATCH')
	dr[, BATCH:=NULL]
	setnames(dr, 'BATCH_LABEL', 'BATCH')
	#	to every batch add one 0 and one 1
	if(1)
	{
		tmp		<- dr[, list(TAXA=c('PRIOR1','PRIOR2'), PANGEA_ID=c('PRIOR1','PRIOR2'), STUDY_ID=c('PRIOR1','PRIOR2'), SANGER_ID=c('PRIOR1','PRIOR2'), UNASS=c(0,1), COHORT=COHORT[1] ), by=c('BATCH','PR')]
		dr		<- rbind(dr, tmp, use.names=TRUE, fill=TRUE)			
	}
	#	average viral load per batch
	tmp		<- dr[, {
				ans		<- NA_character_
				tmp		<- which((COHORT=='BW-Mochudi' & !is.na(RECENTVL)) | abs(RECENTVLDATE-SAMPLEDATE)<1)
				if(length(tmp)<length(COHORT)*.5)
					ans	<- 'No VL measured'
				tmp		<- as.character(cut(median(RECENTVL[tmp]), breaks=c(0, 1e4, 2e4, 4e4, 1e5, Inf), labels=c('<1e4','1e4-2e4','2e4-4e4','4e4-1e5','>1e5')))
				if(!is.na(tmp) & is.na(ans))
					ans	<- tmp
				list(BATCHVL= ans)
			}, by=c('BATCH','PR')]
	dr		<- merge(dr, tmp, by=c('BATCH','PR'))
	#	ART status
	dr[, ART:= as.integer(ARTSTART<SAMPLEDATE)]
	set(dr, dr[, which(is.na(ART))], 'ART', 0L)
	set(dr, dr[, which(ART==0 & everSelfReportArt==1 & SAMPLEDATE<FirstSelfReportArt)], 'ART', 1L)
	set(dr, dr[, which(ART==0 & COHORT=='BW-Mochudi' & CURRENTLYONART=='Y')], 'ART', 1L)
	set(dr, dr[, which(ART==0 & COHORT=='AC_Resistance')], 'ART', 0L)
	set(dr, NULL, 'ART', dr[, factor(ART, levels=c(0L,1L), labels=c('no ART', 'ART started or self reported'))])
	#	recent viral load per individual
	dr[, VL:='No VL measured']
	tmp		<- dr[, which((COHORT=='BW-Mochudi' & !is.na(RECENTVL)) | abs(RECENTVLDATE-SAMPLEDATE)<1)]
	set(dr, tmp, 'VL', dr[tmp, cut(RECENTVL, breaks=c(-1, 1e4, 2e4, 4e4, 1e5, Inf), labels=c('<1e4','1e4-2e4','2e4-4e4','4e4-1e5','>1e5'))])
	#	subtype
	#	keep as is
	#set(dr, dr[, which(COMET_Region1%in%c('B','C'))], 'ST', 'B or C')
	#	nucleotide mutations
	dr[, NT_DIFFc:= as.character(NT_DIFF_SUM)] 
	set(dr, dr[, which(is.na(NT_DIFFc))], 'NT_DIFFc', 'Unassembled')
	set(dr, dr[, which(!NT_DIFFc%in%c('0','Unassembled'))], 'NT_DIFFc', 'at least one mutation')
	#	primer assembled
	dr[, ASSEMBLED:= 'Yes'] 
	set(dr, dr[, which(is.na(NT_DIFF_SUM))], 'ASSEMBLED', 'No')		
	#	primers
	set(dr, NULL, 'PR', dr[, gsub('NT_DIFF_','',PR)])
	#	primers + nt_diff
	set(dr, NULL, 'PR_NTDIFFc', dr[, paste(PR,'-',NT_DIFFc,sep='')])
	#	region, community number, household number
	tmp		<- dr[, which(!grepl('PRIOR',TAXA) & is.na(LOC))]
	set(dr, tmp, 'LOC', dr[tmp, COHORT])
	set(dr, NULL, 'COMM_NUM', dr[, as.character(COMM_NUM)])
	tmp		<- dr[, which(!grepl('PRIOR',TAXA) & is.na(COMM_NUM))]
	set(dr, tmp, 'COMM_NUM', dr[tmp, COHORT])
	#	extraction IDs
	set(dr, NULL, 'EXTRACT_ID', dr[,as.integer(gsub('^0+','',EXTRACT_ID))])
	#	sample date
	dr[, SAMPLEDATEc:= cut(SAMPLEDATE, breaks=seq(2010.25, 2015, 0.25), labels=seq(2010.25, 2015-0.25, 0.25))]
	#	amplicon
	dr[, AMPLICON:='two']
	set(dr, dr[, which(PR=='1R')], 'AMPLICON', 'one')
	set(dr, dr[, which(PR=='3F')], 'AMPLICON', 'three')
	#
	#	exclude AfricaCentre because pre-selected after amplification
	#
	dr		<- subset(dr, COHORT!='AC_Resistance')
	#
	#	re-level so that 'presumably good' factors are the reference
	#
	set(dr, NULL, 'AMPLICON', dr[, relevel(factor(AMPLICON), ref='one')])
	set(dr, NULL, 'PR', dr[, relevel(factor(PR), ref='1R')])
	set(dr, NULL, 'PR_NTDIFFc', dr[, relevel(factor(PR_NTDIFFc), ref='1R-0')])
	set(dr, NULL, 'ASSEMBLED', dr[, relevel(factor(ASSEMBLED), ref='Yes')])
	#	defining the reference level for batches is tricky:
	#	if it s a very good batch, then this could simply be down to high VL
	#	if it s an RCCS batch with average viral load, then it could still be down to community?
	#	OK how about we take a batch that performs as well as a typical UG-MRC run, eg 66% on pol?
	#		--> this is '14683'
	set(dr, NULL, 'BATCH', dr[, relevel(factor(BATCH), ref='14683 (BW-Mochudi)')])
	set(dr, NULL, 'COHORT', dr[, relevel(factor(COHORT), ref='BW-Mochudi')])	
	#set(dr, NULL, 'COHORT', dr[, relevel(factor(COHORT), ref='AC_Resistance')])	
	set(dr, NULL, 'VL', dr[, relevel(factor(VL), ref='>1e5')])
	set(dr, NULL, 'BATCHVL', dr[, relevel(factor(BATCHVL), ref='4e4-1e5')])
	set(dr, NULL, 'ART', dr[, relevel(ART, ref='no ART')])
	set(dr, NULL, 'COMET_1F1R', dr[, relevel(factor(COMET_1F1R), ref='A1')])
	set(dr, NULL, 'COMET_3F4F', dr[, relevel(factor(COMET_3F4F), ref='A1')])
	set(dr, NULL, 'COMET_4F3R', dr[, relevel(factor(COMET_4F3R), ref='A1')])
	set(dr, NULL, 'COMET_CONS', dr[, relevel(factor(COMET_CONS), ref='A1')])
	set(dr, NULL, 'NT_DIFFc', dr[, relevel(factor(NT_DIFFc), ref='0')])
	set(dr, NULL, 'LOC', dr[, relevel(factor(LOC), ref='13')])
	#set(dr, NULL, 'COMM_NUM', dr[, relevel(factor(COMM_NUM), ref='106')])
	set(dr, NULL, 'COMM_NUM', dr[, relevel(factor(COMM_NUM), ref='BW-Mochudi')])
	set(dr, NULL, 'SAMPLEDATEc', dr[, relevel(factor(SAMPLEDATEc), ref='2010.25')])
	#
	#	select data
	#
	#	use <60% vs >80%
	dr[, ASS:= as.numeric(as.character(cut(UNASS, breaks=c(-1, 0.6, 0.8, 2), labels=c('1','0.5','0'))))]	
	dr		<- subset(dr, ASS!=.5)
	dr[, UNASS:= 1-ASS]
	#
	#	prepare data sets for logistic regression
	#	
	drs		<- subset(dr, select=c("ASS", "UNASS", "TAXA", "PR", "GENE", "PR_NTDIFFc", "ASSEMBLED", "AMPLICON", "EXTRACT_ID", "PANGEA_ID", "BATCH", "BATCH2", "ART", "BATCHVL", "VL", "COHORT", "LOC", "COMM_NUM", "COMET_1F1R", "COMET_3F4F", "COMET_4F3R", "COMET_CONS", "NT_DIFFc"))	
	drs1r	<- subset(drs, !grepl('PRIOR',TAXA) & PR=='1R' & GENE=='2F-1R')
	set(drs1r, NULL, 'PR_NTDIFFc', drs1r[,as.character(factor(as.character(PR_NTDIFFc), levels=c('1R-0','1R-Unassembled','1R-at least one mutation'),labels=c('1R no mutation','1R unassembled','1R at least one mutation')))])		
	drs3f	<- subset(drs, !grepl('PRIOR',TAXA) & PR=='3F' & GENE=='2R-4F')
	set(drs3f, NULL, 'PR_NTDIFFc', drs3f[,as.character(factor(as.character(PR_NTDIFFc), levels=c('3F-0','3F-Unassembled','3F-at least one mutation'),labels=c('3F no mutation','3F unassembled','3F at least one mutation')))])		
	drs2f	<- subset(drs, PR=='2F' & GENE=='1R-3F-firsthalf')	
	drs2r	<- subset(drs, PR=='2R' & GENE=='1R-3F-secondhalf')	
	tmp		<- subset(dr, is.na(GENE) | GENE=='1R-3F')
	tmp2	<- dcast.data.table(tmp, BATCH+TAXA+ASS+UNASS~PR, value.var='NT_DIFFc')
	setnames(tmp2, c('2F','2R'), c('p2F','p2R'))
	tmp2[, PR_NTDIFFc:='2F or 2R unassembled']
	set(tmp2, tmp2[, which(p2F=='0' & p2R=='0')], 'PR_NTDIFFc', '0')
	set(tmp2, tmp2[, which(p2F=='at least one mutation' & p2R=='0')], 'PR_NTDIFFc', '2F at least one mutation')
	set(tmp2, tmp2[, which(p2F=='0' & p2R=='at least one mutation')], 'PR_NTDIFFc', '2R at least one mutation, 2F no mutation')
	set(tmp2, tmp2[, which(p2F=='at least one mutation' & p2R=='at least one mutation')], 'PR_NTDIFFc', '2F at least one mutation')
	tmp2[, ASSEMBLED:='Yes']
	set(tmp2, tmp2[, which(PR_NTDIFFc=='2F or 2R unassembled')], 'ASSEMBLED', 'No')
	set(tmp2, NULL, 'PR_NTDIFFc', tmp2[, relevel(factor(PR_NTDIFFc), ref='0')])
	set(tmp2, NULL, 'ASSEMBLED', tmp2[, relevel(factor(ASSEMBLED), ref='Yes')])
	tmp2[, AMPLICON:='two']
	tmp		<- subset(tmp, PR=='2F', select=setdiff(colnames(tmp),c('PR','NT_DIFFc','PR_NTDIFFc','UNASS','ASS','NT_DIFF_SUM','ASSEMBLED','GENE','LEN','ACTG','AMPLICON')))
	drs1r	<- merge(tmp, subset(drs1r,select=c(TAXA,BATCH,ASS,UNASS,AMPLICON,PR_NTDIFFc,ASSEMBLED)), by=c('TAXA','BATCH'))
	drs2a	<- merge(tmp, subset(tmp2,select=c(TAXA,BATCH,ASS,UNASS,AMPLICON,PR_NTDIFFc,ASSEMBLED)), by=c('TAXA','BATCH'))	
	drs3f	<- merge(tmp, subset(drs3f,select=c(TAXA,BATCH,ASS,UNASS,AMPLICON,PR_NTDIFFc,ASSEMBLED)), by=c('TAXA','BATCH'))	
	drs12a	<- rbind(drs2a, drs1r, use.names=TRUE)
	drs23a	<- rbind(drs2a, drs3f, use.names=TRUE)	
	drs2	<- subset(drs2a, select=c("ASS", "UNASS", "TAXA", "EXTRACT_ID", "PANGEA_ID", "BATCH", "BATCH2", "PR_NTDIFFc", "AMPLICON", "ASSEMBLED", "ART", "BATCHVL", "VL", "COHORT", "LOC", "COMM_NUM", "COMET_1F1R", "COMET_3F4F", "COMET_4F3R", "COMET_CONS"))
	drs2t	<- subset(drs2a, select=c("ASS", "UNASS", "TAXA", "EXTRACT_ID", "PANGEA_ID", "BATCH", "BATCH2", "PR_NTDIFFc", "AMPLICON","ASSEMBLED", "ART", "BATCHVL", "VL", "COHORT", "LOC", "COMM_NUM", "COMET_1F1R", "COMET_3F4F", "COMET_4F3R", "COMET_CONS", "SAMPLEDATEc", 'SAMPLEDATE'))
	set(drs2t, NULL, 'SAMPLEDATEc', drs2t[, relevel(factor(SAMPLEDATEc), ref='2011.5')])
	drs2npr	<- subset(drs2, !is.na(LOC))
	drs12	<- subset(drs12a, select=c("ASS", "UNASS", "TAXA", "EXTRACT_ID", "PANGEA_ID", "BATCH", "BATCH2", "PR_NTDIFFc", "AMPLICON", "ASSEMBLED", "ART", "BATCHVL", "VL", "COHORT", "LOC", "COMM_NUM", "COMET_1F1R", "COMET_3F4F", "COMET_4F3R", "COMET_CONS"))
	drs12npr<- subset(drs12, !is.na(LOC))
	drs12t	<- subset(drs12a, select=c("ASS", "UNASS", "TAXA", "EXTRACT_ID", "PANGEA_ID", "BATCH", "BATCH2", "PR_NTDIFFc", "AMPLICON","ASSEMBLED", "ART", "BATCHVL", "VL", "COHORT", "LOC", "COMM_NUM", "COMET_1F1R", "COMET_3F4F", "COMET_4F3R", "COMET_CONS", "SAMPLEDATEc", 'SAMPLEDATE'))
	set(drs12t, NULL, 'SAMPLEDATEc', drs12t[, relevel(factor(SAMPLEDATEc), ref='2011.5')])	
	drs23	<- subset(drs23a, select=c("ASS", "UNASS", "TAXA", "EXTRACT_ID", "PANGEA_ID", "BATCH", "BATCH2", "PR_NTDIFFc", "AMPLICON", "ASSEMBLED", "ART", "BATCHVL", "VL", "COHORT", "LOC", "COMM_NUM", "COMET_1F1R", "COMET_3F4F", "COMET_4F3R", "COMET_CONS"))
	drs23npr<- subset(drs23, !is.na(LOC))	
	#
	#	look at the full region 1R-3F <60% vs >80%
	#	
	m2.1		<- glm(data=drs2, UNASS ~ VL + COHORT + PR_NTDIFFc + ART + BATCHVL + BATCH , family='binomial')
	m12.1		<- glm(data=drs12, UNASS ~ VL + COHORT + PR_NTDIFFc + ART + BATCHVL + BATCH , family='binomial')
	#	sub-analyses to calculate AIC
	m2.1b		<- glm(data=drs2npr, UNASS ~ VL + COHORT + PR_NTDIFFc + ART + BATCHVL + BATCH, family='binomial')
	m2.1c		<- glm(data=subset(drs2npr, COMET_1F1R!='check' & VL!='No VL measured'), UNASS ~ VL + COHORT + PR_NTDIFFc + ART + BATCH, family='binomial')
	#	odds ratios
	m2.1.or		<- cbind(data.table(COEF=names(coef(m2.1))), as.data.table( exp(cbind(OR = coef(m2.1), confint(m2.1))) ) )
	setnames(m2.1.or, c('2.5 %','97.5 %'), c('l95','u95'))
	m12.1.or	<- cbind(data.table(COEF=names(coef(m12.1))), as.data.table( exp(cbind(OR = coef(m12.1), confint(m12.1))) ) )
	setnames(m12.1.or, c('2.5 %','97.5 %'), c('l95','u95'))	
	#subset(m2.1.or, l95>1)
	#subset(m12.1.or, l95>1)
	#	exclude sig plates
	tmp			<- subset(m12.1.or, l95>1 & grepl('BATCH',COEF))[, as.integer(regmatches(COEF, regexpr('[0-9]+',COEF)))]
	m2.1e		<- glm(data=subset(drs2npr,!BATCH2%in%tmp), UNASS ~ VL + COHORT + PR_NTDIFFc + ART, family='binomial')
	summary(m2.1e)
	m2.1e.or	<- cbind(data.table(COEF=names(coef(m2.1e))), as.data.table( exp(cbind(OR = coef(m2.1e), confint(m2.1e))) ) )
	setnames(m2.1e.or, c('2.5 %','97.5 %'), c('l95','u95'))
	#	exclude also those that did not assemble	
	m2.1f		<- glm(data=subset(drs2npr,!BATCH2%in%tmp & !grepl('unassembled',PR_NTDIFFc)), UNASS ~ VL + COHORT + PR_NTDIFFc + ART, family='binomial')
	summary(m2.1f)
	m2.1f.or	<- cbind(data.table(COEF=names(coef(m2.1f))), as.data.table( exp(cbind(OR = coef(m2.1f), confint(m2.1f))) ) )
	setnames(m2.1f.or, c('2.5 %','97.5 %'), c('l95','u95'))
	#subset(m2.1e.or, l95>1)
	#	how many batches w significant odds ratio?
	tmp		<- subset(m12.1.or, grepl('BATCH',COEF) & l95>1)	
	nrow(tmp) / ( nrow(subset(m12.1.or, grepl('BATCH',COEF)))+1 )
	#	21.9%  21/25 (84%) from RCCS		
	# 	where from?
	tmp		<- subset(m12.1.or, l95>1 & grepl('BATCH',COEF))[, as.integer(regmatches(COEF, regexpr('[0-9]+',COEF)))]
	subset(dm, as.integer(gsub('_.*','',SANGER_ID))%in%tmp )
	#	I like this model in the end!
	#	- 2F comes out sig and not 2R as in the individual analysis
	#	- UG-MRC worse than RCCS, batches offer better explanation than study for RCCS
	#	- there is a nice series of batches that fail, from 15886 to 15964

	#	repeat with amplicon 1 included
	tmp			<- subset(m12.1.or, l95>1 & grepl('BATCH',COEF))[, as.integer(regmatches(COEF, regexpr('[0-9]+',COEF)))]
	m12.1e		<- glm(data=subset(drs12npr,!BATCH2%in%tmp), UNASS ~ VL + COHORT + PR_NTDIFFc + ART, family='binomial')
	summary(m12.1e)
	m12.1e.or	<- cbind(data.table(COEF=names(coef(m12.1e))), as.data.table( exp(cbind(OR = coef(m12.1e), confint(m12.1e))) ) )
	setnames(m12.1e.or, c('2.5 %','97.5 %'), c('l95','u95'))	
	#m12.1g		<- glm(data=subset(drs12npr,!BATCH2%in%tmp), UNASS ~ VL + COHORT + AMPLICON:ASSEMBLED + ART, family='binomial')
	#summary(m12.1g)
	#m12.1g.or	<- cbind(data.table(COEF=names(coef(m12.1g))), as.data.table( exp(cbind(OR = coef(m12.1g), confint(m12.1g))) ) )
	#setnames(m12.1g.or, c('2.5 %','97.5 %'), c('l95','u95'))	
	#
	#	repeat with amplicon 3 included
	tmp			<- subset(m12.1.or, l95>1 & grepl('BATCH',COEF))[, as.integer(regmatches(COEF, regexpr('[0-9]+',COEF)))]
	m23.1e		<- glm(data=subset(drs23npr,!BATCH2%in%tmp & !grepl('unassembled',PR_NTDIFFc)), UNASS ~ VL + COHORT + PR_NTDIFFc + ART, family='binomial')
	summary(m23.1e)
	m23.1e.or	<- cbind(data.table(COEF=names(coef(m23.1e))), as.data.table( exp(cbind(OR = coef(m23.1e), confint(m23.1e))) ) )
	setnames(m23.1e.or, c('2.5 %','97.5 %'), c('l95','u95'))	
	#	- see if batch VL can be dropped..
	m2.2		<- glm(data=drs2, UNASS ~ VL + COHORT + PR_NTDIFFc + ART + BATCH, family='binomial')
	m2.2b		<- glm(data=drs2npr, UNASS ~ VL + COHORT + PR_NTDIFFc + ART + BATCH, family='binomial')
	summary(m2.2)
	#	Hmmm
	#	- the batch effect from 15892 to 15964 are not present any more
	m2.3		<- glm(data=drs2npr, UNASS ~ VL + COHORT + PR_NTDIFFc + ART, family='binomial')	
	summary(m2.3)
	#	- without batches, RCCS worse than UG-MRC: batches offer better explanation than study for RCCS
	#	- 2R one mutation also not sig
	m2.4		<- glm(data=drs2npr, UNASS ~ VL + PR_NTDIFFc + ART + COMM_NUM, family='binomial')	
	summary(m2.4)
	#
	#	combined model communities + batches. 
	#	just a few communities with borderline significant coefficients
	#	AIC worse than for batch models, but not terribly worse
	#	the problem is that the prior is deleted because communities are missing, and so there are singularities
	m2.6		<- glm(data=drs2, UNASS ~ VL + COMM_NUM + PR_NTDIFFc + ART + BATCHVL + BATCH, family='binomial')
	summary(m2.6)
	tmp			<- subset(drs2t, !grepl('PRIOR',TAXA))
	tmp2		<- subset(m12.1.or, l95>1 & OR<=3 & grepl('BATCH',COEF))[, as.integer(regmatches(COEF, regexpr('[0-9]+',COEF)))]
	tmp2		<- tmp[, which(BATCH2%in%tmp2)]
	set(tmp, tmp2, 'BATCH', tmp[tmp2,paste(BATCH,' *',sep='')])
	tmp2		<- subset(m12.1.or, l95>1 & OR>3 & OR<=10 & grepl('BATCH',COEF))[, as.integer(regmatches(COEF, regexpr('[0-9]+',COEF)))]
	tmp2		<- tmp[, which(BATCH2%in%tmp2)]
	set(tmp, tmp2, 'BATCH', tmp[tmp2,paste(BATCH,' **',sep='')])
	tmp2		<- subset(m12.1.or, l95>1 & OR>10 & grepl('BATCH',COEF))[, as.integer(regmatches(COEF, regexpr('[0-9]+',COEF)))]
	tmp2		<- tmp[, which(BATCH2%in%tmp2)]
	set(tmp, tmp2, 'BATCH', tmp[tmp2,paste(BATCH,' ***',sep='')])
	tmp2		<- tmp[, which(COHORT=='RCCS')]
	set(tmp, tmp2, 'COMM_NUM', tmp[tmp2,paste('Rakai-',COMM_NUM,sep='')])		
	tmp2		<- tmp[, list(EXTRACT_MED=median(EXTRACT_ID)), by=c('BATCH')]	
	setkey(tmp2, EXTRACT_MED)
	set(tmp2, NULL, 'BATCH3', tmp2[, factor(EXTRACT_MED, levels=EXTRACT_MED, labels=BATCH)])
	tmp			<- tmp[, list(N=length(TAXA), P=mean(UNASS==1), SD= ifelse(all(is.na(SAMPLEDATE)), -1L, as.integer(mean(SAMPLEDATE, na.rm=1)<2012.25))), by=c('BATCH','BATCH2','COMM_NUM')]	
	tmp			<- merge(tmp, tmp2, by='BATCH')
	ggplot( tmp, aes(x=BATCH3, y=COMM_NUM, size=N, colour=P, pch=factor(SD, levels=c(-1,0,1),labels=c('Unknown','No','Yes')))) + 
			geom_point() + 
			scale_colour_gradient(low='blue', high='orange') +
			theme_bw() + theme(axis.text.x=element_text(angle=90, vjust=1)) +
			labs(x='sequencing plate (ordered by sample extraction at UCLH)', y='sampling location', colour='proportion of\nsamples with\n<20% assembled sites\nin 1R-3F', size='number of\nsamples', pch='Average\nsample date\nbefore 2012.25')
	ggsave(file=file.path(wdir, 'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_batch_vs_community.pdf'), w=25, h=10, useDingbats=FALSE)
	#
	#	sub analysis: is there an interaction with extraction ID?
	#	
	if(0)
	{
		setkey(drs2npr, EXTRACT_ID)
		drs2npr[, UNASS_ROLL_EXTRACT:= drs2npr[,rollapply(UNASS, width=20, FUN=mean, align="center", partial=TRUE)]]
		tmp			<- subset(m12.1.or, l95>1 & grepl('BATCH',COEF))[, as.integer(regmatches(COEF, regexpr('[0-9]+',COEF)))]
		ggplot(drs2npr, aes(x=EXTRACT_ID, y=BATCH, fill=UNASS_ROLL_EXTRACT, pch=COHORT, colour=factor(BATCH2%in%tmp, levels=c(FALSE,TRUE),labels=c('N','Y'))))  + 
				geom_point(pch=21, size=2, stroke=0.1) + 
				scale_x_continuous(breaks=seq(0,1e4,200)) +
				scale_colour_manual(values=c('Y'='black','N'='transparent')) +
				scale_fill_gradient(low='blue', high='orange') +
				labs(x='\nUCLH extraction ID', y='Sanger plate\n', colour='significant plate effect', fill='rolling mean\nover 20 consecutive samples of\nsequence with >80% not assembled in 1R3F') +
				theme_bw()
		ggsave(file=file.path(wdir, 'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_batch_vs_extraction.pdf'), w=20, h=20)	
		ggplot(drs2npr, aes(x=EXTRACT_ID, y=UNASS_ROLL_EXTRACT, colour=COHORT)) + 
				geom_line() +
				facet_grid(COHORT~.)
	}
	#
	#	sub analysis: is there an interaction between missingVL and cohort ID?
	#
	dvl			<- copy(drs2npr)
	tmp			<- dvl[, which(VL=='No VL measured')]
	set(dvl, tmp, 'VL', dvl[tmp, paste(VL,COHORT,sep='_')])
	set(dvl, NULL, 'VL', dvl[, relevel(factor(as.character(VL)), ref='>1e5')])
	tmp			<- subset(m12.1.or, l95>1 & grepl('BATCH',COEF))[, as.integer(regmatches(COEF, regexpr('[0-9]+',COEF)))]
	m2.9		<- glm(data=subset(dvl,!BATCH2%in%tmp), UNASS ~ VL + COHORT + PR_NTDIFFc + ART, family='binomial')
	summary(m2.9)
	m2.9.or		<- cbind(data.table(COEF=names(coef(m2.9))), as.data.table( exp(cbind(OR = coef(m2.9), confint(m2.9))) ) )
	#
	#	sub analysis: 	does subtype explain plates and study effects?
	#					exclude no VL
	m2.8		<- glm(data=subset(drs2npr, COMET_1F1R!='check' & VL!='No VL measured'), UNASS ~ VL + COMET_1F1R:COHORT + PR_NTDIFFc + ART, family='binomial')
	#	don t use above: there is a strong study effect, and this should not be attributed to  
	m2.8		<- glm(data=subset(drs2npr, COMET_1F1R!='check' & VL!='No VL measured'), UNASS ~ VL + COHORT + COMET_1F1R + PR_NTDIFFc + ART, family='binomial')
	#	with no VL excluded, there is only two cohorts and the model can adjust for study effect and subtype, 
	#	because all BW sequs are C 
	#	short not sig because the mutations are included
	summary(m2.8)	
	sapply(list(m2.1c, m2.8), AIC)	#AIC comparable but this excludes the odd batches.. ..so what s the point?
	sapply(list(m2.1c, m2.8), BIC)	#BIC prefers m2.8!	
	#	sub on RCCS
	#	TODO: include MRC once we have viral loads
	tmp			<- subset(m12.1.or, l95>1 & grepl('BATCH',COEF))[, as.integer(regmatches(COEF, regexpr('[0-9]+',COEF)))]
	drs12st		<- subset(drs12npr, COHORT%in%c('RCCS','UG-MRC','BW-Mochudi') & COMET_1F1R!='short' & !BATCH2%in%tmp)	
	m2.5.1F1R	<- glm(data=drs12st, UNASS ~ VL + COHORT + ART + PR_NTDIFFc + COMET_1F1R, family='binomial')	
	m2.5.3F4F	<- glm(data=subset(drs12npr, COHORT%in%c('RCCS','UG-MRC','BW-Mochudi') & COMET_3F4F!='short' & !BATCH2%in%tmp), UNASS ~ VL + COHORT + ART + PR_NTDIFFc + COMET_3F4F, family='binomial')
	m2.5.4F3R	<- glm(data=subset(drs12npr, COHORT%in%c('RCCS','UG-MRC','BW-Mochudi') & COMET_4F3R!='short' & !BATCH2%in%tmp), UNASS ~ VL + COHORT + ART + PR_NTDIFFc + COMET_4F3R, family='binomial')
	m2.5.CONS	<- glm(data=subset(drs12npr, COHORT%in%c('RCCS','UG-MRC','BW-Mochudi') & COMET_CONS!='short' & !BATCH2%in%tmp), UNASS ~ VL + COHORT + ART + PR_NTDIFFc + COMET_CONS, family='binomial')
	summary(m2.5.1F1R)
	summary(m2.5.3F4F)	
	summary(m2.5.4F3R)	
	summary(m2.5.CONS)	
	tmp			<- m2.5.1F1R
	m2.5.1F1R.or<- cbind(data.table(COEF=names(coef(tmp))), as.data.table( exp(cbind(OR = coef(tmp), confint(tmp))) ) )
	setnames(m2.5.1F1R.or, c('2.5 %','97.5 %'), c('l95','u95'))		
	tmp			<- m2.5.3F4F
	m2.5.3F4F.or<- cbind(data.table(COEF=names(coef(tmp))), as.data.table( exp(cbind(OR = coef(tmp), confint(tmp))) ) )
	setnames(m2.5.3F4F.or, c('2.5 %','97.5 %'), c('l95','u95'))		
	tmp			<- m2.5.4F3R
	m2.5.4F3R.or<- cbind(data.table(COEF=names(coef(tmp))), as.data.table( exp(cbind(OR = coef(tmp), confint(tmp))) ) )
	setnames(m2.5.4F3R.or, c('2.5 %','97.5 %'), c('l95','u95'))
	tmp			<- m2.5.CONS
	m2.5.CONS.or<- cbind(data.table(COEF=names(coef(tmp))), as.data.table( exp(cbind(OR = coef(tmp), confint(tmp))) ) )
	setnames(m2.5.CONS.or, c('2.5 %','97.5 %'), c('l95','u95'))
	
	
	drs2stvl	<- copy(drs2st)
	tmp			<- drs2stvl[, which(VL=='No VL measured')]
	set(drs2stvl, tmp, 'VL', drs2stvl[tmp, paste(VL,COHORT,sep='_')])
	set(drs2stvl, NULL, 'VL', drs2stvl[, relevel(factor(as.character(VL)), ref='>1e5')])	
	m2.10		<- glm(data=drs2stvl, UNASS ~ VL + COHORT + COMET_CONS + PR_NTDIFFc + ART, family='binomial')
	summary(m2.10)
	tmp			<- m2.10
	m2.10.CONS.or<- cbind(data.table(COEF=names(coef(tmp))), as.data.table( exp(cbind(OR = coef(tmp), confint(tmp))) ) )
	setnames(m2.10.CONS.or, c('2.5 %','97.5 %'), c('l95','u95'))			
	#
	#	sub analysis Rakai: sample date
	#	sample date and batch confounded, should not include batch
	#	sort of overfitting. take out comm num?
	#	cool -- this now recovers the sample date effect!!
	m2.7		<- glm(data=subset(drs2t, COHORT=='RCCS'), UNASS ~ VL + SAMPLEDATEc + PR_NTDIFFc + ART + ST, family='binomial')
	summary(m2.7)
	#
	#	compare AIC of batch model vs community model (must use same data)
	#	the batch models have lower AIC than the comm_num model
	if(0)
	{
		sapply(list(m2.1b, m2.2b, m2.3, m2.4), AIC)	
		sapply(list(m2.1b, m2.2b, m2.3, m2.4), AIC)		
	}	
	#
	#	explore if I can use different link functions (log link to get RR)
	#	... not successful
	if(0)
	{
		p2.5		<- glm(data=subset(drs2npr, ST%in%c('A1','D')), UNASS ~ VL + COHORT + PR_NTDIFFc + ART + ST, family=poisson(link=log))
		r2.5		<- glm(data=subset(drs2npr, ST%in%c('A1','D')), UNASS ~ VL + COHORT + PR_NTDIFFc + ART + ST, family=binomial(link=log), control=glm.control(epsilon=1e-8, maxit=100, trace=TRUE), start=c(-2.12,0.93, 0.96, 0.87, 0.72823123, 0.60499165, -0.02665851,0.46523009, 1.07472518, 0.35393831, 0.94179246, 0.03808863,0))
		#qb2.5		<- glm(data=subset(drs2npr, ST%in%c('A1','D')), UNASS ~ VL + COHORT + PR_NTDIFFc + ART + ST, family=quasibinomial(link=log), start=c(-2.12,0.93, 0.96, 0.87, 0.72823123, 0.60499165, -0.02665851,0.46523009, 1.07472518, 0.35393831, 0.94179246, 0.03808863,0))
		#p2.5		<- glm(data=subset(drs2npr, ST%in%c('A1','D')), UNASS ~ VL + COHORT + PR_NTDIFFc + ART + ST, family=poisson(link=log))
		
		library(ResourceSelection)
		hoslem.test(subset(drs2npr, ST%in%c('A1','D'))[, UNASS], fitted(m2.5))
		hoslem.test(subset(drs2npr, ST%in%c('A1','D'))[, UNASS], fitted(p2.5))
		hoslem.test(subset(drs2npr, ST%in%c('A1','D'))[, UNASS], fitted(r2.5))
		sapply(list(m2.5, p2.5, r2.5), AIC)	#the logit model also has smallest deviance and is easiest to fit
		require(AICcmodavg)
		sapply(list(m2.5, p2.5, r2.5), AICc)	#this is basically the same
		p2.2		<- glm(data=drs2, UNASS ~ VL + COHORT + PR_NTDIFFc + ART + BATCH, family=poisson(link=log))
		#	I don t think r2.2 has converged.. 
		r2.2		<- glm(data=drs2, UNASS ~ VL + COHORT + PR_NTDIFFc + ART + BATCH, family=binomial(link=log), control=glm.control(epsilon=1e-8, maxit=1e3, trace=FALSE),
				start=c(-2.41272163456974, 0.557456598211733, 0.486523777039746, 0.385699510360356, 0.319142195398426, 0.875550188558768, 1.03022591978234, 0.513438150979493, rep(0, 114)))									
	}
	if(0)
	{
		subset(drs2a, !grepl('PRIOR',TAXA))[, list(N_ASS= sum(ASS), N_UNASS=sum(1-ASS), P_UNASS=round(sum(1-ASS)/length(ASS), d=2)), by='BATCH']
		subset(drs2a, !grepl('PRIOR',TAXA))[, list(N_ASS= sum(ASS), N_UNASS=sum(1-ASS), P_UNASS=round(sum(1-ASS)/length(ASS), d=2)), by='COMM_NUM']
		subset(drs2a, !grepl('PRIOR',TAXA))[, list(N_ASS= sum(ASS), N_UNASS=sum(1-ASS), P_UNASS=round(sum(1-ASS)/length(ASS), d=2)), by='COMET_Region1']
		subset(drs2a, !grepl('PRIOR',TAXA))[, table(COMM_NUM, floor(SAMPLEDATE))]		
	}
	#
	#	write odds ratio tables
	#
	tmp		<- copy(m12.1e.or)
	set(tmp, NULL, 'OR', tmp[, round(OR,d=2)])
	set(tmp, NULL, 'CI', tmp[, paste(round(l95,d=2),'-',round(u95,d=2),sep='')])	
	write.csv(subset(tmp, !is.na(OR), select=c(COEF, OR, CI)), row.names=FALSE, file=file.path(wdir, 'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_regression_model12.1e.csv'))
	tmp		<- copy(m2.1e.or)
	set(tmp, NULL, 'OR', tmp[, round(OR,d=2)])
	set(tmp, NULL, 'CI', tmp[, paste(round(l95,d=2),'-',round(u95,d=2),sep='')])	
	write.csv(subset(tmp, !is.na(OR), select=c(COEF, OR, CI)), row.names=FALSE, file=file.path(wdir, 'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_regression_model2.1.csv'))
	tmp		<- copy(m2.5.1F1R.or)
	set(tmp, NULL, 'OR', tmp[, round(OR,d=2)])
	set(tmp, NULL, 'CI', tmp[, paste(round(l95,d=2),'-',round(u95,d=2),sep='')])	
	write.csv(subset(tmp, !is.na(OR), select=c(COEF, OR, CI)), row.names=FALSE, file=file.path(wdir, 'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_regression_model2.5.1F1R.or.csv'))
	tmp		<- copy(m2.5.3F4F.or)
	set(tmp, NULL, 'OR', tmp[, round(OR,d=2)])
	set(tmp, NULL, 'CI', tmp[, paste(round(l95,d=2),'-',round(u95,d=2),sep='')])	
	write.csv(subset(tmp, !is.na(OR), select=c(COEF, OR, CI)), row.names=FALSE, file=file.path(wdir, 'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_regression_model2.5.3F4F.csv'))
	tmp		<- copy(m2.5.4F3R.or)
	set(tmp, NULL, 'OR', tmp[, round(OR,d=2)])
	set(tmp, NULL, 'CI', tmp[, paste(round(l95,d=2),'-',round(u95,d=2),sep='')])	
	write.csv(subset(tmp, !is.na(OR), select=c(COEF, OR, CI)), row.names=FALSE, file=file.path(wdir, 'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_regression_model2.5.4F3R.csv'))
	tmp		<- copy(m2.5.CONS.or)
	set(tmp, NULL, 'OR', tmp[, round(OR,d=2)])
	set(tmp, NULL, 'CI', tmp[, paste(round(l95,d=2),'-',round(u95,d=2),sep='')])	
	write.csv(subset(tmp, !is.na(OR), select=c(COEF, OR, CI)), row.names=FALSE, file=file.path(wdir, 'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_regression_model2.5.CONS.csv'))
	
	#
	#	run batches - subanalysis
	#	
	drs2b	<- copy(drs2a)	
	tmp2	<- unique(subset(drs2b, BATCH2>=15892 & BATCH2<=15964 & !is.na(COMM_NUM), COMM_NUM))
	tmp2	<- unique(subset(merge(drs2b, tmp2, by='COMM_NUM'), select=BATCH2))
	drs2b		<- merge(drs2b, tmp2, by='BATCH2')
	drs2b[, UNASSB:= 0L]
	set(drs2b, drs2b[, which(BATCH2>=15892 & BATCH2<=15964 & !is.na(COMM_NUM))], 'UNASSB', 1L)
	drs2b	<- subset(drs2b, !grepl('PRIOR',TAXA), select=c("UNASSB", "TAXA", "PANGEA_ID", "BATCH", "PR_NTDIFFc", "ART", "BATCHVL", "VL", "COHORT", "LOC", "COMM_NUM", "ST", "SAMPLEDATE", "SAMPLEDATEc"))	
	set(drs2b, NULL, 'SAMPLEDATEc', drs2b[, relevel(factor(SAMPLEDATEc), ref='2011.75')])
	write.csv(drs2b, row.names=FALSE, file=file.path(wdir, 'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_batch_vs_community.csv'))
	s2.1	<- glm(data=subset(drs2b, !is.na(SAMPLEDATE)), UNASSB ~ VL + COMM_NUM + PR_NTDIFFc + ART + ST + SAMPLEDATE, family='binomial')
	summary(s2.1)
	s2.2	<- glm(data=subset(drs2b, !is.na(SAMPLEDATE)), UNASSB ~ VL + PR_NTDIFFc + ART + ST + SAMPLEDATEc, family='binomial')
	summary(s2.2)
	#
	#	do the predictions, use model m2.1
	#
	#	exclude singular batches to get OK prediction
	#
	#	data for predictions
	tmp			<- subset(m12.1.or, l95>1.05 & grepl('BATCH',COEF))[, as.integer(regmatches(COEF, regexpr('[0-9]+',COEF)))] 
	dpr			<- subset(drs2npr, COHORT%in%c('RCCS','UG-MRC','BW-Mochudi') & COMET_CONS!='short' & !BATCH2%in%tmp)	
	#	higher viral load
	tmp		<- subset(dpr, VL!='No VL measured')
	pr0		<- predict(m2.5.CONS, newdata=tmp, type='response', se.fit=TRUE)
	set(tmp, tmp[, which(VL%in%c('<1e4'))], 'VL', '1e4-2e4')
	pr1		<- predict(m2.5.CONS, newdata=tmp, type='response', se.fit=TRUE)
	set(tmp, tmp[, which(VL%in%c('<1e4', '1e4-2e4'))], 'VL', '2e4-4e4')
	pr2		<- predict(m2.5.CONS, newdata=tmp, type='response', se.fit=TRUE)
	set(tmp, tmp[, which(VL%in%c('<1e4', '1e4-2e4', '2e4-4e4'))], 'VL', '4e4-1e5')
	pr3		<- predict(m2.5.CONS, newdata=tmp, type='response', se.fit=TRUE)
	set(tmp, tmp[, which(VL%in%c('<1e4', '1e4-2e4', '2e4-4e4','4e4-1e5'))], 'VL', '>1e5')
	pr4		<- predict(m2.5.CONS, newdata=tmp, type='response', se.fit=TRUE)
	tmp[, H0u_p:= UNASS]
	tmp[, H0all_p:= 1]
	tmp[, H0_p:= pr0$fit]
	tmp[, H0_l:= pr0$fit-1.96*pr0$se.fit]
	tmp[, H0_u:= pr0$fit+1.96*pr0$se.fit]
	tmp[, H1_p:= pr1$fit]
	tmp[, H1_l:= pr1$fit-1.96*pr1$se.fit]
	tmp[, H1_u:= pr1$fit+1.96*pr1$se.fit]
	tmp[, H2_p:= pr2$fit]
	tmp[, H2_l:= pr2$fit-1.96*pr2$se.fit]
	tmp[, H2_u:= pr2$fit+1.96*pr2$se.fit]
	tmp[, H3_p:= pr3$fit]
	tmp[, H3_l:= pr3$fit-1.96*pr3$se.fit]
	tmp[, H3_u:= pr3$fit+1.96*pr3$se.fit]
	tmp[, H4_p:= pr4$fit]
	tmp[, H4_l:= pr4$fit-1.96*pr4$se.fit]
	tmp[, H4_u:= pr4$fit+1.96*pr4$se.fit]
	tmp		<- melt(tmp, id.vars=c('COHORT'), measure.vars=colnames(tmp)[grepl('^H[0-9]+',colnames(tmp))])
	set(tmp, tmp[, which(value<0)], 'value', 0)
	set(tmp, tmp[, which(value>1)], 'value', 1)
	set(tmp, NULL, 'TYPE', tmp[, gsub('_','',regmatches(variable, regexpr('_[a-z]', variable)))])
	set(tmp, NULL, 'variable', tmp[, regmatches(variable, regexpr('[A-Za-z0-9]+', variable))])
	tmp		<- tmp[, list(value=sum(value)), by=c('COHORT','TYPE','variable')]
	tmp2	<- subset(tmp, variable=='H0' & TYPE=='p')
	setnames(tmp2, 'value', 'BL')
	tmp		<- merge(tmp, subset(tmp2, select=c(COHORT, BL)), by='COHORT')
	tmp2	<- subset(tmp, variable=='H0all')
	setnames(tmp2, 'value', 'TOTAL')
	tmp		<- merge(tmp, subset(tmp2, select=c(COHORT, TOTAL)), by='COHORT')
	tmp		<- subset(tmp, variable!='H0u' & variable!='H0all' & TYPE=='p')
	tmp[, REDUCTION_N:= BL-value]
	tmp[, REDUCTION_P:= (BL-value)/TOTAL]
	tmp		<- subset(merge(tmp, data.table(variable=c('H1','H2','H3','H4'), LABEL=c('VL >1e4','VL >2e4','VL >4e4','VL >1e5')), by='variable'), select=c(COHORT, LABEL, REDUCTION_N, REDUCTION_P, TOTAL))
	ans		<- copy(tmp)
	#
	#	primer sites assembled
	tmp		<- copy(dpr)
	pr0		<- predict(m2.5.CONS, newdata=tmp, type='response', se.fit=TRUE)
	set(tmp, tmp[, which(PR_NTDIFFc%in%c('2F or 2R unassembled'))], 'PR_NTDIFFc', '0')
	pr1		<- predict(m2.5.CONS, newdata=tmp, type='response', se.fit=TRUE)
	tmp[, H0u_p:= UNASS]
	tmp[, H0all_p:= 1]
	tmp[, H0_p:= pr0$fit]
	tmp[, H0_l:= pr0$fit-1.96*pr0$se.fit]
	tmp[, H0_u:= pr0$fit+1.96*pr0$se.fit]
	tmp[, H1_p:= pr1$fit]
	tmp[, H1_l:= pr1$fit-1.96*pr1$se.fit]
	tmp[, H1_u:= pr1$fit+1.96*pr1$se.fit]
	tmp		<- melt(tmp, id.vars=c('COHORT'), measure.vars=colnames(tmp)[grepl('^H[0-9]+',colnames(tmp))])
	set(tmp, tmp[, which(value<0)], 'value', 0)
	set(tmp, tmp[, which(value>1)], 'value', 1)
	set(tmp, NULL, 'TYPE', tmp[, gsub('_','',regmatches(variable, regexpr('_[a-z]', variable)))])
	set(tmp, NULL, 'variable', tmp[, regmatches(variable, regexpr('[A-Za-z0-9]+', variable))])
	tmp		<- tmp[, list(value=sum(value)), by=c('COHORT','TYPE','variable')]
	tmp2	<- subset(tmp, variable=='H0' & TYPE=='p')
	setnames(tmp2, 'value', 'BL')
	tmp		<- merge(tmp, subset(tmp2, select=c(COHORT, BL)), by='COHORT')
	tmp2	<- subset(tmp, variable=='H0all')
	setnames(tmp2, 'value', 'TOTAL')
	tmp		<- merge(tmp, subset(tmp2, select=c(COHORT, TOTAL)), by='COHORT')
	tmp		<- subset(tmp, variable!='H0u' & variable!='H0all' & TYPE=='p')
	tmp[, REDUCTION_N:= BL-value]
	tmp[, REDUCTION_P:= (BL-value)/TOTAL]
	tmp		<- subset(merge(tmp, data.table(variable=c('H1'), LABEL=c('Assembled primer sites')), by='variable'), select=c(COHORT, LABEL, REDUCTION_N, REDUCTION_P, TOTAL))
	ans		<- rbind(ans, tmp)
	#
	#	no mutations
	tmp		<- subset(dpr, PR_NTDIFFc!='2F or 2R unassembled')	
	pr0		<- predict(m2.5.CONS, newdata=tmp, type='response', se.fit=TRUE)
	set(tmp, tmp[, which(PR_NTDIFFc%in%c('2F at least one mutation'))], 'PR_NTDIFFc', '0')
	pr1		<- predict(m2.5.CONS, newdata=tmp, type='response', se.fit=TRUE)
	set(tmp, tmp[, which(PR_NTDIFFc%in%c('2F at least one mutation','2R at least one mutation, 2F no mutation'))], 'PR_NTDIFFc', '0')
	pr2		<- predict(m2.5.CONS, newdata=tmp, type='response', se.fit=TRUE)
	tmp[, H0u_p:= UNASS]
	tmp[, H0all_p:= 1]
	tmp[, H0_p:= pr0$fit]
	tmp[, H0_l:= pr0$fit-1.96*pr0$se.fit]
	tmp[, H0_u:= pr0$fit+1.96*pr0$se.fit]
	tmp[, H1_p:= pr1$fit]
	tmp[, H1_l:= pr1$fit-1.96*pr1$se.fit]
	tmp[, H1_u:= pr1$fit+1.96*pr1$se.fit]
	tmp[, H2_p:= pr2$fit]
	tmp[, H2_l:= pr2$fit-1.96*pr2$se.fit]
	tmp[, H2_u:= pr2$fit+1.96*pr2$se.fit]
	tmp		<- melt(tmp, id.vars=c('COHORT'), measure.vars=colnames(tmp)[grepl('^H[0-9]+',colnames(tmp))])
	set(tmp, tmp[, which(value<0)], 'value', 0)
	set(tmp, tmp[, which(value>1)], 'value', 1)
	set(tmp, NULL, 'TYPE', tmp[, gsub('_','',regmatches(variable, regexpr('_[a-z]', variable)))])
	set(tmp, NULL, 'variable', tmp[, regmatches(variable, regexpr('[A-Za-z0-9]+', variable))])
	tmp		<- tmp[, list(value=sum(value)), by=c('COHORT','TYPE','variable')]
	tmp2	<- subset(tmp, variable=='H0' & TYPE=='p')
	setnames(tmp2, 'value', 'BL')
	tmp		<- merge(tmp, subset(tmp2, select=c(COHORT, BL)), by='COHORT')
	tmp2	<- subset(tmp, variable=='H0all')
	setnames(tmp2, 'value', 'TOTAL')
	tmp		<- merge(tmp, subset(tmp2, select=c(COHORT, TOTAL)), by='COHORT')
	tmp		<- subset(tmp, variable!='H0u' & variable!='H0all' & TYPE=='p')
	tmp[, REDUCTION_N:= BL-value]
	tmp[, REDUCTION_P:= (BL-value)/TOTAL]
	tmp		<- subset(merge(tmp, data.table(variable=c('H1','H2'), LABEL=c('PR 2F no mutation','PR 2F, 2R no mutation')), by='variable'), select=c(COHORT, LABEL, REDUCTION_N, REDUCTION_P, TOTAL))
	ans		<- rbind(ans, tmp)
	#
	#	subtype (from sub analysis)
	tmp		<- copy(dpr)
	pr0		<- predict(m2.5.CONS, newdata=tmp, type='response', se.fit=TRUE)
	set(tmp, tmp[, which(COMET_CONS%in%c('D'))], 'COMET_CONS', 'A1')
	pr1		<- predict(m2.5.CONS, newdata=tmp, type='response', se.fit=TRUE)
	set(tmp, tmp[, which(COMET_CONS%in%c('D','C'))], 'COMET_CONS', 'A1')
	pr2		<- predict(m2.5.CONS, newdata=tmp, type='response', se.fit=TRUE)
	set(tmp, tmp[, which(COMET_CONS%in%c('D','C','pot_recombinant'))], 'COMET_CONS', 'A1')
	pr3		<- predict(m2.5.CONS, newdata=tmp, type='response', se.fit=TRUE)
	set(tmp, tmp[, which(COMET_CONS%in%c('D','C','other','pot_recombinant'))], 'COMET_CONS', 'A1')
	pr4		<- predict(m2.5.CONS, newdata=tmp, type='response', se.fit=TRUE)
	tmp[, H0u_p:= UNASS]
	tmp[, H0all_p:= 1]
	tmp[, H0_p:= pr0$fit]
	tmp[, H0_l:= pr0$fit-1.96*pr0$se.fit]
	tmp[, H0_u:= pr0$fit+1.96*pr0$se.fit]
	tmp[, H1_p:= pr1$fit]
	tmp[, H1_l:= pr1$fit-1.96*pr1$se.fit]
	tmp[, H1_u:= pr1$fit+1.96*pr1$se.fit]
	tmp[, H2_p:= pr2$fit]
	tmp[, H2_l:= pr2$fit-1.96*pr2$se.fit]
	tmp[, H2_u:= pr2$fit+1.96*pr2$se.fit]
	tmp[, H3_p:= pr3$fit]
	tmp[, H3_l:= pr3$fit-1.96*pr3$se.fit]
	tmp[, H3_u:= pr3$fit+1.96*pr3$se.fit]
	tmp[, H4_p:= pr4$fit]
	tmp[, H4_l:= pr4$fit-1.96*pr4$se.fit]
	tmp[, H4_u:= pr4$fit+1.96*pr4$se.fit]
	tmp		<- melt(tmp, id.vars=c('COHORT'), measure.vars=colnames(tmp)[grepl('^H[0-9]+',colnames(tmp))])
	set(tmp, tmp[, which(value<0)], 'value', 0)
	set(tmp, tmp[, which(value>1)], 'value', 1)
	set(tmp, NULL, 'TYPE', tmp[, gsub('_','',regmatches(variable, regexpr('_[a-z]', variable)))])
	set(tmp, NULL, 'variable', tmp[, regmatches(variable, regexpr('[A-Za-z0-9]+', variable))])
	tmp		<- tmp[, list(value=sum(value)), by=c('COHORT','TYPE','variable')]
	tmp2	<- subset(tmp, variable=='H0' & TYPE=='p')
	setnames(tmp2, 'value', 'BL')
	tmp		<- merge(tmp, subset(tmp2, select=c(COHORT, BL)), by='COHORT')
	tmp2	<- subset(tmp, variable=='H0all')
	setnames(tmp2, 'value', 'TOTAL')
	tmp		<- merge(tmp, subset(tmp2, select=c(COHORT, TOTAL)), by='COHORT')
	tmp		<- subset(tmp, variable!='H0u' & variable!='H0all' & TYPE=='p')
	tmp[, REDUCTION_N:= BL-value]
	tmp[, REDUCTION_P:= (BL-value)/TOTAL]
	tmp		<- subset(merge(tmp, data.table(variable=c('H1','H2','H3','H4'), LABEL=c('D','D,C','D,C,pot recombinant','D,C,other,pot recomb')), by='variable'), select=c(COHORT, LABEL, REDUCTION_N, REDUCTION_P, TOTAL))
	ans		<- rbind(ans, tmp)
	#
	#	no ART self report 
	tmp		<- copy(dpr)
	pr0		<- predict(m2.5.CONS, newdata=tmp, type='response', se.fit=TRUE)
	set(tmp, tmp[, which(ART%in%c('ART started or self reported'))], 'ART', 'no ART')
	pr1		<- predict(m2.5.CONS, newdata=tmp, type='response', se.fit=TRUE)
	tmp[, H0u_p:= UNASS]
	tmp[, H0all_p:= 1]
	tmp[, H0_p:= pr0$fit]
	tmp[, H0_l:= pr0$fit-1.96*pr0$se.fit]
	tmp[, H0_u:= pr0$fit+1.96*pr0$se.fit]
	tmp[, H1_p:= pr1$fit]
	tmp[, H1_l:= pr1$fit-1.96*pr1$se.fit]
	tmp[, H1_u:= pr1$fit+1.96*pr1$se.fit]
	tmp		<- melt(tmp, id.vars=c('COHORT'), measure.vars=colnames(tmp)[grepl('^H[0-9]+',colnames(tmp))])
	set(tmp, tmp[, which(value<0)], 'value', 0)
	set(tmp, tmp[, which(value>1)], 'value', 1)
	set(tmp, NULL, 'TYPE', tmp[, gsub('_','',regmatches(variable, regexpr('_[a-z]', variable)))])
	set(tmp, NULL, 'variable', tmp[, regmatches(variable, regexpr('[A-Za-z0-9]+', variable))])
	tmp		<- tmp[, list(value=sum(value)), by=c('COHORT','TYPE','variable')]
	tmp2	<- subset(tmp, variable=='H0' & TYPE=='p')
	setnames(tmp2, 'value', 'BL')
	tmp		<- merge(tmp, subset(tmp2, select=c(COHORT, BL)), by='COHORT')
	tmp2	<- subset(tmp, variable=='H0all')
	setnames(tmp2, 'value', 'TOTAL')
	tmp		<- merge(tmp, subset(tmp2, select=c(COHORT, TOTAL)), by='COHORT')
	tmp		<- subset(tmp, variable!='H0u' & variable!='H0all' & TYPE=='p')
	tmp[, REDUCTION_N:= BL-value]
	tmp[, REDUCTION_P:= (BL-value)/TOTAL]
	tmp		<- subset(merge(tmp, data.table(variable=c('H1'), LABEL=c('ART not self-reported or started')), by='variable'), select=c(COHORT, LABEL, REDUCTION_N, REDUCTION_P, TOTAL))
	ans		<- rbind(ans, tmp)
	#
	#	RCCS batches as typical UG-MRC 15034 (UG-MRC)
	#	bad batches as typical from same region 15699 (RCCS)		
	tmp		<- subset(drs2npr, COHORT=='RCCS')
	pr0		<- predict(m2.1, newdata=tmp, type='response', se.fit=TRUE)
	tmp2	<- subset(m12.1.or, l95>1 & grepl('BATCH',COEF))[, gsub('BATCH','',COEF)]
	set(tmp, tmp[, which(BATCH%in%tmp2)], 'BATCH', '15699 (RCCS)')
	pr1		<- predict(m2.1, newdata=tmp, type='response', se.fit=TRUE)	
	tmp2	<- subset(subset(drs2a, !grepl('PRIOR',TAXA) & grepl('RCCS',BATCH))[, list(P_UNASS=round(sum(1-ASS)/length(ASS), d=2)), by='BATCH'], P_UNASS>0.38)[, BATCH]	
	set(tmp, tmp[, which(BATCH%in%tmp2)], 'BATCH', '15034 (UG-MRC)')
	pr2		<- predict(m2.1, newdata=tmp, type='response', se.fit=TRUE)
	tmp[, H0u_p:= UNASS]
	tmp[, H0all_p:= 1]
	tmp[, H0_p:= pr0$fit]
	tmp[, H0_l:= pr0$fit-1.96*pr0$se.fit]
	tmp[, H0_u:= pr0$fit+1.96*pr0$se.fit]
	tmp[, H1_p:= pr1$fit]
	tmp[, H1_l:= pr1$fit-1.96*pr1$se.fit]
	tmp[, H1_u:= pr1$fit+1.96*pr1$se.fit]
	tmp[, H2_p:= pr2$fit]
	tmp[, H2_l:= pr2$fit-1.96*pr2$se.fit]
	tmp[, H2_u:= pr2$fit+1.96*pr2$se.fit]
	tmp		<- melt(tmp, id.vars=c('COHORT'), measure.vars=colnames(tmp)[grepl('^H[0-9]+',colnames(tmp))])
	set(tmp, tmp[, which(value<0)], 'value', 0)
	set(tmp, tmp[, which(value>1)], 'value', 1)
	set(tmp, NULL, 'TYPE', tmp[, gsub('_','',regmatches(variable, regexpr('_[a-z]', variable)))])
	set(tmp, NULL, 'variable', tmp[, regmatches(variable, regexpr('[A-Za-z0-9]+', variable))])
	tmp		<- tmp[, list(value=sum(value)), by=c('COHORT','TYPE','variable')]
	tmp2	<- subset(tmp, variable=='H0' & TYPE=='p')
	setnames(tmp2, 'value', 'BL')
	tmp		<- merge(tmp, subset(tmp2, select=c(COHORT, BL)), by='COHORT')
	tmp2	<- subset(tmp, variable=='H0all')
	setnames(tmp2, 'value', 'TOTAL')
	tmp		<- merge(tmp, subset(tmp2, select=c(COHORT, TOTAL)), by='COHORT')
	tmp		<- subset(tmp, variable!='H0u' & variable!='H0all' & TYPE=='p')
	tmp[, REDUCTION_N:= BL-value]
	tmp[, REDUCTION_P:= (BL-value)/TOTAL]
	tmp		<- subset(merge(tmp, data.table(variable=c('H1','H2'), LABEL=c('no sig plates','as avg UG MRC')), by='variable'), select=c(COHORT, LABEL, REDUCTION_N, REDUCTION_P, TOTAL))
	ans		<- rbind(ans, tmp)	
	#
	#	write csv
	set(ans, NULL, 'REDUCTION_P_L', ans[, paste(round(REDUCTION_P*100,d=1),'pc',sep='')])
	#ans[, paste(unique(LABEL),collapse='", "')]
	#set(ans, NULL, 'LABEL', ans[, factor(LABEL, levels=c("Assembled primer sites", "VL >1e4", "VL >2e4", "VL >4e4", "VL >1e5", "D", "Batch run as typical run from same comm", "PR 2F no mutation", "PR 2F, 2R no mutation", "ART not self-reported or started"))])	
	ans		<- dcast.data.table(ans, LABEL~COHORT, value.var='REDUCTION_P_L')		
	write.csv(ans, row.names=FALSE, file=file.path(wdir, 'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_percentreductions.csv'))
	
	
	#
	#	save
	#
	save(drs2, m2.1, m2.2, m2.3, m2.1.or, file=file.path(wdir, gsub('.rda','_logistic_160905.rda',wfile)))
	tmp		<- copy(m2.1.or)
	set(tmp, NULL, 'COEF', tmp[, gsub('COHORT','',gsub('PR_NTDIFFc','',gsub('VL','viral load ',gsub('BATCH','plate ',COEF))))])
	set(tmp, NULL, 'OR', tmp[, round(OR,d=2)])
	set(tmp, NULL, 'LABEL', tmp[, as.character(OR)])
	set(tmp, NULL, 'l95', tmp[, round(l95,d=2)])
	set(tmp, NULL, 'u95', tmp[, round(u95,d=2)])
	set(tmp, NULL, 'LABEL2', tmp[, paste(l95,'-',u95,sep='')])
	tmp2	<- tmp[, which(l95>1 & l95<=3)]
	set(tmp, tmp2, 'LABEL', tmp[tmp2,paste(LABEL,' *',sep='')])
	tmp2	<- tmp[, which(l95>3)]
	set(tmp, tmp2, 'LABEL', tmp[tmp2,paste(LABEL,' **',sep='')])
	write.csv(subset(tmp, !is.na(OR), select=c(COEF,LABEL,LABEL2)), row.names=FALSE, file=file.path(wdir, gsub('.rda','_logistic_160905.csv',wfile)))
	
	#
	#	old stuff
	#

	#
	#	what are these batches?
	#	
	set(drs2a, NULL, 'BATCH', drs2a[, as.integer(regmatches(BATCH,regexpr('[0-9]+',BATCH)))])
	tmp2	<- subset(m2.1.or, grepl('BATCH',COEF) & grepl('(RCCS)',COEF,fixed=1) & l95>1)[, as.integer(regmatches(COEF,regexpr('[0-9]+',COEF)))]
	db		<- subset(drs2a, BATCH %in% tmp2)
	db[, BATCH_UNASS:='Y']
	tmp		<- subset(drs2a, COHORT=='RCCS' & !BATCH%in%tmp)
	tmp[, BATCH_UNASS:='N']
	db		<- rbind(db,tmp)
	#	compare first locations ( % of samples from... )
	tmp		<- db[, {
				#z	<- subset(db, BATCH_UNASS=='Y')[,table(as.character(COMM_NUM))]
				z	<- table(as.character(COMM_NUM))
				ans	<- as.data.table(binconf(z, sum(z)))
				ans[, TYPE:=names(z)]
				setnames(ans, c('PointEst','Lower','Upper'),c('central','l95','u95'))
				ans	<- melt(ans, id.vars='TYPE')
				set(ans, NULL, 'value', ans[,round(value*100,d=1)])
				ans
			}, by='BATCH_UNASS']		
	dcast.data.table(tmp, TYPE~BATCH_UNASS+variable)
	#
	#	batch models with BATCH VL
	#	
	m.1		<- glm(data=drs, ASS ~ VL + COHORT + PR_NTDIFFc + ART + BATCHVL + BATCH, family='binomial')
	summary(m.1)
	#	batch VL not significant, drop..
	m.2		<- glm(data=drs, ASS ~ BATCH + VL + COHORT + PR_NTDIFFc + ART - 1, family='binomial')
	summary(m.2)
	#	with batches, UG-MRC worse than RCCS + PR_NTDIFFc2F/2R mut significant, but also PR_NTDIFFc2R0 significant (simply knock on)
	m.3		<- glm(data=drs, ASS ~ VL + COHORT + PR_NTDIFFc + ART, family='binomial')
	summary(m.3)
	#	without batches, RCCS worse than UG-MRC: batches offer better explanation than study for RCCS
	
	m.2.or<- cbind(data.table(COEF=names(coef(m.2))), as.data.table( exp(cbind(OR = coef(m.2), confint(m.2))) ) )
	setnames(m.2.or, c('2.5 %','97.5 %'), c('l95','u95'))		
	subset(m.2.or, u95<0.95)
	
	#
	#	batch models with BATCH VL
	#	
	m2f.1	<- glm(data=drs2f, ASS ~ VL + COHORT + ART +  NT_DIFFc + BATCHVL + BATCH, family='binomial')
	summary(m2f.1)
	m2f.1.or<- cbind(data.table(COEF=names(coef(m2f.1))), as.data.table( exp(cbind(OR = coef(m2f.1), confint(m2f.1))) ) )
	setnames(m2f.1.or, c('2.5 %','97.5 %'), c('l95','u95'))		
	m2r.1	<- glm(data=drs2r, ASS ~ VL + COHORT + ART + NT_DIFF_2Rc + BATCHVL + BATCH, family='binomial')
	summary(m2r.1)
	m2r.1.or<- cbind(data.table(COEF=names(coef(m2r.1))), as.data.table( exp(cbind(OR = coef(m2r.1), confint(m2r.1))) ) )
	setnames(m2r.1.or, c('2.5 %','97.5 %'), c('l95','u95'))
	#
	#	batch models without BATCH VL
	#	
	m2f.2	<- glm(data=drs2f, ASS ~ VL + COHORT + ART +  NT_DIFF_2Fc + BATCH, family='binomial')
	summary(m2f.2)
	m2f.2.or<- cbind(data.table(COEF=names(coef(m2f.2))), as.data.table( exp(cbind(OR = coef(m2f.2), confint(m2f.2))) ) )
	setnames(m2f.2.or, c('2.5 %','97.5 %'), c('l95','u95'))		
	m2r.2	<- glm(data=drs2r, ASS ~ VL + COHORT + ART + NT_DIFF_2Rc + BATCH, family='binomial')
	summary(m2r.2)
	m2r.2.or<- cbind(data.table(COEF=names(coef(m2r.2))), as.data.table( exp(cbind(OR = coef(m2r.2), confint(m2r.2))) ) )
	setnames(m2r.2.or, c('2.5 %','97.5 %'), c('l95','u95'))
	#
	#	no BATCHES 
	#	
	m2f.3	<- glm(data=drs2f, ASS ~ VL + COHORT + ART +  NT_DIFF_2Fc, family='binomial')
	summary(m2f.3)
	m2f.3.or<- cbind(data.table(COEF=names(coef(m2f.3))), as.data.table( exp(cbind(OR = coef(m2f.3), confint(m2f.3))) ) )
	setnames(m2f.3.or, c('2.5 %','97.5 %'), c('l95','u95'))		
	m2r.3	<- glm(data=drs2r, ASS ~ VL + COHORT + ART + NT_DIFF_2Rc, family='binomial')
	summary(m2r.3)
	m2r.3.or<- cbind(data.table(COEF=names(coef(m2r.3))), as.data.table( exp(cbind(OR = coef(m2r.3), confint(m2r.3))) ) )
	setnames(m2r.3.or, c('2.5 %','97.5 %'), c('l95','u95'))
	#
	#	no BATCHES and subtype
	#	
	m2f.4	<- glm(data=drs2f, ASS ~ VL + COHORT + ART + ST + NT_DIFF_2Fc, family='binomial')
	summary(m2f.4)
	m2f.4.or<- cbind(data.table(COEF=names(coef(m2f.4))), as.data.table( exp(cbind(OR = coef(m2f.4), confint(m2f.4))) ) )
	setnames(m2f.4.or, c('2.5 %','97.5 %'), c('l95','u95'))		
	m2r.4	<- glm(data=drs2r, ASS ~ VL + COHORT + ART + ST + NT_DIFF_2Rc, family='binomial')
	summary(m2r.4)
	m2r.4.or<- cbind(data.table(COEF=names(coef(m2r.4))), as.data.table( exp(cbind(OR = coef(m2r.4), confint(m2r.4))) ) )
	setnames(m2r.4.or, c('2.5 %','97.5 %'), c('l95','u95'))
	#
	#	NT_DIFF_2Fcat least one mutation 	in models 3,4 and stronger when batches included
	#	NT_DIFF_2Rcat least one mutation	not sig
	#	
	#	NT_DIFF_2FcUnassembled and NT_DIFF_2RcUnassembled always sig
	#
	#	STD	and STB or C					significant for 2F

	#	compare significant batches
	tmp		<- subset(m2f.1.or, u95<0.95)
	tmp[, PR:='2F']
	tmp[, MODEL:= 'adjusting for avg VL in batch']
	tmp2	<- subset(m2r.1.or, u95<0.95)
	tmp2[, PR:='2R']
	tmp2[, MODEL:= 'adjusting for avg VL in batch']
	tmp		<- rbind(tmp, tmp2)	
	tmp2		<- subset(m2f.2.or, u95<0.95)
	tmp2[, PR:='2F']
	tmp2[, MODEL:= 'not adjusting for avg VL in batch']
	tmp		<- rbind(tmp, tmp2)	
	tmp2	<- subset(m2r.2.or, u95<0.95)
	tmp2[, PR:='2R']
	tmp2[, MODEL:= 'not adjusting for avg VL in batch']
	tmp		<- rbind(tmp, tmp2)
	tmp		<- subset(tmp, grepl('BATCH', COEF) & COEF!='BATCHNo matched ID')
	set(tmp, NULL, 'BATCH', tmp[,gsub('BATCH','',COEF)])
	tmp		<- merge(tmp, unique(subset(dr, select=c(BATCH, COHORT))), by='BATCH')
	#
	ggplot(subset(tmp, MODEL=='adjusting for avg VL in batch' & COHORT=='RCCS'), aes(x=BATCH, y=OR, ymin=l95, ymax=u95)) + geom_point() + geom_errorbar() + facet_grid(~PR) + coord_flip()
	ggsave(file=file.path(wdir, 'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_batch_oddsratio_BATCHVL.pdf'), w=7, h=10)
	ggplot(subset(tmp, MODEL=='not adjusting for avg VL in batch' & COHORT=='RCCS'), aes(x=BATCH, y=OR, ymin=l95, ymax=u95)) + geom_point() + geom_errorbar() + facet_grid(~PR) + coord_flip()
	ggsave(file=file.path(wdir, 'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_batch_oddsratio_NOBATCHVL.pdf'), w=7, h=10)
	#
	subset(m2r.4.or, u95<0.95)		
	subset(m2r.1.or, u95<0.95)
	
	#
	#	batch models without BATCH VL and ignoring AC resistance
	#	
	m2f.2b	<- glm(data=subset(drs2f, COHORT!='AC_Resistance'), ASS ~ VL + COHORT + ART +  NT_DIFF_2Fc + BATCH, family='binomial')
	summary(m2f.2b)
	m2f.2b.or<- cbind(data.table(COEF=names(coef(m2f.2b))), as.data.table( exp(cbind(OR = coef(m2f.2b), confint(m2f.2b))) ) )
	setnames(m2f.2b.or, c('2.5 %','97.5 %'), c('l95','u95'))		
	m2r.2b	<- glm(data=subset(drs2r, COHORT!='AC_Resistance'), ASS ~ VL + COHORT + ART + NT_DIFF_2Rc + BATCH, family='binomial')
	summary(m2r.2b)
	m2r.2b.or<- cbind(data.table(COEF=names(coef(m2r.2b))), as.data.table( exp(cbind(OR = coef(m2r.2b), confint(m2r.2b))) ) )
	setnames(m2r.2b.or, c('2.5 %','97.5 %'), c('l95','u95'))
	

	#
	#	Rakai model with communities 
	#	
	m2f.4 <- glm(data=drs2f, ASS ~ VL + ART + LOC + COMM_NUM + NT_DIFF_2Fc, family='binomial')
	summary(m2f.4)
	m2f.4.or<- cbind(data.table(COEF=names(coef(m2f.4))), as.data.table( exp(cbind(OR = coef(m2f.4), confint(m2f.4))) ) )
	setnames(m2f.4.or, c('2.5 %','97.5 %'), c('l95','u95'))	
	subset(m2f.4.or, u95<0.95)	
	m2r.4 <- glm(data=drs2r, ASS ~ VL + ART + LOC + COMM_NUM + NT_DIFF_2Rc, family='binomial')
	summary(m2r.4)
	m2r.4.or<- cbind(data.table(COEF=names(coef(m2r.4))), as.data.table( exp(cbind(OR = coef(m2r.4), confint(m2r.4))) ) )
	setnames(m2r.4.or, c('2.5 %','97.5 %'), c('l95','u95'))	
	subset(m2r.4.or, u95<0.95)		
	

	save(drs2f, drs2r, m2f.1, m2r.1, m2f.1.or, m2r.1.or, m2f.2, m2r.2, m2f.2.or, m2r.2.or, m2f.3, m2r.3, m2f.3.or, m2r.3.or, m2f.4, m2r.4, m2f.4.or, m2r.4.or, file=file.path(wdir, gsub('.rda','_logistic.rda',wfile)))
}
##--------------------------------------------------------------------------------------------------------
##	olli 25.07.16
##--------------------------------------------------------------------------------------------------------
treecomparison.explaingaps.plots.160725<- function()
{
	require(ape)
	require(scales)
	require(data.table)
	require(Hmisc)
	
	wdir			<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/explaingaps'
	wfile			<- 'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.rda'
	#
	#	deal with repeats in global alignment
	#
	load(file.path(wdir, wfile))
	#
	#	re-plot alignment with primers highlighted and all gap columns included too
	#
	min.coverage	<- 600
	min.depth		<- 10	
	#	convert into chunks
	ch				<- lapply(seq_len(nrow(sq)), function(i)
			{
				z	<- gregexpr('1+', paste(as.numeric( !as.character( sq[i,] )%in%c('-','?','n') ), collapse='') )[[1]]
				data.table(PANGEA_ID= rownames(sq)[i], POS=as.integer(z), DEPTH=min.depth, REP=attr(z,"match.length"))
			})
	ch				<- do.call('rbind',ch)	
	#	define SITE
	ch[, SITE:=NA_character_]
	tmp				<- ch[, which(grepl('^R[0-9]+_',PANGEA_ID))]
	set(ch, tmp, 'SITE', 'ZA')
	tmp				<- ch[, which(is.na(SITE) & grepl('PG[0-9]+-[A-Z]+',PANGEA_ID))]	
	set(ch, tmp, 'SITE', ch[tmp, regmatches(PANGEA_ID, regexpr('PG[0-9]+-[A-Z]+',PANGEA_ID))])
	set(ch, NULL, 'SITE', ch[, gsub('PG[0-9]+-','',SITE)])
	ch				<- subset(ch, !is.na(SITE))	
	ch				<- merge(ch, ch[, list(COV=sum(REP)), by='PANGEA_ID'], by='PANGEA_ID')
	#	select min.coverage, select min.depth
	ch			<- subset(ch, COV>=min.coverage & DEPTH>=min.depth)
	#	define chunks
	ch[, POS_NEXT:= POS+REP]	
	ch		<- ch[, list(SITE=SITE, POS=POS, DEPTH=DEPTH, REP=REP, CHUNK=cumsum(as.numeric(c(TRUE, POS[-1]!=POS_NEXT[-length(POS_NEXT)])))), by='PANGEA_ID']
	ch		<- ch[, list(SITE=SITE[1], POS_CH=min(POS), REP_CH=sum(REP), DEPTH_CH= sum(DEPTH*REP)/sum(REP) ), by=c('PANGEA_ID','CHUNK')]
	ch[, DEPTH_MIN:=min.depth]
	set(ch, NULL, 'SITE', ch[, factor(SITE, levels=c('BW', 'ZA', 'UG'), labels=c('Botswana', 'South Africa', 'Uganda'))])
	ch				<- merge(ch, ch[, list(COV=sum(REP_CH)), by='PANGEA_ID'], by='PANGEA_ID')
	ch[, COVP:= COV/ncol(sq)]	
	#	
	require(ggplot2)
	require(viridis)
	ch		<- merge(ch, ch[, list(POS_CHF=min(POS_CH)), by='PANGEA_ID'], by='PANGEA_ID')	
	setkey(ch, SITE, PANGEA_ID)	
	tmp		<- unique(ch)
	setkey(tmp, POS_CHF)	
	tmp[, PLOT:=ceiling(seq_len(nrow(tmp))/1070)]
	ch		<- merge(ch, subset(tmp, select=c(PANGEA_ID, PLOT)), by='PANGEA_ID')
	set(ch, NULL, 'PLOT', ch[, factor(PLOT, levels=c(4,3,2,1), labels=c(4,3,2,1))])
	setkey(ch, POS_CH, SITE)
	set(ch, NULL, 'PANGEA_ID', ch[, factor(PANGEA_ID, levels=unique(PANGEA_ID), labels=unique(PANGEA_ID))])	
	dpani	<- subset(dpan, !is.na(START))[, list(START=START[1], END=START[1]+max(IDX)-1L), by='PR']
	
	ggplot(ch) +
			geom_segment(aes(y=PANGEA_ID, yend=PANGEA_ID, x=POS_CH, xend=POS_CH+REP_CH-1L, colour=SITE)) +  
			geom_rect(data=dpani, aes(xmin=START, xmax=END, ymin=-Inf, ymax=Inf), fill="#3690C0") +
			geom_text(data=dpani, aes(x=START, y=seq_len(length(START))*10+10, label=PR), colour="#3690C0", hjust=-.2, size=2) +
			facet_wrap(~PLOT, scales='free_y', ncol=4) +
			scale_x_continuous(expand=c(0,0), breaks=seq(0,10e3,1e3), minor_breaks=seq(0,10e3,100)) +
			scale_colour_manual(values=c('Botswana'="#1B0C42FF", 'South Africa'="#CF4446FF", 'Uganda'="#781C6DFF")) +						
			labs(x='\nalignment position', y='PANGEA-HIV sequences\n', colour='sampling\nlocation') +
			theme_bw() +
			theme(	axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank(), legend.position='bottom',
					strip.text= element_blank(), strip.background=element_blank()) +
			guides(colour=guide_legend(override.aes=list(size=5)))	
	ggsave(file=file.path(wdir,gsub('.rda','_gapsprimers.pdf',wfile)), w=15, h=15, limitsize = FALSE)
	#
	#	plot Rakai samples by REGA subtype
	#
	setnames(ch, 'PANGEA_ID', 'TAXA')
	tmp		<- unique(subset(dpand, !is.na(PANGEA_ID), select=c(PANGEA_ID, TAXA, STUDY_ID, REGA_GAG_A, REGA_GAG_AS, REGA_GAG_PURE, REGA_GAG_PURES)))
	chr		<- merge(ch, tmp, by='TAXA', all.x=1)
	set(chr, chr[, which(is.na(REGA_GAG_A))], 'REGA_GAG_A', 'No matched ID')	
	chr		<- subset(chr, !is.na(STUDY_ID))
	#	redefine ordering
	chr[, PLOT:=NULL]	
	setkey(chr, REGA_GAG_A, TAXA)	
	tmp		<- unique(chr)
	setkey(tmp, REGA_GAG_A, COVP, TAXA)	
	tmp		<- tmp[, list(TAXA=TAXA, PLOT=REGA_GAG_A, PLOT_ID=seq_along(TAXA)), by='REGA_GAG_A']
	chr		<- merge(chr, subset(tmp, select=c(TAXA, PLOT_ID)), by='TAXA')	
	set(chr, NULL, 'REGA_GAG_A', chr[, factor(REGA_GAG_A, levels=c("A (A1)", "D", "C", "Check bootscan", "D (10_CD)", "CRF 21_A2D", "G", "Sequence error", "Unassigned"))])
	
	ggplot(subset(chr, REGA_GAG_A%in%c("A (A1)", "D", "C", "Check bootscan"))) +
			geom_segment(aes(y=PLOT_ID, yend=PLOT_ID, x=POS_CH, xend=POS_CH+REP_CH-1L, colour=REGA_GAG_A)) +  
			geom_rect(data=dpani, aes(xmin=START, xmax=END, ymin=-Inf, ymax=Inf), fill="black") +			
			facet_wrap(~REGA_GAG_A, scales='free_y', ncol=4) +
			scale_x_continuous(expand=c(0,0), breaks=dpani$START, labels=dpani$PR) +
			scale_y_continuous(expand=c(0,0)) +
			scale_colour_brewer(palette='Set1') +						
			labs(x='\nalignment position', y='Rakai PANGEA-HIV sequences\n', colour='Rega subtype assignment on gag') +
			theme_bw() +
			theme(	legend.position='bottom', strip.text= element_blank(), strip.background=element_blank()) +
			guides(colour=guide_legend(override.aes=list(size=5)))
	ggsave(file=file.path(wdir,gsub('.rda','_gapsprimers_REGAsubtypes.pdf',wfile)), w=15, h=15, limitsize = FALSE)
	#
	#	plot Rakai samples by COMET subtype
	#
	tmp		<- unique(subset(dpand, !is.na(PANGEA_ID), select=c(PANGEA_ID, TAXA, STUDY_ID, COMET_Region1, COMET_Region2, COMET_Region3)))
	chr		<- merge(ch, tmp, by='TAXA', all.x=1)
	chr		<- subset(chr, !is.na(STUDY_ID))
	set(chr, chr[, which(is.na(COMET_Region1))], 'COMET_Region1', 'unassigned')
	set(chr, chr[, which(is.na(COMET_Region2))], 'COMET_Region2', 'unassigned')
	set(chr, chr[, which(is.na(COMET_Region3))], 'COMET_Region3', 'unassigned')
	chr[, COMET_Region123:= paste(COMET_Region1,COMET_Region2,COMET_Region3,sep='-')]
	set(chr, chr[, which(!COMET_Region123%in%c('A1-A1-A1','A2-A2-A2','B-B-B','C-C-C','D-D-D','other-other-other','unassigned-unassigned-unassigned'))], 'COMET_Region123', 'mixed')
	#	Region123 mixes cause and effect: must have good representation of all regions in order to call subtypes well
	chr[, COMET_Region13:= paste(COMET_Region1,COMET_Region3,sep='-')]
	set(chr, chr[, which(!COMET_Region13%in%c('A1-A1','A2-A2','B-B','C-C','D-D','other-other','unassigned-unassigned'))], 'COMET_Region13', 'mixed')
	#	not sure if this should be used..
	#	redefine ordering
	chr[, PLOT:=NULL]	
	setkey(chr, COMET_Region1, TAXA)	
	tmp		<- unique(chr)
	setkey(tmp, COMET_Region1, COVP, TAXA)	
	tmp		<- tmp[, list(TAXA=TAXA, PLOT=COMET_Region1, PLOT_ID=seq_along(TAXA)), by='COMET_Region1']
	chr		<- merge(chr, subset(tmp, select=c(TAXA, PLOT_ID)), by='TAXA')		
	
	ggplot(subset(chr, !COMET_Region1%in%c('check','other'))) +
			geom_segment(aes(y=PLOT_ID, yend=PLOT_ID, x=POS_CH, xend=POS_CH+REP_CH-1L, colour=COMET_Region1)) +  
			geom_rect(data=dpani, aes(xmin=START, xmax=END, ymin=-Inf, ymax=Inf), fill="black") +
			geom_vline(xintercept=c(2200,3000)) +
			facet_wrap(~COMET_Region1, scales='free_y', ncol=4) +
			scale_x_continuous(expand=c(0,0), breaks=dpani$START, labels=dpani$PR) +
			scale_y_continuous(expand=c(0,0)) +
			scale_colour_brewer(palette='Set1') +						
			labs(x='\nalignment position', y='Rakai PANGEA-HIV sequences\n', colour='COMET subtype assignment\non region 1') +
			theme_bw() +
			theme(	legend.position='bottom', strip.text= element_blank(), strip.background=element_blank()) +
			guides(colour=guide_legend(override.aes=list(size=5)))
	ggsave(file=file.path(wdir,gsub('.rda','_gapsprimers_COMETsubtypes.pdf',wfile)), w=15, h=10, limitsize = FALSE)
	#
	#	plot Rakai samples by Sanger processing batch
	#	
	tmp		<- unique(subset(dpand, !is.na(PANGEA_ID), select=c(PANGEA_ID, TAXA, STUDY_ID, SANGER_ID)))
	chr		<- merge(ch, tmp, by='TAXA', all.x=1)
	#chr		<- subset(chr, !is.na(STUDY_ID))
	chr[, BATCH:=NA_character_]
	tmp		<- chr[, which(!is.na(SANGER_ID))]
	set(chr, tmp, 'BATCH', chr[tmp, regmatches(SANGER_ID,regexpr('^[0-9]+', SANGER_ID))])	
	set(chr, chr[, which(is.na(BATCH))], 'BATCH', 'No matched ID')
	#	redefine ordering
	chr[, PLOT:=NULL]	
	setkey(chr, BATCH, TAXA)	
	tmp		<- unique(chr)
	setkey(tmp, BATCH, COVP, TAXA)	
	tmp		<- tmp[, list(TAXA=TAXA, PLOT=BATCH, PLOT_ID=seq_along(TAXA)), by='BATCH']
	chr		<- merge(chr, subset(tmp, select=c(TAXA, PLOT_ID)), by='TAXA')	
	ggplot(chr) +
			geom_segment(aes(y=PLOT_ID, yend=PLOT_ID, x=POS_CH, xend=POS_CH+REP_CH-1L, colour=factor(is.na(STUDY_ID), levels=c(TRUE,FALSE), labels=c('OTHER','Rakai')))) +  
			geom_rect(data=dpani, aes(xmin=START, xmax=END, ymin=-Inf, ymax=Inf), fill="black") +			
			facet_wrap(~BATCH, scales='free_y', ncol=4) +
			scale_x_continuous(expand=c(0,0), breaks=dpani$START, labels=dpani$PR) +
			scale_y_continuous(expand=c(0,0)) +
			#scale_colour_brewer(palette='Set1') +						
			labs(x='\nalignment position', y='Rakai PANGEA-HIV sequences\n', colour='Cohort site') +
			theme_bw() +
			theme(	legend.position='bottom') +
			guides(colour=guide_legend(override.aes=list(size=5)))
	ggsave(file=file.path(wdir,gsub('.rda','_gapsprimers_Batch.pdf',wfile)), w=15, h=100, limitsize = FALSE)	
	#
	#	plot Rakai samples by ART
	#
	#	redefine ordering	
	tmp		<- unique(subset(dpand, !is.na(PANGEA_ID), select=c(PANGEA_ID, TAXA, STUDY_ID, SAMPLEDATE, ARTSTART, selfReportArt, everSelfReportArt, RECENTVL)))
	tmp[, ART:= as.integer(ARTSTART<SAMPLEDATE)]
	set(tmp, tmp[, which(is.na(ART))], 'ART', 0L)
	set(tmp, tmp[, which(ART==0 & everSelfReportArt==1)], 'ART', 2L)
	set(tmp, tmp[, which(ART==0 & RECENTVL<1e4)], 'ART', 3L)
	set(tmp, NULL, 'ART', tmp[, factor(ART, levels=c(0L,1L,2L,3L), labels=c('no ART', 'ART started', 'ART self reported','no ART but VL<1e4'))])
	
	chr		<- merge(ch, tmp, by='TAXA', all.x=1)
	set(chr, chr[, which(is.na(ART))], 'ART', 'No matched ID')	
	chr		<- subset(chr, !is.na(STUDY_ID))
	#	redefine ordering
	chr[, PLOT:=NULL]	
	setkey(chr, TAXA)	
	tmp		<- unique(chr)
	setkey(tmp, ART, COVP, TAXA)	
	tmp		<- tmp[, list(TAXA=TAXA, PLOT=ART, PLOT_ID=seq_along(TAXA)), by='ART']
	chr		<- merge(chr, subset(tmp, select=c(TAXA, PLOT_ID)), by='TAXA')	
	#set(chr, NULL, 'REGA_GAG_A', chr[, factor(REGA_GAG_A, levels=c("A (A1)", "D", "C", "Check bootscan", "D (10_CD)", "CRF 21_A2D", "G", "Sequence error", "Unassigned"))])	
	ggplot(chr) +
			geom_segment(aes(y=PLOT_ID, yend=PLOT_ID, x=POS_CH, xend=POS_CH+REP_CH-1L, colour=ART)) +  
			geom_rect(data=dpani, aes(xmin=START, xmax=END, ymin=-Inf, ymax=Inf), fill="black") +			
			facet_wrap(~ART, scales='free_y', ncol=4) +
			scale_x_continuous(expand=c(0,0), breaks=dpani$START, labels=dpani$PR) +
			scale_y_continuous(expand=c(0,0)) +
			scale_colour_brewer(palette='Set2') +						
			labs(x='\nalignment position', y='Rakai PANGEA-HIV sequences\n', colour='ART status') +
			theme_bw() +
			theme(	legend.position='bottom', strip.text= element_blank(), strip.background=element_blank()) +
			guides(colour=guide_legend(override.aes=list(size=5)))
	ggsave(file=file.path(wdir,gsub('.rda','_gapsprimers_ARTstatus.pdf',wfile)), w=15, h=15, limitsize = FALSE)		
	#
	#	plot Rakai samples by ART
	#
	#	redefine ordering	
	tmp		<- unique(subset(dpand, !is.na(PANGEA_ID), select=c(PANGEA_ID, TAXA, STUDY_ID, SAMPLEDATE, ARTSTART, selfReportArt, everSelfReportArt, RECENTVL)))
	tmp[, VL:= cut(RECENTVL, breaks=c(0, 1e4, 2e4, 4e4, 1e5, Inf), labels=c('<1e4','1e4-2e4','2e4-4e4','4e4-1e5','>1e5'))]	
	set(tmp, tmp[, which(is.na(VL))], 'VL', 'No VL measured')
	chr		<- merge(ch, tmp, by='TAXA', all.x=1)
	set(chr, chr[, which(is.na(VL))], 'VL', 'No matched ID')	
	chr		<- subset(chr, !is.na(STUDY_ID))	
	chr[, PLOT:=NULL]	
	setkey(chr, TAXA)	
	tmp		<- unique(chr)
	setkey(tmp, VL, COVP, TAXA)	
	tmp		<- tmp[, list(TAXA=TAXA, PLOT=VL, PLOT_ID=seq_along(TAXA)), by='VL']
	chr		<- merge(chr, subset(tmp, select=c(TAXA, PLOT_ID)), by='TAXA')	
	ggplot(chr) +
			geom_segment(aes(y=PLOT_ID, yend=PLOT_ID, x=POS_CH, xend=POS_CH+REP_CH-1L, colour=VL)) +  
			geom_rect(data=dpani, aes(xmin=START, xmax=END, ymin=-Inf, ymax=Inf), fill="black") +			
			facet_wrap(~VL, scales='free_y', ncol=6) +
			scale_x_continuous(expand=c(0,0), breaks=dpani$START, labels=dpani$PR) +
			scale_y_continuous(expand=c(0,0)) +
			scale_colour_brewer(palette='Dark2') +						
			labs(x='\nalignment position', y='Rakai PANGEA-HIV sequences\n', colour='Viral load status') +
			theme_bw() +
			theme(	legend.position='bottom', strip.text= element_blank(), strip.background=element_blank()) +
			guides(colour=guide_legend(override.aes=list(size=5)))
	ggsave(file=file.path(wdir,gsub('.rda','_gapsprimers_VLstatus.pdf',wfile)), h=7, w=20, limitsize = FALSE)		
	#
	#	plot Rakai samples by region
	#
	#	redefine ordering	
	tmp		<- unique(subset(dpand, !is.na(PANGEA_ID), select=c(PANGEA_ID, TAXA, STUDY_ID, LOC)))
	set(tmp, NULL, 'LOC', tmp[, factor(LOC)])
	chr		<- merge(ch, tmp, by='TAXA', all.x=1)
	chr		<- subset(chr, !is.na(STUDY_ID))	
	chr[, PLOT:=NULL]	
	setkey(chr, TAXA)	
	tmp		<- unique(chr)
	setkey(tmp, LOC, COVP, TAXA)	
	tmp		<- tmp[, list(TAXA=TAXA, PLOT=LOC, PLOT_ID=seq_along(TAXA)), by='LOC']
	chr		<- merge(chr, subset(tmp, select=c(TAXA, PLOT_ID)), by='TAXA')	
	ggplot(chr) +
			geom_segment(aes(y=PLOT_ID, yend=PLOT_ID, x=POS_CH, xend=POS_CH+REP_CH-1L, colour=LOC)) +  
			geom_rect(data=dpani, aes(xmin=START, xmax=END, ymin=-Inf, ymax=Inf), fill="black") +			
			facet_wrap(~LOC, scales='free_y', ncol=6) +
			scale_x_continuous(expand=c(0,0), breaks=dpani$START, labels=dpani$PR) +
			scale_y_continuous(expand=c(0,0)) +
			scale_colour_brewer(palette='Dark2') +						
			labs(x='\nalignment position', y='Rakai PANGEA-HIV sequences\n', colour='region') +
			theme_bw() +
			theme(	legend.position='bottom', strip.text= element_blank(), strip.background=element_blank()) +
			guides(colour=guide_legend(override.aes=list(size=5)))
	ggsave(file=file.path(wdir,gsub('.rda','_gapsprimers_REGION.pdf',wfile)), h=7, w=20, limitsize = FALSE)	
}
##--------------------------------------------------------------------------------------------------------
##	olli 25.07.16
##--------------------------------------------------------------------------------------------------------
treecomparison.explaingaps.countgaps<- function()
{
	require(ape)
	require(data.table)
	require(big.phylo)
	wdir			<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/explaingaps'
	wfile			<- 'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.rda'
	#
	load(file.path(wdir,wfile))
	#
	#	get summary
	#
	tmp		<- merge(dgd, subset(dm, select=c(TAXA,COHORT)),by='TAXA')
	tmp[, LEN:=END-START+1L]
	#	add the odd RCCS sequence run
	tmp		<- merge(tmp, subset(dm, select=c(TAXA, SANGER_ID)), by='TAXA')
	tmp[, SANGER_PLATE:= as.integer(regmatches(SANGER_ID, regexpr('^[0-9]+', SANGER_ID)))]
	tmp2	<- subset(tmp, SANGER_PLATE>=15892 & SANGER_PLATE<=15964)
	tmp2[, COHORT:='RCCS_run_odd_plates']	
	tmp		<- rbind(tmp, tmp2)
	tmp2	<- subset(tmp, SANGER_PLATE<15892 | SANGER_PLATE>15964)
	tmp2[, COHORT:='RCCS_other_plates']
	tmp		<- rbind(tmp, tmp2)
	#	add all Uganda sequences
	tmp2	<- subset(tmp, COHORT=='UG-MRC'|COHORT=='RCCS')
	tmp2[, COHORT:='Uganda']
	tmp		<- rbind(tmp, tmp2)
	#	look at distribution of gaps by cohort
	dgi		<- subset(tmp, GENE=='GAG+POL+ENV')[, list(P=ecdf(1-UNASS)(seq(0,1,.01)), PC_ASS=seq(0,1,.01)), by=c('COHORT','GENE','LEN')]
	#	                RCCS GAG+POL+ENV 8926 0.2392074    0.2
	dgi		<- subset(tmp, GENE=='GAG+POL+ENV' & COHORT%in%c('BW-Mochudi','UG-MRC','RCCS','AC_Resistance'))[, list(P=ecdf(1-UNASS)(seq(0,1,.01)), PC_ASS=seq(0,1,.01))]
	#					1: 0.206    0.2
	#
	#	make histogram
	tmp2	<- subset(tmp, GENE=='GAG+POL+ENV' & COHORT%in%c('BW-Mochudi','UG-MRC','RCCS','AC_Resistance'))
	set(tmp2, NULL, 'COHORT', tmp2[, factor(COHORT, levels=c('BW-Mochudi','UG-MRC','RCCS','AC_Resistance'),
													labels=c('Mochudi Prevention Project', 'Uganda-MRC', 'Rakai Community Cohort Study', 'South Africa Resistance Cohort'))])				
	tmp2	<- tmp2[, list(P_UNASS=seq(0,1,0.01), CUM= length(UNASS)*ecdf(UNASS)(seq(0,1,0.01))), by='COHORT']
	ggplot(tmp2, aes(x=P_UNASS, fill=COHORT, y=CUM)) + 
			geom_bar(stat='identity', position='stack') +
			scale_x_continuous(labels=percent, expand=c(0,0), limit=c(-0.01,1.01), breaks=seq(0, 1, 0.2)) +
			scale_y_continuous(expand=c(0,0), limits=c(0,4500))+
			scale_fill_manual(values=c('Mochudi Prevention Project'="#33638DFF", 'South Africa Resistance Cohort'="#CF4446FF", 'Rakai Community Cohort Study'="#A8327DFF", 'Uganda-MRC'="#29AF7FFF")) +
			theme_bw() + theme(legend.position='bottom') +
			guides(fill=guide_legend(ncol=2)) + 
			labs(	x='\nx% = proportion of unassembled sites per sequence is at most x%',
					y='PANGEA-HIV sequences (cumulated)\n',
					fill='cohort site')
	ggsave(file.path('~/Dropbox (Infectious Disease)/2016_PANGEA_treecomp/figures','PANGEA_cumulated_by_gaps.pdf'), w=7, h=7)
	#	
	dgi		<- tmp[, list(ASS= c(mean(1-UNASS), quantile(1-UNASS, prob=c(0.25,.5,.75))), STAT=c('mean',paste('qu',c(0.25,.5,.75),sep='')) ), by=c('COHORT','GENE','LEN')]
	dgi[, P:= 100*round(ASS, d=2)]
	dgi[, N:= round(ASS*LEN,d=0)]
	dgi[, LABEL:= paste(N,'/',LEN,' (',P,'%)',sep='')]
	#dcast.data.table(subset(dgi, STAT=='mean'), GENE~COHORT, value.var='LABEL')
	dgi		<- dcast.data.table(subset(dgi, STAT=='mean'), GENE+LEN~COHORT, value.var='P')		
	#
	#	save
	#
	write.csv(dgi, row.names=FALSE, file=file.path(wdir,'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_countgaps.csv'))
	save(dgi, file=file.path(wdir,'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_countgaps.rda'))
}
##--------------------------------------------------------------------------------------------------------
##	olli 25.07.16
##--------------------------------------------------------------------------------------------------------
treecomparison.explaingaps.collect.data<- function()
{
	require(ape)
	require(data.table)
	wdir			<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/explaingaps'
	wfile			<- 'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.rda'
	
	#
	#	deal with repeats in global alignment
	#
	if(0)
	{
		sq				<- read.dna("~/Dropbox (Infectious Disease)/PANGEA_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.fasta", format='fasta')		
		sqi				<- data.table(TAXA=rownames(sq), DUMMY=seq_len(nrow(sq)))
		tmp				<- sqi[, which(duplicated(TAXA))]
		set(sqi, tmp, 'TAXA', sqi[tmp, paste(TAXA,'-R2',sep='')])
		setkey(sqi, DUMMY)
		rownames(sq)	<- sqi[,TAXA]	
		tmp				<- sapply(seq_len(nrow(sq)), function(i) base.freq(sq[i,], all=TRUE, freq=TRUE))
		sqi[, COV:=ncol(sq)-apply( tmp[c('-','?'),], 2, sum	)]	
		sqi[, PNG:= sqi[, factor(grepl('PG',TAXA),levels=c(TRUE,FALSE),labels=c('Y','N'))]]		
		sqi[, SITE:= NA_character_]
		tmp				<- sqi[, which(PNG=='Y')]
		set(sqi, tmp, 'SITE', sqi[tmp, substring(sapply(strsplit(TAXA,'-'),'[[',2),1,2)])
		sqi[, PANGEA_ID:= gsub('-R[0-9]+','',TAXA)]
		#	get gap-column free alignment
		sqp		<- sq[subset(sqi, PNG=='Y' | grepl('HXB2', TAXA))[, TAXA],]
		sqp		<- seq.rmgaps(sqp, rm.only.col.gaps=1, verbose=1)			
		#	save
		write.dna(sq, file=file.path(wdir,'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.fasta'), format='fa')		
		save(sq, sqp, sqi, file=file.path(wdir,wfile))
	}	
	#
	#	find primer coordinates and gene coordinates in PANGEA alignment
	#	
	if(0)
	{		
		load(file.path(wdir,wfile))
		#
		pan1f	<- 'agcc.gggagctctctg'
		pan1r	<- 'cctccaattcc.cctatcatttt'
		pan2f	<- 'gggaagtga.atagc.ggaac'
		pan2r	<- 'ctgccatctgttttccata.tc'	
		pan3f	<- 'ttaaaagaaaaggggggattggg'
		pan3r	<- 'tggc.ytgtaccgtcagcg'
		pan4f	<- 'cctatggcaggaagaagcg'
		pan4r	<- 'ctt.tatgcag..tctgaggg'		#I had to remove one nucleotide here it is: CW ( c. ) instead of (..)
		#	get the forward code of the reverse primers
		dpan	<- rbind( 	data.table(PR='1R', SEQR=strsplit(pan1r,'')[[1]], IDXR=seq_len(nchar(pan1r)), IDX=rev(seq_len(nchar(pan1r))) ),
				data.table(PR='2R', SEQR=strsplit(pan2r,'')[[1]], IDXR=seq_len(nchar(pan2r)), IDX=rev(seq_len(nchar(pan2r)))),
				data.table(PR='3R', SEQR=strsplit(pan3r,'')[[1]], IDXR=seq_len(nchar(pan3r)), IDX=rev(seq_len(nchar(pan3r)))),
				data.table(PR='4R', SEQR=strsplit(pan4r,'')[[1]], IDXR=seq_len(nchar(pan4r)), IDX=rev(seq_len(nchar(pan4r))))	)
		dpan	<- merge(dpan, data.table(SEQR=c('a','t','c','g','.'), SEQ=c('t','a','g','c','.')), by='SEQR')
		setkey(dpan, PR, IDX)		
		#	add the forward primers
		tmp		<- rbind( 	data.table(PR='1F', SEQ=strsplit(pan1f,'')[[1]], IDX=seq_len(nchar(pan1f)) ),
				data.table(PR='2F', SEQ=strsplit(pan2f,'')[[1]], IDX=seq_len(nchar(pan2f)) ),
				data.table(PR='3F', SEQ=strsplit(pan3f,'')[[1]], IDX=seq_len(nchar(pan3f)) ),
				data.table(PR='4F', SEQ=strsplit(pan4f,'')[[1]], IDX=seq_len(nchar(pan4f)) )	)
		dpan	<- rbind( tmp, dpan, use.names=TRUE, fill=TRUE)
		#	check primers exist in HXB2
		load("~/git/hivclust/pkg/data/refseq_hiv1_hxb2.rda")
		hxb2	<- subset(hxb2, !is.na(HXB2.Position))[, paste(as.character(HXB2.K03455),collapse='')]		
		tmp		<- dpan[, list(START=unlist(gregexpr(paste(SEQ, collapse=''),hxb2))), by='PR']
		stopifnot( tmp[, all(START>1)] )
		#	build region table
		tmp		<- dcast.data.table(dpan[, list(SEQ=paste(SEQ,collapse='')), by='PR'], .~PR, value.var='SEQ')		
		dgene	<- rbind( 	data.table(PRIMER='N', GENE='GAG',  LOC='start', SEQ=strsplit("atgggtgcgagagcgtcagtatt",'')[[1]]),
							data.table(PRIMER='N', GENE='GAG',  LOC='end', SEQ=strsplit("aacgacccctcgtcacaataa",'')[[1]]),
							data.table(PRIMER='N', GENE='POL',  LOC='start', SEQ=strsplit("cctcaggtcactctttggca",'')[[1]]),
							data.table(PRIMER='N', GENE='POL',  LOC='end', SEQ=strsplit("gtagacaggatgaggattag",'')[[1]]),
							data.table(PRIMER='N', GENE='ENV',  LOC='start', SEQ=strsplit("atgagagtgaaggagaaatatcag",'')[[1]]),
							data.table(PRIMER='N', GENE='ENV',  LOC='end', SEQ=strsplit("cttggaaaggattttgctataa",'')[[1]]),
							data.table(PRIMER='N', GENE='GAG+POL+ENV',  LOC='start', SEQ=strsplit("atgggtgcgagagcgtcagtatt",'')[[1]]),
							data.table(PRIMER='N', GENE='GAG+POL+ENV',  LOC='end', SEQ=strsplit("cttggaaaggattttgctataa",'')[[1]]),
							data.table(PRIMER='Y', GENE='START-2F',  LOC='start', SEQ=strsplit("atgggtgcgagagcgtcagtatt",'')[[1]]),
							data.table(PRIMER='Y', GENE='START-2F',  LOC='end', SEQ=strsplit(tmp[['2F']],'')[[1]]),							
							data.table(PRIMER='Y', GENE='START-1R',  LOC='start', SEQ=strsplit("atgggtgcgagagcgtcagtatt",'')[[1]]),
							data.table(PRIMER='Y', GENE='START-1R',  LOC='end', SEQ=strsplit(tmp[['1R']],'')[[1]]),														
							data.table(PRIMER='Y', GENE='2F-1R',  LOC='start', SEQ=strsplit(tmp[['2F']],'')[[1]]),
							data.table(PRIMER='Y', GENE='2F-1R',  LOC='end', SEQ=strsplit(tmp[['1R']],'')[[1]]),
							data.table(PRIMER='Y', GENE='1R-3F',  LOC='start', SEQ=strsplit(tmp[['1R']],'')[[1]]),
							data.table(PRIMER='Y', GENE='1R-3F',  LOC='end', SEQ=strsplit(tmp[['3F']],'')[[1]]),
							data.table(PRIMER='Y', GENE='3F-2R',  LOC='start', SEQ=strsplit(tmp[['3F']],'')[[1]]),
							data.table(PRIMER='Y', GENE='3F-2R',  LOC='end', SEQ=strsplit(tmp[['2R']],'')[[1]]),							
							data.table(PRIMER='Y', GENE='3F-4F',  LOC='start', SEQ=strsplit(tmp[['3F']],'')[[1]]),
							data.table(PRIMER='Y', GENE='3F-4F',  LOC='end', SEQ=strsplit(tmp[['4F']],'')[[1]]),														
							data.table(PRIMER='Y', GENE='2R-4F',  LOC='start', SEQ=strsplit(tmp[['2R']],'')[[1]]),
							data.table(PRIMER='Y', GENE='2R-4F',  LOC='end', SEQ=strsplit(tmp[['4F']],'')[[1]]),
							data.table(PRIMER='Y', GENE='4F-3R',  LOC='start', SEQ=strsplit(tmp[['4F']],'')[[1]]),
							data.table(PRIMER='Y', GENE='4F-3R',  LOC='end', SEQ=strsplit(tmp[['3R']],'')[[1]]),
							data.table(PRIMER='Y', GENE='3R-END',  LOC='start', SEQ=strsplit(tmp[['3R']],'')[[1]]),
							data.table(PRIMER='Y', GENE='3R-END',  LOC='end', SEQ=strsplit("cttggaaaggattttgctataa",'')[[1]])						
							)
		#	find positions in alignment
		tmp		<- which(grepl('HXB2',rownames(sqp)))
		sqhxb2	<- paste(as.character(sqp[tmp,]), collapse='')
		tmp		<- dgene[, list(START=unlist(gregexpr(paste(SEQ, collapse='-*'),sqhxb2)), LEN=length(SEQ)), by=c('PRIMER','GENE','LOC')]
		tmp2	<- tmp[, which(grepl('GAG|ENV|POL',GENE) & LOC=='end')]
		set(tmp, tmp2, 'START', tmp[tmp2, START+LEN-1L])
		tmp2	<- tmp[, which(!grepl('GAG|ENV|POL',GENE) & LOC=='end')]
		set(tmp, tmp2, 'START', tmp[tmp2, START-1L])
		tmp2	<- tmp[, which(!grepl('GAG|ENV|POL',GENE) & LOC=='start')]
		set(tmp, tmp2, 'START', tmp[tmp2, START+LEN])
		dgene	<- dcast.data.table(tmp, PRIMER+GENE~LOC, value.var='START')	
		setnames(dgene, colnames(dgene), toupper(colnames(dgene)))
		#	add half regions
		tmp		<- subset(dgene, PRIMER=='Y')
		set(tmp, NULL, 'GENE', tmp[, paste(GENE,'-firsthalf',sep='')])
		set(tmp, NULL, 'END', tmp[, START+(END-START)/2])
		dgene	<- rbind(dgene, tmp)
		set(tmp, NULL, 'GENE', tmp[, gsub('-firsthalf','-secondhalf',GENE)])
		set(tmp, NULL, 'END', tmp[, START+(END-START)*2])
		set(tmp, NULL, 'START', tmp[, START+(END-START)/2])
		dgene	<- rbind(dgene, tmp)
		set(dgene, NULL, 'START', dgene[,round(START,d=0)])
		set(dgene, NULL, 'END', dgene[,round(END,d=0)])
		#
		dgene[, LEN:=END-START+1L]
		set(dgene, NULL, 'GENE', dgene[, factor(GENE, levels=c(	"GAG+POL+ENV","GAG","POL","ENV", 
																"START-2F","START-1R","2F-1R","1R-3F","3F-2R","2R-4F","3F-4F","4F-3R","3R-END",
																"START-2F-firsthalf","START-1R-firsthalf","2F-1R-firsthalf","1R-3F-firsthalf","3F-2R-firsthalf","3F-4F-firsthalf","2R-4F-firsthalf","4F-3R-firsthalf","3R-END-firsthalf",
																"START-2F-secondhalf","START-1R-secondhalf","2F-1R-secondhalf","1R-3F-secondhalf","3F-2R-secondhalf","3F-4F-secondhalf","2R-4F-secondhalf","4F-3R-secondhalf","3R-END-secondhalf"))])
		setkey(dgene, GENE)
		set(dgene, NULL, 'PRIMER', NULL)
		#
		tmp		<- dpan[, list(START=unlist(gregexpr(paste(SEQ, collapse='-*'),sqhxb2))), by=c('PR')]
		dpan	<- merge(dpan, subset(tmp, START>0), by='PR',all.x=1)
		#
		#	add non-ambiguous primers for 1R 2F 2R 3F (none) 4F (none)
		#	pan1r	<- 'cctccaattccYcctatcatttt'. 'y = c or t'. So in forward sense: G or A  
		#	pan2f	<- 'gggaagtgaYatagcWggaac'. 'W= a or t'. 
		#	pan2r	<- 'ctgccatctgttttccataRtc'. 'R= a or g'. So in forward sense: T or C
		tmp		<- subset(dpan, PR=='1R')
		tmp[, PR:='1Ra']
		set(tmp, 12L, 'SEQR', 'c')
		set(tmp, 12L, 'SEQ', 'g')
		dpan	<- rbind(dpan, tmp)
		tmp[, PR:='1Rb']
		set(tmp, 12L, 'SEQR', 't')
		set(tmp, 12L, 'SEQ', 'a')
		dpan	<- rbind(dpan, tmp)
		tmp		<- subset(dpan, PR=='2F')
		tmp[, PR:='2Fa']
		set(tmp, 10L, 'SEQ', 'c')
		set(tmp, 16L, 'SEQ', 'a')
		dpan	<- rbind(dpan, tmp)
		tmp[, PR:='2Fb']
		set(tmp, 10L, 'SEQ', 't')
		set(tmp, 16L, 'SEQ', 't')
		dpan	<- rbind(dpan, tmp)		
		tmp		<- subset(dpan, PR=='2R')
		tmp[, PR:='2Ra']
		set(tmp, 3L, 'SEQR', 'a')
		set(tmp, 3L, 'SEQ', 't')
		dpan	<- rbind(dpan, tmp)
		tmp		<- subset(dpan, PR=='2R')
		tmp[, PR:='2Rb']
		set(tmp, 3L, 'SEQR', 'g')
		set(tmp, 3L, 'SEQ', 'c')
		dpan	<- rbind(dpan, tmp)		
		#
		save(sq, sqp, sqi, dgene, dpan, file=file.path(wdir,wfile))
		
		#	extract START to 1R alignment
		write.dna(sqp[, 1:(1893+23)], file=file.path(wdir, paste(gsub('.rda','',wfile),'_region_1F1R.fasta', sep='')), format='fa', colsep='', nbcol=-1)		
		write.dna(sqp[, 4570:5573], file=file.path(wdir, paste(gsub('.rda','',wfile),'_region_3F4F.fasta', sep='')), format='fa', colsep='', nbcol=-1)		
		write.dna(sqp[, 5555:8903], file=file.path(wdir, paste(gsub('.rda','',wfile),'_region_4F3R.fasta', sep='')), format='fa', colsep='', nbcol=-1)		
		write.dna(sqp[ subset(dgd, GENE=='GAG+POL+ENV' & grepl('UG',TAXA) & LEN*(1-UNASS)>8000)[, TAXA], ], file=file.path(wdir, paste(gsub('.rda','',wfile),'_UGfull.fasta', sep='')), format='fa', colsep='', nbcol=-1)
		write.dna(sqp[ subset(dgd, GENE=='GAG+POL+ENV' & grepl('ZA|BW',TAXA) & LEN*(1-UNASS)>8000)[, TAXA], ], file=file.path(wdir, paste(gsub('.rda','',wfile),'_ZABWfull.fasta', sep='')), format='fa', colsep='', nbcol=-1)		
		write.dna(sqp[, 5555:8903], file=file.path(wdir, paste(gsub('.rda','',wfile),'_region_4F3R.fasta', sep='')), format='fa', colsep='', nbcol=-1)
		#	extract primer alignments and write to file
		subset(dpan, !is.na(START))[, {
					write.dna(sqp[, seq.int(START[1], len=length(START))], file=file.path(wdir, paste(gsub('.rda','',wfile),'_primer_',PR[1],'_start_',START[1],'.fasta', sep='')), format='fa', colsep='', nbcol=-1)
				}, by='PR']	
		#	extract primer alignments and write to file +- 75 bp
		subset(dpan, !is.na(START))[, {
					write.dna(sqp[, seq.int(START[1]-75L, len=length(START)+75L)], file=file.path(wdir, paste(gsub('.rda','',wfile),'_primerplusminus75bp_',PR[1],'_start_',START[1],'.fasta', sep='')), format='fa', colsep='', nbcol=-1)
				}, by='PR']		
		
		
	}
	#
	#	calculate number of mutations in primers and gaps in gene regions
	#
	if(0)
	{
		load(file.path(wdir,wfile))
		#	get differences relative to HXB2 on primer by primer position (IDX)
		dpand		<- subset(dpan, !is.na(START))[, {
					#START	<- rep(1894, 23); z<- '1Ra'
					z		<- PR					
					psq		<- as.character( sqp[!grepl('HXB2',rownames(sqp)), seq.int(START[1], len=length(START))] )
					tmp		<- subset(dpan, PR==z)[, gsub('\\.','n',paste(SEQ, collapse=''))]
					tmp		<- unlist(strsplit(tmp,''))
					z		<- which(tmp=='n')
					if(length(z))
					{
						psq	<- psq[, -z]
						tmp	<- tmp[-z]						
					}
					tmp		<- as.data.table( melt( t( t(psq)==tmp ), varnames=c('TAXA','POS') ) )
					#sqhxb2i	<- which(grepl('HXB2',rownames(psq)))				
					#tmp		<- as.data.table( melt( t( t(psq)==psq[sqhxb2i,] ) ) )
					setnames(tmp, 'value', 'NT_DIFF')
					tmp2	<- as.data.table( melt( psq=='?', varnames=c('TAXA','POS') ) )
					setnames(tmp2, 'value', 'MISS')
					tmp		<- merge(tmp, tmp2, by=c('TAXA','POS'))				
					tmp2	<- as.data.table( melt( psq=='n', varnames=c('TAXA','POS') ) )
					setnames(tmp2, 'value', 'ANY')
					tmp		<- merge(tmp, tmp2, by=c('TAXA','POS'))									
					set(tmp, NULL, 'NT_DIFF', tmp[, as.integer(!NT_DIFF)])
					#	adjust primer position according to z
					while(length(z))
					{
						tmp2<- tmp[, which(POS>=z[1])]
						set(tmp, tmp2, 'POS', tmp[tmp2, POS+1L])
						z	<- z[-1]
					}
					set(tmp, NULL, 'POS', tmp[, paste('PR_',POS,sep='')])					
					set(tmp, tmp[, which(MISS)], 'NT_DIFF', NA_integer_)
					set(tmp, tmp[, which(ANY)], 'NT_DIFF', NA_integer_)
					tmp[, MISS:=NULL]
					tmp[, ANY:=NULL]
					#tmp
					list(TAXA=tmp$TAXA, POS=tmp$POS, NT_DIFF=tmp$NT_DIFF)
				}, by='PR' ]
		
		set(dpand, NULL, 'POS', dpand[, factor(POS, levels=paste('PR_',1:dpand[, length(unique(POS))],sep=''))])
		set(dpand, NULL, 'TAXA', dpand[, as.character(TAXA)])
		set(dpand, NULL, 'POS', dpand[, as.character(POS)])
		#	update 1R
		tmp		<- dcast.data.table(subset(dpand, PR=='1Ra' | PR=='1Rb'), TAXA+POS~PR, value.var='NT_DIFF')
		setnames(tmp, c('1Ra','1Rb'), c('pr1Ra','pr1Rb'))
		tmp[, NT_DIFF:= pr1Ra*pr1Rb]
		tmp[, PR:='1R']
		set(tmp, NULL, c('pr1Ra','pr1Rb'), NULL)
		dpand	<- rbind(subset(dpand, PR!='1R'), tmp)
		#	update 2F
		tmp		<- dcast.data.table(subset(dpand, PR=='2Fa' | PR=='2Fb'), TAXA+POS~PR, value.var='NT_DIFF')
		setnames(tmp, c('2Fa','2Fb'), c('pr2Fa','pr2Fb'))
		tmp[, NT_DIFF:= pr2Fa*pr2Fb]
		tmp[, PR:='2F']
		set(tmp, NULL, c('pr2Fa','pr2Fb'), NULL)
		dpand	<- rbind(subset(dpand, PR!='2F'), tmp)
		#	update 2R
		tmp		<- dcast.data.table(subset(dpand, PR=='2Ra' | PR=='2Rb'), TAXA+POS~PR, value.var='NT_DIFF')
		setnames(tmp, c('2Ra','2Rb'), c('pr2Ra','pr2Rb'))
		tmp[, NT_DIFF:= pr2Ra*pr2Rb]
		tmp[, PR:='2R']
		set(tmp, NULL, c('pr2Ra','pr2Rb'), NULL)
		dpand	<- rbind(subset(dpand, PR!='2R'), tmp)
		#
		#	calculate number of '?' in each of the gene regions
		#
		z		<- as.character( sqp )
		dgd		<- dgene[, {
					#START<- 812; END_B4NXT<- 4341
					tmp			<- z[, seq.int(START, END)]
					tmp			<- apply(tmp=='?', 1, sum)
					list(	TAXA= names(tmp), UNASS= tmp/(END-START+1L) 	)
				}, by=c('GENE','START','END','LEN')]
		#	add number not '?','n','-'
		tmp		<- dgene[, {					
					#START<- 812; END_B4NXT<- 4341
					tmp			<- z[, seq.int(START, END)]
					tmp			<- apply(!(tmp=='?'|tmp=='n'|tmp=='-'), 1, sum)
					list(	TAXA= names(tmp), ACTG= tmp/(END-START+1L) 	)
				}, by=c('GENE','START','END')]
		dgd		<- merge(dgd, tmp, by=c('GENE','START','END','TAXA'))		
		dgd		<- subset(dgd, !grepl('HXB2',TAXA))
		#
		save(sq, sqp, sqi, dgene, dpan, dpand, dgd, file=file.path(wdir,wfile))
	}
	#
	#	get COMET subtypes
	#
	if(0)
	{
		infile		<- "~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/explaingaps/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_region_1F1R_COMETv0.5.txt"
		dc			<- as.data.table(read.table(infile, header=TRUE, sep='\t', stringsAsFactors=FALSE, strip.white=FALSE))
		dc[, COMET_R:='COMET_1F1R']
		dc[, COMET_V:='0.5']
		infile		<- "~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/explaingaps/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_region_3F4F_COMETv0.5.txt"
		tmp			<- as.data.table(read.table(infile, header=TRUE, sep='\t', stringsAsFactors=FALSE, strip.white=FALSE))
		tmp[, COMET_R:='COMET_3F4F']
		tmp[, COMET_V:='0.5']
		dc			<- rbind(dc,tmp)		
		infile		<- "~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/explaingaps/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_region_4F3R_COMETv0.5.txt"
		tmp			<- as.data.table(read.table(infile, header=TRUE, sep='\t', stringsAsFactors=FALSE, strip.white=FALSE))
		tmp[, COMET_R:='COMET_4F3R']
		tmp[, COMET_V:='0.5']
		dc			<- rbind(dc,tmp)
		infile		<- "~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/explaingaps/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_region_1F1R_COMETv2.1.txt"
		tmp			<- as.data.table(read.table(infile, header=TRUE, sep='\t', stringsAsFactors=FALSE, strip.white=FALSE))
		tmp[, COMET_R:='COMET_1F1R']
		tmp[, COMET_V:='2.1']
		dc			<- rbind(dc,tmp,use.name=TRUE,fill=TRUE)
		infile		<- "~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/explaingaps/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_region_3F4F_COMETv2.1.txt"
		tmp			<- as.data.table(read.table(infile, header=TRUE, sep='\t', stringsAsFactors=FALSE, strip.white=FALSE))
		tmp[, COMET_R:='COMET_3F4F']
		tmp[, COMET_V:='2.1']
		dc			<- rbind(dc,tmp,use.name=TRUE,fill=TRUE)		
		infile		<- "~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/explaingaps/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_region_4F3R_COMETv2.1.txt"
		tmp			<- as.data.table(read.table(infile, header=TRUE, sep='\t', stringsAsFactors=FALSE, strip.white=FALSE))
		tmp[, COMET_R:='COMET_4F3R']
		tmp[, COMET_V:='2.1']
		dc			<- rbind(dc,tmp,use.name=TRUE,fill=TRUE)
		infile		<- "~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/explaingaps/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_UGfull_COMETv2.1.txt"
		tmp			<- as.data.table(read.table(infile, header=TRUE, sep='\t', stringsAsFactors=FALSE, strip.white=FALSE))
		tmp[, COMET_R:='COMET_full']
		tmp[, COMET_V:='2.1']
		dc			<- rbind(dc,tmp,use.name=TRUE,fill=TRUE)
		infile		<- "~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/explaingaps/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_ZABWfull_COMETv2.1.txt"
		tmp			<- as.data.table(read.table(infile, header=TRUE, sep='\t', stringsAsFactors=FALSE, strip.white=FALSE))
		tmp[, COMET_R:='COMET_full']
		tmp[, COMET_V:='2.1']
		dc			<- rbind(dc,tmp,use.name=TRUE,fill=TRUE)
		dc[, x:=NULL]
		setnames(dc, c('bootstrap.support','name','subtype','length'), c('COMET_BS','TAXA','COMET_ST','COMET_N'))
		dc			<- subset(dc, !is.na(TAXA))
		dc			<- subset(dc, !grepl('HXB2',TAXA))	
		#	reset length from dgd
		tmp			<- subset(dgd, GENE%in%c('START-1R','3F-4F','4F-3R','GAG+POL+ENV'))
		set(tmp,NULL,'ACTG',tmp[,ACTG*LEN])
		set(tmp,tmp[, which(GENE=='START-1R')],'GENE','COMET_1F1R')
		set(tmp,tmp[, which(GENE=='3F-4F')],'GENE','COMET_3F4F')
		set(tmp,tmp[, which(GENE=='4F-3R')],'GENE','COMET_4F3R')
		set(tmp,tmp[, which(GENE=='GAG+POL+ENV')],'GENE','COMET_full')		
		setnames(tmp,c('GENE','ACTG'),c('COMET_R','COMET_ACTG'))
		dc			<- merge(subset(dc,select=c('TAXA','COMET_ST','COMET_R','COMET_N','COMET_V')),subset(tmp,select=c('TAXA','COMET_R','COMET_ACTG')),by=c('TAXA','COMET_R'))
		if(0)
		{
			set(dc, dc[, which(grepl('unassigned', COMET_ST))], 'COMET_ST', 'unassigned')
			set(dc, dc[, which(grepl('_', COMET_ST))], 'COMET_ST', 'pot_recombinant')
			set(dc, dc[, which(!COMET_ST%in%c('A1','B','C','D','unassigned','pot_recombinant'))], 'COMET_ST', 'other')			
		}
		if(1)	#based on observation that on full genome, all recombinants are 'unassigned'
		{
			set(dc, dc[, which(grepl('unassigned', COMET_ST))], 'COMET_ST', 'pot_recombinant')
			set(dc, dc[, which(grepl('_', COMET_ST))], 'COMET_ST', 'pot_recombinant')
			set(dc, dc[, which(!COMET_ST%in%c('A1','B','C','D','pot_recombinant'))], 'COMET_ST', 'other')						
		}
		#	cross-compare assignments: assignments largely agree, 
		#	except that COMET2.1 now makes more assignments on short sequences, that are mostly assigned A
		dc[, table(COMET_ST,COMET_V,COMET_R)]
		#	use actual length	
		#set(dc, dc[, which(COMET_N<500 & COMET_V=='0.5')], 'COMET_ST', 'short')
		set(dc, dc[, which(COMET_ACTG<500 & COMET_V=='0.5')], 'COMET_ST', 'short')
		set(dc, dc[, which(COMET_ACTG<500 & COMET_V=='2.1')], 'COMET_ST', 'short')
		#	use COMETv2.1
		dc			<- subset(dc, COMET_V=='2.1',select=c('TAXA','COMET_R','COMET_ST','COMET_ACTG'))
		setnames(dc, 'COMET_ACTG','COMET_N')
		tmp			<- dc[, {
					ans	<- 'pot_recombinant'
					z	<- which(!COMET_ST%in%c('short',"unassigned"))
					if(all(COMET_ST[z]=="A1"))
						ans<- 'A1'
					if(all(COMET_ST[z]=="B"))
						ans<- 'B'
					if(all(COMET_ST[z]=="C"))
						ans<- 'C'
					if(all(COMET_ST[z]=="D"))
						ans<- 'D'
					if(all(COMET_ST[z]=="other"))
						ans<- 'other'
					if(length(which(COMET_ST=="unassigned"))>1)
						ans<- 'unassigned'
					if(all(COMET_ST=="short"))
						ans<- 'short'
					list(COMET_R='COMET_CONS',COMET_ST=ans, COMET_N=min(COMET_N))
				}, by='TAXA']
		dc			<- rbind(dc, tmp, use.names=TRUE)
		if(0)
		{
			#	check COMET subtype assignments on ZA|BW full
			tmp			<- merge(dc, subset(dc, COMET_R=='COMET_full' & grepl('ZA|BW',TAXA), TAXA), by='TAXA')
			tmp			<- dcast.data.table(tmp, TAXA~COMET_R, value.var='COMET_ST')
			#	check COMET subtype assignments on UG full
			tmp			<- merge(dc, subset(dc, COMET_R=='COMET_full' & grepl('UG',TAXA), TAXA), by='TAXA')
			tmp			<- dcast.data.table(tmp, TAXA~COMET_R, value.var='COMET_ST')
			tmp[, COMET_CLASS:='agree_pure_subtype']
			set(tmp, tmp[, which(COMET_full=='unassigned' & COMET_CONS%in%c('unassigned','pot_recombinant'))], 'COMET_CLASS', 'agree_pot_recombinant')
			set(tmp, tmp[, which(COMET_full=='unassigned' & !COMET_CONS%in%c('unassigned','pot_recombinant'))], 'COMET_CLASS', 'by_three_regions:_potentially_wrong_pure_subtype')
			set(tmp, tmp[, which(COMET_full!='unassigned' & COMET_CONS%in%c('unassigned','pot_recombinant'))], 'COMET_CLASS', 'by_three_regions:_potentially_wrong_pot_recombinant')
			#	tmp[, table(COMET_CLASS)]
			#                      agree_pot_recombinant                                  agree_pure_subtype by_three_regions:_potentially_wrong_pot_recombinant		by_three_regions:_potentially_wrong_pure_subtype 
			#                                        342                                                 343                                                   2                                                  28 
			#	prop of pure subtypes assignments that may be wrong: 28/(343+28)=0.0754717
			#	prop of pot recombinant and unassigned in 715 from UG: 370/715=0.517
			#	prop of pot recombinant and unassigned in 3628 from UG: (579+104)/3628=0.188
			
		}		
		#	ignore check analysis on 'COMET_full'
		dc			<- subset(dc, COMET_R!='COMET_full')
		tmp2		<- dcast.data.table(dc, TAXA~COMET_R, value.var='COMET_N')
		setnames(tmp2, setdiff(colnames(tmp2),'TAXA'), paste(setdiff(colnames(tmp2),'TAXA'),'_N',sep=''))
		dc			<- dcast.data.table(dc, TAXA~COMET_R, value.var='COMET_ST')
		dc			<- merge(dc, tmp2, by='TAXA')
		#
	}
	#
	#	add meta-data
	#
	if(0)
	{
		dm	<- unique(subset(dgd, select=TAXA))
		dm[, PANGEA_ID:=gsub('-S.*','',TAXA)]		
		#	add RCCS data
		load("~/Dropbox (Infectious Disease)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/RakaiPangeaMetaData.rda")
		rccsData	<- as.data.table(rccsData)		
		rccsData	<- subset(rccsData, select=c('Pangea.id', 'RCCS_studyid', 'date', 'batch', 'birthyr', 'REGION', 'COMM_NUM', 'HH_NUM', 'SEX', 'AGEYRS','firstPosDate', 'arvStartDate', 'selfReportArt', 'everSelfReportArt', 'FirstSelfReportArt', 'recentVL', 'recentVLdate'))
		setnames(rccsData, c('Pangea.id','RCCS_studyid','REGION','date','birthyr','firstPosDate','arvStartDate','recentVL','recentVLdate',"AGEYRS"), c('PANGEA_ID','STUDY_ID','LOC','SAMPLEDATE','DOB','FIRSTPOSDATE','ARTSTART','RECENTVL','RECENTVLDATE',"AGE"))
		rccsData	<- subset(rccsData, !is.na(PANGEA_ID))
		setkey(rccsData, PANGEA_ID)
		rccsData	<- unique(rccsData)	#remove 4 duplicates "K104085" "E106462" "F030186" "F101874"
		rccsData[, COHORT:='RCCS']
		set(rccsData, NULL, 'SAMPLEDATE', rccsData[,hivc.db.Date2numeric(SAMPLEDATE)])	
		set(rccsData, NULL, 'FIRSTPOSDATE', rccsData[,hivc.db.Date2numeric(FIRSTPOSDATE)])
		set(rccsData, NULL, 'ARTSTART', rccsData[,hivc.db.Date2numeric(ARTSTART)])
		set(rccsData, NULL, 'FirstSelfReportArt', rccsData[,hivc.db.Date2numeric(FirstSelfReportArt)])
		set(rccsData, NULL, 'RECENTVLDATE', rccsData[,hivc.db.Date2numeric(RECENTVLDATE)])
		#	add Mochudi data
		load("~/duke/2016_AC/PANGEA_160826/160826_PANGEA_BW_corevariables_n373.rda")
		setnames(bwp, 'PANGEAID', 'PANGEA_ID')
		set(bwp, NULL, 'SEQID', NULL)
		tmp			<- rbind(rccsData, bwp, use.names=TRUE, fill=TRUE) 
		#	add AC data
		load("~/duke/2016_AC/PANGEA_160826/160826_PANGEA_AC_corevariables_n2940.rda")
		setnames(acp, 'PANGEAID', 'PANGEA_ID')
		set(acp, NULL, c('REASONSAMPLING','LATESTARTREGIMEN','CIRCUMCISED','LASTNEGDATE','LASTNUMSEXUALPARTNERS','LATESTARTREGIMENSTARTED'), NULL)		
		tmp			<- rbind(tmp, acp, use.names=TRUE, fill=TRUE)
		dm			<- merge(dm, tmp, by='PANGEA_ID', all.x=1)
		#	add CLASS subtype data
		subtypeSummaryData	<- as.data.table(subtypeSummaryData)
		setnames(subtypeSummaryData, 'RCCS_studyid', 'STUDY_ID')
		dst					<- subset(subtypeSummaryData, select=c(STUDY_ID, gp.class, gag.class, pol.class, vpu.class, env.class, comp.class))
		setkey(dst, STUDY_ID)
		dst					<- unique(dst)	#no duplicates here
		dst					<- melt(dst, id.var='STUDY_ID')
		set(dst, NULL, 'value', dst[, gsub('\\s','',gsub('Complex','',gsub('Subtype', '',value)))])
		set(dst, NULL, 'variable', dst[, paste('SUBTYPE_',toupper(gsub('\\.class','',variable)),sep='')])	
		dst					<- dcast.data.table(dst, STUDY_ID~variable)
		dm					<- merge(dm, dst, by='STUDY_ID',all.x=1)
		#	add COMET subtype data
		dm					<- merge(dm, dc, by='TAXA',all.x=1)
		set(dm, dm[, which(is.na(COMET_1F1R))],'COMET_1F1R_N',0L)
		set(dm, dm[, which(is.na(COMET_1F1R))],'COMET_1F1R','short')
		set(dm, dm[, which(is.na(COMET_3F4F))],'COMET_3F4F_N',0L)
		set(dm, dm[, which(is.na(COMET_3F4F))],'COMET_3F4F','short')
		set(dm, dm[, which(is.na(COMET_4F3R))],'COMET_4F3R_N',0L)
		set(dm, dm[, which(is.na(COMET_4F3R))],'COMET_4F3R','short')
		#	add REGA subtype data
		load("~/Dropbox (Infectious Disease)/PANGEA_alignments/Rega Subtype Analysis/Gag REGA results/regaGag.rda")
		regaGag		<- as.data.table(regaGag)
		regaGag		<- subset(regaGag, select=c(TAXA, assignment, support, pure, pure_support))
		set(regaGag, NULL, 'assignment', regaGag[, gsub('Check the Report|Check the report|Check the bootscan','Check bootscan', gsub('Subtype ','', gsub('HIV-1 ','', assignment)))])
		setnames(regaGag, c('assignment','support','pure','pure_support'),c('REGA_GAG_A','REGA_GAG_AS','REGA_GAG_PURE','REGA_GAG_PURES'))
		dm			<- merge(dm, regaGag, by='TAXA',all.x=1)
		#	add Sanger processing data etc
		dc			<- data.table(read.csv("~/Dropbox (Infectious Disease)/PANGEA_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v150908_Summary.csv"))		
		setnames(dc, c('Sanger.ID','PANGEA.ID','reference.for.mapping','clinical.genome.coverage'), c('SANGER_ID','PANGEA_ID','REF_4_MAPPING','COV'))		
		tmp			<- data.table(read.csv('~/Dropbox (Infectious Disease)/PANGEA_data/PAN_iva_dependencies_9861.txt', sep='\t'))
		setnames(tmp, 'LaneID', 'SANGER_ID')
		set(tmp, NULL, 'SANGER_ID', tmp[, gsub('#','_',SANGER_ID)])
		dc			<- merge(dc, tmp, by='SANGER_ID', all.x=1)
		set(dc, NULL, 'TAXA', dc[, gsub('\\s$','',gsub('^\\s','',as.character(PANGEA_ID)))])
		set(dc, NULL, 'PANGEA_ID', dc[, gsub('-S.*','',TAXA)])		
		#	some of the PANGEA IDs are duplicates and cannot be resolved. Use the coverage.
		tmp			<- as.character(sqp)
		tmp2		<- apply(tmp, 1, function(x) sum(!x%in%c('?','-')))
		tmp2		<- data.table(TAXA=names(tmp2), PANGEA_ID=gsub('-S.*','',names(tmp2)), COV=tmp2)
		z			<- subset(tmp2, grepl('-R2',TAXA))
		#	some manual adjustments to get this to match
		set(tmp2, tmp2[, which(TAXA=='PG14-BW000057-S01150-R2')], 'COV', 8061L)
		set(tmp2, tmp2[, which(TAXA=='PG14-BW000058-S01151-R2')], 'COV', 8090L)
		set(tmp2, tmp2[, which(TAXA=='PG14-BW000059-S01152-R2')], 'COV', 8554L)
		set(tmp2, tmp2[, which(TAXA=='PG14-BW000063-S01156-R2')], 'COV', 7891L)
		set(tmp2, tmp2[, which(TAXA=='PG14-BW000064-S01157-R2')], 'COV', 7872L)		
		set(tmp2, tmp2[, which(TAXA=='PG14-UG002291-S00291-R2')], 'COV', 7392L)
		set(tmp2, tmp2[, which(TAXA=='PG14-UG500523-S02527-R2')], 'COV', 3496L)		
		set(tmp2, tmp2[, which(TAXA=='PG14-UG500526-S02530-R2')], 'COV', 7664L)
		set(tmp2, tmp2[, which(TAXA=='PG14-UG500526-S02530')], 'COV', 7322L)		
		set(tmp2, tmp2[, which(TAXA=='PG14-UG500529-S02533-R2')], 'COV', 7910L)
		set(tmp2, tmp2[, which(TAXA=='PG14-UG500536-S02540-R2')], 'COV', 8210L)
		set(tmp2, tmp2[, which(TAXA=='PG14-UG500541-S02545-R2')], 'COV', 5086L)		
		z			<- subset(tmp2, grepl('-R2',TAXA))
		z			<- merge(subset(dc, COV>0, select=which(colnames(dc)!='TAXA')), subset(z,grepl("PG",TAXA) & COV>0, select=c(PANGEA_ID,TAXA,COV)), by=c('PANGEA_ID','COV'))
		setnames(z, 'TAXA','TAXA_NEW')
		dc			<- merge(dc, subset(z, select=c(SANGER_ID,TAXA_NEW)), by='SANGER_ID', all.x=1)
		tmp2		<- dc[, which(!is.na(TAXA_NEW))]
		set(dc, tmp2, 'TAXA', dc[tmp2,TAXA_NEW])
		dc			<- subset(dc, COV>0)
		setkey(dc, TAXA)
		set(dc, NULL, c('TAXA_NEW','COV','PANGEA_ID'), NULL)				
		dm			<- merge(dm, dc, by='TAXA', all.x=1)
		stopifnot( nrow(subset(dm, is.na(SANGER_ID)))==0 )
		#	complete COHORT
		set(dm, dm[, which(grepl('UG',PANGEA_ID) & is.na(COHORT))], 'COHORT', 'UG-MRC')
		set(dm, dm[, which(grepl('ZA',PANGEA_ID) & is.na(COHORT))], 'COHORT', 'AC_Resistance')
		#	add extraction ID		
		dm[, EXTRACT_ID:= dm[,gsub('-S','',regmatches(TAXA,regexpr('-S[0-9]+', TAXA)))]]		
		#	save
		save(sq, sqp, sqi, dgene, dpan, dpand, dgd, dm, file=file.path(wdir,wfile))		
	}	
}
##--------------------------------------------------------------------------------------------------------
##	olli 27.06.11
treecomparison.bootstrap.mvr.dev<- function(indir=NULL, wdir=NULL)
{
	require(ape)
	require(data.table)
	require(recosystem)
	require(ggplot2)
	#	get master RDA file with all distances
	if(0)	
	{
		wdir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/tree_mvr'	
		indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/simulations'	
		infile	<- '150701_Regional_TRAIN4_SIMULATED.fa'
		#infile	<- '150701_Regional_TRAIN2_SIMULATED.fa'
		#	create tp with IDs -- need this to complete to matrix	
		seq		<- read.dna(file.path(indir, infile), format='fa')
		tp		<- as.data.table( t(combn(rownames(seq),2)) )		
		setnames(tp, c('V1','V2'), c('TAXA1','TAXA2'))
		tmp		<- as.data.table( t(combn(seq_len(nrow(seq)),2)) )
		setnames(tmp, c('V1','V2'), c('ID1','ID2'))
		tp		<- cbind(tp, tmp)
		#
		#	read genetic distances between taxon pairs
		#		
		infiles	<- data.table(FILE=list.files(wdir, pattern='BATCH[0-9]+.rda$', full.names=TRUE))
		infiles[, BATCH:= as.integer(gsub('BATCH','',regmatches(FILE, regexpr('BATCH[0-9]+', FILE))))]
		setkey(infiles, BATCH)
		#	not yet completed
		stopifnot( infiles[, length(setdiff(seq.int(1,400), BATCH))==0]	)
		#	read files
		tmp		<- lapply(infiles[, FILE], function(x)
				{
					load(x)
					ans	<- merge(tp, tpi, by=c('TAXA1','TAXA2'))
					ans
				})
		tp		<- do.call('rbind', tmp)	
		setkey(tp, ID1, ID2)
		tp[, GD_V:= GD_SD*GD_SD]
		#
		#	save tp to file
		#		
		save(tp, seq, file=file.path(wdir, gsub('\\.fa','_tps.rda',infile)))	
	}	
	if(0)
	{
		wdir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/tree_mvr'
		df		<- data.table(FILE=list.files(wdir, pattern='newick$',full.names=1))
		tmp		<- df[, {
					ph2	<- read.tree(FILE)
					list(CH= nrow(unique(data.table(TAXA=ph2$tip.label)))==Ntip(ph2) )
				}, by='FILE']
		
	}
	if(0)
	{
		na.rm.p						<- NA	 
		complete.distance.matrix	<- 0
		seed						<- 123		
		v.mult						<- 1.2
		reco.opts					<- c(dim=750, costp_l1=0, costp_l2=0.001, costq_l1=0, costq_l2=0.001, nthread=1, lrate=0.003, niter=120)
		verbose						<- 1
		
		#wdir	<- '/work/or105/Gates_2014/tree_comparison/mvr'
		wdir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/tree_mvr'
		infile	<- '150701_Regional_TRAIN4_SIMULATED_tps.rda'
		#infile	<- '150701_Regional_TRAIN2_SIMULATED_tps.rda'
		load(file.path(wdir, infile))
		
		loop.rep	<- tp[, unique(REP)]
		loop.gene	<- tp[, unique(GENE)]
		loop.gene	<- "gag+pol+env"
				
		for(gene in loop.gene)
			for(rep in loop.rep)
			{
				#gene	<- "env"; rep<- 1
				#gene	<- "gag+pol+env"; rep<- 1
				tps		<- subset(tp, GENE==gene & REP==rep)
				outfile	<- file.path(wdir, gsub('\\.rda',paste('_GENE_',gene,'_REP_',rep,'_C_',complete.distance.matrix,sep=''), infile))
				tmp		<- seq.mvr.d.and.v(tps, seed=seed, v.mult=v.mult, complete.distance.matrix=complete.distance.matrix, reco.opts=reco.opts, outfile=outfile, verbose=verbose)				
				d		<- tmp$d
				v		<- tmp$v
				tmp		<- NULL
				gc()
				#	write to file
				d		<- as.matrix(d)
				v		<- as.matrix(v)
				d[is.na(d)]	<- -1
				v[is.na(v)]	<- -1
				file.d	<- paste(gsub('\\.rda|\\.newick|\\.tree','',outfile),'_d.phylip',sep='')
				file.v	<- paste(gsub('\\.rda|\\.newick|\\.tree','',outfile),'_v.phylip',sep='')
				seq.write.dna.phylip.triangular(d, file=file.d)
				seq.write.dna.phylip.triangular(v, file=file.v)
				#	call PhyD*
				tmp		<- cmd.phydstar(file.d, outfile=outfile, method='BioNJ', fs=15, binary=TRUE, negative.branch.length=FALSE, lower.triangular=TRUE)
				cat(tmp)
				tmp		<- cmd.phydstar(file.d, outfile=outfile, infile.v=file.v, method='MVR', fs=15, binary=TRUE, negative.branch.length=FALSE, lower.triangular=TRUE)
				system(tmp)
				
				outfile	<- paste(outfile,'_',paste(reco.opts,collapse='_'),'_mvr.newick',sep='')
				options(expressions=5e5)
				write.tree(ph, file=outfile)	
				options(expressions=5e3)
			}
		quit('no')
	}
	if(0)
	{
		tps				<- subset(tp, REP==1 & GENE=='gag+pol+env')		
		tmp				<- dcast.data.table(tps, ID1~ID2, value.var='GD')		
		d				<- cbind(NA_real_, as.matrix(tmp[, -1, with=FALSE]))
		d				<- rbind(d, NA_real_)
		colnames(d)[1]	<- setdiff( as.character(tmp[, ID1]), colnames(d) )
		rownames(d)		<- colnames(d)
		diag(d)			<- 0
		#	complete lower triangular from upper triangular and vice versa
		tmp				<- lower.tri(d) & is.na(d)	
		d[tmp]			<- t(d)[tmp]
		tmp				<- upper.tri(d) & is.na(d)
		d[tmp]			<- t(d)[tmp]
		#	reset names
		tmp				<- subset( tps, select=c(TAXA1, ID1) )
		setnames(tmp, c('TAXA1','ID1'), c('TAXA2','ID2') )
		tmp				<- unique(rbind( tmp, subset( tps, select=c(TAXA2, ID2) ) ))
		setnames(tmp, c('TAXA2','ID2'), c('TAXA','ID') )
		setkey(tmp, ID)
		rownames(d)		<- tmp[, TAXA]
		colnames(d)		<- tmp[, TAXA]
		#	checks	
		stopifnot(ncol(d)==nrow(d))
		stopifnot(length(which(is.na(d)))==2*nrow(subset(tps, is.na(GD))))
		cat('D matrix: proportion of NA entries=',length(which(is.na(d))) / prod(dim(d)))
		#
		#	get variance matrix
		#
		tmp				<- dcast.data.table(tps, ID1~ID2, value.var='GD_V')		
		v				<- cbind(NA_real_, as.matrix(tmp[, -1, with=FALSE]))
		v				<- rbind(v, NA_real_)
		colnames(v)[1]	<- setdiff( as.character(tmp[, ID1]), colnames(v) )
		rownames(v)		<- colnames(v)
		diag(v)			<- 0
		#	complete lower triangular from upper triangular and vice versa
		tmp				<- lower.tri(v) & is.na(v)	
		v[tmp]			<- t(v)[tmp]
		tmp				<- upper.tri(v) & is.na(v)
		v[tmp]			<- t(v)[tmp]
		#	reset names
		tmp				<- subset( tps, select=c(TAXA1, ID1) )
		setnames(tmp, c('TAXA1','ID1'), c('TAXA2','ID2') )
		tmp				<- unique(rbind( tmp, subset( tps, select=c(TAXA2, ID2) ) ))
		setnames(tmp, c('TAXA2','ID2'), c('TAXA','ID') )
		setkey(tmp, ID)
		rownames(v)		<- tmp[, TAXA]
		colnames(v)		<- tmp[, TAXA]
		#	checks	
		stopifnot(ncol(v)==nrow(v))	
		cat('V matrix: proportion of NA entries=',length(which(is.na(v))) / prod(dim(v)))
		#
		#	remove cols/rows that contain nothing else than NAs
		#	
		diag(d)			<- NA_real_
		diag(v)			<- NA_real_		
		tmp				<- apply(d, 1, function(x) !all(is.na(x)))
		cat('\nIn D: found',length(which(!tmp)),'columns / rows with NA only: remove in D and V. ', rownames(d)[!tmp])
		ds				<- d[tmp,tmp]
		vs				<- v[tmp,tmp]	
		tmp				<- apply(vs, 1, function(x) !all(is.na(x)))
		cat('\nIn V: found additional',length(which(!tmp)),'columns / rows with NA only: remove in D and V too. ', rownames(d)[!tmp])
		ds				<- ds[tmp,tmp]
		vs				<- vs[tmp,tmp]			
	}
	if(0)
	{
		wdir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/tree_mvr'
		load(file.path(wdir, '150701_Regional_TRAIN4_SIMULATED_tps.rda'))
		tps		<- subset(tp, REP==1 & GENE=='gag+pol+env', select=c(ID1, ID2, GD))
		#	add upper triangular
		tmp		<- copy(tps)
		set(tmp, NULL, 'ID1', tps[, ID2])
		set(tmp, NULL, 'ID2', tps[, ID1])
		tps		<- rbind(tps, tmp)
		#	add zero diagonal
		tmp		<- tps[, range(ID1)]
		tmp		<- data.table(ID1= seq.int(tmp[1], tmp[2]), ID2= seq.int(tmp[1], tmp[2]), GD=0)
		tps		<- rbind(tps, tmp)
		#	ignore NA entries
		tps.na	<- subset(tps, is.na(GD))
		tps		<- subset(tps, !is.na(GD))
		tps.all	<- rbind(tps, tps.na)
		#	setup matrix completion
		tmp		<- data_memory(tps[,ID1], tps[,ID2], rating=tps[,GD], index1=TRUE)				
		set.seed(123)
		r		<- Reco()
		opts	<- r$tune(tmp, opts=list(dim=c(10, 100, 500, 750), lrate=c(0.01), costp_l1=0, costp_l2=c(0.001, 0.01, 0.1), costq_l1=0, costq_l2=c(0.001, 0.01, 0.1), nthread=1, niter=10))
		#best is  dim=750 costp_l2=0.001 costq_l2=0.001  lrate=0.01 rmse=0.01457322
		opts	<- r$tune(tmp, opts=list(dim=c(500, 750, 1000), lrate=c(0.001, 0.01), costp_l1=0, costp_l2=c(0.0001, 0.001), costq_l1=0, costq_l2=c(0.0001, 0.001), nthread=1, niter=10))
		#best is  dim=1000 costp_l2=0.0001 costq_l2=0.0001  lrate=0.01 rmse=0.01457322
		opts	<- r$tune(tmp, opts=list(dim=c(100, 500), lrate=c(0.003, 0.005), costp_l1=0, costp_l2=c(0.01), costq_l1=0, costq_l2=c(0.01), nthread=1, niter=100))
		
		r$train(tmp, opts=c(dim=500, costp_l1=0, costp_l2=0.01, costq_l1=0, costq_l2=0.01, nthread=1, lrate=0.003, niter=40))
		#rmse 0.0211
		tps.all[, GDp:= r$predict(data_memory(tps.all[,ID1], tps.all[,ID2], index1=TRUE), out_memory())]
		#plot
		ggplot(subset(tps.all, !is.na(GD)), aes(x=GD, y=GDp)) + geom_point(colour='grey80', size=0.5, pch=16) + geom_abline(slope=1, intercept=0)
		ggsave(file=file.path(wdir, 'reco_500_1e-2_3e-3_40.pdf'), w=7, h=7)
		
		r$train(tmp, opts=c(dim=750, costp_l1=0, costp_l2=0.001, costq_l1=0, costq_l2=0.001, nthread=1, lrate=0.003, niter=120))
		#rmse 0.0155
		tps.all[, GDp2:= r$predict(data_memory(tps.all[,ID1], tps.all[,ID2], index1=TRUE), out_memory())]
		#plot
		ggplot(subset(tps.all, !is.na(GD)), aes(x=GD, y=GDp2)) + geom_point(colour='grey80', size=0.5, pch=16) + geom_abline(slope=1, intercept=0)
		ggsave(file=file.path(wdir, 'reco_750_1e-3_3e-3_120.pdf'), w=7, h=7)
		
		#	fill in distance matrix
		tps.all[, GDf:= GD]
		tmp		<- tps.all[, which(is.na(GDf))]
		set(tps.all, tmp, 'GDf', tps.all[tmp, GDp2])
		#	convert to matrix (not necessarily symmetric)
		tmp			<- dcast.data.table( subset(tps.all, select=c(ID1,ID2,GDf)), ID1~ID2, value.var='GDf' )		
		d			<- as.matrix(tmp[, -1, with=FALSE])
		rownames(d)	<- colnames(d)
		#	make symmetric
		d			<- (d+t(d))/2	
		#	some rows/cols may have NAs only -- remove these as the matrix completion problem is ill-specified for these
		tmp			<- subset(tps.na[, list(GDM=length(GD)), by='ID1'], GDM==nrow(d)-1) #subtract one since diagonal is zero
		tmp			<- setdiff(rownames(d), tmp[, as.character(ID1)] )
		d			<- d[tmp, tmp]		
		#
		#	generate variance matrix
		#
		tps				<- subset(tp, REP==1 & GENE=='gag+pol+env', select=c(TAXA1, ID1, TAXA2, ID2, GD_V))	
		tmp				<- dcast.data.table(tps, ID1~ID2, value.var='GD_V')		
		v				<- cbind(NA_real_, as.matrix(tmp[, -1, with=FALSE]))
		v				<- rbind(v, NA_real_)
		colnames(v)[1]	<- setdiff( as.character(tmp[, ID1]), colnames(v) )
		rownames(v)		<- colnames(v)
		diag(v)			<- 0
		#	complete lower triangular from upper triangular and vice versa
		tmp				<- lower.tri(v) & is.na(v)	
		v[tmp]			<- t(v)[tmp]
		tmp				<- upper.tri(v) & is.na(v)
		v[tmp]			<- t(v)[tmp]
		#	set missing variances to large default
		v[is.na(v)]		<- max(v, na.rm=TRUE)*1.2
		v				<- v[rownames(d),colnames(d)]		
		#
		#	reset names
		#
		tmp				<- subset( tps, select=c(TAXA1, ID1) )
		setnames(tmp, c('TAXA1','ID1'), c('TAXA2','ID2') )
		tmp				<- unique(rbind( tmp, subset( tps, select=c(TAXA2, ID2) ) ))
		setnames(tmp, c('TAXA2','ID2'), c('TAXA','ID') )		
		tmp				<- merge(tmp, data.table(ID=as.integer(rownames(d))), by='ID')
		setkey(tmp, ID)
		rownames(d)		<- tmp[, TAXA]
		colnames(d)		<- tmp[, TAXA]		
		rownames(v)		<- tmp[, TAXA]
		colnames(v)		<- tmp[, TAXA]
		#				
		#	run mvr with completed distance and variance matrices
		#
		d				<- as.dist(d)
		v				<- as.dist(v)
		ph				<- mvr(d, v)

	}
	if(0)	#play with basic recosystem example
	{
		#	this is the example		 
		train_set	<- data_file(system.file("dat", "smalltrain.txt", package = "recosystem"))
		test_set	<- data_file(system.file("dat", "smalltest.txt",  package = "recosystem"))
		set.seed(123)
		r			<- Reco()
		opts		<- r$tune(train_set, opts=list(dim=c(10, 20, 30), lrate=c(0.1, 0.2), costp_l1=0, costq_l1=0, nthread=1, niter=10))
		r$train(train_set, opts = c(opts$min, nthread = 1, niter = 20))
		pred_rvec	<- r$predict(test_set, out_memory())
		
		test	 	<- read.table(test_set@source, sep = " ", header = FALSE)
		
		#	this is the same example but from memory
		infile		<- system.file("dat", "smalltrain.txt", package = "recosystem")
		dm			<- as.data.table(read.table(file=infile, sep=' '))
		setnames(dm, c('V1','V2','V3'), c('IDX1','IDX2','D'))				
		infile		<- system.file("dat", "smalltest.txt", package = "recosystem")
		dt			<- as.data.table(read.table(file=infile, sep=' '))
		setnames(dt, c('V1','V2'), c('IDX1','IDX2'))		
		dm.eco		<- data_memory(dm[,IDX1], dm[,IDX2], rating=dm[,D], index1=FALSE)
		dt.eco		<- data_memory(dt[,IDX1], dt[,IDX2], index1=FALSE)		
		set.seed(123)
		r			<- Reco()
		opts2		<- r$tune(dm.eco, opts=list(dim=c(10, 20, 30), lrate=c(0.1, 0.2), costp_l1=0, costq_l1=0, nthread=1, niter=10))
		r$train(dm.eco, opts = c(opts2$min, nthread = 1, niter = 20))
		pred_rvec2	<- r$predict(dt.eco, out_memory())		
		stopifnot( length(which( pred_rvec!=pred_rvec2 ))==0 )	#OK this works
		
		#	this is the same example but from memory and with ordered entries
		infile		<- system.file("dat", "smalltrain.txt", package = "recosystem")
		dm			<- as.data.table(read.table(file=infile, sep=' '))
		setnames(dm, c('V1','V2','V3'), c('IDX1','IDX2','D'))				
		infile		<- system.file("dat", "smalltest.txt", package = "recosystem")
		dt			<- as.data.table(read.table(file=infile, sep=' '))
		setnames(dt, c('V1','V2'), c('IDX1','IDX2'))	
		setkey(dm, IDX1, IDX2)
		setkey(dt, IDX1, IDX2)
		dm.eco		<- data_memory(dm[,IDX1], dm[,IDX2], rating=dm[,D], index1=FALSE)
		dt.eco		<- data_memory(dt[,IDX1], dt[,IDX2], index1=FALSE)		
		set.seed(123)
		r			<- Reco()
		opts3		<- r$tune(dm.eco, opts=list(dim=c(10, 20, 30), lrate=c(0.1, 0.2), costp_l1=0, costq_l1=0, nthread=1, niter=10))
		r$train(dm.eco, opts = c(opts3$min, nthread = 1, niter = 20))
		pred_rvec3	<- r$predict(dt.eco, out_memory())		
		stopifnot( length(which( pred_rvec!=pred_rvec3 ))==0 )	#not identical 
		
		#	reproduce RMSE manually
		infile		<- system.file("dat", "smalltrain.txt", package = "recosystem")
		dm			<- as.data.table(read.table(file=infile, sep=' '))
		setnames(dm, c('V1','V2','V3'), c('IDX1','IDX2','D'))
		dm.eco		<- data_memory(dm[,IDX1], dm[,IDX2], rating=dm[,D], index1=FALSE)
		infile		<- system.file("dat", "smalltest.txt", package = "recosystem")		
		dt			<- copy(dm)
		dt[, D:=NULL]		
		dt.eco		<- data_memory(dt[,IDX1], dt[,IDX2], index1=FALSE)		
		set.seed(123)
		r			<- Reco()
		opts4		<- r$tune(dm.eco, opts=list(dim=c(10, 20, 30), lrate=c(0.1, 0.2), costp_l1=0, costq_l1=0, nthread=1, niter=10))
		r$train(dm.eco, opts = c(opts4$min, nthread = 1, niter = 20))	#RMSE improves with niter
		dm[, PREDICT:= r$predict(dt.eco, out_memory())]
		subset(dm, !is.na(D))[, sqrt(mean((D-PREDICT)*(D-PREDICT)))]	#OK this works
		
		#	reproduce RMSE manually also when entries are ordered and all entries to be predicted?
		#	note: 	diagonal is not automatically considered zero, 
		#			and the matrix is not necessarily symmetric either!		
		infile		<- system.file("dat", "smalltrain.txt", package = "recosystem")
		dm			<- as.data.table(read.table(file=infile, sep=' '))
		setnames(dm, c('V1','V2','V3'), c('IDX1','IDX2','D'))
		setkey(dm, IDX1, IDX2)
		dm.eco		<- data_memory(dm[,IDX1], dm[,IDX2], rating=dm[,D], index1=FALSE)		
		tmp			<- dm[, range(IDX1)]	
		dp			<- as.data.table(expand.grid(IDX1=seq.int(tmp[1],tmp[2]), IDX2=seq.int(tmp[1],tmp[2])))		
		tmp			<- data_memory(dp[,IDX1], dp[,IDX2], index1=FALSE)
		set.seed(123)
		r			<- Reco()
		opts		<- r$tune(dm.eco, opts=list(dim=c(10, 20, 30), lrate=c(0.1, 0.2), costp_l1=0, costq_l1=0, nthread=1, niter=10))
		r$train(dm.eco, opts = c(opts$min, nthread = 1, niter = 20))	
		dp[, PREDICT:= r$predict(tmp, out_memory())]				
		dp			<- merge(dp, dm, by=c('IDX1','IDX2'), all.x=1)
		subset(dp, !is.na(D))[, sqrt(mean((D-PREDICT)*(D-PREDICT)))]	#OK this works too
	}
}
##--------------------------------------------------------------------------------------------------------
##	olli 27.06.11
seq.big.mvr<- function(tps, na.rm.p=NA, mds.args=list('ndim'= 750, type="mspline", "spline.intKnots"=3, "spline.degree"=2), wfile=NA)
{
	require(smacof)
	
	#	select (rep 1 gag+pol+env)
	tps				<- subset(tp, REP==1 & GENE=='gag+pol+env')
	#
	#	get distance matrix
	#
	tmp				<- dcast.data.table(tps, ID1~ID2, value.var='GD')		
	d				<- cbind(NA_real_, as.matrix(tmp[, -1, with=FALSE]))
	d				<- rbind(d, NA_real_)
	colnames(d)[1]	<- setdiff( as.character(tmp[, ID1]), colnames(d) )
	rownames(d)		<- colnames(d)
	diag(d)			<- 0
	#	complete lower triangular from upper triangular and vice versa
	tmp				<- lower.tri(d) & is.na(d)	
	d[tmp]			<- t(d)[tmp]
	tmp				<- upper.tri(d) & is.na(d)
	d[tmp]			<- t(d)[tmp]
	#	reset names
	tmp				<- subset( tps, select=c(TAXA1, ID1) )
	setnames(tmp, c('TAXA1','ID1'), c('TAXA2','ID2') )
	tmp				<- unique(rbind( tmp, subset( tps, select=c(TAXA2, ID2) ) ))
	setnames(tmp, c('TAXA2','ID2'), c('TAXA','ID') )
	setkey(tmp, ID)
	rownames(d)		<- tmp[, TAXA]
	colnames(d)		<- tmp[, TAXA]
	#	checks	
	stopifnot(ncol(d)==nrow(d))
	stopifnot(length(which(is.na(d)))==2*nrow(subset(tps, is.na(GD))))
	cat('D matrix: proportion of NA entries=',length(which(is.na(d))) / prod(dim(d)))
	#
	#	get variance matrix
	#
	tmp				<- dcast.data.table(tps, ID1~ID2, value.var='GD_V')		
	v				<- cbind(NA_real_, as.matrix(tmp[, -1, with=FALSE]))
	v				<- rbind(v, NA_real_)
	colnames(v)[1]	<- setdiff( as.character(tmp[, ID1]), colnames(v) )
	rownames(v)		<- colnames(v)
	diag(v)			<- 0
	#	complete lower triangular from upper triangular and vice versa
	tmp				<- lower.tri(v) & is.na(v)	
	v[tmp]			<- t(v)[tmp]
	tmp				<- upper.tri(v) & is.na(v)
	v[tmp]			<- t(v)[tmp]
	#	reset names
	tmp				<- subset( tps, select=c(TAXA1, ID1) )
	setnames(tmp, c('TAXA1','ID1'), c('TAXA2','ID2') )
	tmp				<- unique(rbind( tmp, subset( tps, select=c(TAXA2, ID2) ) ))
	setnames(tmp, c('TAXA2','ID2'), c('TAXA','ID') )
	setkey(tmp, ID)
	rownames(v)		<- tmp[, TAXA]
	colnames(v)		<- tmp[, TAXA]
	#	checks	
	stopifnot(ncol(v)==nrow(v))	
	cat('V matrix: proportion of NA entries=',length(which(is.na(v))) / prod(dim(v)))
	#
	#	remove cols/rows that contain nothing else than NAs
	#	
	diag(d)			<- NA_real_
	diag(v)			<- NA_real_		
	tmp				<- apply(d, 1, function(x) !all(is.na(x)))
	cat('\nIn D: found',length(which(!tmp)),'columns / rows with NA only: remove in D and V. ', rownames(d)[!tmp])
	ds				<- d[tmp,tmp]
	vs				<- v[tmp,tmp]	
	tmp				<- apply(vs, 1, function(x) !all(is.na(x)))
	cat('\nIn V: found additional',length(which(!tmp)),'columns / rows with NA only: remove in D and V too. ', rownames(d)[!tmp])
	ds				<- ds[tmp,tmp]
	vs				<- vs[tmp,tmp]	
	#
	#	remove cols/rows that contain more than 10% NAs
	#
	if(!is.na(na.rm.p))
	{		
		tmp				<- apply(ds, 1, function(x) length(which(is.na(x))) ) / ncol(ds)
		tmp				<- which(tmp>na.rm.p)
		cat('\nIn D: found',length(tmp),'columns / rows with more than',na.rm.p*100,'% NAs: remove in D and V. ', rownames(ds)[tmp])
		ds				<- ds[-tmp,-tmp]
		vs				<- vs[-tmp,-tmp]			
	}
	#
	#	do MDS to impute missing distances since mvrs fails too often
	#	
	attach(mds.args)
	diag(ds)		<- 0
	#mds.fit		<- mds(ds, ...)
	mds.fit			<- mds(ds, ndim=ndim, type=type, spline.intKnots=spline.intKnots, spline.degree=spline.degree)
	if(!is.na(wfile))
	{
		pdf(file=paste(wfile,'.shephard.pdf',sep=''), width=5, height=5)
		plot(mds.fit, plot.type = "Shepard")	
		dev.off()
	}	
	mds.ifit 		<- inverseMDS(mds.fit$conf) 
	if(!is.na(wfile))
	{
		save(tps, ds, vs, mds.fit, mds.ifit, mds.args, file=paste(wfile,'mds.rda',sep=''))	
	}
	#	unfinished
	
	#rnd.stress 		<- mean( randomstress(n=1591, ndim=500, nrep=100) )	#even just one iteration takes forever
	# njs(ds, fs = 15) # runs fine
	# ph			<- mvrs(ds, vs, fs = 15)			
}
##--------------------------------------------------------------------------------------------------------
##	olli 27.06.11
treecomparison.bootstrap.sd.vs.coverage<- function(indir=NULL, wdir=NULL)
{
	require(ape)
	require(data.table)
	require(ggplot2)
	
	batch.n	<- 3200
	if(is.null(indir) | is.null(wdir))
	{
		indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/simulations'
		wdir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/tree_mvr'		
	}
	
	infile	<- '150701_Regional_TRAIN4_SIMULATED.fa'		
	seq		<- read.dna(file.path(indir, infile), format='fa')
	seqi	<- as.data.table(read.csv(file.path(indir, gsub('\\.fa','_gene.txt',infile)), header=0))
	seqi[, GENE:= regmatches(V2, regexpr('[a-z]+', V2))]
	seqi[, START:= as.integer(gsub('-','',regmatches(V2, regexpr('[0-9]+-', V2))))]
	seqi[, END:= as.integer(gsub('-','',regmatches(V2, regexpr('-[0-9]+', V2))))]
	seqi[, GENE_L:=END-START+1L]
	seqi	<- subset(seqi, select=c(GENE,START,END,GENE_L))
	seqi	<- rbind(seqi, data.table(GENE='gag+pol+env', START=1L, END=seqi[, max(END)], GENE_L=seqi[, max(END)]))
	
	
	load(file.path(wdir,'150701_Regional_TRAIN4_SIMULATED_tps.rda'))
	
	tp		<- merge(tp, subset(seqi, select=c(GENE, GENE_L)), by='GENE')
	set(tp, NULL, 'GENE', tp[, factor(GENE, levels=c('gag','pol','env','gag+pol+env'), labels=c('gag','pol','env','gag+pol+env'))])
	#
	#	do we have higher bootstrap variance if there are more gaps by gene?
	#
	tmp		<- subset(tp, REP==1)
		
	ggplot( subset(tmp, GD_MEAN>0), aes(x=cut(DO/GENE_L, breaks=seq(0,1,0.01), labels=seq(0.01,1,0.01)), y=GD_SD) ) + 	
			geom_boxplot(outlier.shape=NA) +
			coord_cartesian(ylim=c(0,0.18)) +
			scale_x_discrete(breaks=seq(0,1,0.1), labels=paste(100*seq(0,1,0.1),'%',sep='')) +
			scale_y_continuous(expand=c(0,0)) +
			facet_grid(~GENE) + theme_bw() +
			labs(x='\noverlap between taxon pairs\n(% of sequence length)', y='std deviation in genetic distance\n')
	ggsave(file=file.path(wdir, gsub('.fa','_GDSD_by_overlap.pdf',infile)), w=14, h=7)
	
	ggplot( subset(tmp, GD_MEAN>0), aes(x=cut(DO/GENE_L, breaks=seq(0,1,0.01), labels=seq(0.01,1,0.01)), y=GD_SD/GD_MEAN) ) + 	
			geom_boxplot(outlier.shape=NA) +
			coord_cartesian(ylim=c(0,0.6)) +
			scale_x_discrete(breaks=seq(0,1,0.1), labels=paste(100*seq(0,1,0.1),'%',sep='')) +
			scale_y_continuous(expand=c(0,0)) +
			facet_grid(~GENE) + theme_bw() +
			labs(x='\noverlap between taxon pairs\n(% of sequence length)', y='coefficient of variation in genetic distance\nacross bootstrap alignments\n')
	ggsave(file=file.path(wdir, gsub('.fa','_GDCOV_by_overlap.pdf',infile)), w=14, h=7)
}
##--------------------------------------------------------------------------------------------------------
##	olli 27.06.11
treecomparison.gd.dev<- function()
{
	require(ape)
	require(data.table)	
	#
	#	get true trees
	#
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim_internal/freeze_July15'
	tfiles	<- list.files(indir, pattern='newick$', full.names=TRUE)
	tfiles	<- data.table( FILE_T= tfiles[ grepl('SUBSTTREE', tfiles) | grepl('Vill_99', tfiles) | grepl('Vill.*DATEDTREE', tfiles) ] )
	tfiles[, SC:= toupper(gsub('_SUBSTTREE|_DATEDTREE','',gsub('.newick','',basename(FILE_T))))]
	tmp		<- rbind( subset(tfiles, SC=='VILL_99_APR15'), subset(tfiles, SC=='VILL_99_APR15'), subset(tfiles, SC=='VILL_99_APR15') )
	set(tmp, NULL, 'SC', c('150701_VILL_SCENARIO-C','150701_VILL_SCENARIO-D','150701_VILL_SCENARIO-E'))
	tfiles	<- rbind(tfiles, tmp)
	tmp		<- list.files(indir, pattern='newick$', full.names=TRUE)
	tmp		<- data.table( FILE_T= tmp[ grepl('Reg.*DATEDTREE', tmp) ] )
	tmp[, SC:= toupper(gsub('_SUBSTTREE|_DATEDTREE','',gsub('.newick','',basename(FILE_T))))]
	tfiles	<- rbind(tfiles, tmp)
	tfiles[, BRL_T:= 'time']	
	set(tfiles, tfiles[, which(grepl('REG',SC) & grepl('SUBST',FILE_T))], 'BRL_T', 'subst')
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/simulations'
	tmp		<- data.table(FASTA_FILE=list.files(indir, full.names=1, pattern='_TRAIN[0-9]+_SIMULATED.fa$'))	
	tmp[, SC:= toupper(gsub('_SIMULATED.fa','',basename(FASTA_FILE)))]
	tfiles	<- merge(tfiles, tmp, by='SC')
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/tree_mvr'
	tmp		<- data.table(MVR_FILE=list.files(indir, full.names=1, pattern='_SIMULATED_tps.rda$'))	
	tmp[, SC:= toupper(gsub('_SIMULATED_tps.rda','',basename(MVR_FILE)))]
	tfiles	<- merge(tfiles, tmp, by='SC', all.x=1)
	ttrs	<- lapply(tfiles[, FILE_T], function(x)	read.tree(file=x) )
	names(ttrs)	<- tfiles[, SC]	
	tfiles[, IDX_T:=seq_along(ttrs)]
	tfiles[, TAXAN_T:= sapply(ttrs, Ntip)]
	#
	#	read true patristic distances from true tree
	#
	tbrl	<- subset(tfiles, BRL_T=='subst')[, {
				#IDX_T	<- 1
				ph		<- ttrs[[IDX_T]]
				tmp		<- distTips(ph, seq_len(Ntip(ph)), method='patristic', useC=TRUE)
				tmp		<- as.matrix(tmp)
				tmp[upper.tri(tmp, diag=TRUE)]	<- NA_real_
				tmp		<- as.data.table(melt(tmp))								
				setnames(tmp, c('Var1','Var2','value'),c('TAXA1','TAXA2','PD_T'))
				tmp		<- subset(tmp, !is.na(PD_T))
				tmp
			}, by='IDX_T']
	#	z[, table(IDX_T)]
	#	1       3       5       7       9 
	#	1279200 1279200 1279200 1279200 1279200
	#
	#	read true raw genetic distances from sequences
	#	
	tmp		<- subset(tfiles, BRL_T=='subst')[, {
				#FASTA_FILE<- "/Users/Oliver/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/simulations/150701_Regional_TRAIN2_SIMULATED.fa"
				sq		<- read.dna(FASTA_FILE, format='fa')
				tmp		<- dist.dna(sq, model='raw', as.matrix=TRUE, pairwise.deletion=TRUE)
				tmp		<- as.matrix(tmp)
				tmp[upper.tri(tmp, diag=TRUE)]	<- NA_real_
				tmp		<- as.data.table(melt(tmp))								
				setnames(tmp, c('Var1','Var2','value'),c('TAXA1','TAXA2','ALL_GD_RAW_T'))
				tmp		<- subset(tmp, !is.na(ALL_GD_RAW_T))
				tmp				
			}, by='IDX_T']
	tbrl	<- merge(tbrl, tmp, by=c('IDX_T','TAXA1','TAXA2'), all.x=1)
	set(tbrl, NULL, 'TAXA1', tbrl[, as.character(TAXA1)])
	set(tbrl, NULL, 'TAXA2', tbrl[, as.character(TAXA2)])
	#
	#	read genetic distances that I calculated previously
	#
	tmp		<- unique(subset(tfiles, BRL_T=='subst' & !is.na(MVR_FILE), c(IDX_T, MVR_FILE)))
	tmp		<- lapply(seq_len(nrow(tmp)), function(i){
				#MVR_FILE<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/tree_mvr/150701_Regional_TRAIN2_SIMULATED_tps.rda'
				MVR_FILE	<- tmp[i, MVR_FILE]
				load(MVR_FILE)
				z		<- subset(tp, REP==1, select=c(TAXA1, TAXA2, GENE, GD, DO, GD_MEAN, GD_SD))				
				z[, IDX_T:= tmp[i, IDX_T]]
				z
			})
	tmp		<- do.call('rbind',tmp)
	tmp		<- dcast.data.table(melt(tmp, id.vars=c('IDX_T','TAXA1','TAXA2','GENE')), IDX_T+TAXA1+TAXA2~variable+GENE, value.var='value')
	merge(tbrl, tmp, by=c('IDX_T','TAXA1','TAXA2'))
	
	subset(tbrl, IDX_T==3 & TAXA2=='IDPOP_100181|F|DOB_2000.1|2016.31')
	subset(tbrl, IDX_T==3 & grepl('IDPOP_100181',TAXA2))
	
	tmp[, table(IDX_T)]
	
	tbrl	<- merge(tbrl, tmp, by=c('IDX_T','TAXA1','TAXA2'), all=1)
	tbrl	<- subset(tbrl, !is.na(PD_T))
	#
	#	calculate MSEs
	#
	#	MSE from simple raw genetic distance approach
	tmp		<- subset(tbrl, !is.na(ALL_GD_RAW_T) & (is.na(GENE) | GENE=='gag+pol+env'))
	mse		<- tmp[, list( TYPE='ALL_GD_RAW_T', GENE='gag+pol+env', PAIR_N=length(PD_T), MSE=mean((PD_T-ALL_GD_RAW_T)*(PD_T-ALL_GD_RAW_T)) ), by='IDX_T']	
	#	MSE from my previously calculated distances	
	tmp		<- subset(tbrl, !is.na(GENE) & !is.na(GD))
	tmp		<- tmp[, list(		TYPE='GD', PAIR_N=length(PD_T), MSE=mean((PD_T-GD)*(PD_T-GD))	), 	by=c('IDX_T','GENE')]	
	mse		<- rbind(mse, tmp, use.names=TRUE, fill=TRUE)
	
			MSE_GD_MEAN=mean((PD_T-GD_MEAN)*(PD_T-GD_MEAN))
			
	tmp		<- tmp[, list(FULL_GAPS_P=mean(FULL_GAPS_P), GAG_GAPS_P=mean(GAG_GAPS_P), POL_GAPS_P=mean(POL_GAPS_P), ENV_GAPS_P=mean(ENV_GAPS_P)), by='SC']
	tinfo	<- merge(tinfo, tmp, all.x=1, by='SC')	

}

##--------------------------------------------------------------------------------------------------------
##	olli 27.06.11
treecomparison.bootstrap.gd.dev<- function()
{
	require(ape)
	require(data.table)
	
	bsn		<- 1e2
	repn	<- 10
	batch.n	<- 3200
	batch.i	<- 1
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/simulations'	
	infile	<- '150701_Regional_TRAIN4_SIMULATED.fa'	
	outdir	<- indir
	outfile	<- paste(gsub('.fa','',infile),'_GDS_BATCH',batch.i,'.rda',sep='')
	
	seq		<- read.dna(file.path(indir, infile), format='fa')
	#	for this first analysis, to see if the resulting trees are meaningful at all, 
	#	use just gag, pol and full
	seqi	<- as.data.table(read.csv(file.path(indir, gsub('\\.fa','_gene.txt',infile)), header=0))
	seqi[, GENE:= regmatches(V2, regexpr('[a-z]+', V2))]
	seqi[, START:= as.integer(gsub('-','',regmatches(V2, regexpr('[0-9]+-', V2))))]
	seqi[, END:= as.integer(gsub('-','',regmatches(V2, regexpr('-[0-9]+', V2))))]
	seqi	<- subset(seqi, select=c(GENE,START,END))
	#seqi	<- rbind(seqi, data.table(GENE='full', START=1L, END= seqi[, max(END)]))
		
	
	tp		<- as.data.table( t(combn(rownames(seq),2)) )	
	#tp		<- as.data.table( t(combn(sample(50, 10),2)) )	#this won t work. we need the specific ordering of the IDs that is used to create the combinations
	setnames(tp, c('V1','V2'), c('TAXA1','TAXA2'))
	tmp		<- as.data.table( t(combn(seq_len(nrow(seq)),2)) )
	setnames(tmp, c('V1','V2'), c('ID1','ID2'))
	tp		<- cbind(tp, tmp)
	tp[, IDX:= seq_len(nrow(tp))]
	#	subset by batch
	tp[, BATCH:= ceiling(IDX/batch.n)]
	if(!is.na(batch.i))
		tp	<- subset(tp, BATCH==batch.i)
	
	tmp				<- dcast.data.table(tp, TAXA1~TAXA2, value.var='IDX')
	d				<- cbind(NA_real_, as.matrix(tmp[, -1, with=FALSE]))
	d				<- rbind(d, NA_real_)
	colnames(d)[1]	<- setdiff( as.character(tmp[, TAXA1]), colnames(d) )	
	rownames(d)		<- colnames(d)
	diag(d)			<- 0
	#	complete lower triangular from upper triangular and vice versa
	tmp		<- lower.tri(d) & is.na(d)	
	d[tmp]	<- t(d)[tmp]
	tmp		<- upper.tri(d) & is.na(d)
	d[tmp]	<- t(d)[tmp]
	#	
	
	subset(tps, is.na(GD))
	
	stopifnot(ncol(d)==nrow(d))
	cat('proportion of NA entries=',length(which(is.na(d))) / prod(dim(d)))
	
	
	#for each pair, estimate: actual distance, mean distance, variance in distance:
	#	(do this pairwise because otherwise too computationally expensive
	#	by gene
	tpi	<- tp[, 	{
				cat('IDX',IDX, round(IDX/nrow(tp),d=3))
				#TAXA1	<- 'IDPOP_13649|M|DOB_1906.66|2011.23'; TAXA2<- 'IDPOP_27993|F|DOB_1961.29|1991.587'
				#START	<- 1; END<- 1473
				#system.time({
				df.gd	<- seqi[, {
										spc		<- as.character(seq[c(TAXA1,TAXA2), START:END])
										#	use same seed across all bootstrap runs, ie running for every gene is the same as running for the full genome
										#		and running for every taxon pair is the same as running for the whole alignment
										set.seed(42)
										tmp		<- as.data.table(expand.grid(REP=seq_len(repn), BS=seq_len(bsn+1)))
										tmp		<- tmp[, {
													#	take bootstrap sample (except if bsi==1)
													#	the bootstrap is relative to the gene region!											
													spcb	<- copy(spc)
													if(BS>1)
													{
														#	note: bootstrap includes ? columns, which adds uncertainty when genetic distances are evaluated over the overlap columns
														spcb<- spcb[, sample(ncol(spcb), replace=TRUE)]	
													}
													spb		<- as.DNAbin(spcb)
													#	overlap that is not '?'
													do		<- sum(apply( spcb!='?', 2, prod))
													#	count genetic distance on overlap region, ie count gaps '-' as well, on anything that is not '?'
													dn		<- as.numeric( dist.dna(spb, model='N', pairwise.deletion=TRUE) )
													#	add indels to differences, but not when the other sequence is '?'
													tmp		<- which( !apply(spcb=='?', 2, any) )
													if(length(tmp))
														dn	<- dn + as.numeric(dist.dna(spb[, tmp], model='indel'))
													#	DN can be > 0 if DO is 0, because of indels
													list(DN=dn, DO=do)	
												}, by=c('REP','BS')]
										tmp												
									}, by='GENE']
				#	collect distances for genes
				ans		<- subset(df.gd, BS>1)[, list( GD_MEAN=mean(DN/DO), GD_SD=sd(DN/DO) ), by=c('GENE','REP')]				
				tmp		<- subset(df.gd, BS==1)[, list( GD=ifelse(DO==0, NA_real_, DN/DO), DO=DO ), by=c('GENE','REP')]
				ans		<- merge(tmp, ans, by=c('GENE','REP'))
				#	calculate distances for full genome
				tmp		<- subset(df.gd, GENE%in%c('gag','pol','env'))[, list(GENE='gag+pol+env', GD=sum(DN)/sum(DO), DO=sum(DO)), by=c('REP','BS')]
				tmp		<- tmp[, list(GD= GD[BS==1], DO=DO[BS==1], GD_MEAN=mean(GD[BS>1]), GD_SD=sd(GD[BS>1])), by=c('GENE','REP')]
				ans		<- rbind(ans, tmp)
				#})
				ans				
			}, by=c('TAXA1','TAXA2')]
	#	save output to
	save(tpi, file=file.path(outdir,outfile))
}
##--------------------------------------------------------------------------------------------------------
##	olli 27.06.11
treecomparison.create.metadata<- function()
{	
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim_internal/freeze_July15'
	outdir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/simulations'
	infile.prefix	<- '150701_Regional_TRAIN1_'
	
	load( file.path(indir, paste(infile.prefix,'SIMULATED_INTERNAL.R',sep='')) )
	
	tmp			<- subset( df.inds, !is.na(TIME_SEQ), select=c(IDPOP, GENDER, DOB, DOD, DIAG_T, DIAG_CD4, ART1_T, ART1_CD4, TIME_SEQ, RECENT_TR ) )	
	set(tmp, NULL, 'GENDER', tmp[,as.character(GENDER)])
	tmp2		<- tmp[, which(is.na(DIAG_T) & TIME_SEQ<2000)]
	cat(paste('\nSet patient variables to NA for archival seq, n=',length(tmp2)))		
	set(tmp, tmp2, c('DOB','DOD'), NA_real_) 
	set(tmp, tmp2, 'GENDER', NA_character_)
	tmp2		<- tmp[, which(is.na(DIAG_T) & TIME_SEQ>=2000)]
	cat(paste('\nSet patient variables to NA after 2000, n=',length(tmp2)))
	print(tmp[tmp2,])
	set(tmp, tmp2, c('DOB','DOD'), NA_real_) 
	set(tmp, tmp2, 'GENDER', NA_character_)
	set(tmp, NULL, 'GENDER', tmp[,factor(GENDER)])
	
	file			<- paste(outdir, '/', infile.prefix, 'SIMULATED_metadata.csv', sep='')
	cat(paste('\nwrite to file', file))
	write.csv(tmp, file)


}
##--------------------------------------------------------------------------------------------------------
##
##--------------------------------------------------------------------------------------------------------
treecomparison.ana.160627.standardize.MSE<- function()
{	
	require(ggplot2)
	require(data.table)
	require(ape)
	require(scales)	
	require(ggtree)
	require(phangorn)
	
	edir			<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation'
	timetag			<- '160627'
	load(file.path(edir,'submitted_160713_07MSELSD.rda'))
	
	sc		<- copy(sclu.info)
	#
	tmp		<- subset(tinfo, !is.na(IDCLU))[, list(CLU_N=CLU_N[1], MXGPS_CLU= max(GPS), MDGPS_CLU=median(GPS)), by=c('SC','IDCLU')]
	sc		<- merge(sc, tmp, by=c('SC','IDCLU'))	
	set(sc, NULL, 'MODEL', sc[, factor(MODEL, levels=c('V','R'),labels=c('Model: Village','Model: Regional'))])
	set(sc, sc[, which(SC=="VILL_99_APR15")],'SC',"150701_VILL_SCENARIO-C")	
	set(sc, NULL, 'SC', sc[, factor(SC,	levels=c("150701_REGIONAL_TRAIN1", "150701_REGIONAL_TRAIN2", "150701_REGIONAL_TRAIN3", "150701_REGIONAL_TRAIN4","150701_REGIONAL_TRAIN5","150701_VILL_SCENARIO-A","150701_VILL_SCENARIO-B","150701_VILL_SCENARIO-C","150701_VILL_SCENARIO-D","150701_VILL_SCENARIO-E"), 
							labels=c('sc 1','sc 2','sc 3','sc 4','sc 5','sc A','sc B','sc C','sc D','sc E'))])
	set(sc, NULL, 'GAPS', sc[, factor(GAPS, levels=c('none','low','high'),labels=c('none','as for Botswana\nsequences','as for Uganda\nsequences'))])
	set(sc, NULL, 'BEST', sc[, factor(BEST, levels=c('Y','N'),labels=c('best tree','replicate tree'))])									
	set(sc, NULL, 'GENE', sc[, factor(GENE, levels=c('GAG','POL','GAG+POL+ENV'),labels=c('gag','pol','gag+pol+env'))])	
	set(sc, NULL, 'TEAM', sc[, factor(TEAM, levels=sc[, sort(unique(TEAM))],labels=sc[, sort(unique(TEAM))])])
	set(sc, NULL, 'EXT', sc[, factor(EXT, levels=c('~0pc','5pc'),labels=c('~ 0%/year','5%/year'))])
	set(sc, NULL, 'ART', sc[, factor(ART, levels=c('none','fast'),labels=c('none','fast'))])
	sc		<- subset(sc, OTHER=='N')	
	
	
	require(gamlss)
	ggplot(subset(sc, TEAM!='MetaPIGA' & CLU_N<100), aes(x=CLU_N, y=MSE, colour=GENE, pch=TEAM)) + geom_point() + facet_grid(~SC)
	ggplot(subset(sc, TEAM!='MetaPIGA' & SC=='sc 1'), aes(x=CLU_N, y=MSE, colour=GENE, pch=TEAM)) + geom_point()
	ggplot(subset(sc, TEAM!='MetaPIGA'), aes(x=CLU_N, y=MSE, colour=GENE, pch=TEAM)) + geom_point() + coord_cartesian(ylim=c(0,1e4), xlim=c(0,100)) + facet_grid(~SC)
	
	
	#
	#	look reasonable to divive KC by CLU_N*(CLU_N-1)/2
	#
	kc.std.d	<- subset(sc, TEAM!='MetaPIGA' & SC=='sc 1')
	kc.std.m1	<- gamlss(KC~CLU_N-1, data=kc.std.d)
	kc.std.m2	<- gamlss(KC~poly(CLU_N,2, raw=TRUE), data=kc.std.d)	#this allows for a non-zero baseline, which gave much better fit	
	#gamlss(KC~CLU_N+I(CLU_N^2)-1, data=kc.std.d)
	kc.std.m3	<- gamlss(KC~I(CLU_N*(CLU_N-1)/2), data=kc.std.d)
	kc.std.m4	<- gamlss(KC~poly(CLU_N,4, raw=TRUE), data=kc.std.d)
	#kc.std.m4	<- gamlss(KC~I(sqrt(CLU_N*(CLU_N-1)/2))-1, data=kc.std.d)
	kc.std.da	<- subset(sc, TEAM!='MetaPIGA' & SC%in%c('sc 1','sc 2','sc 4'))
	tmp.m1		<- predict(kc.std.m1, data=kc.std.d, newdata=kc.std.da, type='response', se.fit=FALSE)
	tmp.m2		<- predict(kc.std.m2, data=kc.std.d, newdata=kc.std.da,type='response', se.fit=FALSE)
	tmp.m3		<- predict(kc.std.m3, data=kc.std.d, newdata=kc.std.da,type='response', se.fit=FALSE)
	tmp.m4		<- predict(kc.std.m4, data=kc.std.d, newdata=kc.std.da,type='response', se.fit=FALSE)
	kc.std.da[, KC.m1:=tmp.m1]
	kc.std.da[, KC.m2:=tmp.m2]
	kc.std.da[, KC.m3:=tmp.m3]
	kc.std.da[, KC.m4:=tmp.m4]
	kc.std.da	<- melt(kc.std.da, measure.vars=c('KC.m1','KC.m2','KC.m3','KC.m4'))
	set(kc.std.da, NULL, 'variable', kc.std.da[, factor(variable, levels=c('KC.m1','KC.m2','KC.m3','KC.m4'), labels=c('KC~CLU_N-1','KC~poly(CLU_N,2, raw=TRUE)','KC~I(CLU_N*(CLU_N-1)/2)','KC~poly(CLU_N,4, raw=TRUE)'))])
	ggplot(kc.std.da, aes(x=CLU_N)) + geom_point(aes(y=KC,colour=GENE, pch=TEAM)) + 
			geom_line(aes(y=value, linetype=variable, group=variable)) +
			#scale_linetype_manual(values=c('KC~CLU_N-1'='a','KC~poly(CLU_N,2, raw=TRUE)-1'='e','KC~I(CLU_N*(CLU_N-1)/2)-1'='f','KC~poly(CLU_N,4, raw=TRUE)-1'='j')) +
			facet_grid(~SC)
	ggsave(file.path(edir, paste(timetag,'_','dependence_KC_clustersize.pdf',sep='')), w=15, h=7)
}
##--------------------------------------------------------------------------------------------------------
##
##--------------------------------------------------------------------------------------------------------
treecomparison.ana.160627.standardize.KC<- function()
{	
	require(ggplot2)
	require(data.table)
	require(ape)
	require(scales)	
	require(ggtree)
	require(phangorn)
	
	edir			<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation'
	timetag			<- '160713'
	load(paste(edir,'/','submitted_160713_07MSELSD.rda',sep=''))
	
	
	sc		<- copy(sclu.info)
	#
	tmp		<- subset(tinfo, !is.na(IDCLU))[, list(CLU_N=CLU_N[1], MXGPS_CLU= max(GPS), MDGPS_CLU=median(GPS)), by=c('SC','IDCLU')]
	sc		<- merge(sc, tmp, by=c('SC','IDCLU'))	
	set(sc, NULL, 'MODEL', sc[, factor(MODEL, levels=c('V','R'),labels=c('Model: Village','Model: Regional'))])
	set(sc, sc[, which(SC=="VILL_99_APR15")],'SC',"150701_VILL_SCENARIO-C")	
	set(sc, NULL, 'SC', sc[, factor(SC,	levels=c("150701_REGIONAL_TRAIN1", "150701_REGIONAL_TRAIN2", "150701_REGIONAL_TRAIN3", "150701_REGIONAL_TRAIN4","150701_REGIONAL_TRAIN5","150701_VILL_SCENARIO-A","150701_VILL_SCENARIO-B","150701_VILL_SCENARIO-C","150701_VILL_SCENARIO-D","150701_VILL_SCENARIO-E"), 
							labels=c('sc 1','sc 2','sc 3','sc 4','sc 5','sc A','sc B','sc C','sc D','sc E'))])
	set(sc, NULL, 'GAPS', sc[, factor(GAPS, levels=c('none','low','high'),labels=c('none','as for Botswana\nsequences','as for Uganda\nsequences'))])
	set(sc, NULL, 'BEST', sc[, factor(BEST, levels=c('Y','N'),labels=c('best tree','replicate tree'))])									
	set(sc, NULL, 'GENE', sc[, factor(GENE, levels=c('GAG','POL','GAG+POL+ENV'),labels=c('gag','pol','gag+pol+env'))])	
	set(sc, NULL, 'TEAM', sc[, factor(TEAM, levels=sc[, sort(unique(TEAM))],labels=sc[, sort(unique(TEAM))])])
	set(sc, NULL, 'EXT', sc[, factor(EXT, levels=c('~0pc','5pc'),labels=c('~ 0%/year','5%/year'))])
	set(sc, NULL, 'ART', sc[, factor(ART, levels=c('none','fast'),labels=c('none','fast'))])
	sc		<- subset(sc, OTHER=='N')	
	
	
	require(gamlss)
	ggplot(subset(sc, TEAM!='MetaPIGA' & CLU_N<100), aes(x=CLU_N, y=KC, colour=GENE, pch=TEAM)) + geom_point() + facet_grid(~SC)
	ggplot(subset(sc, TEAM!='MetaPIGA' & SC=='sc 1'), aes(x=CLU_N, y=KC, colour=GENE, pch=TEAM)) + geom_point()
	
	#
	#	look reasonable to divive KC by CLU_N*(CLU_N-1)/2
	#
	kc.std.d	<- subset(sc, TEAM!='MetaPIGA' & SC=='sc 1', select=c(SC,TEAM,GENE,CLU_N, KC))
	kc.std.m1	<- gamlss(KC~CLU_N-1, data=kc.std.d)
	kc.std.m2	<- gamlss(KC~poly(CLU_N,2, raw=TRUE), data=kc.std.d)	#this allows for a non-zero baseline, which gave much better fit	
	#gamlss(KC~CLU_N+I(CLU_N^2)-1, data=kc.std.d)
	kc.std.m3	<- gamlss(KC~I(CLU_N*(CLU_N-1)/2), data=kc.std.d)
	kc.std.m4	<- gamlss(KC~poly(CLU_N,4, raw=TRUE), data=kc.std.d)
	#kc.std.m4	<- gamlss(KC~I(sqrt(CLU_N*(CLU_N-1)/2))-1, data=kc.std.d)
	kc.std.da	<- subset(sc, !TEAM%in%c('MetaPIGA','MVR','BioNJ') & SC%in%c('sc 1','sc 2','sc 4'), select=c(SC,TEAM,GENE,CLU_N, KC))
	tmp.m1		<- predict(kc.std.m1, data=kc.std.d, newdata=kc.std.da, type='response', se.fit=FALSE)
	tmp.m2		<- predict(kc.std.m2, data=kc.std.d, newdata=kc.std.da,type='response', se.fit=FALSE)
	tmp.m3		<- predict(kc.std.m3, data=kc.std.d, newdata=kc.std.da,type='response', se.fit=FALSE)
	tmp.m4		<- predict(kc.std.m4, data=kc.std.d, newdata=kc.std.da,type='response', se.fit=FALSE)
	kc.std.da[, KC.m1:=tmp.m1]
	kc.std.da[, KC.m2:=tmp.m2]
	kc.std.da[, KC.m3:=tmp.m3]
	kc.std.da[, KC.m4:=tmp.m4]
	kc.std.da	<- melt(kc.std.da, measure.vars=c('KC.m1','KC.m2','KC.m3','KC.m4'))
	set(kc.std.da, NULL, 'variable', kc.std.da[, factor(variable, levels=c('KC.m1','KC.m2','KC.m3','KC.m4'), labels=c('KC~CLU_N-1','KC~poly(CLU_N,2, raw=TRUE)','KC~I(CLU_N*(CLU_N-1)/2)','KC~poly(CLU_N,4, raw=TRUE)'))])
	ggplot(kc.std.da, aes(x=CLU_N)) + geom_point(aes(y=KC,colour=GENE, pch=TEAM)) + 
			geom_line(aes(y=value, linetype=variable, group=variable)) +
			#scale_linetype_manual(values=c('KC~CLU_N-1'='a','KC~poly(CLU_N,2, raw=TRUE)-1'='e','KC~I(CLU_N*(CLU_N-1)/2)-1'='f','KC~poly(CLU_N,4, raw=TRUE)-1'='j')) +
			facet_grid(~SC)
	ggsave(file.path(edir, paste(timetag,'_','dependence_KC_clustersize.pdf',sep='')), w=15, h=7)
}
##--------------------------------------------------------------------------------------------------------
##	olli 27.06.16
##--------------------------------------------------------------------------------------------------------
treecomparison.submissions.160713<- function()	
{
	load(paste(edir,'/','submitted_160713_QD.rda',sep=''))
	tmp		<- subset(submitted.info, select=c(IDX, TAXA_NJ, NQD))
	tmp2	<- subset(sclu.info, select=c(IDX, IDCLU, NQDC))
	load(paste(edir,'/','submitted_160713.rda',sep=''))
	submitted.info	<- merge(submitted.info, tmp, by='IDX')
	sclu.info		<- merge(sclu.info, tmp2, by=c('IDX','IDCLU'))
	
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/simulations'
	infiles	<- data.table(FILE=list.files(indir, full.names=1, pattern='_TRAIN[0-9]+_SIMULATED.fa$'))
	infiles[, SC:= toupper(gsub('_SIMULATED.fa','',basename(FILE)))]
	tmp		<- infiles[, {
				#FILE<- "/Users/Oliver/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/simulations/150701_Regional_TRAIN4_SIMULATED.fa"
				sq	<- read.dna(FILE, format='fa')
				seqi<- as.data.table(read.csv(gsub('.fa','_gene.txt',FILE), header=FALSE))
				seqi[, GENE:= regmatches(V2, regexpr('[a-z]+', V2))]
				seqi[, START:= as.integer(gsub('-','',regmatches(V2, regexpr('[0-9]+-', V2))))]
				seqi[, END:= as.integer(gsub('-','',regmatches(V2, regexpr('-[0-9]+', V2))))]
				seqi	<- subset(seqi, select=c(GENE,START,END))
				seqi	<- rbind(seqi, data.table(GENE='full', START=1L, END= seqi[, max(END)]))				
				ans		<- seqi[, {
							#START<- 1474; END<- 5466
							z	<- as.character(sq[,START:END])
							tmp	<- apply(z, 2, function(x) all(x%in%c('-','?'))) 							
							z	<- z[, !tmp]
							tmp	<- apply(z, 1, function(x) length(which(x=='?'))) / ncol(z)
							list(TAXA=names(tmp), GAPS_P=tmp)
						}, by='GENE']
				set(ans, NULL, 'GENE', ans[, paste(toupper(GENE),'_GAPS_P',sep='')])
				ans		<- dcast.data.table(ans, TAXA~GENE, value.var='GAPS_P')
				ans
			}, by='SC']
	tmp[, list(FULL_GAPS_P=mean(FULL_GAPS_P), GAG_GAPS_P=mean(GAG_GAPS_P), POL_GAPS_P=mean(POL_GAPS_P), ENV_GAPS_P=mean(ENV_GAPS_P)), by='SC']
	
	tinfo	<- merge(tinfo, tmp, by=c('SC','TAXA'), all.x=1)
	
	save(strs, strs_rtt, ttrs, tinfo, tfiles, tinfo.pairs, submitted.info, sclu.info, lba, file=file.path(edir, 'submitted_160713_RFPDQDTP.rda'))
}
##--------------------------------------------------------------------------------------------------------
##	olli 27.06.16
##--------------------------------------------------------------------------------------------------------
treecomparison.submissions.161123<- function()	
{
	require(data.table)
	require(ape)
	require(adephylo)
	require(phangorn)
	require(parallel)
	#
	#	get true trees
	#
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/simulations_trees'
	tfiles	<- list.files(indir, pattern='newick$', full.names=TRUE)
	
	
	tfiles	<- data.table( FILE_T= tfiles[ grepl('SUBSTTREE', tfiles) | grepl('Vill_99', tfiles) | grepl('.*DATEDTREE', tfiles) ] )
	tfiles[, SC:= toupper(gsub('_SUBSTTREE|_DATEDTREE','',gsub('.newick','',basename(FILE_T))))]
	tmp		<- rbind( subset(tfiles, SC=='VILL_99_APR15'), subset(tfiles, SC=='VILL_99_APR15'), subset(tfiles, SC=='VILL_99_APR15') )
	set(tmp, NULL, 'SC', c('150701_VILL_SCENARIO-C','150701_VILL_SCENARIO-D','150701_VILL_SCENARIO-E'))
	tfiles	<- rbind(tfiles, tmp)
	#tmp		<- list.files(indir, pattern='newick$', full.names=TRUE)
	#tmp		<- data.table( FILE_T= tmp[ grepl('*DATEDTREE', tmp) ] )
	#tmp[, SC:= toupper(gsub('_SUBSTTREE|_DATEDTREE','',gsub('.newick','',basename(FILE_T))))]
	#tfiles	<- rbind(tfiles, tmp)
	tfiles[, BRL_T:= 'time']	
	set(tfiles, NULL, 'SC', tfiles[, gsub('161121_','161121_REGIONAL_',SC)])
	set(tfiles, tfiles[, which(grepl('REG',SC) & grepl('SUBST',FILE_T))], 'BRL_T', 'subst')	
	ttrs	<- lapply(tfiles[, FILE_T], function(x)	read.tree(file=x) )
	names(ttrs)	<- tfiles[, SC]	
	for(z in c('VILL_99_APR15','150701_VILL_SCENARIO-C','150701_VILL_SCENARIO-D','150701_VILL_SCENARIO-E'))
		ttrs[[z]]	<- root(ttrs[[z]], node=Ntip(ttrs[[z]])+2, resolve.root=1)	
	tfiles[, IDX_T:=seq_along(ttrs)]
	tfiles[, TAXAN_T:= sapply(ttrs, Ntip)]
	#	patristic distances on true trees (by time and subst/site)
	tbrl	<- lapply(seq_len(nrow(tfiles)), function(i)
			{
				ph		<- ttrs[[tfiles[i, IDX_T]]]
				tmp		<- cophenetic.phylo(ph)				
				tmp		<- as.data.table(melt(tmp))
				setnames(tmp, c('Var1','Var2','value'),c('TAXA1','TAXA2','PD'))
				tmp		<- subset(tmp, TAXA1!=TAXA2)
				tmp[, IDX_T:= tfiles[i, IDX_T]]
				tmp[, SC:= tfiles[i, SC]]
				tmp[, BRL_T:= tfiles[i, BRL_T]]
				tmp[, TAXAN_T:= tfiles[i, TAXAN_T]]
				tmp
			})
	tbrl	<- do.call('rbind',tbrl)		
	#	info on true trees
	tinfo	<- merge(tfiles, do.call('rbind',lapply(seq_along(ttrs), function(i) data.table(TAXA=ttrs[[i]]$tip.label, IDX_T=i))), by='IDX_T')	
	tinfo[, IDPOP:=NA_character_]
	tmp		<- tinfo[, which(grepl('REGIONAL',SC))]
	set(tinfo, tmp, 'IDPOP', tinfo[tmp,regmatches(TAXA, regexpr('IDPOP_[0-9]+',TAXA))])
	tmp		<- tinfo[, which(!grepl('REGIONAL',SC))]
	set(tinfo, tmp, 'IDPOP', tinfo[tmp, regmatches(TAXA, regexpr('HOUSE[0-9]+-[0-9]+|House[0-9]+-[0-9]+',TAXA))])		
	stopifnot(subset(tinfo, grepl('VILL',SC))[, length(which(substring(TAXA,1,10)!=substring(IDPOP,1,10)))]==0)	
	stopifnot( tinfo[, length(which(is.na(IDPOP)))==0] )	
	set(tinfo, NULL, 'IDPOP', tinfo[,toupper(IDPOP)])
	set(tinfo, NULL, 'TAXA', tinfo[,toupper(TAXA)])
	#	read cluster membership from DATEDCLUTREES	
	tmp		<- list.files(indir, pattern='DATEDCLUTREES', full.names=TRUE)
	tmp		<- rbind( data.table( 	FILE_CLU_T= tmp, 
					SC= gsub('161121_','161121_REGIONAL_',toupper(gsub('_DATEDCLUTREES','',gsub('.newick','',basename(tmp))))),
					BRL_T= 'time'),
			data.table( 	FILE_CLU_T= tmp, 
					SC= gsub('161121_','161121_REGIONAL_',toupper(gsub('_DATEDCLUTREES','',gsub('.newick','',basename(tmp))))),
					BRL_T= 'subst') )
	tfiles	<- merge(tfiles, tmp, by=c('SC','BRL_T'), all=1)	
	tmp		<- subset(tfiles, !is.na(FILE_CLU_T))[, {
				z		<- read.tree(FILE_CLU_T)
				do.call('rbind',lapply(seq_along(z), function(i) data.table(IDCLU=i, TAXA=z[[i]]$tip.label)))				
			}, by=c('SC','BRL_T')]	
	tinfo	<- merge(tinfo, tmp, by=c('SC','BRL_T','TAXA'), all=1)
	tmp		<- subset(tinfo, !is.na(IDCLU))[, list(CLU_N= length(IDPOP)), by=c('SC','BRL_T','IDCLU')]
	tinfo	<- merge(tinfo, tmp, by=c('SC','BRL_T','IDCLU'), all=1)
	#	read sequences and determine %gappiness
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/simulations'
	tmp		<- list.files(indir, pattern='fa$|fasta$', full.names=TRUE)
	tmp		<- data.table( FILE_SEQ_T= tmp, SC= gsub('_GAG|_FULL|_P17','',gsub('161121_','161121_REGIONAL_',toupper(gsub('_SIMULATED','',gsub('\\.fa|\\.fasta','',basename(tmp)))))))
	#z		<- subset(tmp, SC=='VILL_99_APR15')
	#set(z, NULL, 'SC', '150701_VILL_SCENARIO-C')
	#tmp		<- rbind( tmp, z )	
	#tfiles	<- merge(tfiles, tmp, by='SC', all=1)
	#tmp		<- subset(tfiles, !is.na(FILE_SEQ_T) & !grepl('99_Apr15', FILE_T))[, {
	#			cat('\n',FILE_SEQ_T)
	#			z		<- read.dna(FILE_SEQ_T, format='fasta')	
	#			ans		<- sapply(seq_len(nrow(z)), function(i) base.freq(z[i,], all=1))
	#			ans		<- apply(ans[c('n','-','?'),], 2, sum)
	#			list(TAXA=rownames(z), GPS=ans)				
	#		}, by=c('SC','BRL_T')]
	#tinfo	<- merge(tinfo, tmp, by=c('SC','BRL_T','TAXA'), all.x=1)
	#	add % gaps by gene for regional
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/simulations'
	infiles	<- data.table(FILE=list.files(indir, full.names=1, pattern='_TRAIN[0-9]+_SIMULATED.fa$|GTRFIXED.*_SIMULATED.fasta$'))
	infiles[, SC:= gsub('161121_','161121_REGIONAL_',toupper(gsub('_SIMULATED.fa|_SIMULATED.fasta','',basename(FILE))))]
	tmp		<- subset(infiles, grepl('TRAIN|FULL', SC))[, {
				#FILE<- "/Users/Oliver/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/simulations/150701_Regional_TRAIN4_SIMULATED.fa"
				sq	<- read.dna(FILE, format='fa')
				seqi<- as.data.table(read.csv(gsub('\\.fa|\\.fasta','_gene.txt',FILE), header=FALSE))
				seqi[, GENE:= regmatches(V2, regexpr('[a-z]+', V2))]
				seqi[, START:= as.integer(gsub('-','',regmatches(V2, regexpr('[0-9]+-', V2))))]
				seqi[, END:= as.integer(gsub('-','',regmatches(V2, regexpr('-[0-9]+', V2))))]
				seqi	<- subset(seqi, select=c(GENE,START,END))
				seqi	<- rbind(seqi, data.table(GENE='full', START=1L, END= seqi[, max(END)]))				
				ans		<- seqi[, {
							#START<- 1474; END<- 5466
							z	<- as.character(sq[,START:END])
							tmp	<- apply(z, 2, function(x) all(x%in%c('-','?'))) 							
							z	<- z[, !tmp]
							tmp	<- apply(z, 1, function(x) length(which(x=='?'))) / ncol(z)
							list(TAXA=names(tmp), GAPS_P=tmp)
						}, by='GENE']
				set(ans, NULL, 'GENE', ans[, paste(toupper(GENE),'_GAPS_P',sep='')])
				ans		<- dcast.data.table(ans, TAXA~GENE, value.var='GAPS_P')
				ans
			}, by='SC']
	z		<- subset(infiles, !grepl('TRAIN|FULL', SC))[, {
				#FILE<- "/Users/Oliver/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/simulations/161121_GTRFIXED1_GAG_SIMULATED.fasta"
				sq	<- read.dna(FILE, format='fa')
				z	<- as.character(sq)
				zz	<- apply(z, 2, function(x) all(x%in%c('-','?'))) 							
				z	<- z[, !zz]
				ans	<- apply(z, 1, function(x) length(which(x=='?'))) / ncol(z)
				ans	<- data.table(TAXA=names(ans), GAPS_P=ans, GENE=gsub('GTRFIXED[0-9]_','',regmatches(FILE,regexpr('GTRFIXED[0-9]_[A-Z0-9]+',FILE))))				
				ans
			}, by='SC']	
	z	<- dcast.data.table(z, SC+TAXA~GENE, value.var='GAPS_P')
	setnames(z,c('GAG','P17'),c('GAG_GAPS_P','P17_GAPS_P'))
	tmp	<- rbind(tmp, z, use.names=TRUE,fill=TRUE)
	tmp		<- tmp[, list(FULL_GAPS_P=mean(FULL_GAPS_P), GAG_GAPS_P=mean(GAG_GAPS_P), POL_GAPS_P=mean(POL_GAPS_P), ENV_GAPS_P=mean(ENV_GAPS_P), P17_GAPS_P=mean(P17_GAPS_P)), by='SC']
	tinfo	<- merge(tinfo, tmp, all.x=1, by='SC')	
	#
	# to tinfo add actual transmitters
	#
	# check TRAIN1
	load( '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/simulations_trees/150701_Regional_TRAIN1_SIMULATED_INTERNAL.R' )	
	ch				<- subset(tinfo, SC=='150701_REGIONAL_TRAIN1' & BRL_T=='time', TAXA)
	ch[, IDPOP:= as.integer(gsub('IDPOP_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',1)))]
	ch[, GENDER_CH:= sapply(strsplit(TAXA,'|',fixed=1),'[[',2)]
	ch[, DOB_CH:= as.numeric(gsub('DOB_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',3)))]
	ch[, TIME_SEQ_CH:= as.numeric(sapply(strsplit(TAXA,'|',fixed=1),'[[',4))]
	ch				<- merge(subset(df.inds, select=c(IDPOP, GENDER, DOB, TIME_SEQ)), ch, by='IDPOP')	
	stopifnot( ch[, all(abs(DOB-DOB_CH)<=2*.Machine$double.eps)], ch[, all(GENDER==GENDER_CH)], ch[, all(abs(TIME_SEQ-TIME_SEQ_CH)<=0.001)] )
	# OK :-) schedule adding IDPOP_T
	df.trms[, IDPOP:= df.trms[, paste('IDPOP_',IDREC,sep='')]]
	tinfo.add		<- subset(df.trms, select=c(IDPOP, IDTR))	
	df.trms[, IDPOP:= df.trms[, paste('IDPOP_',IDTR,sep='')]]
	tinfo.add		<- merge(tinfo.add, subset(df.trms, select=c(IDPOP, IDREC)), by='IDPOP', all=1)
	set(ch, NULL, 'IDPOP', ch[, paste('IDPOP_',IDPOP,sep='')])
	tinfo.add		<- merge(tinfo.add, subset(ch, select=IDPOP), by='IDPOP')
	tinfo.add[, SC:='150701_REGIONAL_TRAIN1']
	# check TRAIN2
	load( '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/simulations_trees/150701_Regional_TRAIN2_SIMULATED_INTERNAL.R' )	
	ch				<- subset(tinfo, SC=='150701_REGIONAL_TRAIN2' & BRL_T=='time', TAXA)
	ch[, IDPOP:= as.integer(gsub('IDPOP_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',1)))]
	ch[, GENDER_CH:= sapply(strsplit(TAXA,'|',fixed=1),'[[',2)]
	ch[, DOB_CH:= as.numeric(gsub('DOB_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',3)))]
	ch[, TIME_SEQ_CH:= as.numeric(sapply(strsplit(TAXA,'|',fixed=1),'[[',4))]
	ch				<- merge(subset(df.inds, select=c(IDPOP, GENDER, DOB, TIME_SEQ)), ch, by='IDPOP')	
	stopifnot( ch[, all(abs(DOB-DOB_CH)<=2*.Machine$double.eps)], ch[, all(GENDER==GENDER_CH)], ch[, all(abs(TIME_SEQ-TIME_SEQ_CH)<=0.001)] )
	# OK :-) schedule adding IDPOP_T
	df.trms[, IDPOP:= df.trms[, paste('IDPOP_',IDREC,sep='')]]
	tmp				<- subset(df.trms, select=c(IDPOP, IDTR))	
	df.trms[, IDPOP:= df.trms[, paste('IDPOP_',IDTR,sep='')]]
	tmp				<- merge(tmp, subset(df.trms, select=c(IDPOP, IDREC)), by='IDPOP', all=1)
	set(ch, NULL, 'IDPOP', ch[, paste('IDPOP_',IDPOP,sep='')])
	tmp		<- merge(tmp, subset(ch, select=IDPOP), by='IDPOP')
	tmp[, SC:='150701_REGIONAL_TRAIN2']
	tinfo.add		<- rbind(tinfo.add, tmp)
	# check TRAIN4
	ch				<- subset(tinfo, SC=='150701_REGIONAL_TRAIN4' & BRL_T=='time', TAXA)
	ch[, IDPOP:= as.integer(gsub('IDPOP_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',1)))]
	ch[, GENDER_CH:= sapply(strsplit(TAXA,'|',fixed=1),'[[',2)]
	ch[, DOB_CH:= as.numeric(gsub('DOB_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',3)))]
	ch[, TIME_SEQ_CH:= as.numeric(sapply(strsplit(TAXA,'|',fixed=1),'[[',4))]
	ch				<- merge(subset(df.inds, select=c(IDPOP, GENDER, DOB, TIME_SEQ)), ch, by='IDPOP')	
	stopifnot( ch[, all(abs(DOB-DOB_CH)<=2*.Machine$double.eps)], ch[, all(GENDER==GENDER_CH)], ch[, all(abs(TIME_SEQ-TIME_SEQ_CH)<=0.001)] )	
	tmp[, SC:='150701_REGIONAL_TRAIN4']
	tinfo.add		<- rbind(tinfo.add, tmp)
	# check TRAIN3
	load( '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/simulations_trees/150701_Regional_TRAIN3_SIMULATED_INTERNAL.R' )	
	ch				<- subset(tinfo, SC=='150701_REGIONAL_TRAIN3' & BRL_T=='time', TAXA)
	ch[, IDPOP:= as.integer(gsub('IDPOP_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',1)))]
	ch[, GENDER_CH:= sapply(strsplit(TAXA,'|',fixed=1),'[[',2)]
	ch[, DOB_CH:= as.numeric(gsub('DOB_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',3)))]
	ch[, TIME_SEQ_CH:= as.numeric(sapply(strsplit(TAXA,'|',fixed=1),'[[',4))]
	ch				<- merge(subset(df.inds, select=c(IDPOP, GENDER, DOB, TIME_SEQ)), ch, by='IDPOP')	
	stopifnot( ch[, all(abs(DOB-DOB_CH)<=2*.Machine$double.eps)], ch[, all(GENDER==GENDER_CH)], ch[, all(abs(TIME_SEQ-TIME_SEQ_CH)<=0.001)] )
	# OK :-) schedule adding IDPOP_T
	df.trms[, IDPOP:= df.trms[, paste('IDPOP_',IDREC,sep='')]]
	tmp				<- subset(df.trms, select=c(IDPOP, IDTR))	
	df.trms[, IDPOP:= df.trms[, paste('IDPOP_',IDTR,sep='')]]
	tmp				<- merge(tmp, subset(df.trms, select=c(IDPOP, IDREC)), by='IDPOP', all=1)
	set(ch, NULL, 'IDPOP', ch[, paste('IDPOP_',IDPOP,sep='')])
	tmp		<- merge(tmp, subset(ch, select=IDPOP), by='IDPOP')
	tmp[, SC:='150701_REGIONAL_TRAIN3']
	tinfo.add		<- rbind(tinfo.add, tmp)
	# check TRAIN5
	ch				<- subset(tinfo, SC=='150701_REGIONAL_TRAIN5' & BRL_T=='time', TAXA)
	ch[, IDPOP:= as.integer(gsub('IDPOP_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',1)))]
	ch[, GENDER_CH:= sapply(strsplit(TAXA,'|',fixed=1),'[[',2)]
	ch[, DOB_CH:= as.numeric(gsub('DOB_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',3)))]
	ch[, TIME_SEQ_CH:= as.numeric(sapply(strsplit(TAXA,'|',fixed=1),'[[',4))]
	ch				<- merge(subset(df.inds, select=c(IDPOP, GENDER, DOB, TIME_SEQ)), ch, by='IDPOP')	
	stopifnot( ch[, all(abs(DOB-DOB_CH)<=2*.Machine$double.eps)], ch[, all(GENDER==GENDER_CH)], ch[, all(abs(TIME_SEQ-TIME_SEQ_CH)<=0.001)] )
	tmp[, SC:='150701_REGIONAL_TRAIN5']
	tinfo.add		<- rbind(tinfo.add, tmp)
	# check GTRFIXED2
	load( '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/simulations_trees/161121_GTRFIXED2_SIMULATED_INTERNAL.R' )	
	ch				<- subset(tinfo, SC=='161121_REGIONAL_GTRFIXED2' & BRL_T=='time', TAXA)
	ch[, IDPOP:= as.integer(gsub('IDPOP_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',1)))]
	ch[, GENDER_CH:= sapply(strsplit(TAXA,'|',fixed=1),'[[',2)]
	ch[, DOB_CH:= as.numeric(gsub('DOB_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',3)))]
	ch[, TIME_SEQ_CH:= as.numeric(sapply(strsplit(TAXA,'|',fixed=1),'[[',4))]
	ch				<- merge(subset(df.inds, select=c(IDPOP, GENDER, DOB, TIME_SEQ)), ch, by='IDPOP')	
	stopifnot( ch[, all(abs(DOB-DOB_CH)<=2*.Machine$double.eps)], ch[, all(GENDER==GENDER_CH)], ch[, all(abs(TIME_SEQ-TIME_SEQ_CH)<=0.001)] )	
	# OK :-) schedule adding IDPOP_T
	df.trms[, IDPOP:= df.trms[, paste('IDPOP_',IDREC,sep='')]]
	tmp				<- subset(df.trms, select=c(IDPOP, IDTR))	
	df.trms[, IDPOP:= df.trms[, paste('IDPOP_',IDTR,sep='')]]
	tmp				<- merge(tmp, subset(df.trms, select=c(IDPOP, IDREC)), by='IDPOP', all=1)
	set(ch, NULL, 'IDPOP', ch[, paste('IDPOP_',IDPOP,sep='')])
	tmp		<- merge(tmp, subset(ch, select=IDPOP), by='IDPOP')
	tmp[, SC:='161121_REGIONAL_GTRFIXED2']
	tinfo.add		<- rbind(tinfo.add, tmp)	
	# check GTRFIXED3
	load( '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/simulations_trees/161121_GTRFIXED3_SIMULATED_INTERNAL.R' )	
	ch				<- subset(tinfo, SC=='161121_REGIONAL_GTRFIXED3' & BRL_T=='time', TAXA)
	ch[, IDPOP:= as.integer(gsub('IDPOP_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',1)))]
	ch[, GENDER_CH:= sapply(strsplit(TAXA,'|',fixed=1),'[[',2)]
	ch[, DOB_CH:= as.numeric(gsub('DOB_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',3)))]
	ch[, TIME_SEQ_CH:= as.numeric(sapply(strsplit(TAXA,'|',fixed=1),'[[',4))]
	ch				<- merge(subset(df.inds, select=c(IDPOP, GENDER, DOB, TIME_SEQ)), ch, by='IDPOP')	
	stopifnot( ch[, all(abs(DOB-DOB_CH)<=2*.Machine$double.eps)], ch[, all(GENDER==GENDER_CH)], ch[, all(abs(TIME_SEQ-TIME_SEQ_CH)<=0.001)] )	
	# OK :-) schedule adding IDPOP_T
	df.trms[, IDPOP:= df.trms[, paste('IDPOP_',IDREC,sep='')]]
	tmp				<- subset(df.trms, select=c(IDPOP, IDTR))	
	df.trms[, IDPOP:= df.trms[, paste('IDPOP_',IDTR,sep='')]]
	tmp				<- merge(tmp, subset(df.trms, select=c(IDPOP, IDREC)), by='IDPOP', all=1)
	set(ch, NULL, 'IDPOP', ch[, paste('IDPOP_',IDPOP,sep='')])
	tmp		<- merge(tmp, subset(ch, select=IDPOP), by='IDPOP')
	tmp[, SC:='161121_REGIONAL_GTRFIXED3']
	tinfo.add		<- rbind(tinfo.add, tmp)	
	# check GTRFIXED1
	load( '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/simulations_trees/161121_GTRFIXED1_SIMULATED_INTERNAL.R' )	
	ch				<- subset(tinfo, SC=='161121_REGIONAL_GTRFIXED1' & BRL_T=='time', TAXA)
	ch[, IDPOP:= as.integer(gsub('IDPOP_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',1)))]
	ch[, GENDER_CH:= sapply(strsplit(TAXA,'|',fixed=1),'[[',2)]
	ch[, DOB_CH:= as.numeric(gsub('DOB_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',3)))]
	ch[, TIME_SEQ_CH:= as.numeric(sapply(strsplit(TAXA,'|',fixed=1),'[[',4))]
	ch				<- merge(subset(df.inds, select=c(IDPOP, GENDER, DOB, TIME_SEQ)), ch, by='IDPOP')	
	stopifnot( ch[, all(abs(DOB-DOB_CH)<=2*.Machine$double.eps)], ch[, all(GENDER==GENDER_CH)], ch[, all(abs(TIME_SEQ-TIME_SEQ_CH)<=0.001)] )	
	# OK :-) schedule adding IDPOP_T
	df.trms[, IDPOP:= df.trms[, paste('IDPOP_',IDREC,sep='')]]
	tmp				<- subset(df.trms, select=c(IDPOP, IDTR))	
	df.trms[, IDPOP:= df.trms[, paste('IDPOP_',IDTR,sep='')]]
	tmp				<- merge(tmp, subset(df.trms, select=c(IDPOP, IDREC)), by='IDPOP', all=1)
	set(ch, NULL, 'IDPOP', ch[, paste('IDPOP_',IDPOP,sep='')])
	tmp		<- merge(tmp, subset(ch, select=IDPOP), by='IDPOP')
	tmp[, SC:='161121_REGIONAL_GTRFIXED1']
	tinfo.add		<- rbind(tinfo.add, tmp)
	#
	# 	add transmitters for regional to tinfo
	#
	tinfo			<- merge(tinfo, tinfo.add, by=c('IDPOP', 'SC'), all.x=1)
	#
	#	add true node depths to tinfo
	#
	tmp				<- tinfo[,	{
				cat(IDX_T,'\n')
				ph<- ttrs[[IDX_T]]
				list(DEPTH_T=node.depth.edgelength(ph)[seq_len(Ntip(ph))], TAXA=ph$tip.label)
			}, by='IDX_T']
	tinfo			<- merge(tinfo, tmp, by=c('IDX_T','TAXA'), all.x=1)
	#
	# compute closest individual on true trees
	#
	tmp				<- unique(subset(tinfo, select=c(SC, BRL_T, IDX_T)))	
	tmp				<- tmp[, {
				print(IDX_T)
				ph			<- ttrs[[IDX_T]]
				model.reg	<- grepl('REGIONAL',SC)
				treedist.closest.ind(ph, model.reg)
			}, by=c('SC','BRL_T','IDX_T')]
	tinfo			<- merge(tinfo, tmp, by=c('SC','BRL_T','IDX_T','IDPOP'))
	set(tinfo, NULL, 'IDPOP_CL', tinfo[, gsub('IDPOP_','',IDPOP_CL)])	
	#
	#	add if transmitter sampled
	#
	tmp				<- subset(tinfo, grepl('REGIONAL',SC))	
	set(tmp, NULL, 'IDPOP', tmp[,as.integer(gsub('IDPOP_','',IDPOP))])
	setkey(tmp, IDX_T, IDPOP)
	tmp	<- unique(tmp)[, {
				z	<- IDX_T
				list(IDTR_SAMPLED=ifelse(IDTR%in%subset(tmp, IDX_T==z)[['IDPOP']], 'Y', 'N'))
			}, by=c('IDX_T','IDPOP')]
	set(tmp, NULL, 'IDPOP', tmp[, paste('IDPOP_',IDPOP,sep='')])
	tinfo	<- merge(tinfo, tmp, by=c('IDX_T','IDPOP'),all.x=1)
	#
	#	get submitted trees
	#	
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/running_gaps2'
	infiles	<- list.files(indir, pattern='newick$', recursive=1, full.names=1)
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/simple_GTR'
	infiles	<- c(infiles, list.files(indir, pattern='newick$', recursive=1, full.names=1))
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/partiallen'
	infiles	<- c(infiles, list.files(indir, pattern='newick$', recursive=1, full.names=1))
	infiles	<- data.table(FILE=infiles)
	strs	<- lapply(infiles[, FILE], function(x)
			{
				cat('\n',x)
				read.tree(file=x)	
			})
	names(strs)	<- infiles[, FILE]
	#
	#
	#
	submitted.info			<- data.table(FILE=names(strs))
	submitted.info[, IDX:=seq_along(strs)]
	#
	# set team
	#			
	submitted.info[, TEAM:=NA_character_]
	set(submitted.info, submitted.info[, which(grepl('running_gaps',FILE))], 'TEAM', 'RUNGAPS2')
	set(submitted.info, submitted.info[, which(grepl('simple_GTR',FILE))], 'TEAM', 'GTRFIXED')	
	set(submitted.info, submitted.info[, which(grepl('partiallen',FILE))], 'TEAM', 'PLEN')
	stopifnot( submitted.info[, length(which(is.na(TEAM)))==0] )	
	#
	#	scenario
	#	
	submitted.info[, SC:=NA_character_]
	tmp		<- submitted.info[, which(grepl('161121_GTRFIXED[0-9]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, gsub('161121_','161121_REGIONAL_',regmatches(FILE, regexpr('161121_GTRFIXED[0-9]',FILE)))])
	tmp		<- submitted.info[, which(grepl('150701_Regional_TRAIN2', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, regmatches(FILE, regexpr('150701_Regional_TRAIN2',FILE))])
	tmp		<- submitted.info[, which(grepl('161125_Regional_TRAIN1', FILE))]
	set(submitted.info, tmp, 'SC', '150701_REGIONAL_TRAIN1')
	set(submitted.info, NULL, 'SC', submitted.info[, toupper(SC)])
	stopifnot( submitted.info[, length(which(is.na(SC)))==0] )
	#
	#	define running gaps for the running gaps analysis
	#
	submitted.info[, RUNGGAPS_EXCL:=NA_real_]	
	set(submitted.info, tmp, 'RUNGAPS', submitted.info[tmp, as.numeric(gsub('.*TRAIN[0-9]([0-9][0-9]).*','\\1',regmatches(FILE,regexpr('TRAIN[0-9]+',FILE))))/100])	
	#
	#	define running gaps2 selected fraction
	#
	submitted.info[, RUNGAPS_EXCL:=NA_real_]	
	tmp		<- submitted.info[, which(TEAM=='RUNGAPS2')]
	set(submitted.info, tmp, 'RUNGAPS_EXCL', submitted.info[tmp, as.numeric(gsub('.*TRAIN[0-9][0-9][0-9]([0-9][0-9]).*','\\1',FILE))/100])
	#
	#	define partial length
	#
	submitted.info[, PLEN:=NA_real_]	
	tmp		<- submitted.info[, which(TEAM=='PLEN')]
	set(submitted.info, tmp, 'PLEN', submitted.info[tmp, as.numeric(gsub('PL','',regmatches(FILE, regexpr('PL[0-9]+',FILE))))])	
	#
	#	set covariates of scenarios
	#
	tmp		<- data.table(	SC=		c("150701_REGIONAL_TRAIN1","150701_REGIONAL_TRAIN2","161121_REGIONAL_GTRFIXED1","161121_REGIONAL_GTRFIXED2","161121_REGIONAL_GTRFIXED3","150701_REGIONAL_TRAIN3","150701_REGIONAL_TRAIN4" ,"150701_REGIONAL_TRAIN5", "150701_VILL_SCENARIO-A", "150701_VILL_SCENARIO-B", "VILL_99_APR15","150701_VILL_SCENARIO-C", "150701_VILL_SCENARIO-D", "150701_VILL_SCENARIO-E"),
			MODEL=	c('R','R','R','R','R','R','R','R','V','V','V','V','V','V'),
			SEQCOV= c(0.16, 0.16, 0.16, 0.16, 0.16, 0.16, 0.16, 0.16, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6),
			ACUTE=	c('low', 'low', 'low', 'low', 'low',  'high', 'low', 'high', 'high', 'high', 'high', 'high', 'high', 'high'),
			GAPS=	c('none', 'low', 'none', 'low', 'high', 'low', 'high', 'high', 'low', 'high', 'none', 'none', 'low', 'high'), 
			ART=	c('none', 'none', 'none', 'none', 'none', 'none', 'none', 'none', 'none', 'none', 'fast', 'fast', 'fast', 'fast'),
			EXT= 	c('5pc', '5pc', '5pc', '5pc', '5pc', '5pc', '5pc', '5pc', '~0pc', '~0pc', '~0pc', '~0pc', '~0pc', '~0pc')	)
	submitted.info	<- merge(submitted.info, tmp, by='SC')
	#
	#	set which gene used to construct tree (either pol or concatenated gag+pol+env)
	#
	submitted.info[, GENE:=NA_character_]
	set(submitted.info, submitted.info[, which(grepl('FULL', FILE))], 'GENE', 'GAG+POL+ENV')	
	set(submitted.info, submitted.info[, which(grepl('GAG', FILE))], 'GENE', 'GAG')
	set(submitted.info, submitted.info[, which(grepl('P17', FILE))], 'GENE', 'P17')
	set(submitted.info, submitted.info[, which(grepl('GTRFIXED', FILE))], 'GENE', 'GAG+POL+ENV')
	set(submitted.info, submitted.info[, which(grepl('TRAIN1_PL', FILE))], 'GENE', 'GAG+POL+ENV')
	stopifnot(nrow(subset(submitted.info, is.na(GENE)))==0)
	#
	#	best tree for each scenario
	#
	submitted.info[, BEST:='N']
	#
	#	set OTHER (ie old or some preliminary/unknown tree)
	#
	submitted.info[, OTHER:='N']
	#
	#	add BRL_UNITS
	#
	submitted.info[, BRL:='subst']
	#
	#	number taxa in tree
	#
	setkey(submitted.info, IDX)
	submitted.info[, TAXAN:= sapply(strs, Ntip)]	
	#
	#	add index of true tree
	#
	require(phangorn)
	tmp				<- subset(tfiles, select=c('SC','IDX_T','BRL_T'))
	setkey(tmp, SC, BRL_T)
	tmp				<- unique(tmp)
	tmp				<- dcast.data.table(tmp, SC~BRL_T, value.var='IDX_T')
	setnames(tmp, c('subst','time'), c("SUB_IDX_T","TIME_IDX_T"))
	submitted.info	<- merge(submitted.info, tmp, by='SC')
	submitted.info	<- merge(submitted.info, unique(subset(tfiles, select=c('SC','TAXAN_T'))), by='SC')	
	stopifnot(nrow(subset(submitted.info, TAXAN>TAXAN_T))==0)
	#
	#	fix taxa names that teams have changed
	#
	tmp		<- subset(submitted.info, TEAM=='IQTree' & MODEL=='R')[, IDX]
	for(i in tmp)
	{
		strs[[i]]$tip.label		<- sapply(strsplit(strs[[i]]$tip.label,'_'), function(x)	paste(x[1],'_',x[2],'|',x[3],'|',x[4],'_',x[5],'|',x[6],sep='')	)		
	}
	for(i in seq_along(strs))
	{
		strs[[i]]$tip.label		<- toupper(strs[[i]]$tip.label)		
	}
	for(i in seq_along(ttrs))
	{
		ttrs[[i]]$tip.label	<- toupper(ttrs[[i]]$tip.label)
	}		
	tmp2	<- subset(tinfo, BRL_T=='time', select=c(IDPOP,SC,TAXA))
	setkey(tmp2, IDPOP,SC,TAXA)
	tmp2	<- unique(tmp2)
	tmp		<- subset(submitted.info, TEAM=='PhyML' & MODEL=='R')[, IDX]
	for(i in tmp)
	{		
		cat(i,'\n')
		z	<- data.table(IDX=seq_along(strs[[i]]$tip.label), IDPOP=regmatches(strs[[i]]$tip.label, regexpr('IDPOP_[0-9]+',strs[[i]]$tip.label)), SC=subset(submitted.info, IDX==i)[,SC])
		z	<- merge(tmp2, z, by=c('IDPOP','SC'))
		setkey(z, IDX)
		stopifnot(nrow(z)==Ntip(strs[[i]]))
		strs[[i]]$tip.label	<- z[, TAXA]			
	}
	tmp		<- subset(submitted.info, TEAM=='PhyML' & MODEL=='V')[, IDX]
	for(i in tmp)
	{
		cat(i,'\n')
		z	<- data.table(IDX=seq_along(strs[[i]]$tip.label), IDPOP=regmatches(strs[[i]]$tip.label, regexpr('HOUSE[0-9]+-[0-9]+|House[0-9]+-[0-9]+',strs[[i]]$tip.label)), SC=subset(submitted.info, IDX==i)[,SC])
		z	<- merge(tmp2, z, by=c('IDPOP','SC'))
		stopifnot(nrow(z)==length(strs[[i]]$tip.label))
		setkey(z, IDX)
		strs[[i]]$tip.label	<- z[, TAXA]		
	}
	#
	#	check labels and remove labels that do not appear in the observed tree
	#
	tmp	<- submitted.info[, {				
				stree		<- unroot(strs[[IDX]])
				otree		<- unroot(ttrs[[TIME_IDX_T]])				
				z			<- setdiff(otree$tip.label, stree$tip.label)
				list(CHECK= length(z)==abs(diff(c(TAXAN, TAXAN_T))) )
			}, by='IDX']
	tmp	<- merge(subset(tmp, !CHECK), submitted.info, by='IDX')
	for(i in seq_len(nrow(tmp)))
	{
		j			<- tmp[i, IDX]
		cat('\n',j)		
		otree		<- tmp[i, TIME_IDX_T]
		stree		<- unroot(strs[[j]])
		otree		<- unroot(ttrs[[otree]])					
		z			<- merge( data.table(TAXA=stree$tip.label, TYPE='s'), data.table(TAXA=otree$tip.label, TYPE='o'), by='TAXA', all=1)			
		strs[[j]]	<- drop.tip(stree, subset( z, is.na(TYPE.y))[, TAXA])		
	}
	submitted.info[, TAXAN:= sapply(strs, Ntip)]
	#
	#	are trees rooted?
	#
	setkey(submitted.info, IDX)
	submitted.info[, ROOTED:=factor(sapply(strs, is.rooted),levels=c(TRUE,FALSE),labels=c('Y','N'))]
	#
	#	read LSD trees
	#
	indir			<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/LSD_2'
	infiles			<- data.table(FILE=list.files(indir, pattern='LSD.date.newick', full.names=TRUE))
	infiles[, IDX:= as.integer(gsub('IDX_','',regmatches(basename(FILE),regexpr('IDX_[0-9]+',basename(FILE)))))]
	setkey(infiles, IDX)
	setdiff(1:385, infiles[,IDX])
	
	strs_lsd		<- vector('list', submitted.info[, max(IDX)])
	for(i in seq_len(nrow(infiles)))		
	{
		#	i<- 439
		IDX						<- infiles[i,IDX]
		FILE					<- infiles[i,FILE]
		cat('\n',IDX)
		ph						<- read.tree(FILE)
		stopifnot( !is.null(ph) )
		stopifnot( identical(sort(strs_rtt[[IDX]]$tip.label), sort(ph$tip.label)) )
		strs_lsd[[IDX]]			<- ph
		#names(strs_lsd[[IDX]])	<- FILE					
	}
	setkey(submitted.info, IDX)
	submitted.info[, WITH_LSD:= factor(sapply(strs_lsd, is.null), levels=c(TRUE,FALSE), labels=c('N','Y'))]	 	
	#
	#	SAVE so far
	#
	outdir		<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation'	
	save(strs, strs_lsd, ttrs, tinfo, tfiles, tbrl, submitted.info, file=file.path(outdir,'submitted_161123.rda'))	
}
##--------------------------------------------------------------------------------------------------------
##	olli 27.06.16
##--------------------------------------------------------------------------------------------------------
treecomparison.submissions.160627<- function()	
{
	require(data.table)
	require(ape)
	require(adephylo)
	require(phangorn)
	require(parallel)
	#
	#	get true trees
	#
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim_internal/freeze_July15'
	tfiles	<- list.files(indir, pattern='newick$', full.names=TRUE)
	tfiles	<- data.table( FILE_T= tfiles[ grepl('SUBSTTREE', tfiles) | grepl('Vill_99', tfiles) | grepl('Vill.*DATEDTREE', tfiles) ] )
	tfiles[, SC:= toupper(gsub('_SUBSTTREE|_DATEDTREE','',gsub('.newick','',basename(FILE_T))))]
	tmp		<- rbind( subset(tfiles, SC=='VILL_99_APR15'), subset(tfiles, SC=='VILL_99_APR15'), subset(tfiles, SC=='VILL_99_APR15') )
	set(tmp, NULL, 'SC', c('150701_VILL_SCENARIO-C','150701_VILL_SCENARIO-D','150701_VILL_SCENARIO-E'))
	tfiles	<- rbind(tfiles, tmp)
	tmp		<- list.files(indir, pattern='newick$', full.names=TRUE)
	tmp		<- data.table( FILE_T= tmp[ grepl('Reg.*DATEDTREE', tmp) ] )
	tmp[, SC:= toupper(gsub('_SUBSTTREE|_DATEDTREE','',gsub('.newick','',basename(FILE_T))))]
	tfiles	<- rbind(tfiles, tmp)
	tfiles[, BRL_T:= 'time']	
	set(tfiles, tfiles[, which(grepl('REG',SC) & grepl('SUBST',FILE_T))], 'BRL_T', 'subst')	
	ttrs	<- lapply(tfiles[, FILE_T], function(x)	read.tree(file=x) )
	names(ttrs)	<- tfiles[, SC]	
	for(z in c('VILL_99_APR15','150701_VILL_SCENARIO-C','150701_VILL_SCENARIO-D','150701_VILL_SCENARIO-E'))
		ttrs[[z]]	<- root(ttrs[[z]], node=Ntip(ttrs[[z]])+2, resolve.root=1)	
	tfiles[, IDX_T:=seq_along(ttrs)]
	tfiles[, TAXAN_T:= sapply(ttrs, Ntip)]
	#	patristic distances on true trees (by time and subst/site)
	tbrl	<- lapply(seq_len(nrow(tfiles)), function(i)
					{
						ph		<- ttrs[[tfiles[i, IDX_T]]]
						tmp		<- cophenetic.phylo(ph)				
						tmp		<- as.data.table(melt(tmp))
						setnames(tmp, c('Var1','Var2','value'),c('TAXA1','TAXA2','PD'))
						tmp		<- subset(tmp, TAXA1!=TAXA2)
						tmp[, IDX_T:= tfiles[i, IDX_T]]
						tmp[, SC:= tfiles[i, SC]]
						tmp[, BRL_T:= tfiles[i, BRL_T]]
						tmp[, TAXAN_T:= tfiles[i, TAXAN_T]]
						tmp
					})
	tbrl	<- do.call('rbind',tbrl)		
	#	info on true trees
	tinfo	<- merge(tfiles, do.call('rbind',lapply(seq_along(ttrs), function(i) data.table(TAXA=ttrs[[i]]$tip.label, IDX_T=i))), by='IDX_T')	
	tinfo[, IDPOP:=NA_character_]
	tmp		<- tinfo[, which(grepl('REGIONAL',SC))]
	set(tinfo, tmp, 'IDPOP', tinfo[tmp,regmatches(TAXA, regexpr('IDPOP_[0-9]+',TAXA))])
	tmp		<- tinfo[, which(!grepl('REGIONAL',SC))]
	set(tinfo, tmp, 'IDPOP', tinfo[tmp, regmatches(TAXA, regexpr('HOUSE[0-9]+-[0-9]+|House[0-9]+-[0-9]+',TAXA))])		
	stopifnot(subset(tinfo, grepl('VILL',SC))[, length(which(substring(TAXA,1,10)!=substring(IDPOP,1,10)))]==0)	
	stopifnot( tinfo[, length(which(is.na(IDPOP)))==0] )	
	set(tinfo, NULL, 'IDPOP', tinfo[,toupper(IDPOP)])
	set(tinfo, NULL, 'TAXA', tinfo[,toupper(TAXA)])
	#	read cluster membership from DATEDCLUTREES	
	tmp		<- list.files(indir, pattern='DATEDCLUTREES', full.names=TRUE)
	tmp		<- rbind( data.table( 	FILE_CLU_T= tmp, 
									SC= toupper(gsub('_DATEDCLUTREES','',gsub('.newick','',basename(tmp)))),
									BRL_T= 'time'),
					  data.table( 	FILE_CLU_T= tmp, 
									SC= toupper(gsub('_DATEDCLUTREES','',gsub('.newick','',basename(tmp)))),
									BRL_T= 'subst') )
	tfiles	<- merge(tfiles, tmp, by=c('SC','BRL_T'), all=1)	
	tmp		<- subset(tfiles, !is.na(FILE_CLU_T))[, {
				z		<- read.tree(FILE_CLU_T)
				do.call('rbind',lapply(seq_along(z), function(i) data.table(IDCLU=i, TAXA=z[[i]]$tip.label)))				
			}, by=c('SC','BRL_T')]	
	tinfo	<- merge(tinfo, tmp, by=c('SC','BRL_T','TAXA'), all=1)
	tmp		<- subset(tinfo, !is.na(IDCLU))[, list(CLU_N= length(IDPOP)), by=c('SC','BRL_T','IDCLU')]
	tinfo	<- merge(tinfo, tmp, by=c('SC','BRL_T','IDCLU'), all=1)
	#	read sequences and determine %gappiness
	tmp		<- list.files(indir, pattern='fa$|fasta$', full.names=TRUE)
	tmp		<- data.table( FILE_SEQ_T= tmp, SC= toupper(gsub('_SIMULATED','',gsub('.fa','',basename(tmp)))))
	z		<- subset(tmp, SC=='VILL_99_APR15')
	set(z, NULL, 'SC', '150701_VILL_SCENARIO-C')
	tmp		<- rbind( tmp, z )	
	tfiles	<- merge(tfiles, tmp, by='SC', all=1)
	tmp		<- subset(tfiles, !is.na(FILE_SEQ_T))[, {
				z		<- read.dna(FILE_SEQ_T, format='fasta')	
				ans		<- sapply(seq_len(nrow(z)), function(i) base.freq(z[i,], all=1))
				ans		<- apply(ans[c('n','-','?'),], 2, sum)
				list(TAXA=rownames(z), GPS=ans)				
			}, by=c('SC','BRL_T')]
	tinfo	<- merge(tinfo, tmp, by=c('SC','BRL_T','TAXA'), all.x=1)
	#	add % gaps by gene for regional
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/simulations'
	infiles	<- data.table(FILE=list.files(indir, full.names=1, pattern='_TRAIN[0-9]+_SIMULATED.fa$'))
	infiles[, SC:= toupper(gsub('_SIMULATED.fa','',basename(FILE)))]
	tmp		<- infiles[, {
				#FILE<- "/Users/Oliver/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/simulations/150701_Regional_TRAIN4_SIMULATED.fa"
				sq	<- read.dna(FILE, format='fa')
				seqi<- as.data.table(read.csv(gsub('.fa','_gene.txt',FILE), header=FALSE))
				seqi[, GENE:= regmatches(V2, regexpr('[a-z]+', V2))]
				seqi[, START:= as.integer(gsub('-','',regmatches(V2, regexpr('[0-9]+-', V2))))]
				seqi[, END:= as.integer(gsub('-','',regmatches(V2, regexpr('-[0-9]+', V2))))]
				seqi	<- subset(seqi, select=c(GENE,START,END))
				seqi	<- rbind(seqi, data.table(GENE='full', START=1L, END= seqi[, max(END)]))				
				ans		<- seqi[, {
							#START<- 1474; END<- 5466
							z	<- as.character(sq[,START:END])
							tmp	<- apply(z, 2, function(x) all(x%in%c('-','?'))) 							
							z	<- z[, !tmp]
							tmp	<- apply(z, 1, function(x) length(which(x=='?'))) / ncol(z)
							list(TAXA=names(tmp), GAPS_P=tmp)
						}, by='GENE']
				set(ans, NULL, 'GENE', ans[, paste(toupper(GENE),'_GAPS_P',sep='')])
				ans		<- dcast.data.table(ans, TAXA~GENE, value.var='GAPS_P')
				ans
			}, by='SC']
	tmp		<- tmp[, list(FULL_GAPS_P=mean(FULL_GAPS_P), GAG_GAPS_P=mean(GAG_GAPS_P), POL_GAPS_P=mean(POL_GAPS_P), ENV_GAPS_P=mean(ENV_GAPS_P)), by='SC']
	tinfo	<- merge(tinfo, tmp, all.x=1, by='SC')	
	#
	# to tinfo add actual transmitters
	#
	# check TRAIN1
	load( '/Users/Oliver/Dropbox (Infectious Disease)/PANGEAHIVsim_internal/freeze_July15/150701_Regional_TRAIN1_SIMULATED_INTERNAL.R' )	
	ch				<- subset(tinfo, SC=='150701_REGIONAL_TRAIN1' & BRL_T=='time', TAXA)
	ch[, IDPOP:= as.integer(gsub('IDPOP_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',1)))]
	ch[, GENDER_CH:= sapply(strsplit(TAXA,'|',fixed=1),'[[',2)]
	ch[, DOB_CH:= as.numeric(gsub('DOB_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',3)))]
	ch[, TIME_SEQ_CH:= as.numeric(sapply(strsplit(TAXA,'|',fixed=1),'[[',4))]
	ch				<- merge(subset(df.inds, select=c(IDPOP, GENDER, DOB, TIME_SEQ)), ch, by='IDPOP')	
	stopifnot( ch[, all(DOB==DOB_CH)], ch[, all(GENDER==GENDER_CH)], ch[, all(TIME_SEQ==TIME_SEQ_CH)] )
	subset(ch, TIME_SEQ!=TIME_SEQ_CH)
	# OK :-) schedule adding IDPOP_T
	df.trms[, IDPOP:= df.trms[, paste('IDPOP_',IDREC,sep='')]]
	tinfo.add		<- subset(df.trms, select=c(IDPOP, IDTR))	
	df.trms[, IDPOP:= df.trms[, paste('IDPOP_',IDTR,sep='')]]
	tinfo.add		<- merge(tinfo.add, subset(df.trms, select=c(IDPOP, IDREC)), by='IDPOP', all=1)
	set(ch, NULL, 'IDPOP', ch[, paste('IDPOP_',IDPOP,sep='')])
	tinfo.add		<- merge(tinfo.add, subset(ch, select=IDPOP), by='IDPOP')
	tinfo.add[, SC:='150701_REGIONAL_TRAIN1']
	# check TRAIN2
	load( '/Users/Oliver/Dropbox (Infectious Disease)/PANGEAHIVsim_internal/freeze_July15/150701_Regional_TRAIN2_SIMULATED_INTERNAL.R' )	
	ch				<- subset(tinfo, SC=='150701_REGIONAL_TRAIN2' & BRL_T=='time', TAXA)
	ch[, IDPOP:= as.integer(gsub('IDPOP_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',1)))]
	ch[, GENDER_CH:= sapply(strsplit(TAXA,'|',fixed=1),'[[',2)]
	ch[, DOB_CH:= as.numeric(gsub('DOB_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',3)))]
	ch[, TIME_SEQ_CH:= as.numeric(sapply(strsplit(TAXA,'|',fixed=1),'[[',4))]
	ch				<- merge(subset(df.inds, select=c(IDPOP, GENDER, DOB, TIME_SEQ)), ch, by='IDPOP')	
	stopifnot( ch[, all(DOB==DOB_CH)], ch[, all(GENDER==GENDER_CH)], ch[, all(TIME_SEQ==TIME_SEQ_CH)] )
	subset(ch, TIME_SEQ!=TIME_SEQ_CH)
	# OK :-) schedule adding IDPOP_T
	df.trms[, IDPOP:= df.trms[, paste('IDPOP_',IDREC,sep='')]]
	tmp				<- subset(df.trms, select=c(IDPOP, IDTR))	
	df.trms[, IDPOP:= df.trms[, paste('IDPOP_',IDTR,sep='')]]
	tmp				<- merge(tmp, subset(df.trms, select=c(IDPOP, IDREC)), by='IDPOP', all=1)
	set(ch, NULL, 'IDPOP', ch[, paste('IDPOP_',IDPOP,sep='')])
	tmp		<- merge(tmp, subset(ch, select=IDPOP), by='IDPOP')
	tmp[, SC:='150701_REGIONAL_TRAIN2']
	tinfo.add		<- rbind(tinfo.add, tmp)
	# check TRAIN4
	ch				<- subset(tinfo, SC=='150701_REGIONAL_TRAIN4' & BRL_T=='time', TAXA)
	ch[, IDPOP:= as.integer(gsub('IDPOP_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',1)))]
	ch[, GENDER_CH:= sapply(strsplit(TAXA,'|',fixed=1),'[[',2)]
	ch[, DOB_CH:= as.numeric(gsub('DOB_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',3)))]
	ch[, TIME_SEQ_CH:= as.numeric(sapply(strsplit(TAXA,'|',fixed=1),'[[',4))]
	ch				<- merge(subset(df.inds, select=c(IDPOP, GENDER, DOB, TIME_SEQ)), ch, by='IDPOP')	
	stopifnot( ch[, all(DOB==DOB_CH)], ch[, all(GENDER==GENDER_CH)], ch[, all(TIME_SEQ==TIME_SEQ_CH)] )
	subset(ch, TIME_SEQ!=TIME_SEQ_CH)	
	tmp[, SC:='150701_REGIONAL_TRAIN4']
	tinfo.add		<- rbind(tinfo.add, tmp)
	# check TRAIN3
	load( '/Users/Oliver/Dropbox (Infectious Disease)/PANGEAHIVsim_internal/freeze_July15/150701_Regional_TRAIN3_SIMULATED_INTERNAL.R' )	
	ch				<- subset(tinfo, SC=='150701_REGIONAL_TRAIN3' & BRL_T=='time', TAXA)
	ch[, IDPOP:= as.integer(gsub('IDPOP_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',1)))]
	ch[, GENDER_CH:= sapply(strsplit(TAXA,'|',fixed=1),'[[',2)]
	ch[, DOB_CH:= as.numeric(gsub('DOB_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',3)))]
	ch[, TIME_SEQ_CH:= as.numeric(sapply(strsplit(TAXA,'|',fixed=1),'[[',4))]
	ch				<- merge(subset(df.inds, select=c(IDPOP, GENDER, DOB, TIME_SEQ)), ch, by='IDPOP')	
	stopifnot( ch[, all(DOB==DOB_CH)], ch[, all(GENDER==GENDER_CH)], ch[, all(TIME_SEQ==TIME_SEQ_CH)] )
	subset(ch, TIME_SEQ!=TIME_SEQ_CH)
	# OK :-) schedule adding IDPOP_T
	df.trms[, IDPOP:= df.trms[, paste('IDPOP_',IDREC,sep='')]]
	tmp				<- subset(df.trms, select=c(IDPOP, IDTR))	
	df.trms[, IDPOP:= df.trms[, paste('IDPOP_',IDTR,sep='')]]
	tmp				<- merge(tmp, subset(df.trms, select=c(IDPOP, IDREC)), by='IDPOP', all=1)
	set(ch, NULL, 'IDPOP', ch[, paste('IDPOP_',IDPOP,sep='')])
	tmp		<- merge(tmp, subset(ch, select=IDPOP), by='IDPOP')
	tmp[, SC:='150701_REGIONAL_TRAIN3']
	tinfo.add		<- rbind(tinfo.add, tmp)
	# check TRAIN5
	ch				<- subset(tinfo, SC=='150701_REGIONAL_TRAIN5' & BRL_T=='time', TAXA)
	ch[, IDPOP:= as.integer(gsub('IDPOP_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',1)))]
	ch[, GENDER_CH:= sapply(strsplit(TAXA,'|',fixed=1),'[[',2)]
	ch[, DOB_CH:= as.numeric(gsub('DOB_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',3)))]
	ch[, TIME_SEQ_CH:= as.numeric(sapply(strsplit(TAXA,'|',fixed=1),'[[',4))]
	ch				<- merge(subset(df.inds, select=c(IDPOP, GENDER, DOB, TIME_SEQ)), ch, by='IDPOP')	
	stopifnot( ch[, all(DOB==DOB_CH)], ch[, all(GENDER==GENDER_CH)], ch[, all(TIME_SEQ==TIME_SEQ_CH)] )
	subset(ch, TIME_SEQ!=TIME_SEQ_CH)
	tmp[, SC:='150701_REGIONAL_TRAIN5']
	tinfo.add		<- rbind(tinfo.add, tmp)
	#
	# 	add transmitters for regional to tinfo
	#
	tinfo			<- merge(tinfo, tinfo.add, by=c('IDPOP', 'SC'), all.x=1)
	#
	#	add true node depths to tinfo
	#
	tmp				<- tinfo[,	{
									cat(IDX_T,'\n')
									ph<- ttrs[[IDX_T]]
									list(DEPTH_T=node.depth.edgelength(ph)[seq_len(Ntip(ph))], TAXA=ph$tip.label)
								}, by='IDX_T']
	tinfo			<- merge(tinfo, tmp, by=c('IDX_T','TAXA'), all.x=1)
	#
	# compute closest individual on true trees
	#
	tmp				<- unique(subset(tinfo, select=c(SC, BRL_T, IDX_T)))	
	tmp				<- tmp[, {
				print(IDX_T)
				ph			<- ttrs[[IDX_T]]
				model.reg	<- grepl('REGIONAL',SC)
				treedist.closest.ind(ph, model.reg)
			}, by=c('SC','BRL_T','IDX_T')]
	tinfo			<- merge(tinfo, tmp, by=c('SC','BRL_T','IDX_T','IDPOP'))
	set(tinfo, NULL, 'IDPOP_CL', tinfo[, gsub('IDPOP_','',IDPOP_CL)])	
	#
	#	add if transmitter sampled
	#
	tmp				<- subset(tinfo, grepl('REGIONAL',SC))	
	set(tmp, NULL, 'IDPOP', tmp[,as.integer(gsub('IDPOP_','',IDPOP))])
	setkey(tmp, IDX_T, IDPOP)
	tmp	<- unique(tmp)[, {
				z	<- IDX_T
				list(IDTR_SAMPLED=ifelse(IDTR%in%subset(tmp, IDX_T==z)[['IDPOP']], 'Y', 'N'))
			}, by=c('IDX_T','IDPOP')]
	set(tmp, NULL, 'IDPOP', tmp[, paste('IDPOP_',IDPOP,sep='')])
	tinfo	<- merge(tinfo, tmp, by=c('IDX_T','IDPOP'),all.x=1)
	#
	#	get submitted trees
	#	
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/IQTree/IQTree201507'
	infiles	<- list.files(indir, pattern='treefile$', recursive=1, full.names=1)
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/IQTree/IQTree201510'
	infiles	<- c(infiles, list.files(indir, pattern='treefile$', recursive=1, full.names=1))	
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/IQTree/IQTREE_Update_gag'
	infiles	<- c(infiles, list.files(indir, pattern='treefile$', recursive=1, full.names=1))	
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/PhyML'
	infiles	<- c(infiles, list.files(indir, pattern='*tree*', recursive=1, full.names=1))
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/RAxML'
	infiles	<- c(infiles, list.files(indir, pattern='*RAxML_bestTree*', recursive=1, full.names=1))
	infiles	<- c(infiles, list.files(indir, pattern="best_tree", recursive=1, full.names=1))	
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/running_gaps'
	infiles	<- c(infiles, list.files(indir, pattern="newick", recursive=1, full.names=1))	
	infiles	<- data.table(FILE=infiles)
	strs	<- lapply(infiles[, FILE], function(x)
			{
				cat(x)
				read.tree(file=x)	
			})
	names(strs)	<- infiles[, FILE]
	#
	#	add MetaPIGA full genome trees, version 160713. stored as list of newicks
	#
	indir					<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/MetaPIGA_160713'
	tmp						<- list.files(indir, pattern='*txt*', recursive=1, full.names=1)
	tmp.trees				<- lapply(tmp, function(x)
			{
				cat(x)
				read.tree(file=x)	
			})	
	MetaPIGA.trees			<- c(tmp.trees[[1]],tmp.trees[[2]],tmp.trees[[3]],tmp.trees[[4]],tmp.trees[[5]],tmp.trees[[6]])
	tmp						<- sapply( seq_along(tmp), function(i) paste(gsub('.txt','',tmp[i]), '_tree', seq_along(tmp.trees[[i]]), sep='') )
	names(MetaPIGA.trees)	<- unlist(tmp)
	strs					<- c(strs, unclass(MetaPIGA.trees))
	#
	#	to keep old index, add MVR trees now
	#
	options(warn=2)
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/tree_mvr'
	tmp		<- data.table(FILE=list.files(indir, pattern="newick", recursive=1, full.names=1))	
	mvrtrs	<- lapply(tmp[, FILE], function(x)
			{
				cat(x)
				read.tree(file=x)	
			})
	names(mvrtrs)	<- tmp[, FILE]
	strs			<- c(strs, mvrtrs)
	#
	#	add MetaPIGA trees, version 150831. stored as nexus
	#
	#indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/MetaPIGA_150831'
	#tmp		<-  list.files(indir, pattern='*result*', recursive=1, full.names=1)
	#tmp		<- data.table(FILE=tmp)	
	#tmp.trees			<- lapply(tmp[, FILE], function(x)
	#		{
	#			cat(x)
	#			read.nexus(file=x)	
	#		})
	#sapply(tmp.trees, length)
	#MetaPIGA.trees			<- c(lapply(tmp.trees, '[[', 1), lapply(tmp.trees, '[[', 2), lapply(tmp.trees, '[[', 3), lapply(tmp.trees, '[[', 4))
	#names(MetaPIGA.trees)	<- c(sapply(tmp.trees, function(x) paste(names(x)[1],'_use',sep='')), sapply(tmp.trees, function(x) names(x)[2]), sapply(tmp.trees, function(x) names(x)[3]), sapply(tmp.trees, function(x) names(x)[4]))	
	#names(MetaPIGA.trees)	<- gsub("'",'',names(MetaPIGA.trees), fixed=1)	
	#strs					<- c(strs, MetaPIGA.trees)
	#
	#
	#
	submitted.info			<- data.table(FILE=names(strs))
	submitted.info[, IDX:=seq_along(strs)]
	#
	# set team
	#			
	submitted.info[, TEAM:=NA_character_]
	set(submitted.info, submitted.info[, which(grepl('RAXML|RAxML',FILE))], 'TEAM', 'RAXML')
	set(submitted.info, submitted.info[, which(grepl('IQTree',FILE))], 'TEAM', 'IQTree')
	set(submitted.info, submitted.info[, which(grepl('MetaPIGA|Consensus pruning|Best individual of population',FILE))], 'TEAM', 'MetaPIGA')
	set(submitted.info, submitted.info[, which(grepl('PhyML',FILE))], 'TEAM', 'PhyML')
	set(submitted.info, submitted.info[, which(grepl('running_gaps',FILE))], 'TEAM', 'RUNGAPS_ExaML')
	set(submitted.info, submitted.info[, which(grepl('tree_mvr.*MVR_C_0\\.newick',FILE))], 'TEAM', 'MVR')
	set(submitted.info, submitted.info[, which(grepl('tree_mvr.*BioNJ_C_0\\.newick',FILE))], 'TEAM', 'BioNJ')
	stopifnot( submitted.info[, length(which(is.na(TEAM)))==0] )	
	#
	#	scenario
	#	
	submitted.info[, SC:=NA_character_]
	tmp		<- submitted.info[, which(grepl('150701_Regional_TRAIN[0-9]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, regmatches(FILE, regexpr('150701_Regional_TRAIN[0-9]',FILE))])
	tmp		<- submitted.info[, which(grepl('150701_Vill_SCENARIO-[A-Z]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, regmatches(FILE, regexpr('150701_Vill_SCENARIO-[A-Z]',FILE))])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('TRAIN[0-9]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, paste('150701_Regional_',regmatches(FILE, regexpr('TRAIN[0-9]',FILE)),sep='')])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('scenario[A-Z]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, paste('150701_Vill_',regmatches(FILE, regexpr('scenario[A-Z]',FILE)),sep='')])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('150701_regional_train[0-9]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, regmatches(FILE, regexpr('150701_regional_train[0-9]',FILE))])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('150701_vill_scenario-[A-Z]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, regmatches(FILE, regexpr('150701_vill_scenario-[A-Z]',FILE))])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('Vill_99_Apr15', FILE))]
	set(submitted.info, tmp, 'SC', 'Vill_99_Apr15')	
	set(submitted.info, NULL, 'SC', submitted.info[, toupper(SC)])
	tmp		<- submitted.info[, which(grepl('150701_VILL_SCENARIO[A-Z]', SC))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, gsub('150701_VILL_SCENARIO','150701_VILL_SCENARIO-',SC)])
	stopifnot( submitted.info[, length(which(is.na(SC)))==0] )
	#
	#	define running gaps for the running gaps analysis
	#
	submitted.info[, RUNGGAPS:=NA_real_]	
	tmp		<- submitted.info[, which(TEAM=='RUNGAPS_ExaML')]
	set(submitted.info, tmp, 'RUNGAPS', submitted.info[tmp, as.numeric(gsub('TRAIN[0-9]','',regmatches(FILE,regexpr('TRAIN[0-9]+',FILE))))/100])
	#
	#	set covariates of scenarios
	#
	tmp		<- data.table(	SC=		c("150701_REGIONAL_TRAIN1","150701_REGIONAL_TRAIN2","150701_REGIONAL_TRAIN3","150701_REGIONAL_TRAIN4" ,"150701_REGIONAL_TRAIN5", "150701_VILL_SCENARIO-A", "150701_VILL_SCENARIO-B", "VILL_99_APR15","150701_VILL_SCENARIO-C", "150701_VILL_SCENARIO-D", "150701_VILL_SCENARIO-E"),
							MODEL=	c('R','R','R','R','R','V','V','V','V','V','V'),
							SEQCOV= c(0.16, 0.16, 0.16, 0.16, 0.16, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6),
							ACUTE=	c('low', 'low', 'high', 'low', 'high', 'high', 'high', 'high', 'high', 'high', 'high'),
							GAPS=	c('none', 'low', 'low', 'high', 'high', 'low', 'high', 'none', 'none', 'low', 'high'), 
							ART=	c('none', 'none', 'none', 'none', 'none', 'none', 'none', 'fast', 'fast', 'fast', 'fast'),
							EXT= 	c('5pc', '5pc', '5pc', '5pc', '5pc', '~0pc', '~0pc', '~0pc', '~0pc', '~0pc', '~0pc')	)
	submitted.info	<- merge(submitted.info, tmp, by='SC')
	#
	#	set which gene used to construct tree (either pol or concatenated gag+pol+env)
	#
	submitted.info[, GENE:=NA_character_]
	set(submitted.info, submitted.info[, which(TEAM=='RAXML' & grepl('full', FILE))], 'GENE', 'GAG+POL+ENV')
	set(submitted.info, submitted.info[, which(TEAM=='RAXML' & grepl('pol', FILE))], 'GENE', 'POL')
	set(submitted.info, submitted.info[, which(TEAM=='RAXML' & grepl('gag', FILE))], 'GENE', 'GAG')
	stopifnot(nrow(subset(submitted.info, TEAM=='RAXML' & is.na(GENE)))==0)
	set(submitted.info, submitted.info[, which(TEAM=='PhyML' & grepl('gag', FILE))], 'GENE', 'GAG')
	set(submitted.info, submitted.info[, which(TEAM=='PhyML' & grepl('pol', FILE))], 'GENE', 'POL')
	set(submitted.info, submitted.info[, which(TEAM=='PhyML' & grepl('gagpolenv', FILE))], 'GENE', 'GAG+POL+ENV')
	stopifnot(nrow(subset(submitted.info, TEAM=='RAXML' & is.na(GENE)))==0)	
	set(submitted.info, submitted.info[, which(TEAM=='MetaPIGA' & grepl('gag', FILE))], 'GENE', 'GAG')
	set(submitted.info, submitted.info[, which(TEAM=='MetaPIGA' & grepl('all', FILE))], 'GENE', 'GAG+POL+ENV')
	stopifnot(nrow(subset(submitted.info, TEAM=='MetaPIGA' & is.na(GENE)))==0)
	set(submitted.info, submitted.info[, which(TEAM=='IQTree' & grepl('[0-9]_partition', FILE))], 'GENE', 'GAG+POL+ENV')
	set(submitted.info, submitted.info[, which(TEAM=='IQTree' & grepl('[0-9]_pol_partition', FILE))], 'GENE', 'POL')
	set(submitted.info, submitted.info[, which(TEAM=='IQTree' & grepl('[0-9]_gag_partition', FILE))], 'GENE', 'GAG')
	stopifnot(nrow(subset(submitted.info, TEAM=='IQTree' & is.na(GENE)))==0)	
	set(submitted.info, submitted.info[, which(TEAM=='RUNGAPS_ExaML' & grepl('FULL_SIMULATED', FILE))], 'GENE', 'GAG+POL+ENV')
	set(submitted.info, submitted.info[, which(TEAM=='RUNGAPS_ExaML' & grepl('GAG_SIMULATED', FILE))], 'GENE', 'GAG')
	set(submitted.info, submitted.info[, which(TEAM=='RUNGAPS_ExaML' & grepl('GAGPP_SIMULATED', FILE))], 'GENE', 'GAG+PARTIALPOL')
	set(submitted.info, submitted.info[, which(TEAM=='RUNGAPS_ExaML' & grepl('P17_SIMULATED', FILE))], 'GENE', 'P17')
	stopifnot(nrow(subset(submitted.info, TEAM=='RUNGAPS_ExaML' & is.na(GENE)))==0)
	set(submitted.info, submitted.info[, which(TEAM%in%c('MVR','BioNJ') & grepl('gag+pol+env', FILE))], 'GENE', 'GAG+POL+ENV')
	set(submitted.info, submitted.info[, which(TEAM%in%c('MVR','BioNJ') & grepl('pol', FILE))], 'GENE', 'POL')
	set(submitted.info, submitted.info[, which(TEAM%in%c('MVR','BioNJ') & grepl('gag', FILE))], 'GENE', 'GAG')
	set(submitted.info, submitted.info[, which(TEAM%in%c('MVR','BioNJ') & grepl('env', FILE))], 'GENE', 'ENV')
	stopifnot(nrow(subset(submitted.info, TEAM%in%c('MVR','BioNJ') & is.na(GENE)))==0)	
	#
	#	best tree for each scenario
	#
	submitted.info[, BEST:='N']
	set(submitted.info, submitted.info[, which(grepl('RAxML', FILE) & grepl('best_tree', FILE))], 'BEST', 'Y')	
	#	copied from ListOfBestTrees_IQTree150818.txt
	#	there are several best trees for some scenarios
	tmp	<- c( 	'150701_Vill_SCENARIO-A_IQTree150814_partition_12_3_07',	
				'150701_Vill_SCENARIO-A_IQTree150814_partition_12_3_04.',	
				'150701_Vill_SCENARIO-B_IQTree150814_partition_12_3_03.',	
				'Vill_99_Apr15_IQTree150814_partition_123.',			
				'150701_Vill_SCENARIO-D_IQTree150814_partition_12_3.',	
				'150701_Vill_SCENARIO-E_IQTree150814_partition_12_3.',
				'150701_Vill_SCENARIO-A_IQTree150814_pol_partition_12_3.',	
				'150701_Vill_SCENARIO-B_IQTree150814_pol_partition_12_3_05.',
				'Vill_99_Apr15_IQTree150814_pol_partition_12_3_09.',		
				'Vill_99_Apr15_IQTree150814_pol_partition_12_3_10.',		
				'150701_Vill_SCENARIO-D_IQTree150814_pol_partition_12_3_05.',
				'150701_Vill_SCENARIO-D_IQTree150814_pol_partition_12_3_06.',
				'150701_Vill_SCENARIO-D_IQTree150814_pol_partition_12_3_09.',
				'150701_Vill_SCENARIO-E_IQTree150814_pol_partition_12_3_06.',
				'150701_Regional_TRAIN1_IQTree150818_partition_123_03.',	
				'150701_Regional_TRAIN1_IQTree150818_pol_partition_123_05.')
	tmp	<- sapply(tmp, function(x) submitted.info[, which((grepl('IQTree150814/', FILE, fixed=1) | grepl('IQTree150818/', FILE, fixed=1)) & grepl(x, FILE, fixed=1))] )
	set(submitted.info, tmp, 'BEST', 'Y')
	tmp	<- c(	'150701_Regional_TRAIN2_IQTree151019_partition_123_10', 			
				'150701_Regional_TRAIN3_IQTree151019_partition_123_03',	
				'150701_Regional_TRAIN4_IQTree151019_partition_123_10',	
				'150701_Regional_TRAIN5_IQTree151019_partition_123_01',
				'150701_Regional_TRAIN2_IQTree151019_pol_partition_123_08',
				'150701_Regional_TRAIN3_IQTree151019_pol_partition_123_08',	
				'150701_Regional_TRAIN4_IQTree151019_pol_partition_123_05',
				'150701_Regional_TRAIN5_IQTree151019_pol_partition_123_10')
	tmp	<- sapply(tmp, function(x) submitted.info[, which((grepl('IQTree151019', FILE, fixed=1)) & grepl(x, FILE, fixed=1))] )
	set(submitted.info, tmp, 'BEST', 'Y')
	tmp	<- c(	'150701_Regional_TRAIN1_IQTree160530_gag_partition_123_09', 			
				'150701_Regional_TRAIN2_IQTree160530_gag_partition_123_06',	
				'150701_Regional_TRAIN4_IQTree160530_gag_partition_123_02')
	tmp	<- sapply(tmp, function(x) submitted.info[, which((grepl('IQTREE_Update_gag', FILE, fixed=1)) & grepl(x, FILE, fixed=1))] )
	set(submitted.info, tmp, 'BEST', 'Y')
	#	for RUNGAPS_ExaML we only have one tree per gap coverage, so all are 'best'
	set(submitted.info, submitted.info[, which(TEAM=='RUNGAPS_ExaML')], 'BEST', 'Y')
	#	PhyML: read log likelihood	 
	lkl		<- data.table(FILE= list.files('~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/PhyML', pattern='*phyml_stats*', recursive=1, full.names=1))
	lkl		<- lkl[, {
				z		<- readLines(FILE)
				z		<- as.numeric( gsub('\\s','',gsub('Log-likelihood:','',gsub('^\\.','',z[ which( grepl('Log-likelihood', z) ) ]))) )
				list( LOGLKL=z )
			}, by='FILE']
	lkl[, SC:=NA_character_]
	tmp		<- lkl[, which(is.na(SC) & grepl('150701_regional_train[0-9]', FILE))]
	set(lkl, tmp, 'SC', lkl[tmp, regmatches(FILE, regexpr('150701_regional_train[0-9]',FILE))])
	tmp		<- lkl[, which(is.na(SC) & grepl('150701_vill_scenario-[A-Z]', FILE))]
	set(lkl, tmp, 'SC', lkl[tmp, regmatches(FILE, regexpr('150701_vill_scenario-[A-Z]',FILE))])
	set(lkl, NULL, 'SC', lkl[, toupper(SC)])
	lkl[, TEAM:= 'PhyML']
	lkl[, GENE:= NA_character_]
	set(lkl, lkl[, which(TEAM=='PhyML' & grepl('gag', FILE))], 'GENE', 'GAG')
	set(lkl, lkl[, which(TEAM=='PhyML' & grepl('pol', FILE))], 'GENE', 'POL')
	set(lkl, lkl[, which(TEAM=='PhyML' & grepl('gagpolenv', FILE))], 'GENE', 'GAG+POL+ENV')
	setkey(lkl, GENE, SC)
	tmp		<- lkl[, list( FILE_BEST=FILE[ which.max(LOGLKL)[1] ] ), by=c('GENE','SC','TEAM')]
	set(tmp, NULL, 'FILE_BEST', tmp[, gsub('phyml_stats','phyml_tree',FILE_BEST)])	
	set(submitted.info, submitted.info[, which(FILE%in%tmp$FILE_BEST)], 'BEST', 'Y')
	#
	#	set OTHER (ie old or some preliminary/unknown tree)
	#
	submitted.info[, OTHER:='N']
	#	all MetaPIGA trees in 'MetaPIGA_150831' are old
	set(submitted.info, submitted.info[, which(TEAM=='MetaPIGA' & !grepl('MetaPIGA',FILE))], 'OTHER', 'Y')
	#	IQTree did several uploads, use only most recent in main analysis
	set(submitted.info, submitted.info[, which(grepl('150701_Regional_TRAIN1_IQTree150814', FILE))], 'OTHER', 'Y')
	set(submitted.info, submitted.info[, which(TEAM=='IQTree' & MODEL=='R' & !grepl('TRAIN1', SC) & grepl('201507/',FILE,fixed=1))], 'OTHER', 'Y')
	#	RAxML gag_1606 are old
	set(submitted.info, submitted.info[, which(TEAM=='RAXML' & grepl('gag_gene_1606',FILE))], 'OTHER', 'Y')
	#
	#	add BRL_UNITS
	#
	submitted.info[, BRL:='subst']
	#
	#	number taxa in tree
	#
	setkey(submitted.info, IDX)
	submitted.info[, TAXAN:= sapply(strs, Ntip)]	
	#
	#	add index of true tree
	#
	require(phangorn)
	tmp				<- subset(tfiles, select=c('SC','IDX_T','BRL_T'))
	tmp				<- dcast.data.table(tmp, SC~BRL_T, value.var='IDX_T')
	setnames(tmp, c('subst','time'), c("SUB_IDX_T","TIME_IDX_T"))
	submitted.info	<- merge(submitted.info, tmp, by='SC')
	submitted.info	<- merge(submitted.info, unique(subset(tfiles, select=c('SC','TAXAN_T'))), by='SC')	
	stopifnot(nrow(subset(submitted.info, TAXAN>TAXAN_T))==0)
	#
	#	fix taxa names that teams have changed
	#
	tmp		<- subset(submitted.info, TEAM=='IQTree' & MODEL=='R')[, IDX]
	for(i in tmp)
	{
		strs[[i]]$tip.label		<- sapply(strsplit(strs[[i]]$tip.label,'_'), function(x)	paste(x[1],'_',x[2],'|',x[3],'|',x[4],'_',x[5],'|',x[6],sep='')	)		
	}
	for(i in seq_along(strs))
	{
		strs[[i]]$tip.label		<- toupper(strs[[i]]$tip.label)		
	}
	for(i in seq_along(ttrs))
	{
		ttrs[[i]]$tip.label	<- toupper(ttrs[[i]]$tip.label)
	}		
	tmp2	<- subset(tinfo, BRL_T=='time', select=c(IDPOP,SC,TAXA))
	setkey(tmp2, IDPOP,SC,TAXA)
	tmp2	<- unique(tmp2)
	tmp		<- subset(submitted.info, TEAM=='PhyML' & MODEL=='R')[, IDX]
	for(i in tmp)
	{		
		cat(i,'\n')
		z	<- data.table(IDX=seq_along(strs[[i]]$tip.label), IDPOP=regmatches(strs[[i]]$tip.label, regexpr('IDPOP_[0-9]+',strs[[i]]$tip.label)), SC=subset(submitted.info, IDX==i)[,SC])
		z	<- merge(tmp2, z, by=c('IDPOP','SC'))
		setkey(z, IDX)
		stopifnot(nrow(z)==Ntip(strs[[i]]))
		strs[[i]]$tip.label	<- z[, TAXA]			
	}
	tmp		<- subset(submitted.info, TEAM=='PhyML' & MODEL=='V')[, IDX]
	for(i in tmp)
	{
		cat(i,'\n')
		z	<- data.table(IDX=seq_along(strs[[i]]$tip.label), IDPOP=regmatches(strs[[i]]$tip.label, regexpr('HOUSE[0-9]+-[0-9]+|House[0-9]+-[0-9]+',strs[[i]]$tip.label)), SC=subset(submitted.info, IDX==i)[,SC])
		z	<- merge(tmp2, z, by=c('IDPOP','SC'))
		stopifnot(nrow(z)==length(strs[[i]]$tip.label))
		setkey(z, IDX)
		strs[[i]]$tip.label	<- z[, TAXA]		
	}
	#
	#	check labels and remove labels that do not appear in the observed tree
	#
	tmp	<- submitted.info[, {				
				stree		<- unroot(strs[[IDX]])
				otree		<- unroot(ttrs[[TIME_IDX_T]])				
				z			<- setdiff(otree$tip.label, stree$tip.label)
				list(CHECK= length(z)==abs(diff(c(TAXAN, TAXAN_T))) )
			}, by='IDX']
	tmp	<- merge(subset(tmp, !CHECK), submitted.info, by='IDX')
	for(i in seq_len(nrow(tmp)))
	{
		j			<- tmp[i, IDX]
		cat('\n',j)		
		otree		<- tmp[i, TIME_IDX_T]
		stree		<- unroot(strs[[j]])
		otree		<- unroot(ttrs[[otree]])					
		z			<- merge( data.table(TAXA=stree$tip.label, TYPE='s'), data.table(TAXA=otree$tip.label, TYPE='o'), by='TAXA', all=1)
		z[, IDPOP:= gsub('IDPOP_','',regmatches(TAXA, regexpr('IDPOP_[0-9]+',TAXA)))]	
		strs[[j]]	<- drop.tip(stree, subset( z, is.na(TYPE.y))[, TAXA])		
	}
	submitted.info[, TAXAN:= sapply(strs, Ntip)]
	#
	#	are trees rooted?
	#
	setkey(submitted.info, IDX)
	submitted.info[, ROOTED:=factor(sapply(strs, is.rooted),levels=c(TRUE,FALSE),labels=c('Y','N'))]
	#
	#	read LSD trees
	#
	indir			<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/LSD'
	infiles			<- data.table(FILE=list.files(indir, pattern='LSD.date.newick', full.names=TRUE))
	infiles[, IDX:= as.integer(gsub('IDX_','',regmatches(basename(FILE),regexpr('IDX_[0-9]+',basename(FILE)))))]
	setkey(infiles, IDX)
	strs_lsd		<- vector('list', submitted.info[, max(IDX)])
	for(i in seq_len(nrow(infiles)))		
	{
		#	i<- 439
		IDX						<- infiles[i,IDX]
		FILE					<- infiles[i,FILE]
		cat('\n',IDX)
		ph						<- read.tree(FILE)
		stopifnot( !is.null(ph) )
		stopifnot( identical(sort(strs_rtt[[IDX]]$tip.label), sort(ph$tip.label)) )
		strs_lsd[[IDX]]			<- ph
		#names(strs_lsd[[IDX]])	<- FILE					
	}
	setkey(submitted.info, IDX)
	submitted.info[, WITH_LSD:= factor(sapply(strs_lsd, is.null), levels=c(TRUE,FALSE), labels=c('N','Y'))]	 	
	#
	#	SAVE so far
	#
	outdir		<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation'	
	save(strs, strs_lsd, ttrs, tinfo, tfiles, tbrl, submitted.info, file=file.path(outdir,'submitted_160713.rda'))	
}
##--------------------------------------------------------------------------------------------------------
##	olli 27.06.16
##--------------------------------------------------------------------------------------------------------
treecomparison.submissions.160627.addLSDtrees<- function()	
{
	require(data.table)
	require(ape)
	require(phangorn)
	wfile		<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation/submitted_160713_05QD.rda'
	indir		<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/LSD'
	
	load(wfile)
	infiles			<- data.table(FILE=list.files(indir, pattern='LSD.date.newick', full.names=TRUE))
	infiles[, IDX:= as.integer(gsub('IDX_','',regmatches(basename(FILE),regexpr('IDX_[0-9]+',basename(FILE)))))]
	setkey(infiles, IDX)
	
	strs_lsd		<- vector('list', submitted.info[, max(IDX)])
	for(i in seq_len(nrow(infiles)))		
	{
		#	i<- 439
		IDX						<- infiles[i,IDX]
		FILE					<- infiles[i,FILE]
		cat('\n',IDX)
		ph						<- read.tree(FILE)
		stopifnot( !is.null(ph) )
		stopifnot( identical(sort(strs_rtt[[IDX]]$tip.label), sort(ph$tip.label)) )
		strs_lsd[[IDX]]			<- ph
		#names(strs_lsd[[IDX]])	<- FILE					
	}
	setkey(submitted.info, IDX)
	submitted.info[, WITH_LSD:= factor(sapply(strs_lsd, is.null), levels=c(TRUE,FALSE), labels=c('N','Y'))]
	
	save(strs, strs_rtt, strs_lsd, ttrs, tbrl, tinfo, tfiles, submitted.info, sclu.info, lba, file=gsub('_05QD.rda','_06LSD.rda',wfile))
	
	
	str1		<- unroot(strs_rtt[[1]])
	str2		<- unroot(strs_lsd[[1]])				
	z			<- setdiff(str1$tip.label, str2$tip.label)				
	RF.dist(str1, str2, check.labels=TRUE)
	
}
##--------------------------------------------------------------------------------------------------------
##	olli 27.06.16
##--------------------------------------------------------------------------------------------------------
treecomparison.submissions.160627.stuffoncluster<- function(file)	
{
	require(data.table)
	require(ape)
	require(adephylo)
	require(phangorn)
	
	load(file)	
	#
	#	re-root simulated trees with rtt
	#		
	options(show.error.messages = FALSE)		
	readAttempt		<- try(suppressWarnings(load(gsub('.rda','_01rerooted.rda',file))))
	options(show.error.messages = TRUE)			
	if( 0 & inherits(readAttempt, "try-error") )
	{
		options(warn=2)
		strs_rtt	<- lapply(seq_along(strs), function(i)
				{
					cat('\n',i)
					#i	<- 628 ; i<- 241; i<- 571
					ph	<- strs[[i]]
					tmp	<- data.table(TAXA=ph$tip.label)
					set(tmp, NULL, 'T_SEQ', tmp[, as.numeric(regmatches(TAXA, regexpr('[0-9]*\\.[0-9]+$|[0-9]+$', TAXA))) ])					
					#phr	<- rtt(ph, tmp[, as.numeric(T_SEQ)])
					phr	<- rtt(ph, tmp[, T_SEQ], ncpu=1)	#this may drop a tip!
					phr
				})
		names(strs_rtt)	<- names(strs)
		options(warn=0)
		#
		#	ladderize all trees
		#		
		ttrs	<- lapply(ttrs, ladderize)
		strs	<- lapply(strs, ladderize)
		strs_rtt<- lapply(strs_rtt, ladderize)	
		#	save intermediate	
		save(strs, strs_rtt, strs_lsd, ttrs, tinfo, tbrl, tfiles, submitted.info, file=gsub('.rda','_01rerooted.rda',file))
	}	#			
	options(show.error.messages = FALSE)		
	readAttempt		<- try(suppressWarnings(load(gsub('.rda','_03RF.rda',file))))
	options(show.error.messages = TRUE)			
	if( 0 & inherits(readAttempt, "try-error") )
	{		
		#
		#	long branch attraction
		#
		lba		<- submitted.info[, {
					#IDX<-638; SUB_IDX_T<- 1
					cat(IDX,'\n')
					ph	<- strs_rtt[[IDX]]
					tmp	<- data.table(DEPTH=node.depth.edgelength(ph)[seq_len(Ntip(ph))], TAXA=ph$tip.label)
					tmp	<- merge(tmp, unique(subset(tinfo, IDX_T==SUB_IDX_T, c(TAXA,DEPTH_T))), by='TAXA')
					tmp
				}, by=c('MODEL','SC','TEAM','GAPS','GENE','IDX')]
		#
		# compute on true trees the proportion if either transmitter or among recipients
		#
		tmp				<- treedist.closest.ind.obs(tinfo, gd.thresh=Inf)
		setnames(tmp, 'TR_REC_perc_T', 'TR_REC_perc_T_Inf')	
		submitted.info	<- merge(submitted.info, tmp, by='SUB_IDX_T', all.x=1)
		tmp				<- treedist.closest.ind.obs(tinfo, gd.thresh=0.045)
		setnames(tmp, 'TR_REC_perc_T', 'TR_REC_perc_T_45')			
		submitted.info	<- merge(submitted.info, tmp, by='SUB_IDX_T', all.x=1)
		tmp				<- treedist.closest.ind.obs(tinfo, gd.thresh=0.01)
		setnames(tmp, c('TPAIR_PHCL_T','NTPAIR_PHCL_T'), c('TPAIR_PHCL_T_1','NTPAIR_PHCL_T_1'))			
		submitted.info	<- merge(submitted.info, tmp, by='SUB_IDX_T', all.x=1)				
		tinfo.pairs		<- treedist.closest.ind.obs(tinfo, gd.thresh=Inf, rtn.pairs=TRUE)
		setnames(tinfo.pairs, 'TRUE_PAIR','TRUE_PAIR_Inf')
		#
		# compute closest individual on simulated trees and determine proportion if either transmitter or among recipients
		#	
		tmp				<- treedist.closest.ind.reconstructed(submitted.info, tinfo, strs, gd.thresh=Inf)
		setnames(tmp, 'TR_REC_perc', 'TR_REC_perc_Inf')		
		submitted.info	<- merge(submitted.info, tmp, by='IDX', all.x=1)	
		tmp				<- treedist.closest.ind.reconstructed(submitted.info, tinfo, strs, gd.thresh=0.045)
		setnames(tmp, 'TR_REC_perc', 'TR_REC_perc_45')		
		submitted.info	<- merge(submitted.info, tmp, by='IDX', all.x=1)
		tmp				<- treedist.closest.ind.reconstructed(submitted.info, tinfo, strs, gd.thresh=0.01)
		setnames(tmp, c('TPAIR_PHCL', 'NTPAIR_PHCL'), c('TPAIR_PHCL_1', 'NTPAIR_PHCL_1'))		
		submitted.info	<- merge(submitted.info, tmp, by='IDX', all.x=1)			
		tmp				<- treedist.closest.ind.reconstructed.oftruepairs(submitted.info, tinfo.pairs, strs)
		submitted.info	<- merge(submitted.info, tmp, by='IDX', all.x=1)
		#
		#	compute Robinson Fould of complete tree
		#
		tmp				<- treedist.robinsonfould.wrapper(submitted.info, ttrs, strs_rtt)
		submitted.info	<- merge(submitted.info, tmp, by='IDX')
		#	compute Robinson Fould of clusters, then take sum
		tmp				<- treedist.robinsonfouldclusters.wrapper(submitted.info, ttrs, strs_rtt, tinfo)
		sclu.info		<- merge(subset(submitted.info, select=c("IDX","SC","FILE","TEAM","MODEL","SEQCOV","ACUTE","GAPS","ART","EXT","BEST","OTHER","GENE","TAXAN","ROOTED","BRL","SUB_IDX_T","TIME_IDX_T","TAXAN_T")), tmp, by='IDX')	
		#	save intermediate	
		save(strs, strs_rtt, strs_lsd, ttrs, tinfo, tbrl, tfiles, submitted.info, sclu.info, lba, file=gsub('.rda','_03RF.rda',file))
	}
	
	options(show.error.messages = FALSE)		
	readAttempt		<- try(suppressWarnings(load(gsub('.rda','_04PD.rda',file))))
	options(show.error.messages = TRUE)			
	if( 0 & inherits(readAttempt, "try-error") )
	{	
		#
		#	path distance of complete trees
		#
		cat('\nPath distances on rooted trees')
		tmp				<- treedist.pathdifference.wrapper(submitted.info, ttrs, strs_rtt, use.weight=FALSE)
		tmp[, TAXA_NJ:=NULL]
		submitted.info	<- merge(submitted.info, tmp, by='IDX', all.x=1)
		#	path distance of clusters	
		cat('\nPath distances on clusters')
		tmp				<- subset(submitted.info, MODEL=='R')
		tmp				<- treedist.pathdifference.clusters.wrapper(tmp, ttrs, strs_rtt, tinfo, use.weight=FALSE)
		tmp[, TAXA_NC:=NULL]
		sclu.info		<- merge(sclu.info, tmp, by=c('IDX','IDCLU'), all.x=1)
		#	save intermediate	
		save(strs, strs_rtt, strs_lsd, ttrs, tinfo, tbrl, tfiles, submitted.info, sclu.info, lba, file=gsub('.rda','_04PD.rda',file))
	}
	
	options(show.error.messages = FALSE)		
	readAttempt		<- try(suppressWarnings(load(gsub('.rda','_05QD.rda',file))))
	options(show.error.messages = TRUE)			
	if( 0 & inherits(readAttempt, "try-error") )
	{		
		#
		#	quartet distances of complete trees
		#
		cat('\nQuartett distances on rooted trees')
		tmp				<- treedist.quartetdifference.wrapper(submitted.info, ttrs, strs_rtt)
		tmp[, TAXA_NJ:=NULL]
		submitted.info	<- merge(submitted.info, tmp, by='IDX', all.x=1)
		#	quartet distance of clusters
		cat('\nQuartett distances on clusters')
		tmp				<- treedist.quartetdifference.clusters.wrapper(submitted.info, ttrs, strs_rtt, tinfo)
		tmp[, TAXA_NC:=NULL]
		sclu.info		<- merge(sclu.info, tmp, by=c('IDX','IDCLU'), all.x=1)
		save(strs, strs_rtt, strs_lsd, ttrs, tinfo, tbrl, tfiles, submitted.info, sclu.info, lba, file=gsub('.rda','_05QD.rda',file))
	}
	
	options(show.error.messages = FALSE)		
	readAttempt		<- try(suppressWarnings(load(gsub('.rda','_06PDLSD.rda',file))))
	options(show.error.messages = TRUE)			
	if( inherits(readAttempt, "try-error") )
	{	
		#
		#	path distance of complete LSD trees
		#
		cat('\nPath distances on LSD trees')
		tmp				<- subset(submitted.info, WITH_LSD=='Y')
		tmp				<- treedist.pathdifference.wrapper(tmp, ttrs, strs_lsd, use.brl=FALSE, use.weight=TRUE)
		setnames(tmp, c('PD','NPD','NPDSQ'), c('PD_LSD','NPD_LSD','NPDSQ_LSD'))		
		tmp[, TAXA_NJ:=NULL]
		submitted.info	<- merge(submitted.info, tmp, by='IDX', all.x=1)
		#	path distance of LSD clusters
		cat('\nPath distances on LSD clusters')
		tmp				<- subset(submitted.info, WITH_LSD=='Y' & MODEL=='R')
		tmp				<- treedist.pathdifference.clusters.wrapper(tmp, ttrs, strs_lsd, tinfo, use.brl=FALSE, use.weight=TRUE)
		setnames(tmp, c('PD','NPD','NPDSQ'), c('PD_LSD','NPD_LSD','NPDSQ_LSD'))
		tmp[, TAXA_NC:=NULL]
		sclu.info		<- merge(sclu.info, tmp, by=c('IDX','IDCLU'), all.x=1)
		#	save intermediate	
		save(strs, strs_rtt, strs_lsd, ttrs, tinfo, tbrl, tfiles, submitted.info, sclu.info, lba, file=gsub('.rda','_06PDLSD.rda',file))
	}
	
	options(show.error.messages = FALSE)		
	readAttempt		<- try(suppressWarnings(load(gsub('.rda','_07MSELSD.rda',file))))
	options(show.error.messages = TRUE)			
	if( inherits(readAttempt, "try-error") )
	{			
		#	MSE between true time distances and reconstructed patristic distances in LSD tree
		cat('\nMSE of edges on LSD trees')
		tmp				<- subset(submitted.info, WITH_LSD=='Y')
		tmp				<- treedist.MSE.wrapper(tmp, strs_lsd, tbrl, tinfo, use.brl=FALSE)
		setnames(tmp, c('MSE','MAE','MSE_TP','MAE_TP'), c('MSE_LSD','MAE_LSD','MSE_TP_LSD','MAE_TP_LSD'))
		tmp[, TAXA_NJ:=NULL]
		submitted.info	<- merge(submitted.info, tmp, by='IDX', all.x=1)
		cat('\nMSE of edges on LSD clusters')
		tmp				<- subset(submitted.info, WITH_LSD=='Y' & MODEL=='R')
		tmp				<- treedist.MSE.clusters.wrapper(tmp, strs_lsd, tbrl, tinfo, use.brl=FALSE)
		setnames(tmp, c('MSE','MAE','MSE_TP','MAE_TP'), c('MSE_LSD','MAE_LSD','MSE_TP_LSD','MAE_TP_LSD'))
		tmp[, TAXA_NC:=NULL]
		sclu.info		<- merge(sclu.info, tmp, by=c('IDX','IDCLU'), all.x=1)
		#	save intermediate	
		save(strs, strs_rtt, strs_lsd, ttrs, tinfo, tbrl, tfiles, submitted.info, sclu.info, lba, file=gsub('.rda','_07MSELSD.rda',file))
	}
	
	options(show.error.messages = FALSE)		
	readAttempt		<- try(suppressWarnings(load(gsub('.rda','_08MSEGD.rda',file))))
	options(show.error.messages = TRUE)			
	if( inherits(readAttempt, "try-error") )
	{			
		#	MSE between true time distances and reconstructed patristic distances in LSD tree
		cat('\nMSE of edges on SUB trees')
		tmp				<- subset(submitted.info, MODEL=='R')
		tmp				<- treedist.MSE.wrapper(tmp, strs_rtt, tbrl, tinfo, use.brl=TRUE)
		setnames(tmp, c('MSE','MAE','MSE_TP','MAE_TP'), c('MSE_GD','MAE_GD','MSE_TP_GD','MAE_TP_GD'))
		tmp[, TAXA_NJ:=NULL]
		tmp[, EDGE_NJ:=NULL]
		submitted.info	<- merge(submitted.info, tmp, by='IDX', all.x=1)
		cat('\nMSE of edges on SUB clusters')
		tmp				<- subset(submitted.info, MODEL=='R')
		tmp				<- treedist.MSE.clusters.wrapper(tmp, strs_rtt, tbrl, tinfo, use.brl=TRUE)
		setnames(tmp, c('MSE','MAE','MSE_TP','MAE_TP'), c('MSE_GD','MAE_GD','MSE_TP_GD','MAE_TP_GD'))
		tmp[, TAXA_NC:=NULL]
		tmp[, EDGE_NC:=NULL]
		sclu.info		<- merge(sclu.info, tmp, by=c('IDX','IDCLU'), all.x=1)
		#	save intermediate	
		save(strs, strs_rtt, strs_lsd, ttrs, tinfo, tbrl, tfiles, submitted.info, sclu.info, lba, file=gsub('.rda','_07MSELSD.rda',file))
	}
	options(show.error.messages = FALSE)		
	readAttempt		<- try(suppressWarnings(load(gsub('.rda','_09SBRL.rda',file))))
	options(show.error.messages = TRUE)				
	if( inherits(readAttempt, "try-error") )
	{		
		#
		#	sum of branch lengths per tree
		#
		tmp			<- unique(subset(tinfo, select=IDX_T))
		tmp			<- tmp[, {
								ph	<- ttrs[[IDX_T]]
								list(SUM_BRANCHES_T=sum(ph$edge.length))					
							}, by='IDX_T']
		setnames(tmp, 'IDX_T', 'SUB_IDX_T')			
		submitted.info	<- merge(submitted.info, tmp, by='SUB_IDX_T')
		
		tmp		<- submitted.info[, {
					#IDX<-638; SUB_IDX_T<- 1
					#cat(IDX,'\n')
					ph	<- strs[[IDX]]
					list(SUM_BRANCHES=sum(ph$edge.length))					
				}, by=c('IDX')]
		submitted.info	<- merge(submitted.info, tmp, by='IDX')
		#	save intermediate	
		save(strs, strs_rtt, strs_lsd, ttrs, tinfo, tbrl, tfiles, submitted.info, sclu.info, lba, file=gsub('.rda','_09SBRL.rda',file))
	}
	#
	#	ADD other summaries
	#
	#load( file.path(outdir, 'submitted_160704_KC.rda') )
	#sclu.info.kc	<- copy(sclu.info)
	#load( file.path(outdir, 'submitted_160627_QDPD.rda') )
	#sclu.info		<- merge(sclu.info, subset(sclu.info.kc, select=c(IDX, TEAM, GENE, BRL, IDCLU, KC)), by=c('IDX','TEAM','GENE','BRL','IDCLU'))
	#save(strs, strs_rtt, ttrs, tinfo, submitted.info, sclu.info, lba,  file=file.path(outdir,'submitted_160627_QDPDKC.rda'))
}
#
#
#
treecomparison.submissions.151101<- function()	
{
	require(data.table)
	require(ape)
	require(phangorn)
	#
	#	get true trees
	#
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim_internal/freeze_July15'
	tfiles	<- list.files(indir, pattern='newick$', full.names=TRUE)
	tfiles	<- data.table( FILE_T=tfiles[ grepl('SUBSTTREE', tfiles) | grepl('Vill_99', tfiles) | grepl('Vill.*DATEDTREE', tfiles) ] )
	tfiles[, SC:= toupper(gsub('_SUBSTTREE|_DATEDTREE','',gsub('.newick','',basename(FILE_T))))]
	tmp		<- rbind( subset(tfiles, SC=='VILL_99_APR15'), subset(tfiles, SC=='VILL_99_APR15'), subset(tfiles, SC=='VILL_99_APR15') )
	set(tmp, NULL, 'SC', c('150701_VILL_SCENARIO-C','150701_VILL_SCENARIO-D','150701_VILL_SCENARIO-E'))
	tfiles	<- rbind(tfiles, tmp)
	ttrs	<- lapply(tfiles[, FILE_T], function(x)	read.tree(file=x) )
	names(ttrs)	<- tfiles[, SC]
	tfiles[, IDX_T:=seq_along(ttrs)]
	tfiles[, TAXAN_T:= sapply(ttrs, Ntip)]
	#	info on true trees
	tinfo	<- merge(tfiles, do.call('rbind',lapply(seq_along(ttrs), function(i) data.table(TAXA=ttrs[[i]]$tip.label, IDX_T=i))), by='IDX_T')	
	tinfo[, IDPOP:=NA_character_]
	tmp		<- tinfo[, which(grepl('REGIONAL',SC))]
	set(tinfo, tmp, 'IDPOP', tinfo[tmp,regmatches(TAXA, regexpr('IDPOP_[0-9]+',TAXA))])
	tmp		<- tinfo[, which(!grepl('REGIONAL',SC))]
	set(tinfo, tmp, 'IDPOP', tinfo[tmp, regmatches(TAXA, regexpr('HOUSE[0-9]+-[0-9]+|House[0-9]+-[0-9]+',TAXA))])		
	stopifnot(subset(tinfo, grepl('VILL',SC))[, length(which(substring(TAXA,1,10)!=substring(IDPOP,1,10)))]==0)	
	stopifnot( tinfo[, length(which(is.na(IDPOP)))==0] )	
	set(tinfo, NULL, 'IDPOP', tinfo[,toupper(IDPOP)])
	set(tinfo, NULL, 'TAXA', tinfo[,toupper(TAXA)])
	#	read cluster membership from DATEDCLUTREES	
	tmp		<- list.files(indir, pattern='DATEDCLUTREES', full.names=TRUE)
	tmp		<- data.table( FILE_CLU_T= tmp, SC= toupper(gsub('_DATEDCLUTREES','',gsub('.newick','',basename(tmp))))) 
	tfiles	<- merge(tfiles, tmp, by='SC', all=1)	
	tmp		<- subset(tfiles, !is.na(FILE_CLU_T))[, {
				z		<- read.tree(FILE_CLU_T)
				do.call('rbind',lapply(seq_along(z), function(i) data.table(IDCLU=i, TAXA=z[[i]]$tip.label)))				
			}, by='SC']	
	tinfo	<- merge(tinfo, tmp, by=c('SC','TAXA'), all=1)
	tmp		<- subset(tinfo, !is.na(IDCLU))[, list(CLU_N= length(IDPOP)), by=c('SC','IDCLU')]
	tinfo	<- merge(tinfo, tmp, by=c('SC','IDCLU'), all=1)
	#	read sequences and determine %gappiness
	tmp		<- list.files(indir, pattern='fa$|fasta$', full.names=TRUE)
	tmp		<- data.table( FILE_SEQ_T= tmp, SC= toupper(gsub('_SIMULATED','',gsub('.fa','',basename(tmp)))))
	z		<- subset(tmp, SC=='VILL_99_APR15')
	set(z, NULL, 'SC', '150701_VILL_SCENARIO-C')
	tmp		<- rbind( tmp, z )	
	tfiles	<- merge(tfiles, tmp, by='SC', all=1)
	tmp		<- subset(tfiles, !is.na(FILE_SEQ_T))[, {
				z		<- read.dna(FILE_SEQ_T, format='fasta')	
				ans		<- sapply(seq_len(nrow(z)), function(i) base.freq(z[i,], all=1))
				ans		<- apply(ans[c('n','-','?'),], 2, sum)
				list(TAXA=rownames(z), GPS=ans)				
			}, by='SC']
	tinfo	<- merge(tinfo, tmp, by=c('SC','TAXA'), all.x=1)
	#
	#	get submitted trees
	#	
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/IQTree/IQTree201507'
	infiles	<- list.files(indir, pattern='treefile$', recursive=1, full.names=1)
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/IQTree/IQTree201510'
	infiles	<- c(infiles, list.files(indir, pattern='treefile$', recursive=1, full.names=1))	
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/PhyML'
	infiles	<- c(infiles, list.files(indir, pattern='*tree*', recursive=1, full.names=1))
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/RAxML'
	infiles	<- c(infiles, list.files(indir, pattern='*RAxML_bestTree*', recursive=1, full.names=1))
	infiles	<- c(infiles, list.files(indir, pattern="best_tree.newick", recursive=1, full.names=1))
	infiles	<- data.table(FILE=infiles)
	strs	<- lapply(infiles[, FILE], function(x)
			{
				cat(x)
				read.tree(file=x)	
			})
	names(strs)	<- infiles[, FILE]
	
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/MetaPIGA'
	tmp		<-  list.files(indir, pattern='*result*', recursive=1, full.names=1)
	tmp		<- data.table(FILE=tmp)	
	tmp.trees			<- lapply(tmp[, FILE], function(x)
			{
				cat(x)
				read.nexus(file=x)	
			})
	sapply(tmp.trees, length)
	MetaPIGA.trees			<- c(lapply(tmp.trees, '[[', 1), lapply(tmp.trees, '[[', 2), lapply(tmp.trees, '[[', 3), lapply(tmp.trees, '[[', 4))
	names(MetaPIGA.trees)	<- c(sapply(tmp.trees, function(x) paste(names(x)[1],'_use',sep='')), sapply(tmp.trees, function(x) names(x)[2]), sapply(tmp.trees, function(x) names(x)[3]), sapply(tmp.trees, function(x) names(x)[4]))	
	names(MetaPIGA.trees)	<- gsub("'",'',names(MetaPIGA.trees), fixed=1)	
	strs					<- c(strs, MetaPIGA.trees)	
	submitted.info			<- data.table(FILE=names(strs))
	#
	#
	#	
	submitted.info[, IDX:=seq_along(strs)]	
	submitted.info[, TEAM:=NA_character_]
	set(submitted.info, submitted.info[, which(grepl('RAXML|RAxML',FILE))], 'TEAM', 'RAXML')
	set(submitted.info, submitted.info[, which(grepl('IQTree',FILE))], 'TEAM', 'IQTree')
	set(submitted.info, submitted.info[, which(grepl('MetaPIGA|Consensus pruning|Best individual of population',FILE))], 'TEAM', 'MetaPIGA')
	set(submitted.info, submitted.info[, which(grepl('PhyML',FILE))], 'TEAM', 'PhyML')	
	stopifnot( submitted.info[, length(which(is.na(TEAM)))==0] )
	#
	#	scenario
	#	
	submitted.info[, SC:=NA_character_]
	tmp		<- submitted.info[, which(grepl('150701_Regional_TRAIN[0-9]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, regmatches(FILE, regexpr('150701_Regional_TRAIN[0-9]',FILE))])
	tmp		<- submitted.info[, which(grepl('150701_Vill_SCENARIO-[A-Z]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, regmatches(FILE, regexpr('150701_Vill_SCENARIO-[A-Z]',FILE))])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('TRAIN[0-9]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, paste('150701_Regional_',regmatches(FILE, regexpr('TRAIN[0-9]',FILE)),sep='')])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('scenario[A-Z]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, paste('150701_Vill_',regmatches(FILE, regexpr('scenario[A-Z]',FILE)),sep='')])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('150701_regional_train[0-9]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, regmatches(FILE, regexpr('150701_regional_train[0-9]',FILE))])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('150701_vill_scenario-[A-Z]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, regmatches(FILE, regexpr('150701_vill_scenario-[A-Z]',FILE))])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('Vill_99_Apr15', FILE))]
	set(submitted.info, tmp, 'SC', 'Vill_99_Apr15')	
	set(submitted.info, NULL, 'SC', submitted.info[, toupper(SC)])
	tmp		<- submitted.info[, which(grepl('150701_VILL_SCENARIO[A-Z]', SC))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, gsub('150701_VILL_SCENARIO','150701_VILL_SCENARIO-',SC)])
	stopifnot( submitted.info[, length(which(is.na(SC)))==0] )
	#
	#	set covariates of scenarios
	#
	tmp		<- data.table(	SC=		c("150701_REGIONAL_TRAIN1","150701_REGIONAL_TRAIN2","150701_REGIONAL_TRAIN3","150701_REGIONAL_TRAIN4" ,"150701_REGIONAL_TRAIN5", "150701_VILL_SCENARIO-A", "150701_VILL_SCENARIO-B", "VILL_99_APR15","150701_VILL_SCENARIO-C", "150701_VILL_SCENARIO-D", "150701_VILL_SCENARIO-E"),
			MODEL=	c('R','R','R','R','R','V','V','V','V','V','V'),
			SEQCOV= c(0.16, 0.16, 0.16, 0.16, 0.16, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6),
			ACUTE=	c('low', 'low', 'high', 'low', 'high', 'high', 'high', 'high', 'high', 'high', 'high'),
			GAPS=	c('none', 'low', 'low', 'high', 'high', 'low', 'high', 'none', 'none', 'low', 'high'), 
			ART=	c('none', 'none', 'none', 'none', 'none', 'none', 'none', 'fast', 'fast', 'fast', 'fast'),
			EXT= 	c('5pc', '5pc', '5pc', '5pc', '5pc', '~0pc', '~0pc', '~0pc', '~0pc', '~0pc', '~0pc')
	)
	submitted.info	<- merge(submitted.info, tmp, by='SC')
	#
	#	best tree for each scenario
	#
	submitted.info[, BEST:='N']
	set(submitted.info, submitted.info[, which(grepl('RAxML', FILE) & grepl('best_tree', FILE))], 'BEST', 'Y')	
	#	copied from ListOfBestTrees_IQTree150818.txt
	#	there are several best trees for some scenarios
	tmp	<- c( '150701_Vill_SCENARIO-A_IQTree150814_partition_12_3_07',	
			'150701_Vill_SCENARIO-A_IQTree150814_partition_12_3_04.',	
			'150701_Vill_SCENARIO-B_IQTree150814_partition_12_3_03.',	
			'Vill_99_Apr15_IQTree150814_partition_123.',			
			'150701_Vill_SCENARIO-D_IQTree150814_partition_12_3.',	
			'150701_Vill_SCENARIO-E_IQTree150814_partition_12_3.',
			'150701_Vill_SCENARIO-A_IQTree150814_pol_partition_12_3.',	
			'150701_Vill_SCENARIO-B_IQTree150814_pol_partition_12_3_05.',
			'Vill_99_Apr15_IQTree150814_pol_partition_12_3_09.',		
			'Vill_99_Apr15_IQTree150814_pol_partition_12_3_10.',		
			'150701_Vill_SCENARIO-D_IQTree150814_pol_partition_12_3_05.',
			'150701_Vill_SCENARIO-D_IQTree150814_pol_partition_12_3_06.',
			'150701_Vill_SCENARIO-D_IQTree150814_pol_partition_12_3_09.',
			'150701_Vill_SCENARIO-E_IQTree150814_pol_partition_12_3_06.',
			'150701_Regional_TRAIN1_IQTree150818_partition_123_03.',	
			'150701_Regional_TRAIN1_IQTree150818_pol_partition_123_05.')
	tmp	<- sapply(tmp, function(x) submitted.info[, which((grepl('IQTree150814/', FILE, fixed=1) | grepl('IQTree150818/', FILE, fixed=1)) & grepl(x, FILE, fixed=1))] )
	set(submitted.info, tmp, 'BEST', 'Y')
	tmp	<- c('150701_Regional_TRAIN2_IQTree151019_partition_123_10', 			
			'150701_Regional_TRAIN3_IQTree151019_partition_123_03',	
			'150701_Regional_TRAIN4_IQTree151019_partition_123_10',	
			'150701_Regional_TRAIN5_IQTree151019_partition_123_01',
			'150701_Regional_TRAIN2_IQTree151019_pol_partition_123_08',
			'150701_Regional_TRAIN3_IQTree151019_pol_partition_123_08',	
			'150701_Regional_TRAIN4_IQTree151019_pol_partition_123_05',
			'150701_Regional_TRAIN5_IQTree151019_pol_partition_123_10')
	tmp	<- sapply(tmp, function(x) submitted.info[, which((grepl('IQTree151019', FILE, fixed=1)) & grepl(x, FILE, fixed=1))] )
	set(submitted.info, tmp, 'BEST', 'Y')
	#	PhyML no replicates: all files best
	set(submitted.info, submitted.info[, which(TEAM=='PhyML')], 'BEST', 'Y')		
	#
	#	set OTHER (ie old or some preliminary/unknown tree)
	#
	submitted.info[, OTHER:='N']
	#	MetaPIGA tree to be used is first in nexus list (which was tagged with best above)
	set(submitted.info, submitted.info[, which(TEAM=='MetaPIGA' & !grepl('use', FILE))], 'OTHER', 'Y')
	#	IQTree did several uploads, use only most recent in main analysis
	set(submitted.info, submitted.info[, which(grepl('150701_Regional_TRAIN1_IQTree150814', FILE))], 'OTHER', 'Y')
	set(submitted.info, submitted.info[, which(TEAM=='IQTree' & MODEL=='R' & !grepl('TRAIN1', SC) & grepl('201507/',FILE,fixed=1))], 'OTHER', 'Y')
	#
	#	set which gene used to construct tree (either pol or concatenated gag+pol+env)
	#
	submitted.info[, GENE:=NA_character_]
	set(submitted.info, submitted.info[, which(TEAM=='RAXML' & grepl('full', FILE))], 'GENE', 'GAG+POL+ENV')
	set(submitted.info, submitted.info[, which(TEAM=='RAXML' & grepl('pol', FILE))], 'GENE', 'POL')
	stopifnot(nrow(subset(submitted.info, TEAM=='RAXML' & is.na(GENE)))==0)
	set(submitted.info, submitted.info[, which(TEAM=='PhyML')], 'GENE', 'POL')
	set(submitted.info, submitted.info[, which(TEAM=='MetaPIGA')], 'GENE', 'GAG+POL+ENV')	
	set(submitted.info, submitted.info[, which(TEAM=='IQTree' & grepl('[0-9]_partition', FILE))], 'GENE', 'GAG+POL+ENV')
	set(submitted.info, submitted.info[, which(TEAM=='IQTree' & grepl('[0-9]_pol_partition', FILE))], 'GENE', 'POL')
	stopifnot(nrow(subset(submitted.info, TEAM=='IQTree' & is.na(GENE)))==0)
	#
	#	number taxa in tree
	#
	setkey(submitted.info, IDX)
	submitted.info[, TAXAN:= sapply(strs, Ntip)]
	#
	#	are trees rooted?
	#
	setkey(submitted.info, IDX)
	submitted.info[, ROOTED:=factor(sapply(strs, is.rooted),levels=c(TRUE,FALSE),labels=c('Y','N'))]
	#
	#	add index of true tree
	#
	require(phangorn)
	submitted.info	<- merge(submitted.info, subset(tfiles, select=c('SC','IDX_T','TAXAN_T')), by='SC')
	stopifnot(nrow(subset(submitted.info, TAXAN>TAXAN_T))==0)
	#
	#	fix taxa names that teams have changed
	#
	tmp		<- subset(submitted.info, TEAM=='IQTree' & MODEL=='R')[, IDX]
	for(i in tmp)
	{
		strs[[i]]$tip.label	<- sapply(strsplit(strs[[i]]$tip.label,'_'), function(x)	paste(x[1],'_',x[2],'|',x[3],'|',x[4],'_',x[5],'|',x[6],sep='')	)
	}
	for(i in seq_along(strs))
	{
		strs[[i]]$tip.label	<- toupper(strs[[i]]$tip.label)
	}
	for(i in seq_along(ttrs))
	{
		ttrs[[i]]$tip.label	<- toupper(ttrs[[i]]$tip.label)
	}	
	###
	tmp		<- subset(submitted.info, TEAM=='PhyML' & MODEL=='R')[, IDX]
	for(i in tmp)
	{
		
		z	<- data.table(IDX=seq_along(strs[[i]]$tip.label), IDPOP=regmatches(strs[[i]]$tip.label, regexpr('IDPOP_[0-9]+',strs[[i]]$tip.label)), SC=subset(submitted.info, IDX==i)[,SC])
		z	<- merge(subset(tinfo, select=c(IDPOP,SC,TAXA)), z, by=c('IDPOP','SC'))
		setkey(z, IDX)
		strs[[i]]$tip.label	<- z[, TAXA]
	}
	tmp		<- subset(submitted.info, TEAM=='PhyML' & MODEL=='V')[, IDX]
	for(i in tmp)
	{
		
		z	<- data.table(IDX=seq_along(strs[[i]]$tip.label), IDPOP=regmatches(strs[[i]]$tip.label, regexpr('HOUSE[0-9]+-[0-9]+|House[0-9]+-[0-9]+',strs[[i]]$tip.label)), SC=subset(submitted.info, IDX==i)[,SC])
		z	<- merge(subset(tinfo, select=c(IDPOP,SC,TAXA)), z, by=c('IDPOP','SC'))
		stopifnot(nrow(z)==length(strs[[i]]$tip.label))
		setkey(z, IDX)
		strs[[i]]$tip.label	<- z[, TAXA]
	}
	#
	#	compute Robinson Fould of complete tree
	#
	tmp				<- treedist.robinsonfould.wrapper(submitted.info, ttrs, strs)
	submitted.info	<- merge(submitted.info, tmp, by='IDX')
	#	compute Robinson Fould of clusters, then take sum
	tmp				<- treedist.robinsonfouldclusters.wrapper(submitted.info, ttrs, strs, tinfo)
	sclu.info		<- merge(submitted.info, tmp, by='IDX')
	#
	outfile	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation/submitted_151101.rda'
	save(strs, ttrs, tinfo, submitted.info, sclu.info, file=outfile)
}
treecomparison.submissions.151016<- function()	
{
	require(data.table)
	require(ape)
	require(phangorn)
	#
	#	get true trees
	#
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim_internal/freeze_July15'
	tfiles	<- list.files(indir, pattern='newick$', full.names=TRUE)
	tfiles	<- data.table( FILE_T=tfiles[ grepl('SUBSTTREE', tfiles) | grepl('Vill_99', tfiles) | grepl('Vill.*DATEDTREE', tfiles) ] )
	tfiles[, SC:= toupper(gsub('_SUBSTTREE|_DATEDTREE','',gsub('.newick','',basename(FILE_T))))]
	tmp		<- rbind( subset(tfiles, SC=='VILL_99_APR15'), subset(tfiles, SC=='VILL_99_APR15'), subset(tfiles, SC=='VILL_99_APR15') )
	set(tmp, NULL, 'SC', c('150701_VILL_SCENARIO-C','150701_VILL_SCENARIO-D','150701_VILL_SCENARIO-E'))
	tfiles	<- rbind(tfiles, tmp)
	ttrs	<- lapply(tfiles[, FILE_T], function(x)	read.tree(file=x) )
	names(ttrs)	<- tfiles[, SC]
	tfiles[, IDX_T:=seq_along(ttrs)]
	tfiles[, TAXAN_T:= sapply(ttrs, Ntip)]
	#	info on true trees
	tinfo	<- merge(tfiles, do.call('rbind',lapply(seq_along(ttrs), function(i) data.table(TAXA=ttrs[[i]]$tip.label, IDX_T=i))), by='IDX_T')	
	tinfo[, IDPOP:=NA_character_]
	tmp		<- tinfo[, which(grepl('REGIONAL',SC))]
	set(tinfo, tmp, 'IDPOP', tinfo[tmp,regmatches(TAXA, regexpr('IDPOP_[0-9]+',TAXA))])
	tmp		<- tinfo[, which(!grepl('REGIONAL',SC))]
	set(tinfo, tmp, 'IDPOP', tinfo[tmp, regmatches(TAXA, regexpr('HOUSE[0-9]+-[0-9]+|House[0-9]+-[0-9]+',TAXA))])		
	stopifnot(subset(tinfo, grepl('VILL',SC))[, length(which(substring(TAXA,1,10)!=substring(IDPOP,1,10)))]==0)	
	stopifnot( tinfo[, length(which(is.na(IDPOP)))==0] )	
	set(tinfo, NULL, 'IDPOP', tinfo[,toupper(IDPOP)])
	set(tinfo, NULL, 'TAXA', tinfo[,toupper(TAXA)])
	#	read cluster membership from DATEDCLUTREES	
	tmp		<- list.files(indir, pattern='DATEDCLUTREES', full.names=TRUE)
	tmp		<- data.table( FILE_CLU_T= tmp, SC= toupper(gsub('_DATEDCLUTREES','',gsub('.newick','',basename(tmp))))) 
	tfiles	<- merge(tfiles, tmp, by='SC', all=1)	
	tmp		<- subset(tfiles, !is.na(FILE_CLU_T))[, {
														z		<- read.tree(FILE_CLU_T)
														do.call('rbind',lapply(seq_along(z), function(i) data.table(IDCLU=i, TAXA=z[[i]]$tip.label)))				
													}, by='SC']	
	tinfo	<- merge(tinfo, tmp, by=c('SC','TAXA'), all=1)
	tmp		<- subset(tinfo, !is.na(IDCLU))[, list(CLU_N= length(IDPOP)), by=c('SC','IDCLU')]
	tinfo	<- merge(tinfo, tmp, by=c('SC','IDCLU'), all=1)
	#	read sequences and determine %gappiness
	tmp		<- list.files(indir, pattern='fa$|fasta$', full.names=TRUE)
	tmp		<- data.table( FILE_SEQ_T= tmp, SC= toupper(gsub('_SIMULATED','',gsub('.fa','',basename(tmp)))))
	z		<- subset(tmp, SC=='VILL_99_APR15')
	set(z, NULL, 'SC', '150701_VILL_SCENARIO-C')
	tmp		<- rbind( tmp, z )	
	tfiles	<- merge(tfiles, tmp, by='SC', all=1)
	tmp		<- subset(tfiles, !is.na(FILE_SEQ_T))[, {
				z		<- read.dna(FILE_SEQ_T, format='fasta')	
				ans		<- sapply(seq_len(nrow(z)), function(i) base.freq(z[i,], all=1))
				ans		<- apply(ans[c('n','-','?'),], 2, sum)
				list(TAXA=rownames(z), GPS=ans)				
			}, by='SC']
	tinfo	<- merge(tinfo, tmp, by=c('SC','TAXA'), all.x=1)
	#
	#	get submitted trees
	#	
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/IQTree/IQTree201507'
	infiles	<- list.files(indir, pattern='treefile$', recursive=1, full.names=1)
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/PhyML'
	infiles	<- c(infiles, list.files(indir, pattern='*tree*', recursive=1, full.names=1))
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/RAxML'
	infiles	<- c(infiles, list.files(indir, pattern='*RAxML_bestTree*', recursive=1, full.names=1))
	infiles	<- c(infiles, list.files(indir, pattern="best_tree.newick", recursive=1, full.names=1))
	infiles	<- data.table(FILE=infiles)
	strs	<- lapply(infiles[, FILE], function(x)
			{
				cat(x)
				read.tree(file=x)	
			})
	names(strs)	<- infiles[, FILE]
	
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/MetaPIGA'
	tmp		<-  list.files(indir, pattern='*result*', recursive=1, full.names=1)
	tmp		<- data.table(FILE=tmp)	
	tmp.trees			<- lapply(tmp[, FILE], function(x)
			{
				cat(x)
				read.nexus(file=x)	
			})
	sapply(tmp.trees, length)
	MetaPIGA.trees			<- c(lapply(tmp.trees, '[[', 1), lapply(tmp.trees, '[[', 2), lapply(tmp.trees, '[[', 3), lapply(tmp.trees, '[[', 4))
	names(MetaPIGA.trees)	<- c(sapply(tmp.trees, function(x) paste(names(x)[1],'_use',sep='')), sapply(tmp.trees, function(x) names(x)[2]), sapply(tmp.trees, function(x) names(x)[3]), sapply(tmp.trees, function(x) names(x)[4]))	
	names(MetaPIGA.trees)	<- gsub("'",'',names(MetaPIGA.trees), fixed=1)	
	strs					<- c(strs, MetaPIGA.trees)	
	submitted.info			<- data.table(FILE=names(strs))
	#
	#
	#	
	submitted.info[, IDX:=seq_along(strs)]	
	submitted.info[, TEAM:=NA_character_]
	set(submitted.info, submitted.info[, which(grepl('RAXML|RAxML',FILE))], 'TEAM', 'RAXML')
	set(submitted.info, submitted.info[, which(grepl('IQTree',FILE))], 'TEAM', 'IQTree')
	set(submitted.info, submitted.info[, which(grepl('MetaPIGA|Consensus pruning|Best individual of population',FILE))], 'TEAM', 'MetaPIGA')
	set(submitted.info, submitted.info[, which(grepl('PhyML',FILE))], 'TEAM', 'PhyML')	
	stopifnot( submitted.info[, length(which(is.na(TEAM)))==0] )
	#
	#	scenario
	#	
	submitted.info[, SC:=NA_character_]
	tmp		<- submitted.info[, which(grepl('150701_Regional_TRAIN[0-9]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, regmatches(FILE, regexpr('150701_Regional_TRAIN[0-9]',FILE))])
	tmp		<- submitted.info[, which(grepl('150701_Vill_SCENARIO-[A-Z]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, regmatches(FILE, regexpr('150701_Vill_SCENARIO-[A-Z]',FILE))])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('TRAIN[0-9]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, paste('150701_Regional_',regmatches(FILE, regexpr('TRAIN[0-9]',FILE)),sep='')])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('scenario[A-Z]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, paste('150701_Vill_',regmatches(FILE, regexpr('scenario[A-Z]',FILE)),sep='')])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('150701_regional_train[0-9]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, regmatches(FILE, regexpr('150701_regional_train[0-9]',FILE))])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('150701_vill_scenario-[A-Z]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, regmatches(FILE, regexpr('150701_vill_scenario-[A-Z]',FILE))])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('Vill_99_Apr15', FILE))]
	set(submitted.info, tmp, 'SC', 'Vill_99_Apr15')	
	set(submitted.info, NULL, 'SC', submitted.info[, toupper(SC)])
	tmp		<- submitted.info[, which(grepl('150701_VILL_SCENARIO[A-Z]', SC))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, gsub('150701_VILL_SCENARIO','150701_VILL_SCENARIO-',SC)])
	stopifnot( submitted.info[, length(which(is.na(SC)))] )
	#
	#	set covariates of scenarios
	#
	tmp		<- data.table(	SC=		c("150701_REGIONAL_TRAIN1","150701_REGIONAL_TRAIN2","150701_REGIONAL_TRAIN3","150701_REGIONAL_TRAIN4" ,"150701_REGIONAL_TRAIN5", "150701_VILL_SCENARIO-A", "150701_VILL_SCENARIO-B", "VILL_99_APR15","150701_VILL_SCENARIO-C", "150701_VILL_SCENARIO-D", "150701_VILL_SCENARIO-E"),
				MODEL=	c('R','R','R','R','R','V','V','V','V','V','V'),
				SEQCOV= c(0.16, 0.16, 0.16, 0.16, 0.16, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6),
				ACUTE=	c('low', 'low', 'high', 'low', 'high', 'high', 'high', 'high', 'high', 'high', 'high'),
				GAPS=	c('none', 'low', 'low', 'high', 'high', 'low', 'high', 'none', 'none', 'low', 'high'), 
				ART=	c('none', 'none', 'none', 'none', 'none', 'none', 'none', 'fast', 'fast', 'fast', 'fast'),
				EXT= 	c('5pc', '5pc', '5pc', '5pc', '5pc', '~0pc', '~0pc', '~0pc', '~0pc', '~0pc', '~0pc')
				)
	submitted.info	<- merge(submitted.info, tmp, by='SC')
	#
	#	best tree for each scenario
	#
	submitted.info[, BEST:='N']
	set(submitted.info, submitted.info[, which(grepl('RAxML', FILE) & grepl('best_tree', FILE))], 'BEST', 'Y')	
	#	copied from ListOfBestTrees_IQTree150818.txt
	#	there are several best trees for some scenarios
	tmp	<- c( '150701_Vill_SCENARIO-A_IQTree150814_partition_12_3_07',	
	'150701_Vill_SCENARIO-A_IQTree150814_partition_12_3_04.',	
	'150701_Vill_SCENARIO-B_IQTree150814_partition_12_3_03.',	
	'Vill_99_Apr15_IQTree150814_partition_123.',			
	'150701_Vill_SCENARIO-D_IQTree150814_partition_12_3.',	
	'150701_Vill_SCENARIO-E_IQTree150814_partition_12_3.',
	'150701_Vill_SCENARIO-A_IQTree150814_pol_partition_12_3.',	
	'150701_Vill_SCENARIO-B_IQTree150814_pol_partition_12_3_05.',
	'Vill_99_Apr15_IQTree150814_pol_partition_12_3_09.',		
	'Vill_99_Apr15_IQTree150814_pol_partition_12_3_10.',		
	'150701_Vill_SCENARIO-D_IQTree150814_pol_partition_12_3_05.',
	'150701_Vill_SCENARIO-D_IQTree150814_pol_partition_12_3_06.',
	'150701_Vill_SCENARIO-D_IQTree150814_pol_partition_12_3_09.',
	'150701_Vill_SCENARIO-E_IQTree150814_pol_partition_12_3_06.',
	'150701_Regional_TRAIN1_IQTree150818_partition_123_03.',	
	'150701_Regional_TRAIN2_IQTree150814_partition_123_03.',	
	'150701_Regional_TRAIN3_IQTree150814_partition_123_01.',	
	'150701_Regional_TRAIN4_IQTree150814_partition_123_02.',	
	'150701_Regional_TRAIN5_IQTree150814_partition_123.',
	'150701_Regional_TRAIN1_IQTree150818_pol_partition_123_05.',
	'150701_Regional_TRAIN2_IQTree150814_pol_partition_123_10.',
	'150701_Regional_TRAIN3_IQTree150814_pol_partition_123_05.',
	'150701_Regional_TRAIN3_IQTree150814_pol_partition_123_06.',
	'150701_Regional_TRAIN3_IQTree150814_pol_partition_123_08.',
	'150701_Regional_TRAIN4_IQTree150814_pol_partition_123_10.',
	'150701_Regional_TRAIN5_IQTree150814_pol_partition_123_05.')
	tmp	<- sapply(tmp, function(x) submitted.info[, which((grepl('IQTree150814/', FILE, fixed=1) | grepl('IQTree150818/', FILE, fixed=1)) & grepl(x, FILE, fixed=1))] )
	set(submitted.info, tmp, 'BEST', 'Y')
	#	PhyML no replicates: all files best
	set(submitted.info, submitted.info[, which(TEAM=='PhyML')], 'BEST', 'Y')		
	#
	#	set OTHER (ie old or some preliminary/unknown tree)
	#
	submitted.info[, OTHER:='N']
	#	MetaPIGA tree to be used is first in nexus list (which was tagged with best above)
	set(submitted.info, submitted.info[, which(TEAM=='MetaPIGA' & !grepl('use', FILE))], 'OTHER', 'Y')
	#	IQTree did several uploads, use only most recent in main analysis
	set(submitted.info, submitted.info[, which(grepl('150701_Regional_TRAIN1_IQTree150814', FILE))], 'OTHER', 'Y')
	#
	#	set which gene used to construct tree (either pol or concatenated gag+pol+env)
	#
	submitted.info[, GENE:=NA_character_]
	set(submitted.info, submitted.info[, which(TEAM=='RAXML' & grepl('full', FILE))], 'GENE', 'GAG+POL+ENV')
	set(submitted.info, submitted.info[, which(TEAM=='RAXML' & grepl('pol', FILE))], 'GENE', 'POL')
	stopifnot(nrow(subset(submitted.info, TEAM=='RAXML' & is.na(GENE)))==0)
	set(submitted.info, submitted.info[, which(TEAM=='PhyML')], 'GENE', 'POL')
	set(submitted.info, submitted.info[, which(TEAM=='MetaPIGA')], 'GENE', 'GAG+POL+ENV')	
	set(submitted.info, submitted.info[, which(TEAM=='IQTree' & grepl('[0-9]_partition', FILE))], 'GENE', 'GAG+POL+ENV')
	set(submitted.info, submitted.info[, which(TEAM=='IQTree' & grepl('[0-9]_pol_partition', FILE))], 'GENE', 'POL')
	stopifnot(nrow(subset(submitted.info, TEAM=='IQTree' & is.na(GENE)))==0)
	#
	#	number taxa in tree
	#
	setkey(submitted.info, IDX)
	submitted.info[, TAXAN:= sapply(strs, Ntip)]
	#
	#	are trees rooted?
	#
	setkey(submitted.info, IDX)
	submitted.info[, ROOTED:=factor(sapply(strs, is.rooted),levels=c(TRUE,FALSE),labels=c('Y','N'))]
	#
	#	add index of true tree
	#
	require(phangorn)
	submitted.info	<- merge(submitted.info, subset(tfiles, select=c('SC','IDX_T','TAXAN_T')), by='SC')
	stopifnot(nrow(subset(submitted.info, TAXAN>TAXAN_T))==0)
	#
	#	fix taxa names that teams have changed
	#
	tmp		<- subset(submitted.info, TEAM=='IQTree' & MODEL=='R')[, IDX]
	for(i in tmp)
	{
		strs[[i]]$tip.label	<- sapply(strsplit(strs[[i]]$tip.label,'_'), function(x)	paste(x[1],'_',x[2],'|',x[3],'|',x[4],'_',x[5],'|',x[6],sep='')	)
	}
	for(i in seq_along(strs))
	{
		strs[[i]]$tip.label	<- toupper(strs[[i]]$tip.label)
	}
	for(i in seq_along(ttrs))
	{
		ttrs[[i]]$tip.label	<- toupper(ttrs[[i]]$tip.label)
	}	
	###
	tmp		<- subset(submitted.info, TEAM=='PhyML' & MODEL=='R')[, IDX]
	for(i in tmp)
	{
		
		z	<- data.table(IDX=seq_along(strs[[i]]$tip.label), IDPOP=regmatches(strs[[i]]$tip.label, regexpr('IDPOP_[0-9]+',strs[[i]]$tip.label)), SC=subset(submitted.info, IDX==i)[,SC])
		z	<- merge(subset(tinfo, select=c(IDPOP,SC,TAXA)), z, by=c('IDPOP','SC'))
		setkey(z, IDX)
		strs[[i]]$tip.label	<- z[, TAXA]
	}
	tmp		<- subset(submitted.info, TEAM=='PhyML' & MODEL=='V')[, IDX]
	for(i in tmp)
	{
		
		z	<- data.table(IDX=seq_along(strs[[i]]$tip.label), IDPOP=regmatches(strs[[i]]$tip.label, regexpr('HOUSE[0-9]+-[0-9]+|House[0-9]+-[0-9]+',strs[[i]]$tip.label)), SC=subset(submitted.info, IDX==i)[,SC])
		z	<- merge(subset(tinfo, select=c(IDPOP,SC,TAXA)), z, by=c('IDPOP','SC'))
		stopifnot(nrow(z)==length(strs[[i]]$tip.label))
		setkey(z, IDX)
		strs[[i]]$tip.label	<- z[, TAXA]
	}
	#
	#	compute Robinson Fould of complete tree
	#
	tmp				<- treedist.robinsonfould.wrapper(submitted.info, ttrs, strs)
	submitted.info	<- merge(submitted.info, tmp, by='IDX')
	#
	#	compute path differences on complete trees
	#
	setkey(submitted.info, IDX)
	#tmp				<- subset(submitted.info, IDX==463)[1,]
	#IDX<- 1; IDX_T<- 1
	#IDX<- 822; IDX_T<- 11
	tmp				<- submitted.info[, {
				cat('\nAt IDX', IDX)
				stree		<- unroot(strs[[IDX]])
				otree		<- unroot(multi2di(ttrs[[IDX_T]]))				
				if(!is.binary.tree(stree))
				{
					cat('\nFound non-binary tree at IDX',IDX)
					stree	<- multi2di(stree)
				}
				#print(stree)
				#print(otree)
				z			<- setdiff(otree$tip.label, stree$tip.label)
				stopifnot( length(z)==abs(diff(c(Ntip(otree), Ntip(stree)))) )
				if(length(z))
					otree	<- unroot(drop.tip(otree, z))								
				#normalize with choose(n,2)		
				tmp			<- treedist.pathdifference(otree, stree, lambda=0)
				list(PD=tmp['path'], NPD=tmp['path.std'])
			}, by='IDX']
	submitted.info	<- merge(submitted.info, tmp, by='IDX')
	#	compute Robinson Fould of clusters, then take sum
	tmp			<- treedist.robinsonfouldclusters.wrapper(submitted.info, ttrs, strs, tinfo)
	sclu.info	<- merge(submitted.info, tmp, by='IDX')		
	#
	#
	#
	outfile	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation/submitted_151023.rda'
	save(strs, ttrs, tinfo, submitted.info, sclu.info, file=outfile)
}
##--------------------------------------------------------------------------------------------------------
##	olli 03.12.15
##--------------------------------------------------------------------------------------------------------
treecomparison.submissions.update.160430<- function()
{
	edir			<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation'
	timetag			<- '160430'
	#
	# collect results so far
	#	
	file			<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation/submitted_151203.rda'
	load(file)	
	#
	# to tinfo add actual transmitters
	#
	# check TRAIN1
	load( '/Users/Oliver/Dropbox (Infectious Disease)/PANGEAHIVsim_internal/freeze_July15/150701_Regional_TRAIN1_SIMULATED_INTERNAL.R' )	
	ch				<- subset(tinfo, SC=='150701_REGIONAL_TRAIN1' & BRL_T=='time', TAXA)
	ch[, IDPOP:= as.integer(gsub('IDPOP_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',1)))]
	ch[, GENDER_CH:= sapply(strsplit(TAXA,'|',fixed=1),'[[',2)]
	ch[, DOB_CH:= as.numeric(gsub('DOB_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',3)))]
	ch[, TIME_SEQ_CH:= as.numeric(sapply(strsplit(TAXA,'|',fixed=1),'[[',4))]
	ch				<- merge(subset(df.inds, select=c(IDPOP, GENDER, DOB, TIME_SEQ)), ch, by='IDPOP')	
	stopifnot( ch[, all(DOB==DOB_CH)], ch[, all(GENDER==GENDER_CH)], ch[, all(TIME_SEQ==TIME_SEQ_CH)] )
	subset(ch, TIME_SEQ!=TIME_SEQ_CH)
	# OK :-) schedule adding IDPOP_T
	df.trms[, IDPOP:= df.trms[, paste('IDPOP_',IDREC,sep='')]]
	tinfo.add		<- subset(df.trms, select=c(IDPOP, IDTR))	
	df.trms[, IDPOP:= df.trms[, paste('IDPOP_',IDTR,sep='')]]
	tinfo.add		<- merge(tinfo.add, subset(df.trms, select=c(IDPOP, IDREC)), by='IDPOP', all=1)
	set(ch, NULL, 'IDPOP', ch[, paste('IDPOP_',IDPOP,sep='')])
	tinfo.add		<- merge(tinfo.add, subset(ch, select=IDPOP), by='IDPOP')
	tinfo.add[, SC:='150701_REGIONAL_TRAIN1']
	# check TRAIN2
	load( '/Users/Oliver/Dropbox (Infectious Disease)/PANGEAHIVsim_internal/freeze_July15/150701_Regional_TRAIN2_SIMULATED_INTERNAL.R' )	
	ch				<- subset(tinfo, SC=='150701_REGIONAL_TRAIN2' & BRL_T=='time', TAXA)
	ch[, IDPOP:= as.integer(gsub('IDPOP_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',1)))]
	ch[, GENDER_CH:= sapply(strsplit(TAXA,'|',fixed=1),'[[',2)]
	ch[, DOB_CH:= as.numeric(gsub('DOB_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',3)))]
	ch[, TIME_SEQ_CH:= as.numeric(sapply(strsplit(TAXA,'|',fixed=1),'[[',4))]
	ch				<- merge(subset(df.inds, select=c(IDPOP, GENDER, DOB, TIME_SEQ)), ch, by='IDPOP')	
	stopifnot( ch[, all(DOB==DOB_CH)], ch[, all(GENDER==GENDER_CH)], ch[, all(TIME_SEQ==TIME_SEQ_CH)] )
	subset(ch, TIME_SEQ!=TIME_SEQ_CH)
	# OK :-) schedule adding IDPOP_T
	df.trms[, IDPOP:= df.trms[, paste('IDPOP_',IDREC,sep='')]]
	tmp				<- subset(df.trms, select=c(IDPOP, IDTR))	
	df.trms[, IDPOP:= df.trms[, paste('IDPOP_',IDTR,sep='')]]
	tmp				<- merge(tmp, subset(df.trms, select=c(IDPOP, IDREC)), by='IDPOP', all=1)
	set(ch, NULL, 'IDPOP', ch[, paste('IDPOP_',IDPOP,sep='')])
	tmp		<- merge(tmp, subset(ch, select=IDPOP), by='IDPOP')
	tmp[, SC:='150701_REGIONAL_TRAIN2']
	tinfo.add		<- rbind(tinfo.add, tmp)
	# check TRAIN4
	ch				<- subset(tinfo, SC=='150701_REGIONAL_TRAIN4' & BRL_T=='time', TAXA)
	ch[, IDPOP:= as.integer(gsub('IDPOP_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',1)))]
	ch[, GENDER_CH:= sapply(strsplit(TAXA,'|',fixed=1),'[[',2)]
	ch[, DOB_CH:= as.numeric(gsub('DOB_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',3)))]
	ch[, TIME_SEQ_CH:= as.numeric(sapply(strsplit(TAXA,'|',fixed=1),'[[',4))]
	ch				<- merge(subset(df.inds, select=c(IDPOP, GENDER, DOB, TIME_SEQ)), ch, by='IDPOP')	
	stopifnot( ch[, all(DOB==DOB_CH)], ch[, all(GENDER==GENDER_CH)], ch[, all(TIME_SEQ==TIME_SEQ_CH)] )
	subset(ch, TIME_SEQ!=TIME_SEQ_CH)	
	tmp[, SC:='150701_REGIONAL_TRAIN4']
	tinfo.add		<- rbind(tinfo.add, tmp)
	# check TRAIN3
	load( '/Users/Oliver/Dropbox (Infectious Disease)/PANGEAHIVsim_internal/freeze_July15/150701_Regional_TRAIN3_SIMULATED_INTERNAL.R' )	
	ch				<- subset(tinfo, SC=='150701_REGIONAL_TRAIN3' & BRL_T=='time', TAXA)
	ch[, IDPOP:= as.integer(gsub('IDPOP_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',1)))]
	ch[, GENDER_CH:= sapply(strsplit(TAXA,'|',fixed=1),'[[',2)]
	ch[, DOB_CH:= as.numeric(gsub('DOB_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',3)))]
	ch[, TIME_SEQ_CH:= as.numeric(sapply(strsplit(TAXA,'|',fixed=1),'[[',4))]
	ch				<- merge(subset(df.inds, select=c(IDPOP, GENDER, DOB, TIME_SEQ)), ch, by='IDPOP')	
	stopifnot( ch[, all(DOB==DOB_CH)], ch[, all(GENDER==GENDER_CH)], ch[, all(TIME_SEQ==TIME_SEQ_CH)] )
	subset(ch, TIME_SEQ!=TIME_SEQ_CH)
	# OK :-) schedule adding IDPOP_T
	df.trms[, IDPOP:= df.trms[, paste('IDPOP_',IDREC,sep='')]]
	tmp				<- subset(df.trms, select=c(IDPOP, IDTR))	
	df.trms[, IDPOP:= df.trms[, paste('IDPOP_',IDTR,sep='')]]
	tmp				<- merge(tmp, subset(df.trms, select=c(IDPOP, IDREC)), by='IDPOP', all=1)
	set(ch, NULL, 'IDPOP', ch[, paste('IDPOP_',IDPOP,sep='')])
	tmp		<- merge(tmp, subset(ch, select=IDPOP), by='IDPOP')
	tmp[, SC:='150701_REGIONAL_TRAIN3']
	tinfo.add		<- rbind(tinfo.add, tmp)
	# check TRAIN5
	ch				<- subset(tinfo, SC=='150701_REGIONAL_TRAIN5' & BRL_T=='time', TAXA)
	ch[, IDPOP:= as.integer(gsub('IDPOP_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',1)))]
	ch[, GENDER_CH:= sapply(strsplit(TAXA,'|',fixed=1),'[[',2)]
	ch[, DOB_CH:= as.numeric(gsub('DOB_','',sapply(strsplit(TAXA,'|',fixed=1),'[[',3)))]
	ch[, TIME_SEQ_CH:= as.numeric(sapply(strsplit(TAXA,'|',fixed=1),'[[',4))]
	ch				<- merge(subset(df.inds, select=c(IDPOP, GENDER, DOB, TIME_SEQ)), ch, by='IDPOP')	
	stopifnot( ch[, all(DOB==DOB_CH)], ch[, all(GENDER==GENDER_CH)], ch[, all(TIME_SEQ==TIME_SEQ_CH)] )
	subset(ch, TIME_SEQ!=TIME_SEQ_CH)
	tmp[, SC:='150701_REGIONAL_TRAIN5']
	tinfo.add		<- rbind(tinfo.add, tmp)
	# add transmitters for regional to tinfo
	tinfo			<- merge(tinfo, tinfo.add, by=c('IDPOP', 'SC'), all.x=1)
	#
	# compute closest individual on true trees
	#
	tmp				<- unique(subset(tinfo, select=c(SC, BRL_T, IDX_T)))
	tmp				<- tmp[, {
				print(IDX_T)
				ph			<- ttrs[[IDX_T]]
				model.reg	<- grepl('REGIONAL',SC)
				treedist.closest.ind(ph, model.reg)
			}, by=c('SC','BRL_T','IDX_T')]
	tinfo			<- merge(tinfo, tmp, by=c('SC','BRL_T','IDX_T','IDPOP'))
	set(tinfo, NULL, 'IDPOP_CL', tinfo[, gsub('IDPOP_','',IDPOP_CL)])
	#
	# compute closest individual on simulated trees and determine proportion if either transmitter or among recipients
	#	
	sucl			<- subset(submitted.info, MODEL=='R')[, {
				print(IDX)
				#IDX<- 557; SUB_IDX_T<-2; SC<- '150701_REGIONAL_TRAIN2'
				ph			<- strs[[IDX]]
				model.reg	<- grepl('REGIONAL',SC)
				ans			<- treedist.closest.ind(ph, model.reg)
				ans			<- subset(ans, GD<=0.045)
				tmp			<- subset(tinfo, IDX_T==SUB_IDX_T, c(IDPOP, IDTR, IDREC))
				ans			<- merge(ans, tmp, by='IDPOP')
				set(ans, NULL, 'IDPOP_CL', ans[, gsub('IDPOP_','',IDPOP_CL)])
				ans[, IDCL:= as.character(IDTR)]
				tmp			<- ans[, which(!is.na(IDREC))]
				set(ans, tmp, 'IDCL', ans[tmp, paste(IDCL, IDREC, sep=',')])
				if(nrow(ans))
				{
					tmp			<- ans[, list(CLD= IDPOP_CL%in%strsplit(IDCL,',')[[1]]) , by='IDPOP']
					tmp			<- tmp[, mean(CLD)]
				}
				if(!nrow(ans))
					tmp			<- NA_real_
				list(TR_REC_perc= tmp)
			}, by=c('IDX')]
	submitted.info	<- merge(submitted.info, sucl, by='IDX', all.x=1)
	# compute same proportion on true trees
	tmp				<- subset(tinfo, BRL_T=='subst')[, {
				ans		<- data.table(IDPOP=IDPOP, IDCL= as.character(IDTR), IDREC=IDREC, IDPOP_CL=IDPOP_CL, GD=GD)
				ans		<- subset(ans, GD<=0.045)
				tmp		<- ans[, which(!is.na(IDREC))]
				set(ans, tmp, 'IDCL', ans[tmp, paste(IDCL, IDREC, sep=',')])
				tmp			<- ans[, list(CLD= IDPOP_CL%in%strsplit(IDCL,',')[[1]]) , by='IDPOP']
				list(TR_REC_perc_T= tmp[, mean(CLD)])
			}, by='IDX_T']
	setnames(tmp, 'IDX_T','SUB_IDX_T')
	submitted.info	<- merge(submitted.info, tmp, by='SUB_IDX_T', all.x=1)
	#
	# save
	#
	outfile			<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation/submitted_160430.rda'
	save(strs, strs_lsd_brl, strs_lsd_date, ttrs, tinfo, submitted.info, sclu.info, ttdists1, RFttdists, file=outfile)
}
##--------------------------------------------------------------------------------------------------------
##	olli 03.12.15
##--------------------------------------------------------------------------------------------------------
treecomparison.submissions.update.151203<- function()
{
	edir			<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation'
	timetag			<- '151203'
	#
	# collect results so far
	#
	file			<- "~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation/submitted_151119_SRFQD.rda"
	load(file)
	load("~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation/submitted_151203_CCinfowithTTdists.rda")
	# loads myinfo ttdists1 RFttdists
	submitted.info	<- copy(myinfo)
	setnames(submitted.info, c('RF','NRF','NQD','kc0','kc1','kc_x','kc_y','rf_x','rf_y'), c('SB_RF','SB_NRF','SB_NQD','LSD_KC_L0','LSD_KC_L1','LSD_KC_L0_MDSx','LSD_KC_L0_MDSy','LSD_RF_MDSx','LSD_RF_MDSy'))
	setnames(sclu.info, c('NRFC','NQDC'), c('SB_NRFC','SB_NQDC'))	
	#
	# calculate RF on 'strs_lsd_date'
	#
	#	take topology of true dated trees and compare to topology of LSD dated trees with RF
	tmp				<- treedist.robinsonfould.wrapper(submitted.info, ttrs, strs_lsd_date)
	setnames(tmp, c('RF','NRF'), c('LSD_RF','LSD_NRF'))
	submitted.info	<- merge(submitted.info, tmp, by='IDX')
	#	take topology of true dated trees and compare to topology of LSD dated trees with RF
	tmp				<- treedist.robinsonfouldclusters.wrapper(submitted.info, ttrs, strs_lsd_date, tinfo)
	setnames(tmp, c('RFC','NRFC'), c('LSD_RFC','LSD_NRFC'))	
	sclu.info		<- merge(sclu.info, subset(tmp, select=c(IDX, IDCLU, LSD_NRFC)), by=c('IDX','IDCLU'))
	
	strs.new			<- strs
	ttrs.new			<- ttrs
	tinfo.new			<- copy(tinfo)
	submitted.info.new	<- copy(submitted.info)
	sclu.info.new		<- copy(sclu.info)
	outdir		<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation'
	z			<- load(paste(outdir, 'submitted_151119_S.rda', sep='/'))	
	stopifnot( nrow(subset(merge(subset(submitted.info, select=c('FILE','IDX')), subset(submitted.info.new, select=c('FILE','IDX')), by='FILE'), IDX.x!=IDX.y))==0 )	
	submitted.info		<- merge(submitted.info.new, subset(submitted.info, select=c('IDX','NQD','lm_intercept','lm_slope','lm_rsq')), by='IDX')
	strs				<- strs.new
	ttrs				<- ttrs.new
	sclu.info			<- merge(sclu.info.new, subset(sclu.info, select=c('IDX','IDCLU','NQDC')), by=c('IDX','IDCLU'))
	
	
	#
	# save
	#
	outfile	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation/submitted_151203.rda'
	save(strs, strs_lsd_brl, strs_lsd_date, ttrs, tinfo, submitted.info, sclu.info, ttdists1, RFttdists, file=outfile)	
}
##--------------------------------------------------------------------------------------------------------
##	olli 27.06.11
##--------------------------------------------------------------------------------------------------------
treecomparison.ana.160627.sclu<- function()
{	
	require(ggplot2)
	require(data.table)
	require(ape)
	require(scales)	
	require(ggtree)
	require(phangorn)
	
	edir			<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation'
	timetag			<- '160713'	
	load(file.path(edir,'submitted_160713_07MSELSD.rda'))
		
	set(sclu.info, sclu.info[, which(grepl('gag+pol+env',FILE,fixed=1))], 'GENE', 'GAG+POL+ENV')
	
	sc		<- copy(sclu.info)		
	tmp		<- subset(submitted.info, TEAM=='RUNGAPS_ExaML', c(IDX, RUNGAPS, GENE))
	sc		<- merge(sc, tmp, by=c('IDX','GENE'), all.x=1)	
	sc		<- merge(sc, data.table(GENE=c('P17','GAG','GAG+PARTIALPOL','POL','GAG+POL+ENV'), GENE_L=c(396, 1440, 3080, 2843, 6807)), by='GENE')	
	tmp		<- unique( subset(tinfo, BRL_T=='subst' & grepl('REG',SC), select=c(SC, TAXA, ENV_GAPS_P, FULL_GAPS_P,  GAG_GAPS_P, POL_GAPS_P)) )
	tmp		<- tmp[, list(FULL_GAPS_P=mean(FULL_GAPS_P), GAG_GAPS_P=mean(GAG_GAPS_P), POL_GAPS_P=mean(POL_GAPS_P), ENV_GAPS_P=mean(ENV_GAPS_P)), by='SC']
	tmp		<- melt(tmp, id.vars=c('SC'), variable.name='GENE', value.name='GAPS_P')
	set(tmp, NULL, 'GENE', tmp[, gsub('FULL','GAG+POL+ENV',gsub('_GAPS_P','',GENE))])	
	sc		<- merge(sc, tmp, by=c('SC','GENE'), all.x=1)
	#
	tmp		<- subset(tinfo, !is.na(IDCLU))[, list(CLU_N=CLU_N[1], MXGPS_CLU= max(GPS), MDGPS_CLU=median(GPS)), by=c('SC','IDCLU')]
	sc		<- merge(sc, tmp, by=c('SC','IDCLU'))	
	set(sc, NULL, 'MODEL', sc[, factor(MODEL, levels=c('V','R'),labels=c('Model: Village','Model: Regional'))])
	set(sc, sc[, which(SC=="VILL_99_APR15")],'SC',"150701_VILL_SCENARIO-C")	
	set(sc, NULL, 'SC', sc[, factor(SC,	levels=c("150701_REGIONAL_TRAIN1", "150701_REGIONAL_TRAIN2", "150701_REGIONAL_TRAIN3", "150701_REGIONAL_TRAIN4","150701_REGIONAL_TRAIN5","150701_VILL_SCENARIO-A","150701_VILL_SCENARIO-B","150701_VILL_SCENARIO-C","150701_VILL_SCENARIO-D","150701_VILL_SCENARIO-E"), 
							labels=c('sc 1','sc 2','sc 3','sc 4','sc 5','sc A','sc B','sc C','sc D','sc E'))])
	set(sc, NULL, 'GAPS', sc[, factor(GAPS, levels=c('none','low','high'),labels=c('none','as for\nBotswana\nsequences','as for\nUganda\nsequences'))])
	set(sc, NULL, 'BEST', sc[, factor(BEST, levels=c('Y','N'),labels=c('best tree','replicate tree'))])
	set(sc, sc[, which(GENE=='P17')], 'GENE', 'gag (p17)')
	set(sc, sc[, which(GENE=='GAG')], 'GENE', 'gag')
	set(sc, sc[, which(GENE=='GAG+PARTIALPOL')], 'GENE', 'gag + pol (prot,p51)')		
	set(sc, sc[, which(GENE=='POL')], 'GENE', 'pol')
	set(sc, sc[, which(GENE=='GAG+POL+ENV')], 'GENE', 'gag+pol+env')
	set(sc, sc[, which(TEAM=='IQTree')], 'TEAM', 'IQ-TREE')
	set(sc, sc[, which(TEAM=='RAXML')], 'TEAM', 'RAxML')
	set(sc, NULL, 'EXT', sc[, factor(EXT, levels=c('~0pc','5pc'),labels=c('~ 0%/year','5%/year'))])
	set(sc, NULL, 'ART', sc[, factor(ART, levels=c('none','fast'),labels=c('none','fast'))])
	sc		<- subset(sc, OTHER=='N')
	#
	#	add size adjusted KC
	#	
	if(1)
	{
		require(gamlss)
		kc.std.d	<- subset(sc, TEAM!='MetaPIGA' & SC=='sc 1', select=c(SC, IDX, IDCLU, TEAM,GENE,CLU_N, KC))
		kc.std.m3	<- gamlss(KC~I(CLU_N*(CLU_N-1)/2), data=kc.std.d)
		tmp2		<- subset(sc, SC%in%c('sc 1','sc 2','sc 4'), select=c(SC,IDX, IDCLU, TEAM,GENE,CLU_N, KC))	
		tmp2[, KCadj:= KC / predict(kc.std.m3, data=kc.std.d, newdata=tmp2,type='response', se.fit=FALSE)]	
		ggplot(tmp2, aes(x=CLU_N)) + geom_point(aes(y=KCadj,colour=GENE, pch=TEAM)) +	scale_y_log10() + facet_grid(~SC)
		sc			<- merge(sc, subset(tmp2, select=c(IDX, IDCLU, KCadj)), by=c('IDX','IDCLU'), all.x=1)
		#
		tmp		<- melt(subset(sc, IDX==45), measure.var=c('NPD','NPDSQ','NRFC','NQDC','KCadj'))
		ggplot( tmp, aes(x=CLU_N, y=value, colour=GENE, pch=TEAM)) + geom_point() + facet_grid(GENE+TEAM+IDX~variable)
		#
		#	check dependence on size of cluster
		#
		ggplot( melt(sc, measure.var=c('NPD','NPDSQ','NRFC','NQDC','KCadj')), aes(x=CLU_N, y=value, colour=GENE, pch=TEAM)) + geom_point() + facet_grid(variable+TEAM+IDX~GENE, scales='free_y')
		file	<- file.path(edir, paste(timetag,'_','dependence_on_clustersize.pdf',sep=''))
		ggsave(file=file, w=10, h=1000, limitsize = FALSE, useDingbats=FALSE)		
	}
	#
	#
	#	
	sc		<- sc[, list( 	NRFme=mean(NRFC, na.rm=TRUE), 	
							NQDme=mean(NQDC, na.rm=TRUE), 	
							NPDme=mean(NPD, na.rm=TRUE), 	
							NPDSQme=mean(NPDSQ, na.rm=TRUE),		
							KCAme=mean(KCadj, na.rm=TRUE),
							NRFmd=median(NRFC, na.rm=TRUE), 	
							NQDmd=median(NQDC, na.rm=TRUE), 
							NPDmd=median(NPD, na.rm=TRUE), 	
							NPDSQmd=median(NPDSQ, na.rm=TRUE),	
							KCAmd=median(KCadj, na.rm=TRUE)
							), by=c('SC','GENE','GENE_L','TEAM','BEST','IDX','FILE','GAPS','GAPS_P','RUNGAPS','MODEL','TAXAN','TAXAN_T','ROOTED','SEQCOV','ART','ACUTE','EXT','OTHER')]
	sc		<- subset(sc, MODEL=='Model: Regional')
	#
	#	KC distance standardized with TEAMS separated
	#
	tmp		<- subset(sc, ACUTE=='low' & TEAM%in%c('IQ-TREE', 'PhyML', 'RAxML', 'MetaPIGA'))
	set(tmp, NULL, 'TEAM', tmp[, factor(TEAM)])
	set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag','pol','gag+pol+env'))])
	ggplot(tmp, aes(x=GAPS)) +
			geom_jitter(aes(y=KCAme, colour=GENE), position=position_jitter(w=0.8, h = 0), size=2) +			
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 
			scale_y_continuous(expand=c(0,0), limits=c(0,2.7)) + #, breaks=seq(0,1,0.1), minor_breaks=seq(0,1,0.05)) +
			scale_shape_manual(values=c('IQ-TREE'=15, 'PhyML'=12, 'RAxML'=8, 'MetaPIGA'=17)) +
			labs(	x='\nUnassembled sites in simulated sequences', 
					y='incorrectly estimated topologies of subtrees with 4 taxa\n(standardized Quartett distance)\n',
					colour='part of genome used\nfor tree reconstruction',
					pch='algorithm') +
			theme_bw() + theme(legend.position='bottom') + facet_grid(~TEAM)
	file	<- file.path(edir, paste(timetag,'_','KC_clumean_polvsall_by_gaps_taxan1600_Acute10pc_by_TEAM.pdf',sep=''))
	ggsave(file=file, w=12, h=5, useDingbats=FALSE)
	#
	#	KC distance standardized with TEAMS separated
	#
	tmp		<- subset(sc, ACUTE=='low' & GENE=='gag+pol+env' & TEAM%in%c('IQ-TREE', 'PhyML', 'RAxML', 'MetaPIGA'))
	set(tmp, NULL, 'TEAM', tmp[, factor(TEAM)])
	set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag','pol','gag+pol+env'))])
	ggplot(tmp, aes(x=GAPS)) +
			geom_jitter(aes(y=KCAme, pch=TEAM, colour=GENE), position=position_jitter(w=0.8, h = 0), size=2) +			
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 
			scale_y_continuous(expand=c(0,0), limits=c(0,2)) + #, breaks=seq(0,1,0.1), minor_breaks=seq(0,1,0.05)) +
			scale_shape_manual(values=c('IQ-TREE'=15, 'PhyML'=12, 'RAxML'=8, 'MetaPIGA'=17)) +
			labs(	x='\nunassembled PANGEA-HIV sequences', 
					y='standardized Quartett distance\n',
					colour='part of simulated genome\nused or tree reconstruction',
					pch='tree reconstruction\nmethod') +
			theme_bw() + theme(legend.position='bottom') 
	file	<- file.path(edir, paste(timetag,'_','KC_clumean_polvsall_by_gaps_taxan1600_Acute10pc_by_IQTree.pdf',sep=''))	
	ggsave(file=file, w=4.5, h=6, useDingbats=FALSE)
	#
	#	quartett distance standardized by n choose 4 with TEAMS separated
	#
	tmp		<- subset(sc, ACUTE=='low' & TEAM%in%c('IQ-TREE', 'PhyML', 'RAxML', 'MetaPIGA'))
	set(tmp, NULL, 'TEAM', tmp[, factor(TEAM)])
	set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag','pol','gag+pol+env'))])
	ggplot(tmp, aes(x=GAPS)) +
			geom_jitter(aes(y=NQDme, colour=GENE), position=position_jitter(w=0.8, h = 0), size=2) +			
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 
			scale_y_continuous(labels = scales::percent, expand=c(0,0), limits=c(0, 0.45), breaks=seq(0,1,0.1), minor_breaks=seq(0,1,0.05)) +
			scale_shape_manual(values=c('IQ-TREE'=15, 'PhyML'=12, 'RAxML'=8, 'MetaPIGA'=17)) +
			labs(	x='\nUnassembled sites in simulated sequences', 
					y='incorrectly estimated topologies of subtrees with 4 taxa\n(standardized Quartett distance)\n',
					colour='part of genome used\nfor tree reconstruction',
					pch='algorithm') +
			theme_bw() + theme(legend.position='bottom') + facet_grid(~TEAM)
	file	<- file.path(edir, paste(timetag,'_','QD_clumean_polvsall_by_gaps_taxan1600_Acute10pc_by_TEAM.pdf',sep=''))
	ggsave(file=file, w=12, h=5, useDingbats=FALSE)
	#
	#	quartett distance standardized by n choose 4 only IQ-Tree
	#
	tmp		<- subset(sc, ACUTE=='low' & GENE=='gag+pol+env' & TEAM%in%c('IQ-TREE', 'PhyML', 'RAxML', 'MetaPIGA'))
	tmp[, list(NQDme=mean(NQDme)), by=c('SC','GENE')]
	set(tmp, NULL, 'TEAM', tmp[, factor(TEAM)])
	set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag','pol','gag+pol+env'))])
	ggplot(tmp, aes(x=GAPS)) +
			geom_jitter(aes(y=NQDme, pch=TEAM, colour=GENE), position=position_jitter(w=0.8, h = 0), size=2) +			
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 
			scale_y_continuous(labels = scales::percent, expand=c(0,0), limits=c(0, 0.28), breaks=seq(0,1,0.1), minor_breaks=seq(0,1,0.05)) +
			scale_shape_manual(values=c('IQ-TREE'=15, 'PhyML'=12, 'RAxML'=8, 'MetaPIGA'=17)) +
			labs(	x='\nUnassembled sites in simulated sequences', 
					y='proportion among all subtrees with 4 taxa\n',
					colour='part of simulated genome\nused for tree reconstruction',
					pch='tree reconstruction\nmethod') +
			theme_bw() + theme(legend.position='bottom')
	file	<- file.path(edir, paste(timetag,'_','QD_clumean_polvsall_by_gaps_taxan1600_Acute10pc_IQTree.pdf',sep=''))
	ggsave(file=file, w=4.5, h=6, useDingbats=FALSE)	
	#
	#	quartett distance standardized by n choose 4 MVR and BIONJ
	#
	tmp	<- subset(sc, ACUTE=='low' & TEAM%in%c('RAxML','BioNJ','MVR')) 
	ggplot(tmp, aes(x=GAPS)) +
			geom_jitter(aes(y=NQDme, colour=GENE, pch=TEAM), position=position_jitter(w=0.8, h = 0), size=2) +			
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 
			scale_y_continuous(labels = scales::percent, limits=c(0,1), expand=c(0,0)) +
			#scale_shape_manual(values=c('IQ-TREE'=15, 'PhyML'=12, 'RAxML'=8, 'MetaPIGA'=17)) +
			labs(	x='\nUnassembled sites in simulated sequences', 
					y='incorrectly estimated topologies of subtrees with 4 taxa\n(standardized Quartett distance)\n',
					colour='part of genome used\nfor tree reconstruction',
					pch='algorithm') +
			theme_bw() + theme(legend.position='bottom') 
	file	<- file.path(edir, paste(timetag,'_','QD_clumean_polvsall_by_gaps_taxan1600_Acute10pc_MVRBioNJ.pdf',sep=''))
	ggsave(file=file, w=5, h=7, useDingbats=FALSE)		
	#
	#	quartett distance standardized by n choose 4 
	#
	ggplot(tmp, aes(x=GAPS)) +
			geom_jitter(aes(y=NQDme, colour=GENE, pch=TEAM), position=position_jitter(w=0.8, h = 0), size=2) +			
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 
			scale_y_continuous(labels = scales::percent, expand=c(0,0), limits=c(0, 0.45)) +
			scale_shape_manual(values=c('IQ-TREE'=15, 'PhyML'=12, 'RAxML'=8, 'MetaPIGA'=17)) +
			labs(	x='\nUnassembled sites in simulated sequences', 
					y='incorrectly estimated topologies of subtrees with 4 taxa\n(standardized Quartett distance)\n',
					colour='part of genome used\nfor tree reconstruction',
					pch='algorithm') +
			theme_bw() + theme(legend.position='bottom') 
	file	<- file.path(edir, paste(timetag,'_','QD_clumean_polvsall_by_gaps_taxan1600_Acute10pc.pdf',sep=''))
	ggsave(file=file, w=5, h=7, useDingbats=FALSE)	
	#
	#	quartett distance by % gaps on x-axis
	#	
	ggplot(subset(sc, TEAM=='IQ-TREE' & ACUTE=='low'), aes(x=GAPS_P)) +
			#geom_point(aes(y=NQDme, colour=interaction(GENE,GAPS))) +
			geom_boxplot(aes(y=NQDme, colour=GENE), width=0.1, outlier.shape=NA) +
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) +
			scale_x_continuous(labels = scales::percent) +			
			scale_y_continuous(labels = scales::percent, expand=c(0,0), limits=c(0, 0.35)) +
			labs(	x='\nUnassembled sites in simulated sequences', 
					y='incorrectly estimated topologies of subtrees with 4 taxa\n(standardized Quartett distance)\n',
					colour='part of genome used\nfor tree reconstruction',
					pch='algorithm') +
			theme_bw() + theme(legend.position='bottom') + facet_grid(~GAPS, scales='free_x', space='free_x')
	file	<- file.path(edir, paste(timetag,'_','QD_clumean_polvsall_by_gapspc_taxan1600_Acute10pc.pdf',sep=''))
	ggsave(file=file, w=6, h=5, useDingbats=FALSE)
	#
	#	all tree metrics
	#
	tmp		<- subset(sc, ACUTE=='low' & TEAM%in%c('IQ-TREE', 'PhyML', 'RAxML', 'MetaPIGA'))
	set(tmp, NULL, 'TEAM', tmp[, factor(TEAM)])
	set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag','pol','gag+pol+env'))])
	set(tmp, NULL, 'NRFme', tmp[, 100*NRFme])
	set(tmp, NULL, 'NQDme', tmp[, 100*NQDme])
	tmp		<- melt(tmp, measure.vars=c('NRFme','NQDme','NPDSQme','KCAme'))
	set(tmp, tmp[, which(variable=='NRFme')], 'variable','std. Robinson Fould distance\n(%)')
	set(tmp, tmp[, which(variable=='NQDme')], 'variable','std. Quartett distance\n(%)')
	set(tmp, tmp[, which(variable=='NPDSQme')], 'variable','Path distance\n(upper bound is 2)')
	set(tmp, tmp[, which(variable=='KCAme')], 'variable','Kendall-Colijn distance\n(standardized)')
	
	
	ggplot(tmp, aes(x=GAPS)) +
			geom_jitter(aes(y=value, colour=GENE), position=position_jitter(w=0.8, h = 0), size=2) +			
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 
			scale_y_continuous(limits=c(0,NA)) +
			scale_shape_manual(values=c('IQ-TREE'=15, 'PhyML'=12, 'RAxML'=8, 'MetaPIGA'=17)) +
			labs(	x='\nUnassembled sites in simulated sequences', 
					y='Distance between true and reconstructed tree topologies\n',
					colour='part of genome used\nfor tree reconstruction',
					pch='algorithm') +
			theme_bw() + theme(legend.position='bottom') + facet_grid(variable~TEAM, scales='free_y')
	file	<- file.path(edir, paste(timetag,'_','TOPOOTHER_clumean_polvsall_by_gaps_taxan1600_Acute10pc.pdf',sep=''))
	ggsave(file=file, w=15, h=15, useDingbats=FALSE)
	#
	#	all tree metrics incl MVR BioNJ
	#
	tmp		<- subset(sc, ACUTE=='low' & TEAM%in%c('IQ-TREE', 'PhyML', 'RAxML', 'MetaPIGA', 'MVR', 'BioNJ'))
	set(tmp, NULL, 'TEAM', tmp[, factor(TEAM)])
	set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag','pol','gag+pol+env'))])
	set(tmp, NULL, 'NRFme', tmp[, 100*NRFme])
	set(tmp, NULL, 'NQDme', tmp[, 100*NQDme])
	tmp		<- melt(tmp, measure.vars=c('NRFme','NQDme','NPDSQme'))
	set(tmp, tmp[, which(variable=='NRFme')], 'variable','std. Robinson Fould distance\n(%)')
	set(tmp, tmp[, which(variable=='NQDme')], 'variable','std. Quartett distance\n(%)')
	set(tmp, tmp[, which(variable=='NPDSQme')], 'variable','Path distance\n(upper bound is 2)')
	
	ggplot(tmp, aes(x=GAPS)) +
			geom_jitter(aes(y=value, colour=GENE, pch=TEAM), position=position_jitter(w=0.8, h = 0), size=2) +			
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 
			scale_y_continuous(limits=c(0,NA)) +
			scale_shape_manual(values=c('IQ-TREE'=15, 'PhyML'=12, 'RAxML'=8, 'MetaPIGA'=17, 'BioNJ'=7, 'MVR'=9)) +
			labs(	x='\nUnassembled sites in simulated sequences', 
					y='Distance between true and reconstructed tree topologies\n',
					colour='part of genome used\nfor tree reconstruction',
					pch='algorithm') +
			theme_bw() + theme(legend.position='bottom') + facet_grid(variable~TEAM, scales='free_y')
	file	<- file.path(edir, paste(timetag,'_','TOPOOTHER_clumean_polvsall_by_gaps_taxan1600_Acute10pc_withMVR.pdf',sep=''))
	ggsave(file=file, w=20, h=15, useDingbats=FALSE)	
	#
	#	increasing gap coverage with ExaML - missing sites
	#
	tmp		<- subset(sc, TEAM=='RUNGAPS_ExaML' & !grepl('p51', GENE))
	tmp		<- subset(tmp, GENE=='gag' & RUNGAPS==0.02 | GENE=='gag (p17)' & RUNGAPS==0.02 | GENE=='gag+pol+env')
	tmp[, MISSING_P:= (RUNGAPS*GENE_L + (6807-GENE_L))/6807]
	set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag (p17)','gag','gag+pol+env'))])	
	ggplot(tmp, aes(x=MISSING_P)) +
			geom_point(aes(y=NQDme, colour=GENE), size=2, pch=16) +
			scale_colour_manual(values=c('gag (p17)'="#8C510A", 'gag'='red','gag+pol+env'="#3F4788FF")) +
			scale_shape_manual(values=c('Botswana'=23, 'Uganda'=24)) +
			scale_x_continuous(labels = scales::percent, expand=c(0,0), limits=c(0, 1), breaks=seq(0,1,0.1)) +
			scale_y_continuous(labels = scales::percent, expand=c(0,0), limits=c(0, 0.3), breaks=seq(0,1,0.05)) +
			#scale_shape_manual(values=c('IQ-TREE'=15, 'PhyML'=12, 'RAxML'=8, 'MetaPIGA'=17)) +
			labs(	x='\nProportion of missing sites, relative to gag+pol+env genome', 
					y='proportion among all subtrees with 4 taxa\n',
					colour='part of simulated genome\nused for tree reconstruction') +
			theme_bw() + theme(legend.position='bottom') 
	file	<- file.path(edir, paste(timetag,'_','QD_clumean_p17full_by_missingsites_taxan1600_Acute10pc.pdf',sep=''))
	ggsave(file=file, w=5, h=7)
	#
	#	increasing gap coverage with ExaML
	#	
	tmp		<- subset(sc, TEAM=='RUNGAPS_ExaML' & !grepl('p51', GENE))
	set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag (p17)','gag','gag+pol+env'))])
	setkey(tmp, GENE, RUNGAPS)	
	tmp2	<- tmp[, 	list(	RUNGAPS= RUNGAPS, 
						NQDmeSM= predict(loess(NQDme~RUNGAPS, span=5))), 
						by='GENE']
	tmp		<- merge(tmp, tmp2, by=c('GENE','RUNGAPS'))		
	tmp2	<- merge( rbind(data.table(GENE=c('gag (p17)','gag','gag+pol+env'), RUNGAPS=c(0.11, 0.08, 0.17), LOC='Botswana'), data.table(GENE=c('gag (p17)','gag','gag+pol+env'), RUNGAPS=c(0.21, 0.18, 0.47), LOC='Uganda')), tmp2,by=c('GENE','RUNGAPS')) 
	#set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag (p17)','gag','gag + pol (prot,p51)','gag+pol+env'))])
	ggplot(tmp, aes(x=RUNGAPS)) +
			geom_point(aes(y=NQDme, colour=GENE), size=2, pch=8) +
			geom_line(aes(y=NQDmeSM, colour=GENE), size=0.5) +
			geom_point(data=tmp2, aes(y=NQDmeSM, pch=LOC), size=2.5, fill='black') +			
			scale_colour_manual(values=c('gag (p17)'="#8C510A", 'gag'='red','gag+pol+env'="#3F4788FF")) +
			scale_shape_manual(values=c('Botswana'=23, 'Uganda'=24)) +
			scale_x_continuous(labels = scales::percent, expand=c(0,0), limits=c(0, 0.61), breaks=seq(0,1,0.1)) +
			scale_y_continuous(labels = scales::percent, expand=c(0,0), limits=c(0, 0.4)) +
			#scale_shape_manual(values=c('IQ-TREE'=15, 'PhyML'=12, 'RAxML'=8, 'MetaPIGA'=17)) +
			labs(	x='\nUnassembled sites in simulated sequences', 
					y='incorrectly estimated topologies of subtrees with 4 taxa\n(standardized Quartett distance)\n',
					colour='part of genome used\nfor tree reconstruction',
					pch='sampling location'
					) +
			theme_bw() + theme(legend.position='bottom') 
	file	<- file.path(edir, paste(timetag,'_','QD_clumean_p17full_by_rungaps_taxan1600_Acute10pc.pdf',sep=''))
	ggsave(file=file, w=5, h=7)		
}
##--------------------------------------------------------------------------------------------------------
##	olli 27.06.11
##--------------------------------------------------------------------------------------------------------
treecomparison.ana.160627.strs<- function()
{	
	require(ggplot2)
	require(data.table)
	require(ape)
	require(scales)	
	require(ggtree)
	require(phangorn)
	
	edir			<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation'
	#timetag		<- '160627'
	timetag			<- '160713'		
	load(file.path(edir,'submitted_160713_09SBRL.rda'))
	set(submitted.info, submitted.info[, which(grepl('gag+pol+env',FILE,fixed=1))], 'GENE', 'GAG+POL+ENV')
	
	sa		<- copy(submitted.info)	
	sa		<- merge(sa, data.table(GENE=c('P17','GAG','GAG+PARTIALPOL','POL','GAG+POL+ENV'), GENE_L=c(396, 1440, 3080, 2843, 6807)), by='GENE')
	set(sa, NULL, 'MODEL', sa[, factor(MODEL, levels=c('V','R'),labels=c('Model: Village','Model: Regional'))])
	set(sa, sa[, which(SC=="VILL_99_APR15")],'SC',"150701_VILL_SCENARIO-C")	
	set(sa, NULL, 'SC', sa[, factor(SC,	levels=c("150701_REGIONAL_TRAIN1", "150701_REGIONAL_TRAIN2", "150701_REGIONAL_TRAIN3", "150701_REGIONAL_TRAIN4","150701_REGIONAL_TRAIN5","150701_VILL_SCENARIO-A","150701_VILL_SCENARIO-B","150701_VILL_SCENARIO-C","150701_VILL_SCENARIO-D","150701_VILL_SCENARIO-E"), 
										labels=c('sc 1','sc 2','sc 3','sc 4','sc 5','sc A','sc B','sc C','sc D','sc E'))])
	set(sa, NULL, 'GAPS', sa[, factor(GAPS, levels=c('none','low','high'),labels=c('none','as for\nBotswana\nsequences','as for\nUganda\nsequences'))])	
	set(sa, NULL, 'BEST', sa[, factor(BEST, levels=c('Y','N'),labels=c('best tree','replicate tree'))])
	set(sa, sa[, which(GENE=='P17')], 'GENE', 'gag (p17)')
	set(sa, sa[, which(GENE=='GAG')], 'GENE', 'gag')
	set(sa, sa[, which(GENE=='GAG+PARTIALPOL')], 'GENE', 'gag + pol (prot,p51)')		
	set(sa, sa[, which(GENE=='POL')], 'GENE', 'pol')
	set(sa, sa[, which(GENE=='GAG+POL+ENV')], 'GENE', 'gag+pol+env')
	set(sa, sa[, which(TEAM=='IQTree')], 'TEAM', 'IQ-TREE')
	set(sa, sa[, which(TEAM=='RAXML')], 'TEAM', 'RAxML')
	set(sa, NULL, 'EXT', sa[, factor(EXT, levels=c('~0pc','5pc'),labels=c('~ 0%/year','5%/year'))])
	set(sa, NULL, 'ACUTE', sa[, factor(ACUTE, levels=c('low','high'),labels=c('10%','40%'))])
	set(sa, NULL, 'ART', sa[, factor(ART, levels=c('none','fast'),labels=c('none','fast'))])
	sa		<- subset(sa, OTHER=='N')
	#
	#	on full tree
	#	
	
	#
	#	sum of branch lengths on full genome by method
	#
	tmp		<- subset(sa, ACUTE=='10%' & GENE=='gag+pol+env' & TEAM%in%c('PhyML','MetaPIGA','IQ-TREE', 'RAxML'), select=c(IDX, GENE, TEAM, GAPS, SUM_BRANCHES_T, SUM_BRANCHES))
	set(tmp, NULL, 'TEAM', tmp[, factor(TEAM)])
	set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag','pol','gag+pol+env'))])	
	ggplot(tmp, aes(x=GAPS)) +
			geom_jitter(aes(y=sign(SUM_BRANCHES_T-SUM_BRANCHES)*log10(abs(SUM_BRANCHES_T-SUM_BRANCHES)), colour=GENE, pch=TEAM), position=position_jitter(w=0.8, h = 0), size=2) +
			scale_shape_manual(values=c('IQ-TREE'=15, 'PhyML'=12, 'RAxML'=8, 'MetaPIGA'=17)) +			
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) +			
			scale_y_continuous(labels=c(-100,-10,0,10,100), limits=c(-log10(300),log10(300)), expand=c(0,0), breaks=seq(-2,2,1)) +
			labs(	x='\nunassembled sites of PANGEA-HIV sequences', 
					y='sum of branches in true tree -\nsum of branches in reconstructed tree',					
					colour='part of simulated genome\nused for tree reconstruction',
					pch='tree reconstruction\nmethod') +
			theme_bw() + theme(legend.position='bottom')
	file	<- file.path(edir, paste(timetag,'_','sumBranches_by_gaps.pdf',sep=''))
	ggsave(file=file, w=4.5, h=6, useDingbats=FALSE)
	#
	#	true transmission pairs among closest with distance < 1% 
	#	
	tmp		<- subset(sa, ACUTE=='10%' & GENE=='gag+pol+env' & TEAM%in%c('PhyML','MetaPIGA','IQ-TREE', 'RAxML'), select=c(IDX, GENE, TEAM, GAPS, TPAIR_PHCL_1, NTPAIR_PHCL_1))
	set(tmp, NULL, 'TEAM', tmp[, factor(TEAM)])
	set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag','pol','gag+pol+env'))])	
	ggplot(tmp, aes(x=GAPS)) +
			geom_jitter(aes(y=NTPAIR_PHCL_1/(TPAIR_PHCL_1+NTPAIR_PHCL_1), colour=GENE, pch=TEAM), position=position_jitter(w=0.8, h = 0), size=2) +
			scale_shape_manual(values=c('IQ-TREE'=15, 'PhyML'=12, 'RAxML'=8, 'MetaPIGA'=17)) +			
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 
			scale_y_continuous(labels = scales::percent, limits=c(0,1), expand=c(0,0)) +
			labs(	x='\nunassembled sites of PANGEA-HIV sequences', 
					y='no transmission pair\n',					
					colour='part of simulated genome\nused for tree reconstruction',
					pch='tree reconstruction\nmethod') +
			theme_bw() + theme(legend.position='bottom')
	file	<- file.path(edir, paste(timetag,'_','pTransPairAmong1PCDist_by_gaps.pdf',sep=''))
	ggsave(file=file, w=4.5, h=6, useDingbats=FALSE)
	#
	#	long
	#
	tmp		<- subset(sa, ACUTE=='10%' & GENE=='gag+pol+env' & TEAM%in%c('PhyML', 'MetaPIGA','IQ-TREE', 'RAxML'), select=c(IDX, GENE, TEAM, GAPS, TPAIR_PHCL_1, NTPAIR_PHCL_1))
	set(tmp, NULL, 'TEAM', tmp[, factor(TEAM)])
	set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag','pol','gag+pol+env'))])
	tmp		<- melt(tmp, measure.vars=c('TPAIR_PHCL_1','NTPAIR_PHCL_1'))
	set(tmp, tmp[, which(variable=='TPAIR_PHCL_1')], 'variable', 'Yes')
	set(tmp, tmp[, which(variable=='NTPAIR_PHCL_1')], 'variable', 'No')
	ggplot(tmp, aes(x=factor(IDX), fill=variable)) +
			geom_bar(aes(y=value), stat='identity', position='stack') +
			facet_wrap(~TEAM+GAPS, ncol=3, scales='free') +
			scale_fill_manual(values=c('Yes'='red','No'="grey60")) +
			labs(	x='\nreplicate tree reconstructions', 
					y='phylogenetic pairs < 1% subst/site\n',					
					fill='true transmission pair\nin simulation') +
			theme_bw() + theme(legend.position='bottom', axis.text.x=element_blank(), axis.ticks.x=element_blank())
	file	<- file.path(edir, paste(timetag,'_','pTransPairAmong1PCDist_by_gaps_long.pdf',sep=''))
	ggsave(file=file, w=8, h=12, useDingbats=FALSE)
	#
	#	by missing sites
	#
	tmp		<- subset(sa, TEAM=='RUNGAPS_ExaML' & !grepl('p51', GENE))
	tmp		<- subset(tmp, GENE=='gag' & RUNGAPS==0.02 | GENE=='gag (p17)' & RUNGAPS==0.02 | GENE=='gag+pol+env')
	tmp[, MISSING_P:= (RUNGAPS*GENE_L + (6807-GENE_L))/6807]
	set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag (p17)','gag','gag+pol+env'))])
	setkey(tmp, GENE, RUNGAPS)			
	ggplot(tmp, aes(x=MISSING_P)) +
			geom_point(aes(y=NTPAIR_PHCL_1/(TPAIR_PHCL_1+NTPAIR_PHCL_1), colour=GENE), size=2, pch=16) +
			#geom_line(aes(y=1-TR_PAIR_recSM, colour=GENE), size=0.5) +
			#geom_point(data=tmp2, aes(y=1-TR_PAIR_recSM, pch=LOC), size=2.5, fill='black') +			
			scale_colour_manual(values=c('gag (p17)'="#8C510A", 'gag + pol (prot,p51)'='green','gag'='red','gag+pol+env'="#3F4788FF")) +
			scale_shape_manual(values=c('Botswana'=23, 'Uganda'=24)) +
			scale_x_continuous(labels = scales::percent, expand=c(0,0), breaks=seq(0,1,0.1), limits=c(0,1)) +
			scale_y_continuous(labels = scales::percent, expand=c(0,0), limits=c(0,1)) +
			#scale_shape_manual(values=c('IQ-TREE'=15, 'PhyML'=12, 'RAxML'=8, 'MetaPIGA'=17)) +
			labs(	x='\nproportion of missing sites, relative to gag+pol+env genome', 
					y='proportion of sampled transmission pairs\n',
					colour='part of simulated genome\nused for tree reconstruction',
					pch='sampling location') +
			theme_bw() + theme(legend.position='bottom') 
	file	<- file.path(edir, paste(timetag,'_','pTransPairAmong1PCDist_by_missingsites.pdf',sep=''))
	ggsave(file=file, w=4.5, h=6, useDingbats=FALSE)
	#
	#	by rungaps
	#
	tmp		<- subset(sa, TEAM=='RUNGAPS_ExaML' & !grepl('p51', GENE))
	set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag (p17)','gag','gag+pol+env'))])
	setkey(tmp, GENE, RUNGAPS)			
	tmp[, Y:=NTPAIR_PHCL_1/(TPAIR_PHCL_1+NTPAIR_PHCL_1)]
	
	tmp2	<- tmp[, 	list(	RUNGAPS= RUNGAPS, 
								YM= predict(loess(Y~RUNGAPS, span=5, degree=1))), 
			by='GENE']
	tmp		<- merge(tmp, tmp2, by=c('GENE','RUNGAPS'))		
	tmp2	<- merge( rbind(data.table(GENE=c('gag (p17)','gag','gag + pol (prot,p51)','gag+pol+env'), RUNGAPS=c(0.11, 0.08, 0.14, 0.17), LOC='Botswana'), data.table(GENE=c('gag (p17)','gag','gag + pol (prot,p51)','gag+pol+env'), RUNGAPS=c(0.21, 0.18, 0.34, 0.47), LOC='Uganda')), tmp2,by=c('GENE','RUNGAPS')) 
	#set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag (p17)','gag','gag + pol (prot,p51)','gag+pol+env'))])
	ggplot(tmp, aes(x=RUNGAPS)) +
			geom_point(aes(y=Y, colour=GENE), size=2, pch=8) +
			geom_line(aes(y=YM, colour=GENE), size=0.5) +
			geom_point(data=tmp2, aes(y=YM, pch=LOC), size=2.5, fill='black') +			
			scale_colour_manual(values=c('gag (p17)'="#8C510A", 'gag + pol (prot,p51)'='green','gag'='red','gag+pol+env'="#3F4788FF")) +
			scale_shape_manual(values=c('Botswana'=23, 'Uganda'=24)) +
			scale_x_continuous(labels = scales::percent, expand=c(0,0), breaks=seq(0,1,0.1)) +
			scale_y_continuous(labels = scales::percent, expand=c(0,0), limits=c(0,0.8)) +
			#scale_shape_manual(values=c('IQ-TREE'=15, 'PhyML'=12, 'RAxML'=8, 'MetaPIGA'=17)) +
			labs(	x='\nUnassembled sites in simulated sequences', 
					y='no transmission pair\n',
					colour='part of genome used\nfor tree reconstruction',
					pch='sampling location') +
			theme_bw() + theme(legend.position='bottom') 
	file	<- file.path(edir, paste(timetag,'_','pTransPairAmong1PCDist_by_rungaps.pdf',sep=''))
	ggsave(file=file, w=5, h=7, useDingbats=FALSE)

	
	
	#
	#	proportion of recovered transmission pairs 
	#	
	tmp		<- subset(sa, !is.na(TR_PAIR_rec) & ACUTE=='10%' & GENE=='gag+pol+env' & TEAM%in%c('IQ-TREE', 'PhyML', 'RAxML', 'MetaPIGA'))
	set(tmp, NULL, 'TEAM', tmp[, factor(TEAM)])
	set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag','pol','gag+pol+env'))])	
	ggplot(tmp, aes(x=GAPS)) +
			geom_jitter(aes(y=1-TR_PAIR_rec, colour=GENE, pch=TEAM), position=position_jitter(w=0.8, h = 0), size=2) +
			scale_shape_manual(values=c('IQ-TREE'=15, 'PhyML'=12, 'RAxML'=8, 'MetaPIGA'=17)) +			
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 
			scale_y_continuous(labels = scales::percent, limit=c(0,0.25)) +
			labs(	x='\nunassembled sites of PANGEA-HIV sequences', 
					y='proportion of sampled transmission pairs\n',					
					colour='part of simulated genome used\nfor tree reconstruction',
					pch='tree reconstruction method') +
			theme_bw() + theme(legend.position='bottom')
	file	<- file.path(edir, paste(timetag,'_','pTransPairRecovered_by_gaps.pdf',sep=''))
	ggsave(file=file, w=5, h=7, useDingbats=FALSE)
	#		
	#	proportion of recovered transmission pairs RAxML
	#	
	tmp		<- subset(sa, !is.na(TR_PAIR_rec) & ACUTE=='10%' & TEAM=='RAxML')
	tmp[, list(TR_PAIR_nrec_pc= mean(1-TR_PAIR_rec) ), by=c('SC','GENE')]
	#		
	#	proportion of recovered transmission pairs IQTree
	#	
	tmp		<- subset(sa, !is.na(TR_PAIR_rec) & ACUTE=='10%' & TEAM=='IQ-TREE')
	set(tmp, NULL, 'TEAM', tmp[, factor(TEAM)])
	set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag','pol','gag+pol+env'))])	
	ggplot(tmp, aes(x=GAPS)) +
			geom_jitter(aes(y=1-TR_PAIR_rec, colour=GENE), position=position_jitter(w=0.8, h = 0), size=2) +			
			scale_shape_manual(values=c('IQ-TREE'=15, 'PhyML'=12, 'RAxML'=8, 'MetaPIGA'=17)) +			
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 
			scale_y_continuous(labels = scales::percent, limit=c(-0.01,0.47), expand=c(0,0), minor_breaks=seq(0,1,0.05), breaks=seq(0,1,0.1)) +
			labs(	x='\nUnassembled sites in full-genome sequences', 
					y='phylogenetically closest pairs of individuals\nthat are transmission pairs, out of all such pairs\nthat can be identified in the true tree\n',					
					colour='part of genome used\nfor tree reconstruction',
					pch='algorithm') +
			theme_bw() + theme(legend.position='bottom') + facet_grid(~TEAM)
	file	<- file.path(edir, paste(timetag,'_','pTransPairRecovered_by_gaps_IQTree.pdf',sep=''))
	ggsave(file=file, w=4.5, h=6, useDingbats=FALSE)	
	#		
	#	proportion of recovered transmission pairs by TEAM
	#	
	ggplot(tmp, aes(x=GAPS)) +
			geom_jitter(aes(y=1-TR_PAIR_rec, colour=GENE), position=position_jitter(w=0.8, h = 0), size=2) +			
			scale_shape_manual(values=c('IQ-TREE'=15, 'PhyML'=12, 'RAxML'=8, 'MetaPIGA'=17)) +			
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 
			scale_y_continuous(labels = scales::percent, limit=c(-0.01,0.47), expand=c(0,0), minor_breaks=seq(0,1,0.05), breaks=seq(0,1,0.1)) +
			labs(	x='\nUnassembled sites in full-genome sequences', 
					y='phylogenetically closest pairs of individuals\nthat are transmission pairs, out of all such pairs\nthat can be identified in the true tree\n',					
					colour='part of genome used\nfor tree reconstruction',
					pch='algorithm') +
			theme_bw() + theme(legend.position='bottom') + facet_grid(~TEAM)
	file	<- file.path(edir, paste(timetag,'_','pTransPairRecovered_by_gaps_by_TEAM.pdf',sep=''))
	ggsave(file=file, w=12, h=5, useDingbats=FALSE)
	#		
	#	proportion of recovered transmission pairs by TEAM
	#	
	tmp		<- subset(sa, !is.na(TR_PAIR_rec) & ACUTE=='10%' & TEAM%in%c('IQ-TREE', 'PhyML', 'RAxML', 'MetaPIGA','BioNJ','MVR'))
	set(tmp, NULL, 'TEAM', tmp[, factor(TEAM)])
	set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag','pol','gag+pol+env'))])	
	ggplot(tmp, aes(x=GAPS)) +
			geom_jitter(aes(y=TR_PAIR_rec, colour=GENE, pch=TEAM), position=position_jitter(w=0.8, h = 0), size=2) +
			scale_shape_manual(values=c('IQ-TREE'=15, 'PhyML'=12, 'RAxML'=8, 'MetaPIGA'=17,'MVR'=7,'BioNJ'=9)) +			
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 
			scale_y_continuous(labels = scales::percent, limit=c(0,1)) +
			labs(	x='\nUnassembled sites in full-genome sequences', 
					y='phylogenetically closest pairs of individuals\nthat are transmission pairs, out of all such pairs\nthat can be identified in the true tree\n',					
					colour='part of genome used\nfor tree reconstruction',
					pch='algorithm') +
			theme_bw() + theme(legend.position='bottom') + facet_grid(~TEAM)
	file	<- file.path(edir, paste(timetag,'_','pTransPairRecovered_by_gaps_by_TEAM_withBioNJMVR.pdf',sep=''))
	ggsave(file=file, w=15, h=6, useDingbats=FALSE)	
	#
	#	proportion of recovered transmission pairs by rungaps -- overall missing sites
	#	
	#	confirm proportion unassembled
	if(0)
	{
		dre		<- data.table(FA=list.files('~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/running_gaps_simulations', pattern='fa$', full.names=TRUE))	
		dre[, GENE:= sapply(strsplit(basename(FA),'_'),'[[',4)]
		dre[, RUNGAPS:= sapply(strsplit(basename(FA),'_'),'[[',3)]
		set(dre, NULL, 'RUNGAPS', dre[, as.numeric(gsub('TRAIN2','',RUNGAPS))/100])
		set(dre, NULL, 'GENE', dre[, as.character(factor(GENE, levels=c('P17','GAG','GAGPP', 'FULL'), labels=c('gag (p17)', 'gag', 'gag + pol (prot,p51)','gag+pol+env')))])
		dre		<- dre[, {
					sq	<- read.dna(FA, format='fa')
					z	<- apply(as.character(sq),1,function(x) length(which(x=='?'))) / ncol(sq)
					list(RUNGAPS_E=mean(z))
				}, by=c('FA','GENE','RUNGAPS')]
		dre[, FA:=NULL]	
		tmp		<- subset(sa, TEAM=='RUNGAPS_ExaML' & !grepl('p51', GENE))
		tmp		<- merge(dre, tmp, by=c('GENE','RUNGAPS'))
		tmp[, MISSING_P_E:= (RUNGAPS_E*GENE_L + (6807-GENE_L))/6807]
		tmp[, MISSING_P:= (RUNGAPS*GENE_L + (6807-GENE_L))/6807]
	}
	tmp		<- subset(sa, TEAM=='RUNGAPS_ExaML' & !grepl('p51', GENE))
	tmp		<- subset(tmp, GENE=='gag' & RUNGAPS==0.02 | GENE=='gag (p17)' & RUNGAPS==0.02 | GENE=='gag+pol+env')
	tmp[, MISSING_P:= (RUNGAPS*GENE_L + (6807-GENE_L))/6807]
	set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag (p17)','gag','gag+pol+env'))])
	setkey(tmp, GENE, RUNGAPS)			
	ggplot(tmp, aes(x=MISSING_P)) +
			geom_point(aes(y=1-TR_PAIR_rec, colour=GENE), size=2, pch=16) +
			#geom_line(aes(y=1-TR_PAIR_recSM, colour=GENE), size=0.5) +
			#geom_point(data=tmp2, aes(y=1-TR_PAIR_recSM, pch=LOC), size=2.5, fill='black') +			
			scale_colour_manual(values=c('gag (p17)'="#8C510A", 'gag + pol (prot,p51)'='green','gag'='red','gag+pol+env'="#3F4788FF")) +
			scale_shape_manual(values=c('Botswana'=23, 'Uganda'=24)) +
			scale_x_continuous(labels = scales::percent, expand=c(0,0), breaks=seq(0,1,0.1), limits=c(0,1)) +
			scale_y_continuous(labels = scales::percent, expand=c(0,0), limits=c(0,0.15)) +
			#scale_shape_manual(values=c('IQ-TREE'=15, 'PhyML'=12, 'RAxML'=8, 'MetaPIGA'=17)) +
			labs(	x='\nproportion of missing sites, relative to gag+pol+env genome', 
					y='proportion of sampled transmission pairs\n',
					colour='part of simulated genome\nused for tree reconstruction',
					pch='sampling location') +
			theme_bw() + theme(legend.position='bottom') 
	file	<- file.path(edir, paste(timetag,'_','pTransPairRecovered_p17full_by_missingsites_taxan1600_Acute10pc.pdf',sep=''))
	ggsave(file=file, w=5, h=7, useDingbats=FALSE)
	#
	#	proportion of recovered transmission pairs by rungaps
	#
	tmp		<- subset(sa, TEAM=='RUNGAPS_ExaML' & !grepl('p51', GENE))
	set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag (p17)','gag','gag+pol+env'))])
	setkey(tmp, GENE, RUNGAPS)			
	tmp2	<- tmp[, 	list(	RUNGAPS= RUNGAPS, 
								TR_PAIR_recSM= predict(loess(TR_PAIR_rec~RUNGAPS, span=5))), 
						by='GENE']
	tmp		<- merge(tmp, tmp2, by=c('GENE','RUNGAPS'))		
	tmp2	<- merge( rbind(data.table(GENE=c('gag (p17)','gag','gag + pol (prot,p51)','gag+pol+env'), RUNGAPS=c(0.11, 0.08, 0.14, 0.17), LOC='Botswana'), data.table(GENE=c('gag (p17)','gag','gag + pol (prot,p51)','gag+pol+env'), RUNGAPS=c(0.21, 0.18, 0.34, 0.47), LOC='Uganda')), tmp2,by=c('GENE','RUNGAPS')) 
	#set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag (p17)','gag','gag + pol (prot,p51)','gag+pol+env'))])
	ggplot(tmp, aes(x=RUNGAPS)) +
			geom_point(aes(y=1-TR_PAIR_rec, colour=GENE), size=2, pch=8) +
			geom_line(aes(y=1-TR_PAIR_recSM, colour=GENE), size=0.5) +
			geom_point(data=tmp2, aes(y=1-TR_PAIR_recSM, pch=LOC), size=2.5, fill='black') +			
			scale_colour_manual(values=c('gag (p17)'="#8C510A", 'gag + pol (prot,p51)'='green','gag'='red','gag+pol+env'="#3F4788FF")) +
			scale_shape_manual(values=c('Botswana'=23, 'Uganda'=24)) +
			scale_x_continuous(labels = scales::percent, expand=c(0,0), breaks=seq(0,1,0.1)) +
			scale_y_continuous(labels = scales::percent, expand=c(0,0)) +
			#scale_shape_manual(values=c('IQ-TREE'=15, 'PhyML'=12, 'RAxML'=8, 'MetaPIGA'=17)) +
			labs(	x='\nUnassembled sites in simulated sequences', 
					y='phylogenetically closest pairs of individuals\nthat are transmission pairs, out of all such pairs\nthat can be identified in the true tree\n',
					colour='part of genome used\nfor tree reconstruction',
					pch='sampling location') +
			theme_bw() + theme(legend.position='bottom') 
	file	<- file.path(edir, paste(timetag,'_','pTransPairRecovered_p17full_by_rungaps_taxan1600_Acute10pc.pdf',sep=''))
	ggsave(file=file, w=5, h=7, useDingbats=FALSE)
	#
	#	long branches on regional
	#
	lba.su		<- merge(lba, subset(sa, select=c(IDX, TAXAN, RUNGAPS, OTHER)), by='IDX')
	set(lba.su, NULL, 'GAPS', lba.su[, factor(GAPS, levels=c('none','low','high'),labels=c('none','as for\nBotswana\nsequences','as for\nUganda\nsequences'))])		
	set(lba.su, lba.su[, which(GENE=='P17')], 'GENE', 'gag (p17)')
	set(lba.su, lba.su[, which(GENE=='GAG')], 'GENE', 'gag')
	set(lba.su, lba.su[, which(GENE=='GAG+PARTIALPOL')], 'GENE', 'gag + pol (prot,p51)')		
	set(lba.su, lba.su[, which(GENE=='POL')], 'GENE', 'pol')
	set(lba.su, lba.su[, which(GENE=='GAG+POL+ENV')], 'GENE', 'gag+pol+env')
	set(lba.su, lba.su[, which(TEAM=='IQTree')], 'TEAM', 'IQ-TREE')
	set(lba.su, lba.su[, which(TEAM=='RAXML')], 'TEAM', 'RAxML')	
	#	get reference of error without long branches
	lba.su[, ERR:= DEPTH_T-DEPTH]
	#	these are the closest trees in terms of NRF
	#ref.box		<- subset(lba.su, IDX==858)[, quantile(ERR, p=c(0.25, 0.75))+c(-1,1)*3*diff(quantile(ERR, p=c(0.25, 0.75)))]
	#ref.box		<- subset(lba.su, IDX==2)[, quantile(ERR, p=c(0.25, 0.75))+c(-1,1)*3*diff(quantile(ERR, p=c(0.25, 0.75)))]
	#subset(lba.su, IDX==858)[, sd(ERR)*10]
	#	this suggests the following Tukey criterion:
	#ref.box		<- c(-0.04,0.04)
	ref.box		<- c(-0.1,0.1)
	#
	#	severe branch lengths errors by team
	#
	tmp			<- subset(lba.su, OTHER=='N' & TEAM!='RUNGAPS_ExaML')[, list(TAXAN=length(ERR), OUTLIER_P=mean(ERR<ref.box[1] | ERR>ref.box[2])), by=c('MODEL','SC','TEAM','GAPS','GENE','IDX')]
	tmp[, OUTLIER_N:=TAXAN*OUTLIER_P]
	set(tmp, NULL, 'TEAM', tmp[, factor(TEAM)])
	set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag','pol','gag+pol+env'))])		
	ggplot(tmp, aes(x=GAPS, y=OUTLIER_P, colour=GENE)) + 
			geom_jitter(position=position_jitter(w=0.8, h = 0), size=1) +
			scale_shape_manual(values=c('IQ-TREE'=15, 'PhyML'=12, 'RAxML'=8, 'MetaPIGA'=17)) +			
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 			
			scale_y_continuous(labels = scales::percent) +
			facet_grid(TEAM~GENE, scales='free_y')  + theme_bw() + theme(legend.position='bottom') +
			labs(	x='\nUnassembled sites in simulated sequences', 
					y='branch length to root\n50% too small or too large\n(% of all taxa in tree)\n')
	file	<- file.path(edir, paste(timetag,'_','longbranches.pdf',sep=''))
	ggsave(file=file, w=10, h=10, useDingbats=FALSE)
	#
	tmp			<- subset(lba.su, OTHER=='N' & TEAM!='RUNGAPS_ExaML')[, list(TAXAN=length(ERR), OUTLIER_P=mean(ERR< -0.04 | ERR> 0.04)), by=c('MODEL','SC','TEAM','GAPS','GENE','IDX')]
	tmp[, OUTLIER_N:=TAXAN*OUTLIER_P]
	set(tmp, NULL, 'TEAM', tmp[, factor(TEAM)])
	set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag','pol','gag+pol+env'))])		
	ggplot(tmp, aes(x=GAPS, y=OUTLIER_P, colour=GENE)) + 
			geom_jitter(position=position_jitter(w=0.8, h = 0), size=1) +
			scale_shape_manual(values=c('IQ-TREE'=15, 'PhyML'=12, 'RAxML'=8, 'MetaPIGA'=17)) +			
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 			
			scale_y_continuous(labels = scales::percent) +
			facet_grid(TEAM~GENE, scales='free_y')  + theme_bw() + theme(legend.position='bottom') +
			labs(	x='\nUnassembled sites in simulated sequences', 
					y='branch length to root\n20% too small or too large\n(% of all taxa in tree)\n')
	file	<- file.path(edir, paste(timetag,'_','longbranches20pc.pdf',sep=''))
	ggsave(file=file, w=10, h=10, useDingbats=FALSE)
	#
	#	severe branch lengths errors by rungaps
	#
	tmp		<- subset(lba.su, TEAM=='RUNGAPS_ExaML' & !grepl('p51', GENE))
	tmp		<- tmp[, list(TAXAN=length(ERR), OUTLIER_P=mean(ERR< -0.04 | ERR>0.04)), by=c('MODEL','SC','TEAM','GAPS','GENE','RUNGAPS','IDX')]
	tmp[, OUTLIER_N:=TAXAN*OUTLIER_P]	
	set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag (p17)','gag','gag+pol+env'))])
	setkey(tmp, GENE, RUNGAPS)	
	tmp2	<- tmp[, 	list(	RUNGAPS= RUNGAPS, 
								OUTLIER_P_SM= predict(loess(OUTLIER_P~RUNGAPS, span=5))), 
						by='GENE']
	tmp		<- merge(tmp, tmp2, by=c('GENE','RUNGAPS'))		
	tmp2	<- merge( rbind(data.table(GENE=c('gag (p17)','gag','gag + pol (prot,p51)','gag+pol+env'), RUNGAPS=c(0.11, 0.08, 0.14, 0.17), LOC='Botswana'), data.table(GENE=c('gag (p17)','gag','gag + pol (prot,p51)','gag+pol+env'), RUNGAPS=c(0.21, 0.18, 0.34, 0.47), LOC='Uganda')), tmp2,by=c('GENE','RUNGAPS')) 
	#set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag (p17)','gag','gag + pol (prot,p51)','gag+pol+env'))])
	ggplot(tmp, aes(x=RUNGAPS)) +
			geom_point(aes(y=OUTLIER_P, colour=GENE), size=2, pch=8) +
			#geom_line(aes(y=OUTLIER_P_SM, colour=GENE), size=0.5) +
			#geom_point(data=tmp2, aes(y=OUTLIER_P_SM, pch=LOC), size=2.5, fill='black') +			
			scale_colour_manual(values=c('gag (p17)'="#8C510A", 'gag + pol (prot,p51)'='green','gag'='red','gag+pol+env'="#3F4788FF")) +
			scale_shape_manual(values=c('Botswana'=23, 'Uganda'=24)) +
			scale_x_continuous(labels = scales::percent, expand=c(0,0), breaks=seq(0,1,0.1)) +
			scale_y_continuous(labels = scales::percent) +
			#scale_shape_manual(values=c('IQ-TREE'=15, 'PhyML'=12, 'RAxML'=8, 'MetaPIGA'=17)) +
			labs(	x='\nUnassembled sites in simulated sequences', 
					y='branch length to root\n20% too small or too large\n(% of all taxa in tree)\n',
					colour='part of genome used\nfor tree reconstruction',
					pch='sampling location') +
			theme_bw() + theme(legend.position='bottom') 
	file	<- file.path(edir, paste(timetag,'_','longbranches_p17full_by_rungaps_taxan1600_Acute10pc.pdf',sep=''))
	ggsave(file=file, w=5, h=7, useDingbats=FALSE)
	#		
	#	MSE by TEAM
	#	
	tmp		<- subset(sa, !is.na(TR_PAIR_rec) & ACUTE=='10%' & TEAM%in%c('IQ-TREE', 'PhyML', 'MetaPIGA', 'RAxML'))
	set(tmp, NULL, 'TEAM', tmp[, factor(TEAM)])
	set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag','pol','gag+pol+env'))])	
	ggplot(subset(tmp, !is.na(MSE_LSD)), aes(x=GAPS)) +
			geom_jitter(aes(y=sqrt(MSE_LSD), colour=GENE, pch=TEAM), position=position_jitter(w=0.8, h = 0), size=2) +
			scale_shape_manual(values=c('IQ-TREE'=15, 'PhyML'=12, 'RAxML'=8, 'MetaPIGA'=17)) +			
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 
			scale_y_log10(limit=c(1,1e4)) +			
			labs(	x='\nUnassembled sites in full-genome sequences', 
					y='root mean squared error\nin dated branches\n(years)\n',					
					colour='part of genome used\nfor tree reconstruction',
					pch='algorithm') +
			theme_bw() + theme(legend.position='bottom') + facet_grid(~TEAM)
	file	<- file.path(edir, paste(timetag,'_','RMSE_by_gaps_by_TEAM.pdf',sep=''))
	ggsave(file=file, w=12, h=6, useDingbats=FALSE)
	#		
	#	MAE overall by TEAM (similar)
	#	
	tmp		<- subset(sa, !is.na(TR_PAIR_rec) & ACUTE=='10%' & TEAM%in%c('IQ-TREE', 'PhyML', 'MetaPIGA', 'RAxML'))
	set(tmp, NULL, 'TEAM', tmp[, factor(TEAM)])
	set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag','pol','gag+pol+env'))])	
	ggplot(subset(tmp, !is.na(MAE_LSD)), aes(x=GAPS)) +
			geom_jitter(aes(y=MAE_LSD, colour=GENE, pch=TEAM), position=position_jitter(w=0.8, h = 0), size=2) +
			scale_shape_manual(values=c('IQ-TREE'=15, 'PhyML'=12, 'RAxML'=8, 'MetaPIGA'=17)) +			
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 
			scale_y_log10(expand=c(0,0), limit=c(1,1e4), breaks=c(1,10,100,1000,1e4), minor_breaks=c(seq(1,10,1),seq(10,100,10),seq(100,1000,100),seq(1000,10000,1000))) +			
			labs(	x='\nUnassembled sites in full-genome sequences', 
					y='mean absolute error in dated branches\n(years)\n',					
					colour='part of genome used\nfor tree reconstruction',
					pch='algorithm') +
			theme_bw() + theme(legend.position='bottom') + facet_grid(~TEAM)
	file	<- file.path(edir, paste(timetag,'_','MAE_by_gaps_by_TEAM.pdf',sep=''))
	ggsave(file=file, w=12, h=6, useDingbats=FALSE)	
	#
	#	MAE overall pairs by by rungaps
	#
	tmp		<- subset(sa, TEAM=='RUNGAPS_ExaML' & !grepl('p51', GENE))
	set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag (p17)','gag','gag+pol+env'))])
	setkey(tmp, GENE, RUNGAPS)	
	tmp2	<- tmp[, 	list(	RUNGAPS= RUNGAPS, 
					MAE_LSD_SM= predict(loess(MAE_LSD~RUNGAPS, degree=2, span=20))), 
			by='GENE']
	tmp		<- merge(tmp, tmp2, by=c('GENE','RUNGAPS'))		
	tmp2	<- merge( rbind(data.table(GENE=c('gag (p17)','gag','gag + pol (prot,p51)','gag+pol+env'), RUNGAPS=c(0.11, 0.08, 0.14, 0.17), LOC='Botswana'), data.table(GENE=c('gag (p17)','gag','gag + pol (prot,p51)','gag+pol+env'), RUNGAPS=c(0.21, 0.18, 0.34, 0.47), LOC='Uganda')), tmp2,by=c('GENE','RUNGAPS')) 
	#set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag (p17)','gag','gag + pol (prot,p51)','gag+pol+env'))])
	ggplot(tmp, aes(x=RUNGAPS)) +
			geom_point(aes(y=MAE_LSD, colour=GENE), size=2, pch=8) +
			geom_line(aes(y=MAE_LSD_SM, colour=GENE), size=0.5) +
			geom_point(data=tmp2, aes(y=MAE_LSD_SM, pch=LOC), size=2.5, fill='black') +			
			scale_colour_manual(values=c('gag (p17)'="#8C510A", 'gag + pol (prot,p51)'='green','gag'='red','gag+pol+env'="#3F4788FF")) +
			scale_shape_manual(values=c('Botswana'=23, 'Uganda'=24)) +
			scale_x_continuous(labels = scales::percent, expand=c(0,0), breaks=seq(0,1,0.1)) +
			scale_y_continuous(expand=c(0,0), limit=c(0,20), breaks=seq(0,20,5), minor_breaks=seq(0,10,1)) +			
			#scale_shape_manual(values=c('IQ-TREE'=15, 'PhyML'=12, 'RAxML'=8, 'MetaPIGA'=17)) +
			labs(	x='\nUnassembled sites in simulated sequences', 
					y='mean absolute error in dated branches\n(years)\n',
					colour='part of genome used\nfor tree reconstruction',
					pch='sampling location') +
			theme_bw() + theme(legend.position='bottom') 
	file	<- file.path(edir, paste(timetag,'_','MAE_p17full_by_rungaps_taxan1600_Acute10pc.pdf',sep=''))
	ggsave(file=file, w=5, h=7, useDingbats=FALSE)
	#		
	#	MAE of transmission pairs by RAxML
	#		
	tmp		<- subset(sa, !is.na(MAE_TP_LSD) & ACUTE=='10%' & TEAM=='RAxML')
	tmp[, list(MAE=mean(MAE_TP_LSD)), by=c('SC','GENE')]
	#		
	#	MAE of transmission pairs by TEAM
	#		
	tmp		<- subset(sa, !is.na(MAE_TP_LSD) & ACUTE=='10%' & GENE=='gag+pol+env' & TEAM%in%c('IQ-TREE', 'PhyML', 'MetaPIGA', 'RAxML'))
	set(tmp, NULL, 'TEAM', tmp[, factor(TEAM)])
	set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag','pol','gag+pol+env'))])	
	ggplot(subset(tmp, !is.na(MAE_TP_LSD)), aes(x=GAPS)) +
			geom_jitter(aes(y=MAE_TP_LSD, colour=GENE, pch=TEAM), position=position_jitter(w=0.8, h = 0), size=2) +
			scale_shape_manual(values=c('IQ-TREE'=15, 'PhyML'=12, 'RAxML'=8, 'MetaPIGA'=17)) +			
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 
			scale_y_log10(expand=c(0,0), limits=c(1,10), breaks=c(1,1.5,2,3,4,5,10)) +			
			labs(	x='\nunassembled sites of PANGEA-HIV sequences', 
					y='mean absolute error (years)\n',					
					colour='part of simulated genome\nused for tree reconstruction',
					pch='tree reconstruction\nmethod') +
			theme_bw() + theme(legend.position='bottom')
	file	<- file.path(edir, paste(timetag,'_','MAETP_by_gaps_by_TEAM.pdf',sep=''))
	ggsave(file=file, w=4.5, h=6, useDingbats=FALSE)
	#		
	#	MAE_LSD of transmission pairs IQTree
	#		
	tmp		<- subset(sa, !is.na(MAE_TP_LSD) & ACUTE=='10%' & TEAM=='IQ-TREE')
	set(tmp, NULL, 'TEAM', tmp[, factor(TEAM)])
	set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag','pol','gag+pol+env'))])	
	ggplot(subset(tmp, !is.na(MAE_TP_LSD)), aes(x=GAPS)) +
			geom_jitter(aes(y=MAE_TP_LSD, colour=GENE), position=position_jitter(w=0.8, h = 0), size=2) +
			scale_shape_manual(values=c('IQ-TREE'=15, 'PhyML'=12, 'RAxML'=8, 'MetaPIGA'=17)) +			
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 
			scale_y_log10(expand=c(0,0), limit=c(1,32), breaks=c(1,10,100,1000), minor_breaks=c(seq(1,10,1),seq(10,100,10),seq(100,1000,100),seq(1000,10000,1000))) +			
			labs(	x='\nUnassembled sites in full-genome sequences', 
					y='mean absolute error in dated branches\namong sampled transmission pairs\n(years)\n',					
					colour='part of genome used\nfor tree reconstruction',
					pch='algorithm') +
			theme_bw() + theme(legend.position='bottom') + facet_grid(~TEAM)
	file	<- file.path(edir, paste(timetag,'_','MAETP_by_gaps_IQTree.pdf',sep=''))
	ggsave(file=file, w=4.5, h=6, useDingbats=FALSE)
	#		
	#	MAE of transmission pairs IQTree
	#		
	tmp		<- subset(sa, !is.na(MAE_TP_LSD) & ACUTE=='10%' & TEAM=='IQ-TREE')
	set(tmp, NULL, 'TEAM', tmp[, factor(TEAM)])
	set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag','pol','gag+pol+env'))])	
	ggplot(subset(tmp, !is.na(MAE_TP)), aes(x=GAPS)) +
			geom_jitter(aes(y=MAE_TP, colour=GENE), position=position_jitter(w=0.8, h = 0), size=2) +
			scale_shape_manual(values=c('IQ-TREE'=15, 'PhyML'=12, 'RAxML'=8, 'MetaPIGA'=17)) +			
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 
			#scale_y_log10(expand=c(0,0), limit=c(1,32), breaks=c(1,10,100,1000), minor_breaks=c(seq(1,10,1),seq(10,100,10),seq(100,1000,100),seq(1000,10000,1000))) +			
			labs(	x='\nUnassembled sites in full-genome sequences', 
					y='mean absolute error in dated branches\namong sampled transmission pairs\n(years)\n',					
					colour='part of genome used\nfor tree reconstruction',
					pch='algorithm') +
			theme_bw() + theme(legend.position='bottom') + facet_grid(~TEAM)
	file	<- file.path(edir, paste(timetag,'_','MAETP_by_gaps_IQTree.pdf',sep=''))
	ggsave(file=file, w=4.5, h=6, useDingbats=FALSE)		
	
	#
	#	MAE of transmission pairs by by rungaps
	#
	tmp		<- subset(sa, TEAM=='RUNGAPS_ExaML' & !grepl('p51', GENE))
	tmp		<- subset(tmp, GENE=='gag' & RUNGAPS==0.02 | GENE=='gag (p17)' & RUNGAPS==0.02 | GENE=='gag+pol+env')
	tmp[, MISSING_P:= (RUNGAPS*GENE_L + (6807-GENE_L))/6807]
	set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag (p17)','gag','gag+pol+env'))])	
	#set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag (p17)','gag','gag + pol (prot,p51)','gag+pol+env'))])
	ggplot(tmp, aes(x=MISSING_P)) +
			geom_point(aes(y=MAE_TP_LSD, colour=GENE), size=2, pch=16) +
			#geom_line(aes(y=MAE_TP_LSD_SM, colour=GENE), size=0.5) +
			#geom_point(data=tmp2, aes(y=MAE_TP_LSD_SM, pch=LOC), size=2.5, fill='black') +			
			scale_colour_manual(values=c('gag (p17)'="#8C510A", 'gag + pol (prot,p51)'='green','gag'='red','gag+pol+env'="#3F4788FF")) +
			scale_shape_manual(values=c('Botswana'=23, 'Uganda'=24)) +
			scale_x_continuous(labels = scales::percent, expand=c(0,0), breaks=seq(0,1,0.1), limits=c(0,1)) +
			scale_y_continuous(expand=c(0,0), limit=c(0,6), breaks=seq(0,10,1), minor_breaks=seq(0,10,.5)) +			
			#scale_shape_manual(values=c('IQ-TREE'=15, 'PhyML'=12, 'RAxML'=8, 'MetaPIGA'=17)) +
			labs(	x='\nProportion of missing sites relative to gag+pol+env genome', 
					y='mean absolute error (years)\n',
					colour='part of simulated genome used\nfor tree reconstruction') +
			theme_bw() + theme(legend.position='bottom') 
	file	<- file.path(edir, paste(timetag,'_','MAETP_p17full_by_missingsites_taxan1600_Acute10pc.pdf',sep=''))
	ggsave(file=file, w=5, h=7, useDingbats=FALSE)	
	#
	#
	tmp		<- subset(sa, TEAM=='RUNGAPS_ExaML' & !grepl('p51', GENE))
	set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag (p17)','gag','gag+pol+env'))])
	setkey(tmp, GENE, RUNGAPS)	
	tmp2	<- tmp[, 	list(	RUNGAPS= RUNGAPS, 
								MAE_TP_LSD_SM= predict(loess(MAE_TP_LSD~RUNGAPS, degree=2, span=5))), 
						by='GENE']
	tmp		<- merge(tmp, tmp2, by=c('GENE','RUNGAPS'))		
	tmp2	<- merge( rbind(data.table(GENE=c('gag (p17)','gag','gag + pol (prot,p51)','gag+pol+env'), RUNGAPS=c(0.11, 0.08, 0.14, 0.17), LOC='Botswana'), data.table(GENE=c('gag (p17)','gag','gag + pol (prot,p51)','gag+pol+env'), RUNGAPS=c(0.21, 0.18, 0.34, 0.47), LOC='Uganda')), tmp2,by=c('GENE','RUNGAPS')) 
	#set(tmp, NULL, 'GENE', tmp[, factor(GENE, levels=c('gag (p17)','gag','gag + pol (prot,p51)','gag+pol+env'))])
	ggplot(tmp, aes(x=RUNGAPS)) +
			geom_point(aes(y=MAE_TP_LSD, colour=GENE), size=2, pch=8) +
			geom_line(aes(y=MAE_TP_LSD_SM, colour=GENE), size=0.5) +
			geom_point(data=tmp2, aes(y=MAE_TP_LSD_SM, pch=LOC), size=2.5, fill='black') +			
			scale_colour_manual(values=c('gag (p17)'="#8C510A", 'gag + pol (prot,p51)'='green','gag'='red','gag+pol+env'="#3F4788FF")) +
			scale_shape_manual(values=c('Botswana'=23, 'Uganda'=24)) +
			scale_x_continuous(labels = scales::percent, expand=c(0,0), breaks=seq(0,1,0.1)) +
			scale_y_continuous(expand=c(0,0), limit=c(0,7), breaks=seq(0,10,2), minor_breaks=seq(0,10,.5)) +			
			#scale_shape_manual(values=c('IQ-TREE'=15, 'PhyML'=12, 'RAxML'=8, 'MetaPIGA'=17)) +
			labs(	x='\nUnassembled sites in simulated sequences', 
					y='mean absolute error in dated branches\namong sampled transmission pairs\n(years)\n',
					colour='part of genome used\nfor tree reconstruction',
					pch='sampling location') +
			theme_bw() + theme(legend.position='bottom') 
	file	<- file.path(edir, paste(timetag,'_','MAETP_p17full_by_rungaps_taxan1600_Acute10pc.pdf',sep=''))
	ggsave(file=file, w=5, h=7, useDingbats=FALSE)	
	#
	#
	#	plot simulated trees versus true tree
	#
	require(ggtree)	
	setkey(sa, SC, TEAM, GENE, RUNGAPS)
	invisible(subset(sa, TEAM!='RUNGAPS_ExaML')[, 
						{
							#IDX		<- c(531,532,533,534,535,536,537,538,539,540,748,749,750,751,752,753,754,755,756,757,870)
							#TEAM	<- c(1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  0,  0,  0,  0,  0,  0 , 0  ,0  ,0  ,0 )
							#GENE	<- rep('POL', length(IDX))
							tmp		<- lapply(IDX, function(i) strs_rtt[[i]] )				
							if(MODEL[1]!='R')
								tmp[[length(tmp)+1]]	<- ttrs[[TIME_IDX_T[1]]]
							if(MODEL[1]=='R')
								tmp[[length(tmp)+1]]	<- ttrs[[SUB_IDX_T[1]]]
							tmp		<- lapply( c(length(tmp), seq(1,length(tmp)-1)), function(i) tmp[[i]] )
							class(tmp) 	<- "multiPhylo"
							print(c('Simulated phylogeny',paste(TEAM, GENE, IDX, sep='-')))
							names(tmp)	<- c('Simulated phylogeny',paste(TEAM, GENE, IDX, sep='-'))
							p	<- ggtree(tmp, size=0.1) + facet_wrap(~.id, ncol=10, scales='free_x')				
							pdf(file=file.path(edir, paste(timetag,'_strs_rtt_',SC,'.pdf',sep='')), w=40, h=length(IDX)/10*12)
							print(p)
							dev.off()	
							NULL
						}, by=c('SC')])				
	invisible(subset(sa, TEAM=='RUNGAPS_ExaML')[, 
					{
						#IDX		<- c(531,532,533,534,535,536,537,538,539,540,748,749,750,751,752,753,754,755,756,757,870)
						#TEAM	<- c(1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  0,  0,  0,  0,  0,  0 , 0  ,0  ,0  ,0 )
						#GENE	<- rep('POL', length(IDX))
						tmp		<- lapply(IDX, function(i) strs_rtt[[i]] )				
						tmp[[length(tmp)+1]]	<- ttrs[[SUB_IDX_T[1]]]
						tmp		<- lapply( c(length(tmp), seq(1,length(tmp)-1)), function(i) tmp[[i]] )
						class(tmp) 	<- "multiPhylo"
						print(c('Simulated phylogeny',paste(TEAM, GENE, IDX, RUNGAPS, sep='-')))
						names(tmp)	<- c('Simulated phylogeny',paste(GENE, IDX, RUNGAPS, sep='-'))
						p	<- ggtree(tmp, size=0.1) + facet_wrap(~.id, ncol=10)				
						pdf(file=file.path(edir, paste(timetag,'_strs_rtt_rungaps_',GENE,'.pdf',sep='')), w=40, h=length(IDX)/10*12)
						print(p)
						dev.off()	
						NULL
					}, by=c('GENE')])	
	#
	#	plot trees
	#
	tmpdir			<- '~/duke/tmp'
	sc_id			<- 'sc 4'
	gene_id			<- 'gag'
	team_id			<- 'PhyML'
	for(sc_id in c("sc 1","sc 2","sc 4"))
	{
		for(team_id in c("PhyML"))
		{
			for(gene_id in c('gag','gag+pol+env'))
			{
				#gene_id<- 'POL'
				tmp				<- subset(sa, TEAM==team_id & SC==sc_id & GENE==gene_id)
				#if(team_id=="MetaPIGA")
				#	tmp			<- subset(tmp, grepl("best solution_use",FILE))
				phr				<- ttrs[[ tmp[1, SUB_IDX_T] ]]
				phs				<- lapply(tmp[, IDX], function(i) strs_rtt[[i]] )
				#	drop to common tips
				z				<- setdiff(phr$tip.label, phs[[1]]$tip.label)
				#stopifnot( length(z)==abs(diff(c(TAXAN, TAXAN_T))) )
				if(length(z))
					phr			<- drop.tip(phr, z, rooted=TRUE, root.edge=1)	
				phr				<- multi2di(phr,random =FALSE)
				pho				<- treedist.get.tree.100bs(phr, phs, tmpdir)
				
				ggtree(pho, aes(color=bootstrap), size=0.2) +
						scale_colour_continuous(low="grey70", high="#1B0C42FF", guide="none") +
						#scale_colour_continuous(low="#1B0C42FF", high="#D64B40FF", guide="none") + 
						theme_tree2(legend.position='right') + theme(axis.line.x=element_line()) + labs(x='average substions per site', axis.title=element_text(size=2))
				file			<- file.path(edir, paste(timetag, '_AgreeTree100bs_', sc_id, '_', team_id, '_', gene_id, 'std.pdf', sep=''))		
				ggsave(file=file, w=6, h=35)				
				ggtree(pho, layout="fan", aes(color=bootstrap), size=0.2) +
						scale_colour_continuous(low="grey70", high="#1B0C42FF", guide="none") +
						#scale_colour_continuous(low="#1B0C42FF", high="#D64B40FF", guide="none") + 
						theme_tree2(legend.position='right') + theme(axis.text=element_blank())
				file			<- file.path(edir, paste(timetag, '_AgreeTree100bs_', gsub(' ','-',sc_id), '_', team_id, '_', gene_id, 'circular.pdf', sep=''))		
				ggsave(file=file, w=6, h=6)
			}
		}		
	}
}	
##--------------------------------------------------------------------------------------------------------
##	olli 27.06.11
##--------------------------------------------------------------------------------------------------------
treecomparison.ana.160627<- function()
{	
	require(ggplot2)
	require(data.table)
	require(ape)
	require(scales)	
	require(ggtree)
	require(phangorn)
	
	edir			<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation'
	#timetag		<- '160627'
	timetag			<- '160713'
	
	load(paste(edir,'/','submitted_160713_RFPDQDTP.rda',sep=''))	
	#
	#
	#
	sa		<- copy(submitted.info)	
	#
	set(sa, NULL, 'MODEL', sa[, factor(MODEL, levels=c('V','R'),labels=c('Model: Village','Model: Regional'))])
	set(sa, sa[, which(SC=="VILL_99_APR15")],'SC',"150701_VILL_SCENARIO-C")	
	set(sa, NULL, 'SC', sa[, factor(SC,	levels=c("150701_REGIONAL_TRAIN1", "150701_REGIONAL_TRAIN2", "150701_REGIONAL_TRAIN3", "150701_REGIONAL_TRAIN4","150701_REGIONAL_TRAIN5","150701_VILL_SCENARIO-A","150701_VILL_SCENARIO-B","150701_VILL_SCENARIO-C","150701_VILL_SCENARIO-D","150701_VILL_SCENARIO-E"), 
							labels=c('sc 1','sc 2','sc 3','sc 4','sc 5','sc A','sc B','sc C','sc D','sc E'))])
	set(sa, NULL, 'GAPS', sa[, factor(GAPS, levels=c('none','low','high'),labels=c('none','as for Botswana\nsequences','as for Uganda\nsequences'))])
	set(sa, NULL, 'BEST', sa[, factor(BEST, levels=c('Y','N'),labels=c('best tree','replicate tree'))])									
	set(sa, NULL, 'GENE', sa[, factor(GENE, levels=c('GAG','POL','GAG+POL+ENV'),labels=c('gag','pol','gag+pol+env'))])	
	set(sa, NULL, 'TEAM', sa[, factor(TEAM, levels=sa[, sort(unique(TEAM))],labels=sa[, sort(unique(TEAM))])])
	set(sa, NULL, 'EXT', sa[, factor(EXT, levels=c('~0pc','5pc'),labels=c('~ 0%/year','5%/year'))])
	set(sa, NULL, 'ACUTE', sa[, factor(ACUTE, levels=c('low','high'),labels=c('10%','40%'))])
	set(sa, NULL, 'ART', sa[, factor(ART, levels=c('none','fast'),labels=c('none','fast'))])
	sa		<- subset(sa, OTHER=='N')
	#
	#	on full tree
	#		
	
	#	prob closest on 4.5% (confounded by branch lengths)
	#tmp		<- merge(sa, subset(lba.su, OUTLIER_P<0.2, IDX), by='IDX')	
	tmp		<- subset(sa, !is.na(TR_REC_perc_45) & GENE=='gag+pol+env' & ACUTE=='10%' & TEAM%in%c('RAxML','IQ-TREE'))
	ggplot(tmp, aes(x=GAPS)) +
			geom_jitter(aes(y=TR_REC_perc_45, colour=GENE, pch=TEAM), position=position_jitter(w=0.8, h = 0), size=2) +
			scale_shape_manual(values=c('IQTree'=15, 'PhyML'=12, 'RAXML'=8, 'MetaPIGA'=17)) +
			geom_point(aes(y=TR_REC_perc_T_45), colour="black", size=2) +
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 
			scale_y_continuous(labels = scales::percent, limit=c(0,0.3), expand=c(0,0)) +
			labs(	x='\nGappiness of full-genome sequences', 
					y='phylogenetically closest pairs of individual\nthat are true transmission pairs\n',
					colour='part of genome used\nfor tree reconstruction',
					pch='algorithm') +
			theme_bw() + theme(legend.position='bottom')
	file	<- file.path(edir, paste(timetag,'_','pTransRec45_by_gaps.pdf',sep=''))
	ggsave(file=file, w=5, h=7)
	#	prob closest w/o GD criterion	
	tmp		<- subset(sa, !is.na(TR_REC_perc_Inf) & ACUTE=='10%' & TEAM%in%c('RAXML','IQTree'))
	ggplot(tmp, aes(x=GAPS)) +
			geom_jitter(aes(y=TR_REC_perc_Inf, colour=GENE, pch=TEAM), position=position_jitter(w=0.8, h = 0), size=2) +
			scale_shape_manual(values=c('IQTree'=15, 'PhyML'=12, 'RAXML'=8, 'MetaPIGA'=17)) +
			geom_point(aes(y=TR_REC_perc_T_Inf), colour="black", size=2) +
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 
			scale_y_continuous(labels = scales::percent, limit=c(0,0.1), expand=c(0,0)) +
			labs(	x='\nGappiness of full-genome sequences', 
					y='phylogenetically closest pairs of individuals\nthat are true transmission pairs\n',
					colour='part of genome used\nfor tree reconstruction',
					pch='algorithm') +
			theme_bw() + theme(legend.position='bottom')
	file	<- file.path(edir, paste(timetag,'_','pTransRecInf_by_gaps.pdf',sep=''))
	ggsave(file=file, w=5, h=7)
	#	instead, proportion of recovered transmission pairs? 
	#	get list of correct pairs in true phylogeny. how many of these do we see in reconstructed phylogeny?
	tmp		<- subset(sa, !is.na(TR_PAIR_rec) & ACUTE=='10%')
	ggplot(tmp, aes(x=GAPS)) +
			geom_jitter(aes(y=TR_PAIR_rec, colour=GENE, pch=TEAM), position=position_jitter(w=0.8, h = 0), size=2) +
			scale_shape_manual(values=c('IQTree'=15, 'PhyML'=12, 'RAXML'=8, 'MetaPIGA'=17)) +			
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 
			scale_y_continuous(labels = scales::percent, limit=c(0.5,1)) +
			labs(	x='\nGappiness of full-genome sequences', 
					y='phylogenetically closest pairs of individuals\nthat are transmission pairs, out of all such pairs\nthat can be identified in the true tree\n',					
					colour='part of genome used\nfor tree reconstruction',
					pch='algorithm') +
			theme_bw() + theme(legend.position='bottom')
	file	<- file.path(edir, paste(timetag,'_','pTransPairRecovered_by_gaps.pdf',sep=''))
	ggsave(file=file, w=6, h=7)
	


	ggplot(subset(sa, ACUTE=='10%' & TEAM!='MetaPIGA' & MODEL=='Model: Regional'), aes(x=TAXAN)) +
			geom_jitter(aes(y=NRF, colour=GENE, pch=TEAM), position=position_jitter(w=0.8, h = 0), size=2) +			
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 
			scale_y_continuous(labels = scales::percent, expand=c(0,0), limits=c(0, 1)) +
			scale_shape_manual(values=c('IQTree'=15, 'PhyML'=12, 'RAXML'=8)) +
			labs(	x='\nnumber of taxa in simulated tree', 
					y='RF distance\n(standardized)\n',
					colour='part of genome used\nfor tree reconstruction',
					pch='algorithm') +
			theme_bw() + theme(legend.position='bottom') +
			facet_grid(~GAPS)
	file	<- file.path(edir, paste(timetag,'_','RF_fulltree_polvsall_by_gaps_taxan_RegionalAcute10pc.pdf',sep=''))
	ggsave(file=file, w=5, h=7)
	
	ggplot(subset(sa, ACUTE=='40%' & TEAM!='MetaPIGA'), aes(x=TAXAN)) +
			geom_jitter(aes(y=NRF, colour=GENE, pch=TEAM), position=position_jitter(w=0.8, h = 0), size=2) +			
			scale_colour_manual(values=c('pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 
			scale_y_continuous(labels = scales::percent, expand=c(0,0), limits=c(0, 1)) +
			scale_shape_manual(values=c('IQTree'=15, 'PhyML'=12, 'RAXML'=8)) +
			labs(	x='\nnumber of taxa in simulated tree', 
					y='RF distance\n(standardized)\n',
					colour='part of genome used\nfor tree reconstruction',
					pch='algorithm') +
			theme_bw() + theme(legend.position='bottom') +
			facet_grid(~GAPS)
	file	<- file.path(edir, paste(timetag,'_','RF_fulltree_polvsall_by_gaps_taxan_AllAcute40pc.pdf',sep=''))
	ggsave(file=file, w=5, h=7)
	
	ggplot(subset(sa, TEAM!='MetaPIGA'), aes(x=TAXAN)) +
			geom_jitter(aes(y=NQD, colour=GENE, pch=TEAM), position=position_jitter(w=0.8, h = 0), size=2) +			
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 
			scale_y_continuous(labels = scales::percent, expand=c(0,0), limits=c(0, 1)) +
			scale_shape_manual(values=c('IQTree'=15, 'PhyML'=12, 'RAXML'=8)) +
			labs(	x='\nnumber of taxa in simulated tree', 
					y='Quartet distance\n(standardized)\n',
					colour='part of genome used\nfor tree reconstruction',
					pch='algorithm') +
			theme_bw() + theme(legend.position='bottom') +
			facet_grid(~GAPS)
	file	<- file.path(edir, paste(timetag,'_','QD_polvsall_by_gaps_taxan_AllAcute40pc.pdf',sep=''))
	ggsave(file=file, w=5, h=7)
	
	ggplot(subset(sa, TEAM!='MetaPIGA' & ACUTE=='10%')) +
			geom_jitter(aes(x=GAPS, y=PD, colour=GENE, pch=TEAM), position=position_jitter(w=0.8, h = 0), size=2) +	
			#geom_boxplot(aes(x=cut(TAXAN, breaks=seq(800,1601,200)), y=NPD, colour=GENE)) +
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 
			#scale_y_continuous(labels = scales::percent, expand=c(0,0), limits=c(0, 0.03)) +
			scale_shape_manual(values=c('IQTree'=15, 'PhyML'=12, 'RAXML'=8)) +
			labs(	x='\nGappiness of full-genome sequences', 
					y='path distance\n',
					colour='part of genome used\nfor tree reconstruction',
					pch='algorithm') +
			theme_bw() + theme(legend.position='bottom') +
			facet_grid(TEAM~.)
	file	<- file.path(edir, paste(timetag,'_','PD_polvsall_by_gaps_taxan1600_Acute10pc.pdf',sep=''))
	ggsave(file=file, w=5, h=10)	
	
}
##--------------------------------------------------------------------------------------------------------
##	olli 03.12.15
##--------------------------------------------------------------------------------------------------------
treecomparison.ana.160502<- function()
{	
	require(ggplot2)
	require(data.table)
	require(ape)
	require(scales)	
	require(ggtree)
	require(phangorn)

	edir			<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation'
	#timetag			<- '160502'
	#load(paste(edir,'/','submitted_160430.rda',sep=''))
	timetag			<- '160627'	
	load(paste(edir,'/','submitted_160627.rda',sep=''))
	
	#
	#	plot trees
	#
	tmpdir			<- '~/duke/tmp'
	sc_id			<- '150701_REGIONAL_TRAIN5'
	gene_id			<- 'GAG+POL+ENV'
	team_id			<- 'RAXML'
	for(sc_id in c("150701_REGIONAL_TRAIN1","150701_REGIONAL_TRAIN2","150701_REGIONAL_TRAIN3","150701_REGIONAL_TRAIN4","150701_REGIONAL_TRAIN5"))
	{
		for(team_id in c("IQTree","RAXML"))
		{
			for(gene_id in submitted.info[, unique(GENE)])
			{
				#gene_id<- 'POL'
				tmp				<- subset(submitted.info, TEAM==team_id & SC==sc_id & GENE==gene_id)
				if(team_id=="MetaPIGA")
					tmp			<- subset(tmp, grepl("best solution_use",FILE))
				phr				<- ttrs[[ tmp[1, SUB_IDX_T] ]]
				phs				<- lapply(tmp[, IDX], function(i) strs[[i]] )
				#	drop to common tips
				z				<- setdiff(phr$tip.label, phs[[1]]$tip.label)
				#stopifnot( length(z)==abs(diff(c(TAXAN, TAXAN_T))) )
				if(length(z))
					phr			<- drop.tip(phr, z, rooted=TRUE, root.edge=1)	
				phr				<- multi2di(phr,random =FALSE)
				pho				<- treedist.get.tree.100bs(phr, phs, tmpdir)
				
				ggtree(pho, layout="fan", aes(color=bootstrap), size=0.2) +
						scale_colour_continuous(low="grey70", high="#1B0C42FF", guide="none") +
						#scale_colour_continuous(low="#1B0C42FF", high="#D64B40FF", guide="none") + 
						theme_tree2(legend.position='right') + theme(axis.text=element_blank())
				file			<- file.path(edir, paste(timetag, '_AgreeTree100bs_', sc_id, '_', team_id, '_', gene_id, '.pdf', sep=''))		
				ggsave(file=file, w=6, h=6)
			}
		}		
	}
	#	average gappiness
	
	#
	#	plot prob closest
	#
	sa		<- copy(submitted.info)	
	#
	set(sa, NULL, 'MODEL', sa[, factor(MODEL, levels=c('V','R'),labels=c('Model: Village','Model: Regional'))])
	set(sa, sa[, which(SC=="VILL_99_APR15")],'SC',"150701_VILL_SCENARIO-C")	
	set(sa, NULL, 'SC', sa[, factor(SC,	levels=c("150701_REGIONAL_TRAIN1", "150701_REGIONAL_TRAIN2", "150701_REGIONAL_TRAIN3", "150701_REGIONAL_TRAIN4","150701_REGIONAL_TRAIN5","150701_VILL_SCENARIO-A","150701_VILL_SCENARIO-B","150701_VILL_SCENARIO-C","150701_VILL_SCENARIO-D","150701_VILL_SCENARIO-E"), 
							labels=c('sc 1','sc 2','sc 3','sc 4','sc 5','sc A','sc B','sc C','sc D','sc E'))])
	set(sa, NULL, 'GAPS', sa[, factor(GAPS, levels=c('none','low','high'),labels=c('none','as for Botswana\nsequences','as for Uganda\nsequences'))])
	set(sa, NULL, 'BEST', sa[, factor(BEST, levels=c('Y','N'),labels=c('best tree','replicate tree'))])									
	set(sa, NULL, 'GENE', sa[, factor(GENE, levels=c('POL','GAG+POL+ENV'),labels=c('pol','gag+pol+env'))])	
	set(sa, NULL, 'TEAM', sa[, factor(TEAM, levels=sa[, sort(unique(TEAM))],labels=sa[, sort(unique(TEAM))])])
	set(sa, NULL, 'EXT', sa[, factor(EXT, levels=c('~0pc','5pc'),labels=c('~ 0%/year','5%/year'))])
	set(sa, NULL, 'ACUTE', sa[, factor(ACUTE, levels=c('low','high'),labels=c('10%','40%'))])
	set(sa, NULL, 'ART', sa[, factor(ART, levels=c('none','fast'),labels=c('none','fast'))])
	sa		<- subset(sa, OTHER=='N')
	#
	subset(sa, MODEL=='Model: Regional')[, table(TEAM, GENE, SC)]
	#		
	ggplot(subset(sa, !is.na(TR_REC_perc) & ACUTE=='10%'), aes(x=GAPS)) +
		geom_jitter(aes(y=TR_REC_perc, colour=GENE, pch=TEAM), position=position_jitter(w=0.8, h = 0), size=2) +
		scale_shape_manual(values=c('IQTree'=15, 'PhyML'=12, 'RAXML'=8, 'MetaPIGA'=17)) +
		geom_point(aes(y=TR_REC_perc_T), colour="#D64B40FF", size=2) +
		scale_colour_manual(values=c('pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 
		scale_y_continuous(labels = scales::percent, limit=c(0,0.1), expand=c(0,0)) +
		labs(	x='\nGappiness of full-genome sequences', 
				y='phylogenetically closest individual is\ntransmitter or next infected\n',
				colour='part of genome used\nfor tree reconstruction',
				pch='algorithm') +
		theme_bw() + theme(legend.position='bottom')
	file	<- file.path(edir, paste(timetag,'_','pTransRec_by_gaps.pdf',sep=''))
	ggsave(file=file, w=5, h=7)
	#
	#	Quartett distance
	#
	load(paste(edir,'/','submitted_151101_BLQDKC.rda',sep=''))
	sc		<- copy(sclu.info)
	#
	tmp		<- subset(tinfo, !is.na(IDCLU))[, list(CLU_N=CLU_N[1], MXGPS_CLU= max(GPS), MDGPS_CLU=median(GPS)), by=c('SC','IDCLU')]
	sc		<- merge(sc, tmp, by=c('SC','IDCLU'))	
	set(sc, NULL, 'MODEL', sc[, factor(MODEL, levels=c('V','R'),labels=c('Model: Village','Model: Regional'))])
	set(sc, sc[, which(SC=="VILL_99_APR15")],'SC',"150701_VILL_SCENARIO-C")	
	set(sc, NULL, 'SC', sc[, factor(SC,	levels=c("150701_REGIONAL_TRAIN1", "150701_REGIONAL_TRAIN2", "150701_REGIONAL_TRAIN3", "150701_REGIONAL_TRAIN4","150701_REGIONAL_TRAIN5","150701_VILL_SCENARIO-A","150701_VILL_SCENARIO-B","150701_VILL_SCENARIO-C","150701_VILL_SCENARIO-D","150701_VILL_SCENARIO-E"), 
							labels=c('sc 1','sc 2','sc 3','sc 4','sc 5','sc A','sc B','sc C','sc D','sc E'))])
	set(sc, NULL, 'GAPS', sc[, factor(GAPS, levels=c('none','low','high'),labels=c('none','as for Botswana\nsequences','as for Uganda\nsequences'))])
	set(sc, NULL, 'BEST', sc[, factor(BEST, levels=c('Y','N'),labels=c('best tree','replicate tree'))])									
	set(sc, NULL, 'GENE', sc[, factor(GENE, levels=c('POL','GAG+POL+ENV'),labels=c('pol','gag+pol+env'))])	
	set(sc, NULL, 'TEAM', sc[, factor(TEAM, levels=sc[, sort(unique(TEAM))],labels=sc[, sort(unique(TEAM))])])
	set(sc, NULL, 'EXT', sc[, factor(EXT, levels=c('~0pc','5pc'),labels=c('~ 0%/year','5%/year'))])
	set(sc, NULL, 'ART', sc[, factor(ART, levels=c('none','fast'),labels=c('none','fast'))])
	sc		<- subset(sc, OTHER=='N')	
	sc		<- sc[, list( SB_NQD=mean(NQDC, na.rm=TRUE) ), by=c('SC','GENE','TEAM','BEST','IDX','FILE','GAPS','MODEL','TAXAN','TAXAN_T','ROOTED','SEQCOV','ART','ACUTE','EXT','OTHER')]
	sc		<- subset(sc, MODEL=='Model: Regional')
	
	ggplot(subset(sc, ACUTE=='low' & TEAM!='MetaPIGA'), aes(x=GAPS)) +
			geom_jitter(aes(y=SB_NQD, colour=GENE, pch=TEAM), position=position_jitter(w=0.8, h = 0), size=2) +			
			scale_colour_manual(values=c('pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 
			scale_y_continuous(labels = scales::percent, expand=c(0,0), limits=c(0, 0.4)) +
			scale_shape_manual(values=c('IQTree'=15, 'PhyML'=12, 'RAXML'=8)) +
			labs(	x='\nGappiness of full-genome sequences', 
					y='Quartett distance\n(standardized)\n',
					colour='part of genome used\nfor tree reconstruction',
					pch='algorithm') +
			theme_bw() + theme(legend.position='bottom')
	file	<- file.path(edir, paste(timetag,'_','QD_polvsall_by_gaps.pdf',sep=''))
	ggsave(file=file, w=5, h=7)
	
	
}
##--------------------------------------------------------------------------------------------------------
##	olli 03.12.15
##--------------------------------------------------------------------------------------------------------
treecomparison.ana.151203<- function()
{
	require(ggplot2)
	require(gamlss)
	require(scales)
	
	edir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation'
	timetag	<- '151203'
	load(paste(edir,'/','submitted_151203.rda',sep=''))
	#
	#
	#
	sa		<- copy(submitted.info)
	sc		<- copy(sclu.info)
	#
	set(sa, NULL, 'MODEL', sa[, factor(MODEL, levels=c('V','R'),labels=c('Model: Village','Model: Regional'))])
	set(sa, sa[, which(SC=="VILL_99_APR15")],'SC',"150701_VILL_SCENARIO-C")	
	set(sa, NULL, 'SC', sa[, factor(SC,	levels=c("150701_REGIONAL_TRAIN1", "150701_REGIONAL_TRAIN2", "150701_REGIONAL_TRAIN3", "150701_REGIONAL_TRAIN4","150701_REGIONAL_TRAIN5","150701_VILL_SCENARIO-A","150701_VILL_SCENARIO-B","150701_VILL_SCENARIO-C","150701_VILL_SCENARIO-D","150701_VILL_SCENARIO-E"), 
							labels=c('sc 1','sc 2','sc 3','sc 4','sc 5','sc A','sc B','sc C','sc D','sc E'))])
	set(sa, NULL, 'GAPS', sa[, factor(GAPS, levels=c('none','low','high'),labels=c('Gaps: none','Gaps: low','Gaps: high'))])
	set(sa, NULL, 'BEST', sa[, factor(BEST, levels=c('Y','N'),labels=c('best tree','replicate tree'))])									
	set(sa, NULL, 'GENE', sa[, factor(GENE, levels=c('POL','GAG+POL+ENV'),labels=c('pol','gag+pol+env'))])	
	set(sa, NULL, 'TEAM', sa[, factor(TEAM, levels=sa[, sort(unique(TEAM))],labels=sa[, sort(unique(TEAM))])])
	set(sa, NULL, 'EXT', sa[, factor(EXT, levels=c('~0pc','5pc'),labels=c('~ 0%/year','5%/year'))])
	set(sa, NULL, 'ACUTE', sa[, factor(ACUTE, levels=c('low','high'),labels=c('10%','40%'))])
	set(sa, NULL, 'ART', sa[, factor(ART, levels=c('none','fast'),labels=c('none','fast'))])
	sa		<- subset(sa, OTHER=='N')
	#
	tmp		<- subset(tinfo, !is.na(IDCLU))[, list(CLU_N=CLU_N[1], MXGPS_CLU= max(GPS), MDGPS_CLU=median(GPS)), by=c('SC','IDCLU')]
	sc		<- merge(sc, tmp, by=c('SC','IDCLU'))	
	set(sc, NULL, 'MODEL', sc[, factor(MODEL, levels=c('V','R'),labels=c('Model: Village','Model: Regional'))])
	set(sc, sc[, which(SC=="VILL_99_APR15")],'SC',"150701_VILL_SCENARIO-C")	
	set(sc, NULL, 'SC', sc[, factor(SC,	levels=c("150701_REGIONAL_TRAIN1", "150701_REGIONAL_TRAIN2", "150701_REGIONAL_TRAIN3", "150701_REGIONAL_TRAIN4","150701_REGIONAL_TRAIN5","150701_VILL_SCENARIO-A","150701_VILL_SCENARIO-B","150701_VILL_SCENARIO-C","150701_VILL_SCENARIO-D","150701_VILL_SCENARIO-E"), 
							labels=c('sc 1','sc 2','sc 3','sc 4','sc 5','sc A','sc B','sc C','sc D','sc E'))])
	set(sc, NULL, 'GAPS', sc[, factor(GAPS, levels=c('none','low','high'),labels=c('Gaps: none','Gaps: low','Gaps: high'))])
	set(sc, NULL, 'BEST', sc[, factor(BEST, levels=c('Y','N'),labels=c('best tree','replicate tree'))])									
	set(sc, NULL, 'GENE', sc[, factor(GENE, levels=c('POL','GAG+POL+ENV'),labels=c('pol','gag+pol+env'))])	
	set(sc, NULL, 'TEAM', sc[, factor(TEAM, levels=sc[, sort(unique(TEAM))],labels=sc[, sort(unique(TEAM))])])
	set(sc, NULL, 'EXT', sc[, factor(EXT, levels=c('~0pc','5pc'),labels=c('~ 0%/year','5%/year'))])
	set(sc, NULL, 'ART', sc[, factor(ART, levels=c('none','fast'),labels=c('none','fast'))])
	sc		<- subset(sc, OTHER=='N')
	#
	stopifnot(sc[, !any(is.na(SB_NRFC))], sc[, !any(is.na(SB_NQDC))])
	scp		<- sc[, list( SB_NRF=mean(SB_NRFC, na.rm=TRUE), SB_NQD=mean(SB_NQDC, na.rm=TRUE) ), by=c('SC','GENE','TEAM','BEST','IDX','FILE','GAPS','MODEL','TAXAN','TAXAN_T','ROOTED','SEQCOV','ART','ACUTE','EXT','OTHER','TIME_IDX_T','SUB_IDX_T')]
	if('BILL'%in%colnames(sc))
	{
		tmp		<- sc[, list( BILL=mean(BILL, na.rm=TRUE) ), by=c('SC','GENE','TEAM','BEST','IDX','FILE','GAPS','MODEL','TAXAN','TAXAN_T','ROOTED','SEQCOV','ART','ACUTE','EXT','OTHER','TIME_IDX_T','SUB_IDX_T')]
		scp		<- merge(scp, tmp, by=c('SC','GENE','TEAM','BEST','IDX','FILE','GAPS','MODEL','TAXAN','TAXAN_T','ROOTED','SEQCOV','ART','ACUTE','EXT','OTHER','TIME_IDX_T','SUB_IDX_T') )		
	}	
	if('LSD_NRFC'%in%colnames(sc))
	{
		tmp		<- sc[, list( LSD_NRF=mean(LSD_NRFC, na.rm=TRUE) ), by=c('SC','GENE','TEAM','BEST','IDX','FILE','GAPS','MODEL','TAXAN','TAXAN_T','ROOTED','SEQCOV','ART','ACUTE','EXT','OTHER','TIME_IDX_T','SUB_IDX_T')]
		scp		<- merge(scp, tmp, by=c('SC','GENE','TEAM','BEST','IDX','FILE','GAPS','MODEL','TAXAN','TAXAN_T','ROOTED','SEQCOV','ART','ACUTE','EXT','OTHER','TIME_IDX_T','SUB_IDX_T') )		
	}
	sm		<- rbind( subset(sa, grepl('Village',MODEL)), scp, fill=TRUE, use.names=TRUE)
	#
	#	PRIMARY OBJECTIVE
	#
	#	polvsall by gaps
	#	
	if('NRF'%in%colnames(sm))
	{
		ggplot( subset(sm, TEAM!='MetaPIGA'), aes(y=SB_NRF, x=SC, shape=TEAM, fill=GENE, colour=GENE, size=BEST) ) + 
				geom_jitter(position = position_jitter(height=.001, width=0.2)) +			
				scale_size_manual(values=c(3, 1)) +
				scale_shape_manual(values=c(21,23,24)) +
				scale_fill_brewer(palette='Paired') +
				scale_colour_brewer(palette='Paired') +
				facet_wrap(MODEL~GAPS, scales='free_x') +	
				labs(x='\nsimulated data set', y='Robinson-Fould distance\nof estimated trees with subst/site branches\n(standardized)\n', size='', shape='Method', fill='part of genome', colour='part of genome') +
				theme_bw() 
		ggsave(w=10, h=6, file=paste(edir,'/',timetag,'_RF_SUBST_polvsall_by_gaps.pdf',sep=''))	
		ggplot( subset(sm, TEAM!='MetaPIGA'), aes(y=LSD_NRF, x=SC, shape=TEAM, fill=GENE, colour=GENE, size=BEST) ) + 
				geom_jitter(position = position_jitter(height=.001, width=0.2)) +			
				scale_size_manual(values=c(3, 1)) +
				scale_shape_manual(values=c(21,23,24)) +
				scale_fill_brewer(palette='Paired') +
				scale_colour_brewer(palette='Paired') +
				facet_wrap(MODEL~GAPS, scales='free_x') +	
				labs(x='\nsimulated data set', y='Robinson-Fould distance\nof estimated trees with dated branches\n(standardized)\n', size='', shape='Method', fill='part of genome', colour='part of genome') +
				theme_bw() 
		ggsave(w=10, h=6, file=paste(edir,'/',timetag,'_RF_DATED_polvsall_by_gaps.pdf',sep=''))	
	}
	if('LSD_KC_L1'%in%colnames(sa))
	{		
		ggplot( subset(sa, TEAM!='MetaPIGA'), aes(y=LSD_KC_L1/TAXAN/TAXAN, x=SC) ) + 			
				geom_boxplot(aes(colour=GENE), fill='transparent', size=0.5, outlier.shape=NA, alpha=0.3) +
				geom_jitter(aes(shape=TEAM, fill=GENE, colour=GENE, size=BEST), position = position_jitter(height=.001, width=0.2)) +
				scale_size_manual(values=c(3, 1)) +
				scale_shape_manual(values=c(21,23,24)) +
				scale_fill_brewer(palette='Paired') +
				scale_colour_brewer(palette='Paired') +
				facet_wrap(MODEL~GAPS, scales='free_x') +	
				labs(x='\nsimulated data set', y='Kendall-Colijn\nof estimated trees with dated branches\n(lambda=1, /TX^2)\n', size='', shape='Method', fill='part of genome', colour='part of genome') +
				theme_bw() 
		ggsave(w=10, h=6, file=paste(edir,'/',timetag,'_KC1_DATED_polvsall_by_gaps.pdf',sep=''))
		ggplot( subset(sa, TEAM!='MetaPIGA'), aes(y=LSD_KC_L0/TAXAN/TAXAN, x=SC) ) + 			
				geom_boxplot(aes(colour=GENE), fill='transparent', size=0.5, outlier.shape=NA, alpha=0.3) +
				geom_jitter(aes(shape=TEAM, fill=GENE, colour=GENE, size=BEST), position = position_jitter(height=.001, width=0.2)) +
				scale_size_manual(values=c(3, 1)) +
				scale_shape_manual(values=c(21,23,24)) +
				scale_fill_brewer(palette='Paired') +
				scale_colour_brewer(palette='Paired') +
				facet_wrap(MODEL~GAPS, scales='free_x') +	
				labs(x='\nsimulated data set', y='Kendall-Colijn\nof estimated trees with subst/site branches\n(lambda=0, /TX^2)\n', size='', shape='Method', fill='part of genome', colour='part of genome') +
				theme_bw() 
		ggsave(w=10, h=6, file=paste(edir,'/',timetag,'_KC0_SUBST_polvsall_by_gaps.pdf',sep=''))		
	}
	#
	# correlation between KC0 and KC1
	#
	ggplot( subset(sa, TEAM!='MetaPIGA'), aes(x=LSD_KC_L0, y=LSD_KC_L1) ) +			
			geom_point(aes(shape=TEAM, colour=SC, size=BEST)) + geom_abline(slope=1, intercept=0) +
			scale_colour_brewer(palette='Paired') + scale_size_manual(values=c(3, 1)) + scale_shape_manual(values=c(21,23,24)) +
			facet_grid(MODEL~GAPS) +
			labs(y='Kendall-Colijn\nof estimated trees with dated branches\n(lambda=1, branch lengths)\n', x='Kendall-Colijn\nof estimated trees with dated branches\n(lambda=0, topology)\n', size='', shape='Method', colour='simulated data set') +
			theme_bw() 
	ggsave(w=10, h=6, file=paste(edir,'/',timetag,'_KC0_KC1_correlation.pdf',sep=''))
	#
	# KC discrepancy due to taxa size? .. cannot rule out ..
	#
	ggplot( subset(sa, MODEL=="Model: Village" & TEAM!='MetaPIGA'), aes(y=LSD_KC_L0/TAXAN/TAXAN, x=SC, shape=TEAM, fill=TAXAN, colour=TAXAN, size=BEST) ) + 
			geom_jitter(position = position_jitter(height=.001, width=0.2)) +			
			scale_size_manual(values=c(3, 1)) +
			scale_shape_manual(values=c(21,23,24)) +
			scale_colour_gradientn(colours = rainbow(7)) +
			scale_fill_gradientn(colours = rainbow(7)) +						
			facet_grid(GENE~GAPS) +	
			labs(title="Model: Village\n", x='\nsimulated data set', y='Kendall-Colijn\nof estimated trees with dated branches\n(lambda=0, topology)\n', size='', shape='Method', fill='Taxa in subm tree', colour='Taxa in subm tree') +
			theme_bw() 	
	ggsave(w=10, h=6, file=paste(edir,'/',timetag,'_KC0_polvsall_Village_by_TAXAN.pdf',sep=''))
	ggplot( subset(sa, MODEL=="Model: Village" & TEAM!='MetaPIGA'), aes(y=LSD_NRF, x=SC, shape=TEAM, fill=TAXAN, colour=TAXAN, size=BEST) ) + 
			geom_jitter(position = position_jitter(height=.001, width=0.2)) +			
			scale_size_manual(values=c(3, 1)) +
			scale_shape_manual(values=c(21,23,24)) +
			scale_colour_gradientn(colours = rainbow(7)) +
			scale_fill_gradientn(colours = rainbow(7)) +						
			facet_grid(GENE~GAPS) +	
			labs(title="Model: Village\n", x='\nsimulated data set', y='Robinson-Fould\nof estimated trees with dated branches\n', size='', shape='Method', fill='Taxa in subm tree', colour='Taxa in subm tree') +
			theme_bw() 	
	ggsave(w=10, h=6, file=paste(edir,'/',timetag,'_RF_polvsall_Village_by_TAXAN.pdf',sep=''))
	#
	# not an issue for RF because we see same signal in regional where taxa size is constant
	#
	ggplot( subset(sa, MODEL=="Model: Regional" & TEAM!='MetaPIGA'), aes(y=NRF, x=SC, shape=TEAM, fill=TAXAN, colour=TAXAN, size=BEST) ) + 
			geom_jitter(position = position_jitter(height=.001, width=0.2)) +			
			scale_size_manual(values=c(3, 1)) +
			scale_shape_manual(values=c(21,23,24)) +
			#scale_fill_brewer(palette='Paired') +
			#scale_colour_brewer(palette='Paired') +
			facet_wrap(GENE~GAPS, scales='free_x', ncol=3) +	
			labs(title="Model: Regional\n", x='\nsimulated data set', y='Robinson-Fould\n(standardized)\n', size='', shape='Method', fill='Taxa in subm tree', colour='Taxa in subm tree') +
			theme_bw() 
	ggsave(w=10, h=6, file=paste(edir,'/',timetag,'_RF_polvsall_by_TAXAN_R.pdf',sep=''))
	#
	# plot 2D for RF
	#
	ggplot(subset(sa, SC=='sc 5'), aes(x=LSD_RF_MDSx/(2*TAXAN-6), y=LSD_RF_MDSy/(2*TAXAN-6))) + 
			geom_point(aes(colour=TEAM, shape=GENE, size=BEST)) +
			geom_point(x=0, y=0, colour='black') +
			scale_shape_manual(values=c(17,18)) + scale_size_manual(values=c(3, 1)) +
			scale_fill_brewer(palette='Paired') + scale_colour_brewer(palette='Set1') +
			labs(x='', y='') +
			theme_bw()
	

}
##--------------------------------------------------------------------------------------------------------
##	olli 19.11.15
##--------------------------------------------------------------------------------------------------------
treecomparison.ana.151019<- function()
{
	require(ggplot2)
	require(gamlss)
	
	edir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation'	
	if(0)
	{		
		timetag	<- '151101'
		file	<- paste(edir,'/','submitted_',timetag,'.rda',sep='')
		file	<- paste(edir,'/','submitted_',timetag,'_BLQDKC.rda',sep='')	
		load(file)
		sa		<- copy(submitted.info)
		sc		<- copy(sclu.info)
	}
	if(0)
	{
		file	<- paste(edir,'/','submitted_151016_CC.rda',sep='')
		load(file)
		sa		<- copy(myinfo)		
		sc		<- copy(sclu.info)
	}
	if(0)
	{
		timetag			<- '151101'		
		file			<- paste(edir,'/','submitted_',timetag,'_QD.rda',sep='')	
		load(file)
		tmp				<- copy(submitted.info)
		tmp2			<- copy(sclu.info)
		file			<- paste(edir,'/','submitted_',timetag,'_BL.rda',sep='')
		load(file)
		submitted.info	<- merge(submitted.info, subset(tmp, select=c('IDX','NQD')), by='IDX')
		sclu.info[, BILL.x:=NULL]
		setnames(sclu.info, 'BILL.y','BILL')
		sclu.info		<- merge(sclu.info, subset(tmp2, select=c('IDX','IDCLU','NQDC')), by=c('IDX','IDCLU'))
		outfile	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation/submitted_151101_BLQD.rda'
		save(strs, ttrs, tinfo, submitted.info, sclu.info, file=outfile)
	}
	
	set(sa, NULL, 'MODEL', sa[, factor(MODEL, levels=c('V','R'),labels=c('Model: Village','Model: Regional'))])
	set(sa, sa[, which(SC=="VILL_99_APR15")],'SC',"150701_VILL_SCENARIO-C")	
	set(sa, NULL, 'SC', sa[, factor(SC,	levels=c("150701_REGIONAL_TRAIN1", "150701_REGIONAL_TRAIN2", "150701_REGIONAL_TRAIN3", "150701_REGIONAL_TRAIN4","150701_REGIONAL_TRAIN5","150701_VILL_SCENARIO-A","150701_VILL_SCENARIO-B","150701_VILL_SCENARIO-C","150701_VILL_SCENARIO-D","150701_VILL_SCENARIO-E"), 
										labels=c('sc 1','sc 2','sc 3','sc 4','sc 5','sc A','sc B','sc C','sc D','sc E'))])
	set(sa, NULL, 'GAPS', sa[, factor(GAPS, levels=c('none','low','high'),labels=c('Gaps: none','Gaps: low','Gaps: high'))])
	set(sa, NULL, 'BEST', sa[, factor(BEST, levels=c('Y','N'),labels=c('best tree','replicate tree'))])									
	set(sa, NULL, 'GENE', sa[, factor(GENE, levels=c('POL','GAG+POL+ENV'),labels=c('pol','gag+pol+env'))])	
	set(sa, NULL, 'TEAM', sa[, factor(TEAM, levels=sa[, sort(unique(TEAM))],labels=sa[, sort(unique(TEAM))])])
	set(sa, NULL, 'EXT', sa[, factor(EXT, levels=c('~0pc','5pc'),labels=c('~ 0%/year','5%/year'))])
	set(sa, NULL, 'ACUTE', sa[, factor(ACUTE, levels=c('low','high'),labels=c('10%','40%'))])
	set(sa, NULL, 'ART', sa[, factor(ART, levels=c('none','fast'),labels=c('none','fast'))])
	sa		<- subset(sa, OTHER=='N')
	
	tmp		<- subset(tinfo, !is.na(IDCLU))[, list(CLU_N=CLU_N[1], MXGPS_CLU= max(GPS), MDGPS_CLU=median(GPS)), by=c('SC','IDCLU')]
	sc		<- merge(sc, tmp, by=c('SC','IDCLU'))	
	set(sc, NULL, 'MODEL', sc[, factor(MODEL, levels=c('V','R'),labels=c('Model: Village','Model: Regional'))])
	set(sc, sc[, which(SC=="VILL_99_APR15")],'SC',"150701_VILL_SCENARIO-C")	
	set(sc, NULL, 'SC', sc[, factor(SC,	levels=c("150701_REGIONAL_TRAIN1", "150701_REGIONAL_TRAIN2", "150701_REGIONAL_TRAIN3", "150701_REGIONAL_TRAIN4","150701_REGIONAL_TRAIN5","150701_VILL_SCENARIO-A","150701_VILL_SCENARIO-B","150701_VILL_SCENARIO-C","150701_VILL_SCENARIO-D","150701_VILL_SCENARIO-E"), 
							labels=c('sc 1','sc 2','sc 3','sc 4','sc 5','sc A','sc B','sc C','sc D','sc E'))])
	set(sc, NULL, 'GAPS', sc[, factor(GAPS, levels=c('none','low','high'),labels=c('Gaps: none','Gaps: low','Gaps: high'))])
	set(sc, NULL, 'BEST', sc[, factor(BEST, levels=c('Y','N'),labels=c('best tree','replicate tree'))])									
	set(sc, NULL, 'GENE', sc[, factor(GENE, levels=c('POL','GAG+POL+ENV'),labels=c('pol','gag+pol+env'))])	
	set(sc, NULL, 'TEAM', sc[, factor(TEAM, levels=sc[, sort(unique(TEAM))],labels=sc[, sort(unique(TEAM))])])
	set(sc, NULL, 'EXT', sc[, factor(EXT, levels=c('~0pc','5pc'),labels=c('~ 0%/year','5%/year'))])
	set(sc, NULL, 'ART', sc[, factor(ART, levels=c('none','fast'),labels=c('none','fast'))])
	sc		<- subset(sc, OTHER=='N')
	stopifnot(sc[, !any(is.na(NRFC))], sc[, !any(is.na(NQDC))])
	scp		<- sc[, list( NRF=mean(NRFC, na.rm=TRUE), NQD=mean(NQDC, na.rm=TRUE) ), by=c('SC','GENE','TEAM','BEST','IDX','FILE','GAPS','MODEL','TAXAN','TAXAN_T','ROOTED','SEQCOV','ART','ACUTE','EXT','OTHER','TIME_IDX_T','SUB_IDX_T')]
	if('BILL'%in%colnames(sc))
	{
		tmp		<- sc[, list( BILL=mean(BILL, na.rm=TRUE) ), by=c('SC','GENE','TEAM','BEST','IDX','FILE','GAPS','MODEL','TAXAN','TAXAN_T','ROOTED','SEQCOV','ART','ACUTE','EXT','OTHER','TIME_IDX_T','SUB_IDX_T')]
		scp		<- merge(scp, tmp, by=c('SC','GENE','TEAM','BEST','IDX','FILE','GAPS','MODEL','TAXAN','TAXAN_T','ROOTED','SEQCOV','ART','ACUTE','EXT','OTHER','TIME_IDX_T','SUB_IDX_T') )				
	}	 
	sm		<- rbind( subset(sa, grepl('Village',MODEL)), scp, fill=TRUE, use.names=TRUE)
	#
	#	polvsall by gaps
	#	-->
	#	all leads to improvements throughout
	#	with gaps, the topology without branch lengths is increasingly difficult to estimate
	#	regional overall more difficult!
	if('NRF'%in%colnames(sm))
	{
		ggplot( subset(sm, TEAM!='MetaPIGA'), aes(y=NRF, x=SC, shape=TEAM, fill=GENE, colour=GENE, size=BEST) ) + 
				geom_jitter(position = position_jitter(height=.001, width=0.2)) +			
				scale_size_manual(values=c(3, 1)) +
				scale_shape_manual(values=c(21,23,24)) +
				scale_fill_brewer(palette='Paired') +
				scale_colour_brewer(palette='Paired') +
				facet_wrap(MODEL~GAPS, scales='free_x') +	
				labs(x='\nsimulated data set', y='Robinson-Fould\n(standardized)\n', size='', shape='Method', fill='part of genome', colour='part of genome') +
				theme_bw() 
		ggsave(w=10, h=6, file=paste(edir,'/',timetag,'_RF_polvsall_by_gaps.pdf',sep=''))		
	}
	if('BILL'%in%colnames(sm))
	{
		ggplot( subset(sm, TEAM!='MetaPIGA'), aes(y=BILL, x=SC, shape=TEAM, fill=GENE, colour=GENE, size=BEST) ) + 
				geom_jitter(position = position_jitter(height=.001, width=0.2)) +			
				scale_size_manual(values=c(3, 1)) +
				scale_shape_manual(values=c(21,23,24)) +
				scale_fill_brewer(palette='Paired') +
				scale_colour_brewer(palette='Paired') +
				facet_wrap(MODEL~GAPS, scales='free') +	
				labs(x='\nsimulated data set', y='Geodesic\n(raw)\n', size='', shape='Method', fill='part of genome', colour='part of genome') +
				theme_bw() 
		ggsave(w=10, h=6, file=paste(edir,'/',timetag,'_BL_polvsall_by_gaps.pdf',sep=''))		
	}	
	if('NQD'%in%colnames(sm))
	{
		ggplot( subset(sm, TEAM!='MetaPIGA'), aes(y=NQD, x=SC, shape=TEAM, fill=GENE, colour=GENE, size=BEST) ) + 
				geom_jitter(position = position_jitter(height=.001, width=0.2)) +			
				scale_size_manual(values=c(3, 1)) +
				scale_shape_manual(values=c(21,23,24)) +
				scale_fill_brewer(palette='Paired') +
				scale_colour_brewer(palette='Paired') +
				facet_wrap(MODEL~GAPS, scales='free_x') +	
				labs(x='\nsimulated data set', y='Quartett\n(standardized)\n', size='', shape='Method', fill='part of genome', colour='part of genome') +
				theme_bw() 
		ggsave(w=10, h=6, file=paste(edir,'/',timetag,'_QD_polvsall_by_gaps.pdf',sep=''))		
	}
	if('KC1'%in%colnames(sa))
	{		
		ggplot( subset(sa, TEAM!='MetaPIGA'), aes(y=KC1/TAXAN/TAXAN, x=SC) ) + 			
				geom_boxplot(aes(colour=GENE), fill='transparent', size=0.5, outlier.shape=NA, alpha=0.3) +
				geom_jitter(aes(shape=TEAM, fill=GENE, colour=GENE, size=BEST), position = position_jitter(height=.001, width=0.2)) +
				scale_size_manual(values=c(3, 1)) +
				scale_shape_manual(values=c(21,23,24)) +
				scale_fill_brewer(palette='Paired') +
				scale_colour_brewer(palette='Paired') +
				facet_wrap(MODEL~GAPS, scales='free_x') +	
				labs(x='\nsimulated data set', y='Kendall-Colijn\n(lambda=0, rtt if unrooted, /TX^2)\n', size='', shape='Method', fill='part of genome', colour='part of genome') +
				theme_bw() 
		ggsave(w=10, h=6, file=paste(edir,'/',timetag,'_KC1_polvsall_by_gaps.pdf',sep=''))
		ggplot( subset(sa, TEAM!='MetaPIGA'), aes(y=KC2/TAXAN/TAXAN, x=SC) ) + 			
				geom_boxplot(aes(colour=GENE), fill='transparent', size=0.5, outlier.shape=NA, alpha=0.3) +
				geom_jitter(aes(shape=TEAM, fill=GENE, colour=GENE, size=BEST), position = position_jitter(height=.001, width=0.2)) +
				scale_size_manual(values=c(3, 1)) +
				scale_shape_manual(values=c(21,23,24)) +
				scale_fill_brewer(palette='Paired') +
				scale_colour_brewer(palette='Paired') +
				facet_wrap(MODEL~GAPS, scales='free_x') +	
				labs(x='\nsimulated data set', y='Kendall-Colijn\n(lambda=0, rtt all, /TX^2)\n', size='', shape='Method', fill='part of genome', colour='part of genome') +
				theme_bw() 
		ggsave(w=10, h=6, file=paste(edir,'/',timetag,'_KC2_polvsall_by_gaps.pdf',sep=''))		
	}
	if('NPD'%in%colnames(sa))
	{
		ggplot( subset(sa, TEAM!='MetaPIGA'), aes(y=NPD, x=SC, shape=TEAM, fill=GENE, colour=GENE, size=BEST) ) + 
				geom_jitter(position = position_jitter(height=.001, width=0.2)) +			
				scale_size_manual(values=c(3, 1)) +
				scale_shape_manual(values=c(21,23,24)) +
				scale_fill_brewer(palette='Paired') +
				scale_colour_brewer(palette='Paired') +
				facet_wrap(MODEL~GAPS, scales='free') +	
				labs(x='\nsimulated data set', y='Path difference\n(standardized)\n', size='', shape='Method', fill='part of genome', colour='part of genome') +
				theme_bw() 
		ggsave(w=10, h=6, file=paste(edir,'/151023_PD_polvsall_by_gaps.pdf',sep=''))
	}
	#	RF may be confounded by size of data set when evaluating the extent that regional is more difficult
	#	-->
	#	hard to extrapolate how standardized RF grows with size of data set,
	#	but regression extrapolation suggests there is an effect
	if('NRF'%in%colnames(sm))
	{
		
		z		<- subset(sa, TEAM!='MetaPIGA' & !grepl('Reg',MODEL))
		mo		<- gamlss(NRF~TAXAN+GENE+GAPS, sigma.formula=~TAXAN+GENE+GAPS, family=BE(mu.link='cauchit'), data=z)
		tmp		<- subset(sa, TEAM!='MetaPIGA')
		tmp[, NRFP:=predict(mo, data=z, newdata=subset(tmp, select=c(TAXAN,GENE,GAPS)), what='mu',type='response')]
		ggplot( tmp, aes(x=TAXAN) ) + 
				geom_jitter(aes(y=NRF, shape=TEAM, colour=EXT, fill=EXT,  size=BEST), position = position_jitter(height=.001, width=20), alpha=0.7) +
				geom_line(aes(y=NRFP), colour='black', size=0.5) +
				scale_size_manual(values=c(3, 1)) +
				scale_shape_manual(values=c(21,23,24)) +
				scale_fill_brewer(palette='Set1') + scale_colour_brewer(palette='Set1') +
				scale_y_continuous(breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.1)) +
				facet_grid(GAPS~GENE) +
				labs(x='\nsize of simulated data set', y='Robinson-Fould\n(standardized)\n', size='', shape='Method', fill='trms/outside', colour='trms/outside') +
				theme_bw()
		ggsave(w=10, h=8, file=paste(edir,'/',timetag,'_RF_trmsoutside.pdf',sep=''))
	#
	#	
	#
		ggplot( subset(sa, MODEL=="Model: Village" & TEAM!='MetaPIGA'), aes(y=NRF, x=SC, shape=TEAM, fill=TAXAN, colour=TAXAN, size=BEST) ) + 
				geom_jitter(position = position_jitter(height=.001, width=0.2)) +			
				scale_size_manual(values=c(3, 1)) +
				scale_shape_manual(values=c(21,23,24)) +
				#scale_fill_brewer(palette='Paired') +
				#scale_colour_brewer(palette='Paired') +
				facet_wrap(GENE~GAPS, scales='free_x', ncol=3) +	
				labs(title="Model: Village\n", x='\nsimulated data set', y='Robinson-Fould\n(standardized)\n', size='', shape='Method', fill='Taxa in subm tree', colour='Taxa in subm tree') +
				theme_bw() 
		ggsave(w=10, h=6, file=paste(edir,'/',timetag,'_RF_polvsall_by_TAXAN_V.pdf',sep=''))
		ggplot( subset(sa, MODEL=="Model: Regional" & TEAM!='MetaPIGA'), aes(y=NRF, x=SC, shape=TEAM, fill=TAXAN, colour=TAXAN, size=BEST) ) + 
				geom_jitter(position = position_jitter(height=.001, width=0.2)) +			
				scale_size_manual(values=c(3, 1)) +
				scale_shape_manual(values=c(21,23,24)) +
				#scale_fill_brewer(palette='Paired') +
				#scale_colour_brewer(palette='Paired') +
				facet_wrap(GENE~GAPS, scales='free_x', ncol=3) +	
				labs(title="Model: Regional\n", x='\nsimulated data set', y='Robinson-Fould\n(standardized)\n', size='', shape='Method', fill='Taxa in subm tree', colour='Taxa in subm tree') +
				theme_bw() 
		ggsave(w=10, h=6, file=paste(edir,'/',timetag,'_RF_polvsall_by_TAXAN_R.pdf',sep=''))
	}
	
	
	#	teams
	#	--> 
	#	MetaPIGA fairly bad in terms of RF
	#	PhyML IQTree RAxML similar in terms of RF,
	#	but PhyML behaved poorly when many gaps
	if('NRF'%in%colnames(sm))
	{
		ggplot( sm, aes(y=NRF, x=SC, shape=TEAM, colour=TEAM, fill=TEAM, size=BEST) ) + 
				geom_jitter(position = position_jitter(height=.001, width=0.2), alpha=0.7) +
				scale_size_manual(values=c(4, 1)) +
				scale_shape_manual(values=c(21,22,23,24)) +
				scale_fill_brewer(palette='Set1') + scale_colour_brewer(palette='Set1') +
				scale_y_continuous(breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.1)) +
				facet_grid(MODEL~GENE) +			
				labs(x='\nsimulated data set', y='Robinson-Fould\n(standardized)\n', size='', shape='Method', fill='Method', colour='Method') +
				theme_bw() 
		ggsave(w=10, h=7, file=paste(edir,'/',timetag,'_RF_team_by_scenarioandgene.pdf',sep=''))
	}
	if('NQD'%in%colnames(sm))
	{
		ggplot( sm, aes(y=NQD, x=SC, shape=TEAM, colour=TEAM, fill=TEAM, size=BEST) ) + 
				geom_jitter(position = position_jitter(height = .001, width=0.2), alpha=0.7) +
				scale_size_manual(values=c(4, 1)) +
				scale_shape_manual(values=c(21,22,23,24)) +
				scale_fill_brewer(palette='Set1') + scale_colour_brewer(palette='Set1') +
				scale_y_continuous(breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.1)) +
				facet_grid(MODEL~GENE) +			
				labs(x='\nsimulated data set', y='Quartett\n(standardized)\n', size='', shape='Method', fill='Method', colour='Method') +
				theme_bw() 
		ggsave(w=10, h=7, file=paste(edir,'/',timetag,'_QD_team_by_scenarioandgene.pdf',sep=''))
	}
	if('BILL'%in%colnames(sm))
	{
		ggplot( sm, aes(y=BILL, x=SC, shape=TEAM, colour=TEAM, fill=TEAM, size=BEST) ) + 
				geom_jitter(position = position_jitter(height = .001, width=0.2), alpha=0.7) +
				scale_size_manual(values=c(4, 1)) +
				scale_shape_manual(values=c(21,22,23,24)) +
				scale_fill_brewer(palette='Set1') + scale_colour_brewer(palette='Set1') +
				#scale_y_continuous(breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.1)) +
				facet_grid(MODEL~GENE, scales='free') +			
				labs(x='\nsimulated data set', y='Geodesic\n(raw)\n', size='', shape='Method', fill='Method', colour='Method') +
				theme_bw() 
		ggsave(w=10, h=7, file=paste(edir,'/',timetag,'_BL_team_by_scenarioandgene.pdf',sep=''))
	}
	if('KC1'%in%colnames(sa))
	{
		ggplot( sa, aes(y=KC1/TAXAN/TAXAN, x=SC, shape=TEAM, colour=TEAM, fill=TEAM, size=BEST) ) + 
				geom_jitter(position = position_jitter(height=.001, width=0.2), alpha=0.7) +
				scale_size_manual(values=c(4, 1)) +
				scale_shape_manual(values=c(21,22,23,24)) +
				scale_fill_brewer(palette='Set1') + scale_colour_brewer(palette='Set1') +
				#scale_y_continuous(breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.1)) +
				facet_grid(MODEL~GENE) +			
				labs(x='\nsimulated data set', y='Kendall-Colijn\n(lambda=0, rtt if unrooted, /TX^2)\n', size='', shape='Method', fill='Method', colour='Method') +
				theme_bw() 
		ggsave(w=10, h=7, file=paste(edir,'/',timetag,'_KC1_team_by_scenarioandgene.pdf',sep=''))
		ggplot( sa, aes(y=KC2/TAXAN/TAXAN, x=SC, shape=TEAM, colour=TEAM, fill=TEAM, size=BEST) ) + 
				geom_jitter(position = position_jitter(height=.001, width=0.2), alpha=0.7) +
				scale_size_manual(values=c(4, 1)) +
				scale_shape_manual(values=c(21,22,23,24)) +
				scale_fill_brewer(palette='Set1') + scale_colour_brewer(palette='Set1') +
				#scale_y_continuous(breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.1)) +
				facet_grid(MODEL~GENE) +			
				labs(x='\nsimulated data set', y='Kendall-Colijn\n(lambda=0, rtt all, /TX^2)\n', size='', shape='Method', fill='Method', colour='Method') +
				theme_bw() 
		ggsave(w=10, h=7, file=paste(edir,'/',timetag,'_KC2_team_by_scenarioandgene.pdf',sep=''))
	}
		
	#	taxa excluded:	plot cluster RF as a function of cluster size
	#	-->
	#	excluding taxa did not lead to noticeably lower RFs
	if('NRFC'%in%colnames(sc))
	{
		tmp		<- sc[, list(NRF=median(NRFC, na.rm=TRUE)), by=c('SC','GENE','TEAM','BEST','IDX','FILE','GAPS','MODEL')]
		tmp		<- subset(tmp, TEAM!='MetaPIGA')	
		ggplot( tmp, aes(y=NRF, x=SC, shape=TEAM, fill=GENE, colour=GENE, size=BEST) ) + 
				geom_jitter(position = position_jitter(height = .001, width=0.2)) +			
				scale_size_manual(values=c(3, 1)) +
				scale_shape_manual(values=c(21,22,23,24)) +
				scale_fill_brewer(palette='Paired') +
				scale_colour_brewer(palette='Paired') +
				facet_wrap(MODEL~GAPS, scales='free_x') +	
				labs(x='\nsimulated data set', y='median Robinson-Fould\n(standardized)\n', size='', shape='Method', fill='part of genome', colour='part of genome') +
				theme_bw() 
		ggsave(w=10, h=6, file=paste(edir,'/',timetag,'_RFCLU_polvsall_by_gaps.pdf',sep=''))
	}	
	
	
	ggplot( subset(sc, GENE=='gag+pol+env'), aes(y=NRFC, x=TAXA_NC, size=BEST, shape=TEAM, fill=TAXA_NC<CLU_N, colour=TAXA_NC<CLU_N)) +
			geom_jitter(position = position_jitter(height=.001, width=0.1), alpha=0.7) +
			scale_size_manual(values=c(3, 1)) +
			scale_shape_manual(values=c(21,22,24), guide=FALSE) +
			scale_fill_brewer(palette='Set1') +
			scale_colour_brewer(palette='Set1') +
			scale_y_continuous(breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.1)) +			
			#coord_trans(x='log10') +
			scale_x_log10(breaks=c(1,2,3,4,5,6,8,10,20,50,100,200,300), minor_breaks=NULL) +
			labs(x='\nsize of sampled transmission cluster', y='Robinson-Fould\n(standardized per transmission cluster)\n', size='', shape='Method', fill='taxa excluded\nprior to reconstruction', colour='taxa excluded\nprior to reconstruction') +
			facet_grid(SC~TEAM) +
			theme_bw() 
	ggsave(w=16, h=8, file=paste(edir,'/',timetag,'_RFCLU_iftaxaexcludedbeforetreereconstruction.pdf',sep=''))
	
		
	#	effect of acute in terms of RF? --> Yes
	if('NRF'%in%colnames(sm))
	{
		ggplot( subset(sm, TEAM!='MetaPIGA' & TEAM!='PhyML' & grepl('Reg',MODEL) & !grepl('none',GAPS)), aes(y=NRF, x=ACUTE, shape=TEAM, fill=ACUTE, colour=ACUTE) ) + 
				geom_jitter(aes(size=BEST), position = position_jitter(height=.001, width=0.2), alpha=0.8) +
				geom_boxplot(outlier.shape=NA, colour='black', alpha=0.3) +
				scale_size_manual(values=c(3, 1)) +
				scale_shape_manual(values=c(21,24), guide=FALSE) +
				scale_fill_brewer(palette='Set1', guide=FALSE) +	scale_colour_brewer(palette='Set1', guide=FALSE) +
				facet_grid(GAPS~TEAM+GENE, scales='free_x') +	
				labs(x='\ntransmissions from those in acute infection', y='Robinson-Fould\n(standardized)\n', size='') +
				theme_bw() 
		ggsave(w=10, h=8, file=paste(edir,'/',timetag,'_RF_impactAcute.pdf',sep=''))
	}
	if('BILL'%in%colnames(sm))
	{
		ggplot( subset(sm, TEAM!='MetaPIGA' & TEAM!='PhyML' & grepl('Reg',MODEL) & !grepl('none',GAPS)), aes(y=BILL, x=ACUTE, shape=TEAM, fill=ACUTE, colour=ACUTE) ) + 
				geom_jitter(aes(size=BEST), position = position_jitter(height=.001, width=0.2), alpha=0.8) +
				geom_boxplot(outlier.shape=NA, colour='black', alpha=0.3) +
				scale_size_manual(values=c(3, 1)) +
				scale_shape_manual(values=c(21,24), guide=FALSE) +
				scale_fill_brewer(palette='Set1', guide=FALSE) +	scale_colour_brewer(palette='Set1', guide=FALSE) +
				facet_grid(GAPS~TEAM+GENE, scales='free_x') +	
				labs(x='\ntransmissions from those in acute infection', y='Geodesic\n(raw)\n', size='') +
				theme_bw() 
		ggsave(w=10, h=8, file=paste(edir,'/',timetag,'_BL_impactAcute.pdf',sep=''))
	}	
	if('NQD'%in%colnames(sm))
	{
		ggplot( subset(sm, TEAM!='MetaPIGA' & TEAM!='PhyML' & grepl('Reg',MODEL) & !grepl('none',GAPS)), aes(y=NQD, x=ACUTE, shape=TEAM, fill=ACUTE, colour=ACUTE) ) + 
				geom_jitter(aes(size=BEST), position = position_jitter(height=.001, width=0.2), alpha=0.8) +
				geom_boxplot(outlier.shape=NA, colour='black', alpha=0.3) +
				scale_size_manual(values=c(3, 1)) +
				scale_shape_manual(values=c(21,24), guide=FALSE) +
				scale_fill_brewer(palette='Set1', guide=FALSE) +	scale_colour_brewer(palette='Set1', guide=FALSE) +
				facet_grid(GAPS~TEAM+GENE, scales='free_x') +	
				labs(x='\ntransmissions from those in acute infection', y='Quartett\n(standardized)\n', size='') +
				theme_bw() 
		ggsave(w=10, h=8, file=paste(edir,'/',timetag,'_QD_impactAcute.pdf',sep=''))
	}
	if('KC1'%in%colnames(sa))
	{
		ggplot( subset(sa, TEAM!='MetaPIGA' & TEAM!='PhyML' & grepl('Reg',MODEL) & !grepl('none',GAPS)), aes(y=KC1/TAXAN/TAXAN, x=ACUTE, shape=TEAM, fill=ACUTE, colour=ACUTE) ) + 
				geom_jitter(aes(size=BEST), position = position_jitter(height=.001, width=0.2), alpha=0.8) +
				geom_boxplot(outlier.shape=NA, colour='black', alpha=0.3) +
				scale_size_manual(values=c(3, 1)) +
				scale_shape_manual(values=c(21,24), guide=FALSE) +
				scale_fill_brewer(palette='Set1', guide=FALSE) +	scale_colour_brewer(palette='Set1', guide=FALSE) +
				facet_grid(GAPS~TEAM+GENE, scales='free_x') +	
				labs(x='\ntransmissions from those in acute infection', y='Kendall-Colijn\n(lambda=0, rtt if unrooted, /TX^2)\n', size='') +
				theme_bw() 
		ggsave(w=10, h=8, file=paste(edir,'/',timetag,'_KC1_impactAcute.pdf',sep=''))
		ggplot( subset(sa, TEAM!='MetaPIGA' & TEAM!='PhyML' & grepl('Reg',MODEL) & !grepl('none',GAPS)), aes(y=KC2/TAXAN/TAXAN, x=ACUTE, shape=TEAM, fill=ACUTE, colour=ACUTE) ) + 
				geom_jitter(aes(size=BEST), position = position_jitter(height=.001, width=0.2), alpha=0.8) +
				geom_boxplot(outlier.shape=NA, colour='black', alpha=0.3) +
				scale_size_manual(values=c(3, 1)) +
				scale_shape_manual(values=c(21,24), guide=FALSE) +
				scale_fill_brewer(palette='Set1', guide=FALSE) +	scale_colour_brewer(palette='Set1', guide=FALSE) +
				facet_grid(GAPS~TEAM+GENE, scales='free_x') +	
				labs(x='\ntransmissions from those in acute infection', y='Kendall-Colijn\n(lambda=0, rtt all, /TX^2)\n', size='') +
				theme_bw() 
		ggsave(w=10, h=8, file=paste(edir,'/',timetag,'_KC2_impactAcute.pdf',sep=''))
	}
	
	#	effect of ART roll out in terms of RF? --> No
	if('NRF'%in%colnames(sa))
	{
		ggplot( subset(sa, TEAM!='MetaPIGA' & grepl('Vill',MODEL) & !grepl('none',GAPS)), aes(y=NRF, x=ART, shape=TEAM, fill=ART, colour=ART) ) + 
			geom_jitter(aes(size=BEST), position = position_jitter(height=.001, width=0.2), alpha=0.8) +
			geom_boxplot(outlier.shape=NA, colour='black', alpha=0.3) +
			scale_size_manual(values=c(3, 1)) +
			scale_shape_manual(values=c(21,23,24), guide=FALSE) +
			scale_fill_brewer(palette='Set2', guide=FALSE) +	scale_colour_brewer(palette='Set2', guide=FALSE) +
			facet_grid(GAPS~TEAM+GENE, scales='free_x') +	
			labs(x='\nART roll-out', y='Robinson-Fould\n(standardized)\n', size='') +
			theme_bw()
		ggsave(w=10, h=8, file=paste(edir,'/',timetag,'_RF_impactART.pdf',sep=''))
	}
	if('NQD'%in%colnames(sa))
	{
		ggplot( subset(sa, TEAM!='MetaPIGA' & grepl('Vill',MODEL) & !grepl('none',GAPS)), aes(y=NQD, x=ART, shape=TEAM, fill=ART, colour=ART) ) + 
				geom_jitter(aes(size=BEST), position = position_jitter(height=.001, width=0.2), alpha=0.8) +
				geom_boxplot(outlier.shape=NA, colour='black', alpha=0.3) +
				scale_size_manual(values=c(3, 1)) +
				scale_shape_manual(values=c(21,23,24), guide=FALSE) +
				scale_fill_brewer(palette='Set2', guide=FALSE) +	scale_colour_brewer(palette='Set2', guide=FALSE) +
				facet_grid(GAPS~TEAM+GENE, scales='free_x') +	
				labs(x='\nART roll-out', y='Quartett\n(standardized)\n', size='') +
				theme_bw()
		ggsave(w=10, h=8, file=paste(edir,'/',timetag,'_QD_impactART.pdf',sep=''))
	}
	if('BILL'%in%colnames(sm))
	{
		ggplot( subset(sm, TEAM!='MetaPIGA' & grepl('Vill',MODEL) & !grepl('none',GAPS)), aes(y=BILL, x=ART, shape=TEAM, fill=ART, colour=ART) ) + 
				geom_jitter(aes(size=BEST), position = position_jitter(height=.001, width=0.2), alpha=0.8) +
				geom_boxplot(outlier.shape=NA, colour='black', alpha=0.3) +
				scale_size_manual(values=c(3, 1)) +
				scale_shape_manual(values=c(21,23,24), guide=FALSE) +
				scale_fill_brewer(palette='Set2', guide=FALSE) +	scale_colour_brewer(palette='Set2', guide=FALSE) +
				facet_grid(GAPS~TEAM+GENE, scales='free_x') +	
				labs(x='\nART roll-out', y='Geodesic\n(raw)\n', size='') +
				theme_bw()
		ggsave(w=10, h=8, file=paste(edir,'/',timetag,'_BL_impactART.pdf',sep=''))
	}
	if('KC1'%in%colnames(sa))
	{
		ggplot( subset(sa, TEAM!='MetaPIGA' & grepl('Vill',MODEL) & !grepl('none',GAPS)), aes(y=KC1/TAXAN/TAXAN, x=ART, shape=TEAM, fill=ART, colour=ART) ) + 
				geom_jitter(aes(size=BEST), position = position_jitter(height=.001, width=0.2), alpha=0.8) +
				geom_boxplot(outlier.shape=NA, colour='black', alpha=0.3) +
				scale_size_manual(values=c(3, 1)) +
				scale_shape_manual(values=c(21,23,24), guide=FALSE) +
				scale_fill_brewer(palette='Set2', guide=FALSE) +	scale_colour_brewer(palette='Set2', guide=FALSE) +
				facet_grid(GAPS~TEAM+GENE, scales='free_x') +	
				labs(x='\nART roll-out', y='Kendall-Colijn\n(lambda=0, rtt if unrooted, /TX^2)\n', size='') +
				theme_bw()
		ggsave(w=10, h=8, file=paste(edir,'/',timetag,'_KC1_impactART.pdf',sep=''))
		ggplot( subset(sa, TEAM!='MetaPIGA' & grepl('Vill',MODEL) & !grepl('none',GAPS)), aes(y=KC2/TAXAN/TAXAN, x=ART, shape=TEAM, fill=ART, colour=ART) ) + 
				geom_jitter(aes(size=BEST), position = position_jitter(height=.001, width=0.2), alpha=0.8) +
				geom_boxplot(outlier.shape=NA, colour='black', alpha=0.3) +
				scale_size_manual(values=c(3, 1)) +
				scale_shape_manual(values=c(21,23,24), guide=FALSE) +
				scale_fill_brewer(palette='Set2', guide=FALSE) +	scale_colour_brewer(palette='Set2', guide=FALSE) +
				facet_grid(GAPS~TEAM+GENE, scales='free_x') +	
				labs(x='\nART roll-out', y='Kendall-Colijn\n(lambda=0, rtt all, /TX^2)\n', size='') +
				theme_bw()
		ggsave(w=10, h=8, file=paste(edir,'/',timetag,'_KC2_impactART.pdf',sep=''))
	}
}
treecomparison.submissions.150930<- function()	
{
	require(data.table)
	require(ape)
	
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/IQTree/IQTree201507'
	infiles	<- list.files(indir, pattern='treefile$', recursive=1, full.names=1)
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/PhyML'
	infiles	<- c(infiles, list.files(indir, pattern='*tree*', recursive=1, full.names=1))
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/RAxML'
	infiles	<- c(infiles, list.files(indir, pattern='*RAxML_bestTree*', recursive=1, full.names=1))	
	infiles	<- data.table(FILE=infiles)
	submitted.trees			<- lapply(infiles[, FILE], function(x)
			{
				cat(x)
				read.tree(file=x)	
			})
	names(submitted.trees)	<- infiles[, FILE]
	
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/MetaPIGA'
	tmp		<-  list.files(indir, pattern='*result*', recursive=1, full.names=1)
	tmp		<- data.table(FILE=tmp)
	
	tmp.trees			<- lapply(tmp[, FILE], function(x)
			{
				cat(x)
				read.nexus(file=x)	
			})
	sapply(tmp.trees, length)
	MetaPIGA.trees			<- c(lapply(tmp.trees, '[[', 1), lapply(tmp.trees, '[[', 2), lapply(tmp.trees, '[[', 3), lapply(tmp.trees, '[[', 4))
	names(MetaPIGA.trees)	<- c(sapply(tmp.trees, function(x) names(x)[1]), sapply(tmp.trees, function(x) names(x)[2]), sapply(tmp.trees, function(x) names(x)[3]), sapply(tmp.trees, function(x) names(x)[4]))	
	submitted.trees			<- c(submitted.trees, MetaPIGA.trees)	
	submitted.info			<- data.table(FILE=names(submitted.trees))
	
	submitted.info[, TEAM:=NA_character_]
	set(submitted.info, submitted.info[, which(grepl('RAXML',FILE))], 'TEAM', 'RAXML')
	set(submitted.info, submitted.info[, which(grepl('IQTree',FILE))], 'TEAM', 'IQTree')
	set(submitted.info, submitted.info[, which(grepl('MetaPIGA',FILE))], 'TEAM', 'MetaPIGA')
	set(submitted.info, submitted.info[, which(grepl('PhyML',FILE))], 'TEAM', 'PhyML')
	submitted.info[, SC:=NA_character_]
	tmp		<- submitted.info[, which(grepl('150701_Regional_TRAIN[0-9]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, regmatches(FILE, regexpr('150701_Regional_TRAIN[0-9]',FILE))])
	tmp		<- submitted.info[, which(grepl('150701_Vill_SCENARIO-[A-Z]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, regmatches(FILE, regexpr('150701_Vill_SCENARIO-[A-Z]',FILE))])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('TRAIN[0-9]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, paste('150701_Regional_',regmatches(FILE, regexpr('TRAIN[0-9]',FILE)),sep='')])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('scenario[A-Z]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, paste('150701_Vill_',regmatches(FILE, regexpr('scenario[A-Z]',FILE)),sep='')])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('150701_regional_train[0-9]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, regmatches(FILE, regexpr('150701_regional_train[0-9]',FILE))])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('150701_vill_scenario-[A-Z]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, regmatches(FILE, regexpr('150701_vill_scenario-[A-Z]',FILE))])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('Vill_99_Apr15', FILE))]
	set(submitted.info, tmp, 'SC', 'Vill_99_Apr15')
	
	set(submitted.info, NULL, 'SC', submitted.info[, toupper(SC)])
	tmp		<- submitted.info[, which(grepl('150701_VILL_SCENARIO[A-Z]', SC))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, gsub('150701_VILL_SCENARIO','150701_VILL_SCENARIO-',SC)])
	
	outfile	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/submitted_150911.rda'
	save(submitted.trees, submitted.info, file=outfile)
}