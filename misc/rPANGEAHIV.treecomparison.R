##--------------------------------------------------------------------------------------------------------
##	olli 23.10.15
##--------------------------------------------------------------------------------------------------------
#
#	compute path differences on complete trees
#
treedist.pathdifference.wrapper<- function(submitted.info, ttrs, strs)
{
	require(ape)
	require(phangorn)
	require(data.table)	
	#tmp			<- subset(submitted.info, IDX==463)[1,]
	#IDX<- 463
	#IDX_T<-7
	tmp				<- submitted.info[, {
				cat('\nAt IDX', IDX)
				#IDX<- 241; TIME_IDX_T<- 6; IDX<- 1; TIME_IDX_T<- 1 
				stree		<- strs_rtt[[IDX]]
				otree		<- ttrs[[TIME_IDX_T]]								
				z			<- setdiff(otree$tip.label, stree$tip.label)
				stopifnot( length(z)==abs(diff(c(TAXAN, TAXAN_T))) )
				if(length(z))
					otree	<- drop.tip(otree, z)				
				z			<- path.dist(otree, stree)
				list(PD=z, NPD=z/choose(Ntip(otree),2), NPDSQ=z/sqrt(choose(Ntip(otree),2)))
			}, by='IDX']
	tmp
}
#--------------------------------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------------------------------
treedist.pathdifference.clusters.wrapper<- function(submitted.info, ttrs, strs, tinfo)
{
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
								sclu	<- drop.tip(stree, setdiff(stree$tip.label,TAXA))
								oclu	<- drop.tip(otree, union( setdiff(otree$tip.label, stree$tip.label), setdiff(otree$tip.label,TAXA)))
								z		<- path.dist(oclu, sclu)
								list(PD=z, NPD=z/choose(Ntip(oclu),2), NPDSQ=z/sqrt(choose(Ntip(oclu),2)))								
							}, by='IDCLU']	
				}
				if(!nrow(z))
					ans		<- data.table(PD=NA_real_, NPD=NA_real_, NPDSQ=NA_real_)
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
	min.depth		<- 5
	
	infile			<- '/Users/Oliver/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/readlengths/bam_stats_150218.rda'
	load(infile)
	infile			<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/alignments_160110/PANGEA_HIV_n5003_Imperial_v160110_GlobalAlignment.rda'
	load(infile)	#loads sqi, sq
	outdir			<- '/Users/Oliver/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/data'
	
	
	
	#	add SANGER_ID
	infile.s	<- "~/Dropbox (Infectious Disease)/pangea_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v150908_Summary.csv"
	si			<- as.data.table(read.csv(infile.s, stringsAsFactors=FALSE))
	setnames(si, colnames(si), toupper(gsub('.','_',colnames(si),fixed=1))) 
	set(si, NULL, 'PANGEA_ID', si[, gsub(' ','',PANGEA_ID)])
	setnames(si, 'CLINICAL_GENOME_COVERAGE', 'COV')
	tmp			<- subset(si, select=c(PANGEA_ID, SANGER_ID))
	set(tmp, NULL, 'PANGEA_ID', tmp[,gsub('-','_',PANGEA_ID)])
	#	of the duplicate PANGEA_IDs, consider only those with larger coverage	
	sqi			<- merge(sqi, tmp, by='PANGEA_ID', all.x=1, allow.cartesian=TRUE)
	tmp			<- sqi[, list(SANGER_ID=SANGER_ID[which.max(COV)]), by='PANGEA_ID']
	sqi			<- merge(sqi, tmp, by=c('PANGEA_ID','SANGER_ID'))
	#	select min.coverage
	tmp			<- subset(sqi, !is.na(SITE) & COV>=min.coverage, c(PANGEA_ID, SANGER_ID, COV, SITE))
	tmp2		<- tmp[, which(is.na(SANGER_ID))]
	set(tmp, tmp2, 'SANGER_ID', tmp[tmp2, PANGEA_ID])
	setnames(bam.cov, c('FILE_ID','COV'), c('SANGER_ID','DEPTH'))
	bam.cov		<- merge(bam.cov, tmp, by='SANGER_ID')
	#	select min.depth
	bam.ch		<- subset(bam.cov, DEPTH>=min.depth)
	#	define chunks
	bam.ch[, POS_NEXT:= POS+REP]	
	bam.ch		<- bam.ch[, list(SITE=SITE, POS=POS, DEPTH=DEPTH, REP=REP, CHUNK=cumsum(as.numeric(c(TRUE, POS[-1]!=POS_NEXT[-length(POS_NEXT)])))), by='SANGER_ID']
	bam.ch		<- bam.ch[, list(SITE=SITE[1], POS_CH=min(POS), REP_CH=sum(REP), DEPTH_CH= sum(DEPTH*REP)/sum(REP) ), by=c('SANGER_ID','CHUNK')]
	bam.ch		<- merge(bam.ch, bam.ch[, list(COV=sum(REP_CH)), by='SANGER_ID'], by='SANGER_ID')
	bam.ch[, COVP:= COV/ncol(sq)]	
	bam.ch[, DEPTH_MIN:=min.depth]
	set(bam.ch, NULL, 'SITE', bam.ch[, factor(SITE, levels=c('BW', 'ZA', 'UG'), labels=c('Botswana', 'South Africa', 'Uganda'))])
	#
	unique(bam.ch)[, mean(COVP)]
	#	0.6400059
	unique(bam.ch)[, mean(COV)]
	#	6425.659
	
	#	plot chunks
	require(viridis)
	setkey(bam.ch, POS_CH, SITE)
	set(bam.ch, NULL, 'SANGER_ID', bam.ch[, factor(SANGER_ID, levels=unique(SANGER_ID), labels=unique(SANGER_ID))])
	ggplot(bam.ch, aes(x=SANGER_ID, xend=SANGER_ID, y=POS_CH, yend=POS_CH+REP_CH-1L, colour=SITE)) +
			scale_y_continuous(expand=c(0,0),limits=c(0,10e3), breaks=seq(0,10e3,500), minor_breaks=seq(0,10e3,100)) +
			scale_colour_manual(values=c('Botswana'="#1B0C42FF", 'South Africa'="#CF4446FF", 'Uganda'="#781C6DFF")) +			
			geom_segment() + theme_bw() + labs(y='genome position', x='Sequenced individual') + coord_flip() +
			theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank() )
	ggsave(file=file.path(outdir,'PANGEA_HIV_n5003_Imperial_v160110_selectedchunks.pdf'), w=10, h=60, limitsize = FALSE)
	#
	#	without bam files, just from alignment (min.depth is 10, but we can delete all-gap cols)
	#
	outdir			<- '/Users/Oliver/Dropbox (Infectious Disease)/2016_PANGEA_treecomp/figures'
	min.coverage	<- 600
	min.depth		<- 10
	#with.gaps		<- 1
	#outfile		<- '150623_PANGEAGlobal2681_C5_wgaps.pdf'
	with.gaps		<- 0
	outfile			<- '160110_PANGEAGlobal5003_C10.pdf'
	infile			<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/alignments_160110/PANGEA_HIV_n5003_Imperial_v160110_GlobalAlignment.rda'
	load(infile)	#loads sqi, sq
	#	sq has already LTR trimmed and sites outside ref compendium dropped
	#	drop references
	
	# 	drop gap only columns
	if(!with.gaps)
	{
		tmp			<- apply( as.character(sq), 2, function(x) !all(x%in%c('?','-','n')) ) 
		sq			<- sq[, tmp]		
	}
	#	convert into chunks
	ch				<- lapply(seq_len(nrow(sq)), function(i)
			{
				z	<- gregexpr('1+', paste(as.numeric( !as.character( sq[i,] )%in%c('-','?','n') ), collapse='') )[[1]]
				data.table(PANGEA_ID= rownames(sq)[i], POS=as.integer(z), DEPTH=min.depth, REP=attr(z,"match.length"))
			})
	ch				<- do.call('rbind',ch)	
	#	define SITE
	tmp				<- ch[, which(grepl('^R[0-9]+_',PANGEA_ID))]
	set(ch, tmp, 'SITE', 'ZA')
	tmp				<- ch[, which(is.na(SITE) & grepl('PG[0-9]+_[A-Z]+',PANGEA_ID))]	
	set(ch, tmp, 'SITE', ch[tmp, regmatches(PANGEA_ID, regexpr('PG[0-9]+_[A-Z]+',PANGEA_ID))])
	set(ch, NULL, 'SITE', ch[, gsub('PG[0-9]+_','',SITE)])
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
	setkey(ch, PANGEA_ID)	
	dcast.data.table(unique(ch)[, list(P=seq(0,1,0.1), Q=quantile(COV, p=seq(0,1,0.1))), by='SITE'], P~SITE, value.var='Q')
	dcast.data.table(unique(ch)[, list(Q=seq(1e3,8e3,1e3), P=ecdf(COV)(seq(1e3,8e3,1e3))), by='SITE'], Q~SITE, value.var='P')
	#	proportion of non-gaps in alignment
	#unique(ch)[, sum(COV)] / (nrow(unique(ch))*ncol(so))
	nrow(unique(ch))	
	#	C10: 	4253
	#	    Botswana 	South Africa       	Uganda 
    #     	344          813         		3096 
	unique(ch)[, mean(COVP)]
	#	C10:	0.5848233
	unique(ch)[, mean(COVP), by='SITE']
	#	1:     Botswana 0.7327631
	#	2:       Uganda 0.4991186
	#	3: South Africa 0.7605369
	unique(ch)[, mean(COV)]
	#	C10:	5702.612
	unique(ch)[, mean(COV), by='SITE']
	#	           SITE       V1
	#	1:     Botswana 7356.942
	#	2:       Uganda 5011.151
	#	3: South Africa 7635.791
	#	plot chunks
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
	ggplot(ch, aes(y=PANGEA_ID, yend=PANGEA_ID, x=POS_CH, xend=POS_CH+REP_CH-1L, colour=SITE)) +
			scale_x_continuous(expand=c(0,0), breaks=seq(0,10e3,1e3), minor_breaks=seq(0,10e3,100)) +
			scale_colour_manual(values=c('Botswana'="#1B0C42FF", 'South Africa'="#CF4446FF", 'Uganda'="#781C6DFF")) +			
			geom_segment() + theme_bw() + 
			facet_wrap(~PLOT, scales='free_y', ncol=4) +
			labs(x='\nalignment position', y='PANGEA-HIV sequences\n', colour='sampling\nlocation') + 			
			theme(	axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank(), legend.position='bottom',
					strip.text= element_blank(), strip.background=element_blank()) +
			guides(colour=guide_legend(override.aes=list(size=5)))	
	ggsave(file=file.path(outdir,outfile), w=15, h=15, limitsize = FALSE)
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
				stopifnot( length(z)==abs(diff(c(TAXAN, TAXAN_T))) )
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
					ans		<- list(TR_REC_perc_T= ans[, mean(as.numeric(TRUE_PAIR))]  )	
				}
				ans
			}, by='IDX_T']
	setnames(tmp, 'IDX_T','SUB_IDX_T')
	tmp
}
##--------------------------------------------------------------------------------------------------------
##	olli 13.07.16
##--------------------------------------------------------------------------------------------------------
treedist.closest.ind.reconstructed<- function(submitted.info, tinfo, gd.thresh)
{
	sucl			<- subset(submitted.info, MODEL=='R')[, {
				print(IDX)
				#IDX<- 557; SUB_IDX_T<-2; SC<- '150701_REGIONAL_TRAIN2'
				ph			<- strs[[IDX]]
				model.reg	<- grepl('REGIONAL',SC)
				gds			<- treedist.closest.ind(ph, model.reg)
				gds			<- subset(gds, IDPOP_CL_GD<=gd.thresh)
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
					ans		<- ans[, list(IN=any(IN)), by='CL_PH_PAIR']
					ans		<- ans[, mean(as.numeric(IN))] 
				}
				if(!nrow(gds))
					ans		<- NA_real_				
				list(TR_REC_perc= ans  )				
			}, by=c('IDX')]
	sucl
}
##--------------------------------------------------------------------------------------------------------
##	olli 13.07.16
##--------------------------------------------------------------------------------------------------------
treedist.closest.ind.reconstructed.oftruepairs<- function(submitted.info, tinfo.pairs)
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
	require(phangorn)
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
				stopifnot( length(z)==abs(diff(c(TAXAN, TAXAN_T))) )
				if(length(z))
					otree	<- unroot(drop.tip(otree, z))				
				#https://groups.google.com/forum/#!topic/raxml/JgvxgknTeqw
				#normalize with 2n-6		
				rf			<- RF.dist(otree, stree, check.labels=TRUE)
				list(RF=rf, NRF=rf/(2*Ntip(otree)-6))
			}, by='IDX']
	tmp
}
##--------------------------------------------------------------------------------------------------------
##	olli 03.12.15
##--------------------------------------------------------------------------------------------------------
treedist.robinsonfouldclusters.wrapper<- function(submitted.info, ttrs, strs, tinfo)
{
	require(phangorn)
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
				stopifnot( length(z)==abs(diff(c(TAXAN, TAXAN_T))) )
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
treecomparison.explaingaps.evaluate.160725<- function()
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
	tmp		<- unique(subset(dpand, !is.na(PANGEA_ID), select=c(PANGEA_ID, TAXA, RCCS_studyid, REGA_GAG_A, REGA_GAG_AS, REGA_GAG_PURE, REGA_GAG_PURES)))
	chr		<- merge(ch, tmp, by='TAXA', all.x=1)
	set(chr, chr[, which(is.na(REGA_GAG_A))], 'REGA_GAG_A', 'No matched ID')	
	chr		<- subset(chr, !is.na(RCCS_studyid))
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
	tmp		<- unique(subset(dpand, !is.na(PANGEA_ID), select=c(PANGEA_ID, TAXA, RCCS_studyid, COMET_Region1, COMET_Region2, COMET_Region3)))
	chr		<- merge(ch, tmp, by='TAXA', all.x=1)
	chr		<- subset(chr, !is.na(RCCS_studyid))
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
	#	plot Rakai samples by ART
	#
	#	redefine ordering	
	tmp		<- unique(subset(dpand, !is.na(PANGEA_ID), select=c(PANGEA_ID, TAXA, RCCS_studyid, date, arvStartDate, selfReportArt, everSelfReportArt, recentVL)))
	tmp[, ART:= as.integer(arvStartDate<date)]
	set(tmp, tmp[, which(is.na(ART))], 'ART', 0L)
	set(tmp, tmp[, which(ART==0 & everSelfReportArt==1)], 'ART', 2L)
	set(tmp, tmp[, which(ART==0 & recentVL<1e4)], 'ART', 3L)
	set(tmp, NULL, 'ART', tmp[, factor(ART, levels=c(0L,1L,2L,3L), labels=c('no ART', 'ART started', 'ART self reported','no ART but VL<1e4'))])
	
	chr		<- merge(ch, tmp, by='TAXA', all.x=1)
	set(chr, chr[, which(is.na(ART))], 'ART', 'No matched ID')	
	chr		<- subset(chr, !is.na(RCCS_studyid))
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
	tmp		<- unique(subset(dpand, !is.na(PANGEA_ID), select=c(PANGEA_ID, TAXA, RCCS_studyid, date, arvStartDate, selfReportArt, everSelfReportArt, recentVL)))
	tmp[, VL:= cut(recentVL, breaks=c(0, 1e4, 2e4, 4e4, 1e5, Inf), labels=c('<1e4','1e4-2e4','2e4-4e4','4e4-1e5','>1e5'))]	
	set(tmp, tmp[, which(is.na(VL))], 'VL', 'No VL measured')
	chr		<- merge(ch, tmp, by='TAXA', all.x=1)
	set(chr, chr[, which(is.na(VL))], 'VL', 'No matched ID')	
	chr		<- subset(chr, !is.na(RCCS_studyid))	
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
	tmp		<- unique(subset(dpand, !is.na(PANGEA_ID), select=c(PANGEA_ID, TAXA, RCCS_studyid, REGION)))
	set(tmp, NULL, 'REGION', tmp[, factor(REGION)])
	chr		<- merge(ch, tmp, by='TAXA', all.x=1)
	chr		<- subset(chr, !is.na(RCCS_studyid))	
	chr[, PLOT:=NULL]	
	setkey(chr, TAXA)	
	tmp		<- unique(chr)
	setkey(tmp, REGION, COVP, TAXA)	
	tmp		<- tmp[, list(TAXA=TAXA, PLOT=REGION, PLOT_ID=seq_along(TAXA)), by='REGION']
	chr		<- merge(chr, subset(tmp, select=c(TAXA, PLOT_ID)), by='TAXA')	
	ggplot(chr) +
			geom_segment(aes(y=PLOT_ID, yend=PLOT_ID, x=POS_CH, xend=POS_CH+REP_CH-1L, colour=REGION)) +  
			geom_rect(data=dpani, aes(xmin=START, xmax=END, ymin=-Inf, ymax=Inf), fill="black") +			
			facet_wrap(~REGION, scales='free_y', ncol=6) +
			scale_x_continuous(expand=c(0,0), breaks=dpani$START, labels=dpani$PR) +
			scale_y_continuous(expand=c(0,0)) +
			scale_colour_brewer(palette='Dark2') +						
			labs(x='\nalignment position', y='Rakai PANGEA-HIV sequences\n', colour='region') +
			theme_bw() +
			theme(	legend.position='bottom', strip.text= element_blank(), strip.background=element_blank()) +
			guides(colour=guide_legend(override.aes=list(size=5)))
	ggsave(file=file.path(wdir,gsub('.rda','_gapsprimers_REGION.pdf',wfile)), h=7, w=20, limitsize = FALSE)
	#
	#	select taxa that failed after 1R despite high viral load + not on ART
	#	select taxa that were successful after 1R despite low viral load + not on ART
	#
	set(dpand, NULL, 'date', dpand[,hivc.db.Date2numeric(date)])	
	set(dpand, NULL, 'firstPosDate', dpand[,hivc.db.Date2numeric(firstPosDate)])
	set(dpand, NULL, 'arvStartDate', dpand[,hivc.db.Date2numeric(arvStartDate)])
	set(dpand, NULL, 'FirstSelfReportArt', dpand[,hivc.db.Date2numeric(FirstSelfReportArt)])
	set(dpand, NULL, 'recentVLdate', dpand[,hivc.db.Date2numeric(recentVLdate)])
	#	viral load and coverage
	dsmut	<- subset(dpand,	PR=='1R' & POS=='PR_1' &	
								UNASS_TO_NEXTPRIMER_P<0.2 &
								recentVL<2e4 & abs(date-recentVLdate)<.5, select=c(TAXA, UNASS_TO_NEXTPRIMER_P))
	dsmut[, TYPE:= 'Rakai_VL_low_Coverage_high']		
	tmp		<- subset(dpand,	PR=='1R' & 	POS=='PR_1' &
					UNASS_TO_NEXTPRIMER_P>0.95 &
					(is.na(arvStartDate) | date<arvStartDate) &
					(!everSelfReportArt | everSelfReportArt & date<FirstSelfReportArt) &
					recentVL>4e4 & abs(date-recentVLdate)<.5, select=c(TAXA, UNASS_TO_NEXTPRIMER_P))
	tmp[, TYPE:= 'Rakai_VL_high_Coverage_none']			
	dsmut	<- rbind(dsmut, tmp)
	#	coverage
	tmp		<- subset(dpand,	PR=='1R' & POS=='PR_1' & !is.na(RCCS_studyid) &
					UNASS_TO_NEXTPRIMER_P<0.2, select=c(TAXA, UNASS_TO_NEXTPRIMER_P))
	tmp[, TYPE:= 'Rakai_Coverage_high']			
	dsmut	<- rbind(dsmut, tmp)
	tmp		<- subset(dpand,	PR=='1R' & 	POS=='PR_1' &
					UNASS_TO_NEXTPRIMER_P>0.95, select=c(TAXA, UNASS_TO_NEXTPRIMER_P))
	tmp[, TYPE:= 'Rakai_Coverage_none']			
	dsmut	<- rbind(dsmut, tmp)
	#	population
	tmp		<- subset(dpand,	PR=='1R' & 	POS=='PR_1' & !is.na(PANGEA_ID), select=c(TAXA, UNASS_TO_NEXTPRIMER_P))
	tmp[, TYPE:= 'PANGEA_All']
	dsmut	<- rbind(dsmut, tmp)
	tmp		<- subset(dpand,	PR=='1R' & 	POS=='PR_1' & !is.na(RCCS_studyid), select=c(TAXA, UNASS_TO_NEXTPRIMER_P))
	tmp[, TYPE:= 'Rakai_All']
	dsmut	<- rbind(dsmut, tmp)
	#	subtypes
	tmp		<- subset(dpand,	PR=='1R' & POS=='PR_1' & COMET_Region1=='A1' &
					UNASS_TO_NEXTPRIMER_P>0.95, select=c(TAXA, UNASS_TO_NEXTPRIMER_P))
	tmp[, TYPE:= 'Rakai_A1_Coverage_none']			
	dsmut	<- rbind(dsmut, tmp)
	tmp		<- subset(dpand,	PR=='1R' & POS=='PR_1' & COMET_Region1=='C' &
					UNASS_TO_NEXTPRIMER_P>0.95, select=c(TAXA, UNASS_TO_NEXTPRIMER_P))
	tmp[, TYPE:= 'Rakai_C_Coverage_none']			
	dsmut	<- rbind(dsmut, tmp)
	tmp		<- subset(dpand,	PR=='1R' & POS=='PR_1' & COMET_Region1=='D' &
					UNASS_TO_NEXTPRIMER_P>0.95, select=c(TAXA, UNASS_TO_NEXTPRIMER_P))
	tmp[, TYPE:= 'Rakai_D_Coverage_none']			
	dsmut	<- rbind(dsmut, tmp)
	
	dsmut	<- merge(subset(dpand, PR=='2R' | PR=='2F', select=c(TAXA, PR, POS, NT_DIFF, REGION, COMM_NUM, HH_NUM, SEX, AGEYRS, COMET_Region1)), dsmut, by=c('TAXA'), allow.cartesian=TRUE)	
	set(dsmut, NULL, 'POS', dsmut[,gsub('PR_','',POS)])
	set(dsmut, NULL, 'PR', dsmut[, paste('PR_',PR, sep='')])
	setnames(dsmut, 'UNASS_TO_NEXTPRIMER_P', 'PR1R_UNASS_TO_NEXTPRIMER_P')
	#dcast.data.table(dsmut, TYPE+TAXA+PR1R_UNASS_TO_NEXTPRIMER_P+REGION+COMM_NUM+HH_NUM+SEX+AGEYRS+COMET_Region1  ~  PR+POS, value.var='NT_DIFF')
	
	tmp		<- dsmut[,{
				z	<- round(as.numeric(binconf(length(which(NT_DIFF==1)), length(which(!is.na(NT_DIFF))))), d=3)
				list(EST=c('central','l95','u95'), VAL= z)				
			}, by=c('TYPE','PR','POS')]
	tmp		<- dcast.data.table(tmp, PR+POS+TYPE~EST, value.var='VAL')
	set(tmp, NULL, 'POS', tmp[, as.integer(POS)])
	set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, levels=c("PANGEA_All", "Rakai_All", "Rakai_Coverage_high", "Rakai_Coverage_none", "Rakai_VL_high_Coverage_none", "Rakai_VL_low_Coverage_high", 'Rakai_A1_Coverage_none', 'Rakai_C_Coverage_none', 'Rakai_D_Coverage_none'))])
	
	
	ggplot(subset(tmp, TYPE%in%c('PANGEA_All','Rakai_Coverage_high','Rakai_A1_Coverage_none','Rakai_C_Coverage_none','Rakai_D_Coverage_none')), aes(x=POS, fill=TYPE)) + 
			geom_bar(aes(y=central), stat='identity', width=0.7, position=position_dodge(0.8)) + 
			facet_grid(PR~.) + 
			geom_linerange(aes(ymin= l95, ymax=u95), position=position_dodge(0.8)) +
			theme_bw() + theme(legend.position='bottom') +
			scale_x_continuous(breaks=tmp[, seq_len(max(POS))]) +
			scale_y_continuous(labels=percent, expand=c(0,0), limits=c(0,1)) +
			labs(x='\nNucleotide position in primer', y='PANGEA sequences with mutation from HXB2\n', fill='selected sequences') 
	ggsave(file=file.path(wdir,gsub('.rda','_gapsprimers_2F2R_eval1.pdf',wfile)), h=10, w=15, limitsize = FALSE)
	
	
	ggplot(subset(tmp, TYPE%in%c('PANGEA_All','Rakai_Coverage_none','Rakai_VL_high_Coverage_none','Rakai_VL_low_Coverage_high')), aes(x=POS, fill=TYPE)) + 
			geom_bar(aes(y=central), stat='identity', width=0.7, position=position_dodge(0.8)) + 
			facet_grid(PR~.) + 
			geom_linerange(aes(ymin= l95, ymax=u95), position=position_dodge(0.8)) +
			theme_bw() + theme(legend.position='bottom') +
			scale_x_continuous(breaks=tmp[, seq_len(max(POS))]) +
			scale_y_continuous(labels=percent, expand=c(0,0), limits=c(0,1)) +
			labs(x='\nNucleotide position in primer', y='PANGEA sequences with mutation from HXB2\n', fill='selected sequences') 
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
treecomparison.explaingaps.collect.data<- function()
{
	require(ape)
	require(data.table)
	wdir			<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/explaingaps'
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
		#	save
		write.dna(sq, file=file.path(wdir,'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.fasta'), format='fa')
		wfile			<- 'PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.rda'
		save(sq, sqi, file=file.path(wdir,wfile))
	}	
	#
	#	find primer coordinates
	#	
	if(0)
	{
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
		tmp		<- dpan[, list(START=unlist(gregexpr(paste(SEQ, collapse=''),hxb2))), by='PR']
		stopifnot( tmp[, all(START>1)] )
		#	find primers in global alignment
		sqhxb2i	<- which(grepl('HXB2',rownames(sq)))
		sqhxb2	<- paste(as.character(sq[sqhxb2i,]), collapse='')
		tmp		<- dpan[, list(START=unlist(gregexpr(paste(SEQ, collapse='-*'),sqhxb2))), by='PR']
		#	1F and 4R are not part of PANGEA alignment
		tmp		<- subset(tmp, START>1)
		dpan	<- merge(dpan, tmp, by='PR', all.x=1)
		save(sq, sqi, dpan, file=file.path(wdir,wfile))
	}
	#
	#	evaluate primer alignments and add meta data (rccs VL, AGE, SUBTYPE)
	#
	if(0)
	{
		load(file.path(wdir,wfile))
		#	extract primer alignments and write to file
		subset(dpan, !is.na(START))[, {
					write.dna(sq[, seq.int(START[1], len=length(START))], file=file.path(wdir, paste(gsub('.rda','',wfile),'_primer_',PR[1],'_start_',START[1],'.fasta', sep='')), format='fa')
				}, by='PR']	
		#	get differences relative to HXB2 on primer by primer position (IDX)
		dpand		<- subset(dpan, !is.na(START))[, {
					#START	<- rep(1936, 23)
					psq		<- as.character( sq[, seq.int(START[1], len=length(START))] )
					sqhxb2i	<- which(grepl('HXB2',rownames(psq)))				
					tmp		<- as.data.table( melt( t( t(psq)==psq[sqhxb2i,] ) ) )
					setnames(tmp, c('Var1','Var2','value'), c('TAXA','POS','NT_DIFF'))
					tmp2	<- as.data.table( melt( psq=='?' ) )
					setnames(tmp2, c('Var1','Var2','value'), c('TAXA','POS','MISS'))
					tmp		<- merge(tmp, tmp2, by=c('TAXA','POS'))				
					set(tmp, NULL, 'NT_DIFF', tmp[, as.integer(!NT_DIFF)])
					set(tmp, NULL, 'POS', tmp[, paste('PR_',POS,sep='')])
					set(tmp, NULL, 'POS', tmp[, factor(POS, levels=paste('PR_',1:tmp[, length(unique(POS))],sep=''))])
					set(tmp, tmp[, which(MISS)], 'NT_DIFF', NA_integer_)
					tmp[, MISS:=NULL]
					tmp
				}, by='PR' ]
		#	calculate number of '?' until next primer start and for the next 250 sites
		tmp			<- rbind( 	data.table(PR='1F', END_B4NXT=subset(dpan, PR=='2F')[1, START-1L]),
								data.table(PR='2F', END_B4NXT=subset(dpan, PR=='1R')[ , START[1]+max(IDX)]),								
								data.table(PR='1R', END_B4NXT=subset(dpan, PR=='3F')[1, START-1L]),
								data.table(PR='3F', END_B4NXT=subset(dpan, PR=='2R')[, START[1]+max(IDX)]),
								data.table(PR='2R', END_B4NXT=subset(dpan, PR=='4F')[1, START[1]-1L]),
								data.table(PR='4F', END_B4NXT=subset(dpan, PR=='3R')[, START[1]+max(IDX)])		)	
		dpan		<- merge(dpan, tmp, by='PR', all.x=1)
		dun			<- subset(dpan, !is.na(START) & !is.na(END_B4NXT) & IDX==1)[, {
					#START<- 812; END_B4NXT<- 4341
					tmp			<- as.character( sq[, seq.int(START, END_B4NXT)] )
					tmp			<- apply(tmp=='?',1,sum)
					if(START<END_B4NXT)
					{
						tmp2	<- as.character( sq[, seq.int(START-20, START-1)] )
						tmp2	<- apply(tmp2=='?',1,sum)
					}
					if(START>END_B4NXT)
					{
						tmp2	<- as.character( sq[, seq.int(START+1, START+20)] )
						tmp2	<- apply(tmp2=='?',1,sum)
					}	
					if(START<END_B4NXT)
					{
						tmp3	<- as.character( sq[, seq.int(START, START+250)] )
						tmp3	<- apply(tmp3=='?',1,sum)
					}
					if(START>END_B4NXT)
					{
						tmp3	<- as.character( sq[, seq.int(START-250, START)] )
						tmp3	<- apply(tmp3=='?',1,sum)
					}	
					list(TAXA= names(tmp), UNASS_TO_NEXTPRIMER_P= tmp/length(seq.int(START, END_B4NXT)), UNASS_INDIR_250_P= tmp3/250, DROPPED_B4_P=tmp2/20  )
				}, by='PR']
		dpand		<- merge(dpand, dun, by=c('PR','TAXA'), all.x=1)
		dpand[, PANGEA_ID:= NA_character_]
		tmp			<- dpand[, which(grepl('^PG[0-9]+',TAXA))]
		set(dpand, tmp, 'PANGEA_ID', dpand[tmp, regmatches(TAXA,regexpr('PG[0-9]+-[A-Za-z0-9]+', TAXA))])	
		#	add RCCS data
		load("~/Dropbox (Infectious Disease)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/RakaiPangeaMetaData.rda")
		rccsData			<- as.data.table(rccsData)
		subtypeSummaryData	<- as.data.table(subtypeSummaryData)
		rccsData	<- subset(rccsData, select=c('Pangea.id', 'RCCS_studyid', 'date', 'batch', 'birthyr', 'REGION', 'COMM_NUM', 'HH_NUM', 'SEX', 'AGEYRS','firstPosDate', 'arvStartDate', 'selfReportArt', 'everSelfReportArt', 'FirstSelfReportArt', 'recentVL', 'recentVLdate'))
		setnames(rccsData, 'Pangea.id', 'PANGEA_ID')
		rccsData	<- subset(rccsData, !is.na(PANGEA_ID))
		setkey(rccsData, PANGEA_ID)
		rccsData	<- unique(rccsData)	#remove 4 duplicates "K104085" "E106462" "F030186" "F101874"	
		dpand		<- merge(dpand, rccsData, by='PANGEA_ID', all.x=1)
		#	add CLASS subtype data
		dst			<- subset(subtypeSummaryData, select=c(RCCS_studyid, gp.class, gag.class, pol.class, vpu.class, env.class, comp.class))
		setkey(dst, RCCS_studyid)
		dst			<- unique(dst)	#no duplicates here
		dst			<- melt(dst, id.var='RCCS_studyid')
		set(dst, NULL, 'value', dst[, gsub('\\s','',gsub('Complex','',gsub('Subtype', '',value)))])
		set(dst, NULL, 'variable', dst[, paste('SUBTYPE_',toupper(gsub('\\.class','',variable)),sep='')])	
		dst			<- dcast.data.table(dst, RCCS_studyid~variable)
		dpand		<- merge(dpand, dst, by='RCCS_studyid',all.x=1)
		#	add REGA subtype data
		load("~/Dropbox (Infectious Disease)/PANGEA_alignments/Rega Subtype Analysis/Gag REGA results/regaGag.rda")
		regaGag		<- as.data.table(regaGag)
		regaGag		<- subset(regaGag, select=c(TAXA, assignment, support, pure, pure_support))
		set(regaGag, NULL, 'assignment', regaGag[, gsub('Check the Report|Check the report|Check the bootscan','Check bootscan', gsub('Subtype ','', gsub('HIV-1 ','', assignment)))])
		setnames(regaGag, c('assignment','support','pure','pure_support'),c('REGA_GAG_A','REGA_GAG_AS','REGA_GAG_PURE','REGA_GAG_PURES'))
		dpand		<- merge(dpand, regaGag, by='TAXA',all.x=1)
		#	add COMET subtype data
		indir		<- '~/Dropbox (Infectious Disease)/PANGEA_alignments/Rakai Data for IqTree'
		infiles		<- data.table(FILE=list.files(indir, full.names=1, pattern='.rda$', recursive=1))	
		infiles[, GENE:= regmatches(FILE,regexpr('Region[0-9]',FILE))]
		dc			<- infiles[, {
					load(FILE)
					dc			<- as.data.table(summaryData)
					setnames(dc, 'seqid', 'TAXA')
					tmp			<- dc[, which(grepl('PG[0-9]+-[A-Z]+', TAXA))]
					dc[, SITE:= NA_character_]
					set(dc, tmp, 'SITE', dc[tmp,regmatches(TAXA,regexpr('PG[0-9]+-[A-Z]+', TAXA))])					
					dc			<- subset(dc, dataSource=='Pangea')
					set(dc, NULL, 'TAXA', dc[, regmatches(TAXA,regexpr('PG[0-9]+.*', TAXA))])
					dc			<- subset(dc, select=c(TAXA, CometSubtype))
					set(dc, NULL, 'CometSubtype', dc[, as.character(CometSubtype)])
					dc
				}, by='GENE']
		set(dc, dc[, which(grepl('unassigned', CometSubtype))], 'CometSubtype', 'unassigned')
		set(dc, dc[, which(grepl('check', CometSubtype))], 'CometSubtype', 'check')
		set(dc, dc[, which(!CometSubtype%in%c('A1','A2','B','C','D','unassigned','check'))], 'CometSubtype', 'other')
		set(dc, NULL, 'GENE', dc[, paste('COMET_',GENE,sep='')])
		dc			<- dcast.data.table(dc, TAXA~GENE, value.var='CometSubtype')
		dpand		<- merge(dpand, dc, by='TAXA',all.x=1)		
		#	save
		save(sq, sqi, dpan, dpand, file=file.path(wdir,wfile))
	}
	#
	#
	#
	dcast.data.table(tmp, TAXA~POS, value.var='NT_DIFF')
	subset(dpand, PR=='1R')
	dcast.data.table(subset(dpand, PR=='1R'), TAXA~POS, value.var='NT_DIFF')
}

##--------------------------------------------------------------------------------------------------------
##	olli 27.06.11
treecomparison.bootstrap.sd.vs.coverage<- function(indir=NULL, wdir=NULL)
{
	require(ape)
	require(data.table)
	require(ggplot2)
	bsn		<- 1e3
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
	seqi	<- subset(seqi, select=c(GENE,START,END))
	tp		<- as.data.table( t(combn(rownames(seq),2)) )
	setnames(tp, c('V1','V2'), c('TAXA1','TAXA2'))
	tp[, IDX:= seq_len(nrow(tp))]
	#	subset by batch
	tp[, BATCH:= ceiling(IDX/batch.n)]
	#
	#	calculate overlap in genomic coverage for each taxon pair
	#
	seqi	<- rbind(seqi, data.table(GENE='full', START=1L, END= seqi[, max(END)]))
	tmp		<- tp[, 	{
				#cat('IDX',IDX, round(IDX/nrow(tp),d=3))
				#TAXA1	<- 'IDPOP_13649|M|DOB_1906.66|2011.23'; TAXA2<- 'IDPOP_27993|F|DOB_1961.29|1991.587'
				#START<- 1; END<- 1473				
				df.gd	<- seqi[, {
							spc		<- as.character(seq[c(TAXA1,TAXA2), START:END])
							do		<- sum(apply( spc!='?', 2, prod))
							list(DO=do, GENE_LEN=END-START+1L)								
						}, by='GENE']
				list(GENE= df.gd$GENE, DO= df.gd$DO, GENE_LEN= df.gd$GENE_LEN)
			}, by=c('TAXA1','TAXA2')]
	tp		<- merge(tp, tmp, by=c('TAXA1','TAXA2'))
	#	save
	save(tp, seqi, seq, file=file.path(wdir, gsub('.fa','_save160726.rda',infile)))
	
	stop()
	#
	#	read genetic distances between taxon pairs
	#	
	infiles	<- data.table(FILE=list.files(wdir, pattern='BATCH[0-9]+.rda$', full.names=TRUE))
	infiles[, BATCH:= as.integer(gsub('BATCH','',regmatches(FILE, regexpr('BATCH[0-9]+', FILE))))]
	setkey(infiles, BATCH)
	tmp		<- infiles[, {
				load(FILE)
				tpi
			}, by='FILE']
	tmp[, FILE:=NULL]
	tp		<- merge(tp, tmp, by=c('TAXA1','TAXA2'))
	
	#	dev code
	setkey(tp, TAXA1, TAXA2)
	tmp		<- unique(tp)[, 	{
				#TAXA1	<- 'IDPOP_13649|M|DOB_1906.66|2011.23'; TAXA2<- 'IDPOP_27993|F|DOB_1961.29|1991.587'
				#START<- 1; END<- 1473				
				df.gd	<- seqi[, {
							spc		<- as.character(seq[c(TAXA1,TAXA2), START:END])
							do		<- sum(apply( spc!='?', 2, prod))
							list(DO=do, GENE_LEN=END-START+1L)								
						}, by='GENE']
				list(GENE= df.gd$GENE, DO= df.gd$DO, GENE_LEN= df.gd$GENE_LEN)
			}, by=c('TAXA1','TAXA2')]
	set(tmp, tmp[, which(GENE=='full')], 'GENE', 'gag+pol+env')
	tp		<- merge(tp, tmp, by=c('TAXA1','TAXA2','GENE'))
	save(tp, seqi, seq, file=file.path(wdir, gsub('.fa','_save160726beta.rda',infile)))
	#
	#	do we have higher bootstrap variance if there are more gaps by gene?
	#
	set(tp, NULL, 'GENE', tp[, factor(GENE, levels=c('gag','pol','env','gag+pol+env'), labels=c('gag','pol','env','gag+pol+env'))])	
	ggplot( subset(tp, GENE!='gag+pol+env' & GD_MEAN>0), aes(x=cut(DO/GENE_LEN, breaks=seq(0,1,0.01)), y=GD_SD) ) + 	
			geom_boxplot(outlier.shape=NA) +
			coord_cartesian(ylim=c(0,0.18)) +
			scale_y_continuous(expand=c(0,0)) +
			facet_grid(~GENE) + theme_bw() +
			labs(x='\noverlap between taxon pairs\n(% of sequence length)', y='std deviation in genetic distance\n')
	ggsave(file=file.path(wdir, gsub('.fa','_GDSD_by_overlap.pdf',infile)), w=10, h=5)
	
	ggplot( subset(tp, GENE!='gag+pol+env' & GD_MEAN>0), aes(x=cut(DO/GENE_LEN, breaks=seq(0,1,0.01)), y=GD_SD/GD_MEAN) ) + 	
			geom_boxplot(outlier.shape=NA) +
			coord_cartesian(ylim=c(0,0.6)) +
			scale_y_continuous(expand=c(0,0)) +
			facet_grid(~GENE) + theme_bw() +
			labs(x='\noverlap between taxon pairs\n(% of sequence length)', y='coefficient of variation in genetic distance\nacross bootstrap alignments\n')
	ggsave(file=file.path(wdir, gsub('.fa','_GDCOV_by_overlap.pdf',infile)), w=10, h=5)
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
	setnames(tp, c('V1','V2'), c('TAXA1','TAXA2'))
	tp[, IDX:= seq_len(nrow(tp))]
	#	subset by batch
	tp[, BATCH:= ceiling(IDX/batch.n)]
	if(!is.na(batch.i))
		tp	<- subset(tp, BATCH==batch.i)
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
treecomparison.ana.160627.standardize.KC<- function()
{	
	require(ggplot2)
	require(data.table)
	require(ape)
	require(scales)	
	require(ggtree)
	require(phangorn)
	
	edir			<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation'
	timetag			<- '160627'
	load(paste(edir,'/','submitted_160627_QDPDKC.rda',sep=''))
	
	
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
##	olli 27.06.16
##--------------------------------------------------------------------------------------------------------
treecomparison.submissions.160627<- function()	
{
	require(data.table)
	require(ape)
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
	MetaPIGA.trees			<- vector('list', sum(sapply(tmp.trees, length)))
	for(i in seq_along(tmp.trees))				
		for(j in seq_along(tmp.trees[[i]]))
			MetaPIGA.trees[[j+(i-1)*length(tmp.trees[[1]])]]	<- tmp.trees[[i]][[j]]	
	tmp						<- sapply( seq_along(tmp), function(i) paste(gsub('.txt','',tmp[i]), '_tree', seq_along(tmp.trees[[i]]), sep='') )
	names(MetaPIGA.trees)	<- as.vector(tmp)
	strs					<- c(strs, MetaPIGA.trees)
	#
	#	add MetaPIGA full genome trees, version 150831. stored as nexus
	#
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/MetaPIGA_150831'
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
	set(submitted.info, submitted.info[, which(TEAM=='MetaPIGA')], 'GENE', 'GAG+POL+ENV')	
	set(submitted.info, submitted.info[, which(TEAM=='IQTree' & grepl('[0-9]_partition', FILE))], 'GENE', 'GAG+POL+ENV')
	set(submitted.info, submitted.info[, which(TEAM=='IQTree' & grepl('[0-9]_pol_partition', FILE))], 'GENE', 'POL')
	set(submitted.info, submitted.info[, which(TEAM=='IQTree' & grepl('[0-9]_gag_partition', FILE))], 'GENE', 'GAG')
	stopifnot(nrow(subset(submitted.info, TEAM=='IQTree' & is.na(GENE)))==0)	
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
	###
	tmp		<- subset(submitted.info, TEAM=='PhyML' & MODEL=='R')[, IDX]
	for(i in tmp)
	{		
		cat(i,'\n')
		z	<- data.table(IDX=seq_along(strs[[i]]$tip.label), IDPOP=regmatches(strs[[i]]$tip.label, regexpr('IDPOP_[0-9]+',strs[[i]]$tip.label)), SC=subset(submitted.info, IDX==i)[,SC])
		z	<- merge(unique(subset(tinfo, BRL_T=='time', select=c(IDPOP,SC,TAXA))), z, by=c('IDPOP','SC'))
		setkey(z, IDX)
		stopifnot(nrow(z)==Ntip(strs[[i]]))
		strs[[i]]$tip.label	<- z[, TAXA]			
	}
	tmp		<- subset(submitted.info, TEAM=='PhyML' & MODEL=='V')[, IDX]
	for(i in tmp)
	{
		
		z	<- data.table(IDX=seq_along(strs[[i]]$tip.label), IDPOP=regmatches(strs[[i]]$tip.label, regexpr('HOUSE[0-9]+-[0-9]+|House[0-9]+-[0-9]+',strs[[i]]$tip.label)), SC=subset(submitted.info, IDX==i)[,SC])
		z	<- merge(unique(subset(tinfo, BRL_T=='time', select=c(IDPOP,SC,TAXA))), z, by=c('IDPOP','SC'))
		stopifnot(nrow(z)==length(strs[[i]]$tip.label))
		setkey(z, IDX)
		strs[[i]]$tip.label	<- z[, TAXA]		
	}
	#
	#	re-root simulated trees with rtt
	#		
	options(warn=2)
	strs_rtt	<- lapply(seq_along(strs), function(i)
			{
				cat('\n',i)
				#i	<- 628 ; i<- 241; i<- 571
				ph	<- strs[[i]]
				tmp	<- data.table(TAXA=ph$tip.label)	
				tmp[, T_SEQ:= tmp[, regmatches(TAXA, regexpr('[0-9]*\\.[0-9]+$|[0-9]+$', TAXA)) ]]
				#phr	<- rtt(ph, tmp[, as.numeric(T_SEQ)])
				phr	<- rtt(ph, tmp[, as.numeric(T_SEQ)], ncpu=4)
				phr
			})
	names(strs_rtt)	<- names(strs)
	
	#outdir		<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation'	
	#save(strs, strs_rtt, ttrs, tinfo, tfiles, file=file.path(outdir,'submitted_160713.rda'))
	
	#
	#	ladderize all trees
	#		
	ttrs	<- lapply(ttrs, ladderize)
	strs	<- lapply(strs, ladderize)
	strs_rtt<- lapply(strs_rtt, ladderize)	
	#
	#	plot simulated trees versus true tree
	#
	require(ggtree)	
	setkey(submitted.info, SC, TEAM, GENE)
	invisible(subset(submitted.info, OTHER=='N')[, 
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
				pdf(file=file.path(outdir, paste('strs_rtt_160727_',SC,'.pdf',sep='')), w=40, h=length(IDX)/10*12)
				print(p)
				dev.off()	
				NULL
			}, by=c('SC')])
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
	# compute on true trees the proportion if either transmitter or among recipients
	#
	tmp				<- treedist.closest.ind.obs(tinfo, gd.thresh=Inf)
	setnames(tmp, 'TR_REC_perc_T', 'TR_REC_perc_T_Inf')	
	submitted.info	<- merge(submitted.info, tmp, by='SUB_IDX_T', all.x=1)
	tmp				<- treedist.closest.ind.obs(tinfo, gd.thresh=0.045)
	setnames(tmp, 'TR_REC_perc_T', 'TR_REC_perc_T_45')	
	submitted.info	<- merge(submitted.info, tmp, by='SUB_IDX_T', all.x=1)		
	tinfo.pairs		<- treedist.closest.ind.obs(tinfo, gd.thresh=Inf, rtn.pairs=TRUE)
	setnames(tinfo.pairs, 'TRUE_PAIR','TRUE_PAIR_Inf')
	#
	# compute closest individual on simulated trees and determine proportion if either transmitter or among recipients
	#	
	tmp				<- treedist.closest.ind.reconstructed(submitted.info, tinfo, gd.thresh=Inf)
	setnames(tmp, 'TR_REC_perc', 'TR_REC_perc_Inf')		
	submitted.info	<- merge(submitted.info, tmp, by='IDX', all.x=1)	
	tmp				<- treedist.closest.ind.reconstructed(submitted.info, tinfo, gd.thresh=0.045)
	setnames(tmp, 'TR_REC_perc', 'TR_REC_perc_45')		
	submitted.info	<- merge(submitted.info, tmp, by='IDX', all.x=1)	
	tmp				<- treedist.closest.ind.reconstructed.oftruepairs(submitted.info, tinfo.pairs)
	submitted.info	<- merge(submitted.info, tmp, by='IDX', all.x=1)
	#
	#	compute Robinson Fould of complete tree 
	#
	tmp				<- treedist.robinsonfould.wrapper(submitted.info, ttrs, strs_rtt)
	submitted.info	<- merge(submitted.info, tmp, by='IDX')
	#	compute Robinson Fould of clusters, then take sum
	tmp				<- treedist.robinsonfouldclusters.wrapper(submitted.info, ttrs, strs, tinfo)
	sclu.info		<- merge(subset(submitted.info, select=c("IDX","SC","FILE","TEAM","MODEL","SEQCOV","ACUTE","GAPS","ART","EXT","BEST","OTHER","GENE","TAXAN","ROOTED","BRL","SUB_IDX_T","TIME_IDX_T","TAXAN_T")), tmp, by='IDX')
	#
	#	path distance of complete trees SSS
	#
	tmp				<- treedist.pathdifference.wrapper(submitted.info, ttrs, strs)
	submitted.info	<- merge(submitted.info, tmp, by='IDX')
	#set(submitted.info, NULL, c("SUB_IDX_T.y","SC.y","FILE.y","TEAM.y","MODEL.y","SEQCOV.y","ACUTE.y","GAPS.y","ART.y", "EXT.y", "GENE.y", "BEST.y","OTHER.y","TAXAN.y","ROOTED.y","BRL.y","TIME_IDX_T.y","TAXAN_T.y","TR_REC_perc_T_Inf.y","TR_REC_perc_T_45.y","TR_REC_perc_Inf.y","TR_REC_perc_45.y","RF.y","NRF.y"), NULL)	
	#setnames(submitted.info, c("SUB_IDX_T.x", "SC.x", "FILE.x", "TEAM.x", "MODEL.x", "SEQCOV.x", "ACUTE.x", "GAPS.x", "ART.x", "EXT.x", "GENE.x", "BEST.x", "OTHER.x", "TAXAN.x", "ROOTED.x", "BRL.x", "TIME_IDX_T.x", "TAXAN_T.x", "TR_REC_perc_T_Inf.x", "TR_REC_perc_T_45.x", "TR_REC_perc_Inf.x", "TR_REC_perc_45.x", "RF.x", "NRF.x"),
	#		c("SUB_IDX_T", "SC", "FILE", "TEAM", "MODEL", "SEQCOV", "ACUTE", "GAPS", "ART", "EXT", "GENE", "BEST", "OTHER", "TAXAN", "ROOTED", "BRL", "TIME_IDX_T", "TAXAN_T", "TR_REC_perc_T_Inf", "TR_REC_perc_T_45", "TR_REC_perc_Inf", "TR_REC_perc_45", "RF", "NRF"))
	#
	#	path distance of clusters
	#	
	tmp				<- treedist.pathdifference.clusters.wrapper(submitted.info, ttrs, strs, tinfo)	
	sclu.info		<- merge(sclu.info, tmp, by=c('IDX','IDCLU'))	
	#
	#	SAVE
	#
	outdir		<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation'	
	save(strs, strs_rtt, ttrs, tinfo, tfiles, tinfo.pairs, submitted.info, sclu.info, lba,  file=file.path(outdir,'submitted_160713.rda'))
	#
	#	ADD other summaries
	#
	load( file.path(outdir, 'submitted_160704_KC.rda') )
	sclu.info.kc	<- copy(sclu.info)
	load( file.path(outdir, 'submitted_160627_QDPD.rda') )
	sclu.info		<- merge(sclu.info, subset(sclu.info.kc, select=c(IDX, TEAM, GENE, BRL, IDCLU, KC)), by=c('IDX','TEAM','GENE','BRL','IDCLU'))
	save(strs, strs_rtt, ttrs, tinfo, submitted.info, sclu.info, lba,  file=file.path(outdir,'submitted_160627_QDPDKC.rda'))
}

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
				stopifnot( length(z)==abs(diff(c(TAXAN, TAXAN_T))) )
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
treecomparison.ana.160627<- function()
{	
	require(ggplot2)
	require(data.table)
	require(ape)
	require(scales)	
	require(ggtree)
	require(phangorn)
	
	edir			<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation'
	timetag			<- '160627'
	timetag			<- '160713'
	#load(paste(edir,'/','submitted_160627_QDPDKC.rda',sep=''))
	load(paste(edir,'/','submitted_160713.rda',sep=''))
	#
	#	long branches on regional
	#
	lba			<- merge(lba, subset(submitted.info, select=c(IDX, OTHER)), by='IDX')
	#	get reference of error without long branches
	lba[, ERR:= DEPTH_T-DEPTH]
	#	these are the closest trees in terms of NRF
	ref.box		<- subset(lba, IDX==858)[, quantile(ERR, p=c(0.25, 0.75))+c(-1,1)*3*diff(quantile(ERR, p=c(0.25, 0.75)))]
	ref.box		<- subset(lba, IDX==2)[, quantile(ERR, p=c(0.25, 0.75))+c(-1,1)*3*diff(quantile(ERR, p=c(0.25, 0.75)))]
	subset(lba, IDX==858)[, sd(ERR)*10]
	#	this suggests the following Tukey criterion:
	ref.box		<- c(-0.04,0.04)
	ref.box		<- c(-0.1,0.1)
	lba.su		<- subset(lba, OTHER=='N')[, list(TAXAN=length(ERR), OUTLIER_P=mean(ERR<ref.box[1] | ERR>ref.box[2])), by=c('MODEL','SC','TEAM','GAPS','GENE','IDX')]
	lba.su[, OUTLIER_N:=TAXAN*OUTLIER_P]
	ggplot(lba.su, aes(x=factor(GAPS, levels=c('none','low','high'), labels=c('none','low','high')), y=OUTLIER_P, colour=GENE)) + 
			geom_jitter(position=position_jitter(w=0.8, h = 0), size=1) + 
			scale_y_continuous(labels = scales::percent, limits=c(0, 1)) +
			facet_grid(~TEAM+GENE)  + theme_bw() +
			labs(	x='\ngaps', 
					y='branch length to root\n50% too small or too large\n(% of all taxa in tree)\n')
	file	<- file.path(edir, paste(timetag,'_','longbranches.pdf',sep=''))
	ggsave(file=file, w=15, h=7)
	#
	lba.su	<- subset(lba, OTHER=='N')[, list(TAXAN=length(ERR), OUTLIER_P=mean(ERR< -0.04 | ERR>0.04)), by=c('MODEL','SC','TEAM','GAPS','GENE','IDX')]
	ggplot(lba.su, aes(x=factor(GAPS, levels=c('none','low','high'), labels=c('none','low','high')), y=OUTLIER_P, colour=GENE)) + 
			geom_jitter(position=position_jitter(w=0.8, h = 0), size=1) + 
			scale_y_continuous(labels = scales::percent, limits=c(0, 1)) +
			facet_grid(~TEAM+GENE)  + theme_bw() +
			labs(	x='\ngaps', 
					y='branch length to root\n20% too small or too large\n(% of all taxa in tree)\n')
	file	<- file.path(edir, paste(timetag,'_','longbranches_20pc.pdf',sep=''))
	ggsave(file=file, w=15, h=7)
	
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
	tmp		<- merge(sa, subset(lba.su, OUTLIER_P<0.2, IDX), by='IDX')
	tmp		<- copy(sa)
	tmp		<- subset(tmp, !is.na(TR_REC_perc_45) & ACUTE=='10%' & TEAM%in%c('RAXML','IQTree'))
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
	#
	#	on clusters
	#
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
	#
	#	add size adjusted KC
	#	
	require(gamlss)
	kc.std.d	<- subset(sc, TEAM!='MetaPIGA' & SC=='sc 1')
	kc.std.m3	<- gamlss(KC~I(CLU_N*(CLU_N-1)/2), data=kc.std.d)
	sc			<- subset(sc, TEAM!='MetaPIGA' & SC%in%c('sc 1','sc 2','sc 4'))	
	sc[, KCadj:= KC / predict(kc.std.m3, data=kc.std.d, newdata=sc,type='response', se.fit=FALSE)]	
	ggplot(sc, aes(x=CLU_N)) + geom_point(aes(y=KCadj,colour=GENE, pch=TEAM)) +	scale_y_log10() + facet_grid(~SC)	
	tmp		<- melt(subset(sc, IDX==45), measure.var=c('NPD','NPDSQ','NRFC','NQDC','KCadj'))
	ggplot( tmp, aes(x=CLU_N, y=value, colour=GENE, pch=TEAM)) + geom_point() + facet_grid(GENE+TEAM+IDX~variable)
	#
	#	check dependence on size of cluster
	#
	ggplot( melt(sc, measure.var=c('NPD','NPDSQ','NRFC','NQDC','KCadj')), aes(x=CLU_N, y=value, colour=GENE, pch=TEAM)) + geom_point() + facet_grid(variable+TEAM+IDX~GENE, scales='free_y')
	file	<- file.path(edir, paste(timetag,'_','dependence_on_clustersize.pdf',sep=''))
	ggsave(file=file, w=10, h=1000, limitsize = FALSE)
	#
	#
	#	
	sc		<- sc[, list( 	NRFme=mean(NRFC, na.rm=TRUE), 	NQDme=mean(NQDC, na.rm=TRUE), 	NPDme=mean(NPD, na.rm=TRUE), 	NPDSQme=mean(NPDSQ, na.rm=TRUE),		KCAme=mean(KCadj, na.rm=TRUE),
						  	NRFmd=median(NRFC, na.rm=TRUE), 	NQDmd=median(NQDC, na.rm=TRUE), NPDmd=median(NPD, na.rm=TRUE), 	NPDSQmd=median(NPDSQ, na.rm=TRUE),	KCAmd=median(KCadj, na.rm=TRUE)
							), by=c('SC','GENE','TEAM','BEST','IDX','FILE','GAPS','MODEL','TAXAN','TAXAN_T','ROOTED','SEQCOV','ART','ACUTE','EXT','OTHER')]
	sc		<- subset(sc, MODEL=='Model: Regional')
	
	ggplot(subset(sc, ACUTE=='low' & TEAM!='MetaPIGA'), aes(x=GAPS)) +
			geom_jitter(aes(y=NQDme, colour=GENE, pch=TEAM), position=position_jitter(w=0.8, h = 0), size=2) +			
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 
			scale_y_continuous(labels = scales::percent, expand=c(0,0), limits=c(0, 0.4)) +
			scale_shape_manual(values=c('IQTree'=15, 'PhyML'=12, 'RAXML'=8)) +
			labs(	x='\nGappiness of full-genome sequences', 
					y='Quartett distance\n(standardized)\n',
					colour='part of genome used\nfor tree reconstruction',
					pch='algorithm') +
			theme_bw() + theme(legend.position='bottom') 
	file	<- file.path(edir, paste(timetag,'_','QD_clumean_polvsall_by_gaps_taxan1600_Acute10pc.pdf',sep=''))
	ggsave(file=file, w=5, h=7)
	
	ggplot(subset(sc, ACUTE=='low' & TEAM!='MetaPIGA'), aes(x=GAPS)) +
			geom_jitter(aes(y=NRFme, colour=GENE, pch=TEAM), position=position_jitter(w=0.8, h = 0), size=2) +			
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 
			scale_y_continuous(labels = scales::percent, expand=c(0,0), limits=c(0, 0.6)) +
			scale_shape_manual(values=c('IQTree'=15, 'PhyML'=12, 'RAXML'=8)) +
			labs(	x='\nGappiness of full-genome sequences', 
					y='Robinson Fould\n(standardized)\n',
					colour='part of genome used\nfor tree reconstruction',
					pch='algorithm') +
			theme_bw() + theme(legend.position='bottom') 
	file	<- file.path(edir, paste(timetag,'_','RF_clumean_polvsall_by_gaps_taxan1600_Acute10pc.pdf',sep=''))
	ggsave(file=file, w=5, h=7)
	
	ggplot(subset(sc, ACUTE=='low' & TEAM!='MetaPIGA'), aes(x=GAPS)) +
			geom_jitter(aes(y=NPDSQme, colour=GENE, pch=TEAM), position=position_jitter(w=0.8, h = 0), size=2) +			
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 
			scale_y_continuous(expand=c(0,0), limits=c(0, 2)) +
			scale_shape_manual(values=c('IQTree'=15, 'PhyML'=12, 'RAXML'=8)) +
			labs(	x='\nGappiness of full-genome sequences', 
					y='Path distance\n(upper bound is 2)\n',
					colour='part of genome used\nfor tree reconstruction',
					pch='algorithm') +
			theme_bw() + theme(legend.position='bottom') 
	file	<- file.path(edir, paste(timetag,'_','PD_clumean_polvsall_by_gaps_taxan1600_Acute10pc.pdf',sep=''))
	ggsave(file=file, w=5, h=7)
	
	ggplot(subset(sc, ACUTE=='low' & TEAM!='MetaPIGA'), aes(x=GAPS)) +
			geom_jitter(aes(y=KCAme, colour=GENE, pch=TEAM), position=position_jitter(w=0.8, h = 0), size=2) +			
			scale_colour_manual(values=c('gag'='red','pol'="grey60", 'gag+pol+env'="#3F4788FF")) + 
			scale_y_continuous(expand=c(0,0), limits=c(0, 2)) +
			scale_shape_manual(values=c('IQTree'=15, 'PhyML'=12, 'RAXML'=8)) +
			labs(	x='\nGappiness of full-genome sequences', 
					y='KC distance\n(adjusted for dependence on sub-tree size)\n',
					colour='part of genome used\nfor tree reconstruction',
					pch='algorithm') +
			theme_bw() + theme(legend.position='bottom') 
	file	<- file.path(edir, paste(timetag,'_','KC_clumean_polvsall_by_gaps_taxan1600_Acute10pc.pdf',sep=''))
	ggsave(file=file, w=5, h=7)
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