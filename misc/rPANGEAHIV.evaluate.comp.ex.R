##--------------------------------------------------------------------------------------------------------
##	evaluate results
##	olli 12.08.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate<- function()
{
	require(RColorBrewer)
	dfa		<- project.PANGEA.TEST.pipeline.Aug2015.evaluate.read()
	#	check for updated submissions, and keep 
	dfa		<- project.PANGEA.TEST.pipeline.Aug2015.keep.most.recent.submission(dfa, format='%d.%m.%Y')
	#	save submissions
	outdir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim_internal/documents/external/2015_08_results/results'
	save(dfa, file=paste(outdir,'/submissions.R',sep=''))
	load(paste(outdir,'/submissions.R',sep=''))
	#	add objective legend
	dfa		<- merge(dfa, data.table(USED_GENES=c('pol','all'), USED_GENES_L=c('pol gene','pol+gag+env\ngenome') ), by='USED_GENES')
	set(dfa, NULL, 'TEAM', dfa[, factor(TEAM)])
	tmp		<- data.table( 	OBJ=	c('OBJ_i','OBJ_ii','OBJ_iii','OBJ_iv','OBJ_v','OBJ_vi'),
			OBJ_L=	c('Incidence\nTrend', '%Incidence', 'Incidence\nreduction', '%Acute Ctgr\n(baseline)', '%Acute\n(baseline)', '%Acute\n(endpoint)'))
	set(tmp, NULL, 'OBJ_L2', tmp[, factor(OBJ_L, levels=OBJ_L, labels=OBJ_L)])
	set(tmp, NULL, 'OBJ_L', tmp[, factor(OBJ_L, levels=rev(OBJ_L), labels=rev(OBJ_L))])
	dfa		<- merge(dfa, tmp, by='OBJ')
	#	add data legend
	dfa[, DATA_T2:='NA_character_']
	set(dfa, dfa[, which(DATA_T=='seq')], 'DATA_T2', 'using\nsequences')
	set(dfa, dfa[, which(DATA_T=='phy')], 'DATA_T2', 'using\ntrue tree')
	set(dfa, NULL, 'DATA_T2', dfa[, factor(DATA_T2, levels=rev(c('using\nsequences','using\ntrue tree')), labels=rev(c('using\nsequences','using\ntrue tree')))])		
	#	add scenario type
	set(dfa, NULL, 'DATA_T', dfa[, factor(DATA_T, levels=c('seq','phy'), labels=c('seq','phy'))])
	set(dfa, NULL, 'INT_T', dfa[, factor(INT_T, levels=c('fast','slow','none'), labels=c('fast','slow','none'))])
	set(dfa, NULL, 'AC_T', dfa[, factor(AC_T, levels=c('low','high'), labels=c('low','high'))])
	set(dfa, NULL, 'IMPRT', dfa[, factor(IMPRT*100, levels=c(0,2,5,20), labels=paste(c(0,2,5,20),'%',sep=''))])
	set(dfa, NULL, 'SMPL_C', dfa[, factor(SMPL_C*100, levels=c(8, 16, 30, 60), labels=paste(c(8, 16, 30, 60),'%',sep=''))])
	set(dfa, NULL, 'SMPL_D', dfa[, factor(SMPL_D, levels=c(5,3), labels=c(5,3))])	
	set(dfa, dfa[, which(SMPL_M=='overs')], 'SMPL_M', 'much')
	set(dfa, dfa[, which(SMPL_M=='extrs')], 'SMPL_M', 'extreme')
	set(dfa, dfa[, which(is.na(SMPL_M))], 'SMPL_M', 'extreme')
	set(dfa, NULL, 'SMPL_M', dfa[, factor(SMPL_M, levels=c('much','extreme'), labels=c('much','extreme'))])	
	tmp		<- unique(subset( dfa, select=c(DATAT_L, SC_RND, DATA_T, SC, AC_T, INT_T, IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D) ))
	setkey(tmp, DATAT_L, AC_T, INT_T, DATA_T, IMPRT, SMPL_C, SMPL_D, SMPL_M)
	tmp[, SCENARIO_L:= paste('%AC=',AC_T,' ARTup=',INT_T,' EXT=',IMPRT,'\n',DATA_T,' ',SMPL_N,' ',SMPL_C,' ',SMPL_D,' ',SMPL_M, ' (',SC_RND,')',sep='')]
	dfa		<- merge(dfa, subset(tmp, select=c(SC_RND, SCENARIO_L)), by='SC_RND')
	#	add intervention legend
	dfa[, INT_L:= dfa[, paste('ART scale up\n',as.character(INT_T),sep='')]]
	setkey(dfa, INT_T)
	set(dfa, NULL, 'INT_L', dfa[, factor(INT_L, levels=dfa[, unique(INT_L)], labels=dfa[, unique(INT_L)])])
	#	add %Acute legend
	dfa[, AC_L:= dfa[, paste('%Acute\n',as.character(AC_T),sep='')]]
	setkey(dfa, AC_T)
	set(dfa, NULL, 'AC_L', dfa[, factor(AC_L, levels=dfa[, unique(AC_L)], labels=dfa[, unique(AC_L)])])
	#	get data.table of data sets ~ all primary and secondary objectives
	dfd		<- subset(dfa, select=c(SC_RND, DATA_T, DATAT_L, AC_T, INT_T, IMPRT, SMPL_N, SMPL_C,  SMPL_M, SMPL_D))
	setkey(dfd, DATAT_L, INT_T, AC_T, DATA_T, IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D)
	dfd		<- unique(dfd)
	#	Primary Objectives, on sequences 
	tmp		<- data.table(expand.grid(ANA='Pr_Seq', SC_RND=subset(dfd, DATA_T=='seq')[, SC_RND], stringsAsFactors=FALSE))
	dfr		<- copy(tmp)
	#	Primary Objectives, on trees
	tmp		<- subset(dfd, DATA_T=='phy' & SMPL_D==5 & (DATAT_L=='Regional' & IMPRT=='5%' & SMPL_M=='much' & SMPL_N==1600 | DATAT_L=='Village' & SMPL_C=='30%')  )
	dfr		<- rbind(dfr, data.table(expand.grid(ANA='Pr_Phy', SC_RND=tmp[, SC_RND], stringsAsFactors=FALSE)))
	#	Secondary: sequence coverage
	tmp		<- subset(dfd, DATA_T=='phy' & SMPL_D==5 & INT_T!='none' & (DATAT_L=='Regional' & IMPRT=='5%' & SMPL_M=='much' | DATAT_L=='Village')  )
	dfr		<- rbind(dfr, data.table(expand.grid(ANA='Sc_SeqCoverage_Phy', SC_RND=tmp[, SC_RND], stringsAsFactors=FALSE)))
	#	Secondary: imports
	tmp		<- subset(dfd, DATA_T=='phy' & SMPL_D==5 & INT_T!='none' & AC_T=='high' & (DATAT_L=='Regional' & SMPL_M=='much' & SMPL_C=='8%')  )
	dfr		<- rbind(dfr, data.table(expand.grid(ANA='Sc_Imports_Phy', SC_RND=tmp[, SC_RND], stringsAsFactors=FALSE)))
	#	Secondary: focussed sampling
	tmp		<- subset(dfd, DATA_T=='phy' & SMPL_D==5 & INT_T!='none' & AC_T=='low' & SMPL_C=='8%' & IMPRT=='5%' & DATAT_L=='Regional'  )
	dfr		<- rbind(dfr, data.table(expand.grid(ANA='Sc_SmplFc_Phy', SC_RND=tmp[, SC_RND], stringsAsFactors=FALSE)))
	#	Secondary: sampling duration
	tmp		<- subset(dfd, DATA_T=='phy' & INT_T!='none' & DATAT_L=='Regional' & SMPL_M=='much' & SMPL_C=='8%' & IMPRT=='5%' & INT_T=='fast')
	dfr		<- rbind(dfr, data.table(expand.grid(ANA='Sc_SmplD_Phy', SC_RND=tmp[, SC_RND], stringsAsFactors=FALSE)))
	#	merge with dfa
	dfr		<- dcast.data.table(dfr, SC_RND~ANA, value.var='SC_RND')
	set(dfr, NULL, 'Pr_Phy', dfr[, as.numeric(!is.na(Pr_Phy))])	
	set(dfr, NULL, 'Pr_Seq', dfr[, as.numeric(!is.na(Pr_Seq))])	
	set(dfr, NULL, 'Sc_Imports_Phy', dfr[, as.numeric(!is.na(Sc_Imports_Phy))])
	set(dfr, NULL, 'Sc_SeqCoverage_Phy', dfr[, as.numeric(!is.na(Sc_SeqCoverage_Phy))])
	set(dfr, NULL, 'Sc_SmplD_Phy', dfr[, as.numeric(!is.na(Sc_SmplD_Phy))])
	set(dfr, NULL, 'Sc_SmplFc_Phy', dfr[, as.numeric(!is.na(Sc_SmplFc_Phy))])	
	dfr		<- merge(unique(subset(dfa, select=c(DATAT_L,SC_RND))), dfr, by='SC_RND')
	dfa		<- merge(dfa, dfr, by=c('SC_RND','DATAT_L'))
	#
	#	reset Vancouver to UBC
	#
	set( dfa, dfa[, which(TEAM=='Vancouver')], 'TEAM', 'British Columbia' )	
	#
	#	set team color
	#
	tmp				<- c('Cambridge','Cambridge/Imperial','ETH Zurich','Imperial','British Columbia','True','Cambridge/Imperial (chronos)','Cambridge/Imperial (lsd)','Cambridge/Imperial (mh15)','Cambridge/Imperial (mh30)','Cambridge/Imperial (merged)')
	TEAM_CL			<- c( brewer.pal(5, 'Set1'), 'black', brewer.pal(5, 'Set2') )
	names(TEAM_CL)	<- tmp	
	#	write info on scenario IDs to table
	dfi<- subset( dfa, select=c(SC_RND, DATAT_L, DATA_T, AC_T, INT_T, IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D) )
	setkey(dfi, DATAT_L, DATA_T, AC_T, INT_T, IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D)
	dfi	<- unique(dfi)
	tmp	<- subset(dfa, OBJ!='OBJ_iv' & TEAM!='True' & !grepl('(',TEAM,fixed=1))[, list(SUB_N=length(central)), by=c('SC_RND','TEAM')]
	tmp	<- dcast.data.table(tmp, SC_RND~TEAM, value.var='SUB_N', fill=0L)
	dfi	<- merge(dfi, tmp, by='SC_RND')
	setkey(dfi, DATAT_L, DATA_T, AC_T, INT_T, IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D)	
	file<- paste(outdir,'/SC_RND_info.csv',sep='')
	write.csv(dfi, file=file, row.names=FALSE)
	#
	melt(dfi, measure.vars=c('Cambridge','Cambridge/Imperial','ETH Zurich','Imperial','British Columbia'))[, list(SUM=sum(value), COMPLETE=sum(value)/(5*33)), by='variable']
	#
	project.PANGEA.TEST.pipeline.Aug2015.evaluate.primary.sortedIncidenceRatio(dfa, outdir)
	project.PANGEA.TEST.pipeline.Aug2015.evaluate.primary.IncidenceTrends(dfa, outdir)
	project.PANGEA.TEST.pipeline.Aug2015.evaluate.primary.IncidenceRatioCovariates(dfa, outdir)
	project.PANGEA.TEST.pipeline.Aug2015.evaluate.primary.sortedPCIncidence(dfa, outdir)
	project.PANGEA.TEST.pipeline.Aug2015.evaluate.primary.PCIncidenceCovariates(dfa, outdir)
	#
	project.PANGEA.TEST.pipeline.Aug2015.evaluate.overallnumbers(dfa, outdir)
	#
	#	no results on pol submitted, focus on full genome only. get sample sizes per objective.
	project.PANGEA.TEST.pipeline.Aug2015.evaluate.samplesize(dfr, dfa, outdir)
	#
	tmp	<- subset(dfa, select=c(TEAM, DATA_T, DATAT_L, SIM_SCENARIO))
	setkey(tmp, TEAM, DATA_T, DATAT_L, SIM_SCENARIO)
	tmp	<- unique(tmp)
	tmp[, table(TEAM, DATA_T, DATAT_L )]
	#	for each primary objective
	#	compare results across teams
	project.PANGEA.TEST.pipeline.Aug2015.evaluate.primary.incidence(dfa, outdir, onSeq=1)
	project.PANGEA.TEST.pipeline.Aug2015.evaluate.primary.incidence(dfa, outdir, onSeq=0)
	project.PANGEA.TEST.pipeline.Aug2015.evaluate.primary.acute(dfa, outdir, onSeq=1)
	project.PANGEA.TEST.pipeline.Aug2015.evaluate.primary.acute(dfa, outdir, onSeq=0)
	#	for each secondary objective
	#	compare results across teams	
	project.PANGEA.TEST.pipeline.Aug2015.evaluate.secondary.imports(dfa, outdir)
	project.PANGEA.TEST.pipeline.Aug2015.evaluate.secondary.focussedsampling(dfa, outdir)
	project.PANGEA.TEST.pipeline.Aug2015.evaluate.secondary.sduration(dfa, outdir)
	project.PANGEA.TEST.pipeline.Aug2015.evaluate.secondary.seqcoverage(dfa, outdir)
	#	for each team
	#	all results
	invisible(sapply(setdiff(dfa[, unique(TEAM)],'True'), function(x)
					{		
						#x	<- 'Imperial'
						df		<- subset(dfa, (TEAM=='True' | TEAM==x) & USED_GENES=='all')
						set(df, df[, which(TEAM==x)], 'TEAM', 'estimate')
						set(df, df[, which(TEAM=='True')], 'TEAM', 'true value')
						set(df, NULL, 'TEAM', df[, factor(TEAM, levels=c('estimate','true value'), labels=c('estimate','true value'))])
						ggplot(df, aes(y=SCENARIO_L, x=central, xmin=lower95, xmax=upper95, colour=TEAM, pch=TEAM)) + 
								geom_errorbarh(height=0.3) + geom_point(size=3) + 
								scale_colour_manual(values = c("red","black")) +
								scale_shape_manual(values = c(13,18), guide = FALSE) +
								labs(x='', y='', title= paste('TEAM',x,'\n'), colour='')  +
								facet_grid(DATAT_L~OBJ_L2, scales='free', space='free_y') +
								theme_bw() + theme(legend.position='bottom')
						ggsave(file=paste(outdir,'/res_obj_TEAM_',gsub(' ','_',gsub('\\/|\\(|\\)','',x)),'.pdf',sep=''), w=14, h=0.5*df[, length(unique(SCENARIO_L))])
						#	results using seq data
						df		<- subset(dfa, (TEAM=='True' | TEAM==x) & USED_GENES=='all' & DATA_T=='seq')
						set(df, df[, which(TEAM==x)], 'TEAM', 'estimate')
						set(df, df[, which(TEAM=='True')], 'TEAM', 'true value')
						set(df, NULL, 'TEAM', df[, factor(TEAM, levels=c('estimate','true value'), labels=c('estimate','true value'))])
						ggplot(df, aes(y=SCENARIO_L, x=central, xmin=lower95, xmax=upper95, colour=TEAM, pch=TEAM)) + 
								geom_errorbarh(height=0.3) + geom_point(size=3) + 
								scale_colour_manual(values = c("red","black")) +
								scale_shape_manual(values = c(13,18), guide = FALSE) +
								labs(x='', y='', title= paste('TEAM',x,'\n'), colour='')  +
								facet_grid(DATAT_L~OBJ_L2, scales='free', space='free_y') +
								theme_bw() + theme(legend.position='bottom')
						ggsave(file=paste(outdir,'/res_objonseq_TEAM_',gsub(' ','_',gsub('\\/|\\(|\\)','',x)),'.pdf',sep=''), w=14, h=0.7*df[, length(unique(SCENARIO_L))])	
					}))
}
##--------------------------------------------------------------------------------------------------------
##	evaluate results
##	olli 12.08.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.overallnumbers<- function(dfa, outdir)
{
	#	count total submissions primary vs secondary
	tmp		<- subset(dfa, TEAM!='True' & !grepl('(', TEAM, fixed=1))
	tmp		<- tmp[, list(	Village=length(which(grepl('Vill',SIM_SCENARIO))), Regional=length(which(grepl('Regional',SIM_SCENARIO)))), by=c('TEAM','OBJ_L','USED_GENES_L')]	
	tmp		<- melt(tmp, measure.vars=c('Village','Regional'))	
	ggplot(tmp, aes(x=OBJ_L, y=value, fill=TEAM)) + geom_bar(stat='identity') +
			facet_grid(USED_GENES_L~variable) +			
			guides(fill=guide_legend(ncol=2)) +
			scale_fill_manual(values=TEAM_CL) +
			labs(x='', y='submissions\n(#)', title='Total scenarios submitted\n(using sequence data or true trees)\n') +
			theme_bw()+ theme(legend.position='bottom') + coord_flip()	
	ggsave(file=paste(outdir,'/res_scenarios_total.pdf',sep=''), w=10, h=8)
	
	#	count all submissions for primary objectives
	tmp		<- subset(dfa, TEAM!='True' & !grepl('(', TEAM, fixed=1) & DATA_T=='seq')
	tmp		<- tmp[, list(	Village=length(which(grepl('Vill',SIM_SCENARIO))), Regional=length(which(grepl('Regional',SIM_SCENARIO)))), by=c('TEAM','OBJ_L')]	
	tmp		<- melt(tmp, measure.vars=c('Village','Regional'))	
	ggplot(tmp, aes(x=OBJ_L, y=value, fill=TEAM)) + geom_bar(stat='identity') +
			facet_grid(~variable) +
			labs(x='', y='submissions\n(#)', title='Total scenarios submitted\n(using sequence data)\n') +
			scale_fill_manual(values=TEAM_CL) +
			guides(fill=guide_legend(ncol=2)) +
			theme_bw() + theme(legend.position='bottom') + coord_flip()
	ggsave(file=paste(outdir,'/res_scenarios_total_seqonly.pdf',sep=''), w=10, h=5)
	
	#	count complete submissions for primary objectives
	tmp		<- subset(dfa, TEAM!='True' & !grepl('(', TEAM, fixed=1) & DATA_T=='seq')
	tmp		<- tmp[, list(	Village=as.numeric(length(setdiff(c('01','02','03','04'),SC_RND))==0), Regional=as.numeric(length(setdiff(c('A','B','C','D'),SC))==0)), by=c('TEAM','OBJ_L','USED_GENES_L')]	
	tmp		<- melt(tmp, measure.vars=c('Village','Regional'))	
	ggplot(tmp, aes(x=OBJ_L, y=value, fill=TEAM)) + geom_bar(stat='identity') +
			facet_grid(USED_GENES_L~variable) +
			scale_y_continuous(breaks=seq(1,10,1), minor_breaks=NULL) +
			scale_fill_manual(values=TEAM_CL) +
			labs(x='', y='complete set of 4 submissions\n(#)', title='Complete submissions to evalute primary objectives\n(either village or regional)') +
			guides(fill=guide_legend(ncol=2)) +
			theme_bw() + theme(legend.position='bottom') + coord_flip()
	ggsave(file=paste(outdir,'/res_scenarios_total_seqonlycomplete.pdf',sep=''), w=10, h=7)
	
}
##--------------------------------------------------------------------------------------------------------
##	olli 12.08.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.keep.most.recent.submission<- function(dfa, format='%d.%m.%Y')
{
	
	tmp		<- dfa[, list(SUB_N=length(unique(SUBMISSION_DATE)), SUB_DATE=unique(SUBMISSION_DATE)), by=c('TEAM','DATAT_L','OBJ')]
	tmp		<- subset(tmp, SUB_N>1)
	set(tmp, NULL, 'SUB_DATE', tmp[,as.Date(SUB_DATE, format=format)])
	#	for each objective, determine submissions that are to be discarded
	tmp		<- tmp[, list(SUB_DATE=SUB_DATE[SUB_DATE!=max(SUB_DATE)]), by=c('TEAM','DATAT_L','OBJ')]
	for(i in seq_len(nrow(tmp)))
		set(dfa, dfa[, which(TEAM==tmp$TEAM[i] & DATAT_L==tmp$DATAT_L[i] & OBJ==tmp$OBJ[i] & SUBMISSION_DATE==tmp[, as.character(tmp$SUB_DATE[i], format=format)])], 'central', NA_real_)
	subset(dfa, !is.na(central))	
}
##--------------------------------------------------------------------------------------------------------
##	olli 12.08.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.primary.incidence<- function(dfa, outdir, onSeq=1)
{
	require(gridExtra)
	#,'OBJ_v','OBJ_vi'
	#	compare objectives with / without seq data, village + regional	
	if(onSeq)
	{
		df		<- subset(dfa, Pr_Seq==1)
		title	<- '\nPrimary objective\nIncidence from sequence data\n'
	}		
	if(!onSeq)
	{
		df	<- subset(dfa, Pr_Phy==1)
		title	<- '\nPrimary objective\nIncidence when phylogeny known\n'
	}			
	tmp	<- subset(df, !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_ii'))	
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')		
	p1		<- ggplot(subset(tmp, TEAM!='True'), aes(y=gsub('\n',':',paste(INT_L,'  ',AC_L,sep='')), x=central)) +
			#geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) + 
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			coord_cartesian(xlim=c(0, 10)) +
			scale_x_continuous(breaks=seq(1,9,1)) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM~DATAT_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\n%Incidence', y='')
	tmp	<- subset(df, !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_iii'))	
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')		
	p2		<- ggplot(subset(tmp, TEAM!='True'), aes(y=gsub('\n',':',paste(INT_L,'  ',AC_L,sep='')), x=central)) +
			#geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) + 
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			geom_vline(xintercept=1, colour='grey50', size=0.8) +
			coord_cartesian(xlim=c(0, 2)) +
			scale_colour_manual(values=TEAM_CL) +
			scale_x_continuous(breaks=seq(0.2,1.8,0.4), minor_breaks=seq(0.2,1.8,0.2)) +
			facet_grid(TEAM~DATAT_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\nIncidence reduction', y='')
	pdf(file=paste(outdir,'/res_acrossTEAM_PrimaryIncidence_onSeq',onSeq,'.pdf',sep=''), width = 15, height = 8)
	print(grid.arrange(p1, p2, nrow=1, main=title))
	dev.off() 
}
##--------------------------------------------------------------------------------------------------------
##	olli 12.08.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.primary.acute2<- function(dfa, outdir, onSeq=1)
{
	require(gridExtra)
	if(onSeq)
	{
		df		<- subset(dfa, Pr_Seq==1)
		title	<- '\nPrimary objective\n%Acute from sequence data\n'
	}		
	if(!onSeq)
	{
		df	<- subset(dfa, Pr_Phy==1)
		title	<- '\nSecondary objective\n%Acute when phylogeny known\n'
	}			
	tmp		<- subset(df, !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_v'))	
	tmp3	<- subset(tmp, TEAM=='True', select=c(SC_RND, central))
	setnames(tmp3, 'central', 'central_true')
	tmp		<- subset(merge(tmp, tmp3, by='SC_RND', all.x=TRUE), TEAM!='True')
	
	ggplot(tmp, aes(x=central_true, y=central, colour=TEAM, pch= gsub('\n',':',paste(INT_L,'  ',AC_L,sep='')))) +
			geom_point() +
			geom_abline(slope=1, intercept=0) +
			scale_colour_manual(values=TEAM_CL) +
			facet_wrap(~DATAT_L, scales='free') +
			theme_bw() + theme(panel.margin.x= unit(1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3), pch=guide_legend(ncol=3)) +
			labs(x= '\ntrue % Acute at baseline', y='\nestimated % Acute at baseline', pch='')	
}
##--------------------------------------------------------------------------------------------------------
##	olli 25.11.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.primary.IncidenceTrends<- function(dfa, outdir)
{
	tmp		<- subset(dfa, USED_GENES=='all' & TEAM!='True' & !grepl('(', TEAM,fixed=1) & OBJ=='OBJ_iii')
	tmp2	<- subset(dfa, TEAM=='True' & OBJ=='OBJ_iii')	
	setkey(tmp2, central)
	set(tmp, NULL, 'SC_RND', factor(tmp[,as.character(SC_RND)], levels=tmp2[, SC_RND], labels=tmp2[, SC_RND]))
	set(tmp2, NULL, 'SC_RND', factor(tmp2[,as.character(SC_RND)], levels=tmp2[, SC_RND], labels=tmp2[, SC_RND]))
	setnames(tmp2, 'central', 'central_t')
	tmp2[, TEAM:=NULL]
	tmp2	<- merge( data.table(expand.grid(SC_RND=tmp[, unique(SC_RND)], TEAM=tmp[, unique(TEAM)])), tmp2, by='SC_RND')
	tmp		<- merge(tmp, subset(tmp2, select=c(SC_RND, TEAM, central_t)), by=c('SC_RND','TEAM'), all=1)	
	tmp2	<- subset(dfa, USED_GENES=='all' & TEAM!='True' & !grepl('(', TEAM,fixed=1) & OBJ=='OBJ_i', c(SC_RND, TEAM, central))
	setnames(tmp2, 'central', 'OBJ_i')
	set(tmp2, NULL, 'SC_RND', factor(tmp2[,as.character(SC_RND)], levels=tmp[, levels(SC_RND)], labels=tmp[, levels(SC_RND)]))
	tmp		<- merge(tmp, tmp2, by=c('SC_RND','TEAM'), all=1)
	set(tmp, tmp[, which(is.na(OBJ_i))], 'OBJ_i', 2)
	
	
	cnts	<- as.data.table(melt(tmp[,  table(central_t<0.75, as.character(OBJ_i), as.character(TEAM)) ]))
	setnames(cnts, c('Var1','Var2','Var3','value'), c('TRUE_TREND','OBJ_i','TEAM','CNT'))
	set(cnts, NULL, 'OBJ_i', cnts[, factor(OBJ_i, levels=c(-1,0,1,2), labels=c('declining','stable','increasing','missing'))])
	set(cnts, NULL, 'TRUE_TREND', cnts[, factor(TRUE_TREND, levels=c(FALSE,TRUE), labels=c('g75pc','l75pc'))])
	cnts	<- subset(cnts, TRUE_TREND=='l75pc')
	tmp2	<- cnts[, list(CNT= round(100*CNT[OBJ_i=='declining']/sum(CNT),d=0), OBJ_i='TPR' ), by=c('TRUE_TREND','TEAM')]
	cnts	<- rbind(cnts, tmp2, use.names=TRUE)	
	cnts	<- dcast.data.table( subset(cnts, TRUE_TREND=='l75pc'), TRUE_TREND+OBJ_i~TEAM, value.var='CNT')
	ans		<- copy(cnts)	
	cnts	<- as.data.table(melt(tmp[,  table(central_t>0.85, as.character(OBJ_i), as.character(TEAM)) ]))
	setnames(cnts, c('Var1','Var2','Var3','value'), c('TRUE_TREND','OBJ_i','TEAM','CNT'))
	set(cnts, NULL, 'OBJ_i', cnts[, factor(OBJ_i, levels=c(-1,0,1,2), labels=c('declining','stable','increasing','missing'))])
	set(cnts, NULL, 'TRUE_TREND', cnts[, factor(TRUE_TREND, levels=c(FALSE,TRUE), labels=c('l85pc','g85pc'))])
	cnts	<- subset(cnts, TRUE_TREND=='g85pc')
	tmp2	<- cnts[, list(CNT= round(100*CNT[OBJ_i=='declining']/sum(CNT),d=0), OBJ_i='FPR' ), by=c('TRUE_TREND','TEAM')]
	cnts	<- rbind(cnts, tmp2, use.names=TRUE)	
	cnts	<- dcast.data.table( subset(cnts, TRUE_TREND=='g85pc'), TRUE_TREND+OBJ_i~TEAM, value.var='CNT')	
	ans		<- rbind(ans, cnts)
	write.csv(ans, file=paste(outdir, '/', 'res_acrossTEAM_Primary_IncidenceTrend.csv', sep=''), row.names=FALSE)
}
##--------------------------------------------------------------------------------------------------------
##	olli 26.11.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.primary.CovariatesPCIncidence<- function(dfa, outdir)
{
	tmp		<- subset(dfa, USED_GENES=='all' & TEAM!='True' & !grepl('(', TEAM,fixed=1) & OBJ=='OBJ_ii')
	tmp2	<- subset(dfa, TEAM=='True' & OBJ=='OBJ_ii')	
	setkey(tmp2, central)
	set(tmp, NULL, 'SC_RND', factor(tmp[,as.character(SC_RND)], levels=tmp2[, SC_RND], labels=tmp2[, SC_RND]))
	set(tmp2, NULL, 'SC_RND', factor(tmp2[,as.character(SC_RND)], levels=tmp2[, SC_RND], labels=tmp2[, SC_RND]))
	setnames(tmp2, 'central', 'central_t')
	tmp2[, TEAM:=NULL]
	tmp2	<- merge( data.table(expand.grid(SC_RND=tmp[, unique(SC_RND)], TEAM=tmp[, unique(TEAM)])), tmp2, by='SC_RND')
	tmp		<- merge(tmp, subset(tmp2, select=c(SC_RND, TEAM, central_t)), by=c('SC_RND','TEAM'), all=1)	
	tmp2	<- subset(dfa, USED_GENES=='all' & TEAM!='True' & !grepl('(', TEAM,fixed=1) & OBJ=='OBJ_i', c(SC_RND, TEAM, central))
	setnames(tmp2, 'central', 'OBJ_i')
	set(tmp2, NULL, 'SC_RND', factor(tmp2[,as.character(SC_RND)], levels=tmp[, levels(SC_RND)], labels=tmp[, levels(SC_RND)]))
	tmp		<- merge(tmp, tmp2, by=c('SC_RND','TEAM'), all=1)
	tmp		<- subset(tmp, !is.na(central), c(SC_RND, TEAM, DATAT_L, central, lower95, upper95, DATA_T, AC_T, INT_T, IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D, central_t))
	
	tmp[, RES:= log(central)-log(central_t)]
	tmp[, OUTLIER:= central>20]
	dfc		<- subset(tmp, !OUTLIER)[, list( COR=cor(central_t, central), COR_LOG=cor(log(central_t), log(central)), 											 
					MAE=mean( abs(central-central_t) ), MAE_L=mean( abs(log(central)-log(central_t)) ), 
					MSE=mean( (central-central_t)^2 ), MSE_L=mean( (log(central)-log(central_t))^2 ),
					ARME=mean(central-central_t), ARME_L=mean(log(central)-log(central_t))), by='TEAM']
	
	require(gamlss)
	bw.AIC			<- vector('list', tmp[, length(unique(TEAM))])
	names(bw.AIC)	<- tmp[, unique(TEAM)]
	bw.BIC			<- vector('list', tmp[, length(unique(TEAM))])
	names(bw.BIC)	<- tmp[, unique(TEAM)]
	for(x in names(bw.AIC))
	{
		cat('\nprocess',x)
		df			<- subset(tmp, TEAM==x)	
		set(df, NULL, c('lower95','upper95'), NULL)
		if(x!='Cambridge')
			mnoa	<- gamlss(RES~DATAT_L+DATA_T+AC_T+INT_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D-1, data=df, family=NO, trace=FALSE)
		if(x=='Cambridge')
			mnoa	<- gamlss(RES~AC_T+INT_T+IMPRT+SMPL_N-1, data=df, family=NO, trace=FALSE)
		bw.AIC[[x]]	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE)		
		bw.BIC[[x]] <- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE, k=log(nrow(df)))				
	}
	dfcoeff	<- tmp[, 	{
				z			<- names(coef(bw.AIC[[as.character(TEAM)]]))
				z2			<- names(coef(bw.BIC[[as.character(TEAM)]]))
				length(z)	<- max(length(z),length(z2))
				length(z2)	<- max(length(z),length(z2))
				list( BWAIC=z, BWBIC=z2 )	
			}, by='TEAM']	
}
##--------------------------------------------------------------------------------------------------------
##	olli 26.11.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.primary.sortedPCIncidence<- function(dfa, outdir)
{
	tmp		<- subset(dfa, USED_GENES=='all' & TEAM!='True' & !grepl('(', TEAM,fixed=1) & OBJ=='OBJ_ii')
	tmp2	<- subset(dfa, TEAM=='True' & OBJ=='OBJ_ii')	
	setkey(tmp2, central)
	set(tmp, NULL, 'SC_RND', factor(tmp[,as.character(SC_RND)], levels=tmp2[, SC_RND], labels=tmp2[, SC_RND]))
	set(tmp2, NULL, 'SC_RND', factor(tmp2[,as.character(SC_RND)], levels=tmp2[, SC_RND], labels=tmp2[, SC_RND]))
	setnames(tmp2, 'central', 'central_t')
	tmp2[, TEAM:=NULL]
	tmp2	<- merge( data.table(expand.grid(SC_RND=tmp[, unique(SC_RND)], TEAM=tmp[, unique(TEAM)])), tmp2, by='SC_RND')
	tmp		<- merge(tmp, subset(tmp2, select=c(SC_RND, TEAM, central_t)), by=c('SC_RND','TEAM'), all=1)	
	tmp2	<- subset(dfa, USED_GENES=='all' & TEAM!='True' & !grepl('(', TEAM,fixed=1) & OBJ=='OBJ_i', c(SC_RND, TEAM, central))
	setnames(tmp2, 'central', 'OBJ_i')
	set(tmp2, NULL, 'SC_RND', factor(tmp2[,as.character(SC_RND)], levels=tmp[, levels(SC_RND)], labels=tmp[, levels(SC_RND)]))
	tmp		<- merge(tmp, tmp2, by=c('SC_RND','TEAM'), all=1)
	
	ggplot( tmp, aes(x=SC_RND, group=TEAM, colour=gsub('using\n','',DATA_T2))) +
			geom_point(aes(y=central, shape=factor(OBJ_i, levels=c(-1,0,1), labels=c('declining','stable','increasing'))), position=position_dodge(width = 0.90)) + 
			geom_errorbar(aes(ymin=lower95, ymax=upper95), na.rm=TRUE, position=position_dodge(width = 0.90)) +
			geom_point(aes(y=central_t), colour='black', pch=18) +
			coord_cartesian(ylim=c(0,20)) +
			scale_colour_brewer(palette='Set1') +
			#scale_y_continuous(breaks=seq(0,3,0.5), minor_breaks=seq(0,3,0.1)) +
			#scale_size_manual(values = seq(1.5, by=0.5, length.out=4)) +
			scale_shape_manual(values = c(0, 1, 5)) +
			facet_wrap(~TEAM, ncol=3) +
			labs(x='\nPANGEA data set', y='Estimated and true % Incidence\n', colour='Estimates based on', shape='Estimated trend in incidence during intervention') +
			theme_bw() + 
			theme(legend.position='bottom', axis.text.x=element_text(size=5.5), panel.grid.major.y=element_line(colour='grey70', size=1), panel.grid.minor.y=element_line(colour='grey70', size=0.4))
	ggsave(file=paste(outdir,'/res_acrossTEAM_Primary_PcIncidence.pdf',sep=''), width=10, height=7)
}
##--------------------------------------------------------------------------------------------------------
##	olli 25.11.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.primary.sortedIncidenceRatio<- function(dfa, outdir)
{
	tmp		<- subset(dfa, USED_GENES=='all' & TEAM!='True' & !grepl('(', TEAM,fixed=1) & OBJ=='OBJ_iii')
	tmp2	<- subset(dfa, TEAM=='True' & OBJ=='OBJ_iii')	
	setkey(tmp2, central)
	set(tmp, NULL, 'SC_RND', factor(tmp[,as.character(SC_RND)], levels=tmp2[, SC_RND], labels=tmp2[, SC_RND]))
	set(tmp2, NULL, 'SC_RND', factor(tmp2[,as.character(SC_RND)], levels=tmp2[, SC_RND], labels=tmp2[, SC_RND]))
	setnames(tmp2, 'central', 'central_t')
	tmp2[, TEAM:=NULL]
	tmp2	<- merge( data.table(expand.grid(SC_RND=tmp[, unique(SC_RND)], TEAM=tmp[, unique(TEAM)])), tmp2, by='SC_RND')
	tmp		<- merge(tmp, subset(tmp2, select=c(SC_RND, TEAM, central_t)), by=c('SC_RND','TEAM'), all=1)	
	tmp2	<- subset(dfa, USED_GENES=='all' & TEAM!='True' & !grepl('(', TEAM,fixed=1) & OBJ=='OBJ_i', c(SC_RND, TEAM, central))
	setnames(tmp2, 'central', 'OBJ_i')
	set(tmp2, NULL, 'SC_RND', factor(tmp2[,as.character(SC_RND)], levels=tmp[, levels(SC_RND)], labels=tmp[, levels(SC_RND)]))
	tmp		<- merge(tmp, tmp2, by=c('SC_RND','TEAM'), all=1)
	
	ggplot( tmp, aes(x=SC_RND, group=TEAM, colour=gsub('using\n','',DATA_T2))) +
			geom_point(aes(y=central, shape=factor(OBJ_i, levels=c(-1,0,1), labels=c('declining','stable','increasing'))), position=position_dodge(width = 0.90)) + 
			geom_errorbar(aes(ymin=lower95, ymax=upper95), na.rm=TRUE, position=position_dodge(width = 0.90)) +
			geom_point(aes(y=central_t), colour='black', pch=18) +
			coord_cartesian(ylim=c(0,1.99)) +
			scale_colour_brewer(palette='Set1') +
			scale_y_continuous(breaks=seq(0,3,0.5), minor_breaks=seq(0,3,0.1)) +
			#scale_size_manual(values = seq(1.5, by=0.5, length.out=4)) +
			scale_shape_manual(values = c(0, 1, 5)) +
			facet_wrap(~TEAM, ncol=3) +
			labs(x='\nPANGEA data set', y='Estimated and true incidence ratio\n', colour='Estimates based on', shape='Estimated trend in incidence during intervention') +
			theme_bw() + 
			theme(legend.position='bottom', axis.text.x=element_text(size=5.5), panel.grid.major.y=element_line(colour='grey70', size=1), panel.grid.minor.y=element_line(colour='grey70', size=0.4))
	ggsave(file=paste(outdir,'/res_acrossTEAM_Primary_IncidenceRatio.pdf',sep=''), width=10, height=7)
	
	ggplot( tmp, aes(x=SC_RND, group=TEAM, colour=gsub('using\n','',DATA_T2))) +
			geom_point(aes(y=central, shape=factor(OBJ_i, levels=c(-1,0,1), labels=c('declining','stable','increasing'))), position=position_dodge(width = 0.90)) + 
			geom_errorbar(aes(ymin=lower95, ymax=upper95), na.rm=TRUE, position=position_dodge(width = 0.90)) +
			geom_point(aes(y=central_t), colour='black', pch=18) +
			coord_cartesian(ylim=c(1.99,8)) +
			scale_colour_brewer(palette='Set1') +
			scale_y_continuous(breaks=seq(0,8,2), minor_breaks=seq(0,8,0.5)) +
			#scale_size_manual(values = seq(1.5, by=0.5, length.out=4)) +
			scale_shape_manual(values = c(0, 1, 5)) +
			facet_wrap(~TEAM, ncol=3) +
			labs(x='\nPANGEA data set', y='Outliers in estimated incidence ratio\n', colour='Estimates based on', shape='Estimated trend in incidence during intervention') +
			theme_bw() + 
			theme(legend.position='bottom', axis.text.x=element_text(size=5.5), panel.grid.major.y=element_line(colour='grey70', size=1), panel.grid.minor.y=element_line(colour='grey70', size=0.4))
	ggsave(file=paste(outdir,'/res_acrossTEAM_Primary_IncidenceRatio_Outliers.pdf',sep=''), width=10, height=7)
}
##--------------------------------------------------------------------------------------------------------
##	olli 27.11.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.primary.AcutePercentvsSamplingCoverage<- function(dfa, outdir)
{
	tmp		<- subset(dfa, USED_GENES=='all' & !grepl('(', TEAM,fixed=1) & OBJ%in%c('OBJ_v','OBJ_vi'), c(SC_RND, TEAM, DATAT_L, central, lower95, upper95, DATA_T, AC_T, INT_T, IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D, OBJ))
	tmp2	<- subset(tmp, TEAM=='True', c(SC_RND, OBJ, central))
	setnames(tmp2, c('central'), c('central_t'))
	tmp		<- merge(tmp, tmp2, by=c('SC_RND','OBJ'))
	tmp		<- subset(tmp, TEAM!='True')
	ggplot(tmp, aes(x=SC_RND, y=abs(central-central_t),  group=OBJ)) +
			geom_point(aes(colour=factor(OBJ, levels=c('OBJ_v','OBJ_vi'), labels=c('just before intervention','after intervention'))), position=position_dodge(width = 0.6)) +
			geom_boxplot(aes(fill=factor(OBJ, levels=c('OBJ_v','OBJ_vi'), labels=c('just before intervention','after intervention'))), colour='black', outlier.shape=NA, alpha=0.1) +
			scale_y_continuous(expand=c(0,0), breaks=seq(0,30,10), minor_breaks=seq(0,30,5)) +
			scale_colour_brewer(palette='Set1') +
			scale_fill_brewer(palette='Set1') +
			facet_grid(TEAM~SMPL_C, scales='free_x', space='free_x') +
			theme_bw() + theme(legend.position='bottom', panel.grid.major.y=element_line(colour='grey70', size=0.4), panel.grid.minor.y=element_line(colour='grey70', size=0.4)) +
			labs(x="\nPANGEA data set", y='absolute error between\nestimated and true proportion of early transmissions\n', colour='Proportion early transmissions', fill='Proportion early transmissions')
	ggsave(file=paste(outdir,'/res_acrossTEAM_Primary_AcutePercent_SamplingCoverage.pdf',sep=''), width=10, height=7)
}
##--------------------------------------------------------------------------------------------------------
##	olli 26.11.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.primary.sortedAcutePercent<- function(dfa, outdir)
{
	tmp		<- subset(dfa, USED_GENES=='all' & TEAM!='True' & !grepl('(', TEAM,fixed=1) & OBJ=='OBJ_v')
	tmp2	<- subset(dfa, TEAM=='True' & OBJ=='OBJ_v')	
	setkey(tmp2, central)
	set(tmp, NULL, 'SC_RND', factor(tmp[,as.character(SC_RND)], levels=tmp2[, SC_RND], labels=tmp2[, SC_RND]))
	set(tmp2, NULL, 'SC_RND', factor(tmp2[,as.character(SC_RND)], levels=tmp2[, SC_RND], labels=tmp2[, SC_RND]))
	setnames(tmp2, 'central', 'central_t')
	tmp2[, TEAM:=NULL]
	tmp2	<- merge( data.table(expand.grid(SC_RND=tmp[, unique(SC_RND)], TEAM=tmp[, unique(TEAM)])), tmp2, by='SC_RND')
	tmp		<- merge(tmp, subset(tmp2, select=c(SC_RND, TEAM, DATAT_L, central_t)), by=c('SC_RND','TEAM','DATAT_L'), all=1)	
	tmp2	<- subset(dfa, USED_GENES=='all' & TEAM!='True' & !grepl('(', TEAM,fixed=1) & OBJ=='OBJ_iv', c(SC_RND, TEAM, central))
	setnames(tmp2, 'central', 'OBJ_iv')
	set(tmp2, NULL, 'SC_RND', factor(tmp2[,as.character(SC_RND)], levels=tmp[, levels(SC_RND)], labels=tmp[, levels(SC_RND)]))
	tmp		<- merge(tmp, tmp2, by=c('SC_RND','TEAM'), all=1)
	
	tmp[, SMPL_C_BEFORE:= NA_character_]
	set(tmp, tmp[, which(DATAT_L=='Village'  & !is.na(SMPL_C))], 'SMPL_C_BEFORE', '<1%')
	set(tmp, tmp[, which(DATAT_L=='Regional' & SMPL_C=='8%' & SMPL_M=='much')], 'SMPL_C_BEFORE', '4%')
	set(tmp, tmp[, which(DATAT_L=='Regional' & SMPL_C=='16%' & SMPL_M=='much')], 'SMPL_C_BEFORE', '8%')
	set(tmp, tmp[, which(DATAT_L=='Regional' & SMPL_M=='extreme')], 'SMPL_C_BEFORE', '1%')
	tmp2	<- tmp[, which(!is.na(SMPL_C))]
	set(tmp, tmp2, 'SMPL_C_BEFORE', tmp[tmp2, paste(SMPL_C_BEFORE,' just before intervention\n',SMPL_C,' after intervention',sep='')])
	
	ggplot( tmp, aes(x=SC_RND, group=TEAM, colour=gsub('using\n','',DATA_T2))) +
			geom_point(aes(y=central, shape=SMPL_C_BEFORE, size=SMPL_C_BEFORE), position=position_dodge(width = 0.90)) +
			geom_errorbar(aes(ymin=lower95, ymax=upper95), na.rm=TRUE, position=position_dodge(width = 0.90)) +
			geom_point(aes(y=central_t), colour='black', pch=18) +
			coord_cartesian(ylim=c(0,50)) +
			scale_colour_brewer(palette='Set1') +
			scale_y_continuous(breaks=seq(0,50,10), minor_breaks=seq(0,50,5)) +
			scale_size_manual(values = c(1.5, 4, 1.5, 1.5, 4)) +
			scale_shape_manual(values = c(0, 1, 5, 6, 7)) +
			facet_wrap(DATAT_L~TEAM, ncol=4, scale='free') +
			labs(x='\nPANGEA data set', y='Estimated and true % early transmissions\njust before the intervention\n', colour='Estimates based on', shape='Sampling coverage', size='Sampling coverage') +
			theme_bw() + 
			theme(legend.position='bottom', axis.text.x=element_text(size=8), panel.grid.major.y=element_line(colour='grey70', size=1), panel.grid.minor.y=element_line(colour='grey70', size=0.4)) +
			guides(shape=guide_legend(ncol=3), size=guide_legend(ncol=3))
	ggsave(file=paste(outdir,'/res_acrossTEAM_Primary_PCearlyjustbefore.pdf',sep=''), width=10, height=7)
		
	ggplot( tmp, aes(x=central_t, y=central, colour=gsub('using\n','',DATA_T2))) +
			geom_point(aes(shape=SMPL_C, size=SMPL_C)) +
			geom_errorbar(aes(ymin=lower95, ymax=upper95), na.rm=TRUE) +
			geom_abline(intercept=0, slope=1) +
			coord_cartesian(ylim=c(0,50)) +
			scale_colour_brewer(palette='Set1') +
			scale_size_manual(values = seq(1.5, by=1, length.out=4)) +
			scale_shape_manual(values = c(0, 1, 5, 6)) +			
			facet_grid(DATAT_L~TEAM) +
			theme_bw() +
			theme(legend.position='bottom', panel.grid.major=element_line(colour='grey70', size=0.4), panel.grid.minor=element_line(colour='grey70', size=0.4)) +
			labs(x='\nTrue % early transmissions\njust before the intervention', y='Estimated % early transmissions\njust before the intervention\n', colour='Estimates based on', shape='Sampling coverage', size='Sampling coverage') 
	ggsave(file=paste(outdir,'/res_acrossTEAM_Primary_PCearlyjustbefore2.pdf',sep=''), width=10, height=7)
	
	tmp		<- subset(dfa, USED_GENES=='all' & TEAM!='True' & !grepl('(', TEAM,fixed=1) & OBJ=='OBJ_vi')
	tmp2	<- subset(dfa, TEAM=='True' & OBJ=='OBJ_vi')	
	setkey(tmp2, central)
	set(tmp, NULL, 'SC_RND', factor(tmp[,as.character(SC_RND)], levels=tmp2[, SC_RND], labels=tmp2[, SC_RND]))
	set(tmp2, NULL, 'SC_RND', factor(tmp2[,as.character(SC_RND)], levels=tmp2[, SC_RND], labels=tmp2[, SC_RND]))
	setnames(tmp2, 'central', 'central_t')
	tmp2[, TEAM:=NULL]
	tmp2	<- merge( data.table(expand.grid(SC_RND=tmp[, unique(SC_RND)], TEAM=tmp[, unique(TEAM)])), tmp2, by='SC_RND')
	tmp		<- merge(tmp, subset(tmp2, select=c(SC_RND, TEAM, DATAT_L, central_t)), by=c('SC_RND','TEAM','DATAT_L'), all=1)	
	tmp2	<- subset(dfa, USED_GENES=='all' & TEAM!='True' & !grepl('(', TEAM,fixed=1) & OBJ=='OBJ_iv', c(SC_RND, TEAM, central))
	setnames(tmp2, 'central', 'OBJ_iv')
	set(tmp2, NULL, 'SC_RND', factor(tmp2[,as.character(SC_RND)], levels=tmp[, levels(SC_RND)], labels=tmp[, levels(SC_RND)]))
	tmp		<- merge(tmp, tmp2, by=c('SC_RND','TEAM'), all=1)
	
	tmp[, SMPL_C_BEFORE:= NA_character_]
	set(tmp, tmp[, which(DATAT_L=='Village'  & !is.na(SMPL_C))], 'SMPL_C_BEFORE', '<1%')
	set(tmp, tmp[, which(DATAT_L=='Regional' & SMPL_C=='8%' & SMPL_M=='much')], 'SMPL_C_BEFORE', '4%')
	set(tmp, tmp[, which(DATAT_L=='Regional' & SMPL_C=='16%' & SMPL_M=='much')], 'SMPL_C_BEFORE', '8%')
	set(tmp, tmp[, which(DATAT_L=='Regional' & SMPL_M=='extreme')], 'SMPL_C_BEFORE', '1%')
	tmp2	<- tmp[, which(!is.na(SMPL_C))]
	set(tmp, tmp2, 'SMPL_C_BEFORE', tmp[tmp2, paste(SMPL_C_BEFORE,' just before intervention\n',SMPL_C,' after intervention',sep='')])
	
	
	ggplot( tmp, aes(x=SC_RND, group=TEAM, colour=gsub('using\n','',DATA_T2))) +
			#geom_point(aes(y=central, shape=factor(OBJ_iv, levels=c(-1,0,1), labels=c('<15%','15-30%','>30%'))), position=position_dodge(width = 0.90)) +
			geom_point(aes(y=central, shape=SMPL_C_BEFORE, size=SMPL_C_BEFORE), position=position_dodge(width = 0.90)) +
			geom_errorbar(aes(ymin=lower95, ymax=upper95), na.rm=TRUE, position=position_dodge(width = 0.90)) +
			geom_point(aes(y=central_t), colour='black', pch=18) +
			coord_cartesian(ylim=c(0,50)) +
			scale_colour_brewer(palette='Set1') +
			scale_y_continuous(breaks=seq(0,50,10), minor_breaks=seq(0,50,5)) +
			scale_size_manual(values = c(1.5, 4, 1.5, 1.5, 4)) +
			scale_shape_manual(values = c(0, 1, 5, 6, 7)) +
			facet_wrap(DATAT_L~TEAM, ncol=4, scale='free') +
			labs(x='\nPANGEA data set', y='Estimated and true % early transmissions\nafter the intervention\n', colour='Estimates based on', shape='Sampling coverage', size='Sampling coverage') +
			theme_bw() + 
			theme(legend.position='bottom', axis.text.x=element_text(size=8), panel.grid.major.y=element_line(colour='grey70', size=1), panel.grid.minor.y=element_line(colour='grey70', size=0.4)) +
			guides(shape=guide_legend(ncol=3), size=guide_legend(ncol=3))
	ggsave(file=paste(outdir,'/res_acrossTEAM_Primary_PCearlyafter.pdf',sep=''), width=10, height=7)
	
}
##--------------------------------------------------------------------------------------------------------
##	olli 25.11.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.primary.CovariatesIncidenceRatio<- function(dfa, outdir)
{
	tmp		<- subset(dfa, USED_GENES=='all' & TEAM!='True' & !grepl('(', TEAM,fixed=1) & OBJ=='OBJ_iii')
	tmp2	<- subset(dfa, TEAM=='True' & OBJ=='OBJ_iii')	
	setkey(tmp2, central)
	set(tmp, NULL, 'SC_RND', factor(tmp[,as.character(SC_RND)], levels=tmp2[, SC_RND], labels=tmp2[, SC_RND]))
	set(tmp2, NULL, 'SC_RND', factor(tmp2[,as.character(SC_RND)], levels=tmp2[, SC_RND], labels=tmp2[, SC_RND]))
	setnames(tmp2, 'central', 'central_t')
	tmp2[, TEAM:=NULL]
	tmp2	<- merge( data.table(expand.grid(SC_RND=tmp[, unique(SC_RND)], TEAM=tmp[, unique(TEAM)])), tmp2, by='SC_RND')
	tmp		<- merge(tmp, subset(tmp2, select=c(SC_RND, TEAM, central_t)), by=c('SC_RND','TEAM'), all=1)	
	tmp2	<- subset(dfa, USED_GENES=='all' & TEAM!='True' & !grepl('(', TEAM,fixed=1) & OBJ=='OBJ_i', c(SC_RND, TEAM, central))
	setnames(tmp2, 'central', 'OBJ_i')
	set(tmp2, NULL, 'SC_RND', factor(tmp2[,as.character(SC_RND)], levels=tmp[, levels(SC_RND)], labels=tmp[, levels(SC_RND)]))
	tmp		<- merge(tmp, tmp2, by=c('SC_RND','TEAM'), all=1)
	tmp		<- subset(tmp, !is.na(central), c(SC_RND, TEAM, DATAT_L, central, lower95, upper95, DATA_T, AC_T, INT_T, IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D, central_t))
	
	tmp[, RES:= central-central_t]
	tmp[, OUTLIER:= central>2]
	dfc		<- subset(tmp, !OUTLIER)[, list( COR=cor(central_t, central), COR_LOG=cor(log(central_t), log(central)), 											 
											 MAE=mean( abs(central-central_t) ), MAE_L=mean( abs(log(central)-log(central_t)) ), 
											 MSE=mean( (central-central_t)^2 ), MSE_L=mean( (log(central)-log(central_t))^2 ),
											 ARME= mean( central-central_t ), 
											 HAME_R=1/mean(1/(central/central_t)), ARME_L=mean(log(central)-log(central_t))), by='TEAM']
	
	require(gamlss)
	bw.AIC			<- vector('list', tmp[, length(unique(TEAM))])
	names(bw.AIC)	<- tmp[, unique(TEAM)]
	bw.BIC			<- vector('list', tmp[, length(unique(TEAM))])
	names(bw.BIC)	<- tmp[, unique(TEAM)]
	for(x in names(bw.AIC))
	{
		cat('\nprocess',x)
		df			<- subset(tmp, TEAM==x)	
		set(df, NULL, c('lower95','upper95'), NULL)
		if(x!='Cambridge')
			mnoa	<- gamlss(RES~DATAT_L+DATA_T+AC_T+INT_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D-1, data=df, family=NO, trace=FALSE)
		if(x=='Cambridge')
			mnoa	<- gamlss(RES~AC_T+INT_T+IMPRT+SMPL_N-1, data=df, family=NO, trace=FALSE)
		bw.AIC[[x]]	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE)		
		bw.BIC[[x]] <- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE, k=log(nrow(df)))				
	}
	dfcoeff	<- tmp[, 	{
							z			<- names(coef(bw.AIC[[as.character(TEAM)]]))
							z2			<- names(coef(bw.BIC[[as.character(TEAM)]]))
							length(z)	<- max(length(z),length(z2))
							length(z2)	<- max(length(z),length(z2))
							list( BWAIC=z, BWBIC=z2 )	
						}, by='TEAM']	
}
##--------------------------------------------------------------------------------------------------------
##	olli 26.11.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.primary.CovariatesAcutePercent<- function(dfa, outdir)
{
	tmp		<- subset(dfa, USED_GENES=='all' & TEAM!='True' & !grepl('(', TEAM,fixed=1) & OBJ=='OBJ_v')
	tmp2	<- subset(dfa, TEAM=='True' & OBJ=='OBJ_v')	
	setkey(tmp2, central)
	set(tmp, NULL, 'SC_RND', factor(tmp[,as.character(SC_RND)], levels=tmp2[, SC_RND], labels=tmp2[, SC_RND]))
	set(tmp2, NULL, 'SC_RND', factor(tmp2[,as.character(SC_RND)], levels=tmp2[, SC_RND], labels=tmp2[, SC_RND]))
	setnames(tmp2, 'central', 'central_t')
	tmp2[, TEAM:=NULL]
	tmp2	<- merge( data.table(expand.grid(SC_RND=tmp[, unique(SC_RND)], TEAM=tmp[, unique(TEAM)])), tmp2, by='SC_RND')
	tmp		<- merge(tmp, subset(tmp2, select=c(SC_RND, TEAM, DATAT_L, central_t)), by=c('SC_RND','DATAT_L','TEAM'), all=1)	
	tmp2	<- subset(dfa, USED_GENES=='all' & TEAM!='True' & !grepl('(', TEAM,fixed=1) & OBJ=='OBJ_iv', c(SC_RND, TEAM, central))
	setnames(tmp2, 'central', 'OBJ_iv')
	set(tmp2, NULL, 'SC_RND', factor(tmp2[,as.character(SC_RND)], levels=tmp[, levels(SC_RND)], labels=tmp[, levels(SC_RND)]))
	tmp		<- merge(tmp, tmp2, by=c('SC_RND','TEAM'), all=1)
	tmp		<- subset(tmp, !is.na(central), c(SC_RND, TEAM, DATAT_L, central, lower95, upper95, DATA_T, AC_T, INT_T, IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D, central_t))
	
	tmp[, RES:= central-central_t]
	tmp[, OUTLIER:= central>60]
	dfcb	<- subset(tmp, !OUTLIER)[, list( COR=cor(central_t, central),  											 
					MAE=mean( abs(central-central_t) ), 
					MSE=mean( (central-central_t)^2 ),
					ARME=mean(central-central_t)), by=c('DATAT_L','TEAM')]
	dfcbo	<- subset(tmp, !OUTLIER)[, list( COR=cor(central_t, central),  											 
					MAE=mean( abs(central-central_t) ), 
					MSE=mean( (central-central_t)^2 ),
					ARME=mean(central-central_t)), by=c('TEAM')]
	
	require(gamlss)
	bw.AIC			<- vector('list', tmp[, length(unique(TEAM))])
	names(bw.AIC)	<- tmp[, unique(TEAM)]
	bw.BIC			<- vector('list', tmp[, length(unique(TEAM))])
	names(bw.BIC)	<- tmp[, unique(TEAM)]
	for(x in names(bw.AIC))
	{
		cat('\nprocess',x)
		df			<- subset(tmp, TEAM==x)	
		set(df, NULL, c('lower95','upper95'), NULL)		
		#stepGAIC.VR(mnoa, direction='forward', what='mu', trace=1, scope=~DATAT_L+DATA_T+AC_T+INT_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D)
		if(x!='Cambridge')	# exclude IMPRT SMPL_C
			mnoa	<- gamlss(RES~DATAT_L+DATA_T+AC_T+INT_T+SMPL_N+SMPL_M+SMPL_D-1, data=df, family=NO, trace=FALSE)			
		if(x=='Cambridge')
			mnoa	<- gamlss(RES~AC_T+INT_T+IMPRT+SMPL_N-1, data=df, family=NO, trace=FALSE)
		bw.AIC[[x]]	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE)		
		bw.BIC[[x]] <- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE, k=log(nrow(df)))				
	}
	dfcoeffb	<- tmp[, 	{
				z			<- summary(bw.AIC[[as.character(TEAM)]])[names(coef(bw.AIC[[as.character(TEAM)]])),'Pr(>|t|)']
				z2			<- summary(bw.BIC[[as.character(TEAM)]])[names(coef(bw.BIC[[as.character(TEAM)]])),'Pr(>|t|)']
				z			<- sort(z)
				z2			<- sort(z2)
				length(z)	<- max(length(z),length(z2))
				length(z2)	<- max(length(z),length(z2))
				list( BWAIC=names(z), BWAICp=z, BWBIC=names(z2), BWBICp=z2 )	
			}, by='TEAM']	
	
	
	tmp		<- subset(dfa, USED_GENES=='all' & TEAM!='True' & !grepl('(', TEAM,fixed=1) & OBJ=='OBJ_vi')
	tmp2	<- subset(dfa, TEAM=='True' & OBJ=='OBJ_vi')	
	setkey(tmp2, central)
	set(tmp, NULL, 'SC_RND', factor(tmp[,as.character(SC_RND)], levels=tmp2[, SC_RND], labels=tmp2[, SC_RND]))
	set(tmp2, NULL, 'SC_RND', factor(tmp2[,as.character(SC_RND)], levels=tmp2[, SC_RND], labels=tmp2[, SC_RND]))
	setnames(tmp2, 'central', 'central_t')
	tmp2[, TEAM:=NULL]
	tmp2	<- merge( data.table(expand.grid(SC_RND=tmp[, unique(SC_RND)], TEAM=tmp[, unique(TEAM)])), tmp2, by='SC_RND')
	tmp		<- merge(tmp, subset(tmp2, select=c(SC_RND, TEAM, DATAT_L, central_t)), by=c('SC_RND','DATAT_L','TEAM'), all=1)		
	tmp		<- subset(tmp, !is.na(central), c(SC_RND, TEAM, DATAT_L, central, lower95, upper95, DATA_T, AC_T, INT_T, IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D, central_t))
	
	tmp[, RES:= central-central_t]
	tmp[, OUTLIER:= central>60]
	dfca	<- subset(tmp, !OUTLIER)[, list( COR=cor(central_t, central),  											 
					MAE=mean( abs(central-central_t) ), 
					MSE=mean( (central-central_t)^2 ),
					ARME=mean(central-central_t)), by=c('DATAT_L','TEAM')]
	dfcao	<- subset(tmp, !OUTLIER)[, list( COR=cor(central_t, central),  											 
					MAE=mean( abs(central-central_t) ), 
					MSE=mean( (central-central_t)^2 ),
					ARME=mean(central-central_t)), by=c('TEAM')]
	
	bw.AICR			<- vector('list', tmp[, length(unique(TEAM))])
	names(bw.AICR)	<- tmp[, unique(TEAM)]
	bw.BICR			<- vector('list', tmp[, length(unique(TEAM))])
	names(bw.BICR)	<- tmp[, unique(TEAM)]
	bw.AICV			<- vector('list', tmp[, length(unique(TEAM))])
	names(bw.AICV)	<- tmp[, unique(TEAM)]
	bw.BICV			<- vector('list', tmp[, length(unique(TEAM))])
	names(bw.BICV)	<- tmp[, unique(TEAM)]	
	for(x in names(bw.AICR))
	{
		cat('\nprocess',x)
		df			<- subset(tmp, DATAT_L=='Village' & TEAM==x)	
		set(df, NULL, c('lower95','upper95'), NULL)		
		#stepGAIC.VR(mnoa, direction='forward', what='mu', trace=1, scope=~DATAT_L+DATA_T+AC_T+INT_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D)
		if(x!='Cambridge')	# exclude IMPRT SMPL_C
			mnoa	<- gamlss(RES~DATA_T+AC_T+INT_T+SMPL_N+SMPL_D-1, data=df, family=NO, trace=FALSE)			
		if(x=='Cambridge')
			mnoa	<- gamlss(RES~AC_T+INT_T+IMPRT+SMPL_N-1, data=df, family=NO, trace=FALSE)
		bw.AICV[[x]]<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE)		
		bw.BICV[[x]]<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE, k=log(nrow(df)))
		
		df			<- subset(tmp, DATAT_L=='Regional' & TEAM==x)	
		set(df, NULL, c('lower95','upper95'), NULL)		
		#stepGAIC.VR(mnoa, direction='forward', what='mu', trace=1, scope=~DATAT_L+DATA_T+AC_T+INT_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D)
		if(x!='Cambridge')	# exclude IMPRT SMPL_C
			mnoa	<- gamlss(RES~DATA_T+AC_T+INT_T+SMPL_N+SMPL_M+SMPL_D-1, data=df, family=NO, trace=FALSE)			
		if(x=='Cambridge')
			mnoa	<- gamlss(RES~AC_T+INT_T+IMPRT+SMPL_N-1, data=df, family=NO, trace=FALSE)
		bw.AICR[[x]]<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE)		
		bw.BICR[[x]]<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE, k=log(nrow(df)))
	}
	dfcoeffa	<- tmp[, 	{
				cat(TEAM, DATAT_L)
				z <- z2	<- numeric(0)
				if(DATAT_L=='Village' & !is.null(names(coef(bw.AICV[[as.character(TEAM)]]))))
				{
					z		<- summary(bw.AICV[[as.character(TEAM)]])[,'Pr(>|t|)']
					z		<- z[ intersect(names(z), names(coef(bw.AICV[[as.character(TEAM)]]))) ]
				}					
				if(DATAT_L=='Village' & !is.null(names(coef(bw.BICV[[as.character(TEAM)]]))))
				{
					z2		<- summary(bw.BICV[[as.character(TEAM)]])[,'Pr(>|t|)']
					z2		<- z2[ intersect(names(z2), names(coef(bw.BICV[[as.character(TEAM)]]))) ]					
				}									
				if(DATAT_L=='Regional'& !is.null(names(coef(bw.AICR[[as.character(TEAM)]]))))
				{
					z		<- summary(bw.AICR[[as.character(TEAM)]])[,'Pr(>|t|)']
					z		<- z[ intersect(names(z), names(coef(bw.AICR[[as.character(TEAM)]]))) ]
				}					
				if(DATAT_L=='Regional'& !is.null(names(coef(bw.BICR[[as.character(TEAM)]]))))
				{
					z2		<- summary(bw.BICR[[as.character(TEAM)]])[,'Pr(>|t|)']
					z2		<- z2[ intersect(names(z2), names(coef(bw.BICR[[as.character(TEAM)]]))) ]
				}									
				z			<- sort(z)
				z2			<- sort(z2)
				length(z)	<- max(length(z),length(z2))
				length(z2)	<- max(length(z),length(z2))
				list( BWAIC=names(z), BWAICp=z, BWBIC=names(z2), BWBICp=z2 )	
			}, by=c('TEAM','DATAT_L')]	
	dfcoeffa	<- subset(dfcoeffa, BWAICp<0.05 | BWBIC<.05)
}
##--------------------------------------------------------------------------------------------------------
##	olli 02.12.15
##--------------------------------------------------------------------------------------------------------
grubbs.flag <- function(x) 
{
	outliers 	<- NULL
	test 		<- x
	gr 			<- grubbs.test(test)
	pv 			<- gr$p.value
	# throw an error if there are too few values for the Grubb's test
	if (length(test) < 3 ) 
		stop("Grubb's test requires > 2 input values")
	while(pv < 0.05) 
	{		
		if( grepl('lowest',gr$alternative) )
		{
			outliers	<- c(outliers, min(test))
			test		<- test[-which.min(test)]
		}
		if( grepl('highest',gr$alternative) )
		{
			outliers	<- c(outliers, max(test))
			test		<- test[-which.max(test)]
		}		
		# stop if all but two values are flagged as outliers
		if (length(test) < 3 ) 
		{
			warning("All but two values flagged as outliers")
			break
		}
		gr <- grubbs.test(test)
		pv <- gr$p.value
	}
	return(data.frame(X=x, OUTLIER=(x %in% outliers)))
}
##--------------------------------------------------------------------------------------------------------
##	olli 05.12.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.secondary.TreeSeq<- function(dfa, outdir)
{
	require(exactRankTests)
	
	dfo		<- subset(dfa, OBJ%in%c('OBJ_ii','OBJ_iii','OBJ_v','OBJ_vi') & USED_GENES=='all' & TEAM!='True' & !grepl('(', TEAM,fixed=1), c(SC_RND, TEAM, DATAT_L, central, AC_T, INT_T,DATA_T,IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D, OBJ, Pr_Phy, Pr_Seq))
	dfo		<- rbind( subset(dfo, Pr_Phy==1), subset(dfo, Pr_Seq==1) ) 	
	tmp		<- dcast.data.table( subset(dfo, DATAT_L=='Village'), TEAM+DATAT_L+AC_T+INT_T+SMPL_C+SMPL_M+SMPL_D+OBJ~DATA_T, value.var='SC_RND')	
	dfo		<- dcast.data.table(subset(dfo, DATAT_L=='Regional'), TEAM+DATAT_L+AC_T+INT_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D+OBJ~DATA_T, value.var='SC_RND')
	dfo		<- rbind(dfo, tmp, use.names=TRUE, fill=TRUE)
	dfo		<- subset(dfo, !is.na(seq) & !is.na(phy))
	
	tmp		<- subset(dfa, OBJ%in%c('OBJ_ii','OBJ_iii','OBJ_v','OBJ_vi') & USED_GENES=='all' & TEAM!='True' & !grepl('(', TEAM,fixed=1), c(SC_RND, TEAM, OBJ, central, lower95, upper95))
	setnames(tmp, c('SC_RND','central','lower95','upper95'), c('seq','SEQ_central','SEQ_lower95','SEQ_upper95'))
	dfo		<- merge(dfo, tmp, by=c('TEAM','OBJ','seq'), all.x=1)
	tmp		<- subset(dfa, OBJ%in%c('OBJ_ii','OBJ_iii','OBJ_v','OBJ_vi') & USED_GENES=='all' & TEAM!='True' & !grepl('(', TEAM,fixed=1), c(SC_RND, TEAM, OBJ, central, lower95, upper95))
	setnames(tmp, c('SC_RND','central','lower95','upper95'), c('phy','PHY_central','PHY_lower95','PHY_upper95'))
	dfo		<- merge(dfo, tmp, by=c('TEAM','OBJ','phy'), all.x=1)
	tmp		<- subset(dfa, OBJ%in%c('OBJ_ii','OBJ_iii','OBJ_v','OBJ_vi') & TEAM=='True', c(SC_RND, OBJ, central))	
	setnames(tmp, c('SC_RND','central'), c('seq','SEQ_True'))
	dfo		<- merge(dfo, tmp, by=c('OBJ','seq'), all.x=1)
	tmp		<- subset(dfa, OBJ%in%c('OBJ_ii','OBJ_iii','OBJ_v','OBJ_vi') & TEAM=='True', c(SC_RND, OBJ, central))	
	setnames(tmp, c('SC_RND','central'), c('phy','PHY_True'))
	dfo		<- merge(dfo, tmp, by=c('OBJ','phy'), all.x=1)
	#	set residuals for incidence and incidence reduction
	tmp		<- dfo[, which(OBJ%in%c('OBJ_ii','OBJ_iii','OBJ_v','OBJ_vi'))]
	set(dfo, tmp, 'SEQ_central', dfo[tmp, log(SEQ_central)-log(SEQ_True)])
	set(dfo, tmp, 'SEQ_lower95', dfo[tmp, log(SEQ_lower95)-log(SEQ_True)])
	set(dfo, tmp, 'SEQ_upper95', dfo[tmp, log(SEQ_upper95)-log(SEQ_True)])
	set(dfo, tmp, 'PHY_central', dfo[tmp, log(PHY_central)-log(PHY_True)])
	set(dfo, tmp, 'PHY_lower95', dfo[tmp, log(PHY_lower95)-log(PHY_True)])
	set(dfo, tmp, 'PHY_upper95', dfo[tmp, log(PHY_upper95)-log(PHY_True)])
	#	set residuals for 
	if(0)
	{
		tmp		<- dfo[, which(OBJ%in%c('OBJ_v','OBJ_vi'))]
		set(dfo, tmp, 'SEQ_central', dfo[tmp, SEQ_central-SEQ_True])
		set(dfo, tmp, 'SEQ_lower95', dfo[tmp, SEQ_lower95-SEQ_True])
		set(dfo, tmp, 'SEQ_upper95', dfo[tmp, SEQ_upper95-SEQ_True])
		set(dfo, tmp, 'PHY_central', dfo[tmp, PHY_central-PHY_True])
		set(dfo, tmp, 'PHY_lower95', dfo[tmp, PHY_lower95-PHY_True])
		set(dfo, tmp, 'PHY_upper95', dfo[tmp, PHY_upper95-PHY_True])		
	}
	
	set(dfo, NULL, 'OBJ', dfo[, factor(OBJ, levels=c('OBJ_ii','OBJ_iii','OBJ_v','OBJ_vi'), labels=c('Incidence\nafter intervention','Incidence reduction\nduring intervention','Proportion of early transmissions\njust before intervention','Proportion of early transmissions\nafter intervention'))])
	
	ggplot(dfo, aes(x=seq, y=SEQ_central-PHY_central, ymin=SEQ_lower95-PHY_central, ymax=SEQ_upper95-PHY_central, colour=(SEQ_lower95>PHY_central | SEQ_upper95<PHY_central) )) +
			geom_point() + geom_errorbar() + geom_boxplot(aes(group=OBJ), outlier.shape=NA, fill='transparent', colour='black') +
			scale_colour_brewer(palette='Dark2', guide=FALSE) +
			scale_y_continuous(breaks=seq(-4,4,2), minor_breaks=seq(-4,4,0.5)) +
			facet_grid(TEAM~OBJ) + labs(x='\nPANGEA data set\n(top: data set containing sequences\nbottom: paired data set containing trees)',y='pairwise difference in error of estimates\nobtained from sequences and true trees\n') +
			theme_bw() +
			theme( panel.grid.major.y=element_line(colour='grey70', size=0.4), panel.grid.minor.y=element_line(colour='grey70', size=0.4))
	ggsave(file=paste(outdir,'/res_acrossTEAM_Secondary_TreeSeq.pdf',sep=''), width=10, height=10)
	#	
	
	dfo[,{
				z	<- t.test(SEQ_central, PHY_central, paired=TRUE)
				list(TYPE='TTEST',NV=length(SEQ_central), PV=z$p.value, CIL=z$conf.int[1], CIU=z$conf.int[2])
			}, by=c('TEAM')]
	dfo[,{
				z	<- perm.test(1000*SEQ_central,1000*PHY_central,paired = TRUE,exact = TRUE)				
				list(TYPE='EXACTRANK',PV=z$p.value)
			}, by=c('TEAM')]
	#	estimates on trees are sig different from those on sequences *only for Imperial*	
}
##--------------------------------------------------------------------------------------------------------
##	olli 07.12.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.secondary.Pol<- function(dfa, outdir)
{
	dfpol	<- subset(dfa, OBJ%in%c('OBJ_ii','OBJ_iii','OBJ_v','OBJ_vi') & USED_GENES=='pol' & TEAM!='True' & !grepl('(', TEAM,fixed=1), c(SC_RND, TEAM, DATAT_L, central, lower95, upper95, DATA_T,IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D, OBJ))
	setnames(dfpol, c('central','lower95','upper95'), c('POL_central','POL_lower95','POL_upper95'))
	tmp		<- merge(subset(dfa, USED_GENES=='all'), subset(dfpol, select=c('SC_RND','TEAM','OBJ')), by=c('SC_RND','TEAM','OBJ'))
	setnames(tmp, c('central','lower95','upper95'), c('ALL_central','ALL_lower95','ALL_upper95'))
	dfpol	<- merge(dfpol, subset(tmp, select=c('SC_RND','TEAM','OBJ','ALL_central','ALL_lower95','ALL_upper95')), by=c('SC_RND','TEAM','OBJ'))
	
	tmp		<- subset(dfa, TEAM=='True', c(SC_RND, TEAM, DATAT_L, central, DATA_T,IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D, OBJ))
	setnames(tmp, c('central'), c('TRUE_central'))
	dfpol	<- merge(dfpol, subset(tmp, select=c('SC_RND','OBJ','TRUE_central')), by=c('SC_RND','OBJ'))
	set(dfpol, NULL, 'POL_central', dfpol[, log(POL_central)-log(TRUE_central)] )
	set(dfpol, NULL, 'POL_upper95', dfpol[, log(POL_upper95)-log(TRUE_central)] )
	set(dfpol, NULL, 'POL_lower95', dfpol[, log(POL_lower95)-log(TRUE_central)] )
	set(dfpol, NULL, 'ALL_central', dfpol[, log(ALL_central)-log(TRUE_central)] )
	set(dfpol, NULL, 'ALL_upper95', dfpol[, log(ALL_upper95)-log(TRUE_central)] )
	set(dfpol, NULL, 'ALL_lower95', dfpol[, log(ALL_lower95)-log(TRUE_central)] )
}
##--------------------------------------------------------------------------------------------------------
##	olli 06.12.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.secondary.Errors.v151206<- function(dfa, outdir)
{
	require(plsdepot)
	require(pls)
	require(outliers)
	
	dfo		<- subset(dfa, OBJ%in%c('OBJ_ii','OBJ_iii','OBJ_v','OBJ_vi') & USED_GENES=='all' & TEAM!='True' & !grepl('(', TEAM,fixed=1), c(SC_RND, TEAM, DATAT_L, central, DATA_T,IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D, OBJ))
	dfo		<- dcast.data.table(dfo, SC_RND+TEAM+DATAT_L+DATA_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D~OBJ, value.var='central')			
	tmp		<- subset(dfa, OBJ%in%c('OBJ_ii','OBJ_iii','OBJ_v','OBJ_vi') & TEAM=='True', c(SC_RND, TEAM, DATAT_L, central, DATA_T,IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D, OBJ))
	set(tmp, NULL, 'OBJ', tmp[, paste(OBJ,'_t',sep='')])
	tmp		<- dcast.data.table(tmp, SC_RND~OBJ, value.var='central')
	dfo		<- merge(dfo, tmp, by='SC_RND')	
	set(dfo, NULL, 'DATAT_L', dfo[, as.numeric(factor(DATAT_L, levels=c("Regional","Village"), labels=c("Regional","Village")))])
	set(dfo, NULL, 'SMPL_M', dfo[, as.numeric(SMPL_M)])
	set(dfo, NULL, 'SMPL_D', dfo[, as.numeric(SMPL_D)])	
	set(dfo, NULL, 'DATA_T', dfo[, as.numeric(DATA_T)])
	if(1)
	{
		set(dfo, dfo[, which(IMPRT!='20%')], 'IMPRT', '<=5%')	
		set(dfo, NULL, 'IMPRT', dfo[, as.numeric(factor(as.character(IMPRT), levels=c('<=5%','20%'), labels=c('<=5%','20%')))])	
		set(dfo, dfo[, which(SMPL_C%in%c('8%','30%'))], 'SMPL_C', '1x')
		set(dfo, dfo[, which(SMPL_C%in%c('16%','60%'))], 'SMPL_C', '2x')
		set(dfo, NULL, 'SMPL_C', dfo[, as.numeric(factor(as.character(SMPL_C), levels=c('1x','2x'), labels=c('1x','2x')))])		
	}
	if(0)
	{
		set(dfo, NULL, 'SMPL_C', dfo[, as.numeric(gsub('%','',as.character(SMPL_C)))/100])
		set(dfo, NULL, 'IMPRT', dfo[, as.numeric(gsub('%','',as.character(IMPRT)))/100])		
	}
	setnames(dfo, c('OBJ_ii_t','OBJ_iii_t','OBJ_v_t','OBJ_vi_t'), c('INC_t','INCR_t','ACS_t','ACE_t'))
	#dfo[, R_ii_1:= OBJ_ii-INC_t]
	dfo[, R_ii:= log(OBJ_ii)-log(INC_t)]
	#dfo[, R_iii_1:= OBJ_iii-INCR_t]
	dfo[, R_iii:= log(OBJ_iii)-log(INCR_t)]	
	dfo[, R_v:= OBJ_v-ACS_t]
	#dfo[, R_v_2:= log(OBJ_v)-log(ACS_t)]
	dfo[, R_vi:= OBJ_vi-ACE_t]
	#dfo[, R_vi_2:= log(OBJ_vi)-log(ACE_t)]
	dfo			<- subset(dfo, select=c(SC_RND, TEAM, DATAT_L, DATA_T,IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D, INC_t, INCR_t, ACS_t, ACE_t, R_ii, R_iii, R_v, R_vi))
	dfo			<- melt(dfo, measure.vars=c('R_ii','R_iii','R_v','R_vi'), variable.name='OBJ', value.name='RESID')
	set(dfo, NULL, 'OBJ', dfo[, factor(OBJ, levels=c('R_ii','R_iii','R_v','R_vi'), labels=c('Incidence\nafter intervention','Incidence reduction\nduring intervention','Proportion of early transmissions\njust before intervention','Proportion of early transmissions\nafter intervention'))])
	dfo			<- subset(dfo, !is.na(RESID))
	#	restrict to not Cambridge
	dfo				<- subset(dfo, TEAM!='Cambridge')
	
	#
	#	do PLS on errors
	#	
	plsms2			<- vector('list', dfo[, length(unique(TEAM))])
	plsms3			<- vector('list', dfo[, length(unique(TEAM))])
	plsms5			<- vector('list', dfo[, length(unique(TEAM))])
	plsms6			<- vector('list', dfo[, length(unique(TEAM))])
	names(plsms2)	<- names(plsms3)	<- names(plsms5)	<- names(plsms6)	<- dfo[, unique(TEAM)]
	for(x in names(plsms2))
	{
		df			<- subset(dfo, TEAM==x & OBJ=='Incidence\nafter intervention')	
		plsms2[[x]]	<- plsr(RESID~DATAT_L+IMPRT+SMPL_C+SMPL_M+INC_t, data=df, validation='LOO')
		df			<- subset(dfo, TEAM==x & OBJ=='Incidence reduction\nduring intervention')	
		plsms3[[x]]	<- plsr(RESID~DATAT_L+IMPRT+SMPL_C+SMPL_M+INCR_t, data=df, validation='LOO')
		df			<- subset(dfo, TEAM==x & OBJ=='Proportion of early transmissions\njust before intervention')	
		plsms5[[x]]	<- plsr(RESID~DATAT_L+IMPRT+SMPL_C+SMPL_M+ACS_t, data=df, validation='LOO')
		df			<- subset(dfo, TEAM==x & OBJ=='Proportion of early transmissions\nafter intervention')	
		plsms6[[x]]	<- plsr(RESID~DATAT_L+IMPRT+SMPL_C+SMPL_M+ACE_t, data=df, validation='LOO')
	}		
	# get variance explained
	dfr		<- as.data.table(melt(sapply(plsms2, function(x) drop((R2(x, intercept=FALSE, estimate='train'))$val) ), varnames=c('LATENT','TEAM'), value.name='R2'))
	dfr[, OBJ:='Incidence\nafter intervention']
	tmp		<- as.data.table(melt(sapply(plsms3, function(x) drop((R2(x, intercept=FALSE, estimate='train'))$val) ), varnames=c('LATENT','TEAM'), value.name='R2'))
	tmp[, OBJ:='Incidence reduction\nduring intervention']
	dfr		<- rbind(dfr, tmp)
	tmp		<- as.data.table(melt(sapply(plsms5, function(x) drop((R2(x, intercept=FALSE, estimate='train'))$val) ), varnames=c('LATENT','TEAM'), value.name='R2'))
	tmp[, OBJ:='Proportion of early transmissions\njust before intervention']
	dfr		<- rbind(dfr, tmp)	
	tmp		<- as.data.table(melt(sapply(plsms6, function(x) drop((R2(x, intercept=FALSE, estimate='train'))$val) ), varnames=c('LATENT','TEAM'), value.name='R2'))
	tmp[, OBJ:='Proportion of early transmissions\nafter intervention']
	dfr		<- rbind(dfr, tmp)
	set(dfr, NULL, 'LATENT', dfr[, paste('LV',gsub(' comps','',LATENT),sep='')])
	dfr		<- merge(dfr, dfr[, list(LATENT=LATENT, R2each= c(R2[1],diff(R2))), by=c('TEAM','OBJ')], by=c('TEAM','OBJ','LATENT'))
	#oh wow! the variance explained by the various latent variables is quite different from each other!!	
	dfr		<- dfr[order(TEAM, OBJ, -R2each)]
	dfr		<- merge(dfr, dfr[, list(LATENT=LATENT, LATENTO= paste('lv',seq_along(LATENT),sep='')), by=c('TEAM','OBJ')], by=c('TEAM','OBJ','LATENT'))
	set(dfr, NULL, 'LATENTO', dfr[, factor(LATENTO, levels=paste('lv',1:10,sep=''), labels=paste('lv',1:10,sep=''))])
	
	#collect key X variables with high latent variable weights
	dfl		<- do.call('rbind',lapply(paste('Comp',1:5), function(x){
						do.call('rbind', 	list(	data.table(TEAM=tmp[, as.character(unique(TEAM))])[, list( X=names(plsms2[[TEAM]]$loading.weights[, 'Comp 1']), LATENT=x, LOAD=plsms2[[TEAM]]$loading.weights[, x], OBJ='Incidence\nafter intervention'), by='TEAM'],
										data.table(TEAM=tmp[, as.character(unique(TEAM))])[, list( X=names(plsms3[[TEAM]]$loading.weights[, 'Comp 1']), LATENT=x, LOAD=plsms3[[TEAM]]$loading.weights[, x], OBJ='Incidence reduction\nduring intervention'), by='TEAM'],
										data.table(TEAM=tmp[, as.character(unique(TEAM))])[, list( X=names(plsms5[[TEAM]]$loading.weights[, 'Comp 1']), LATENT=x, LOAD=plsms5[[TEAM]]$loading.weights[, x], OBJ='Proportion of early transmissions\njust before intervention'), by='TEAM'],
										data.table(TEAM=tmp[, as.character(unique(TEAM))])[, list( X=names(plsms6[[TEAM]]$loading.weights[, 'Comp 1']), LATENT=x, LOAD=plsms6[[TEAM]]$loading.weights[, x], OBJ='Proportion of early transmissions\nafter intervention'), by='TEAM'] 	))						
					}))
	set(dfl, NULL, 'LATENT', dfl[, gsub('Comp ','LV',LATENT)])		
	dfl		<- subset(dfl, !is.nan(LOAD))
	# use variable influence projection (Wold et al 1993)
	# https://books.google.co.uk/books?id=QhHdGt8TG80C&pg=PA2&lpg=PA2&dq=PLS+contribution+of+each+variable&source=bl&ots=vWaqNtCTYz&sig=RT9STQ3SzXk1tU2ZNYlRycgxIQ8&hl=en&sa=X&ved=0ahUKEwiH1eHmirvJAhXK7hoKHcBXC_YQ6AEILzAC#v=onepage&q=PLS%20contribution%20of%20each%20variable&f=false
	dfl		<- merge(dfl, dfl[, list(X=X, LOADstd=LOAD^2/sum(LOAD^2)), by=c('TEAM','OBJ','LATENT')], by=c('TEAM','OBJ','LATENT','X'))
	dfl		<- merge(dfl, dfr, by=c('TEAM','OBJ','LATENT'))
	set(dfl, NULL, 'X', dfl[, factor(X, 	levels=c("INC_t","INCR_t","ACS_t","ACE_t","DATAT_L","DATA_T","IMPRT","SMPL_N","SMPL_C","SMPL_M","SMPL_D"),
							labels=c('True % incidence','True incidence ratio','True % early transmissions just before intervention','True % early transmissions after intervention','Simulation model','Data provided', 'Frequency of viral introductions','Sequences (#)','Sequence coverage','Proportion of sequences from after intervention start','Sampling duration after intervention start'))])
	set(dfl, NULL, 'OBJ', dfl[, factor(OBJ, levels=c('Incidence\nafter intervention','Incidence reduction\nduring intervention','Proportion of early transmissions\njust before intervention','Proportion of early transmissions\nafter intervention'), labels=c('Incidence\nafter intervention','Incidence reduction\nduring intervention','Proportion of early transmissions\njust before intervention','Proportion of early transmissions\nafter intervention'))])
	setkey(dfl, TEAM, OBJ, LATENTO)
	dfl		<- merge(dfl, dfl[, list( LATENTO=LATENTO, LOADcm=cumsum(LOADstd*R2each)), by=c('TEAM','OBJ','X')], by=c('TEAM','OBJ','LATENTO','X'))
	#do barplot	
	ggplot(subset(dfl, LATENTO%in%paste('lv',1:4,sep='')), aes(x=LATENTO, y=100*LOADcm, fill=X, alpha= factor(LATENTO%in%paste('lv',1:4,sep=''),levels=c(TRUE,FALSE),labels=c(1,0)) )) + 
			geom_bar(stat="identity", colour='black') +
			scale_fill_manual(values=c("#762A83","#9970AB","#C2A5CF","#E7D4E8", "#80B1D3", "#FDB462", "#A6DBA0","#5AAE61","#1B7837")) +
			scale_alpha_manual(values=c(1,0.2)) +
			scale_x_discrete(labels=c("1","1-2","1-3","1-4")) +
			scale_y_continuous(expand=c(0,0), limit=c(0,100), breaks=seq(0,100,20), minor_breaks=seq(0,100,10)) +
			facet_grid(TEAM~OBJ) +
			theme_bw() + theme(panel.margin = unit(0.8, "lines"), legend.position='bottom',panel.grid.major.y=element_line(colour='grey70', size=1), panel.grid.minor.y=element_line(colour='grey70', size=0.4)) +
			labs(x='\nfirst n PLS latent factors',fill='Variable',y='variance in error explained\n(%)\n') +
			guides(fill=guide_legend(ncol=3), alpha=FALSE)
	ggsave(file=paste(outdir,'/res_acrossTEAM_Secondary_Errors_PLSbyLatentFactors_v2.pdf',sep=''), width=12, height=12)	
	#do influence plot
	tmp2	<- dfl[, length(unique(X))]
	tmp2	<- dfl[, list(VIP= sqrt(tmp2*sum(R2each*LOADstd)) ), by=c('TEAM','OBJ','X')]
	ggplot(tmp2, aes(x=X, y=VIP, fill=X))	+ geom_bar(stat="identity", colour='black') + 
			scale_y_continuous(expand=c(0,0), limit=c(0,2.9), breaks=seq(0.5,5,0.5), minor_breaks=seq(0,5,0.1)) +
			theme_bw() + coord_flip() + facet_grid(TEAM~OBJ) + 
			theme(legend.position='bottom', axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major.x=element_line(colour='grey70', size=1), panel.grid.minor.x=element_line(colour='grey70', size=0.4)) +
			labs(x='', y='\nVariable influence projection', fill='Variable') +
			guides(fill=guide_legend(ncol=3))
	ggsave(file=paste(outdir,'/res_acrossTEAM_Secondary_Outliers_PLSvip_v2.pdf',sep=''), width=12, height=12)
	#
	# do the same without true values
	#
	plsms2			<- vector('list', dfo[, length(unique(TEAM))])
	plsms3			<- vector('list', dfo[, length(unique(TEAM))])
	plsms5			<- vector('list', dfo[, length(unique(TEAM))])
	plsms6			<- vector('list', dfo[, length(unique(TEAM))])
	names(plsms2)	<- names(plsms3)	<- names(plsms5)	<- names(plsms6)	<- dfo[, unique(TEAM)]
	for(x in names(plsms2))
	{
		df			<- subset(dfo, TEAM==x & OBJ=='Incidence\nafter intervention')	
		plsms2[[x]]	<- plsr(RESID~DATAT_L+IMPRT+SMPL_C+SMPL_M+SMPL_D, data=df, validation='LOO')
		df			<- subset(dfo, TEAM==x & OBJ=='Incidence reduction\nduring intervention')	
		plsms3[[x]]	<- plsr(RESID~DATAT_L+IMPRT+SMPL_C+SMPL_M+SMPL_D, data=df, validation='LOO')
		df			<- subset(dfo, TEAM==x & OBJ=='Proportion of early transmissions\njust before intervention')	
		plsms5[[x]]	<- plsr(RESID~DATAT_L+IMPRT+SMPL_C+SMPL_M+SMPL_D, data=df, validation='LOO')
		df			<- subset(dfo, TEAM==x & OBJ=='Proportion of early transmissions\nafter intervention')	
		plsms6[[x]]	<- plsr(RESID~DATAT_L+IMPRT+SMPL_C+SMPL_M+SMPL_D, data=df, validation='LOO')
	}		
	# get variance explained
	dfr		<- as.data.table(melt(sapply(plsms2, function(x) drop((R2(x, intercept=FALSE, estimate='train'))$val) ), varnames=c('LATENT','TEAM'), value.name='R2'))
	dfr[, OBJ:='Incidence\nafter intervention']
	tmp		<- as.data.table(melt(sapply(plsms3, function(x) drop((R2(x, intercept=FALSE, estimate='train'))$val) ), varnames=c('LATENT','TEAM'), value.name='R2'))
	tmp[, OBJ:='Incidence reduction\nduring intervention']
	dfr		<- rbind(dfr, tmp)
	tmp		<- as.data.table(melt(sapply(plsms5, function(x) drop((R2(x, intercept=FALSE, estimate='train'))$val) ), varnames=c('LATENT','TEAM'), value.name='R2'))
	tmp[, OBJ:='Proportion of early transmissions\njust before intervention']
	dfr		<- rbind(dfr, tmp)	
	tmp		<- as.data.table(melt(sapply(plsms6, function(x) drop((R2(x, intercept=FALSE, estimate='train'))$val) ), varnames=c('LATENT','TEAM'), value.name='R2'))
	tmp[, OBJ:='Proportion of early transmissions\nafter intervention']
	dfr		<- rbind(dfr, tmp)
	set(dfr, NULL, 'LATENT', dfr[, paste('LV',gsub(' comps','',LATENT),sep='')])
	dfr		<- merge(dfr, dfr[, list(LATENT=LATENT, R2each= c(R2[1],diff(R2))), by=c('TEAM','OBJ')], by=c('TEAM','OBJ','LATENT'))
	#oh wow! the variance explained by the various latent variables is quite different from each other!!	
	dfr		<- dfr[order(TEAM, OBJ, -R2each)]
	dfr		<- merge(dfr, dfr[, list(LATENT=LATENT, LATENTO= paste('lv',seq_along(LATENT),sep='')), by=c('TEAM','OBJ')], by=c('TEAM','OBJ','LATENT'))
	set(dfr, NULL, 'LATENTO', dfr[, factor(LATENTO, levels=paste('lv',1:10,sep=''), labels=paste('lv',1:10,sep=''))])	
	#collect key X variables with high latent variable weights
	dfl		<- do.call('rbind',lapply(paste('Comp',1:5), function(x){
						do.call('rbind', 	list(	data.table(TEAM=tmp[, as.character(unique(TEAM))])[, list( X=names(plsms2[[TEAM]]$loading.weights[, 'Comp 1']), LATENT=x, LOAD=plsms2[[TEAM]]$loading.weights[, x], OBJ='Incidence\nafter intervention'), by='TEAM'],
										data.table(TEAM=tmp[, as.character(unique(TEAM))])[, list( X=names(plsms3[[TEAM]]$loading.weights[, 'Comp 1']), LATENT=x, LOAD=plsms3[[TEAM]]$loading.weights[, x], OBJ='Incidence reduction\nduring intervention'), by='TEAM'],
										data.table(TEAM=tmp[, as.character(unique(TEAM))])[, list( X=names(plsms5[[TEAM]]$loading.weights[, 'Comp 1']), LATENT=x, LOAD=plsms5[[TEAM]]$loading.weights[, x], OBJ='Proportion of early transmissions\njust before intervention'), by='TEAM'],
										data.table(TEAM=tmp[, as.character(unique(TEAM))])[, list( X=names(plsms6[[TEAM]]$loading.weights[, 'Comp 1']), LATENT=x, LOAD=plsms6[[TEAM]]$loading.weights[, x], OBJ='Proportion of early transmissions\nafter intervention'), by='TEAM'] 	))						
					}))
	set(dfl, NULL, 'LATENT', dfl[, gsub('Comp ','LV',LATENT)])		
	dfl		<- subset(dfl, !is.nan(LOAD))
	# use variable influence projection (Wold et al 1993)
	# https://books.google.co.uk/books?id=QhHdGt8TG80C&pg=PA2&lpg=PA2&dq=PLS+contribution+of+each+variable&source=bl&ots=vWaqNtCTYz&sig=RT9STQ3SzXk1tU2ZNYlRycgxIQ8&hl=en&sa=X&ved=0ahUKEwiH1eHmirvJAhXK7hoKHcBXC_YQ6AEILzAC#v=onepage&q=PLS%20contribution%20of%20each%20variable&f=false
	dfl		<- merge(dfl, dfl[, list(X=X, LOADstd=LOAD^2/sum(LOAD^2)), by=c('TEAM','OBJ','LATENT')], by=c('TEAM','OBJ','LATENT','X'))
	dfl		<- merge(dfl, dfr, by=c('TEAM','OBJ','LATENT'))
	set(dfl, NULL, 'X', dfl[, factor(X, 	levels=c("DATAT_L","DATA_T","INC_t","INCR_t","ACS_t","ACE_t","IMPRT","SMPL_N","SMPL_C","SMPL_M","SMPL_D"),
							labels=c('Simulation model','Data provided', 'True % incidence','True incidence ratio','True % early transmissions just before intervention','True % early transmissions after intervention','Viral introductions','Sequences (#)','Sequence coverage','Sequences from after intervention start','Sampling duration after intervention start'))])
	setkey(dfl, TEAM, OBJ, LATENTO)
	dfl		<- merge(dfl, dfl[, list( LATENTO=LATENTO, LOADcm=cumsum(LOADstd*R2each)), by=c('TEAM','OBJ','X')], by=c('TEAM','OBJ','LATENTO','X'))
	#do barplot	
	ggplot(subset(dfl, LATENTO%in%paste('lv',1:4,sep='')), aes(x=LATENTO, y=100*LOADcm, fill=X, alpha= factor(LATENTO%in%paste('lv',1:4,sep=''),levels=c(TRUE,FALSE),labels=c(1,0)) )) + 
			geom_bar(stat="identity", colour='black') +
			#scale_fill_brewer(palette='Spectral') +
			scale_alpha_manual(values=c(1,0.2)) +
			scale_x_discrete(labels=c("1","1-2","1-3","1-4")) +
			scale_y_continuous(expand=c(0,0), limit=c(0,100), breaks=seq(0,100,20), minor_breaks=seq(0,100,10)) +
			facet_grid(TEAM~OBJ) +
			theme_bw() + theme(panel.margin = unit(0.8, "lines"), legend.position='bottom',panel.grid.major.y=element_line(colour='grey70', size=1), panel.grid.minor.y=element_line(colour='grey70', size=0.4)) +
			labs(x='\nfirst n PLS latent factors',fill='Variable',y='variance in error explained\n(%)\n') +
			guides(fill=guide_legend(ncol=3), alpha=FALSE)
	ggsave(file=paste(outdir,'/res_acrossTEAM_Secondary_Errors_PLSbyLatentFactors_v3.pdf',sep=''), width=12, height=12)	
	#do influence plot
	tmp2	<- dfl[, length(unique(X))]
	tmp2	<- dfl[, list(VIP= sqrt(tmp2*sum(R2each*LOADstd)) ), by=c('TEAM','OBJ','X')]
	ggplot(tmp2, aes(x=X, y=VIP, fill=X))	+ geom_bar(stat="identity", colour='black') + 
			scale_y_continuous(expand=c(0,0), limit=c(0,2.9), breaks=seq(0.5,5,0.5), minor_breaks=seq(0,5,0.1)) +
			theme_bw() + coord_flip() + facet_grid(TEAM~OBJ) + 
			theme(legend.position='bottom', axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major.x=element_line(colour='grey70', size=1), panel.grid.minor.x=element_line(colour='grey70', size=0.4)) +
			labs(x='', y='\nVariable influence projection', fill='Variable') +
			guides(fill=guide_legend(ncol=3))
	ggsave(file=paste(outdir,'/res_acrossTEAM_Secondary_Outliers_PLSvip_v3.pdf',sep=''), width=12, height=12)



	#
	# Do stepwise model selection in regression 
	#	
	require(gamlss)
	bw.AIC2			<- vector('list', dfo[, length(unique(TEAM))])
	bw.BIC2			<- vector('list', dfo[, length(unique(TEAM))])
	bw.AIC3			<- vector('list', dfo[, length(unique(TEAM))])
	bw.BIC3			<- vector('list', dfo[, length(unique(TEAM))])
	bw.AIC5			<- vector('list', dfo[, length(unique(TEAM))])
	bw.BIC5			<- vector('list', dfo[, length(unique(TEAM))])
	bw.AIC6			<- vector('list', dfo[, length(unique(TEAM))])
	bw.BIC6			<- vector('list', dfo[, length(unique(TEAM))])	
	names(bw.AIC6)	<- names(bw.AIC5)	<- names(bw.AIC3)	<-names(bw.AIC2)	<- dfo[, unique(TEAM)]
	names(bw.BIC6)	<- names(bw.BIC5)	<- names(bw.BIC3)	<-names(bw.BIC2)	<- dfo[, unique(TEAM)]	
	for(x in names(bw.AIC2))
	{
		cat('\nprocess',x)
		tmp				<- subset(dfo, TEAM==x & OBJ=='Incidence\nafter intervention')			
		mnoa			<- gamlss(RESID~DATA_T+DATAT_L+IMPRT+SMPL_C+SMPL_M+SMPL_D+INC_t, data=tmp, family=NO, trace=FALSE)
		bw.AIC2[[x]]	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE)		
		bw.BIC2[[x]] 	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE, k=log(nrow(tmp)))
		tmp				<- subset(dfo, TEAM==x & OBJ=='Incidence reduction\nduring intervention')			
		mnoa			<- gamlss(RESID~DATA_T+DATAT_L+IMPRT+SMPL_C+SMPL_M+SMPL_D+INCR_t, data=tmp, family=NO, trace=FALSE)
		bw.AIC3[[x]]	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE)		
		bw.BIC3[[x]] 	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE, k=log(nrow(tmp)))
		tmp				<- subset(dfo, TEAM==x & OBJ=='Proportion of early transmissions\njust before intervention')			
		mnoa			<- gamlss(RESID~DATA_T+DATAT_L+IMPRT+SMPL_C+SMPL_M+SMPL_D+ACS_t, data=tmp, family=NO, trace=FALSE)
		bw.AIC5[[x]]	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE)		
		bw.BIC5[[x]] 	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE, k=log(nrow(tmp)))
		tmp				<- subset(dfo, TEAM==x & OBJ=='Proportion of early transmissions\nafter intervention')			
		mnoa			<- gamlss(RESID~DATA_T+DATAT_L+IMPRT+SMPL_C+SMPL_M+SMPL_D+ACE_t, data=tmp, family=NO, trace=FALSE)
		bw.AIC6[[x]]	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE)		
		bw.BIC6[[x]] 	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE, k=log(nrow(tmp)))
	}
	dforc			<- dfo[, 	{
				cat(as.character(TEAM), as.character(OBJ))
				zz		<- zy	<- NULL
				if(OBJ=='Incidence\nafter intervention')
				{
					zz	<- bw.AIC2[[as.character(TEAM)]]
					zy	<- bw.BIC2[[as.character(TEAM)]]
				}
				if(OBJ=='Incidence reduction\nduring intervention')
				{
					zz	<- bw.AIC3[[as.character(TEAM)]]
					zy	<- bw.BIC3[[as.character(TEAM)]]
				}
				if(OBJ=='Proportion of early transmissions\njust before intervention')
				{
					zz	<- bw.AIC5[[as.character(TEAM)]]
					zy	<- bw.BIC5[[as.character(TEAM)]]
				}
				if(OBJ=='Proportion of early transmissions\nafter intervention')
				{
					zz	<- bw.AIC6[[as.character(TEAM)]]
					zy	<- bw.BIC6[[as.character(TEAM)]]
				}		
				z 		<- z2	<- numeric(0)
				if(!is.null(names(coef(zz))))
				{
					z		<- summary(zz)[,'Pr(>|t|)']
					z		<- z[ intersect(names(z), names(coef(zz))) ]					
				}								
				if(!is.null(names(coef(zy))))
				{
					z2		<- summary(zy)[,'Pr(>|t|)']
					z2		<- z2[ intersect(names(z2), names(coef(zy))) ]					
				}		
				z			<- sort(z)
				z2			<- sort(z2)
				length(z)	<- max(length(z),length(z2))
				length(z2)	<- max(length(z),length(z2))				
				list( AICn=as.character(names(z)), AICp=z, BICn=as.character(names(z2)), BICp=z2 )	
			}, by=c('TEAM','OBJ')]	
	dforc			<- subset(dforc, AICp<0.05 | BICp<0.05)
	dforc			<- dcast.data.table(dforc, OBJ+BICn~TEAM, value.var='BICp')
	dforc			<- subset(dforc, !grepl('Intercept',BICn) & nchar(BICn)>0)
	setkey(dforc, OBJ, BICn)
	write.csv(dforc, row.names=FALSE, file=paste(outdir,'/res_acrossTEAM_Secondary_Errors.csv',sep=''))
}
##--------------------------------------------------------------------------------------------------------
##	olli 02.12.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.secondary.Errors<- function(dfa, outdir)
{
	require(plsdepot)
	require(pls)
	require(outliers)
	
	dfo		<- subset(dfa, OBJ%in%c('OBJ_ii','OBJ_iii','OBJ_v','OBJ_vi') & USED_GENES=='all' & TEAM!='True' & !grepl('(', TEAM,fixed=1), c(SC_RND, TEAM, DATAT_L, central, DATA_T,IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D, OBJ))
	dfo		<- dcast.data.table(dfo, SC_RND+TEAM+DATAT_L+DATA_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D~OBJ, value.var='central')			
	tmp		<- subset(dfa, OBJ%in%c('OBJ_ii','OBJ_iii','OBJ_v','OBJ_vi') & TEAM=='True', c(SC_RND, TEAM, DATAT_L, central, DATA_T,IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D, OBJ))
	set(tmp, NULL, 'OBJ', tmp[, paste(OBJ,'_t',sep='')])
	tmp		<- dcast.data.table(tmp, SC_RND~OBJ, value.var='central')
	dfo		<- merge(dfo, tmp, by='SC_RND')	
	set(dfo, NULL, 'DATAT_L', dfo[, as.numeric(factor(DATAT_L))])
	set(dfo, NULL, 'SMPL_M', dfo[, as.numeric(SMPL_M)])
	set(dfo, NULL, 'SMPL_D', dfo[, as.numeric(SMPL_D)])	
	set(dfo, NULL, 'DATA_T', dfo[, as.numeric(DATA_T)])
	set(dfo, NULL, 'SMPL_C', dfo[, as.numeric(gsub('%','',as.character(SMPL_C)))/100])
	set(dfo, NULL, 'IMPRT', dfo[, as.numeric(gsub('%','',as.character(IMPRT)))/100])
	setnames(dfo, c('OBJ_ii_t','OBJ_iii_t','OBJ_v_t','OBJ_vi_t'), c('INC_t','INCR_t','ACS_t','ACE_t'))
	#dfo[, R_ii_1:= OBJ_ii-INC_t]
	dfo[, R_ii:= log(OBJ_ii)-log(INC_t)]
	#dfo[, R_iii_1:= OBJ_iii-INCR_t]
	dfo[, R_iii:= log(OBJ_iii)-log(INCR_t)]	
	dfo[, R_v:= OBJ_v-ACS_t]
	#dfo[, R_v_2:= log(OBJ_v)-log(ACS_t)]
	dfo[, R_vi:= OBJ_vi-ACE_t]
	#dfo[, R_vi_2:= log(OBJ_vi)-log(ACE_t)]
	dfo			<- subset(dfo, select=c(SC_RND, TEAM, DATAT_L, DATA_T,IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D, INC_t, INCR_t, ACS_t, ACE_t, R_ii, R_iii, R_v, R_vi))
	dfo			<- melt(dfo, measure.vars=c('R_ii','R_iii','R_v','R_vi'), variable.name='OBJ', value.name='RESID')
	set(dfo, NULL, 'OBJ', dfo[, factor(OBJ, levels=c('R_ii','R_iii','R_v','R_vi'), labels=c('Incidence\nafter intervention','Incidence reduction\nduring intervention','Proportion of early transmissions\njust before intervention','Proportion of early transmissions\nafter intervention'))])
	dfo			<- subset(dfo, !is.na(RESID))
	#	restrict to not Cambridge
	dfo				<- subset(dfo, TEAM!='Cambridge')
	
	#
	#	do PLS on errors
	#	
	plsms2			<- vector('list', dfo[, length(unique(TEAM))])
	plsms3			<- vector('list', dfo[, length(unique(TEAM))])
	plsms5			<- vector('list', dfo[, length(unique(TEAM))])
	plsms6			<- vector('list', dfo[, length(unique(TEAM))])
	names(plsms2)	<- names(plsms3)	<- names(plsms5)	<- names(plsms6)	<- dfo[, unique(TEAM)]
	for(x in names(plsms2))
	{
		df			<- subset(dfo, TEAM==x & OBJ=='Incidence\nafter intervention')	
		plsms2[[x]]	<- plsr(RESID~DATAT_L+DATA_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D+INC_t+INCR_t+ACE_t, data=df, validation='LOO')
		df			<- subset(dfo, TEAM==x & OBJ=='Incidence reduction\nduring intervention')	
		plsms3[[x]]	<- plsr(RESID~DATAT_L+DATA_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D+INC_t+INCR_t+ACE_t, data=df, validation='LOO')
		df			<- subset(dfo, TEAM==x & OBJ=='Proportion of early transmissions\njust before intervention')	
		plsms5[[x]]	<- plsr(RESID~DATAT_L+DATA_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D+INC_t+INCR_t+ACE_t, data=df, validation='LOO')
		df			<- subset(dfo, TEAM==x & OBJ=='Proportion of early transmissions\nafter intervention')	
		plsms6[[x]]	<- plsr(RESID~DATAT_L+DATA_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D+INC_t+INCR_t+ACE_t, data=df, validation='LOO')
	}		
	# get variance explained
	dfr		<- as.data.table(melt(sapply(plsms2, function(x) drop((R2(x, intercept=FALSE, estimate='train'))$val) ), varnames=c('LATENT','TEAM'), value.name='R2'))
	dfr[, OBJ:='Incidence\nafter intervention']
	tmp		<- as.data.table(melt(sapply(plsms3, function(x) drop((R2(x, intercept=FALSE, estimate='train'))$val) ), varnames=c('LATENT','TEAM'), value.name='R2'))
	tmp[, OBJ:='Incidence reduction\nduring intervention']
	dfr		<- rbind(dfr, tmp)
	tmp		<- as.data.table(melt(sapply(plsms5, function(x) drop((R2(x, intercept=FALSE, estimate='train'))$val) ), varnames=c('LATENT','TEAM'), value.name='R2'))
	tmp[, OBJ:='Proportion of early transmissions\njust before intervention']
	dfr		<- rbind(dfr, tmp)	
	tmp		<- as.data.table(melt(sapply(plsms6, function(x) drop((R2(x, intercept=FALSE, estimate='train'))$val) ), varnames=c('LATENT','TEAM'), value.name='R2'))
	tmp[, OBJ:='Proportion of early transmissions\nafter intervention']
	dfr		<- rbind(dfr, tmp)
	set(dfr, NULL, 'LATENT', dfr[, paste('LV',gsub(' comps','',LATENT),sep='')])
	dfr		<- merge(dfr, dfr[, list(LATENT=LATENT, R2each= c(R2[1],diff(R2))), by=c('TEAM','OBJ')], by=c('TEAM','OBJ','LATENT'))
	#oh wow! the variance explained by the various latent variables is quite different from each other!!	
	dfr		<- dfr[order(TEAM, OBJ, -R2each)]
	dfr		<- merge(dfr, dfr[, list(LATENT=LATENT, LATENTO= paste('lv',seq_along(LATENT),sep='')), by=c('TEAM','OBJ')], by=c('TEAM','OBJ','LATENT'))
	set(dfr, NULL, 'LATENTO', dfr[, factor(LATENTO, levels=paste('lv',1:10,sep=''), labels=paste('lv',1:10,sep=''))])
	
	#collect key X variables with high latent variable weights
	dfl		<- do.call('rbind',lapply(paste('Comp',1:10), function(x){
						do.call('rbind', 	list(	data.table(TEAM=tmp[, as.character(unique(TEAM))])[, list( X=names(plsms2[[TEAM]]$loading.weights[, 'Comp 1']), LATENT=x, LOAD=plsms2[[TEAM]]$loading.weights[, x], OBJ='Incidence\nafter intervention'), by='TEAM'],
										data.table(TEAM=tmp[, as.character(unique(TEAM))])[, list( X=names(plsms3[[TEAM]]$loading.weights[, 'Comp 1']), LATENT=x, LOAD=plsms3[[TEAM]]$loading.weights[, x], OBJ='Incidence reduction\nduring intervention'), by='TEAM'],
										data.table(TEAM=tmp[, as.character(unique(TEAM))])[, list( X=names(plsms5[[TEAM]]$loading.weights[, 'Comp 1']), LATENT=x, LOAD=plsms5[[TEAM]]$loading.weights[, x], OBJ='Proportion of early transmissions\njust before intervention'), by='TEAM'],
										data.table(TEAM=tmp[, as.character(unique(TEAM))])[, list( X=names(plsms6[[TEAM]]$loading.weights[, 'Comp 1']), LATENT=x, LOAD=plsms6[[TEAM]]$loading.weights[, x], OBJ='Proportion of early transmissions\nafter intervention'), by='TEAM'] 	))						
					}))
	set(dfl, NULL, 'LATENT', dfl[, gsub('Comp ','LV',LATENT)])		
	dfl		<- subset(dfl, !is.nan(LOAD))
	# use variable influence projection (Wold et al 1993)
	# https://books.google.co.uk/books?id=QhHdGt8TG80C&pg=PA2&lpg=PA2&dq=PLS+contribution+of+each+variable&source=bl&ots=vWaqNtCTYz&sig=RT9STQ3SzXk1tU2ZNYlRycgxIQ8&hl=en&sa=X&ved=0ahUKEwiH1eHmirvJAhXK7hoKHcBXC_YQ6AEILzAC#v=onepage&q=PLS%20contribution%20of%20each%20variable&f=false
	dfl		<- merge(dfl, dfl[, list(X=X, LOADstd=LOAD^2/sum(LOAD^2)), by=c('TEAM','OBJ','LATENT')], by=c('TEAM','OBJ','LATENT','X'))
	dfl		<- merge(dfl, dfr, by=c('TEAM','OBJ','LATENT'))
	set(dfl, NULL, 'X', dfl[, factor(X, 	levels=c("DATAT_L","DATA_T","ACE_t","INC_t","INCR_t","IMPRT","SMPL_N","SMPL_C","SMPL_M","SMPL_D"),
							labels=c('Simulation model','Data provided','%Acute', '%Incidence','Incidence ratio','Viral introductions','Sequences (#)','Sequence coverage','Sequences from after intervention start','Sampling duration after intervention start'))])
	setkey(dfl, TEAM, OBJ, LATENTO)
	dfl		<- merge(dfl, dfl[, list( LATENTO=LATENTO, LOADcm=cumsum(LOADstd*R2each)), by=c('TEAM','OBJ','X')], by=c('TEAM','OBJ','LATENTO','X'))
	#do barplot	
	ggplot(subset(dfl, LATENTO%in%paste('lv',1:4,sep='')), aes(x=LATENTO, y=100*LOADcm, fill=X, alpha= factor(LATENTO%in%paste('lv',1:4,sep=''),levels=c(TRUE,FALSE),labels=c(1,0)) )) + 
			geom_bar(stat="identity", colour='black') +
			#scale_fill_brewer(palette='Spectral') +
			scale_alpha_manual(values=c(1,0.2)) +
			scale_x_discrete(labels=c("1","1-2","1-3","1-4")) +
			scale_y_continuous(expand=c(0,0), limit=c(0,100), breaks=seq(0,100,20), minor_breaks=seq(0,100,10)) +
			facet_grid(TEAM~OBJ) +
			theme_bw() + theme(panel.margin = unit(0.8, "lines"), legend.position='bottom',panel.grid.major.y=element_line(colour='grey70', size=1), panel.grid.minor.y=element_line(colour='grey70', size=0.4)) +
			labs(x='\nfirst n PLS latent factors',fill='Variable',y='variance in error explained\n(%)\n') +
			guides(fill=guide_legend(ncol=3), alpha=FALSE)
	ggsave(file=paste(outdir,'/res_acrossTEAM_Secondary_Errors_PLSbyLatentFactors.pdf',sep=''), width=12, height=12)
	
	#do influence plot
	tmp2	<- dfl[, length(unique(X))]
	tmp2	<- dfl[, list(VIP= sqrt(tmp2*sum(R2each*LOADstd)) ), by=c('TEAM','OBJ','X')]
	ggplot(tmp2, aes(x=X, y=VIP, fill=X))	+ geom_bar(stat="identity", colour='black') + 
			scale_y_continuous(expand=c(0,0), limit=c(0,2.9), breaks=seq(0.5,5,0.5), minor_breaks=seq(0,5,0.1)) +
			theme_bw() + coord_flip() + facet_grid(TEAM~OBJ) + 
			theme(legend.position='bottom', axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major.x=element_line(colour='grey70', size=1), panel.grid.minor.x=element_line(colour='grey70', size=0.4)) +
			labs(x='', y='\nVariable influence projection', fill='Variable') +
			guides(fill=guide_legend(ncol=3))
	ggsave(file=paste(outdir,'/res_acrossTEAM_Secondary_Outliers_PLSvip.pdf',sep=''), width=12, height=12)
	
	#
	# Do stepwise model selection in regression 
	#	
	require(gamlss)
	bw.AIC2			<- vector('list', dfo[, length(unique(TEAM))])
	bw.BIC2			<- vector('list', dfo[, length(unique(TEAM))])
	bw.AIC3			<- vector('list', dfo[, length(unique(TEAM))])
	bw.BIC3			<- vector('list', dfo[, length(unique(TEAM))])
	bw.AIC5			<- vector('list', dfo[, length(unique(TEAM))])
	bw.BIC5			<- vector('list', dfo[, length(unique(TEAM))])
	bw.AIC6			<- vector('list', dfo[, length(unique(TEAM))])
	bw.BIC6			<- vector('list', dfo[, length(unique(TEAM))])	
	names(bw.AIC6)	<- names(bw.AIC5)	<- names(bw.AIC3)	<-names(bw.AIC2)	<- dfo[, unique(TEAM)]
	names(bw.BIC6)	<- names(bw.BIC5)	<- names(bw.BIC3)	<-names(bw.BIC2)	<- dfo[, unique(TEAM)]	
	for(x in names(bw.AIC2))
	{
		cat('\nprocess',x)
		tmp				<- subset(dfo, TEAM==x & OBJ=='Incidence\nafter intervention')			
		mnoa			<- gamlss(RESID~DATAT_L+DATA_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D+INC_t+INCR_t+ACE_t, data=tmp, family=NO, trace=FALSE)
		bw.AIC2[[x]]	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE)		
		bw.BIC2[[x]] 	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE, k=log(nrow(tmp)))
		tmp				<- subset(dfo, TEAM==x & OBJ=='Incidence reduction\nduring intervention')			
		mnoa			<- gamlss(RESID~DATAT_L+DATA_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D+INC_t+INCR_t+ACE_t, data=tmp, family=NO, trace=FALSE)
		bw.AIC3[[x]]	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE)		
		bw.BIC3[[x]] 	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE, k=log(nrow(tmp)))
		tmp				<- subset(dfo, TEAM==x & OBJ=='Proportion of early transmissions\njust before intervention')			
		mnoa			<- gamlss(RESID~DATAT_L+DATA_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D+INC_t+INCR_t+ACE_t, data=tmp, family=NO, trace=FALSE)
		bw.AIC5[[x]]	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE)		
		bw.BIC5[[x]] 	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE, k=log(nrow(tmp)))
		tmp				<- subset(dfo, TEAM==x & OBJ=='Proportion of early transmissions\nafter intervention')			
		mnoa			<- gamlss(RESID~DATAT_L+DATA_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D+INC_t+INCR_t+ACE_t, data=tmp, family=NO, trace=FALSE)
		bw.AIC6[[x]]	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE)		
		bw.BIC6[[x]] 	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE, k=log(nrow(tmp)))
	}
	dforc			<- dfo[, 	{
				cat(as.character(TEAM), as.character(OBJ))
				zz		<- zy	<- NULL
				if(OBJ=='Incidence\nafter intervention')
				{
					zz	<- bw.AIC2[[as.character(TEAM)]]
					zy	<- bw.BIC2[[as.character(TEAM)]]
				}
				if(OBJ=='Incidence reduction\nduring intervention')
				{
					zz	<- bw.AIC3[[as.character(TEAM)]]
					zy	<- bw.BIC3[[as.character(TEAM)]]
				}
				if(OBJ=='Proportion of early transmissions\njust before intervention')
				{
					zz	<- bw.AIC5[[as.character(TEAM)]]
					zy	<- bw.BIC5[[as.character(TEAM)]]
				}
				if(OBJ=='Proportion of early transmissions\nafter intervention')
				{
					zz	<- bw.AIC6[[as.character(TEAM)]]
					zy	<- bw.BIC6[[as.character(TEAM)]]
				}		
				z 		<- z2	<- numeric(0)
				if(!is.null(names(coef(zz))))
				{
					z		<- summary(zz)[,'Pr(>|t|)']
					z		<- z[ intersect(names(z), names(coef(zz))) ]					
				}								
				if(!is.null(names(coef(zy))))
				{
					z2		<- summary(zy)[,'Pr(>|t|)']
					z2		<- z2[ intersect(names(z2), names(coef(zy))) ]					
				}		
				z			<- sort(z)
				z2			<- sort(z2)
				length(z)	<- max(length(z),length(z2))
				length(z2)	<- max(length(z),length(z2))				
				list( AICn=as.character(names(z)), AICp=z, BICn=as.character(names(z2)), BICp=z2 )	
			}, by=c('TEAM','OBJ')]	
	dforc			<- subset(dforc, AICp<0.05 | BICp<0.05)
	dforc			<- dcast.data.table(dforc, OBJ+BICn~TEAM, value.var='BICp')
	dforc			<- subset(dforc, !grepl('Intercept',BICn) & nchar(BICn)>0)
	setkey(dforc, OBJ, BICn)
	write.csv(dforc, row.names=FALSE, file=paste(outdir,'/res_acrossTEAM_Secondary_Errors.csv',sep=''))
}
##--------------------------------------------------------------------------------------------------------
##	olli 06.12.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.secondary.Outliers.v151206<- function(dfa, outdir)
{
	require(plsdepot)
	require(pls)
	require(outliers)
	
	dfo		<- subset(dfa, OBJ%in%c('OBJ_ii','OBJ_iii','OBJ_v','OBJ_vi') & USED_GENES=='all' & TEAM!='True' & !grepl('(', TEAM,fixed=1), c(SC_RND, TEAM, DATAT_L, central, DATA_T,IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D, OBJ))
	dfo		<- dcast.data.table(dfo, SC_RND+TEAM+DATAT_L+DATA_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D~OBJ, value.var='central')			
	tmp		<- subset(dfa, OBJ%in%c('OBJ_ii','OBJ_iii','OBJ_v','OBJ_vi') & TEAM=='True', c(SC_RND, TEAM, DATAT_L, central, DATA_T,IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D, OBJ))
	set(tmp, NULL, 'OBJ', tmp[, paste(OBJ,'_t',sep='')])
	tmp		<- dcast.data.table(tmp, SC_RND~OBJ, value.var='central')
	dfo		<- merge(dfo, tmp, by='SC_RND')	
	set(dfo, NULL, 'DATAT_L', dfo[, as.numeric(factor(DATAT_L, levels=c("Regional","Village"), labels=c("Regional","Village")))])
	set(dfo, NULL, 'SMPL_M', dfo[, as.numeric(SMPL_M)])
	set(dfo, NULL, 'SMPL_D', dfo[, as.numeric(SMPL_D)])	
	set(dfo, NULL, 'DATA_T', dfo[, as.numeric(DATA_T)])
	if(1)
	{
		set(dfo, dfo[, which(IMPRT!='20%')], 'IMPRT', '<=5%')	
		set(dfo, NULL, 'IMPRT', dfo[, as.numeric(factor(as.character(IMPRT), levels=c('<=5%','20%'), labels=c('<=5%','20%')))])	
		set(dfo, dfo[, which(SMPL_C%in%c('8%','30%'))], 'SMPL_C', '1x')
		set(dfo, dfo[, which(SMPL_C%in%c('16%','60%'))], 'SMPL_C', '2x')
		set(dfo, NULL, 'SMPL_C', dfo[, as.numeric(factor(as.character(SMPL_C), levels=c('1x','2x'), labels=c('1x','2x')))])		
	}
	if(0)
	{
		set(dfo, NULL, 'SMPL_C', dfo[, as.numeric(gsub('%','',as.character(SMPL_C)))/100])
		set(dfo, NULL, 'IMPRT', dfo[, as.numeric(gsub('%','',as.character(IMPRT)))/100])		
	}
	setnames(dfo, c('OBJ_ii_t','OBJ_iii_t','OBJ_v_t','OBJ_vi_t'), c('INC_t','INCR_t','ACS_t','ACE_t'))	
	dfo[, R_ii_1:= OBJ_ii-INC_t]
	dfo[, R_ii_2:= log(OBJ_ii)-log(INC_t)]
	dfo[, R_iii_1:= OBJ_iii-INCR_t]
	dfo[, R_iii_2:= log(OBJ_iii)-log(INCR_t)]	
	dfo[, R_v_1:= OBJ_v-ACS_t]
	dfo[, R_v_2:= log(OBJ_v)-log(ACS_t)]
	dfo[, R_vi_1:= OBJ_vi-ACE_t]
	dfo[, R_vi_2:= log(OBJ_vi)-log(ACE_t)]	
	#
	#	Calculate outliers
	#
	set(dfo, NULL, c('R_ii_1','R_iii_1','R_v_2','R_vi_2'), NULL)
	setnames(dfo, c('R_ii_2','R_iii_2','R_v_1','R_vi_1'), c('R_ii','R_iii','R_v','R_vi'))	
	dfo			<- subset(dfo, select=c(SC_RND, TEAM, DATAT_L, DATA_T,IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D, INC_t, INCR_t, ACS_t, ACE_t, R_ii, R_iii, R_v, R_vi))
	dfo			<- melt(dfo, measure.vars=c('R_ii','R_iii','R_v','R_vi'), variable.name='OBJ', value.name='RESID')
	dfo			<- subset(dfo, !is.na(RESID))
	tmp			<- dfo[, {
				z	<- 1.5*diff(quantile(RESID, p=c(0.25,0.75)))
				z	<- c(quantile(RESID, p=0.25)-z, quantile(RESID, p=0.25)+z)				
				list(TEAM=TEAM, SC_RND=SC_RND, OU_GR= as.data.table(grubbs.flag( RESID ))[, OUTLIER], OU_TK= RESID<z[1] | RESID>z[2] )
			}, by='OBJ']
	dfo			<- merge(dfo, tmp, by=c('TEAM','OBJ','SC_RND'))
	dfo[, OUTLIER:=0]
	set(dfo, dfo[, which(OU_TK)], 'OUTLIER', 1)
	set(dfo, dfo[, which(OU_GR)], 'OUTLIER', 2)
	set(dfo, NULL, 'OUTLIER', dfo[, factor(OUTLIER, levels=c(0,1,2), labels=c('No','Mild','Extreme'))])
	set(dfo, NULL, 'OBJ', dfo[, factor(OBJ, levels=c('R_ii','R_iii','R_v','R_vi'), labels=c('Incidence\nafter intervention','Incidence reduction\nduring intervention','Proportion of early transmissions\njust before intervention','Proportion of early transmissions\nafter intervention'))])
	#
	#	do PLS on OU_TK as in http://link.springer.com/article/10.1007/s00439-003-0921-9#/page-1 or http://bib.oxfordjournals.org/content/8/1/32.full
	#
	set(dfo, NULL, 'OU_TK', dfo[, as.numeric(OU_TK)])
	set(dfo, NULL, 'OU_GR', dfo[, as.numeric(OU_GR)])
	#	restrict to not Cambridge
	dfo				<- subset(dfo, TEAM!='Cambridge')
	#
	plsms2			<- vector('list', dfo[, length(unique(TEAM))])
	plsms3			<- vector('list', dfo[, length(unique(TEAM))])
	plsms5			<- vector('list', dfo[, length(unique(TEAM))])
	plsms6			<- vector('list', dfo[, length(unique(TEAM))])
	names(plsms2)	<- names(plsms3)	<- names(plsms5)	<- names(plsms6)	<- dfo[, unique(TEAM)]
	for(x in names(plsms2))
	{
		df			<- subset(dfo, TEAM==x & OBJ=='Incidence\nafter intervention')	
		plsms2[[x]]	<- plsr(OU_TK~DATAT_L+IMPRT+SMPL_C+SMPL_M+INC_t, data=df, validation='LOO', scale=FALSE)
		df			<- subset(dfo, TEAM==x & OBJ=='Incidence reduction\nduring intervention')	
		plsms3[[x]]	<- plsr(OU_TK~DATAT_L+IMPRT+SMPL_C+SMPL_M+INCR_t, data=df, validation='LOO', scale=FALSE)
		df			<- subset(dfo, TEAM==x & OBJ=='Proportion of early transmissions\njust before intervention')	
		if(x=='ETH Zurich')
			plsms5[[x]]	<- plsr(OU_TK~DATAT_L+IMPRT+SMPL_C+SMPL_M+ACS_t, data=df, validation='LOO', scale=TRUE)
		if(x!='ETH Zurich')
			plsms5[[x]]	<- plsr(OU_TK~DATAT_L+IMPRT+SMPL_C+SMPL_M+ACS_t, data=df, validation='LOO', scale=FALSE)
		df			<- subset(dfo, TEAM==x & OBJ=='Proportion of early transmissions\nafter intervention')	
		plsms6[[x]]	<- plsr(OU_TK~DATAT_L+IMPRT+SMPL_C+SMPL_M+ACE_t, data=df, validation='LOO', scale=FALSE)
	}		
	# get variance explained
	dfr		<- as.data.table(melt(sapply(plsms2, function(x) drop((R2(x, intercept=FALSE, estimate='train'))$val) ), varnames=c('LATENT','TEAM'), value.name='R2'))
	dfr[, OBJ:='Incidence\nafter intervention']
	tmp		<- as.data.table(melt(sapply(plsms3, function(x) drop((R2(x, intercept=FALSE, estimate='train'))$val) ), varnames=c('LATENT','TEAM'), value.name='R2'))
	tmp[, OBJ:='Incidence reduction\nduring intervention']
	dfr		<- rbind(dfr, tmp)
	tmp		<- as.data.table(melt(sapply(plsms5, function(x) drop((R2(x, intercept=FALSE, estimate='train'))$val) ), varnames=c('LATENT','TEAM'), value.name='R2'))
	tmp[, OBJ:='Proportion of early transmissions\njust before intervention']
	dfr		<- rbind(dfr, tmp)	
	tmp		<- as.data.table(melt(sapply(plsms6, function(x) drop((R2(x, intercept=FALSE, estimate='train'))$val) ), varnames=c('LATENT','TEAM'), value.name='R2'))
	tmp[, OBJ:='Proportion of early transmissions\nafter intervention']
	dfr		<- rbind(dfr, tmp)
	set(dfr, NULL, 'LATENT', dfr[, paste('LV',gsub(' comps','',LATENT),sep='')])
	dfr		<- merge(dfr, dfr[, list(LATENT=LATENT, R2each= c(R2[1],diff(R2))), by=c('TEAM','OBJ')], by=c('TEAM','OBJ','LATENT'))
	#oh wow! the variance explained by the various latent variables is quite different from each other!!	
	dfr		<- dfr[order(TEAM, OBJ, -R2each)]
	dfr		<- merge(dfr, dfr[, list(LATENT=LATENT, LATENTO= paste('lv',seq_along(LATENT),sep='')), by=c('TEAM','OBJ')], by=c('TEAM','OBJ','LATENT'))
	set(dfr, NULL, 'LATENTO', dfr[, factor(LATENTO, levels=paste('lv',1:10,sep=''), labels=paste('lv',1:10,sep=''))])
	#collect coefficients of first x latent factors
	tmp		<- do.call('rbind',lapply(paste(1:5,'comps'), function(x){
					do.call('rbind', 	list(	data.table(TEAM=tmp[, as.character(unique(TEAM))])[, list( X=names(plsms2[[TEAM]]$coefficients[,, '1 comps']), LATENT=x, CBETA=plsms2[[TEAM]]$coefficients[,, x], OBJ='Incidence\nafter intervention'), by='TEAM'],
												data.table(TEAM=tmp[, as.character(unique(TEAM))])[, list( X=names(plsms3[[TEAM]]$coefficients[,, '1 comps']), LATENT=x, CBETA=plsms3[[TEAM]]$coefficients[,, x], OBJ='Incidence reduction\nduring intervention'), by='TEAM'],
												data.table(TEAM=tmp[, as.character(unique(TEAM))])[, list( X=names(plsms5[[TEAM]]$coefficients[,, '1 comps']), LATENT=x, CBETA=plsms5[[TEAM]]$coefficients[,, x], OBJ='Proportion of early transmissions\njust before intervention'), by='TEAM'],
												data.table(TEAM=tmp[, as.character(unique(TEAM))])[, list( X=names(plsms6[[TEAM]]$coefficients[,, '1 comps']), LATENT=x, CBETA=plsms6[[TEAM]]$coefficients[,, x], OBJ='Proportion of early transmissions\nafter intervention'), by='TEAM'] 	))						
				}))
	tmp		<- subset(tmp, !is.nan(CBETA))
	set(tmp, NULL, 'LATENT', tmp[, paste('LV',gsub(' comps','',LATENT),sep='')])		
	dfr		<- merge(tmp, dfr, by=c('TEAM','OBJ','LATENT'))	
	#collect key X variables with high latent variable weights
	dfl		<- do.call('rbind',lapply(paste('Comp',1:5), function(x){
						do.call('rbind', 	list(	data.table(TEAM=tmp[, as.character(unique(TEAM))])[, list( X=names(plsms2[[TEAM]]$loading.weights[, 'Comp 1']), LATENT=x, LOAD=plsms2[[TEAM]]$loading.weights[, x], OBJ='Incidence\nafter intervention'), by='TEAM'],
													data.table(TEAM=tmp[, as.character(unique(TEAM))])[, list( X=names(plsms3[[TEAM]]$loading.weights[, 'Comp 1']), LATENT=x, LOAD=plsms3[[TEAM]]$loading.weights[, x], OBJ='Incidence reduction\nduring intervention'), by='TEAM'],
													data.table(TEAM=tmp[, as.character(unique(TEAM))])[, list( X=names(plsms5[[TEAM]]$loading.weights[, 'Comp 1']), LATENT=x, LOAD=plsms5[[TEAM]]$loading.weights[, x], OBJ='Proportion of early transmissions\njust before intervention'), by='TEAM'],
													data.table(TEAM=tmp[, as.character(unique(TEAM))])[, list( X=names(plsms6[[TEAM]]$loading.weights[, 'Comp 1']), LATENT=x, LOAD=plsms6[[TEAM]]$loading.weights[, x], OBJ='Proportion of early transmissions\nafter intervention'), by='TEAM'] 	))						
					}))
	set(dfl, NULL, 'LATENT', dfl[, gsub('Comp ','LV',LATENT)])		
	dfl		<- subset(dfl, !is.nan(LOAD))
	# use variable influence projection (Wold et al 1993)
	# https://books.google.co.uk/books?id=QhHdGt8TG80C&pg=PA2&lpg=PA2&dq=PLS+contribution+of+each+variable&source=bl&ots=vWaqNtCTYz&sig=RT9STQ3SzXk1tU2ZNYlRycgxIQ8&hl=en&sa=X&ved=0ahUKEwiH1eHmirvJAhXK7hoKHcBXC_YQ6AEILzAC#v=onepage&q=PLS%20contribution%20of%20each%20variable&f=false
	dfl		<- merge(dfl, dfl[, list(X=X, LOADstd=LOAD^2/sum(LOAD^2)), by=c('TEAM','OBJ','LATENT')], by=c('TEAM','OBJ','LATENT','X'))
	dfl		<- merge(dfl, dfr, by=c('TEAM','OBJ','LATENT','X'), all=1)
	set(dfl, NULL, 'X', dfl[, factor(X, 	levels=c("INC_t","INCR_t","ACS_t","ACE_t","DATAT_L","DATA_T","IMPRT","SMPL_N","SMPL_C","SMPL_M","SMPL_D"),
							labels=c('True % incidence','True incidence ratio','True % early transmissions just before intervention','True % early transmissions after intervention',
										'Village simulation model vs. Regional model','Data provided', 'Frequency of viral introductions 20%/year vs. <=5%/year','Sequences (#)','High sequence coverage (80% for Village, 16% for Regional) vs. lower coverage (40% for Village, 8% for Regional)','Proportion of sequences from after intervention start >80% vs. 50%','Sampling duration after intervention start'))])
	set(dfl, NULL, 'OBJ', dfl[, factor(OBJ, levels=c('Incidence\nafter intervention','Incidence reduction\nduring intervention','Proportion of early transmissions\njust before intervention','Proportion of early transmissions\nafter intervention'), labels=c('Incidence\nafter intervention','Incidence reduction\nduring intervention','Proportion of early transmissions\njust before intervention','Proportion of early transmissions\nafter intervention'))])
	setkey(dfl, TEAM, OBJ, LATENTO)
	dfl		<- merge(dfl, dfl[, list( LATENTO=LATENTO, LOADcm=cumsum(LOADstd*R2each)), by=c('TEAM','OBJ','X')], by=c('TEAM','OBJ','LATENTO','X'))
	dfl[, CBS:= factor(sign(CBETA), levels=c(-1,0,1),labels=c('+','','-'))]
	setkey(dfl, OBJ, TEAM, LATENTO, X)
	dfl		<- merge(dfl, dfl[, list(X=X, POS= cumsum(LOADcm)-0.5*LOADcm), by=c('OBJ','TEAM','LATENTO')],by=c('OBJ','TEAM','LATENTO','X')) 
	
	
	#do barplot	
	ggplot(subset(dfl, LATENTO%in%paste('lv',1:4,sep='')), aes(x=LATENTO, y=100*LOADcm, fill=X)) + 
			geom_bar(stat="identity", colour='black') +
			scale_fill_manual(values=c("#762A83","#9970AB","#C2A5CF","#E7D4E8", "#80B1D3", "#FDB462", "#A6DBA0","#5AAE61","#1B7837")) +
			#scale_fill_brewer(palette='PRGn') +
			scale_x_discrete(labels=c('1','1-2','1-3','1-4')) +
			scale_y_continuous(expand=c(0,0), limit=c(0,100), breaks=seq(0,100,20), minor_breaks=seq(0,100,10)) +
			facet_grid(TEAM~OBJ) +
			theme_bw() + theme(panel.margin = unit(0.8, "lines"), legend.position='bottom', panel.grid.major.y=element_line(colour='grey70', size=1), panel.grid.minor.y=element_line(colour='grey70', size=0.4)) +
			labs(x='\nfirst n PLS latent factors',fill='Strong error\npredictor',y='variance in outlier presence explained\n(%)\n') +
			guides(fill=guide_legend(ncol=3))
	ggsave(file=paste(outdir,'/res_acrossTEAM_Secondary_Outliers_PLSbyLatentFactors_v2.pdf',sep=''), width=12, height=12)
	#	simpler barplot
	ggplot(subset(dfl, LATENTO=='lv4'), aes(x=TEAM)) +
			geom_bar(aes(y=100*LOADcm, fill=X), stat="identity", colour='black') +
			scale_fill_manual(values=c("#762A83","#9970AB","#C2A5CF","#E7D4E8", "#80B1D3", "#FDB462", "#A6DBA0","#5AAE61","#1B7837")) +
			geom_text(aes(y=100*POS, label=CBS)) +
			#scale_fill_brewer(palette='PRGn') +
			#scale_x_discrete(labels=c('1','1-2','1-3','1-4')) +
			scale_y_continuous(expand=c(0,0), limit=c(0,100), breaks=seq(0,100,20), minor_breaks=seq(0,100,10)) +
			facet_grid(~OBJ) +
			coord_flip() +
			theme_bw() + theme(panel.margin = unit(0.8, "lines"), legend.position='bottom', panel.grid.major.x=element_line(colour='grey70', size=1), panel.grid.minor.x=element_line(colour='grey70', size=0.4)) +
			labs(x='',fill='Strong error\npredictor',y='\nvariance in outlier presence explained\n(%)') +
			guides(fill=guide_legend(ncol=2))	
	ggsave(file=paste(outdir,'/res_acrossTEAM_Secondary_Outliers_PLSbyLatentFactors_v3.pdf',sep=''), width=12, height=5)
	
	#do influence plot
	tmp2	<- dfl[, length(unique(X))]
	tmp2	<- dfl[, list(VIP= sqrt(tmp2*sum(R2each*LOADstd)) ), by=c('TEAM','OBJ','X')]
	ggplot(tmp2, aes(x=X, y=VIP, fill=X))	+ geom_bar(stat="identity", colour='black') + 
			scale_y_continuous(expand=c(0,0), limit=c(0,2.9), breaks=seq(0.5,5,0.5), minor_breaks=seq(0,5,0.1)) +
			theme_bw() + coord_flip() + facet_grid(TEAM~OBJ) + 
			theme(legend.position='bottom', axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major.x=element_line(colour='grey70', size=1), panel.grid.minor.x=element_line(colour='grey70', size=0.4)) +
			labs(x='', y='\nVariable influence projection', fill='Variable') +
			guides(fill=guide_legend(ncol=3))
	ggsave(file=paste(outdir,'/res_acrossTEAM_Secondary_Outliers_PLSbyLatentFactors_v2.pdf',sep=''), width=12, height=5)
	#
	#	repeat without true values
	#
	plsms2			<- vector('list', dfo[, length(unique(TEAM))])
	plsms3			<- vector('list', dfo[, length(unique(TEAM))])
	plsms5			<- vector('list', dfo[, length(unique(TEAM))])
	plsms6			<- vector('list', dfo[, length(unique(TEAM))])
	names(plsms2)	<- names(plsms3)	<- names(plsms5)	<- names(plsms6)	<- dfo[, unique(TEAM)]
	for(x in names(plsms2))
	{
		df			<- subset(dfo, TEAM==x & OBJ=='Incidence\nafter intervention')	
		plsms2[[x]]	<- plsr(OU_TK~DATAT_L+IMPRT+SMPL_C+SMPL_M+SMPL_D, data=df, validation='LOO')
		df			<- subset(dfo, TEAM==x & OBJ=='Incidence reduction\nduring intervention')	
		plsms3[[x]]	<- plsr(OU_TK~DATAT_L+IMPRT+SMPL_C+SMPL_M+SMPL_D, data=df, validation='LOO')
		df			<- subset(dfo, TEAM==x & OBJ=='Proportion of early transmissions\njust before intervention')	
		plsms5[[x]]	<- plsr(OU_TK~DATAT_L+IMPRT+SMPL_C+SMPL_M+SMPL_D, data=df, validation='LOO')
		df			<- subset(dfo, TEAM==x & OBJ=='Proportion of early transmissions\nafter intervention')	
		plsms6[[x]]	<- plsr(OU_TK~DATAT_L+IMPRT+SMPL_C+SMPL_M+SMPL_D, data=df, validation='LOO')
	}		
	# get variance explained
	dfr		<- as.data.table(melt(sapply(plsms2, function(x) drop((R2(x, intercept=FALSE, estimate='train'))$val) ), varnames=c('LATENT','TEAM'), value.name='R2'))
	dfr[, OBJ:='Incidence\nafter intervention']
	tmp		<- as.data.table(melt(sapply(plsms3, function(x) drop((R2(x, intercept=FALSE, estimate='train'))$val) ), varnames=c('LATENT','TEAM'), value.name='R2'))
	tmp[, OBJ:='Incidence reduction\nduring intervention']
	dfr		<- rbind(dfr, tmp)
	tmp		<- as.data.table(melt(sapply(plsms5, function(x) drop((R2(x, intercept=FALSE, estimate='train'))$val) ), varnames=c('LATENT','TEAM'), value.name='R2'))
	tmp[, OBJ:='Proportion of early transmissions\njust before intervention']
	dfr		<- rbind(dfr, tmp)	
	tmp		<- as.data.table(melt(sapply(plsms6, function(x) drop((R2(x, intercept=FALSE, estimate='train'))$val) ), varnames=c('LATENT','TEAM'), value.name='R2'))
	tmp[, OBJ:='Proportion of early transmissions\nafter intervention']
	dfr		<- rbind(dfr, tmp)
	set(dfr, NULL, 'LATENT', dfr[, paste('LV',gsub(' comps','',LATENT),sep='')])
	dfr		<- merge(dfr, dfr[, list(LATENT=LATENT, R2each= c(R2[1],diff(R2))), by=c('TEAM','OBJ')], by=c('TEAM','OBJ','LATENT'))
	#oh wow! the variance explained by the various latent variables is quite different from each other!!	
	dfr		<- dfr[order(TEAM, OBJ, -R2each)]
	dfr		<- merge(dfr, dfr[, list(LATENT=LATENT, LATENTO= paste('lv',seq_along(LATENT),sep='')), by=c('TEAM','OBJ')], by=c('TEAM','OBJ','LATENT'))
	set(dfr, NULL, 'LATENTO', dfr[, factor(LATENTO, levels=paste('lv',1:10,sep=''), labels=paste('lv',1:10,sep=''))])
	
	#collect key X variables with high latent variable weights
	dfl		<- do.call('rbind',lapply(paste('Comp',1:5), function(x){
						do.call('rbind', 	list(	data.table(TEAM=tmp[, as.character(unique(TEAM))])[, list( X=names(plsms2[[TEAM]]$loading.weights[, 'Comp 1']), LATENT=x, LOAD=plsms2[[TEAM]]$loading.weights[, x], OBJ='Incidence\nafter intervention'), by='TEAM'],
										data.table(TEAM=tmp[, as.character(unique(TEAM))])[, list( X=names(plsms3[[TEAM]]$loading.weights[, 'Comp 1']), LATENT=x, LOAD=plsms3[[TEAM]]$loading.weights[, x], OBJ='Incidence reduction\nduring intervention'), by='TEAM'],
										data.table(TEAM=tmp[, as.character(unique(TEAM))])[, list( X=names(plsms5[[TEAM]]$loading.weights[, 'Comp 1']), LATENT=x, LOAD=plsms5[[TEAM]]$loading.weights[, x], OBJ='Proportion of early transmissions\njust before intervention'), by='TEAM'],
										data.table(TEAM=tmp[, as.character(unique(TEAM))])[, list( X=names(plsms6[[TEAM]]$loading.weights[, 'Comp 1']), LATENT=x, LOAD=plsms6[[TEAM]]$loading.weights[, x], OBJ='Proportion of early transmissions\nafter intervention'), by='TEAM'] 	))						
					}))
	set(dfl, NULL, 'LATENT', dfl[, gsub('Comp ','LV',LATENT)])		
	dfl		<- subset(dfl, !is.nan(LOAD))
	# use variable influence projection (Wold et al 1993)
	# https://books.google.co.uk/books?id=QhHdGt8TG80C&pg=PA2&lpg=PA2&dq=PLS+contribution+of+each+variable&source=bl&ots=vWaqNtCTYz&sig=RT9STQ3SzXk1tU2ZNYlRycgxIQ8&hl=en&sa=X&ved=0ahUKEwiH1eHmirvJAhXK7hoKHcBXC_YQ6AEILzAC#v=onepage&q=PLS%20contribution%20of%20each%20variable&f=false
	dfl		<- merge(dfl, dfl[, list(X=X, LOADstd=LOAD^2/sum(LOAD^2)), by=c('TEAM','OBJ','LATENT')], by=c('TEAM','OBJ','LATENT','X'))
	dfl		<- merge(dfl, dfr, by=c('TEAM','OBJ','LATENT'))
	set(dfl, NULL, 'X', dfl[, factor(X, 	levels=c("DATAT_L","DATA_T","INC_t","INCR_t","ACS_t","ACE_t","IMPRT","SMPL_N","SMPL_C","SMPL_M","SMPL_D"),
							labels=c('Simulation model','Data provided', 'True % incidence','True incidence ratio','True % early transmissions just before intervention','True % early transmissions after intervention','Viral introductions','Sequences (#)','Sequence coverage','Sequences from after intervention start','Sampling duration after intervention start'))])
	setkey(dfl, TEAM, OBJ, LATENTO)
	dfl		<- merge(dfl, dfl[, list( LATENTO=LATENTO, LOADcm=cumsum(LOADstd*R2each)), by=c('TEAM','OBJ','X')], by=c('TEAM','OBJ','LATENTO','X'))
	#do barplot	
	ggplot(subset(dfl, LATENTO%in%paste('lv',1:4,sep='')), aes(x=LATENTO, y=100*LOADcm, fill=X)) + 
			geom_bar(stat="identity", colour='black') +
			#scale_fill_brewer(palette='Spectral') +
			scale_x_discrete(labels=c('1','1-2','1-3','1-4')) +
			scale_y_continuous(expand=c(0,0), limit=c(0,100), breaks=seq(0,100,20), minor_breaks=seq(0,100,10)) +
			facet_grid(TEAM~OBJ) +
			theme_bw() + theme(panel.margin = unit(0.8, "lines"), legend.position='bottom', panel.grid.major.y=element_line(colour='grey70', size=1), panel.grid.minor.y=element_line(colour='grey70', size=0.4)) +
			labs(x='\nfirst n PLS latent factors',fill='Predictor',y='variance in outlier status explained\n(%)\n') +
			guides(fill=guide_legend(ncol=3))
	ggsave(file=paste(outdir,'/res_acrossTEAM_Secondary_Outliers_PLSbyLatentFactors_v3.pdf',sep=''), width=12, height=12)
	
	#do influence plot
	tmp2	<- dfl[, length(unique(X))]
	tmp2	<- dfl[, list(VIP= sqrt(tmp2*sum(R2each*LOADstd)) ), by=c('TEAM','OBJ','X')]
	ggplot(tmp2, aes(x=X, y=VIP, fill=X))	+ geom_bar(stat="identity", colour='black') + 
			scale_y_continuous(expand=c(0,0), limit=c(0,2.9), breaks=seq(0.5,5,0.5), minor_breaks=seq(0,5,0.1)) +
			theme_bw() + coord_flip() + facet_grid(TEAM~OBJ) + 
			theme(legend.position='bottom', axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major.x=element_line(colour='grey70', size=1), panel.grid.minor.x=element_line(colour='grey70', size=0.4)) +
			labs(x='', y='\nVariable influence projection', fill='Variable') +
			guides(fill=guide_legend(ncol=3))
	ggsave(file=paste(outdir,'/res_acrossTEAM_Secondary_Outliers_PLSvip_v2.pdf',sep=''), width=12, height=12)


	#
	# Do stepwise model selection in regression 
	#	
	require(gamlss)
	bw.AIC2			<- vector('list', dfo[, length(unique(TEAM))])
	bw.BIC2			<- vector('list', dfo[, length(unique(TEAM))])
	bw.AIC3			<- vector('list', dfo[, length(unique(TEAM))])
	bw.BIC3			<- vector('list', dfo[, length(unique(TEAM))])
	bw.AIC5			<- vector('list', dfo[, length(unique(TEAM))])
	bw.BIC5			<- vector('list', dfo[, length(unique(TEAM))])
	bw.AIC6			<- vector('list', dfo[, length(unique(TEAM))])
	bw.BIC6			<- vector('list', dfo[, length(unique(TEAM))])	
	names(bw.AIC6)	<- names(bw.AIC5)	<- names(bw.AIC3)	<-names(bw.AIC2)	<- dfo[, unique(TEAM)]
	names(bw.BIC6)	<- names(bw.BIC5)	<- names(bw.BIC3)	<-names(bw.BIC2)	<- dfo[, unique(TEAM)]	
	for(x in names(bw.AIC2))
	{
		cat('\nprocess',x)
		tmp				<- subset(dfo, TEAM==x & OBJ=='Incidence\nafter intervention')			
		mnoa			<- gamlss(OU_TK~DATAT_L+DATA_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D+INC_t+INCR_t+ACE_t, data=tmp, family=BI, trace=FALSE)
		bw.AIC2[[x]]	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE)		
		bw.BIC2[[x]] 	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE, k=log(nrow(tmp)))
		tmp				<- subset(dfo, TEAM==x & OBJ=='Incidence reduction\nduring intervention')			
		mnoa			<- gamlss(OU_TK~DATAT_L+DATA_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D+INC_t+INCR_t+ACE_t, data=tmp, family=BI, trace=FALSE)
		bw.AIC3[[x]]	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE)		
		bw.BIC3[[x]] 	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE, k=log(nrow(tmp)))
		tmp				<- subset(dfo, TEAM==x & OBJ=='Proportion of early transmissions\njust before intervention')			
		mnoa			<- gamlss(OU_TK~DATAT_L+DATA_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D+INC_t+INCR_t+ACE_t, data=tmp, family=BI, trace=FALSE)
		bw.AIC5[[x]]	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE)		
		bw.BIC5[[x]] 	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE, k=log(nrow(tmp)))
		tmp				<- subset(dfo, TEAM==x & OBJ=='Proportion of early transmissions\nafter intervention')			
		mnoa			<- gamlss(OU_TK~DATAT_L+DATA_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D+INC_t+INCR_t+ACE_t, data=tmp, family=BI, trace=FALSE)
		bw.AIC6[[x]]	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE)		
		bw.BIC6[[x]] 	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE, k=log(nrow(tmp)))
	}
	dfoc			<- dfo[, 	{
				cat(as.character(TEAM), as.character(OBJ))
				zz		<- zy	<- NULL
				if(OBJ=='Incidence\nafter intervention')
				{
					zz	<- bw.AIC2[[as.character(TEAM)]]
					zy	<- bw.BIC2[[as.character(TEAM)]]
				}
				if(OBJ=='Incidence reduction\nduring intervention')
				{
					zz	<- bw.AIC3[[as.character(TEAM)]]
					zy	<- bw.BIC3[[as.character(TEAM)]]
				}
				if(OBJ=='Proportion of early transmissions\njust before intervention')
				{
					zz	<- bw.AIC5[[as.character(TEAM)]]
					zy	<- bw.BIC5[[as.character(TEAM)]]
				}
				if(OBJ=='Proportion of early transmissions\nafter intervention')
				{
					zz	<- bw.AIC6[[as.character(TEAM)]]
					zy	<- bw.BIC6[[as.character(TEAM)]]
				}		
				z 		<- z2	<- numeric(0)
				if(!is.null(names(coef(zz))))
				{
					z		<- summary(zz)[,'Pr(>|t|)']
					z		<- z[ intersect(names(z), names(coef(zz))) ]					
				}								
				if(!is.null(names(coef(zy))))
				{
					z2		<- summary(zy)[,'Pr(>|t|)']
					z2		<- z2[ intersect(names(z2), names(coef(zy))) ]					
				}		
				z			<- sort(z)
				z2			<- sort(z2)
				length(z)	<- max(length(z),length(z2))
				length(z2)	<- max(length(z),length(z2))				
				list( AICn=as.character(names(z)), AICp=z, BICn=as.character(names(z2)), BICp=z2 )	
			}, by=c('TEAM','OBJ')]	
	subset(dfoc, AICp<0.05 | BICp<0.05)
	
	
	subset(dfo, OU_TK & TEAM=='Cambridge/Imperial')
	subset(dfo, OU_TK & TEAM=='Imperial')
	subset(dfo, OU_TK & TEAM=='British Columbia')
	dfo[, table(OU_GR,OU_TK)]
	dfo[, table(OU_TK,OBJ,as.character(TEAM))]
	dfo[, table(OU_TK,OBJ)]
	subset(dfo, OU_TK & OBJ%in%c('R_ii','R_vi'))[, table(SC_RND,as.character(TEAM))]
	subset(dfo, OU_TK & OBJ%in%c('R_ii','R_vi'))[, table(SC_RND,OBJ)]
	subset(dfo, OU_TK)
	subset(dfo, OU_GR)
	
}
##--------------------------------------------------------------------------------------------------------
##	olli 02.12.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.secondary.Outliers<- function(dfa, outdir)
{
	require(plsdepot)
	require(pls)
	require(outliers)
	
	dfo		<- subset(dfa, OBJ%in%c('OBJ_ii','OBJ_iii','OBJ_v','OBJ_vi') & USED_GENES=='all' & TEAM!='True' & !grepl('(', TEAM,fixed=1), c(SC_RND, TEAM, DATAT_L, central, DATA_T,IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D, OBJ))
	dfo		<- dcast.data.table(dfo, SC_RND+TEAM+DATAT_L+DATA_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D~OBJ, value.var='central')			
	tmp		<- subset(dfa, OBJ%in%c('OBJ_ii','OBJ_iii','OBJ_v','OBJ_vi') & TEAM=='True', c(SC_RND, TEAM, DATAT_L, central, DATA_T,IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D, OBJ))
	set(tmp, NULL, 'OBJ', tmp[, paste(OBJ,'_t',sep='')])
	tmp		<- dcast.data.table(tmp, SC_RND~OBJ, value.var='central')
	dfo		<- merge(dfo, tmp, by='SC_RND')	
	set(dfo, NULL, 'DATAT_L', dfo[, as.numeric(factor(DATAT_L))])
	set(dfo, NULL, 'SMPL_M', dfo[, as.numeric(SMPL_M)])
	set(dfo, NULL, 'SMPL_D', dfo[, as.numeric(SMPL_D)])	
	set(dfo, NULL, 'DATA_T', dfo[, as.numeric(DATA_T)])
	set(dfo, NULL, 'SMPL_C', dfo[, as.numeric(gsub('%','',as.character(SMPL_C)))/100])
	set(dfo, NULL, 'IMPRT', dfo[, as.numeric(gsub('%','',as.character(IMPRT)))/100])
	setnames(dfo, c('OBJ_ii_t','OBJ_iii_t','OBJ_v_t','OBJ_vi_t'), c('INC_t','INCR_t','ACS_t','ACE_t'))
	dfo[, R_ii_1:= OBJ_ii-INC_t]
	dfo[, R_ii_2:= log(OBJ_ii)-log(INC_t)]
	dfo[, R_iii_1:= OBJ_iii-INCR_t]
	dfo[, R_iii_2:= log(OBJ_iii)-log(INCR_t)]	
	dfo[, R_v_1:= OBJ_v-ACS_t]
	dfo[, R_v_2:= log(OBJ_v)-log(ACS_t)]
	dfo[, R_vi_1:= OBJ_vi-ACE_t]
	dfo[, R_vi_2:= log(OBJ_vi)-log(ACE_t)]
	
	tmp2		<- melt(dfo, measure.vars=c('R_ii_1','R_ii_2'))	
	ggplot(tmp2, aes(x=value, colour=TEAM)) + geom_density() + facet_grid(~variable) + coord_cartesian(xlim=c(-10,10))
	#	use R_ii_2
	tmp2		<- melt(dfo, measure.vars=c('R_iii_1','R_iii_2'))	
	ggplot(tmp2, aes(x=value)) + geom_histogram(aes(colour=TEAM), fill='transparent') + geom_density(colour='black') + facet_grid(TEAM~variable, scales='free_x')
	#	use R_iii_2	
	tmp2		<- melt(dfo, measure.vars=c('R_v_1','R_v_2'))	
	ggplot(tmp2, aes(x=value)) + geom_histogram(aes(colour=TEAM), fill='transparent') + geom_density(colour='black') + facet_grid(TEAM~variable, scales='free_x')
	#	use R_v_1
	tmp2		<- melt(dfo, measure.vars=c('R_vi_1','R_vi_2'))	
	ggplot(tmp2, aes(x=value)) + geom_histogram(aes(colour=TEAM), fill='transparent') + geom_density(colour='black') + facet_grid(TEAM~variable, scales='free_x')
	#	use R_vi_1
	
	#
	#	Calculate outliers
	#
	set(dfo, NULL, c('R_ii_1','R_iii_1','R_v_2','R_vi_2'), NULL)
	setnames(dfo, c('R_ii_2','R_iii_2','R_v_1','R_vi_1'), c('R_ii','R_iii','R_v','R_vi'))	
	dfo			<- subset(dfo, select=c(SC_RND, TEAM, DATAT_L, DATA_T,IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D, INC_t, INCR_t, ACS_t, ACE_t, R_ii, R_iii, R_v, R_vi))
	dfo			<- melt(dfo, measure.vars=c('R_ii','R_iii','R_v','R_vi'), variable.name='OBJ', value.name='RESID')
	dfo			<- subset(dfo, !is.na(RESID))
	tmp			<- dfo[, {
				z	<- 1.5*diff(quantile(RESID, p=c(0.25,0.75)))
				z	<- c(quantile(RESID, p=0.25)-z, quantile(RESID, p=0.25)+z)				
				list(TEAM=TEAM, SC_RND=SC_RND, OU_GR= as.data.table(grubbs.flag( RESID ))[, OUTLIER], OU_TK= RESID<z[1] | RESID>z[2] )
			}, by='OBJ']
	dfo			<- merge(dfo, tmp, by=c('TEAM','OBJ','SC_RND'))
	dfo[, OUTLIER:=0]
	set(dfo, dfo[, which(OU_TK)], 'OUTLIER', 1)
	set(dfo, dfo[, which(OU_GR)], 'OUTLIER', 2)
	set(dfo, NULL, 'OUTLIER', dfo[, factor(OUTLIER, levels=c(0,1,2), labels=c('No','Mild','Extreme'))])
	set(dfo, NULL, 'OBJ', dfo[, factor(OBJ, levels=c('R_ii','R_iii','R_v','R_vi'), labels=c('Incidence\nafter intervention','Incidence reduction\nduring intervention','Proportion of early transmissions\njust before intervention','Proportion of early transmissions\nafter intervention'))])
	ggplot(dfo, aes(y=TEAM, x=RESID, label=SC_RND, colour=OUTLIER, size=OUTLIER)) + 
			geom_text(position=position_jitter(width=0, height=0.3)) + 
			scale_colour_brewer(palette='Dark2') +
			scale_size_manual(values=c(1,2,3)) +
			theme_bw() + facet_grid(~OBJ, scales='free') +
			labs(x='\nerror in phylogenetic estimate',y='',size='Outlier',colour='Outlier')
	ggsave(file=paste(outdir,'/res_acrossTEAM_Secondary_Outliers.pdf',sep=''), width=12, height=5)
	
	
	#
	#	do PLS on OU_TK as in http://link.springer.com/article/10.1007/s00439-003-0921-9#/page-1 or http://bib.oxfordjournals.org/content/8/1/32.full
	#
	set(dfo, NULL, 'OU_TK', dfo[, as.numeric(OU_TK)])
	set(dfo, NULL, 'OU_GR', dfo[, as.numeric(OU_GR)])
	#	restrict to not Cambridge
	dfo				<- subset(dfo, TEAM!='Cambridge')
	
	plsms2			<- vector('list', dfo[, length(unique(TEAM))])
	plsms3			<- vector('list', dfo[, length(unique(TEAM))])
	plsms5			<- vector('list', dfo[, length(unique(TEAM))])
	plsms6			<- vector('list', dfo[, length(unique(TEAM))])
	names(plsms2)	<- names(plsms3)	<- names(plsms5)	<- names(plsms6)	<- dfo[, unique(TEAM)]
	for(x in names(plsms2))
	{
		df			<- subset(dfo, TEAM==x & OBJ=='Incidence\nafter intervention')	
		plsms2[[x]]	<- plsr(OU_TK~DATAT_L+DATA_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D+INC_t+INCR_t+ACE_t, data=df, validation='LOO')
		df			<- subset(dfo, TEAM==x & OBJ=='Incidence reduction\nduring intervention')	
		plsms3[[x]]	<- plsr(OU_TK~DATAT_L+DATA_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D+INC_t+INCR_t+ACE_t, data=df, validation='LOO')
		df			<- subset(dfo, TEAM==x & OBJ=='Proportion of early transmissions\njust before intervention')	
		plsms5[[x]]	<- plsr(OU_TK~DATAT_L+DATA_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D+INC_t+INCR_t+ACE_t, data=df, validation='LOO')
		df			<- subset(dfo, TEAM==x & OBJ=='Proportion of early transmissions\nafter intervention')	
		plsms6[[x]]	<- plsr(OU_TK~DATAT_L+DATA_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D+INC_t+INCR_t+ACE_t, data=df, validation='LOO')
	}		
	# get variance explained
	dfr		<- as.data.table(melt(sapply(plsms2, function(x) drop((R2(x, intercept=FALSE, estimate='train'))$val) ), varnames=c('LATENT','TEAM'), value.name='R2'))
	dfr[, OBJ:='Incidence\nafter intervention']
	tmp		<- as.data.table(melt(sapply(plsms3, function(x) drop((R2(x, intercept=FALSE, estimate='train'))$val) ), varnames=c('LATENT','TEAM'), value.name='R2'))
	tmp[, OBJ:='Incidence reduction\nduring intervention']
	dfr		<- rbind(dfr, tmp)
	tmp		<- as.data.table(melt(sapply(plsms5, function(x) drop((R2(x, intercept=FALSE, estimate='train'))$val) ), varnames=c('LATENT','TEAM'), value.name='R2'))
	tmp[, OBJ:='Proportion of early transmissions\njust before intervention']
	dfr		<- rbind(dfr, tmp)	
	tmp		<- as.data.table(melt(sapply(plsms6, function(x) drop((R2(x, intercept=FALSE, estimate='train'))$val) ), varnames=c('LATENT','TEAM'), value.name='R2'))
	tmp[, OBJ:='Proportion of early transmissions\nafter intervention']
	dfr		<- rbind(dfr, tmp)
	set(dfr, NULL, 'LATENT', dfr[, paste('LV',gsub(' comps','',LATENT),sep='')])
	dfr		<- merge(dfr, dfr[, list(LATENT=LATENT, R2each= c(R2[1],diff(R2))), by=c('TEAM','OBJ')], by=c('TEAM','OBJ','LATENT'))
	#oh wow! the variance explained by the various latent variables is quite different from each other!!	
	dfr		<- dfr[order(TEAM, OBJ, -R2each)]
	dfr		<- merge(dfr, dfr[, list(LATENT=LATENT, LATENTO= paste('lv',seq_along(LATENT),sep='')), by=c('TEAM','OBJ')], by=c('TEAM','OBJ','LATENT'))
	set(dfr, NULL, 'LATENTO', dfr[, factor(LATENTO, levels=paste('lv',1:10,sep=''), labels=paste('lv',1:10,sep=''))])
	
	#collect key X variables with high latent variable weights
	dfl		<- do.call('rbind',lapply(paste('Comp',1:10), function(x){
						do.call('rbind', 	list(	data.table(TEAM=tmp[, as.character(unique(TEAM))])[, list( X=names(plsms2[[TEAM]]$loading.weights[, 'Comp 1']), LATENT=x, LOAD=plsms2[[TEAM]]$loading.weights[, x], OBJ='Incidence\nafter intervention'), by='TEAM'],
										data.table(TEAM=tmp[, as.character(unique(TEAM))])[, list( X=names(plsms3[[TEAM]]$loading.weights[, 'Comp 1']), LATENT=x, LOAD=plsms3[[TEAM]]$loading.weights[, x], OBJ='Incidence reduction\nduring intervention'), by='TEAM'],
										data.table(TEAM=tmp[, as.character(unique(TEAM))])[, list( X=names(plsms5[[TEAM]]$loading.weights[, 'Comp 1']), LATENT=x, LOAD=plsms5[[TEAM]]$loading.weights[, x], OBJ='Proportion of early transmissions\njust before intervention'), by='TEAM'],
										data.table(TEAM=tmp[, as.character(unique(TEAM))])[, list( X=names(plsms6[[TEAM]]$loading.weights[, 'Comp 1']), LATENT=x, LOAD=plsms6[[TEAM]]$loading.weights[, x], OBJ='Proportion of early transmissions\nafter intervention'), by='TEAM'] 	))						
					}))
	set(dfl, NULL, 'LATENT', dfl[, gsub('Comp ','LV',LATENT)])		
	dfl		<- subset(dfl, !is.nan(LOAD))
	# use variable influence projection (Wold et al 1993)
	# https://books.google.co.uk/books?id=QhHdGt8TG80C&pg=PA2&lpg=PA2&dq=PLS+contribution+of+each+variable&source=bl&ots=vWaqNtCTYz&sig=RT9STQ3SzXk1tU2ZNYlRycgxIQ8&hl=en&sa=X&ved=0ahUKEwiH1eHmirvJAhXK7hoKHcBXC_YQ6AEILzAC#v=onepage&q=PLS%20contribution%20of%20each%20variable&f=false
	dfl		<- merge(dfl, dfl[, list(X=X, LOADstd=LOAD^2/sum(LOAD^2)), by=c('TEAM','OBJ','LATENT')], by=c('TEAM','OBJ','LATENT','X'))
	dfl		<- merge(dfl, dfr, by=c('TEAM','OBJ','LATENT'))
	set(dfl, NULL, 'X', dfl[, factor(X, 	levels=c("DATAT_L","DATA_T","ACE_t","INC_t","INCR_t","IMPRT","SMPL_N","SMPL_C","SMPL_M","SMPL_D"),
							labels=c('Simulation model','Data provided','%Acute', '%Incidence','Incidence ratio','Viral introductions','Sequences (#)','Sequence coverage','Sequences from after intervention start','Sampling duration after intervention start'))])
	setkey(dfl, TEAM, OBJ, LATENTO)
	dfl		<- merge(dfl, dfl[, list( LATENTO=LATENTO, LOADcm=cumsum(LOADstd*R2each)), by=c('TEAM','OBJ','X')], by=c('TEAM','OBJ','LATENTO','X'))
	#do barplot	
	ggplot(subset(dfl, LATENTO%in%paste('lv',1:4,sep='')), aes(x=LATENTO, y=100*LOADcm, fill=X)) + 
			geom_bar(stat="identity", colour='black') +
			#scale_fill_brewer(palette='Spectral') +
			scale_x_discrete(labels=c('1','1-2','1-3','1-4')) +
			scale_y_continuous(expand=c(0,0), limit=c(0,100), breaks=seq(0,100,20), minor_breaks=seq(0,100,10)) +
			facet_grid(TEAM~OBJ) +
			theme_bw() + theme(panel.margin = unit(0.8, "lines"), legend.position='bottom', panel.grid.major.y=element_line(colour='grey70', size=1), panel.grid.minor.y=element_line(colour='grey70', size=0.4)) +
			labs(x='\nfirst n PLS latent factors',fill='Variable',y='variance in outlier status explained\n(%)\n') +
			guides(fill=guide_legend(ncol=3))
	ggsave(file=paste(outdir,'/res_acrossTEAM_Secondary_Outliers_PLSbyLatentFactors.pdf',sep=''), width=12, height=12)
	
	#do influence plot
	tmp2	<- dfl[, length(unique(X))]
	tmp2	<- dfl[, list(VIP= sqrt(tmp2*sum(R2each*LOADstd)) ), by=c('TEAM','OBJ','X')]
	ggplot(tmp2, aes(x=X, y=VIP, fill=X))	+ geom_bar(stat="identity", colour='black') + 
			scale_y_continuous(expand=c(0,0), limit=c(0,2.9), breaks=seq(0.5,5,0.5), minor_breaks=seq(0,5,0.1)) +
			theme_bw() + coord_flip() + facet_grid(TEAM~OBJ) + 
			theme(legend.position='bottom', axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major.x=element_line(colour='grey70', size=1), panel.grid.minor.x=element_line(colour='grey70', size=0.4)) +
			labs(x='', y='\nVariable influence projection', fill='Variable') +
			guides(fill=guide_legend(ncol=3))
	ggsave(file=paste(outdir,'/res_acrossTEAM_Secondary_Outliers_PLSvip.pdf',sep=''), width=12, height=12)
	
	
	#
	# Do stepwise model selection in regression 
	#	
	require(gamlss)
	bw.AIC2			<- vector('list', dfo[, length(unique(TEAM))])
	bw.BIC2			<- vector('list', dfo[, length(unique(TEAM))])
	bw.AIC3			<- vector('list', dfo[, length(unique(TEAM))])
	bw.BIC3			<- vector('list', dfo[, length(unique(TEAM))])
	bw.AIC5			<- vector('list', dfo[, length(unique(TEAM))])
	bw.BIC5			<- vector('list', dfo[, length(unique(TEAM))])
	bw.AIC6			<- vector('list', dfo[, length(unique(TEAM))])
	bw.BIC6			<- vector('list', dfo[, length(unique(TEAM))])	
	names(bw.AIC6)	<- names(bw.AIC5)	<- names(bw.AIC3)	<-names(bw.AIC2)	<- dfo[, unique(TEAM)]
	names(bw.BIC6)	<- names(bw.BIC5)	<- names(bw.BIC3)	<-names(bw.BIC2)	<- dfo[, unique(TEAM)]	
	for(x in names(bw.AIC2))
	{
		cat('\nprocess',x)
		tmp				<- subset(dfo, TEAM==x & OBJ=='Incidence\nafter intervention')			
		mnoa			<- gamlss(OU_TK~DATAT_L+DATA_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D+INC_t+INCR_t+ACE_t, data=tmp, family=BI, trace=FALSE)
		bw.AIC2[[x]]	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE)		
		bw.BIC2[[x]] 	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE, k=log(nrow(tmp)))
		tmp				<- subset(dfo, TEAM==x & OBJ=='Incidence reduction\nduring intervention')			
		mnoa			<- gamlss(OU_TK~DATAT_L+DATA_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D+INC_t+INCR_t+ACE_t, data=tmp, family=BI, trace=FALSE)
		bw.AIC3[[x]]	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE)		
		bw.BIC3[[x]] 	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE, k=log(nrow(tmp)))
		tmp				<- subset(dfo, TEAM==x & OBJ=='Proportion of early transmissions\njust before intervention')			
		mnoa			<- gamlss(OU_TK~DATAT_L+DATA_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D+INC_t+INCR_t+ACE_t, data=tmp, family=BI, trace=FALSE)
		bw.AIC5[[x]]	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE)		
		bw.BIC5[[x]] 	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE, k=log(nrow(tmp)))
		tmp				<- subset(dfo, TEAM==x & OBJ=='Proportion of early transmissions\nafter intervention')			
		mnoa			<- gamlss(OU_TK~DATAT_L+DATA_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D+INC_t+INCR_t+ACE_t, data=tmp, family=BI, trace=FALSE)
		bw.AIC6[[x]]	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE)		
		bw.BIC6[[x]] 	<- stepGAIC.VR(mnoa, direction='backward', what='mu', trace=FALSE, k=log(nrow(tmp)))
	}
	dfoc			<- dfo[, 	{
				cat(as.character(TEAM), as.character(OBJ))
				zz		<- zy	<- NULL
				if(OBJ=='Incidence\nafter intervention')
				{
					zz	<- bw.AIC2[[as.character(TEAM)]]
					zy	<- bw.BIC2[[as.character(TEAM)]]
				}
				if(OBJ=='Incidence reduction\nduring intervention')
				{
					zz	<- bw.AIC3[[as.character(TEAM)]]
					zy	<- bw.BIC3[[as.character(TEAM)]]
				}
				if(OBJ=='Proportion of early transmissions\njust before intervention')
				{
					zz	<- bw.AIC5[[as.character(TEAM)]]
					zy	<- bw.BIC5[[as.character(TEAM)]]
				}
				if(OBJ=='Proportion of early transmissions\nafter intervention')
				{
					zz	<- bw.AIC6[[as.character(TEAM)]]
					zy	<- bw.BIC6[[as.character(TEAM)]]
				}		
				z 		<- z2	<- numeric(0)
				if(!is.null(names(coef(zz))))
				{
					z		<- summary(zz)[,'Pr(>|t|)']
					z		<- z[ intersect(names(z), names(coef(zz))) ]					
				}								
				if(!is.null(names(coef(zy))))
				{
					z2		<- summary(zy)[,'Pr(>|t|)']
					z2		<- z2[ intersect(names(z2), names(coef(zy))) ]					
				}		
				z			<- sort(z)
				z2			<- sort(z2)
				length(z)	<- max(length(z),length(z2))
				length(z2)	<- max(length(z),length(z2))				
				list( AICn=as.character(names(z)), AICp=z, BICn=as.character(names(z2)), BICp=z2 )	
			}, by=c('TEAM','OBJ')]	
	subset(dfoc, AICp<0.05 | BICp<0.05)
	
	
	subset(dfo, OU_TK & TEAM=='Cambridge/Imperial')
	subset(dfo, OU_TK & TEAM=='Imperial')
	subset(dfo, OU_TK & TEAM=='British Columbia')
	dfo[, table(OU_GR,OU_TK)]
	dfo[, table(OU_TK,OBJ,as.character(TEAM))]
	dfo[, table(OU_TK,OBJ)]
	subset(dfo, OU_TK & OBJ%in%c('R_ii','R_vi'))[, table(SC_RND,as.character(TEAM))]
	subset(dfo, OU_TK & OBJ%in%c('R_ii','R_vi'))[, table(SC_RND,OBJ)]
	subset(dfo, OU_TK)
	subset(dfo, OU_GR)
	
}
##--------------------------------------------------------------------------------------------------------
##	olli 26.11.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.secondary.PLSIncidence<- function(dfa, outdir)
{
	require(plsdepot)
	require(pls)
	tmp		<- subset(dfa, USED_GENES=='all' & TEAM!='True' & !grepl('(', TEAM,fixed=1) & OBJ=='OBJ_ii')
	tmp2	<- dcast.data.table(subset(dfa, TEAM=='True', c(SC_RND, OBJ, central)), SC_RND~OBJ, value.var='central')
	setnames(tmp2, c('OBJ_ii','OBJ_iii','OBJ_vi'), c('INC_t','INCR_t','ACE_t'))
	tmp		<- merge(tmp, tmp2, by='SC_RND')
	tmp		<- subset(tmp, !is.na(central), c(SC_RND, TEAM, DATAT_L, central, lower95, upper95, DATA_T,IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D, INC_t, INCR_t, ACE_t))	
	set(tmp, NULL, 'DATAT_L', tmp[, as.numeric(factor(DATAT_L))])
	set(tmp, NULL, 'SMPL_M', tmp[, as.numeric(SMPL_M)])
	set(tmp, NULL, 'SMPL_D', tmp[, as.numeric(SMPL_D)])	
	set(tmp, NULL, 'DATA_T', tmp[, as.numeric(DATA_T)])
	set(tmp, NULL, 'SMPL_C', tmp[, as.numeric(gsub('%','',as.character(SMPL_C)))/100])
	set(tmp, NULL, 'IMPRT', tmp[, as.numeric(gsub('%','',as.character(IMPRT)))/100])	
	tmp[, RES:= central-INC_t]
	
	plsms			<- vector('list', tmp[, length(unique(TEAM))])
	names(plsms)	<- tmp[, unique(TEAM)]
	for(x in names(plsms))
	{
		df			<- subset(tmp, TEAM==x)	
		plsms[[x]]	<- plsr(RES~DATAT_L+DATA_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D+INC_t+INCR_t+ACE_t, data=df, validation='LOO')		
	}
	sapply(plsms, explvar)	
	#first component explains 99% of variance in X
	tmp2			<- sapply(plsms, function(x) drop((R2(x, intercept=FALSE, estimate='train'))$val) )
	tmp2$Cambridge	<- c(tmp2$Cambridge, {x<- rep(NA, 8); names(x)<- paste(3:10,'comps'); x})	
	dfr				<- as.data.table(melt(t(do.call('rbind',tmp2)), varnames=c('LATENT','TEAM'), value.name='R2'))
	set(dfr, NULL, 'LATENT', dfr[, paste('LV',gsub(' comps','',LATENT),sep='')])
	dfr		<- merge(dfr, dfr[, list(LATENT=LATENT, R2each= c(R2[1],diff(R2))), by='TEAM'], by=c('TEAM','LATENT'))
	#oh wow! the variance explained by the various latent variables is quite different from each other!!	
	dfr		<- dfr[order(TEAM, -R2each)]
	dfr		<- merge(dfr, dfr[, list(LATENT=LATENT, LATENTO= paste('lv',seq_along(LATENT),sep='')), by=c('TEAM')], by=c('TEAM','LATENT'))
	set(dfr, NULL, 'LATENTO', dfr[, factor(LATENTO, levels=paste('lv',1:10,sep=''), labels=paste('lv',1:10,sep=''))])
	
	#collect key X variables with high latent variable weights
	dfl		<- do.call('rbind',lapply(paste('Comp',1:10), function(x){
						tmp2	<- NA_real_ 
						if(x%in%colnames(plsms[[TEAM]]$loading.weights))
							tmp2	<- 	plsms[[TEAM]]$loading.weights[, x]
						data.table(TEAM=tmp[, as.character(unique(TEAM))])[, list( X=names(plsms[[TEAM]]$loading.weights[, 'Comp 1']), LATENT=x, LOAD=tmp2), by='TEAM']
					}))
	set(dfl, NULL, 'LATENT', dfl[, gsub('Comp ','LV',LATENT)])		
	# use variable influence projection (Wold et al 1993)
	# https://books.google.co.uk/books?id=QhHdGt8TG80C&pg=PA2&lpg=PA2&dq=PLS+contribution+of+each+variable&source=bl&ots=vWaqNtCTYz&sig=RT9STQ3SzXk1tU2ZNYlRycgxIQ8&hl=en&sa=X&ved=0ahUKEwiH1eHmirvJAhXK7hoKHcBXC_YQ6AEILzAC#v=onepage&q=PLS%20contribution%20of%20each%20variable&f=false
	dfl		<- merge(dfl, dfl[, list(X=X, LOADstd=LOAD^2/sum(LOAD^2)), by=c('TEAM','LATENT')], by=c('TEAM','LATENT','X'))
	dfl		<- merge(dfl, dfr, by=c('TEAM','LATENT'))
	set(dfl, NULL, 'X', dfl[, factor(X, 	levels=c("DATAT_L","DATA_T","ACE_t","INC_t","INCR_t","IMPRT","SMPL_N","SMPL_C","SMPL_M","SMPL_D"),
							labels=c('Simulation model','Data provided','%Acute', '%Incidence','Incidence ratio','Viral introductions','Sequences (#)','Sequence coverage','Sequences from after intervention start','Sampling duration after intervention start'))])
	setkey(dfl, TEAM, LATENTO)
	dfl		<- merge(dfl, dfl[, list( LATENTO=LATENTO, LOADcm=cumsum(LOADstd*R2each)), by=c('TEAM','X')], by=c('TEAM','LATENTO','X'))
	#do barplot
	#does not show total variance explained :-(
	ggplot(dfl, aes(x=LATENTO, y=100*LOADcm, fill=X)) + 
			geom_bar(stat="identity", colour='black') +
			#scale_fill_brewer(palette='Spectral') +
			scale_y_continuous(expand=c(0,0), limit=c(0,100), breaks=seq(0,100,20), minor_breaks=seq(0,100,10)) +
			facet_grid(~TEAM) +
			theme_bw() + theme(legend.position='bottom', axis.text.x=element_blank(),panel.grid.major.y=element_line(colour='grey70', size=1), panel.grid.minor.y=element_line(colour='grey70', size=0.4)) +
			labs(x='\nfirst n PLS latent factors',fill='Variable',y='variance in error explained\n(%)') +
			guides(fill=guide_legend(ncol=3))
	ggsave(file=paste(outdir,'/res_acrossTEAM_Secondary_IncEnd_PLSbyLatentFactors.pdf',sep=''), width=12, height=6)
	
	#do influence plot
	tmp2	<- dfl[, length(unique(X))]
	tmp2	<- dfl[, list(VIP= sqrt(tmp2*sum(R2each*LOADstd)) ), by=c('TEAM','X')]
	ggplot(tmp2, aes(x=X, y=VIP, fill=X))	+ geom_bar(stat="identity", colour='black') + 
			scale_y_continuous(expand=c(0,0), limit=c(0,2.9), breaks=seq(0.5,5,0.5), minor_breaks=seq(0,5,0.1)) +
			theme_bw() + coord_flip() + facet_grid(~TEAM) + 
			theme(legend.position='bottom', axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major.x=element_line(colour='grey70', size=1), panel.grid.minor.x=element_line(colour='grey70', size=0.4)) +
			labs(x='', y='\nVariable influence projection', fill='Variable') +
			guides(fill=guide_legend(ncol=3))
	ggsave(file=paste(outdir,'/res_acrossTEAM_Secondary_IncEnd_PLSvip.pdf',sep=''), width=10, height=5)
	
	# TRAIN
	data("gasoline")
	gasTrain	<- gasoline[1:50, ]
	gasTest 	<- gasoline[51:60, ]
	gas1 		<- plsr(octane ~ NIR, ncomp = 10, data = gasTrain, validation = "LOO")
	plot(gas1, ncomp = 2, asp = 1, line = TRUE)
	plot(RMSEP(gas1), legendpos = "topright")
	plot(gas1, plottype = "scores", comps = 1:3)
	plot(gas1, plottype = "coef", ncomp = 1:3, legendpos = "bottomleft", labels = "numbers", xlab = "nm")
	plot(gas1, plottype = "loadings", comps = 1:5, legendpos = "topleft", labels = "numbers", xlab = "nm")
	plot(gas1, plottype = "correlation", comps = 1:4)
	data("yarn")
	y1 			<- plsr(density ~ NIR, ncomp = 10, data = yarn, validation = "LOO")
	plot(y1, plottype = "scores", comps = 1:3)
	
	df			<- subset(tmp, TEAM=='Cambridge/Imperial')	
	pls1		<- plsr(RES~DATAT_L+DATA_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D+INC_t+INCR_t+ACE_t, data=df, validation='LOO')
	plot(RMSEP(pls1), legendpos = "topright")	#4 components?
	explvar(pls1)								#first component explains 99.9% of variance in X
	plot(pls1, plottype = "loadings", comps = 1:3, labels=c("Model","Data\nprovided","Viral\nintro","SampleN", "SampleC","SampleM","SampleD","Inc","INCR","pcAcute"), legendpos='topright')
	#	the first component has -1 loading at SampleN
	plot(pls1, plottype = "coef", comps = 1, labels=c("Model","Data\nprovided","Viral\nintro","SampleN", "SampleC","SampleM","SampleD","Inc","INCR","pcAcute"), legendpos='topright')
	#	the first component has a dip at SampleN
	plot(gas1, plottype = "correlation", comps = 1:4)
	pls1$coefficients		#for every C: Y~beta*var in component C 
	pls1$scores				#obs x component
	pls1$loadings			#var x component
	pls1$loading.weights	#var x component
	
	#https://books.google.co.uk/books?id=Qsn6yjRXOaMC&pg=PA143&lpg=PA143&dq=what+the+loadings+and+loading+weights+in+pls&source=bl&ots=cB2qb_rUZ1&sig=lIWjCXvw95pejNTLY0c5RO58dm0&hl=en&sa=X&ved=0ahUKEwj9ttSj7rrJAhWGaRQKHZ9nB4oQ6AEIKDAB#v=onepage&q=what%20the%20loadings%20and%20loading%20weights%20in%20pls&f=false
	#Multivariate Data Analysis - in Practice: An Introduction to Multivariate ... By Kim H. Esbensen, Dominique Guyot, Frank Westad, Lars P. Houmolle
	#P loadings are very much like PCA
	#W loading weights. w1 characterizes the first PLS component direction in X space
	#If P and W are similar, then the dominant structures in X happen to be those with max correlation to Y
	
	pls2	<- plsreg1(subset(df, select=c(DATAT_L, DATA_T, IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D, INC_t, INCR_t, ACE_t)), subset(df, select=RES), comps=10, crosval=TRUE)
	#	this seems quite different! The first component explains 63% of variance
}
##--------------------------------------------------------------------------------------------------------
##	olli 26.11.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.secondary.PLSAcute<- function(dfa, outdir)
{
	require(plsdepot)
	require(pls)
	tmp		<- subset(dfa, USED_GENES=='all' & TEAM!='True' & !grepl('(', TEAM,fixed=1) & OBJ=='OBJ_vi')
	tmp2	<- dcast.data.table(subset(dfa, TEAM=='True', c(SC_RND, OBJ, central)), SC_RND~OBJ, value.var='central')
	setnames(tmp2, c('OBJ_ii','OBJ_iii','OBJ_vi'), c('INC_t','INCR_t','ACE_t'))
	tmp		<- merge(tmp, tmp2, by='SC_RND')
	tmp		<- subset(tmp, !is.na(central), c(SC_RND, TEAM, DATAT_L, central, lower95, upper95, DATA_T,IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D, INC_t, INCR_t, ACE_t))	
	set(tmp, NULL, 'DATAT_L', tmp[, as.numeric(factor(DATAT_L))])
	set(tmp, NULL, 'SMPL_M', tmp[, as.numeric(SMPL_M)])
	set(tmp, NULL, 'SMPL_D', tmp[, as.numeric(SMPL_D)])	
	set(tmp, NULL, 'DATA_T', tmp[, as.numeric(DATA_T)])
	set(tmp, NULL, 'SMPL_C', tmp[, as.numeric(gsub('%','',as.character(SMPL_C)))/100])
	set(tmp, NULL, 'IMPRT', tmp[, as.numeric(gsub('%','',as.character(IMPRT)))/100])	
	tmp[, RES:= central-ACE_t]
	
	plsms			<- vector('list', tmp[, length(unique(TEAM))])
	names(plsms)	<- tmp[, unique(TEAM)]
	for(x in names(plsms))
	{
		df			<- subset(tmp, TEAM==x)	
		plsms[[x]]	<- plsr(RES~DATAT_L+DATA_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D+INC_t+INCR_t+ACE_t, data=df, validation='LOO')		
	}
	sapply(plsms, explvar)	
	#first component explains 99% of variance in X
	dfr		<- as.data.table(melt(sapply(plsms, function(x) drop((R2(x, intercept=FALSE, estimate='train'))$val) ), varnames=c('LATENT','TEAM'), value.name='R2'))
	set(dfr, NULL, 'LATENT', dfr[, paste('LV',gsub(' comps','',LATENT),sep='')])
	dfr		<- merge(dfr, dfr[, list(LATENT=LATENT, R2each= c(R2[1],diff(R2))), by='TEAM'], by=c('TEAM','LATENT'))
	#oh wow! the variance explained by the various latent variables is quite different from each other!!	
	dfr		<- dfr[order(TEAM, -R2each)]
	dfr		<- merge(dfr, dfr[, list(LATENT=LATENT, LATENTO= paste('lv',seq_along(LATENT),sep='')), by=c('TEAM')], by=c('TEAM','LATENT'))
	set(dfr, NULL, 'LATENTO', dfr[, factor(LATENTO, levels=paste('lv',1:10,sep=''), labels=paste('lv',1:10,sep=''))])
	
	#collect key X variables with high latent variable weights
	dfl		<- do.call('rbind',lapply(paste('Comp',1:10), function(x){
							data.table(TEAM=tmp[, as.character(unique(TEAM))])[, list( X=names(plsms[[TEAM]]$loading.weights[, 'Comp 1']), LATENT=x, LOAD=plsms[[TEAM]]$loading.weights[, x]), by='TEAM']
						}))
	set(dfl, NULL, 'LATENT', dfl[, gsub('Comp ','LV',LATENT)])		
	# use variable influence projection (Wold et al 1993)
	# https://books.google.co.uk/books?id=QhHdGt8TG80C&pg=PA2&lpg=PA2&dq=PLS+contribution+of+each+variable&source=bl&ots=vWaqNtCTYz&sig=RT9STQ3SzXk1tU2ZNYlRycgxIQ8&hl=en&sa=X&ved=0ahUKEwiH1eHmirvJAhXK7hoKHcBXC_YQ6AEILzAC#v=onepage&q=PLS%20contribution%20of%20each%20variable&f=false
	dfl		<- merge(dfl, dfl[, list(X=X, LOADstd=LOAD^2/sum(LOAD^2)), by=c('TEAM','LATENT')], by=c('TEAM','LATENT','X'))
	dfl		<- merge(dfl, dfr, by=c('TEAM','LATENT'))
	set(dfl, NULL, 'X', dfl[, factor(X, 	levels=c("DATAT_L","DATA_T","ACE_t","INC_t","INCR_t","IMPRT","SMPL_N","SMPL_C","SMPL_M","SMPL_D"),
							labels=c('Simulation model','Data provided','%Acute', '%Incidence','Incidence ratio','Viral introductions','Sequences (#)','Sequence coverage','Sequences from after intervention start','Sampling duration after intervention start'))])
	setkey(dfl, TEAM, LATENTO)
	dfl		<- merge(dfl, dfl[, list( LATENTO=LATENTO, LOADcm=cumsum(LOADstd*R2each)), by=c('TEAM','X')], by=c('TEAM','LATENTO','X'))
	#do barplot
	#does not show total variance explained :-(
	ggplot(dfl, aes(x=LATENTO, y=100*LOADcm, fill=X)) + 
			geom_bar(stat="identity", colour='black') +
			#scale_fill_brewer(palette='Spectral') +
			scale_y_continuous(expand=c(0,0), limit=c(0,100), breaks=seq(0,100,20), minor_breaks=seq(0,100,10)) +
			facet_grid(~TEAM) +
			theme_bw() + theme(legend.position='bottom', axis.text.x=element_blank(),panel.grid.major.y=element_line(colour='grey70', size=1), panel.grid.minor.y=element_line(colour='grey70', size=0.4)) +
			labs(x='\nfirst n PLS latent factors',fill='Variable',y='variance in error explained\n(%)') +
			guides(fill=guide_legend(ncol=3))
	ggsave(file=paste(outdir,'/res_acrossTEAM_Secondary_AcuteEnd_PLSbyLatentFactors.pdf',sep=''), width=10, height=6)
	
	#do influence plot
	tmp2	<- dfl[, length(unique(X))]
	tmp2	<- dfl[, list(VIP= sqrt(tmp2*sum(R2each*LOADstd)) ), by=c('TEAM','X')]
	ggplot(tmp2, aes(x=X, y=VIP, fill=X))	+ geom_bar(stat="identity", colour='black') + 
			scale_y_continuous(expand=c(0,0), limit=c(0,2.9), breaks=seq(0.5,5,0.5), minor_breaks=seq(0,5,0.1)) +
			theme_bw() + coord_flip() + facet_grid(~TEAM) + 
			theme(legend.position='bottom', axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major.x=element_line(colour='grey70', size=1), panel.grid.minor.x=element_line(colour='grey70', size=0.4)) +
			labs(x='', y='\nVariable influence projection', fill='Variable') +
			guides(fill=guide_legend(ncol=3))
	ggsave(file=paste(outdir,'/res_acrossTEAM_Secondary_AcuteEnd_PLSvip.pdf',sep=''), width=10, height=5)
	
	# TRAIN
	data("gasoline")
	gasTrain	<- gasoline[1:50, ]
	gasTest 	<- gasoline[51:60, ]
	gas1 		<- plsr(octane ~ NIR, ncomp = 10, data = gasTrain, validation = "LOO")
	plot(gas1, ncomp = 2, asp = 1, line = TRUE)
	plot(RMSEP(gas1), legendpos = "topright")
	plot(gas1, plottype = "scores", comps = 1:3)
	plot(gas1, plottype = "coef", ncomp = 1:3, legendpos = "bottomleft", labels = "numbers", xlab = "nm")
	plot(gas1, plottype = "loadings", comps = 1:5, legendpos = "topleft", labels = "numbers", xlab = "nm")
	plot(gas1, plottype = "correlation", comps = 1:4)
	data("yarn")
	y1 			<- plsr(density ~ NIR, ncomp = 10, data = yarn, validation = "LOO")
	plot(y1, plottype = "scores", comps = 1:3)
	
	df			<- subset(tmp, TEAM=='Cambridge/Imperial')	
	pls1		<- plsr(RES~DATAT_L+DATA_T+IMPRT+SMPL_N+SMPL_C+SMPL_M+SMPL_D+INC_t+INCR_t+ACE_t, data=df, validation='LOO')
	plot(RMSEP(pls1), legendpos = "topright")	#4 components?
	explvar(pls1)								#first component explains 99.9% of variance in X
	plot(pls1, plottype = "loadings", comps = 1:3, labels=c("Model","Data\nprovided","Viral\nintro","SampleN", "SampleC","SampleM","SampleD","Inc","INCR","pcAcute"), legendpos='topright')
	#	the first component has -1 loading at SampleN
	plot(pls1, plottype = "coef", comps = 1, labels=c("Model","Data\nprovided","Viral\nintro","SampleN", "SampleC","SampleM","SampleD","Inc","INCR","pcAcute"), legendpos='topright')
	#	the first component has a dip at SampleN
	plot(gas1, plottype = "correlation", comps = 1:4)
	pls1$coefficients		#for every C: Y~beta*var in component C 
	pls1$scores				#obs x component
	pls1$loadings			#var x component
	pls1$loading.weights	#var x component
	
	#https://books.google.co.uk/books?id=Qsn6yjRXOaMC&pg=PA143&lpg=PA143&dq=what+the+loadings+and+loading+weights+in+pls&source=bl&ots=cB2qb_rUZ1&sig=lIWjCXvw95pejNTLY0c5RO58dm0&hl=en&sa=X&ved=0ahUKEwj9ttSj7rrJAhWGaRQKHZ9nB4oQ6AEIKDAB#v=onepage&q=what%20the%20loadings%20and%20loading%20weights%20in%20pls&f=false
	#Multivariate Data Analysis - in Practice: An Introduction to Multivariate ... By Kim H. Esbensen, Dominique Guyot, Frank Westad, Lars P. Houmolle
	#P loadings are very much like PCA
	#W loading weights. w1 characterizes the first PLS component direction in X space
	#If P and W are similar, then the dominant structures in X happen to be those with max correlation to Y

	pls2	<- plsreg1(subset(df, select=c(DATAT_L, DATA_T, IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D, INC_t, INCR_t, ACE_t)), subset(df, select=RES), comps=10, crosval=TRUE)
	#	this seems quite different! The first component explains 63% of variance
}
##--------------------------------------------------------------------------------------------------------
##	olli 12.08.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.primary.acute<- function(dfa, outdir, onSeq=1)
{
	require(gridExtra)
	if(onSeq)
	{
		df		<- subset(dfa, Pr_Seq==1)
		title	<- '\nPrimary objective\n%Acute from sequence data\n'
	}		
	if(!onSeq)
	{
		df	<- subset(dfa, Pr_Phy==1)
		title	<- '\nSecondary objective\n%Acute when phylogeny known\n'
	}			
	tmp	<- subset(df, !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_v'))	
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')		
	p1	<- ggplot(subset(tmp, TEAM!='True'), aes(y=gsub('\n',':',paste(INT_L,'  ',AC_L,sep='')), x=central)) +
			#geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) + 
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			coord_cartesian(xlim=c(0, 50)) +
			scale_x_continuous(breaks=seq(10,40,10), minor_breaks=seq(0,50,2)) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM~DATAT_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\n% transmissions from individuals in their first 3 months of infection\nat baseline', y='')
	tmp	<- subset(df, !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_vi'))	
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')		
	p2		<- ggplot(subset(tmp, TEAM!='True'), aes(y=gsub('\n',':',paste(INT_L,'  ',AC_L,sep='')), x=central)) +
			#geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) + 
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +			
			coord_cartesian(xlim=c(0, 50)) +
			scale_x_continuous(breaks=seq(10,40,10), minor_breaks=seq(0,50,2)) +			
			scale_colour_manual(values=TEAM_CL) +			
			facet_grid(TEAM~DATAT_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\n% transmissions from individuals in their first 3 months of infection\nat end of intervention', y='')
	pdf(file=paste(outdir,'/res_acrossTEAM_PrimaryAcute_onSeq',onSeq,'.pdf',sep=''), width = 15, height = 8)
	print(grid.arrange(p1, p2, nrow=1, main=title))
	dev.off() 
}
##--------------------------------------------------------------------------------------------------------
##	olli 12.08.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.secondary.seqcoverage<- function(dfa, outdir)
{
	title	<- '\nSecondary objective\nImpact of sequence coverage\n'
	
	tmp		<- subset(dfa, Sc_SeqCoverage_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_ii'))
	set(tmp, NULL, 'INT_L',tmp[,gsub('\n',': ',INT_L)])
	tmp2	<- as.data.table(expand.grid(SMPL_C=tmp[, unique(SMPL_C)], central=1, DATAT_L=tmp[, unique(DATAT_L)], INT_L=tmp[, unique(INT_L)], AC_L=tmp[,unique(AC_L)]))	
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, INT_L, IMPRT)	
	p1		<- ggplot(subset(tmp, TEAM!='True'), aes(y=SMPL_C, x=central)) +
			geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) +
			coord_cartesian(xlim=c(0, 8)) +
			scale_x_continuous(breaks=seq(1,7,1)) +			
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM~DATAT_L+INT_L+AC_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\n% Incidence', y='sequence coverage\n')
	
	tmp	<- subset(dfa, Sc_SeqCoverage_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_iii'))
	set(tmp, NULL, 'INT_L',tmp[,gsub('\n',': ',INT_L)])
	tmp2	<- as.data.table(expand.grid(SMPL_C=tmp[, unique(SMPL_C)], central=1, DATAT_L=tmp[, unique(DATAT_L)], INT_L=tmp[, unique(INT_L)], AC_L=tmp[,unique(AC_L)]))
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, INT_L, IMPRT)	
	p2	<- ggplot(subset(tmp, TEAM!='True'), aes(y=SMPL_C, x=central)) +
			geom_point(data=tmp2, size=1, colour='transparent') +
			geom_vline(xintercept=1, colour='grey50', size=0.8) +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) +
			coord_cartesian(xlim=c(0, 2)) +
			scale_colour_manual(values=TEAM_CL) +
			scale_x_continuous(breaks=seq(0.5,1.5,0.5), minor_breaks=seq(0.25,1.75,0.25)) +
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +			
			facet_grid(TEAM~DATAT_L+INT_L+AC_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\nIncidence reduction', y='sequence coverage\n')
	
	
	tmp	<- subset(dfa, Sc_SeqCoverage_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_v'))
	set(tmp, NULL, 'INT_L',tmp[,gsub('\n',': ',INT_L)])
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, INT_L, IMPRT)	
	p3	<- ggplot(subset(tmp, TEAM!='True'), aes(y=SMPL_C, x=central)) +
			geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) + 
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			coord_cartesian(xlim=c(0, 50)) +
			scale_x_continuous(breaks=seq(10,40,10), minor_breaks=seq(0,50,2)) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM~DATAT_L+INT_L+AC_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\n% trms from individuals in their first 3 months of infection\nat baseline', y='sequence coverage\n')
	
	
	tmp	<- subset(dfa, Sc_SeqCoverage_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_vi'))
	set(tmp, NULL, 'INT_L',tmp[,gsub('\n',': ',INT_L)])
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, INT_L, IMPRT)	
	p4	<- ggplot(subset(tmp, TEAM!='True'), aes(y=SMPL_C, x=central)) +
			geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) + 
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			coord_cartesian(xlim=c(0, 50)) +
			scale_x_continuous(breaks=seq(10,40,10), minor_breaks=seq(0,50,2)) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM~DATAT_L+INT_L+AC_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\n% trms from individuals in their first 3 months of infection\nat end of intervention', y='sequence coverage\n')
	
	pdf(file=paste(outdir,'/res_acrossTEAM_Secondary_SeqCoverage.pdf',sep=''), width = 15, height = 30)
	print(grid.arrange(p1, p2, p3, p4, nrow=4, main=title))
	dev.off()
}
##--------------------------------------------------------------------------------------------------------
##	olli 12.08.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.secondary.imports<- function(dfa, outdir)
{
	title	<- '\nSecondary objective\nImpact of 20% / year transmissions from outside\n'
	
	tmp	<- subset(dfa, Sc_Imports_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_ii'))
	set(tmp, NULL, 'INT_L',tmp[,gsub('\n',': ',INT_L)])
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, INT_L, IMPRT)	
	p1	<- ggplot(subset(tmp, TEAM!='True'), aes(y=paste(IMPRT,' / year',sep=''), x=central)) +
			#geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) +
			coord_cartesian(xlim=c(0, 8)) +
			scale_x_continuous(breaks=seq(1,7,1)) +			
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM~DATAT_L+INT_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\n% Incidence', y='transmissions from outside\n')
	
	tmp	<- subset(dfa, Sc_Imports_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_iii'))
	set(tmp, NULL, 'INT_L',tmp[,gsub('\n',': ',INT_L)])
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, INT_L, IMPRT)	
	p2	<- ggplot(subset(tmp, TEAM!='True'), aes(y=paste(IMPRT,' / year',sep=''), x=central)) +
			geom_vline(xintercept=1, colour='grey50', size=0.8) +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) +
			coord_cartesian(xlim=c(0, 4)) +
			scale_colour_manual(values=TEAM_CL) +
			scale_x_continuous(breaks=seq(0.5,3.5,0.5), minor_breaks=seq(0.25,3.75,0.25)) +
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +			
			facet_grid(TEAM~DATAT_L+INT_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\nIncidence reduction', y='transmissions from outside\n')
	
	
	tmp	<- subset(dfa, Sc_Imports_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_v'))
	set(tmp, NULL, 'INT_L',tmp[,gsub('\n',': ',INT_L)])
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, INT_L, IMPRT)	
	p3	<- ggplot(subset(tmp, TEAM!='True'), aes(y=paste(IMPRT,' / year',sep=''), x=central)) +
			#geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) + 
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			coord_cartesian(xlim=c(0, 50)) +
			scale_x_continuous(breaks=seq(10,40,10), minor_breaks=seq(0,50,2)) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM~DATAT_L+INT_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\n% trms from individuals in their first 3 months of infection\nat baseline', y='transmissions from outside\n')
	
	
	tmp	<- subset(dfa, Sc_Imports_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_vi'))
	set(tmp, NULL, 'INT_L',tmp[,gsub('\n',': ',INT_L)])
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, INT_L, IMPRT)	
	p4	<- ggplot(subset(tmp, TEAM!='True'), aes(y=paste(IMPRT,' / year',sep=''), x=central)) +
			#geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) + 
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			coord_cartesian(xlim=c(0, 50)) +
			scale_x_continuous(breaks=seq(10,40,10), minor_breaks=seq(0,50,2)) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM~DATAT_L+INT_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\n% trms from individuals in their first 3 months of infection\nat end of intervention', y='transmissions from outside\n')
	
	pdf(file=paste(outdir,'/res_acrossTEAM_Secondary_Imports.pdf',sep=''), width = 12, height = 12)
	print(grid.arrange(p1, p2, p3, p4, nrow=2, main=title))
	dev.off()
}
##--------------------------------------------------------------------------------------------------------
##	olli 12.08.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.secondary.focussedsampling<- function(dfa, outdir)
{
	title	<- '\nSecondary objective\nImpact of 50% vs 85% of samples obtained after intervention start\n'
	
	tmp	<- subset(dfa, Sc_SmplFc_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_ii'))
	set(tmp, NULL, 'INT_L',tmp[,gsub('\n',': ',INT_L)])
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, INT_L, SMPL_M)	
	p1	<- ggplot(subset(tmp, TEAM!='True'), aes(y=factor(as.character(SMPL_M),levels=c('much','extreme'),labels=c('50%','85%')), x=central)) +
			#geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) +
			coord_cartesian(xlim=c(0, 8)) +
			scale_x_continuous(breaks=seq(1,7,1)) +			
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM~DATAT_L+INT_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\n% Incidence', y='% samples obtained after\nintervention start\n')
	
	tmp	<- subset(dfa, Sc_SmplFc_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_iii'))
	set(tmp, NULL, 'INT_L',tmp[,gsub('\n',': ',INT_L)])
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, INT_L, SMPL_M)	
	p2	<- ggplot(subset(tmp, TEAM!='True'), aes(y=factor(as.character(SMPL_M),levels=c('much','extreme'),labels=c('50%','85%')), x=central)) +
			geom_vline(xintercept=1, colour='grey50', size=0.8) +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) +
			coord_cartesian(xlim=c(0, 4)) +
			scale_colour_manual(values=TEAM_CL) +
			scale_x_continuous(breaks=seq(0.5,3.5,0.5), minor_breaks=seq(0.25,3.75,0.25)) +
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +			
			facet_grid(TEAM~DATAT_L+INT_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\nIncidence reduction', y='% samples obtained after\nintervention start\n')
	
	
	tmp	<- subset(dfa, Sc_SmplFc_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_v'))
	set(tmp, NULL, 'INT_L',tmp[,gsub('\n',': ',INT_L)])
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, INT_L, SMPL_M)	
	p3	<- ggplot(subset(tmp, TEAM!='True'), aes(y=factor(as.character(SMPL_M),levels=c('much','extreme'),labels=c('50%','85%')), x=central)) +
			#geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) + 
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			coord_cartesian(xlim=c(0, 50)) +
			scale_x_continuous(breaks=seq(10,40,10), minor_breaks=seq(0,50,2)) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM~DATAT_L+INT_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\n% trms from individuals in their first 3 months of infection\nat baseline', y='% samples obtained after\nintervention start\n')
	
	
	tmp	<- subset(dfa, Sc_SmplFc_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_vi'))
	set(tmp, NULL, 'INT_L',tmp[,gsub('\n',': ',INT_L)])
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, INT_L, SMPL_M)	
	p4	<- ggplot(subset(tmp, TEAM!='True'), aes(y=factor(as.character(SMPL_M),levels=c('much','extreme'),labels=c('50%','85%')), x=central)) +
			#geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) + 
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			coord_cartesian(xlim=c(0, 50)) +
			scale_x_continuous(breaks=seq(10,40,10), minor_breaks=seq(0,50,2)) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM~DATAT_L+INT_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\n% trms from individuals in their first 3 months of infection\nat end of intervention', y='% samples obtained after\nintervention start\n')
	
	
	pdf(file=paste(outdir,'/res_acrossTEAM_Secondary_SFocus.pdf',sep=''), width = 12, height = 12)
	print(grid.arrange(p1, p2, p3, p4, nrow=2, main=title))
	dev.off()
}
##--------------------------------------------------------------------------------------------------------
##	olli 12.08.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.secondary.sduration<- function(dfa, outdir)
{
	title	<- '\nSecondary objective\nImpact of sampling duration 3 yrs vs 5 yrs\n'
	
	tmp	<- subset(dfa, Sc_SmplD_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_ii'))
	set(tmp, NULL, 'AC_L',tmp[,gsub('\n',': ',AC_L)])
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, AC_L, IMPRT)	
	p1	<- ggplot(subset(tmp, TEAM!='True'), aes(y=paste(SMPL_D,' yrs',sep=''), x=central)) +
			#geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) +
			coord_cartesian(xlim=c(0, 8)) +
			scale_x_continuous(breaks=seq(1,7,1)) +			
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM~DATAT_L+AC_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\n% Incidence', y='duration of intensified sampling\nafter intervention start\n')
	
	tmp	<- subset(dfa, Sc_SmplD_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_iii'))
	set(tmp, NULL, 'AC_L',tmp[,gsub('\n',': ',AC_L)])
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, AC_L, IMPRT)	
	p2	<- ggplot(subset(tmp, TEAM!='True'), aes(y=paste(SMPL_D,' yrs',sep=''), x=central)) +
			geom_vline(xintercept=1, colour='grey50', size=0.8) +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) +
			coord_cartesian(xlim=c(0, 4)) +
			scale_colour_manual(values=TEAM_CL) +
			scale_x_continuous(breaks=seq(0.5,3.5,0.5), minor_breaks=seq(0.25,3.75,0.25)) +
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +			
			facet_grid(TEAM~DATAT_L+AC_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\nIncidence reduction', y='duration of intensified sampling\nafter intervention start\n')
	
	
	tmp	<- subset(dfa, Sc_SmplD_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_v'))
	set(tmp, NULL, 'AC_L',tmp[,gsub('\n',': ',AC_L)])
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, AC_L, IMPRT)	
	p3	<- ggplot(subset(tmp, TEAM!='True'), aes(y=paste(SMPL_D,' yrs',sep=''), x=central)) +
			#geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) + 
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			coord_cartesian(xlim=c(0, 50)) +
			scale_x_continuous(breaks=seq(10,40,10), minor_breaks=seq(0,50,2)) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM~DATAT_L+AC_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\n% trms from individuals in their first 3 months of infection\nat baseline', y='duration of intensified sampling\nafter intervention start\n')
	
	
	tmp	<- subset(dfa, Sc_SmplD_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_vi'))
	set(tmp, NULL, 'AC_L',tmp[,gsub('\n',': ',AC_L)])
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, AC_L, IMPRT)	
	p4	<- ggplot(subset(tmp, TEAM!='True'), aes(y=paste(SMPL_D,' yrs',sep=''), x=central)) +
			#geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) + 
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			coord_cartesian(xlim=c(0, 50)) +
			scale_x_continuous(breaks=seq(10,40,10), minor_breaks=seq(0,50,2)) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM~DATAT_L+AC_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\n% trms from individuals in their first 3 months of infection\nat end of intervention', y='duration of intensified sampling\nafter intervention start\n')
	
	pdf(file=paste(outdir,'/res_acrossTEAM_Secondary_SDuration.pdf',sep=''), width = 12, height = 12)
	print(grid.arrange(p1, p2, p3, p4, nrow=2, main=title))
	dev.off()
}
##--------------------------------------------------------------------------------------------------------
##	evaluate results
##	olli 08.05.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.read<- function()
{
	#	read truth for regional simus	
	indir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/PANGEAHIVsim_internal/documents/external/2015_05_results'	
	file	<- paste(indir, '/answers_Regional_Feb2015_rFormat.csv', sep='')
	df		<- read.submission.Feb2015(file, verbose=0, reset.OBJiv.conservative=1)
	#	read truth for village simus
	file	<- paste(indir, '/answers_Village_Feb2015-yr43_rFormat.csv', sep='')
	tmp		<- read.submission.Feb2015(file, verbose=0, reset.OBJiv.conservative=1)
	set(tmp, NULL, 'TEAM', 'True')
	df		<- rbind(df, tmp)
	#	read submissions from May 2015
	tmp		<- list.files(indir, pattern='csv$')
	tmp		<- tmp[!grepl('answers',tmp)]
	#	read Eriks multiple submissions from May 2015
	tmp2	<- data.table(FILE=tmp[grepl('cambImp',tmp)])
	tmp2[, RUN:= tmp2[,  sapply( strsplit(FILE,'_'), function(x) rev(x)[1] )]]
	set(tmp2, NULL, 'RUN', tmp2[, substr(RUN, 1, nchar(RUN)-4)])
	set(tmp2, NULL, 'RUN', tmp2[, gsub('results0','',RUN)])
	dfs		<- do.call('rbind',lapply(seq_len(nrow(tmp2)), function(i)
					{
						z	<- read.submission.Feb2015( paste(indir, '/', tmp2[i, FILE], sep=''), verbose=0, reset.OBJiv.conservative=1 )
						set(z, NULL, 'TEAM', z[, paste(TEAM, ' (', tmp2[i, RUN], ')', sep='')])
						z
					}))
	tmp		<- tmp[!grepl('cambImp',tmp)]
	tmp		<- do.call('rbind',lapply(tmp, function(x) read.submission.Feb2015(paste(indir,'/',x,sep=''), verbose=0, reset.OBJiv.conservative=1)))
	dfs		<- rbind(dfs, tmp)
	#	read submissions from August 2015
	indir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/PANGEAHIVsim_internal/documents/external/2015_08_results'
	tmp		<- list.files(indir, pattern='csv$')
	tmp		<- tmp[!grepl('answers',tmp)]
	tmp2	<- tmp[grepl('Vancouver',tmp)]
	stopifnot(length(tmp2)==1)
	tmp2	<- read.submission.May2015(paste(indir,'/',tmp2,sep=''), verbose=0)
	dfs		<- rbind(dfs, tmp2)
	tmp2	<- tmp[!grepl('Vancouver',tmp)]
	tmp2	<- do.call('rbind',lapply(tmp2, function(x) read.submission.Aug2015(paste(indir,'/',x,sep=''), verbose=0, reset.OBJiv.conservative=1)))
	dfs		<- rbind(dfs, tmp2)
	# 	change team name
	set(dfs, dfs[, which(TEAM=='Colijn')],'TEAM','Imperial')
	#	construct Erik's gold submission
	#	for regional tree, use mergedTab
	tmp		<- subset(dfs, grepl('merged', TEAM) & grepl('Regional',SIM_SCENARIO))	
	tmp[, TEAM:='Cambridge/Imperial']	
	#tmp		<- subset(dfs, grepl('mh30', TEAM) & grepl('Regional',SIM_SCENARIO))	
	#tmp[, TEAM:='Cambridge/Imperial']
	#tmp2	<- subset(dfs, grepl('mh15', TEAM) & grepl('Regional',SIM_SCENARIO))	
	#tmp2[, TEAM:='Cambridge/Imperial']
	#tmp		<- merge(tmp, tmp2, by=c('TEAM','SUBMISSION_DATE','SIM_SCENARIO','USED_GENES','OBJ'), all=1)
	#tmp2	<- tmp[, which(is.na(central.x))]
	#set(tmp, tmp2, 'central.x', tmp[tmp2, central.y])
	#set(tmp, tmp2, 'lower95.x', tmp[tmp2, lower95.y])
	#set(tmp, tmp2, 'upper95.x', tmp[tmp2, upper95.y])
	#setnames(tmp, c('central.x', 'lower95.x', 'upper95.x'), c('central', 'lower95', 'upper95'))
	#set(tmp, NULL, c('central.y', 'lower95.y', 'upper95.y'), NULL)
	dfs		<- rbind(dfs, tmp)
	#	for village tree, use mh30 where available and mh15 where mh30 not available
	tmp		<- subset(dfs, grepl('mh30', TEAM) & grepl('Vill',SIM_SCENARIO))	
	tmp[, TEAM:='Cambridge/Imperial']
	tmp2	<- subset(dfs, grepl('mh15', TEAM) & grepl('Vill',SIM_SCENARIO))	
	tmp2[, TEAM:='Cambridge/Imperial']
	tmp		<- merge(tmp, tmp2, by=c('TEAM','SUBMISSION_DATE','SIM_SCENARIO','USED_GENES','OBJ'), all=1)
	tmp2	<- tmp[, which(is.na(central.x))]
	set(tmp, tmp2, 'central.x', tmp[tmp2, central.y])
	set(tmp, tmp2, 'lower95.x', tmp[tmp2, lower95.y])
	set(tmp, tmp2, 'upper95.x', tmp[tmp2, upper95.y])
	setnames(tmp, c('central.x', 'lower95.x', 'upper95.x'), c('central', 'lower95', 'upper95'))
	set(tmp, NULL, c('central.y', 'lower95.y', 'upper95.y'), NULL)
	dfs		<- rbind(dfs, tmp)
	#	for village seq, use LSD
	tmp		<- subset(dfs, grepl('lsd', TEAM) & grepl('Vill',SIM_SCENARIO))	
	tmp[, TEAM:='Cambridge/Imperial']
	dfs		<- rbind(dfs, tmp)
	#	define data types (seq or phylo)
	dfa		<- rbind(dfs, df)
	dfa[, DATA_T:=NA_character_]
	set(dfa, dfa[, which(grepl('Vill_0[1-4]', SIM_SCENARIO))], 'DATA_T', 'seq')
	set(dfa, dfa[, which(!grepl('Vill_0[1-4]', SIM_SCENARIO))], 'DATA_T', 'phy')	
	set(dfa, dfa[, which(grepl('FirstObj', SIM_SCENARIO))], 'DATA_T', 'seq')
	set(dfa, dfa[, which(grepl('SecondObj', SIM_SCENARIO))], 'DATA_T', 'phy')
	stopifnot(!any(is.na(dfa[, DATA_T])))
	#	define randomized scenario IDs
	dfa[, SC_RND:=NA_character_]
	tmp		<- dfa[, which(grepl('Regional',SIM_SCENARIO))]
	set(dfa, tmp, 'SC_RND', dfa[tmp, substring(regmatches(SIM_SCENARIO,regexpr('sc[A-Z]',SIM_SCENARIO)),3)])
	tmp		<- dfa[, which(grepl('Vill',SIM_SCENARIO))]
	set(dfa, tmp, 'SC_RND', dfa[tmp, substring(regmatches(SIM_SCENARIO,regexpr('Vill_[0-9]+',SIM_SCENARIO)),6)])
	stopifnot(!any(is.na(dfa[, SC_RND])))
	
	#	describe regional simulations in terms of fast/low intervention high/low acute	
	set.seed(42)
	dfi			<- data.table(FILE=list.files('/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2014/2014_Gates/methods_comparison_pipeline/FINAL', '.*zip$', full.names=FALSE))	
	dfi[, SC:= sapply(strsplit(FILE, '_'),'[[',3)]
	dfi[, CONFIG:= sapply(strsplit(SC, '-'),'[[',2)]
	set(dfi, NULL, 'SC', dfi[, sapply(strsplit(SC, '-'),'[[',1)])
	dfi[, DATAT:= sapply(strsplit(FILE, '_'),'[[',5)]
	set(dfi, NULL, 'DATAT', dfi[, gsub('.zip','',DATAT,fixed=T)])	
	set(dfi, NULL, 'OBJECTIVE', 'SecondObj')
	set(dfi, dfi[,which(CONFIG=='sq')],'OBJECTIVE', 'FirstObj')
	dfi			<- merge(dfi,dfi[, list(FILE=FILE, DUMMY=sample(length(FILE),length(FILE))), by='OBJECTIVE'],by=c('OBJECTIVE','FILE'))
	tmp			<- dfi[, which(OBJECTIVE=='SecondObj')]
	set(dfi, tmp, 'DUMMY', dfi[tmp, DUMMY] + dfi[OBJECTIVE=='FirstObj', max(DUMMY)])	
	setkey(dfi, DUMMY)
	dfi[, SC_RND:= toupper(letters[seq_len(nrow(dfi))])]
	dfi			<- subset(dfi, select=c(SC, SC_RND, CONFIG))
	set(dfi, NULL, 'SC', dfi[, substring(SC, 3)])
	dfi			<- merge( dfi, data.table(SC= c('A','B','C','D','E','F'), AC_T=c('low','high','low','high','low','high'), INT_T=c('fast','fast','slow','slow','none','none')), by='SC' )
	tmp			<- data.table(	CONFIG=	c('sq','s2x','y3','mFP85','ph','tr20'),
			IMPRT=	c(.05, .05, .05, .05, .05, .2),
			SMPL_N=	c(1600, 3200, 1280, 1600, 1600, 1600),
			SMPL_C= c(0.08, 0.16, 0.08, 0.08, 0.08, 0.08),
			SMPL_M=	c('overs', 'overs', 'overs', 'extrs', 'overs', 'overs'),
			SMPL_D= c(5, 5, 3, 5, 5, 5))
	dfi			<- merge( dfi, tmp, by='CONFIG')					
	set(dfi, NULL, c('CONFIG'), NULL)
	#	add info for village
	tmp			<- data.table(	SC_RND= c('03','02','01','04','05','08','06','07','11','09','12','10','00'),
			AC_T=	c('low','low','high','high','low','low','high','high','low','low','high','high','low'),
			INT_T=	c('fast','slow','fast','slow','fast','slow','fast','slow','fast','slow','fast','slow','none'),
			#SMPL_C=	c(0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.25, 0.25, 0.25, 0.25, 0.25),
			SMPL_C=	c(0.3, 0.3, 0.3, 0.3, 0.6, 0.6, 0.6, 0.6, 0.3, 0.3, 0.3, 0.3, 0.3),
			SMPL_D= 5,
			SMPL_N= c(777, 857, 957, 1040, 1469, 1630, 1831, 1996, 638, 686, 956, 1012, 872),
			IMPRT=	c(0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0, 0, 0, 0, 0.02))	
	dfi			<- rbind(dfi, tmp, fill=TRUE,use.names=TRUE)
	#	merge to dfa
	cat(paste('\nnumber of rows before merge with dfi, n=', nrow(dfa)))
	dfa			<- merge(dfa, dfi, by='SC_RND')
	cat(paste('\nnumber of rows before merge with dfi, n=', nrow(dfa)))
	#	set answers to numerical
	set(dfa, dfa[, which(OBJ%in%c('OBJ_i','OBJ_iv'))], c('lower95','upper95'), NA_character_)
	set(dfa, dfa[, which(central=='decreasing')], c('central'), '-1')
	set(dfa, dfa[, which(central=='stable')], c('central'), '0')
	set(dfa, dfa[, which(central=='increasing')], c('central'), '1')	
	set(dfa, dfa[, which(central=='<15%')], c('central'), '-1')
	set(dfa, dfa[, which(central=='15%-30%')], c('central'), '0')
	set(dfa, dfa[, which(central=='>30%')], c('central'), '1')
	set(dfa, NULL, 'central', dfa[, as.numeric(central)])
	set(dfa, NULL, 'lower95', dfa[, as.numeric(lower95)])
	set(dfa, NULL, 'upper95', dfa[, as.numeric(upper95)])	
	#	add simulation type
	dfa[, DATAT_L:='NA_character_']
	set(dfa, dfa[, which(grepl('Vill',SIM_SCENARIO))], 'DATAT_L','Village')
	set(dfa, dfa[, which(grepl('Regional',SIM_SCENARIO))], 'DATAT_L','Regional')
	#	transform %incidence 0.01 to 1%
	tmp		<- dfa[, which(OBJ=='OBJ_ii')]
	set(dfa, tmp, 'central', dfa[tmp, 100*central])
	set(dfa, tmp, 'lower95', dfa[tmp, 100*lower95])
	set(dfa, tmp, 'upper95', dfa[tmp, 100*upper95])	
	#	transform %Acute 0.01 to 1%
	tmp		<- dfa[, which(OBJ%in%c('OBJ_v','OBJ_vi'))]
	set(dfa, tmp, 'central', dfa[tmp, 100*central])
	set(dfa, tmp, 'lower95', dfa[tmp, 100*lower95])
	set(dfa, tmp, 'upper95', dfa[tmp, 100*upper95])	
	#	fix submission dates
	set(dfa, NULL, 'SUBMISSION_DATE', dfa[, gsub('\\.15','\\.2015',SUBMISSION_DATE)])
	dfa
}
##--------------------------------------------------------------------------------------------------------
##	evaluate results
##	olli 12.08.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.samplesize<- function(dfr, dfa, outdir)
{
	tmp		<- dfr[, lapply(.SD, sum ), .SDcol=c('Pr_Phy','Pr_Seq','Sc_Imports_Phy','Sc_SeqCoverage_Phy','Sc_SmplD_Phy','Sc_SmplFc_Phy'), by='DATAT_L']
	tmp2	<- dfa[, unique(TEAM[TEAM!='True' & !grepl('(',TEAM,fixed=T)])]
	tmp		<- melt(tmp, id.vars='DATAT_L')
	set(tmp, NULL, 'value', tmp[,value]*length(tmp2))
	tmp		<- dcast.data.table(tmp, DATAT_L~variable, value.var='value')	
	tmp[, N:='TOTAL']
	tmp2	<- unique(subset(dfa, select=c(DATAT_L,OBJ)))
	tmp2[, N:='TOTAL']
	tmp		<- merge(tmp, tmp2, by=c('DATAT_L','N'))	
	tmp2	<- subset(dfa, USED_GENES=='all' & TEAM!='True' & !grepl('(',TEAM,fixed=T))[, lapply(.SD, sum ), by=c('DATAT_L','OBJ','TEAM'), .SDcol=c('Pr_Phy','Pr_Seq','Sc_Imports_Phy','Sc_SeqCoverage_Phy','Sc_SmplD_Phy','Sc_SmplFc_Phy')]
	tmp2[, N:='SUBMITTED']	
	tmp		<- rbind(tmp, tmp2, use.names=T, fill=T)	
	tmp2	<- dcast.data.table(melt(tmp, id.vars=c('DATAT_L','OBJ','N','TEAM')), DATAT_L+OBJ~variable+N,fun.aggregate=sum,value.var='value')
	file	<- paste(outdir,'/SampleSizesByAnalysis_PolGagEnv','.csv',sep='')
	cat('\nWrite to',file)
	write.csv(tmp2, file=file, row.names=FALSE)
	
	set(tmp, tmp[, which(is.na(TEAM))], 'TEAM', 'True')
	z	<- melt(subset(tmp, OBJ=='OBJ_i'), id.vars=c('DATAT_L','N','OBJ','TEAM'))
	set(z, NULL, 'variable', z[, factor(variable, 	levels=c('Pr_Seq','Pr_Phy','Sc_Imports_Phy','Sc_SeqCoverage_Phy','Sc_SmplD_Phy','Sc_SmplFc_Phy'), 
							labels=c('Primary Objective\non sequences','Primary Objective\non true tree','Secondary\nimports','Secondary\nsequence coverage','Secondary\nsampling duration','Secondary\nfocussed sampling'))])
	ggplot(z, aes(x=factor(N, levels=c('SUBMITTED','TOTAL'), labels=c('submitted','if submissions\nhad been\ncomplete')), fill=TEAM, y=value)) + geom_bar(stat='identity') + facet_grid(DATAT_L~variable) +
			scale_fill_manual(values=TEAM_CL) +			
			theme_bw() + theme(legend.position='bottom') + labs(x='', y='total')
	file	<- paste(outdir,'/SampleSizesByAnalysis_PolGagEnv.pdf',sep='')
	cat('\nPlot to',file)
	ggsave(file=file, w=12, h=4)
	NULL	
}