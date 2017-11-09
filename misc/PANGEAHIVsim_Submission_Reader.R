read.example.201502<- function()
{
	file	<- '~/Dropbox\ (Infectious\ Disease)/PANGEAHIVsim/201502/PANGEAHIVsim_Submission_Example.csv'
	file	<- '/Users/Oliver/Downloads/reg5_results0.csv'
	file	<- '~/Dropbox (SPH Imperial College)/PANGEAHIVsim_internal/documents/external/2015_05_results/PANGEAHIVsim_Submission_ColijnV2.csv'
	file	<- '~/Dropbox (SPH Imperial College)/PANGEAHIVsim_internal/documents/external/2015_05_results/cambImp_regional_mh15.csv'	
	file	<- '~/Dropbox (SPH Imperial College)/PANGEAHIVsim_internal/documents/external/2015_05_results/PANGEAHIVsim_Submission_Vancouver.csv'
	file	<- '~/Dropbox (SPH Imperial College)/PANGEAHIVsim_internal/documents/external/2015_05_results/PANGEAHIVsim_Submission_BD.csv'
	file	<- '~/Dropbox (SPH Imperial College)/PANGEAHIVsim_internal/documents/external/2015_05_results/PANGEAHIVsim_Submission_ETHZurich_2.csvERROR'
	file	<- '~/Dropbox (SPH Imperial College)/PANGEAHIVsim_internal/documents/external/2015_05_results/PANGEAHIVsim_Submission_ETHZurich.csv'
	ans2	<- read.submission.Feb2015(file, verbose=0)
}

read.example.201505<- function()
{
	file	<- '~/Dropbox\ (Infectious\ Disease)/PANGEAHIVsim/201502/PANGEAHIVsim_Submission_May2015_Example.csv'
	ans2	<- read.submission.May2015(file, verbose=0)
}

read.example.201508<- function()
{
	file	<- '~/Dropbox (SPH Imperial College)/PANGEAHIVsim_internal/documents/external/2015_08_results/previous_got_updated/PANGEAHIVsim_Submission_ETHZurich_201508.csv'
	file	<- '~/Dropbox (SPH Imperial College)/PANGEAHIVsim_internal/documents/external/2015_08_results/PANGEAHIVsim_Submission_RegSeqs_Colijn.csv'
	file	<- '~/Dropbox (SPH Imperial College)/PANGEAHIVsim_internal/documents/external/2015_08_results/PANGEAHIVsim_Submission_ETHZurich_20150818.csv'
	ans2	<- read.submission.Aug2015(file, verbose=0)
	file	<- '~/Dropbox (SPH Imperial College)/PANGEAHIVsim_internal/documents/external/2015_08_results/PANGEAHIVsim_Submission_Vancouver_Aug2015.csv'
	ans2	<- read.submission.May2015(file, verbose=0)
}

read.submission.Aug2015<- function(file, verbose=1, warn.all=0, reset.OBJiv.conservative=1)
{
	require(data.table)
	cat('\nThis is read.submission version 15-08-12.')
	if(verbose)
	{
		cat('\nReminder of objectives\n(use verbose=0 to suppress this message)')
		
		cat('\nOBJ_i\tDuring the evaluation period, was incidence stable, declining or increasing?\n\t\tAnswer: "stable", "declining", "increasing", or "NA"')
		cat('\nOBJ_ii\tWhat is the annual % incidence in the last year of the evaluation period?\n\t\tAnswer: numerical or "NA"')
		cat('\nOBJ_iii\tComparing the last year of the evaluation period to the year preceding the evaluation period, what is the ratio in annual % incidence?\n\t\tAnswer: numerical or "NA"')
		
		cat('\nOBJ_iv\tWas the proportion of transmissions that originated from individuals in early HIV infection in the year preceding the evaluation period below 10%, between 10-30%, or above 30%?\n\t\tAnswer: "<10%", "10-30%", ">30%" or "NA"')
		cat('\nOBJ_v\tWhat is the proportion of transmissions that originated from individuals in early HIV infection in the year preceding the evaluation period?\n\t\tAnswer: numerical or "NA"')
		cat('\nOBJ_vi\tWhat is the proportion of transmissions that originated from individuals in early HIV infection in the last year of the evaluation period?\n\t\tAnswer: numerical or "NA"')
		
		cat('\n\nPlease use:\nscenario names SIM_SCENARIO that correspond to the available file names\nUSED_GENES\teither "pol" or "all"')
	}
	cat('\nSkip first 22 rows that should contain comment rows starting with #.')
	cat(paste('\nreading', file))
	df	<- read.csv(file, stringsAsFactors=FALSE, comment.char="#", blank.lines.skip=TRUE, skip=22)
	df	<- as.data.table(df)	
	#	check column names
	df.colnm	<- c( "TEAM", "SUBMISSION_DATE", "SIM_SCENARIO", "USED_GENES", "OBJ_i", "OBJ_ii", "OBJ_iii", "OBJ_iv", "OBJ_v", "OBJ_vi", "ESTIMATE" )
	tmp			<- setdiff( df.colnm, names(df) )
	if(length(tmp))
		stop(paste('Found missing columns, ', paste(tmp, collapse=','), sep=''))
	tmp			<- setdiff( names(df), df.colnm )
	if(length(tmp))
		warning(paste('Ignore extra columns, ', paste(tmp, collapse=','), sep=''))
	df			<- df[, df.colnm, with=0]
	#
	#	check columns
	#
	#	check TEAM
	if(df[, length(unique(TEAM))>1])
		stop('Found more than one TEAM name.')
	#	check SUBMISSION_DATE
	if(df[, length(unique(SUBMISSION_DATE))>1])
		stop('Found more than one SUBMISSION_DATE.')
	#	check SIM_SCENARIO
	df.sc		<- data.table(SIM_SCENARIO= c( 	paste('Vill_',sprintf("%02d",0:12),'_Feb2015_5yr',sep=''), 
					paste('Vill_',sprintf("%02d",c(0:1,4,6:7,9:12)),'_Feb2015',sep=''),
					paste('Vill_',sprintf("%02d",c(2:3,5,8)),'_Feb2015_3yr',sep=''),
					paste('150129_PANGEAsim_Regional_FirstObj_sc',LETTERS[seq(1,4)],'_SIMULATED_SEQ',sep=''),
					paste('REGIONAL_sc',LETTERS[seq(1,20)],sep=''),
					paste('150129_PANGEAsim_Regional_SecondObj_sc',LETTERS[seq(5,20)],'_SIMULATED_DATEDTREE',sep='')												
			)) 
	tmp			<- unique( setdiff( df[, SIM_SCENARIO], df.sc[, SIM_SCENARIO] ) )
	if(length(tmp))
		stop(paste('Found invalid scenarios',paste(tmp,collapse=',')))
	#	if scenarios are coded as REGIONAL, transform back to original SCENARIO code
	tmp			<- df[, which(grepl('REGIONAL_sc[A-D]', SIM_SCENARIO))]	
	if(length(tmp))
	{
		if(verbose)
			cat('Found SIM_SCENARIO REGIONAL_sc[A-D]. Re-code as 150129_PANGEAsim_Regional_FirstObj_sc[A-D]_SIMULATED_SEQ. n=', length(tmp))
		set(df, tmp, 'SIM_SCENARIO', df[tmp,paste(gsub('REGIONAL_','150129_PANGEAsim_Regional_FirstObj_',SIM_SCENARIO),'_SIMULATED_SEQ',sep='')])
	}
	tmp			<- df[, which(grepl('REGIONAL_sc[E-T]', SIM_SCENARIO))]	
	if(length(tmp))
	{
		if(verbose)
			cat('Found SIM_SCENARIO REGIONAL_sc[E-T]. Re-code as 150129_PANGEAsim_Regional_SecondObj_sc[E-T]_SIMULATED_DATEDTREE. n=', length(tmp))
		set(df, tmp, 'SIM_SCENARIO', df[tmp,paste(gsub('REGIONAL_','150129_PANGEAsim_Regional_SecondObj_',SIM_SCENARIO),'_SIMULATED_DATEDTREE',sep='')])
	}	
	df			<- merge(df, df.sc, by='SIM_SCENARIO')		
	#	check USED_GENES
	set(df, NULL, 'USED_GENES', df[, tolower(gsub('\\s','',as.character(USED_GENES)))])
	tmp			<- df[, which( !USED_GENES%in%c('pol','all') )]
	if(length(tmp))
		stop(paste('Found invalid USED_GENES', paste(df[tmp,USED_GENES], collapse=',')))
	tmp			<- df[, which( is.na(USED_GENES) )]
	if(length(tmp))
		stop(paste('Found missing USED_GENES in rows', paste(tmp, collapse=',')))	
	#	check OBJ_i
	set(df, NULL, 'OBJ_i', df[, tolower(gsub('\\s','',as.character(OBJ_i)))])
	set(df, df[, which(OBJ_i=='declining')], 'OBJ_i', 'decreasing')
	tmp			<- df[, which( !is.na(OBJ_i) & !OBJ_i%in%c('stable','increasing','decreasing') )]
	if(length(tmp))
		stop(paste('Found invalid answer to OBJ_i', paste(df[tmp,OBJ_i], collapse=',')))
	tmp			<- df[, which( is.na(OBJ_i) )]
	if(length(tmp))
		cat(paste('\nFound missing answer to OBJ_i in rows', paste(tmp, collapse=',')))
	#
	#	check OBJ_iv and update OBJ_iv: consider <15% 15-30% >30% rather than <10% etc if still old coding
	#
	set(df, NULL, 'OBJ_iv', df[, tolower(gsub('\\s','',as.character(OBJ_iv)))])
	tmp			<- df[, which( OBJ_iv%in%c('<10%','10-30%','10%-30%') )]
	if(length(tmp))
	{
		if(verbose)
			cat(paste('Found old coding to OBJ_iv', paste(df[tmp,OBJ_iv], collapse=',')))
		tmp2		<- df[, which(is.na(OBJ_iv))]
		stopifnot( df[, !any(!is.na(OBJ_iv) & is.na(OBJ_v))] )
		set(df, df[, which(!is.na(OBJ_v))], 'OBJ_iv', NA_character_)
		#	if OBJ_iv missing and OBJ_v provided, set automatically
		tmp			<- df[, which(is.na(OBJ_iv) & !is.na(OBJ_v))]
		if(length(tmp))
			set(df, tmp, 'OBJ_iv', df[tmp, cut(OBJ_v, breaks=c(-Inf,0.15,0.3,Inf), labels=c('<15%','15%-30%','>30%'))])
		set(df, tmp2, 'OBJ_iv', NA_character_)		
	}
	if(!length(tmp))
		reset.OBJiv.conservative	<- 0		#nothing reset
	tmp			<- df[, which( !is.na(OBJ_iv) & !OBJ_iv%in%c('<15%','15-30%','15%-30%','>30%') )]	
	if(length(tmp))
		stop(paste('Found invalid answers to OBJ_iv', paste(df[tmp,OBJ_iv], collapse=',')))
	#	set OBJ_iv 15-30% to 15%-30%
	tmp			<- df[, which( OBJ_iv=='15-30%' )]
	set(df, tmp, 'OBJ_iv', '15%-30%')
	tmp			<- df[, which( is.na(OBJ_iv) )]
	if(length(tmp))
		cat(paste('\nFound missing answer to OBJ_iv in rows', paste(tmp, collapse=',')))			
	#	check OBJ_ii OBJ_v OBJ_vi
	for(x in c('OBJ_ii','OBJ_v','OBJ_vi'))
	{		
		tmp	<- which(df[[x]]>1 | df[[x]]<0)
		if(length(tmp))
			stop(paste('Expect value in [0-1]. Found invalid numerical entry for',x,'in rows',paste(tmp, collapse=',')))
		tmp	<- which(is.na(df[[x]]))
		if(verbose & length(tmp))
			cat(paste('\nFound missing value for',x,'in rows',paste(tmp, collapse=',')))
	}
	#	check OBJ_ii OBJ_iii OBJ_v OBJ_vi
	for(x in c('OBJ_iii'))
	{		
		tmp	<- which(df[[x]]<=0)
		if(length(tmp))
			stop(paste('Expect value in [0-1]. Found invalid numerical entry for',x,'in rows',paste(tmp, collapse=',')))
		tmp	<- which(is.na(df[[x]]))
		if(verbose & length(tmp))
			cat(paste('\nFound missing value for',x,'in rows',paste(tmp, collapse=',')))
	}
	#	check for duplicates
	tmp			<- df[, list(CH= length(OBJ_i)), by=c('TEAM','SIM_SCENARIO','SUBMISSION_DATE','USED_GENES','ESTIMATE')]
	tmp			<- subset(tmp, CH!=1)
	if(nrow(tmp))
	{
		print(tmp)
		stop('Found duplicate submissions')
	}
	#
	#	check consistency of OBJ_i and OBJ_iii
	#
	tmp			<- df[, which(OBJ_i=='stable' & abs(1-OBJ_iii)>0.1)]
	if(length(tmp))
		warning(paste('\nPlease check: Potential conflict between OBJ_i stable and OBJ_iii for scenario',paste(df[tmp,SIM_SCENARIO], collapse=', ')))
	tmp			<- df[, which(OBJ_i=='increasing' & OBJ_iii<1)]
	if(length(tmp))
		warning(paste('\nPlease check: Potential conflict between OBJ_i increasing and OBJ_iii for scenario',paste(df[tmp,SIM_SCENARIO], collapse=', ')))
	tmp			<- df[, which(OBJ_i=='decreasing' & OBJ_iii>1)]
	if(length(tmp))
		warning(paste('\nPlease check: Potential conflict between OBJ_i decreasing and OBJ_iii for scenario',paste(df[tmp,SIM_SCENARIO], collapse=', ')))	
	#	check ESTIMATE
	set(df, NULL, 'ESTIMATE', df[, tolower(gsub('\\s','',as.character(ESTIMATE)))])
	tmp			<- df[, which( !ESTIMATE%in%c('central','lower95%','upper95%') )]
	if(length(tmp))
		stop(paste('Expect either "central", "lower95%", "upper95%". Found invalid ESTIMATE in rows', paste(tmp, collapse=',')))
	tmp			<- df[, which( is.na(ESTIMATE) )]
	if(length(tmp))
		stop(paste('Found missing ESTIMATE in rows', paste(tmp, collapse=',')))		
	#
	#	check if at least one estimate provided per row
	#
	tmp			<- subset(df, ESTIMATE=='central')[, {
				list(CH= all(is.na(c(OBJ_i, OBJ_ii, OBJ_iii, OBJ_iv, OBJ_v, OBJ_vi))))
			}, by=c('TEAM','SIM_SCENARIO','SUBMISSION_DATE','USED_GENES')]
	tmp			<- tmp[, which(CH)]
	if(length(tmp))
		stop('Found rows with no submitted answer, ',paste(tmp, collapse=','))
	#
	#	remove NA rows confidence intervals
	#
	tmp			<- df[, {
				list(CH= all(is.na(c(OBJ_i, OBJ_ii, OBJ_iii, OBJ_iv, OBJ_v, OBJ_vi))))
			}, by=c('TEAM','SIM_SCENARIO','SUBMISSION_DATE','USED_GENES','ESTIMATE')]
	if(any(tmp[,CH]))
	{
		cat(paste('\nFound rows with no estimate. Removing rows',paste(tmp[,which(CH)], collapse=',')))
		df		<- merge(df, subset(tmp, CH==FALSE), by=c('TEAM','SIM_SCENARIO','SUBMISSION_DATE','USED_GENES','ESTIMATE'))
		df[, CH:=NULL]
	}
	#
	#	check if central estimate provided
	#
	tmp			<- df[, {
				lapply( .SD, function(x) all(is.na(ESTIMATE[ x ])) | 'central'%in%ESTIMATE[ !is.na(x) ] )			
			}, by=c('TEAM','SIM_SCENARIO','SUBMISSION_DATE','USED_GENES'), .SDcols=c('ESTIMATE','OBJ_i','OBJ_ii','OBJ_iii','OBJ_iv','OBJ_v','OBJ_vi')]
	for(x in c('OBJ_i','OBJ_ii','OBJ_iii','OBJ_iv','OBJ_v','OBJ_vi'))
		if(!all(tmp[[x]]))
			stop('Found scenarios with no central answer for',x)
	#
	#	check if lower and upper estimate when one is not missing
	#
	tmp			<- df[, {
				lapply( .SD, function(x) all(is.na(ESTIMATE[ x ])) | !'lower95%'%in%ESTIMATE[ !is.na(x) ]  |  ('lower95%'%in%ESTIMATE[ !is.na(x) ]&'upper95%'%in%ESTIMATE[ !is.na(x) ])	)			
			}, by=c('TEAM','SIM_SCENARIO','SUBMISSION_DATE','USED_GENES'), .SDcols=c('ESTIMATE','OBJ_i','OBJ_ii','OBJ_iii','OBJ_iv','OBJ_v','OBJ_vi')]
	for(x in c('OBJ_i','OBJ_ii','OBJ_iii','OBJ_iv','OBJ_v','OBJ_vi'))
		if(!all(tmp[[x]]))
			stop('Found scenarios with lower95% but no upper95%, or vice versa',x)
	#
	#	check that lower < upper estimate when not missing
	#
	tmp			<- df[, {
				lapply( .SD, function(x) all(is.na(ESTIMATE[ x ])) | !'lower95%'%in%ESTIMATE[ !is.na(x) ]  |  ('lower95%'%in%ESTIMATE[ !is.na(x) ]&'upper95%'%in%ESTIMATE[ !is.na(x) ]&x[ESTIMATE=='lower95%']<=x[ESTIMATE=='upper95%'])	)			
			}, by=c('TEAM','SIM_SCENARIO','SUBMISSION_DATE','USED_GENES'), .SDcols=c('ESTIMATE','OBJ_i','OBJ_ii','OBJ_iii','OBJ_iv','OBJ_v','OBJ_vi')]
	for(x in c('OBJ_ii','OBJ_iii','OBJ_v','OBJ_vi'))
		if(!all(tmp[[x]]))
			stop('Found scenarios with lower95% > upper95%',x)	
	#
	#	check first objective for regional simulations
	#
	df.sc		<- data.table(SIM_SCENARIO= paste('150129_PANGEAsim_Regional_FirstObj_sc',LETTERS[seq(1,4)],'_SIMULATED_SEQ',sep='') ) 
	tmp			<- unique( setdiff( df.sc[, SIM_SCENARIO], df[, SIM_SCENARIO] ) )
	if(warn.all & length(tmp))
		warning(paste('Found missing scenarios for first objective?',paste(tmp,collapse=',')))
	#
	#	check second objective for regional simulations
	#
	df.sc		<- data.table(SIM_SCENARIO= paste('150129_PANGEAsim_Regional_SecondObj_sc',LETTERS[seq(5,20)],'_SIMULATED_DATEDTREE',sep='') ) 
	tmp			<- unique( setdiff( df.sc[, SIM_SCENARIO], df[, SIM_SCENARIO] ) )
	if(warn.all & length(tmp) & length(tmp)<nrow(df.sc))
		warning(paste('Found only a subset of submissions for second objective? Not submitted:',paste(tmp,collapse=',')))	
	#
	cat('\nPassed checks.\n')
	#
	#	re-format
	#
	for(x in c('OBJ_ii','OBJ_iii','OBJ_v','OBJ_vi'))
	{		
		set(df, NULL, x, as.character(df[[x]]))
	}
	df	<- melt(df, measure.vars=which(grepl('OBJ',names(df))), variable.name="OBJ")
	set(df, NULL, 'ESTIMATE', df[, gsub('%','',ESTIMATE)])
	df	<- dcast.data.table(df, TEAM+SUBMISSION_DATE+SIM_SCENARIO+USED_GENES+OBJ~ESTIMATE,  value.var='value')
	df	<- subset( df, !is.na(central) )	
	#	be conservative on %Acute OBJ_i: set to NA if CIs overlap
	if(reset.OBJiv.conservative & all(c('upper95','lower95')%in%colnames(df)))
	{
		df	<- merge(df, subset(df, OBJ=='OBJ_v'), by=c('TEAM','SUBMISSION_DATE','SIM_SCENARIO','USED_GENES'), all.x=1)
		tmp	<- c( df[, which(OBJ.x=='OBJ_iv' & central.x=='<15%' & upper95.y>0.15)],
				df[, which(OBJ.x=='OBJ_iv' & central.x=='15%-30%' & (upper95.y>0.3 | lower95.y<0.15))],
				df[, which(OBJ.x=='OBJ_iv' & central.x=='>30%' & lower95.y<0.3)] )
		cat(paste('\nSetting OBJ_iv to NA because confidence intervals of OBJ_v overlap boundaries, n=', length(tmp)))  
		set(df, tmp, 'central.x', NA_character_)	
		set(df, NULL, c('OBJ.y','central.y','upper95.y','lower95.y'), NULL)
		setnames(df, c('OBJ.x','central.x','lower95.x','upper95.x'), c('OBJ','central','lower95','upper95'))	
	}
	setkey(df,  TEAM, SUBMISSION_DATE, SIM_SCENARIO, USED_GENES)
	cat(paste('\nFound submissions for unique scenarios, n=', df[, length(unique(SIM_SCENARIO))]))
	cat(paste('\nFound submissions for unique USED_GENES, n=', df[, length(unique(USED_GENES))]))
	cat(paste('\nFound submissions for unique SIM_SCENARIO x USED_GENES, n=', nrow(unique(df))))
	cat(paste('\nFound submissions for unique objectives, n=', df[, length(unique(OBJ))]))
	cat(paste('\nFound total estimates, n=', nrow(df)))
	if(!any('lower95'==names(df)))
		df[, lower95:=NA_real_]
	if(!any('upper95'==names(df)))
		df[, upper95:=NA_real_]	
	cat(paste('\nFound total estimates with confidence intervals, n=', nrow(subset(df, !is.na(lower95)))))
	#
	df
}

read.submission.Feb2015<- function(file, verbose=1, warn.all=0, reset.OBJiv.conservative=1)
{
	require(data.table)
	cat('\nThis is read.submission version 15-05-10.')
	if(verbose)
	{
		cat('\nReminder of objectives\n(use verbose=0 to suppress this message)')
		
		cat('\nOBJ_i\tDuring the evaluation period, was incidence stable, declining or increasing?\n\t\tAnswer: "stable", "declining", "increasing", or "NA"')
		cat('\nOBJ_ii\tWhat is the annual % incidence in the last year of the evaluation period?\n\t\tAnswer: numerical or "NA"')
		cat('\nOBJ_iii\tComparing the last year of the evaluation period to the year preceding the evaluation period, what is the ratio in annual % incidence?\n\t\tAnswer: numerical or "NA"')
		
		cat('\nOBJ_iv\tWas the proportion of transmissions that originated from individuals in early HIV infection in the year preceding the evaluation period below 10%, between 10-30%, or above 30%?\n\t\tAnswer: "<10%", "10-30%", ">30%" or "NA"')
		cat('\nOBJ_v\tWhat is the proportion of transmissions that originated from individuals in early HIV infection in the year preceding the evaluation period?\n\t\tAnswer: numerical or "NA"')
		cat('\nOBJ_vi\tWhat is the proportion of transmissions that originated from individuals in early HIV infection in the last year of the evaluation period?\n\t\tAnswer: numerical or "NA"')
		
		cat('\n\nPlease use:\nscenario names SIM_SCENARIO that correspond to the available file names\nUSED_GENES\teither "pol" or "all"')
	}
	cat('\nSkip first 22 rows that should contain comment rows starting with #.')
	cat(paste('\nreading', file))
	df	<- read.csv(file, stringsAsFactors=FALSE, comment.char="#", blank.lines.skip=TRUE, skip=22)
	df	<- as.data.table(df)	
	#	check column names
	df.colnm	<- c( "TEAM", "SUBMISSION_DATE", "SIM_SCENARIO", "USED_GENES", "OBJ_i", "OBJ_ii", "OBJ_iii", "OBJ_iv", "OBJ_v", "OBJ_vi", "ESTIMATE" )
	tmp			<- setdiff( df.colnm, names(df) )
	if(length(tmp))
		stop(paste('Found missing columns, ', paste(tmp, collapse=','), sep=''))
	tmp			<- setdiff( names(df), df.colnm )
	if(length(tmp))
		warning(paste('Ignore extra columns, ', paste(tmp, collapse=','), sep=''))
	df			<- df[, df.colnm, with=0]
	#
	#	check columns
	#
	#	check TEAM
	if(df[, length(unique(TEAM))>1])
		stop('Found more than one TEAM name.')
	#	check SUBMISSION_DATE
	if(df[, length(unique(SUBMISSION_DATE))>1])
		stop('Found more than one SUBMISSION_DATE.')
	#	check SIM_SCENARIO
	df.sc		<- data.table(SIM_SCENARIO= c( 	paste('Vill_',sprintf("%02d",0:12),'_Feb2015_5yr',sep=''), 
					paste('Vill_',sprintf("%02d",c(0:1,4,6:7,9:12)),'_Feb2015',sep=''),
					paste('Vill_',sprintf("%02d",c(2:3,5,8)),'_Feb2015_3yr',sep=''),
					paste('150129_PANGEAsim_Regional_FirstObj_sc',LETTERS[seq(1,4)],'_SIMULATED_SEQ',sep=''),
					paste('150129_PANGEAsim_Regional_SecondObj_sc',LETTERS[seq(5,20)],'_SIMULATED_DATEDTREE',sep='')												
			)) 
	tmp			<- unique( setdiff( df[, SIM_SCENARIO], df.sc[, SIM_SCENARIO] ) )
	if(length(tmp))
		stop(paste('Found invalid scenarios',paste(tmp,collapse=',')))
	df			<- merge(df, df.sc, by='SIM_SCENARIO')		
	#	check USED_GENES
	set(df, NULL, 'USED_GENES', df[, tolower(gsub('\\s','',as.character(USED_GENES)))])
	tmp			<- df[, which( !USED_GENES%in%c('pol','all') )]
	if(length(tmp))
		stop(paste('Found invalid USED_GENES', paste(df[tmp,USED_GENES], collapse=',')))
	tmp			<- df[, which( is.na(USED_GENES) )]
	if(length(tmp))
		stop(paste('Found missing USED_GENES in rows', paste(tmp, collapse=',')))	
	#	check OBJ_i
	set(df, NULL, 'OBJ_i', df[, tolower(gsub('\\s','',as.character(OBJ_i)))])
	set(df, df[, which(OBJ_i=='declining')], 'OBJ_i', 'decreasing')
	tmp			<- df[, which( !is.na(OBJ_i) & !OBJ_i%in%c('stable','increasing','decreasing') )]
	if(length(tmp))
		stop(paste('Found invalid answer to OBJ_i', paste(df[tmp,OBJ_i], collapse=',')))
	tmp			<- df[, which( is.na(OBJ_i) )]
	if(length(tmp))
		cat(paste('\nFound missing answer to OBJ_i in rows', paste(tmp, collapse=',')))		
	#	check OBJ_iv
	set(df, NULL, 'OBJ_iv', df[, tolower(gsub('\\s','',as.character(OBJ_iv)))])
	tmp			<- df[, which( !is.na(OBJ_iv) & !OBJ_iv%in%c('<10%','10-30%','10%-30%','>30%') )]
	if(length(tmp))
		stop(paste('Found invalid answer to OBJ_iv', paste(df[tmp,OBJ_iv], collapse=',')))
	tmp			<- df[, which( OBJ_iv=='10-30%' )]
	set(df, tmp, 'OBJ_iv', '10%-30%')
	tmp			<- df[, which( is.na(OBJ_iv) )]
	if(length(tmp))
		cat(paste('\nFound missing answer to OBJ_iv in rows', paste(tmp, collapse=',')))			
	#	check OBJ_ii OBJ_v OBJ_vi
	for(x in c('OBJ_ii','OBJ_v','OBJ_vi'))
	{		
		tmp	<- which(df[[x]]>1 | df[[x]]<0)
		if(length(tmp))
			stop(paste('Expect value in [0-1]. Found invalid numerical entry for',x,'in rows',paste(tmp, collapse=',')))
		tmp	<- which(is.na(df[[x]]))
		if(verbose & length(tmp))
			cat(paste('\nFound missing value for',x,'in rows',paste(tmp, collapse=',')))
	}
	#	check OBJ_ii OBJ_iii OBJ_v OBJ_vi
	for(x in c('OBJ_iii'))
	{		
		tmp	<- which(df[[x]]<=0)
		if(length(tmp))
			stop(paste('Expect value in [0-1]. Found invalid numerical entry for',x,'in rows',paste(tmp, collapse=',')))
		tmp	<- which(is.na(df[[x]]))
		if(verbose & length(tmp))
			cat(paste('\nFound missing value for',x,'in rows',paste(tmp, collapse=',')))
	}
	#	check for duplicates
	tmp			<- df[, list(CH= length(OBJ_i)), by=c('TEAM','SIM_SCENARIO','SUBMISSION_DATE','USED_GENES','ESTIMATE')]
	tmp			<- subset(tmp, CH!=1)
	if(nrow(tmp))
	{
		print(tmp)
		stop('Found duplicate submissions')
	}
	#
	#	check consistency of OBJ_i and OBJ_iii
	#
	tmp			<- df[, which(OBJ_i=='stable' & abs(1-OBJ_iii)>0.1)]
	if(length(tmp))
		warning(paste('\nPlease check: Potential conflict between OBJ_i stable and OBJ_iii for scenario',paste(df[tmp,SIM_SCENARIO], collapse=', ')))
	tmp			<- df[, which(OBJ_i=='increasing' & OBJ_iii<1)]
	if(length(tmp))
		warning(paste('\nPlease check: Potential conflict between OBJ_i increasing and OBJ_iii for scenario',paste(df[tmp,SIM_SCENARIO], collapse=', ')))
	tmp			<- df[, which(OBJ_i=='decreasing' & OBJ_iii>1)]
	if(length(tmp))
		warning(paste('\nPlease check: Potential conflict between OBJ_i decreasing and OBJ_iii for scenario',paste(df[tmp,SIM_SCENARIO], collapse=', ')))	
	#
	#	updated OBJ_iv: consider <15% 15-30% >30%
	#
	tmp2		<- df[, which(is.na(OBJ_iv))]
	stopifnot( df[, !any(!is.na(OBJ_iv) & is.na(OBJ_v))] )
	set(df, df[, which(!is.na(OBJ_v))], 'OBJ_iv', NA_character_)
	#	if OBJ_iv missing and OBJ_v provided, set automatically
	tmp			<- df[, which(is.na(OBJ_iv) & !is.na(OBJ_v))]
	if(length(tmp))
		set(df, tmp, 'OBJ_iv', df[tmp, cut(OBJ_v, breaks=c(-Inf,0.15,0.3,Inf), labels=c('<15%','15%-30%','>30%'))])
	set(df, tmp2, 'OBJ_iv', NA_character_)
	#	check ESTIMATE
	set(df, NULL, 'ESTIMATE', df[, tolower(gsub('\\s','',as.character(ESTIMATE)))])
	tmp			<- df[, which( !ESTIMATE%in%c('central','lower95%','upper95%') )]
	if(length(tmp))
		stop(paste('Expect either "central", "lower95%", "upper95%". Found invalid ESTIMATE in rows', paste(tmp, collapse=',')))
	tmp			<- df[, which( is.na(ESTIMATE) )]
	if(length(tmp))
		stop(paste('Found missing ESTIMATE in rows', paste(tmp, collapse=',')))		
	#
	#	check if at least one estimate provided per row
	#
	tmp			<- subset(df, ESTIMATE=='central')[, {
				list(CH= all(is.na(c(OBJ_i, OBJ_ii, OBJ_iii, OBJ_iv, OBJ_v, OBJ_vi))))
			}, by=c('TEAM','SIM_SCENARIO','SUBMISSION_DATE','USED_GENES')]
	tmp			<- tmp[, which(CH)]
	if(length(tmp))
		stop('Found rows with no submitted answer, ',paste(tmp, collapse=','))
	#
	#	remove NA rows confidence intervals
	#
	tmp			<- df[, {
				list(CH= all(is.na(c(OBJ_i, OBJ_ii, OBJ_iii, OBJ_iv, OBJ_v, OBJ_vi))))
			}, by=c('TEAM','SIM_SCENARIO','SUBMISSION_DATE','USED_GENES','ESTIMATE')]
	if(any(tmp[,CH]))
	{
		cat(paste('\nFound rows with no estimate. Removing rows',paste(tmp[,which(CH)], collapse=',')))
		df		<- merge(df, subset(tmp, CH==FALSE), by=c('TEAM','SIM_SCENARIO','SUBMISSION_DATE','USED_GENES','ESTIMATE'))
		df[, CH:=NULL]
	}
	#
	#	check if central estimate provided
	#
	tmp			<- df[, {
				lapply( .SD, function(x) all(is.na(ESTIMATE[ x ])) | 'central'%in%ESTIMATE[ !is.na(x) ] )			
			}, by=c('TEAM','SIM_SCENARIO','SUBMISSION_DATE','USED_GENES'), .SDcols=c('ESTIMATE','OBJ_i','OBJ_ii','OBJ_iii','OBJ_iv','OBJ_v','OBJ_vi')]
	for(x in c('OBJ_i','OBJ_ii','OBJ_iii','OBJ_iv','OBJ_v','OBJ_vi'))
		if(!all(tmp[[x]]))
			stop('Found scenarios with no central answer for',x)
	#
	#	check if lower and upper estimate when one is not missing
	#
	tmp			<- df[, {
				lapply( .SD, function(x) all(is.na(ESTIMATE[ x ])) | !'lower95%'%in%ESTIMATE[ !is.na(x) ]  |  ('lower95%'%in%ESTIMATE[ !is.na(x) ]&'upper95%'%in%ESTIMATE[ !is.na(x) ])	)			
			}, by=c('TEAM','SIM_SCENARIO','SUBMISSION_DATE','USED_GENES'), .SDcols=c('ESTIMATE','OBJ_i','OBJ_ii','OBJ_iii','OBJ_iv','OBJ_v','OBJ_vi')]
	for(x in c('OBJ_i','OBJ_ii','OBJ_iii','OBJ_iv','OBJ_v','OBJ_vi'))
		if(!all(tmp[[x]]))
			stop('Found scenarios with lower95% but no upper95%, or vice versa',x)
	#
	#	check that lower < upper estimate when not missing
	#
	tmp			<- df[, {
				lapply( .SD, function(x) all(is.na(ESTIMATE[ x ])) | !'lower95%'%in%ESTIMATE[ !is.na(x) ]  |  ('lower95%'%in%ESTIMATE[ !is.na(x) ]&'upper95%'%in%ESTIMATE[ !is.na(x) ]&x[ESTIMATE=='lower95%']<=x[ESTIMATE=='upper95%'])	)			
			}, by=c('TEAM','SIM_SCENARIO','SUBMISSION_DATE','USED_GENES'), .SDcols=c('ESTIMATE','OBJ_i','OBJ_ii','OBJ_iii','OBJ_iv','OBJ_v','OBJ_vi')]
	for(x in c('OBJ_ii','OBJ_iii','OBJ_v','OBJ_vi'))
		if(!all(tmp[[x]]))
			stop('Found scenarios with lower95% > upper95%',x)	
	#
	#	check first objective for regional simulations
	#
	df.sc		<- data.table(SIM_SCENARIO= paste('150129_PANGEAsim_Regional_FirstObj_sc',LETTERS[seq(1,4)],'_SIMULATED_SEQ',sep='') ) 
	tmp			<- unique( setdiff( df.sc[, SIM_SCENARIO], df[, SIM_SCENARIO] ) )
	if(warn.all & length(tmp))
		warning(paste('Found missing scenarios for first objective?',paste(tmp,collapse=',')))
	#
	#	check second objective for regional simulations
	#
	df.sc		<- data.table(SIM_SCENARIO= paste('150129_PANGEAsim_Regional_SecondObj_sc',LETTERS[seq(5,20)],'_SIMULATED_DATEDTREE',sep='') ) 
	tmp			<- unique( setdiff( df.sc[, SIM_SCENARIO], df[, SIM_SCENARIO] ) )
	if(warn.all & length(tmp) & length(tmp)<nrow(df.sc))
		warning(paste('Found only a subset of submissions for second objective? Not submitted:',paste(tmp,collapse=',')))	
	#
	cat('\nPassed checks.\n')
	#
	#	re-format
	#
	for(x in c('OBJ_ii','OBJ_iii','OBJ_v','OBJ_vi'))
	{		
		set(df, NULL, x, as.character(df[[x]]))
	}
	df	<- melt(df, measure.vars=which(grepl('OBJ',names(df))), variable.name="OBJ")
	set(df, NULL, 'ESTIMATE', df[, gsub('%','',ESTIMATE)])
	df	<- dcast.data.table(df, TEAM+SUBMISSION_DATE+SIM_SCENARIO+USED_GENES+OBJ~ESTIMATE,  value.var='value')
	df	<- subset( df, !is.na(central) )	
	#	be conservative on %Acute OBJ_i: set to NA if CIs overlap
	if(reset.OBJiv.conservative & all(c('upper95','lower95')%in%colnames(df)))
	{
		df	<- merge(df, subset(df, OBJ=='OBJ_v'), by=c('TEAM','SUBMISSION_DATE','SIM_SCENARIO','USED_GENES'), all.x=1)
		tmp	<- c( df[, which(OBJ.x=='OBJ_iv' & central.x=='<15%' & upper95.y>0.15)],
				df[, which(OBJ.x=='OBJ_iv' & central.x=='15%-30%' & (upper95.y>0.3 | lower95.y<0.15))],
				df[, which(OBJ.x=='OBJ_iv' & central.x=='>30%' & lower95.y<0.3)] )
		cat(paste('\nSetting OBJ_iv to NA because confidence intervals of OBJ_v overlap boundaries, n=', length(tmp)))  
		set(df, tmp, 'central.x', NA_character_)	
		set(df, NULL, c('OBJ.y','central.y','upper95.y','lower95.y'), NULL)
		setnames(df, c('OBJ.x','central.x','lower95.x','upper95.x'), c('OBJ','central','lower95','upper95'))	
	}
	setkey(df,  TEAM, SUBMISSION_DATE, SIM_SCENARIO, USED_GENES)
	cat(paste('\nFound submissions for unique scenarios, n=', df[, length(unique(SIM_SCENARIO))]))
	cat(paste('\nFound submissions for unique USED_GENES, n=', df[, length(unique(USED_GENES))]))
	cat(paste('\nFound submissions for unique SIM_SCENARIO x USED_GENES, n=', nrow(unique(df))))
	cat(paste('\nFound submissions for unique objectives, n=', df[, length(unique(OBJ))]))
	cat(paste('\nFound total estimates, n=', nrow(df)))
	if(!any('lower95'==names(df)))
		df[, lower95:=NA_real_]
	if(!any('upper95'==names(df)))
		df[, upper95:=NA_real_]	
	cat(paste('\nFound total estimates with confidence intervals, n=', nrow(subset(df, !is.na(lower95)))))
	#
	df
}

read.submission.May2015<- function(file, verbose=1, warn.all=0)
{
	require(data.table)
	cat('\nThis is read.submission version 15-08-12.')
	if(verbose)
	{
		cat('\nReminder of objectives\n(use verbose=0 to suppress this message)')
		
		cat('\nOBJ_i\tDuring the evaluation period, was incidence stable, declining or increasing?\n\t\tAnswer: "stable", "declining", "increasing", or "NA"')
		cat('\nOBJ_ii\tWhat is the annual % incidence in the last year of the evaluation period?\n\t\tAnswer: numerical or "NA"')
		cat('\nOBJ_iii\tComparing the last year of the evaluation period to the year preceding the evaluation period, what is the ratio in annual % incidence?\n\t\tAnswer: numerical or "NA"')
		
		cat('\nOBJ_iv\tWas the proportion of transmissions that originated from individuals in early HIV infection in the year preceding the evaluation period below 10%, between 10-30%, or above 30%?\n\t\tAnswer: "<10%", "10-30%", ">30%" or "NA"')
		cat('\nOBJ_v\tWhat is the proportion of transmissions that originated from individuals in early HIV infection in the year preceding the evaluation period?\n\t\tAnswer: numerical or "NA"')
		cat('\nOBJ_vi\tWhat is the proportion of transmissions that originated from individuals in early HIV infection in the last year of the evaluation period?\n\t\tAnswer: numerical or "NA"')
		
		cat('\n\nPlease use:\nscenario names SIM_SCENARIO that correspond to the available file names\nUSED_GENES\teither "pol" or "all"')
		cat('\nSkip first 22 rows that should contain comment rows starting with #.')
	}	
	cat(paste('\nreading', file))
	df	<- read.csv(file, stringsAsFactors=FALSE, comment.char="#", blank.lines.skip=TRUE, skip=22)
	df	<- as.data.table(df)	
	#	check column names
	df.colnm	<- c( "TEAM", "SUBMISSION_DATE", "SIM_SCENARIO", "USED_GENES", "OBJ_i", "OBJ_ii", "OBJ_iii", "OBJ_iv_3m", "OBJ_v_3m", "OBJ_vi_3m", "OBJ_iv_6m", "OBJ_v_6m", "OBJ_vi_6m", "OBJ_iv_12m", "OBJ_v_12m", "OBJ_vi_12m", "ESTIMATE" )
	tmp			<- setdiff( df.colnm, names(df) )
	if(length(tmp))
		stop(paste('Found missing columns, ', paste(tmp, collapse=','), sep=''))
	tmp			<- setdiff( names(df), df.colnm )
	if(length(tmp))
		warning(paste('Ignore extra columns, ', paste(tmp, collapse=','), sep=''))
	df			<- df[, df.colnm, with=0]
	#
	#	check columns
	#
	#	check TEAM
	if(df[, length(unique(TEAM))>1])
		stop('Found more than one TEAM name.')
	#	check SUBMISSION_DATE
	if(df[, length(unique(SUBMISSION_DATE))>1])
		stop('Found more than one SUBMISSION_DATE.')
	#	check SIM_SCENARIO
	df.sc		<- data.table(SIM_SCENARIO= c( 	paste('Vill_',sprintf("%02d",0:12),'_Feb2015_5yr',sep=''), 
												paste('Vill_',sprintf("%02d",c(0:1,4,6:7,9:12)),'_Feb2015',sep=''),
												paste('Vill_',sprintf("%02d",c(2:3,5,8)),'_Feb2015_3yr',sep=''),
												paste('150129_PANGEAsim_Regional_FirstObj_sc',LETTERS[seq(1,4)],'_SIMULATED_SEQ',sep=''),
												paste('150129_PANGEAsim_Regional_SecondObj_sc',LETTERS[seq(5,20)],'_SIMULATED_DATEDTREE',sep='')												
												)) 
	tmp			<- unique( setdiff( df[, SIM_SCENARIO], df.sc[, SIM_SCENARIO] ) )
	if(length(tmp))
		stop(paste('Found invalid scenarios',paste(tmp,collapse=',')))
	df			<- merge(df, df.sc, by='SIM_SCENARIO')		
	#	check USED_GENES
	set(df, NULL, 'USED_GENES', df[, tolower(gsub('\\s','',as.character(USED_GENES)))])
	tmp			<- df[, which( !USED_GENES%in%c('pol','all') )]
	if(length(tmp))
		stop(paste('Found invalid USED_GENES', paste(df[tmp,USED_GENES], collapse=',')))
	tmp			<- df[, which( is.na(USED_GENES) )]
	if(length(tmp))
		stop(paste('Found missing USED_GENES in rows', paste(tmp, collapse=',')))	
	#	check OBJ_i
	set(df, NULL, 'OBJ_i', df[, tolower(gsub('\\s','',as.character(OBJ_i)))])
	set(df, df[, which(OBJ_i=='declining')], 'OBJ_i', 'decreasing')
	tmp			<- df[, which( !is.na(OBJ_i) & !OBJ_i%in%c('stable','increasing','decreasing') )]
	if(length(tmp))
		stop(paste('Found invalid answer to OBJ_i', paste(df[tmp,OBJ_i], collapse=',')))
	tmp			<- df[, which( is.na(OBJ_i) )]
	if(length(tmp))
		cat(paste('\nFound missing answer to OBJ_i in rows', paste(tmp, collapse=',')))		
	#	check OBJ_iv_3m
	set(df, NULL, 'OBJ_iv_3m', df[, tolower(gsub('\\s','',as.character(OBJ_iv_3m)))])
	tmp			<- df[, which( !is.na(OBJ_iv_3m) & !OBJ_iv_3m%in%c('<15%','15-30%','15%-30%','>30%') )]
	if(length(tmp))
		stop(paste('Found invalid answer to OBJ_iv_3m', paste(df[tmp,OBJ_iv_3m], collapse=',')))
	#	check OBJ_iv_6m
	set(df, NULL, 'OBJ_iv_6m', df[, tolower(gsub('\\s','',as.character(OBJ_iv_6m)))])
	tmp			<- df[, which( !is.na(OBJ_iv_6m) & !OBJ_iv_6m%in%c('<15%','15-30%','15%-30%','>30%') )]
	if(length(tmp))
		stop(paste('Found invalid answer to OBJ_iv_6m', paste(df[tmp,OBJ_iv_6m], collapse=',')))
	#	check OBJ_iv_12m
	set(df, NULL, 'OBJ_iv_12m', df[, tolower(gsub('\\s','',as.character(OBJ_iv_12m)))])
	tmp			<- df[, which( !is.na(OBJ_iv_12m) & !OBJ_iv_12m%in%c('<15%','15-30%','15%-30%','>30%') )]
	if(length(tmp))
		stop(paste('Found invalid answer to OBJ_iv_12m', paste(df[tmp,OBJ_iv_12m], collapse=',')))
	#	
	set(df, df[, which( OBJ_iv_3m=='15-30%' )], 'OBJ_iv_3m', '15%-30%')
	set(df, df[, which( OBJ_iv_6m=='15-30%' )], 'OBJ_iv_6m', '15%-30%')
	set(df, df[, which( OBJ_iv_12m=='15-30%' )], 'OBJ_iv_12m', '15%-30%')
	tmp			<- df[, which( is.na(OBJ_iv_3m) )]
	if(length(tmp))
		cat(paste('\nFound missing answer to OBJ_iv_3m in rows', paste(tmp, collapse=',')))
	tmp			<- df[, which( is.na(OBJ_iv_6m) )]
	if(length(tmp))
		cat(paste('\nFound missing answer to OBJ_iv_6m in rows', paste(tmp, collapse=',')))			
	tmp			<- df[, which( is.na(OBJ_iv_12m) )]
	if(length(tmp))
		cat(paste('\nFound missing answer to OBJ_iv_12m in rows', paste(tmp, collapse=',')))			
	#	check OBJ_ii OBJ_v_Xm OBJ_vi_Xm
	for(x in c('OBJ_ii','OBJ_v_3m','OBJ_vi_3m','OBJ_v_6m','OBJ_vi_6m','OBJ_v_12m','OBJ_vi_12m'))
	{		
		tmp	<- which(df[[x]]>1 | df[[x]]<0)
		if(length(tmp))
			stop(paste('Expect value in [0-1]. Found invalid numerical entry for',x,'in rows',paste(tmp, collapse=',')))
		tmp	<- which(is.na(df[[x]]))
		if(verbose & length(tmp))
			cat(paste('\nFound missing value for',x,'in rows',paste(tmp, collapse=',')))
	}
	#	check OBJ_ii OBJ_iii OBJ_v OBJ_vi
	for(x in c('OBJ_iii'))
	{		
		tmp	<- which(df[[x]]<=0)
		if(length(tmp))
			stop(paste('Expect value in [0-1]. Found invalid numerical entry for',x,'in rows',paste(tmp, collapse=',')))
		tmp	<- which(is.na(df[[x]]))
		if(verbose & length(tmp))
			cat(paste('\nFound missing value for',x,'in rows',paste(tmp, collapse=',')))
	}
	#	check for duplicates
	tmp			<- df[, list(CH= length(OBJ_i)), by=c('TEAM','SIM_SCENARIO','SUBMISSION_DATE','USED_GENES','ESTIMATE')]
	tmp			<- subset(tmp, CH!=1)
	if(nrow(tmp))
	{
		print(tmp)
		stop('Found duplicate submissions')
	}
	#
	#	check consistency of OBJ_i and OBJ_iii
	#
	tmp			<- df[, which(OBJ_i=='stable' & abs(1-OBJ_iii)>0.1)]
	if(length(tmp))
		warning(paste('\nPlease check: Potential conflict between OBJ_i stable and OBJ_iii for scenario',paste(df[tmp,SIM_SCENARIO], collapse=', ')))
	tmp			<- df[, which(OBJ_i=='increasing' & OBJ_iii<1)]
	if(length(tmp))
		warning(paste('\nPlease check: Potential conflict between OBJ_i increasing and OBJ_iii for scenario',paste(df[tmp,SIM_SCENARIO], collapse=', ')))
	tmp			<- df[, which(OBJ_i=='decreasing' & OBJ_iii>1)]
	if(length(tmp))
		warning(paste('\nPlease check: Potential conflict between OBJ_i decreasing and OBJ_iii for scenario',paste(df[tmp,SIM_SCENARIO], collapse=', ')))	
	#
	#	updated OBJ_iv: consider <15% 15-30% >30%
	#
	#tmp2		<- df[, which(is.na(OBJ_iv))]
	#stopifnot( df[, !any(!is.na(OBJ_iv) & is.na(OBJ_v))] )
	#set(df, df[, which(!is.na(OBJ_v))], 'OBJ_iv', NA_character_)
	#	if OBJ_iv missing and OBJ_v provided, set automatically
	#tmp			<- df[, which(is.na(OBJ_iv) & !is.na(OBJ_v))]
	#if(length(tmp))
	#	set(df, tmp, 'OBJ_iv', df[tmp, cut(OBJ_v, breaks=c(-Inf,0.15,0.3,Inf), labels=c('<15%','15%-30%','>30%'))])
	#set(df, tmp2, 'OBJ_iv', NA_character_)
	#	check ESTIMATE
	set(df, NULL, 'ESTIMATE', df[, tolower(gsub('\\s','',as.character(ESTIMATE)))])
	tmp			<- df[, which( !ESTIMATE%in%c('central','lower95%','upper95%') )]
	if(length(tmp))
		stop(paste('Expect either "central", "lower95%", "upper95%". Found invalid ESTIMATE in rows', paste(tmp, collapse=',')))
	tmp			<- df[, which( is.na(ESTIMATE) )]
	if(length(tmp))
		stop(paste('Found missing ESTIMATE in rows', paste(tmp, collapse=',')))		
	#
	#	check if at least one estimate provided per row
	#
	tmp			<- subset(df, ESTIMATE=='central')[, {
							list(CH= all(is.na(c(OBJ_i, OBJ_ii, OBJ_iii, OBJ_iv_3m, OBJ_v_3m, OBJ_vi_3m, OBJ_iv_6m, OBJ_v_6m, OBJ_vi_6m, OBJ_iv_12m, OBJ_v_12m, OBJ_vi_12m))))
						}, by=c('TEAM','SIM_SCENARIO','SUBMISSION_DATE','USED_GENES')]
	tmp			<- tmp[, which(CH)]
	if(length(tmp))
		stop('Found rows with no submitted answer, ',paste(tmp, collapse=','))
	#
	#	remove NA rows confidence intervals
	#
	tmp			<- df[, {
							list(CH= all(is.na(c(OBJ_i, OBJ_ii, OBJ_iii, OBJ_iv_3m, OBJ_v_3m, OBJ_vi_3m, OBJ_iv_6m, OBJ_v_6m, OBJ_vi_6m, OBJ_iv_12m, OBJ_v_12m, OBJ_vi_12m))))
						}, by=c('TEAM','SIM_SCENARIO','SUBMISSION_DATE','USED_GENES','ESTIMATE')]
	if(any(tmp[,CH]))
	{
		cat(paste('\nFound rows with no estimate. Removing rows',paste(tmp[,which(CH)], collapse=',')))
		df		<- merge(df, subset(tmp, CH==FALSE), by=c('TEAM','SIM_SCENARIO','SUBMISSION_DATE','USED_GENES','ESTIMATE'))
		df[, CH:=NULL]
	}
	#
	#	check if central estimate provided
	#
	tmp			<- df[, {
				lapply( .SD, function(x) all(is.na(ESTIMATE[ x ])) | 'central'%in%ESTIMATE[ !is.na(x) ] )			
			}, by=c('TEAM','SIM_SCENARIO','SUBMISSION_DATE','USED_GENES'), .SDcols=c('ESTIMATE','OBJ_i','OBJ_ii','OBJ_iii','OBJ_iv_3m','OBJ_v_3m','OBJ_vi_3m','OBJ_iv_6m','OBJ_v_6m','OBJ_vi_6m','OBJ_iv_12m','OBJ_v_12m','OBJ_vi_12m')]
	for(x in c('OBJ_i','OBJ_ii','OBJ_iii','OBJ_iv_3m','OBJ_v_3m','OBJ_vi_3m','OBJ_iv_6m','OBJ_v_6m','OBJ_vi_6m','OBJ_iv_12m','OBJ_v_12m','OBJ_vi_12m'))
		if(!all(tmp[[x]]))
			stop('Found scenarios with no central answer for',x)
	#
	#	check if lower and upper estimate when one is not missing
	#
	tmp			<- df[, {
				lapply( .SD, function(x) all(is.na(ESTIMATE[ x ])) | !'lower95%'%in%ESTIMATE[ !is.na(x) ]  |  ('lower95%'%in%ESTIMATE[ !is.na(x) ]&'upper95%'%in%ESTIMATE[ !is.na(x) ])	)			
			}, by=c('TEAM','SIM_SCENARIO','SUBMISSION_DATE','USED_GENES'), .SDcols=c('ESTIMATE','OBJ_i','OBJ_ii','OBJ_iii','OBJ_iv_3m','OBJ_v_3m','OBJ_vi_3m','OBJ_iv_6m','OBJ_v_6m','OBJ_vi_6m','OBJ_iv_12m','OBJ_v_12m','OBJ_vi_12m')]
	for(x in c('OBJ_i','OBJ_ii','OBJ_iii','OBJ_iv_3m','OBJ_v_3m','OBJ_vi_3m','OBJ_iv_6m','OBJ_v_6m','OBJ_vi_6m','OBJ_iv_12m','OBJ_v_12m','OBJ_vi_12m'))
		if(!all(tmp[[x]]))
			stop('Found scenarios with lower95% but no upper95%, or vice versa',x)
	#
	#	check that lower < upper estimate when not missing
	#
	tmp			<- df[, {
				lapply( .SD, function(x) all(is.na(ESTIMATE[ x ])) | !'lower95%'%in%ESTIMATE[ !is.na(x) ]  |  ('lower95%'%in%ESTIMATE[ !is.na(x) ]&'upper95%'%in%ESTIMATE[ !is.na(x) ]&x[ESTIMATE=='lower95%']<=x[ESTIMATE=='upper95%'])	)			
			}, by=c('TEAM','SIM_SCENARIO','SUBMISSION_DATE','USED_GENES'), .SDcols=c('ESTIMATE','OBJ_i','OBJ_ii','OBJ_iii','OBJ_iv_3m','OBJ_v_3m','OBJ_vi_3m','OBJ_iv_6m','OBJ_v_6m','OBJ_vi_6m','OBJ_iv_12m','OBJ_v_12m','OBJ_vi_12m')]
	for(x in c('OBJ_i','OBJ_ii','OBJ_iii','OBJ_iv_3m','OBJ_v_3m','OBJ_vi_3m','OBJ_iv_6m','OBJ_v_6m','OBJ_vi_6m','OBJ_iv_12m','OBJ_v_12m','OBJ_vi_12m'))
		if(!all(tmp[[x]]))
			stop('Found scenarios with lower95% > upper95%',x)
	#
	#	check that 3m<6m
	#
	tmp			<- df[, which(!is.na(OBJ_v_3m) & !is.na(OBJ_v_6m) & ESTIMATE=='central')]
	if(length(tmp))
		invisible(df[tmp,][,{
					if(OBJ_v_3m>OBJ_v_6m)
						stop( '\nFound OBJ_v_3m>OBJ_v_6m', paste(TEAM, SIM_SCENARIO, SUBMISSION_DATE, USED_GENES, collapse=''))
				}, by=c('TEAM','SIM_SCENARIO','SUBMISSION_DATE','USED_GENES')])
	#
	#	check that 6m<12m
	#
	tmp			<- df[, which(!is.na(OBJ_v_6m) & !is.na(OBJ_v_12m) & ESTIMATE=='central')]
	if(length(tmp))
		invisible(df[tmp,][,{
							if(OBJ_v_6m>OBJ_v_12m)
								stop( '\nFound OBJ_v_6m>OBJ_v_12m', paste(TEAM, SIM_SCENARIO, SUBMISSION_DATE, USED_GENES, collapse=''))
						}, by=c('TEAM','SIM_SCENARIO','SUBMISSION_DATE','USED_GENES')])
	#
	#	check first objective for regional simulations
	#
	df.sc		<- data.table(SIM_SCENARIO= paste('150129_PANGEAsim_Regional_FirstObj_sc',LETTERS[seq(1,4)],'_SIMULATED_SEQ',sep='') ) 
	tmp			<- unique( setdiff( df.sc[, SIM_SCENARIO], df[, SIM_SCENARIO] ) )
	if(warn.all & length(tmp))
		warning(paste('Found missing scenarios for first objective?',paste(tmp,collapse=',')))
	#
	#	check second objective for regional simulations
	#
	df.sc		<- data.table(SIM_SCENARIO= paste('150129_PANGEAsim_Regional_SecondObj_sc',LETTERS[seq(5,20)],'_SIMULATED_DATEDTREE',sep='') ) 
	tmp			<- unique( setdiff( df.sc[, SIM_SCENARIO], df[, SIM_SCENARIO] ) )
	if(warn.all & length(tmp) & length(tmp)<nrow(df.sc))
		warning(paste('Found only a subset of submissions for second objective? Not submitted:',paste(tmp,collapse=',')))	
	#
	cat('\nPassed checks.\n')
	#
	#	re-format
	#
	for(x in c('OBJ_ii','OBJ_iii','OBJ_v_3m','OBJ_vi_3m','OBJ_v_6m','OBJ_vi_6m','OBJ_v_12m','OBJ_vi_12m'))
	{		
		set(df, NULL, x, as.character(df[[x]]))
	}
	#	set 'OBJ_v_3m' to 'OBJ_v' etc
	tmp	<- colnames(df)[ grepl('3m', colnames(df))]
	setnames(df, tmp, gsub('_*3m','',tmp))
	df	<- data.table:::melt.data.table(df, measure.vars=which(grepl('OBJ',names(df))), variable.name="OBJ")
	set(df, NULL, 'ESTIMATE', df[, gsub('%','',ESTIMATE)])
	df	<- dcast.data.table(df, TEAM+SUBMISSION_DATE+SIM_SCENARIO+USED_GENES+OBJ~ESTIMATE,  value.var='value')
	df	<- subset( df, !is.na(central) )	
	
	setkey(df,  TEAM, SUBMISSION_DATE, SIM_SCENARIO, USED_GENES)
	cat(paste('\nFound submissions for unique scenarios, n=', df[, length(unique(SIM_SCENARIO))]))
	cat(paste('\nFound submissions for unique USED_GENES, n=', df[, length(unique(USED_GENES))]))
	cat(paste('\nFound submissions for unique SIM_SCENARIO x USED_GENES, n=', nrow(unique(df))))
	cat(paste('\nFound submissions for unique objectives, n=', df[, length(unique(OBJ))]))
	cat(paste('\nFound total estimates, n=', nrow(df)))
	if(!any('lower95'==names(df)))
		df[, lower95:=NA_real_]
	if(!any('upper95'==names(df)))
		df[, upper95:=NA_real_]	
	cat(paste('\nFound total estimates with confidence intervals, n=', nrow(subset(df, !is.na(lower95)))))
	#
	df
}

