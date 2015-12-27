#	re-name the following:
outdir			<- '/Users/Oliver/duke/2015_various'
\dontrun{
##--------------------------------------------------------------------------------------------------------
##	first example  
##--------------------------------------------------------------------------------------------------------
#	input arguments for the simulation
pipeline.args	<- sim.regional.args( 	yr.start=1985, 
					yr.end=2020, 
					seed=42, 
					s.MODEL='Fixed2Prop', 
					report.prop.recent=1.0,
					s.PREV.max.n=1600, 
					s.INTERVENTION.prop=0.5, 
					s.INTERVENTION.start=2015, 
					s.INTERVENTION.mul=NA, 
					s.ARCHIVAL.n=50,
					epi.model='HPTN071', 
					epi.acute='high', 
					epi.intervention='fast', 
					epi.dt=1/48, 
					epi.import=0.05, 
					root.edge.fixed=0,
					v.N0tau=1, 
					v.r=2.851904, 
					v.T50=-2,
					wher.mu=log(0.00447743)-0.5^2/2, 
					wher.sigma=0.5, 
					bwerm.mu=log(0.002239075)-0.3^2/2, 
					bwerm.sigma=0.3, 
					er.gamma=4,
					dbg.GTRparam=0, 
					dbg.rER=0, 
					index.starttime.mode='fix1970', 
					startseq.mode='one', 
					seqtime.mode='AtDiag')									
#	produce UNIX script to generate the simulation
cat(sim.regional(outdir, pipeline.args=pipeline.args))
#	now run this script from the command line
	
##--------------------------------------------------------------------------------------------------------
##	The sequence data sets of the PANGEA-HIV methods comparison exercise (Primary Objectives)  
##--------------------------------------------------------------------------------------------------------
pipeline.args	<- sim.regional.args( 	yr.start=1985, yr.end=2020, seed=42, 
					s.MODEL='Fixed2Prop', report.prop.recent=1.0, 
					s.PREV.max.n=1600, s.INTERVENTION.prop=0.5, s.INTERVENTION.start=2015, s.INTERVENTION.mul= NA, s.ARCHIVAL.n=50,
					epi.model='HPTN071', epi.acute=NA, epi.intervention= NA, epi.dt=1/48, epi.import=0.05, root.edge.fixed=0,
					v.N0tau=1, v.r=2.851904, v.T50=-2,
					wher.mu=log(0.00447743)-0.5^2/2, wher.sigma=0.5, bwerm.mu=log(0.002239075)-0.3^2/2, bwerm.sigma=0.3, er.gamma=4,
					dbg.GTRparam=0, dbg.rER=0, index.starttime.mode='fix1970', startseq.mode='one', seqtime.mode='AtDiag')								
pipeline.vary	<- data.table(	label= c('D','C','A','B'),																										
				epi.acute= c('low','low','high','high'),
				epi.intervention= c('fast','slow','fast','slow')
								)		
invisible(pipeline.vary[, {									
			set(pipeline.args, which( pipeline.args$stat=='epi.acute' ), 'v', as.character(epi.acute))
			set(pipeline.args, which( pipeline.args$stat=='epi.intervention' ), 'v', as.character(epi.intervention))												
			tmpdir			<- paste(outdir,'-Dataset',label,sep='')
			dir.create(tmpdir, showWarnings=FALSE)																														
			file			<- sim.regional(tmpdir, pipeline.args=pipeline.args)
			#system(file)	
			}])
##--------------------------------------------------------------------------------------------------------
##	The tree data sets of the PANGEA-HIV methods comparison exercise (Secondary Objectives)  
##--------------------------------------------------------------------------------------------------------
pipeline.args	<- sim.regional.args( 	yr.start=1985, yr.end=NA, seed=NA, s.MODEL='Fixed2Prop', report.prop.recent=1.0,
					s.PREV.max.n=NA, s.INTERVENTION.prop=NA, s.INTERVENTION.start=2015, s.INTERVENTION.mul= NA, s.ARCHIVAL.n=50,
					epi.model='HPTN071', epi.dt=1/48, epi.import=NA, root.edge.fixed=0,
					v.N0tau=1, v.r=2.851904, v.T50=-2,
					wher.mu=log(0.00447743)-0.5^2/2, wher.sigma=0.5, bwerm.mu=log(0.002239075)-0.3^2/2, bwerm.sigma=0.3, er.gamma=4,
					dbg.GTRparam=0, dbg.rER=0, index.starttime.mode='fix1970', startseq.mode='one', seqtime.mode='AtDiag')								
pipeline.vary	<- data.table(	label=					c('O',	 'F',	'T',	'S',	'Q',	'I', 	'G',	'J',	'K',	'R',	'N',	'M',	'L',	'P',	'E',	'H'),
				epi.acute=		c('low', 'high','low',	'low',	'low',	'low',	'low',	'high',	'high',	'low',	'low',	'high',	'high',	'high',	'high',	'high'),
				epi.intervention=	c('fast','fast','fast',	'fast',	'slow',	'fast',	'slow',	'fast',	'slow',	'slow',	'none',	'none',	'fast',	'fast',	'slow',	'slow'),
				yr.end=			c(2018,	 2018,  2020,  	2020,   2020,    2020, 	2020,	2020,	2020,	2020,	2020,	2020,	2020,	2020,	2020,	2020),
				epi.import=		c(0.05,  0.05,  0.05,   0.05,   0.05,    0.05,	0.05,	0.05,	0.05,	0.05,	0.05,	0.05,	0.05,	0.2,	0.2,	0.05),
				s.PREV.max.n=		c(1280,  1280,  1600,  	1600,   1600,    3200, 	3200,	3200,	3200,	1600,	1600,	1600,	1600,	1600,	1600,	1600),
				s.INTERVENTION.prop=	c(0.375, 0.375, 0.5,   	0.85,   0.85,    0.5,	0.5,	0.5,	0.5,	0.5,	0.5,	0.5,	0.5,	0.5,	0.5,	0.5),
				seed=                   c(17,    17,    5,     	13,     13,      11, 	11,		11,		11,		5,		5,		5,		5,		7,		7,		5))
invisible(pipeline.vary[, {	
			set(pipeline.args, which( pipeline.args$stat=='epi.acute' ), 'v', as.character(epi.acute))
			set(pipeline.args, which( pipeline.args$stat=='epi.intervention' ), 'v', as.character(epi.intervention))																	
			set(pipeline.args, which( pipeline.args$stat=='yr.end' ), 'v', as.character(yr.end))
			set(pipeline.args, which( pipeline.args$stat=='epi.import' ), 'v', as.character(epi.import))
			set(pipeline.args, which( pipeline.args$stat=='s.PREV.max.n' ), 'v', as.character(s.PREV.max.n))											
			set(pipeline.args, which( pipeline.args$stat=='s.INTERVENTION.prop' ), 'v', as.character(s.INTERVENTION.prop))
			set(pipeline.args, which( pipeline.args$stat=='s.seed' ), 'v', as.character(seed))
			tmpdir			<- paste(outdir,'-Dataset',label,sep='')
			dir.create(tmpdir, showWarnings=FALSE)																														
			file			<- sim.regional(tmpdir, pipeline.args=pipeline.args)
			#system(file)				
			}])
}
