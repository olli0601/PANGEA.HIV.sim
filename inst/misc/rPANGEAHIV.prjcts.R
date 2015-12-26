


prog.hello<- function()	
{
	print('hello')
}
##--------------------------------------------------------------------------------------------------------
##	select between host sequences
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.DATA.SSAfg.select.between.host<- function()
{
	label.sep				<- '|'
	label.idx.country.id	<- 2
	label.idx.label			<- 3
	label.idx.ctime			<- 4		
	#read hand aligned sequences
	DATA		<<- '~/git/HPTN071sim/data_rootseq'
	infile.gag	<- 'PANGEA_SSAfg_gag_140806_n556_final.fasta'
	infile.pol	<- 'PANGEA_SSAfg_pol_140806_n556_final.fasta'
	infile.env	<- 'PANGEA_SSAfg_env_140806_n556_final.fasta'
	seq.gag		<- read.dna( paste(DATA,'/',infile.gag, sep=''), format='fasta' )
	seq.pol		<- read.dna( paste(DATA,'/',infile.pol, sep=''), format='fasta' )
	seq.env		<- read.dna( paste(DATA,'/',infile.env, sep=''), format='fasta' )
	#check seq names are the same
	stopifnot( all(rownames(seq.gag)==rownames(seq.pol)), all(rownames(seq.gag)==rownames(seq.env)) )	
	#ensure unique seq names
	tmp			<- paste( rownames(seq.gag), seq_len(nrow(seq.gag)), sep=label.sep )
	rownames(seq.gag)	<- rownames(seq.pol)	<- rownames(seq.env)	<- tmp	
	tmp			<- strsplit( rownames(seq.gag), label.sep, fixed=1 )
	label.info	<- data.table( SEQ_NAME=sapply(tmp, function(x) paste(x, collapse=label.sep)), COUNTRY_ID= sapply(tmp, '[[', label.idx.country.id), LABEL= sapply(tmp, '[[', label.idx.label), CALENDAR_TIME= sapply(tmp, '[[', label.idx.ctime))
	label.info	<- label.info[, {
							list(SEQ_NAME=SEQ_NAME, COUNTRY_ID=COUNTRY_ID, LABEL_UNIQUE=paste(LABEL,seq_along(SEQ_NAME),sep=label.sep), CALENDAR_TIME=CALENDAR_TIME)
						}, by=LABEL]
	setkey(label.info, SEQ_NAME)
	tmp			<- label.info[rownames(seq.gag), ][, paste('C',COUNTRY_ID,LABEL_UNIQUE,CALENDAR_TIME,sep=label.sep) ]
	rownames(seq.gag)	<- rownames(seq.pol)	<- rownames(seq.env)	<- tmp
	set(label.info, NULL, 'SEQ_NAME', label.info[, paste('C',COUNTRY_ID,LABEL_UNIQUE,CALENDAR_TIME,sep=label.sep) ])	
	#concatenate gag pol env
	seq.c		<- do.call('cbind', list(seq.gag, seq.pol, seq.env))
	#keep only first sequence with same label	
	tmp			<- label.info[, as.integer(sapply( strsplit( LABEL_UNIQUE, label.sep, fixed=1 ), '[[', 2 ))]
	label.info[, LABEL_NO:= tmp]
	label.info[, SELECT:= FALSE]	
	#search for further within host seqs: break up label by '.' or '_'
	setkey(label.info, LABEL)
	tmp			<- label.info[, which( grepl('\\.', LABEL) | grepl('_', LABEL) ) ]	
	#select manually not within host seqs
	include		<- c( "00BW3876_9", "00BW3886_8", "00BW3891_6", "00BW3970_2", "00BW5031_1","02ET_288","4403bmLwk4_fl7",                   
		"702010141_CH141.w12","702010293_CH293.w8a","702010432_CH432.w4","702010440_CH440.w4","703010085_CH085.w4a",
		"703010131_CH131_TF", "703010167_CH167.w8", "703010200_CH200_TFa", "703010228_CH228_TFa", "703010256_CH256.w96", 
		"703010269_CH269.w24", "704010042_CH042.mo6", "705010067_CH067_TF", "705010162_CH162.mo6", "705010185_CH185.mo6", 
		"705010198_CH198_TF", "705010534_CH534.w12", "706010164_CH164.mo6", "707010457_CH457.w8", "89SM_145", "90SE_364", 
		"93MW_965", "96BWMO1_5", "BD16_10", "BD22_11", "BD39_8", "BD9_11", "C.703010159.S.0dps.fl", "C.704010042.S.0dps.fl", 
		"C.705010162.e.wg2", "C.CAP210.w02.0dps.1_00_F4","C.CAP239.w02.0dps.1_02_F32", "C.CAP45.w02.0dps.1_05_T1", "CAP174_4w", "CAP206_8w_F1",	              
		"CAP210_5w", "CAP228_8w_F2", "CAP229_7w", "CAP239_5w_F1", "CAP244_8w_F1", "CAP248_9w", "CAP255_8w_F1", "CAP256_6w",                    
		"CAP257_7w_F1", "CAP30_5w_F4", "CAP45_5w_F1", "CAP61_8w_F3", "CAP63_5w_F4", "CAP65_6w", "CAP84_3w_F2", "CAP85_5w_F1",                  
		"CAP88_5w_F2", "CAP8_3w_F2", "C_ZA_1069MB", "C_ZA_1184MB", "C_ZA_1189MB", "C_ZA_J112MA", "TV001_patent", "TV002_patent",
		"ZM246F_flA1", "ZM247F_flA1", "ZM249M_flC1", "pZAC_R3714")
	include		<- c(include, label.info[ !grepl('\\.', LABEL) & !grepl('_', LABEL),   ][, LABEL])
	#these are excluded based on the '.' and '_' search:
	#c("4403bmLwk4_fl11","702010293_CH293.w8b","703010200_CH200_TFb","703010200_CH200_TFc", "703010228_CH228_TFb",
	#	"703010256_CH256_TF", "704010042_CH042_TF", "705010162_CH162_TF", "705010185_CH185_TF", "706010164_CH164_TF",
	#	"C.CAP210.w02.0dps.1_00_T11", "C.CAP210.w02.0dps.1_00_T2B", "C.CAP210.w02.0dps.1_00_T3", "C.CAP210.w02.0dps.1_00_T36",
	#	"C.CAP210.w02.0dps.1_00_T3C", "C.CAP210.w02.0dps.1_00_T4", "C.CAP210.w02.0dps.1_00_T43", "C.CAP210.w02.0dps.1_00_T5",
	#	"C.CAP210.w02.0dps.1_00_T6", "C.CAP210.w05.21dps.2_00", "C.CAP210.w12.70dps.2_05_T13", "C.CAP210.w12.70dps.2_05_T13C",
	#	"C.CAP210.w12.70dps.2_05_T2", "C.CAP210.w12.70dps.2_05_T39w", "C.CAP210.w12.70dps.2_05_T42", "C.CAP210.w12.70dps.2_05_T5",
	#	"C.CAP210.w12.70dps.2_05_T8", "C.CAP210.w26.168dps.3_10_T13B", "C.CAP210.w26.168dps.3_10_T20w", "C.CAP210.w26.168dps.3_10_T23B", 
	#	"C.CAP210.w26.168dps.3_10_T24B", "C.CAP210.w26.168dps.3_10_T24C", "C.CAP210.w26.168dps.3_10_T28", "C.CAP210.w26.168dps.3_10_T40",
	#	"C.CAP210.w26.168dps.3_10_T42", "C.CAP210.w26.168dps.3_10_T43B", "C.CAP210.w26.168dps.3_10_T47", "C.CAP210.w26.168dps.3_10_T49B",
	#	"C.CAP239.w02.0dps.1_02_T8", "C.CAP239.w05.21dps.2_00", "C.CAP239.w05.21dps.2_00_T11", "C.CAP239.w05.21dps.2_00_T17", 
	#	"C.CAP239.w05.21dps.2_00_T18", "C.CAP239.w05.21dps.2_00_T19", "C.CAP239.w05.21dps.2_00_T21", "C.CAP239.w05.21dps.2_00_T3",
	#	"C.CAP239.w05.21dps.2_00_T49", "C.CAP239.w05.21dps.2_00_T8", "C.CAP239.w11.63dps.2_05_T37", "C.CAP239.w11.63dps.2_05_T47",
	#	"C.CAP239.w117.805dps.4_21_T44", "C.CAP239.w117.805dps.4_21_T46", "C.CAP239.w117.805dps.4_21_T50", "C.CAP239.w22.140dps.3_09_F1", 
	#	"C.CAP239.w22.140dps.3_09_T10", "C.CAP239.w22.140dps.3_09_T17", "C.CAP239.w22.140dps.3_09_T20", "C.CAP239.w22.140dps.3_09_T36", 
	#	"C.CAP239.w22.140dps.3_09_T39", "C.CAP239.w22.140dps.3_09_W16", "C.CAP45.w02.0dps.1_05_F3", "C.CAP45.w02.0dps.1_05_T2", "C.CAP45.w05.21dps.2_00", 
	#	"C.CAP45.w05.21dps.2_00_T11", "C.CAP45.w05.21dps.2_00_T12", "C.CAP45.w05.21dps.2_00_T14", "C.CAP45.w05.21dps.2_00_T5", "C.CAP45.w05.21dps.2_00_T9", 
	#	"C.CAP45.w12.70dps.2_05_F1", "C.CAP45.w12.70dps.2_05_T11", "C.CAP45.w12.70dps.2_05_T13b", "C.CAP45.w12.70dps.2_05_T18",
	#	"C.CAP45.w65.455dps.4_17_T14", "C.CAP45.w65.455dps.4_17_T14B", "ZM246F_flA10", "ZM246F_flA2", "ZM246F_flA6", "ZM246F_flB1", 
	#	"ZM246F_flC12", "ZM246F_flC3" , "ZM246F_flC5", "ZM246F_flC7", "ZM246F_flD5", "ZM247F_flA12", "ZM247F_flA2", "ZM247F_flB8", 
	#	"ZM247F_flB9", "ZM247F_flE10", "ZM247F_flE11", "ZM247F_flE3", "ZM247F_flF10", "ZM247F_flF7", "ZM247F_flG11", "ZM247F_flH1",
	#	"ZM249M_flC5", "ZM249M_flE10", "ZM249M_flE8", "ZM249M_flF1", "TV001_patent", "TV001_patent", "TV002_patent", "chimeric_MJ4")
	exclude	<- c(	'C|ZA|C_ZA_1184MB|1|2000','C|ZA|C_ZA_1189MB|1|2000','C|ZA|C_ZA_J112MA|1|2000','C|ZA|C_ZA_1069MB|1|2000',
					'C|ZA|TV001_patent|1|1998','C|ZA|TV002_patent|1|1998','C|ZA|03ZASK212B1|1|2003',
					'C|ZA|C.CAP239.w02.0dps.1_02_F32|1|2005',"C|ZA|C.CAP210.w02.0dps.1_00_F4|1|2005","C|ZA|C.CAP45.w02.0dps.1_05_T1|1|2005",
					'C|MW|C.703010159.S.0dps.fl|1|2007','C|ZA|C.704010042.S.0dps.fl|1|2007','C|ZA|705010162_CH162.mo6|1|2007',
					'C|BW|96BW15C05|1|1996','C|BW|96BW15B03|1|1996',
					'C|BW|96BW01B22|1|1996','C|BW|96BW01B03|1|1996',"C|BW|96BW11B01|1|1996","C|BW|96BW1104|1|1996",
					'C|BW|96BW16B01|1|1996','C|BW|96BW1626|1|1996', "C|BW|96BW0502|1|1996", "C|BW|96BW06H51|1|1996","C|BW|96BW06|1|1996",
					"C|BW|96BW0407|1|1996","C|BW|96BW0402|1|1996","C|BW|96BW0410|1|1996","C|BW|96BW0408|1|1996")
	set(label.info, label.info[, which(  LABEL_NO==1L & LABEL%in%include )], 'SELECT', TRUE)
	set(label.info, label.info[, which(  SEQ_NAME%in%exclude )], 'SELECT', FALSE)
	#select
	tmp			<- subset(label.info, SELECT)[, SEQ_NAME]
	seq.c		<- seq.c[tmp, ]
	#save the whole lot
	tmp			<- rownames(seq.c)
	seq.gag		<- seq.gag[tmp, ]
	seq.pol		<- seq.pol[tmp, ]
	seq.env		<- seq.env[tmp, ]
	seq			<- seq.c			#need 'seq' because expected for 3SEQ
	outdir		<- '~/duke/2014_Gates/methods_comparison_rootseqsim/140811'
	outfile		<- 'PANGEA_SSAfgBwh_140811_n415_final.R'
	file		<- paste(outdir, '/', outfile, sep='')
	save(seq, seq.gag, seq.pol, seq.env, file=file)
	#check for recombinants
}
##--------------------------------------------------------------------------------------------------------
##	run 3SEQ
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.3SEQ.SSAfg.rm.recombinants<- function()
{
	require(XML)
	require(ape)
	require(r3SEQ)
	#DATA			<<- "/work/or105/Gates_2014"
	DATA			<<- '/Users/Oliver/duke/2014_Gates'
	indir			<- paste(DATA,'methods_comparison_rootseqsim/140811',sep='/')	
	infile			<- 'PANGEA_SSAfgBwh_140811_n415_final.R'
	outfile			<- 'PANGEA_SSAfgBwhRc-_140811_n390.R'
	#
	#	run 3SEQ
	r3seq.pipe.run.3seq(indir, infile, batch.n=5, hpc.walltime=1, hpc.q=NA, hpc.mem="500mb", hpc.nproc=1)
	#	parse 3SEQ output
	argv			<<-	r3seq.cmd.process.3SEQ.output(indir, infile, '', resume=1, verbose=1) 
	argv			<<- unlist(strsplit(argv,' '))
	df.recomb		<- r3seq.prog.process.3SEQ.output()	
	#	subset( df.recomb, adjp<1e-4 & min_rec_length>500)[, hist(log10(adjp), breaks=100)]
	df.recomb		<- subset( df.recomb, adjp<1e-7 & min_rec_length>500)
	cat(paste('\nfound potential recombinants, n=',nrow(df.recomb)))
	#
	file		<- paste(indir, '/', infile, sep='')
	load(file)
	tmp			<- setdiff(rownames(seq), df.recomb[, child])
	seq.gag		<- seq.gag[tmp,]
	seq.pol		<- seq.pol[tmp,]
	seq.env		<- seq.env[tmp,]
	seq			<- seq[tmp,]
	file		<- paste(indir, '/', outfile, sep='')
	save(seq, seq.gag, seq.pol, seq.env, file=file)
}
##--------------------------------------------------------------------------------------------------------
##	run BEAST XML file
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.BEAST.SSAfg.runXML<- function()
{
	#DATA		<<- "/work/or105/Gates_2014"
	DATA		<<- '/Users/Oliver/duke/2014_Gates'	
	indir		<- paste(DATA,'methods_comparison_rootseqsim/140907',sep='/')
	#search for XML files in indir
	infiles		<- list.files(indir, pattern=paste(".xml$",sep=''))
	insignat	<- ''	
	hpc.ncpu	<- 8
	
	for(infile in infiles)
	{
		infile		<- substr(infile, 1, nchar(infile)-4) 		
		cmd			<- hivc.cmd.beast.runxml(indir, infile, insignat, prog.beast=PR.BEAST, prog.beast.opt=" -beagle -working", hpc.tmpdir.prefix="beast", hpc.ncpu=hpc.ncpu)
		tmp			<- paste(infile,'.timetrees',sep='')	
		cmd			<- paste(cmd, hivc.cmd.beast.read.nexus(indir, tmp, indir, tree.id=NA, method.node.stat='any.node'), sep='\n')
		cmd			<- paste(cmd, hivc.cmd.beast.run.treeannotator(indir, infile, insignat, prog.beastmcc=PR.BEASTMCC, beastmcc.burnin=500, beastmcc.heights="median"), sep='\n')
		cat(cmd)	
		cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.q="pqeph", hpc.nproc=hpc.ncpu, hpc.walltime=91, hpc.mem="3700mb")		
		outdir		<- indir
		outfile		<- paste("b2m.",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='')					
		hivc.cmd.hpccaller(outdir, outfile, cmd)		
	}
}
##--------------------------------------------------------------------------------------------------------
##	run BEAST XML file
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.BEAST.SSApg.runXML<- function()
{
	#DATA		<<- "/work/or105/Gates_2014"
	DATA		<<- '/Users/Oliver/duke/2014_Gates'	
	indir		<- paste(DATA,'methods_comparison_rootseqsim/140902',sep='/')
	#search for XML files in indir
	infiles		<- list.files(indir, pattern=paste(".xml$",sep=''))
	insignat	<- ''	
	hpc.ncpu	<- 8
	
	for(infile in infiles)
	{
		infile		<- substr(infile, 1, nchar(infile)-4) 		
		cmd			<- hivc.cmd.beast.runxml(indir, infile, insignat, prog.beast=PR.BEAST, prog.beast.opt=" -beagle -working", hpc.tmpdir.prefix="beast", hpc.ncpu=hpc.ncpu)
		tmp			<- paste(infile,'.timetrees',sep='')	
		cmd			<- paste(cmd, hivc.cmd.beast.read.nexus(indir, tmp, indir, tree.id=NA, method.node.stat='any.node'), sep='\n')
		cmd			<- paste(cmd, hivc.cmd.beast.run.treeannotator(indir, infile, insignat, prog.beastmcc=PR.BEASTMCC, beastmcc.burnin=500, beastmcc.heights="median"), sep='\n')
		cat(cmd)	
		cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.q="pqeph", hpc.nproc=hpc.ncpu, hpc.walltime=791, hpc.mem="3700mb")		
		outdir		<- indir
		outfile		<- paste("bpg.",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='')					
		hivc.cmd.hpccaller(outdir, outfile, cmd)		
	}
}
##--------------------------------------------------------------------------------------------------------
##	run BEAST XML file
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.BEAST.SSApg.runXML<- function()
{
	#DATA		<<- "/work/or105/Gates_2014"
	DATA		<<- '/Users/Oliver/duke/2014_Gates'	
	indir		<- paste(DATA,'methods_comparison_rootseqsim/140830',sep='/')
	#search for XML files in indir
	infiles		<- list.files(indir, pattern=paste(".xml$",sep=''))
	insignat	<- ''	
	hpc.ncpu	<- 8
	
	for(infile in infiles)
	{
		infile		<- substr(infile, 1, nchar(infile)-4) 		
		cmd			<- hivc.cmd.beast.runxml(indir, infile, insignat, prog.beast=PR.BEAST, prog.beast.opt=" -beagle -working", hpc.tmpdir.prefix="beast", hpc.ncpu=hpc.ncpu)
		tmp			<- paste(infile,'.timetrees',sep='')	
		cmd			<- paste(cmd, hivc.cmd.beast.read.nexus(indir, tmp, indir, tree.id=NA, method.node.stat='any.node'), sep='\n')
		cmd			<- paste(cmd, hivc.cmd.beast.run.treeannotator(indir, infile, insignat, prog.beastmcc=PR.BEASTMCC, beastmcc.burnin=500, beastmcc.heights="median"), sep='\n')
		cat(cmd)	
		cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.q="pqeph", hpc.nproc=hpc.ncpu, hpc.walltime=91, hpc.mem="3700mb")		
		outdir		<- indir
		outfile		<- paste("b2m.",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='')					
		hivc.cmd.hpccaller(outdir, outfile, cmd)		
	}
}
##--------------------------------------------------------------------------------------------------------
##	create BEAST XML file for each gene
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.BEAST.SSApg.createXML<- function()
{
	require(hivclust)
	require(XML)
	require(ape)
	require(r3SEQ)
	#DATA			<<- "/work/or105/Gates_2014"
	DATA			<<- '/Users/Oliver/duke/2014_Gates'
	s.seed			<- 42
	if(1)	
	{		
		#
		#	to define sequences for each BEAST run
		#	compute NJ tree and define clusters
		#	
		#	v09: estimate frequencies
		#	v10: empirical frequencies -- chain gets 'stuck' in better lkl for particular freq for A/T
		#		 center evol rate on previous estimates (Wawer, Delatorre) Normal( 1955, with 95%CI 1935-1975-> sigma=12 )
		#		 evol rate prior U(0.0001-0.006)
		infile.beast.gag<- '/Users/Oliver/git/HPTN071sim/data_rootseq/BEAST_template_v10gag.xml'
		infile.beast.pol<- '/Users/Oliver/git/HPTN071sim/data_rootseq/BEAST_template_v10pol.xml'
		infile.beast.env<- '/Users/Oliver/git/HPTN071sim/data_rootseq/BEAST_template_v10env.xml'
		indir			<- paste(DATA,'methods_comparison_rootseqsim/140907',sep='/')	
		infile			<- 'PANGEA_SSAfgBwhRc-_140811_n390.R'
		file			<- paste(indir, '/', infile, sep='')
		load(file)		
		#	remove sequences without calendar time 		
		label.sep				<- '|'
		label.idx.ctime			<- 5		
		tmp						<- sapply( strsplit( rownames(seq.gag), label.sep, fixed=1 ), '[[', label.idx.ctime )
		tmp						<- rownames(seq.gag)[ which(is.na(as.numeric(tmp))) ]
		cat(paste('\nExclude sequences with no calendar date, ', paste(tmp, collapse=' ')))
		tmp						<- setdiff(rownames(seq.gag), tmp)		
		seq.gag					<- seq.gag[tmp,]
		#	exclude last 2 nucleotides in gag to avoid incomple AA
		seq.gag					<- seq.gag[,1:1440]
		seq.pol					<- seq.pol[tmp,]
		seq.env					<- seq.env[tmp,]
		seq						<- seq[tmp,]
		#	get NJ tree and plot
		tmp				<- dist.dna( seq )
		seq.ph			<- nj(tmp)				
		file			<- paste( indir, '/', substr(infile,1,nchar(infile)-2), '_njtree.pdf', sep='' )	
		pdf(file=file, w=10, h=80)
		plot(seq.ph, show.tip=TRUE)
		dev.off()			
		#
		#	get 3 sequence pools of equal size
		#
		set.seed(s.seed)
		pool.n			<- 3
		tmp				<- hivc.clu.brdist.stats(seq.ph, eval.dist.btw="leaf", stat.fun=hivc.clu.min.transmission.cascade)
		thresh.brl		<- 0.055
		clustering		<- hivc.clu.clusterbythresh(seq.ph, thresh.brl=thresh.brl, dist.brl=tmp, retval="all")
		#	allocate clustering tips into 3 distinct clusters
		seq.clumem		<- data.table( PH_NODE_ID=seq_len(Ntip(seq.ph)), CLU_ID=clustering$clu.mem[ seq_len(Ntip(seq.ph)) ] )
		setkey(seq.clumem, CLU_ID)		
		tmp				<- which(!is.na(seq.clumem[, CLU_ID]))
		tmp				<- seq.clumem[tmp,][, list(CLU_N=-length(PH_NODE_ID)), by='CLU_ID']
		setkey(tmp, CLU_N)
		set(tmp, NULL, 'POOL_ID', tmp[, cumsum(-CLU_N)]) 		
		set(tmp, NULL, 'POOL_ID', tmp[, ceiling( POOL_ID / max(POOL_ID) * pool.n ) ] )
		seq.clumem		<- merge(seq.clumem, subset(tmp, select=c(CLU_ID, POOL_ID)), by='CLU_ID', all.x=TRUE)
		#	allocate non-clustering tips into 3 distinct clusters
		tmp				<- subset(seq.clumem,!is.na(POOL_ID))[, list(NOCLU_N= ceiling( nrow(seq.clumem) / pool.n ) - length(PH_NODE_ID)), by='POOL_ID']		
		set(tmp, 1L, 'NOCLU_N', tmp[1,NOCLU_N] - ( tmp[, sum(NOCLU_N)] - ( Ntip(seq.ph) - nrow(subset(seq.clumem,!is.na(POOL_ID))) )) )		
		set(seq.clumem, seq.clumem[, which(is.na(POOL_ID))], 'POOL_ID',  rep(tmp[,POOL_ID], tmp[,NOCLU_N]) )
		seq.clumem[, table(POOL_ID)]	
		#
		#	for each sequence pool, set up BEAST run
		#
		verbose				<- 1
		bxml.template.gag	<- xmlTreeParse(infile.beast.gag, useInternalNodes=TRUE, addFinalizer = TRUE)
		bxml.template.pol	<- xmlTreeParse(infile.beast.pol, useInternalNodes=TRUE, addFinalizer = TRUE)
		bxml.template.env	<- xmlTreeParse(infile.beast.env, useInternalNodes=TRUE, addFinalizer = TRUE)
		for(pool.id in seq_len(pool.n))
		{			
			pool.seqnames	<- seq.ph$tip.label[ subset(seq.clumem, POOL_ID==pool.id)[, PH_NODE_ID] ]
			#
			#
			#	
			cat(paste('\ncreate GAG BEAST XML file for seqs=',paste(pool.seqnames, collapse=' ')))
			pool.infile		<- paste(  substr(infile,1,nchar(infile)-2),'_geneGAG_pool',pool.id, sep='' )
			#	write XML file with new sequences
			bxml			<- newXMLDoc(addFinalizer=T)
			bxml.beast		<- newXMLNode("beast", doc=bxml, addFinalizer=T)
			tmp				<- newXMLCommentNode(text=paste("Generated by HIVCLUST from template",infile.beast.gag), parent=bxml.beast, doc=bxml, addFinalizer=T)
			#	add new set of GAG sequences into GAG alignment
			tmp				<- seq.gag[pool.seqnames,]
			bxml			<- hivc.beast.add.seq(bxml, tmp, df=NULL, beast.label.datepos= 5, beast.label.sep= '|', beast.date.direction= "forwards", beast.date.units= "years", beast.alignment.id="GAG.alignment", beast.alignment.dataType= "nucleotide", verbose=1)
			#	copy from template	
			bt.beast		<- getNodeSet(bxml.template.gag, "//beast")[[1]]
			dummy			<- sapply(seq.int( 1, xmlSize(bt.beast) ), function(i)
					{
						if( class(bt.beast[[i]])[1]=="XMLInternalCommentNode" )
							dummy<- newXMLCommentNode(text=xmlValue(bt.beast[[i]]), parent=bxml.beast, doc=bxml, addFinalizer=T)
						else
							dummy<- addChildren( bxml.beast, xmlClone( bt.beast[[i]], addFinalizer=T, doc=bxml ) )
					})
			#	change gmrf dimensions	
			tmp			<- getNodeSet(bxml, "//*[@id='skyride.logPopSize']")
			if(length(tmp)!=1)	stop("unexpected number of *[@id='skyride.logPopSize'")
			tmp			<- tmp[[1]]
			xmlAttrs(tmp)["dimension"]	<-	length(pool.seqnames)-1  
			tmp			<- getNodeSet(bxml, "//*[@id='skyride.groupSize']")
			if(length(tmp)!=1)	stop("unexpected number of *[@id='skyride.groupSize'")
			tmp			<- tmp[[1]]
			xmlAttrs(tmp)["dimension"]	<-	length(pool.seqnames)-1			
			#	change outfile name 
			bxml.onodes	<- getNodeSet(bxml, "//*[@fileName]")
			tmp			<- sapply(bxml.onodes, function(x) xmlGetAttr(x,"fileName"))
			tmp			<- gsub("(time).","time",tmp,fixed=1)
			tmp			<- gsub("(subst).","subst",tmp,fixed=1)	
			tmp			<- sapply(strsplit(tmp,'.',fixed=1), function(x)	paste(pool.infile, '.', tail(x,1), sep=''))
			dummy		<- sapply(seq_along(bxml.onodes), function(i){		xmlAttrs(bxml.onodes[[i]])["fileName"]<- tmp[i]		})
			#	write to file
			file		<- paste(indir,'/',pool.infile,".xml", sep='')
			if(verbose)	cat(paste("\nwrite xml file to",file))
			saveXML(bxml, file=file)
			#
			#	POL
			#
			cat(paste('\ncreate POL BEAST XML file for seqs=',paste(pool.seqnames, collapse=' ')))
			pool.infile		<- paste(  substr(infile,1,nchar(infile)-2),'_genePOL_pool',pool.id, sep='' )
			#	write XML file with new sequences
			bxml			<- newXMLDoc(addFinalizer=T)
			bxml.beast		<- newXMLNode("beast", doc=bxml, addFinalizer=T)
			tmp				<- newXMLCommentNode(text=paste("Generated by HIVCLUST from template",infile.beast.pol), parent=bxml.beast, doc=bxml, addFinalizer=T)
			#	add new set of GAG sequences into GAG alignment
			tmp				<- seq.pol[pool.seqnames,]
			bxml			<- hivc.beast.add.seq(bxml, tmp, df=NULL, beast.label.datepos= 5, beast.label.sep= '|', beast.date.direction= "forwards", beast.date.units= "years", beast.alignment.id="POL.alignment", beast.alignment.dataType= "nucleotide", verbose=1)
			#	copy from template	
			bt.beast		<- getNodeSet(bxml.template.pol, "//beast")[[1]]
			dummy			<- sapply(seq.int( 1, xmlSize(bt.beast) ), function(i)
					{
						if( class(bt.beast[[i]])[1]=="XMLInternalCommentNode" )
							dummy<- newXMLCommentNode(text=xmlValue(bt.beast[[i]]), parent=bxml.beast, doc=bxml, addFinalizer=T)
						else
							dummy<- addChildren( bxml.beast, xmlClone( bt.beast[[i]], addFinalizer=T, doc=bxml ) )
					})
			#	change gmrf dimensions	
			tmp			<- getNodeSet(bxml, "//*[@id='skyride.logPopSize']")
			if(length(tmp)!=1)	stop("unexpected number of *[@id='skyride.logPopSize'")
			tmp			<- tmp[[1]]
			xmlAttrs(tmp)["dimension"]	<-	length(pool.seqnames)-1  
			tmp			<- getNodeSet(bxml, "//*[@id='skyride.groupSize']")
			if(length(tmp)!=1)	stop("unexpected number of *[@id='skyride.groupSize'")
			tmp			<- tmp[[1]]
			xmlAttrs(tmp)["dimension"]	<-	length(pool.seqnames)-1			
			#	change outfile name 
			bxml.onodes	<- getNodeSet(bxml, "//*[@fileName]")
			tmp			<- sapply(bxml.onodes, function(x) xmlGetAttr(x,"fileName"))
			tmp			<- gsub("(time).","time",tmp,fixed=1)
			tmp			<- gsub("(subst).","subst",tmp,fixed=1)	
			tmp			<- sapply(strsplit(tmp,'.',fixed=1), function(x)	paste(pool.infile, '.', tail(x,1), sep=''))
			dummy		<- sapply(seq_along(bxml.onodes), function(i){		xmlAttrs(bxml.onodes[[i]])["fileName"]<- tmp[i]		})
			#	write to file
			file		<- paste(indir,'/',pool.infile,".xml", sep='')
			if(verbose)	cat(paste("\nwrite xml file to",file))
			saveXML(bxml, file=file)
			#
			#	ENV
			#
			cat(paste('\ncreate ENV BEAST XML file for seqs=',paste(pool.seqnames, collapse=' ')))
			pool.infile		<- paste(  substr(infile,1,nchar(infile)-2),'_geneENV_pool',pool.id, sep='' )
			#	write XML file with new sequences
			bxml			<- newXMLDoc(addFinalizer=T)
			bxml.beast		<- newXMLNode("beast", doc=bxml, addFinalizer=T)
			tmp				<- newXMLCommentNode(text=paste("Generated by HIVCLUST from template",infile.beast.env), parent=bxml.beast, doc=bxml, addFinalizer=T)
			#	add new set of GAG sequences into GAG alignment
			tmp				<- seq.env[pool.seqnames,]
			bxml			<- hivc.beast.add.seq(bxml, tmp, df=NULL, beast.label.datepos= 5, beast.label.sep= '|', beast.date.direction= "forwards", beast.date.units= "years", beast.alignment.id="ENV.alignment", beast.alignment.dataType= "nucleotide", verbose=1)
			#	copy from template	
			bt.beast		<- getNodeSet(bxml.template.env, "//beast")[[1]]
			dummy			<- sapply(seq.int( 1, xmlSize(bt.beast) ), function(i)
					{
						if( class(bt.beast[[i]])[1]=="XMLInternalCommentNode" )
							dummy<- newXMLCommentNode(text=xmlValue(bt.beast[[i]]), parent=bxml.beast, doc=bxml, addFinalizer=T)
						else
							dummy<- addChildren( bxml.beast, xmlClone( bt.beast[[i]], addFinalizer=T, doc=bxml ) )
					})
			#	change gmrf dimensions	
			tmp			<- getNodeSet(bxml, "//*[@id='skyride.logPopSize']")
			if(length(tmp)!=1)	stop("unexpected number of *[@id='skyride.logPopSize'")
			tmp			<- tmp[[1]]
			xmlAttrs(tmp)["dimension"]	<-	length(pool.seqnames)-1  
			tmp			<- getNodeSet(bxml, "//*[@id='skyride.groupSize']")
			if(length(tmp)!=1)	stop("unexpected number of *[@id='skyride.groupSize'")
			tmp			<- tmp[[1]]
			xmlAttrs(tmp)["dimension"]	<-	length(pool.seqnames)-1			
			#	change outfile name 
			bxml.onodes	<- getNodeSet(bxml, "//*[@fileName]")
			tmp			<- sapply(bxml.onodes, function(x) xmlGetAttr(x,"fileName"))
			tmp			<- gsub("(time).","time",tmp,fixed=1)
			tmp			<- gsub("(subst).","subst",tmp,fixed=1)	
			tmp			<- sapply(strsplit(tmp,'.',fixed=1), function(x)	paste(pool.infile, '.', tail(x,1), sep=''))
			dummy		<- sapply(seq_along(bxml.onodes), function(i){		xmlAttrs(bxml.onodes[[i]])["fileName"]<- tmp[i]		})
			#	write to file
			file		<- paste(indir,'/',pool.infile,".xml", sep='')
			if(verbose)	cat(paste("\nwrite xml file to",file))
			saveXML(bxml, file=file)			
		}		
	}
	#
}
##--------------------------------------------------------------------------------------------------------
##	create BEAST XML file
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.BEAST.SSAfg.createXML<- function()
{
	require(hivclust)
	require(XML)
	require(ape)
	require(r3SEQ)
	#DATA			<<- "/work/or105/Gates_2014"
	DATA			<<- '/Users/Oliver/duke/2014_Gates'
	s.seed			<- 42
	if(1)	
	{		
		#
		#	to define sequences for each BEAST run
		#	compute NJ tree and define clusters
		#
		infile.beast	<- '/Users/Oliver/git/HPTN071sim/data_rootseq/BEAST_template_v09.xml'
		indir			<- paste(DATA,'methods_comparison_rootseqsim/140830',sep='/')	
		infile			<- 'PANGEA_SSAfgBwhRc-_140811_n390.R'
		file			<- paste(indir, '/', infile, sep='')
		load(file)		
		#	remove sequences without calendar time 		
		label.sep				<- '|'
		label.idx.ctime			<- 5		
		tmp						<- sapply( strsplit( rownames(seq.gag), label.sep, fixed=1 ), '[[', label.idx.ctime )
		tmp						<- rownames(seq.gag)[ which(is.na(as.numeric(tmp))) ]
		cat(paste('\nExclude sequences with no calendar date, ', paste(tmp, collapse=' ')))
		tmp						<- setdiff(rownames(seq.gag), tmp)		
		seq.gag					<- seq.gag[tmp,]
		#	exclude last 2 nucleotides in gag to avoid incomple AA
		seq.gag					<- seq.gag[,1:1440]
		seq.pol					<- seq.pol[tmp,]
		seq.env					<- seq.env[tmp,]
		seq						<- seq[tmp,]
		#	get NJ tree and plot
		tmp				<- dist.dna( seq )
		seq.ph			<- nj(tmp)				
		file			<- paste( indir, '/', substr(infile,1,nchar(infile)-2), '_njtree.pdf', sep='' )	
		pdf(file=file, w=10, h=80)
		plot(seq.ph, show.tip=TRUE)
		dev.off()			
		#
		#	get 3 sequence pools of equal size
		#
		set.seed(s.seed)
		pool.n			<- 3
		tmp				<- hivc.clu.brdist.stats(seq.ph, eval.dist.btw="leaf", stat.fun=hivc.clu.min.transmission.cascade)
		thresh.brl		<- 0.055
		clustering		<- hivc.clu.clusterbythresh(seq.ph, thresh.brl=thresh.brl, dist.brl=tmp, retval="all")
		#	allocate clustering tips into 3 distinct clusters
		seq.clumem		<- data.table( PH_NODE_ID=seq_len(Ntip(seq.ph)), CLU_ID=clustering$clu.mem[ seq_len(Ntip(seq.ph)) ] )
		setkey(seq.clumem, CLU_ID)		
		tmp				<- which(!is.na(seq.clumem[, CLU_ID]))
		tmp				<- seq.clumem[tmp,][, list(CLU_N=-length(PH_NODE_ID)), by='CLU_ID']
		setkey(tmp, CLU_N)
		set(tmp, NULL, 'POOL_ID', tmp[, cumsum(-CLU_N)]) 		
		set(tmp, NULL, 'POOL_ID', tmp[, ceiling( POOL_ID / max(POOL_ID) * pool.n ) ] )
		seq.clumem		<- merge(seq.clumem, subset(tmp, select=c(CLU_ID, POOL_ID)), by='CLU_ID', all.x=TRUE)
		#	allocate non-clustering tips into 3 distinct clusters
		tmp				<- subset(seq.clumem,!is.na(POOL_ID))[, list(NOCLU_N= ceiling( nrow(seq.clumem) / pool.n ) - length(PH_NODE_ID)), by='POOL_ID']		
		set(tmp, 1L, 'NOCLU_N', tmp[1,NOCLU_N] - ( tmp[, sum(NOCLU_N)] - ( Ntip(seq.ph) - nrow(subset(seq.clumem,!is.na(POOL_ID))) )) )		
		set(seq.clumem, seq.clumem[, which(is.na(POOL_ID))], 'POOL_ID',  rep(tmp[,POOL_ID], tmp[,NOCLU_N]) )
		seq.clumem[, table(POOL_ID)]	
		#
		#	for each sequence pool, set up BEAST run
		#
		verbose			<- 1
		bxml.template	<- xmlTreeParse(infile.beast, useInternalNodes=TRUE, addFinalizer = TRUE)
		for(pool.id in seq_len(pool.n))
		{
			pool.infile		<- paste(  substr(infile,1,nchar(infile)-2),'_pool',pool.id, sep='' )
			pool.seqnames	<- seq.ph$tip.label[ subset(seq.clumem, POOL_ID==pool.id)[, PH_NODE_ID] ]
			#
			if(0)
			{
				cat(paste('\ncreate FASTA file for seqs=',paste(pool.seqnames, collapse=' ')))
				file			<- paste(indir,'/',pool.infile,"_GAG.fa", sep='')
				write.dna(seq.gag[pool.seqnames,], format='fasta', file=file)
				file			<- paste(indir,'/',pool.infile,"_POL.fa", sep='')
				write.dna(seq.pol[pool.seqnames,], format='fasta', file=file)
				file			<- paste(indir,'/',pool.infile,"_ENV.fa", sep='')
				write.dna(seq.env[pool.seqnames,], format='fasta', file=file)
				tmp				<- do.call('cbind',list( seq.gag[pool.seqnames,],seq.pol[pool.seqnames,],seq.env[pool.seqnames,] ))
				file			<- paste(indir,'/',pool.infile,".fa", sep='')
				write.dna(tmp, format='fasta', file=file)				
			}
			#
			cat(paste('\ncreate BEAST XML file for seqs=',paste(pool.seqnames, collapse=' ')))
			#	write XML file with new sequences
			bxml			<- newXMLDoc(addFinalizer=T)
			bxml.beast		<- newXMLNode("beast", doc=bxml, addFinalizer=T)
			newXMLCommentNode(text=paste("Generated by HIVCLUST from template",infile.beast), parent=bxml.beast, doc=bxml, addFinalizer=T)
			#	add new set of ENV sequences into alignment ID 1
			tmp				<- seq.env[pool.seqnames,]
			bxml			<- hivc.beast.add.seq(bxml, tmp, df=NULL, beast.label.datepos= 5, beast.label.sep= '|', beast.date.direction= "forwards", beast.date.units= "years", beast.alignment.id="alignment1", beast.alignment.dataType= "nucleotide", verbose=1)
			#	add new set of GAG sequences into alignment ID 2
			tmp				<- seq.gag[pool.seqnames,]
			bxml			<- hivc.beast.add.seq(bxml, tmp, df=NULL, beast.label.datepos= 5, beast.label.sep= '|', beast.date.direction= "forwards", beast.date.units= "years", beast.alignment.id="alignment2", beast.alignment.dataType= "nucleotide", verbose=1)
			#	add new set of POL sequences into alignment ID 3
			tmp				<- seq.pol[pool.seqnames,]
			bxml			<- hivc.beast.add.seq(bxml, tmp, df=NULL, beast.label.datepos= 5, beast.label.sep= '|', beast.date.direction= "forwards", beast.date.units= "years", beast.alignment.id="alignment3", beast.alignment.dataType= "nucleotide", verbose=1)
			#	copy from template	
			bt.beast		<- getNodeSet(bxml.template, "//beast")[[1]]
			dummy			<- sapply(seq.int( 1, xmlSize(bt.beast) ), function(i)
					{
						if( class(bt.beast[[i]])[1]=="XMLInternalCommentNode" )
							dummy<- newXMLCommentNode(text=xmlValue(bt.beast[[i]]), parent=bxml.beast, doc=bxml, addFinalizer=T)
						else
							dummy<- addChildren( bxml.beast, xmlClone( bt.beast[[i]], addFinalizer=T, doc=bxml ) )
					})
			#	change gmrf dimensions	
			tmp			<- getNodeSet(bxml, "//*[@id='skyride.logPopSize']")
			if(length(tmp)!=1)	stop("unexpected number of *[@id='skyride.logPopSize'")
			tmp			<- tmp[[1]]
			xmlAttrs(tmp)["dimension"]	<-	length(pool.seqnames)-1  
			tmp			<- getNodeSet(bxml, "//*[@id='skyride.groupSize']")
			if(length(tmp)!=1)	stop("unexpected number of *[@id='skyride.groupSize'")
			tmp			<- tmp[[1]]
			xmlAttrs(tmp)["dimension"]	<-	length(pool.seqnames)-1			
			#	change outfile name 
			bxml.onodes	<- getNodeSet(bxml, "//*[@fileName]")
			tmp			<- sapply(bxml.onodes, function(x) xmlGetAttr(x,"fileName"))
			tmp			<- gsub("(time).","time",tmp,fixed=1)
			tmp			<- gsub("(subst).","subst",tmp,fixed=1)	
			tmp			<- sapply(strsplit(tmp,'.',fixed=1), function(x)	paste(pool.infile, '.', tail(x,1), sep=''))
			dummy		<- sapply(seq_along(bxml.onodes), function(i){		xmlAttrs(bxml.onodes[[i]])["fileName"]<- tmp[i]		})
			#	write to file
			file		<- paste(indir,'/',pool.infile,".xml", sep='')
			if(verbose)	cat(paste("\nwrite xml file to",file))
			saveXML(bxml, file=file)
		}		
	}
	#
}
##--------------------------------------------------------------------------------------------------------
##	process BEAST log file and extract GTR parameters
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.BEAST.SSAfg.getGTR<- function()
{
	tree.id.burnin		<- 2e7
	tree.id.labelsep	<- '|'
	dir.name			<- '/Users/Oliver/duke/2014_Gates'  	
	indir				<- paste(dir.name,'methods_comparison_rootseqsim/140813',sep='/')
	outdir				<- indir
	#	search for BEAST output
	infiles				<- list.files(indir)
	infiles				<- infiles[ sapply(infiles, function(x) grepl('pool[0-9].log$',x) ) ]
	#	collect log variables
	log.df	<- lapply(seq_along(infiles), function(i)
			{
				infile	<- infiles[i]
				cat(paste('\nprocess file', infile))
				file	<- paste(indir, '/', infile, sep='')
				df		<- as.data.table(read.delim(file, comment.char='#'))
				cat(paste('\nignore logs for\n',paste(colnames(df)[ !grepl('state|POL|GAG|ENV|ucld|meanRate',colnames(df)) ], collapse=', ') ))	
				df		<- subset(df, select=which(grepl('state|POL|GAG|ENV|ucld|meanRate',colnames(df))))
				log.df	<- c( paste('frequencies',1:4,sep=''), 'mu', 'alpha', '\\.at', '\\.ac', '\\.cg', '\\.ag', '\\.gt', 'treeLikelihood' )
				log.df	<- lapply( log.df, function(x)
						{
							tmp		<- melt( subset(df, select=which(grepl(paste('state|',x,sep=''),colnames(df)))), id='state', value.name=x)
							tmp[, GENE:= tmp[,substr(variable, 1, 3)]]
							tmp[, CODON_POS:=tmp[, regmatches(variable, regexpr('CP[1-3]',variable))]]
							subset(tmp, select=which(colnames(tmp)!='variable'))				
						})
				tmp		<- log.df[[1]]
				for(j in seq_along(log.df)[-1])
					tmp	<- merge(tmp, log.df[[j]], by=c('state','GENE','CODON_POS'))
				log.df	<- tmp	
				log.df	<- merge(log.df, subset(df, select=which(grepl('state|ucld|meanRate',colnames(df)))), by='state')
				log.df[, FILE:= regmatches(infile, regexpr('pool[0-9]+',infile))]
				log.df
			})
	log.df	<- do.call('rbind',log.df)	
	setnames(log.df, colnames(log.df), gsub('\\.','',colnames(log.df),fixed=TRUE))
	setnames(log.df, paste('frequencies',1:4,sep=''), c('a','c','g','t') )
	log.df	<- subset(log.df, state>tree.id.burnin)
	#	check mean of relative rates
	#tmp	<- log.df[, list(mmu=mean(mu), n=length(mu)), by=c('FILE','GENE', 'state')]
	
	file	<- paste( substr(infiles[1],1,nchar(infiles[1])-9),'log.R',sep='' )
	file	<- paste( outdir, '/', file, sep='')
	cat(paste('\nsave to file', file))
	save(file=file, log.df)
	# 
	tmp		<- copy(log.df)
	tmp		<- melt(tmp, id=c('state','GENE','CODON_POS','FILE'))
	ggplot( tmp, aes(x=value, fill=FILE)) + geom_histogram() + facet_grid(GENE+CODON_POS~variable, scales='free')
	file	<- paste( substr(file,1,nchar(file)-1),'pdf',sep='' )
	ggsave(file, h=15, w=20)
}
##--------------------------------------------------------------------------------------------------------
##	process BEAST log file and extract GTR parameters
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.BEAST.SSApg.getGTR<- function()
{
	require(XML)
	tree.id.burnin		<- 2e7
	tree.id.labelsep	<- '|'
	dir.name			<- '/Users/Oliver/duke/2014_Gates'  	
	indir				<- paste(dir.name,'methods_comparison_rootseqsim/140907',sep='/')
	outdir				<- indir
	#	compute frequencies from BEAST xml
	infiles				<- list.files(indir)
	infiles				<- infiles[ sapply(infiles, function(x) grepl('pool[0-9].xml$',x) ) ]
	freq	<- lapply(seq_along(infiles), function(i)
			{
				infile	<- infiles[i]
				cat(paste('\nprocess file', infile))
				file	<- paste(indir, '/', infile, sep='')
				bxml		<- xmlTreeParse(file, useInternalNodes=TRUE, addFinalizer = TRUE)
				bseq		<- hivc.beast.get.sequences(bxml, verbose=1)					
				tmp			<- tolower(do.call('rbind',strsplit(bseq[, SEQ],'')))
				bseq.CP1	<- as.DNAbin( tmp[, seq.int(1, ncol(tmp), by=3)] )
				bseq.CP2	<- as.DNAbin( tmp[, seq.int(2, ncol(tmp), by=3)] )
				bseq.CP3	<- as.DNAbin( tmp[, seq.int(3, ncol(tmp), by=3)] )
				bseq[1, substr(ALIGNMENT_ID,1,3)]				
				tmp			<- list(	data.table(FREQ=c('a','c','g','t'), VALUE=base.freq( bseq.CP1 ), CODON_POS='CP1', GENE=bseq[1, substr(ALIGNMENT_ID,1,3)]),
									 	data.table(FREQ=c('a','c','g','t'), VALUE=base.freq( bseq.CP2 ), CODON_POS='CP2', GENE=bseq[1, substr(ALIGNMENT_ID,1,3)]),
										data.table(FREQ=c('a','c','g','t'), VALUE=base.freq( bseq.CP3 ), CODON_POS='CP3', GENE=bseq[1, substr(ALIGNMENT_ID,1,3)])	)
				tmp			<- do.call('rbind', tmp)
				tmp[, FILE:=regmatches(infile, regexpr('pool[0-9]+',infile))]
				tmp
			})
	freq				<- do.call('rbind', freq)
	freq				<- dcast.data.table(freq, CODON_POS+GENE+FILE~FREQ, value.var='VALUE')
	#	search for BEAST output
	infiles				<- list.files(indir)
	infiles				<- infiles[ sapply(infiles, function(x) grepl('pool[0-9].log$',x) ) ]
	#	collect log variables
	log.df	<- lapply(seq_along(infiles), function(i)
			{
				infile	<- infiles[i]
				cat(paste('\nprocess file', infile))
				file	<- paste(indir, '/', infile, sep='')
				df		<- as.data.table(read.delim(file, comment.char='#'))
				cat(paste('\nignore logs for\n',paste(colnames(df)[ !grepl('state|POL|GAG|ENV|ucld|meanRate|coefficientOfVariation|treeModel.rootHeight',colnames(df)) ], collapse=', ') ))	
				df		<- subset(df, select=which(grepl('state|POL|GAG|ENV|ucld|meanRate|coefficientOfVariation|treeModel.rootHeight',colnames(df))))
				#log.df	<- c( paste('frequencies',1:4,sep=''), 'mu', 'alpha', '\\.at', '\\.ac', '\\.cg', '\\.ag', '\\.gt' )
				log.df	<- c( 'mu', 'alpha', '\\.at', '\\.ac', '\\.cg', '\\.ag', '\\.gt' )
				log.df	<- lapply( log.df, function(x)
						{
							tmp		<- melt( subset(df, select=which(grepl(paste('state|',x,sep=''),colnames(df)))), id='state', value.name=x)
							tmp[, GENE:= tmp[,substr(variable, 1, 3)]]
							tmp[, CODON_POS:=tmp[, regmatches(variable, regexpr('CP[1-3]',variable))]]
							subset(tmp, select=which(colnames(tmp)!='variable'))				
						})
				tmp		<- log.df[[1]]
				for(j in seq_along(log.df)[-1])
					tmp	<- merge(tmp, log.df[[j]], by=c('state','GENE','CODON_POS'))
				log.df	<- tmp	
				log.df	<- merge(log.df, subset(df, select=which(grepl('state|ucld|meanRate|coefficientOfVariation|treeModel.rootHeight',colnames(df)))), by='state')
				log.df[, FILE:= regmatches(infile, regexpr('pool[0-9]+',infile))]
				log.df
			})
	log.df	<- do.call('rbind',log.df)	
	setnames(log.df, colnames(log.df), gsub('\\.','',colnames(log.df),fixed=TRUE))
	#	handle empirical freq
	if(any(grepl('frequencies',colnames(log.df))))
		setnames(log.df, paste('frequencies',1:4,sep=''), c('a','c','g','t') )
	if(all(!grepl('frequencies',colnames(log.df))))
	{
		log.df	<- merge(log.df, freq, by=c('GENE','CODON_POS','FILE'))
	}
	log.df	<- subset(log.df, state>tree.id.burnin)
	#	check mean of relative rates
	tmp	<- log.df[, list(mmu=mean(mu), n=length(mu)), by=c('FILE','GENE', 'state')]
	stopifnot( nrow(subset(tmp, abs(mmu-1)>2*EPS))==0 )
	#
	file	<- paste( substr(infiles[1],1,nchar(infiles[1])-9),'log.R',sep='' )
	file	<- "PANGEA_SSAfgBwhRc-_140907_n390_log.R"
	file	<- paste( outdir, '/', file, sep='')
	cat(paste('\nsave to file', file))
	save(file=file, log.df)
	# 
	tmp		<- copy(log.df)
	tmp		<- melt(tmp, id=c('state','GENE','CODON_POS','FILE'))
	ggplot( tmp, aes(x=value, fill=FILE)) + geom_histogram() + facet_grid(GENE+CODON_POS~variable, scales='free')
	file	<- paste( substr(file,1,nchar(file)-1),'pdf',sep='' )
	ggsave(file, h=15, w=20)
}
##--------------------------------------------------------------------------------------------------------
##	get ancestral sequences from BEAST XML
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.BEAST.SSApg.getancestralseq.from.output<- function()
{
	tree.id.burnin		<- 2e7
	tree.id.labelsep	<- '|'
	ancseq.excl.timediff<- 3
	
	dir.name			<- '/Users/Oliver/duke/2014_Gates'  	
	indir				<- paste(dir.name,'methods_comparison_rootseqsim/140907',sep='/')
	ancseq.label.prefix	<- 'PANGEA_SSApgBwhRc-_140907_n390'
	outdir				<- indir
	#	search for BEAST output
	files				<- list.files(indir)
	files				<- files[ sapply(files, function(x) grepl('pool[0-9].R$',x) ) ]	
	if(!length(files))	stop('cannot find files matching criteria')
	
	#	load and process BEAST PARSER output
	#	sampling times are different for each gene, as they come from different trees
	anc.seq				<- lapply(files, function(file)
			{
				cat(paste('\nProcess file=', file  ))
				load( paste(indir, file, sep='/') )	#	expect tree, node.stat
				#	compute gag pol env ancestral sequences		
				anc.seq	<- PANGEA.RootSeqSim.get.ancestral.seq.pg(tree, node.stat, tree.id.sep='_', tree.id.idx.mcmcit=2, tree.id.burnin=tree.id.burnin, label.sep=tree.id.labelsep, label.idx.ctime=5)								
				anc.seq[, POOL:= regmatches(file, regexpr('pool[0-9]+', file)) ]
				set(anc.seq, NULL, 'BEAST_MCMC_IT', NULL )
				anc.seq
			})
	anc.seq				<- do.call('rbind',anc.seq)
	set(anc.seq, NULL, 'POOL', anc.seq[, factor(POOL)])
	set(anc.seq, NULL, 'GENE', anc.seq[, factor(GENE)])
	set(anc.seq, NULL, 'TREE_ID', anc.seq[, factor(TREE_ID)])
	anc.seq[, LABEL:=NULL]
	#	check we have exactly 3 genes for every inner node	
	tmp	<- anc.seq[, list(n=length(GENE)), by=c('POOL','TREE_ID','NODE_ID')]
	stopifnot(tmp[, all(n==3)])
	# prelim save
	file				<- paste( outdir, '/', substr(files[1],1,nchar(files[1])-7), 'AncSeq_Raw.R',sep='' )
	save(anc.seq, file=file)	
	#	concatenate sequences in time order (so we don t loose too many sequences)
	anc.seq.gag			<- subset(anc.seq, GENE=='GAG')
	setnames(anc.seq.gag, c('SEQ','CALENDAR_TIME','POOL'), c('GAG','GAG_CALENDAR_TIME','GAG_POOL'))	
	anc.seq.pol			<- subset(anc.seq, GENE=='POL')
	setnames(anc.seq.pol, c('SEQ','CALENDAR_TIME','POOL'), c('POL','POL_CALENDAR_TIME','POL_POOL'))
	anc.seq.env			<- subset(anc.seq, GENE=='ENV')
	setnames(anc.seq.env, c('SEQ','CALENDAR_TIME','POOL'), c('ENV','ENV_CALENDAR_TIME','ENV_POOL'))
	setkey(anc.seq.gag, GAG_CALENDAR_TIME)
	setkey(anc.seq.pol, POL_CALENDAR_TIME)	
	setkey(anc.seq.env, ENV_CALENDAR_TIME)
	#	
	anc.seq				<- cbind( subset(anc.seq.gag, select=c(GAG, GAG_CALENDAR_TIME)), subset(anc.seq.pol, select=c(POL, POL_CALENDAR_TIME)) )
	anc.seq				<- cbind( anc.seq, subset(anc.seq.env, select=c(ENV, ENV_CALENDAR_TIME) ))
	cat(paste('\nFound starting sequences, n=', nrow(anc.seq)))
	#	exclude concatenated genes if TIME_SEQ difference too large
	anc.seq[, d.gp:= abs(GAG_CALENDAR_TIME-POL_CALENDAR_TIME)]
	anc.seq[, d.ge:= abs(GAG_CALENDAR_TIME-ENV_CALENDAR_TIME)]
	anc.seq[, d.pe:= abs(POL_CALENDAR_TIME-ENV_CALENDAR_TIME)]
	anc.seq				<- subset(anc.seq, d.gp<=ancseq.excl.timediff & d.ge<=ancseq.excl.timediff & d.pe<=ancseq.excl.timediff)
	cat(paste('\nKeep starting sequences with sufficiently close TIME_SEQ, n=', nrow(anc.seq)))
	anc.seq[, d.gp:=NULL]
	anc.seq[, d.ge:=NULL]
	anc.seq[, d.pe:=NULL]
	#	finalize
	anc.seq[, CALENDAR_TIME:= (GAG_CALENDAR_TIME+POL_CALENDAR_TIME+ENV_CALENDAR_TIME)/3]
	anc.seq[, GAG_CALENDAR_TIME:=NULL]
	anc.seq[, POL_CALENDAR_TIME:=NULL]
	anc.seq[, ENV_CALENDAR_TIME:=NULL]
	anc.seq[, LABEL:= paste(ancseq.label.prefix, tree.id.labelsep, 'SEQ_', seq_len(nrow(anc.seq)), tree.id.labelsep, round(CALENDAR_TIME, d=4), sep='')]
	anc.seq		<- subset(anc.seq, CALENDAR_TIME>1935)
	anc.seq.gag	<- anc.seq.pol	<- anc.seq.env	<- NULL
	gc()	
	#
	#	plot distribution of node times
	#
	tmp			<- subset(anc.seq, select=c(CALENDAR_TIME))
	ggplot(tmp, aes(x=CALENDAR_TIME)) + geom_histogram(binwidth=1) + scale_x_continuous(breaks=seq(1900,2020, 5), name='estimated inner node time\n(year)')
	file		<- paste( outdir, '/', substr(files[1],1,nchar(files[1])-7), 'AncSeq_Times.pdf',sep='' )
	ggsave(file=file, w=10, h=6)	
	#
	#	return DNAbin
	#
	tmp			<- c(seq(1, nrow(anc.seq), 5e4), nrow(anc.seq)+1)
	anc.seq.gag	<- lapply(seq_along(tmp)[-length(tmp)], function(i)
			{
				cat(paste('\nProcess GAG up to',tmp[i+1]-1))
				anc.seq.gag				<- tolower(do.call('rbind',strsplit(anc.seq[seq.int(tmp[i], tmp[i+1]-1), GAG],'')))
				rownames(anc.seq.gag)	<- anc.seq[seq.int(tmp[i], tmp[i+1]-1), LABEL]
				anc.seq.gag				<- as.DNAbin(anc.seq.gag)				
			})
	anc.seq.gag	<- do.call('rbind', anc.seq.gag)
	anc.seq.pol	<- lapply(seq_along(tmp)[-length(tmp)], function(i)
			{
				cat(paste('\nProcess POL up to',tmp[i+1]-1))
				anc.seq.pol				<- tolower(do.call('rbind',strsplit(anc.seq[seq.int(tmp[i], tmp[i+1]-1), POL],'')))
				rownames(anc.seq.pol)	<- anc.seq[seq.int(tmp[i], tmp[i+1]-1), LABEL]
				anc.seq.pol				<- as.DNAbin(anc.seq.pol)				
			})
	anc.seq.pol	<- do.call('rbind', anc.seq.pol)
	anc.seq.env	<- lapply(seq_along(tmp)[-length(tmp)], function(i)
			{
				cat(paste('\nProcess ENV up to',tmp[i+1]-1))
				anc.seq.env				<- tolower(do.call('rbind',strsplit(anc.seq[seq.int(tmp[i], tmp[i+1]-1), ENV],'')))
				rownames(anc.seq.env)	<- anc.seq[seq.int(tmp[i], tmp[i+1]-1), LABEL]
				anc.seq.env				<- as.DNAbin(anc.seq.env)				
			})
	anc.seq.env	<- do.call('rbind', anc.seq.env)
	#
	#	return info data.table
	#
	set( anc.seq, NULL, 'GAG', NULL )
	set( anc.seq, NULL, 'POL', NULL )
	set( anc.seq, NULL, 'ENV', NULL )
	anc.seq.info		<- anc.seq
	#	save
	file				<- "/Users/Oliver/duke/2014_Gates/methods_comparison_rootseqsim/140907/PANGEA_SSAfgBwhRc-_140907_n390_AncSeq.R"		
	cat(paste('\nwrite Ancestral Sequences to ',file))
	save(file=file, anc.seq.gag, anc.seq.pol, anc.seq.env, anc.seq.info)
}
##--------------------------------------------------------------------------------------------------------
##	get anecestral sequences from BEAST XML
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.BEAST.SSAfg.getancestralseq.from.output<- function()
{
	tree.id.burnin		<- 2e7
	tree.id.labelsep	<- '|'
	dir.name			<- '/Users/Oliver/duke/2014_Gates'  	
	indir				<- paste(dir.name,'methods_comparison_rootseqsim/140813',sep='/')
	outdir				<- indir
	#	search for BEAST output
	files				<- list.files(indir)
	files				<- files[ sapply(files, function(x) grepl('pool[0-9].R$',x) ) ]	
	if(!length(files))	stop('cannot find files matching criteria')
	
	#	load and process BEAST PARSER output
	anc.seq				<- lapply(files, function(file)
			{
				cat(paste('\nProcess file=', file  ))
				load( paste(indir, file, sep='/') )	#	expect tree, node.stat
				#	compute gag pol env ancestral sequences		
				anc.seq	<- PANGEA.RootSeqSim.get.ancestral.seq(tree, node.stat, tree.id.sep='_', tree.id.idx.mcmcit=2, tree.id.burnin=tree.id.burnin, label.sep=tree.id.labelsep, label.idx.ctime=5)				
				set(anc.seq, NULL, 'LABEL', anc.seq[, paste( substr(file,1,nchar(file)-2), LABEL, sep=tree.id.labelsep )] )				
				set(anc.seq, NULL, 'TREE_ID', NULL )
				set(anc.seq, NULL, 'NODE_ID', NULL )
				set(anc.seq, NULL, 'BEAST_MCMC_IT', NULL )
				anc.seq
			})
	anc.seq				<- do.call('rbind',anc.seq)
	#
	#	return DNAbin
	#
	anc.seq.gag				<- tolower(do.call('rbind',strsplit(anc.seq[, GAG],'')))
	rownames(anc.seq.gag)	<- anc.seq[, LABEL]
	anc.seq.gag				<- as.DNAbin(anc.seq.gag)		
	anc.seq.pol				<- tolower(do.call('rbind',strsplit(anc.seq[, POL],'')))
	rownames(anc.seq.pol)	<- anc.seq[, LABEL]
	anc.seq.pol				<- as.DNAbin(anc.seq.pol)		
	anc.seq.env				<- tolower(do.call('rbind',strsplit(anc.seq[, ENV],'')))
	rownames(anc.seq.env)	<- anc.seq[, LABEL]
	anc.seq.env				<- as.DNAbin(anc.seq.env)	
	
	set( anc.seq, NULL, 'GAG', NULL )
	set( anc.seq, NULL, 'POL', NULL )
	set( anc.seq, NULL, 'ENV', NULL )
	anc.seq.info			<- anc.seq
	#anc.seq					<- cbind(anc.seq.gag, anc.seq.pol, anc.seq.env)
	#
	outfile				<- paste( substr(files[1],1,nchar(files[1])-7), 'AncSeq.R',sep='' )
	file				<- paste(outdir, outfile, sep='/')
	cat(paste('\nwrite Ancestral Sequences to ',file))
	save(file=file, anc.seq.gag, anc.seq.pol, anc.seq.env, anc.seq.info)
}
##--------------------------------------------------------------------------------------------------------
##	check simulated sequences: create NJ tree and estimate R2
##	olli 14.09.14
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.SSApg.NJR2<- function()
{
	require(phytools)
	tree.id.labelsep		<- '|'
	tree.id.label.idx.ctime	<- 4 
	indir		<- '/Users/Oliver/git/HPTN071sim/tmp140914/140716_RUN001_INTERNAL'  
	outdir		<- '/Users/Oliver/git/HPTN071sim/tmp140914/140716_RUN001_INTERNAL'
	#indir		<- '/Users/Oliver/git/HPTN071sim/tmp140912/140911_DSPS_RUN002_INTERNAL'  
	#outdir		<- '/Users/Oliver/git/HPTN071sim/tmp140912/140911_DSPS_RUN002_INTERNAL'
	
	infile		<- list.files(indir, '.*INTERNAL.R$', full.names=FALSE)
	stopifnot(length(infile)==1)
	#	load simulated data
	file			<- paste(indir, '/', infile, sep='')
	cat(paste('\nLoading file', file))
	load(file)		#expect "df.epi"    "df.trms"   "df.inds"   "df.sample" "df.seq"
	#	load outgroup sequences
	file			<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfg_HXB2outgroup.R')
	cat(paste('\nLoading outgroup seq from file', file))
	load(file)		#expect "outgroup.seq.gag" "outgroup.seq.pol" "outgroup.seq.env"
	
	tmp				<- tolower(do.call('rbind',strsplit(df.seq[, GAG],'')))
	rownames(tmp)	<- df.seq[, LABEL]
	df.seq.gag		<- as.DNAbin(tmp)
	tmp				<- tolower(do.call('rbind',strsplit(df.seq[, POL],'')))
	rownames(tmp)	<- df.seq[, LABEL]
	df.seq.pol		<- as.DNAbin(tmp)	
	tmp				<- tolower(do.call('rbind',strsplit(df.seq[, ENV],'')))
	rownames(tmp)	<- df.seq[, LABEL]
	df.seq.env		<- as.DNAbin(tmp)
	#
	#	get R2 for df.seq.pol
	#
	seq				<- df.seq.pol
	seq				<- rbind(seq, outgroup.seq.pol[, seq_len(ncol(seq))])
	#	get NJ tree	
	tmp				<- dist.dna(seq)
	nj				<- nj(tmp)
	tmp				<- which(nj$tip.label=="HXB2")
	nj				<- reroot(nj, tmp, nj$edge.length[which(nj$edge[,2]==tmp)])
	nj				<- ladderize(nj)		
	file			<- paste( outdir, '/', substr(infile,1,nchar(infile)-20),'INFO_simu_NJpol.pdf', sep='' )	
	pdf(file=file, w=10, h=150)
	plot(nj, show.tip=TRUE, cex=0.5)
	add.scale.bar()
	dev.off()			
	#	get root to tip divergence
	nj				<- drop.tip(nj,'HXB2')
	tmp				<- node.depth.edgelength(nj)
	nj.info			<- data.table(LABEL=nj$tip.label, ROOT2TIP=tmp[seq_len(Ntip(nj))] )
	set(nj.info, NULL, 'CALENDAR_TIME', nj.info[, as.numeric(sapply(strsplit(LABEL, tree.id.labelsep, fixed=TRUE),'[[',tree.id.label.idx.ctime))] )
	tmp				<- lm(ROOT2TIP~CALENDAR_TIME, data=nj.info)		 
	set( nj.info, NULL, 'ROOT2TIP_LM', predict(tmp, type='response') ) 	
	tmp2			<- c( R2=round(summary(tmp)$r.squared,d=3), SLOPE= as.numeric(round(coef(tmp)['CALENDAR_TIME'],d=4)), TMRCA=as.numeric(round( -coef(tmp)['(Intercept)']/coef(tmp)['CALENDAR_TIME'], d=1 )) )
	ggplot(nj.info, aes(x=CALENDAR_TIME, y=ROOT2TIP)) + geom_point(alpha=0.5) + geom_line(aes(y=ROOT2TIP_LM)) +
			#scale_x_continuous(breaks=seq(1980,2020,2)) +						
			labs(x='Sequence sampling date', y='root-to-tip divergence\n(HIV-1 pol sequences)') +
			annotate("text", x=nj.info[, min(CALENDAR_TIME)], y=nj.info[, 0.9*max(ROOT2TIP)], label=paste("R2=", tmp2['R2'],'\nSlope=',tmp2['SLOPE'],'\nTMRCA=',tmp2['TMRCA'], sep=''), hjust = 0, size = 4) +
			theme(legend.position=c(0,1), legend.justification=c(0,1))		
	file			<- paste( outdir, '/', substr(infile,1,nchar(infile)-20),'INFO_simu_NJpolR2.pdf', sep='' )
	ggsave(file=file, w=10, h=6)
	#
	#	get R2 for concatenated genome
	#
	seq				<- cbind(df.seq.gag,df.seq.pol,df.seq.env)
	tmp				<- cbind(outgroup.seq.gag[,1:ncol(df.seq.gag)], outgroup.seq.pol, outgroup.seq.env)
	seq				<- rbind(seq,tmp)
	#	get NJ tree	
	tmp				<- dist.dna(seq)
	nj				<- nj(tmp)
	tmp				<- which(nj$tip.label=="HXB2")
	nj				<- reroot(nj, tmp, nj$edge.length[which(nj$edge[,2]==tmp)])
	nj				<- ladderize(nj)		
	file			<- paste( outdir, '/', substr(infile,1,nchar(infile)-20),'INFO_simu_NJ.pdf', sep='' )	
	pdf(file=file, w=10, h=150)
	plot(nj, show.tip=TRUE, cex=0.5)
	add.scale.bar()
	dev.off()		
	#	get root to tip divergence
	nj				<- drop.tip(nj,'HXB2')
	tmp				<- node.depth.edgelength(nj)
	nj.info			<- data.table(LABEL=nj$tip.label, ROOT2TIP=tmp[seq_len(Ntip(nj))] )
	set(nj.info, NULL, 'CALENDAR_TIME', nj.info[, as.numeric(sapply(strsplit(LABEL, tree.id.labelsep, fixed=TRUE),'[[',tree.id.label.idx.ctime))] )
	tmp				<- lm(ROOT2TIP~CALENDAR_TIME, data=nj.info)		 
	set( nj.info, NULL, 'ROOT2TIP_LM', predict(tmp, type='response') ) 	
	tmp2			<- c( R2=round(summary(tmp)$r.squared,d=3), SLOPE= as.numeric(round(coef(tmp)['CALENDAR_TIME'],d=4)), TMRCA=as.numeric(round( -coef(tmp)['(Intercept)']/coef(tmp)['CALENDAR_TIME'], d=1 )) )
	ggplot(nj.info, aes(x=CALENDAR_TIME, y=ROOT2TIP)) + geom_point(alpha=0.5) + geom_line(aes(y=ROOT2TIP_LM)) +
			#scale_x_continuous(breaks=seq(1980,2020,2)) +						
			labs(x='Sequence sampling date', y='root-to-tip divergence\n(HIV-1 concatenated sequences)') +
			annotate("text", x=nj.info[, min(CALENDAR_TIME)], y=nj.info[, 0.9*max(ROOT2TIP)], label=paste("R2=", tmp2['R2'],'\nSlope=',tmp2['SLOPE'],'\nTMRCA=',tmp2['TMRCA'], sep=''), hjust = 0, size = 4) +
			theme(legend.position=c(0,1), legend.justification=c(0,1))		
	file			<- paste( outdir, '/', substr(infile,1,nchar(infile)-20),'INFO_simu_NJR2.pdf', sep='' )
	ggsave(file=file, w=10, h=6)
}
##--------------------------------------------------------------------------------------------------------
##	check simulated sequences: create ExaML tree and estimate R2
##	olli 14.09.14
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline<- function()
{	
	indir			<- system.file(package="rPANGEAHIVsim", "misc")
	indir			<- ifelse(indir=='','/Users/Oliver/git/HPTN071sim/raw_trchain',indir)
	infile.ind		<- '140716_RUN001'
	infile.trm		<- '140716_RUN001'	
	
	pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
													s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.2, 
													epi.model='HPTN071', epi.dt=1/48, epi.import=0.1,
													v.N0tau=3.58e4, v.r=2, v.T50=-1,
													wher.mu=NA, wher.sigma=NA,
													bwerm.mu=1.5, bwerm.sigma=0.12,
													startseq.backdate=NA )	
	pipeline.vary	<- data.table(wher.mu=c(0.005, 0.004, 0.003), wher.sigma=c(0.8, 0.7, 0.6), label=c('-5','-4','-3'))				
	
	dummy			<- pipeline.vary[, {				
											set(pipeline.args, which( pipeline.args$stat=='wher.mu' ), 'v', as.character(wher.mu))
											set(pipeline.args, which( pipeline.args$stat=='wher.sigma' ), 'v', as.character(wher.sigma))
											print(pipeline.args)
											#	re-name the following:
											tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140914'
											tmpdir			<- paste(tmpdir,label,sep='')
											dir.create(tmpdir, showWarnings=FALSE)						
											#						
											file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
											file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
											file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
											system(file)
										}, by='label']										
	
}
##--------------------------------------------------------------------------------------------------------
##	check simulated sequences: create ExaML tree and estimate R2
##	olli 14.09.14
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.SSApg.BEAST<- function()
{
	require(phytools)
	require(hivclust)
	require(XML)
	tree.id.labelsep		<- '|'
	tree.id.label.idx.ctime	<- 4 
	indir		<- '/Users/Oliver/duke/2014_Gates/methods_comparison_pipeline/140914'
	#indir		<- '/Users/Oliver/git/HPTN071sim/tmp140914/140716_RUN001_INTERNAL'  
	outdir		<- indir
	infiles		<- list.files(indir, '.*INTERNAL.R$', full.names=FALSE)
	#stopifnot(length(infiles)==1)
	#	read BEAST template files	
	infile.beast.pol	<- system.file(package="rPANGEAHIVsim", "misc",'BEAST_template_vTESTpol.xml')
	bxml.template.pol	<- xmlTreeParse(infile.beast.pol, useInternalNodes=TRUE, addFinalizer = TRUE)
	#infile.beast.gag	<- '/Users/Oliver/git/HPTN071sim/data_rootseq/BEAST_template_v10gag.xml'
	#infile.beast.env	<- '/Users/Oliver/git/HPTN071sim/data_rootseq/BEAST_template_v10env.xml'
	#bxml.template.gag	<- xmlTreeParse(infile.beast.gag, useInternalNodes=TRUE, addFinalizer = TRUE)	
	#bxml.template.env	<- xmlTreeParse(infile.beast.env, useInternalNodes=TRUE, addFinalizer = TRUE)	
	#
	#	run ExaML 
	#
	for(i in seq_along(infiles))
	{
		infile			<- infiles[i]
		#	load simulated data
		file			<- paste(indir, '/', infile, sep='')
		cat(paste('\nLoading file', file))
		load(file)		#expect "df.epi"    "df.trms"   "df.inds"   "df.sample" "df.seq"
		#	load outgroup sequences
		file			<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfg_HXB2outgroup.R')
		cat(paste('\nLoading outgroup seq from file', file))
		load(file)		#expect "outgroup.seq.gag" "outgroup.seq.pol" "outgroup.seq.env"
		#	concatenate sequences
		tmp				<- tolower(do.call('rbind',strsplit(df.seq[, GAG],'')))
		rownames(tmp)	<- df.seq[, LABEL]
		df.seq.gag		<- as.DNAbin(tmp)
		tmp				<- tolower(do.call('rbind',strsplit(df.seq[, POL],'')))
		rownames(tmp)	<- df.seq[, LABEL]
		df.seq.pol		<- as.DNAbin(tmp)	
		tmp				<- tolower(do.call('rbind',strsplit(df.seq[, ENV],'')))
		rownames(tmp)	<- df.seq[, LABEL]
		df.seq.env		<- as.DNAbin(tmp)
		seq				<- cbind(df.seq.gag,df.seq.pol,df.seq.env)
		tmp				<- cbind(outgroup.seq.gag[,1:ncol(df.seq.gag)], outgroup.seq.pol, outgroup.seq.env)
		seq				<- rbind(seq,tmp)
		#	get 100 'divergent' sequences from different clusters
		tmp				<- dist.dna( seq )
		seq.ph			<- nj(tmp)				
		tmp				<- hivc.clu.brdist.stats(seq.ph, eval.dist.btw="leaf", stat.fun=hivc.clu.min.transmission.cascade)
		thresh.brl		<- 0.045
		clustering		<- hivc.clu.clusterbythresh(seq.ph, thresh.brl=thresh.brl, dist.brl=tmp, retval="all")
		cat(paste('\nFound clusters, n=', length(clustering$clu.idx)))
		#	Take 1 sequence from each cluster
		seq.select		<- data.table( PH_NODE_ID=seq_len(Ntip(seq.ph)), CLU_ID=clustering$clu.mem[ seq_len(Ntip(seq.ph)) ] )
		seq.select		<- subset(seq.select, !is.na(CLU_ID))[, list(LABEL= seq.ph$tip.label[PH_NODE_ID[1]]), by='CLU_ID']
		seq.select		<- merge(df.seq, seq.select, by='LABEL')
		#
		#	create BEAST XML
		#
		#
		#	POL
		#
		cat(paste('\ncreate POL BEAST XML file for seqs=',paste( seq.select[,LABEL], collapse=' ')))
		pool.infile		<- paste(  substr(infile,1,nchar(infile)-21),'_TEST_pol', sep='' )
		#	write XML file with new sequences
		bxml			<- newXMLDoc(addFinalizer=T)
		bxml.beast		<- newXMLNode("beast", doc=bxml, addFinalizer=T)
		tmp				<- newXMLCommentNode(text=paste("Generated by HIVCLUST from template",infile.beast.pol), parent=bxml.beast, doc=bxml, addFinalizer=T)
		#	add new set of GAG sequences into GAG alignment
		tmp				<- tolower(do.call('rbind',strsplit(seq.select[, POL],'')))
		rownames(tmp)	<- seq.select[, LABEL]
		tmp				<- as.DNAbin(tmp)
		bxml			<- hivc.beast.add.seq(bxml, tmp, df=NULL, beast.label.datepos=4, beast.label.sep= '|', beast.date.direction= "forwards", beast.date.units= "years", beast.alignment.id="POL.alignment", beast.alignment.dataType= "nucleotide", verbose=1)
		#	copy from template	
		bt.beast		<- getNodeSet(bxml.template.pol, "//beast")[[1]]
		dummy			<- sapply(seq.int( 1, xmlSize(bt.beast) ), function(i)
				{
					if( class(bt.beast[[i]])[1]=="XMLInternalCommentNode" )
						dummy<- newXMLCommentNode(text=xmlValue(bt.beast[[i]]), parent=bxml.beast, doc=bxml, addFinalizer=T)
					else
						dummy<- addChildren( bxml.beast, xmlClone( bt.beast[[i]], addFinalizer=T, doc=bxml ) )
				})
		#	change gmrf dimensions	
		tmp			<- getNodeSet(bxml, "//*[@id='skyride.logPopSize']")
		if(length(tmp)!=1)	stop("unexpected number of *[@id='skyride.logPopSize'")
		tmp			<- tmp[[1]]
		xmlAttrs(tmp)["dimension"]	<-	nrow(seq.select)-1  
		tmp			<- getNodeSet(bxml, "//*[@id='skyride.groupSize']")
		if(length(tmp)!=1)	stop("unexpected number of *[@id='skyride.groupSize'")
		tmp			<- tmp[[1]]
		xmlAttrs(tmp)["dimension"]	<-	nrow(seq.select)-1			
		#	change outfile name 
		bxml.onodes	<- getNodeSet(bxml, "//*[@fileName]")
		tmp			<- sapply(bxml.onodes, function(x) xmlGetAttr(x,"fileName"))
		tmp			<- gsub("(time).","time",tmp,fixed=1)
		tmp			<- gsub("(subst).","subst",tmp,fixed=1)	
		tmp			<- sapply(strsplit(tmp,'.',fixed=1), function(x)	paste(pool.infile, '.', tail(x,1), sep=''))
		dummy		<- sapply(seq_along(bxml.onodes), function(i){		xmlAttrs(bxml.onodes[[i]])["fileName"]<- tmp[i]		})
		#	write to file
		file		<- paste(indir,'/',pool.infile,".xml", sep='')
		cat(paste("\nwrite xml file to",file))
		saveXML(bxml, file=file)
	}
}
##--------------------------------------------------------------------------------------------------------
##	check simulated sequences: create ExaML tree and estimate R2
##	olli 14.09.14
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.SSApg.ExaMLR2<- function()
{
	require(phytools)
	require(hivclust)
	tree.id.labelsep		<- '|'
	tree.id.label.idx.ctime	<- 4 
	indir		<- '/Users/Oliver/git/HPTN071sim/tmp140910/140716_RUN001_INTERNAL'  
	outdir		<- '/Users/Oliver/git/HPTN071sim/tmp140910/140716_RUN001_INTERNAL'
	infiles		<- list.files(indir, '.*INTERNAL.R$', full.names=FALSE)
	stopifnot(length(infiles)==1)
	#
	#	run ExaML 
	#
	for(i in seq_along(infiles))
	{
		infile		<- infiles[i]
		#	load simulated data
		file			<- paste(indir, '/', infile, sep='')
		cat(paste('\nLoading file', file))
		load(file)		#expect "df.epi"    "df.trms"   "df.inds"   "df.sample" "df.seq"
		#	load outgroup sequences
		file			<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfg_HXB2outgroup.R')
		cat(paste('\nLoading outgroup seq from file', file))
		load(file)		#expect "outgroup.seq.gag" "outgroup.seq.pol" "outgroup.seq.env"
		
		tmp				<- tolower(do.call('rbind',strsplit(df.seq[, GAG],'')))
		rownames(tmp)	<- df.seq[, LABEL]
		df.seq.gag		<- as.DNAbin(tmp)
		tmp				<- tolower(do.call('rbind',strsplit(df.seq[, POL],'')))
		rownames(tmp)	<- df.seq[, LABEL]
		df.seq.pol		<- as.DNAbin(tmp)	
		tmp				<- tolower(do.call('rbind',strsplit(df.seq[, ENV],'')))
		rownames(tmp)	<- df.seq[, LABEL]
		df.seq.env		<- as.DNAbin(tmp)
		#
		#	run ExaML on pol
		#
		seq				<- df.seq.pol
		seq				<- rbind(seq, outgroup.seq.pol[, seq_len(ncol(seq))])
		infile.seq.sig	<- "Sun_Sep_14_12:59:06_2013"
		infile.seq		<- paste(substr(infile,1,nchar(infile)-20),'INFO_simu_polseq',sep='')
		file			<- paste( outdir, '/', infile.seq,'_',gsub('/',':',infile.seq.sig),'.R', sep='' )
		save(seq, file=file)
		#	run ExaML
		cmd				<- hivc.cmd.examl.bootstrap.on.one.machine(indir, infile.seq, infile.seq.sig, infile.seq.sig, bs.from=0, bs.to=0, verbose=1)
		cmd				<- hivc.cmd.hpcwrapper(cmd, hpc.walltime=21, hpc.q= NA, hpc.mem="450mb", hpc.nproc=1)
		cmd.hpccaller(outdir, paste("exa",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.'), cmd)
		Sys.sleep(1)	
		#
		#	run ExaML on concatenated
		#
		seq				<- cbind(df.seq.gag,df.seq.pol,df.seq.env)
		tmp				<- cbind(outgroup.seq.gag[,1:ncol(df.seq.gag)], outgroup.seq.pol, outgroup.seq.env)
		seq				<- rbind(seq,tmp)
		infile.seq.sig	<- "Sun_Sep_14_12:59:06_2013"
		infile.seq		<- paste(substr(infile,1,nchar(infile)-20),'INFO_simu_concseq',sep='')
		file			<- paste( outdir, '/', infile.seq,'_',gsub('/',':',infile.seq.sig),'.R', sep='' )
		save(seq, file=file)
		#	run ExaML
		cmd				<- hivc.cmd.examl.bootstrap.on.one.machine(indir, infile.seq, infile.seq.sig, infile.seq.sig, bs.from=0, bs.to=0, verbose=1)
		cmd				<- hivc.cmd.hpcwrapper(cmd, hpc.walltime=21, hpc.q= NA, hpc.mem="450mb", hpc.nproc=1)
		cmd.hpccaller(outdir, paste("exa",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.'), cmd)
		Sys.sleep(1)	
	}	
	#
	#	evaluate R2 for pol
	#
	infiles			<- list.files(indir, '^ExaML_result.*polseq.*finaltree.000$', full.names=FALSE)
	i				<- 1
	infile			<- infiles[i]
	file			<- paste(indir,'/',infile,sep='')
	ph				<- read.tree(file)
	
	tmp				<- which(ph$tip.label=="HXB2")
	ph				<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])
	ph				<- ladderize(ph)		
	file			<- paste( outdir, '/', substr(infile,1,nchar(infile)-20),'INFO_simu_ExaMLpol.pdf', sep='' )	
	pdf(file=file, w=10, h=150)
	plot(ph, show.tip=TRUE, cex=0.5)
	add.scale.bar()
	dev.off()			
	#	get root to tip divergence
	ph				<- drop.tip(ph,'HXB2')
	tmp				<- node.depth.edgelength(ph)
	ph.info			<- data.table(LABEL=ph$tip.label, ROOT2TIP=tmp[seq_len(Ntip(ph))] )
	set(ph.info, NULL, 'CALENDAR_TIME', ph.info[, as.numeric(sapply(strsplit(LABEL, tree.id.labelsep, fixed=TRUE),'[[',tree.id.label.idx.ctime))] )
	tmp				<- lm(ROOT2TIP~CALENDAR_TIME, data=ph.info)		 
	set( ph.info, NULL, 'ROOT2TIP_LM', predict(tmp, type='response') ) 	
	tmp2			<- c( R2=round(summary(tmp)$r.squared,d=3), SLOPE= as.numeric(round(coef(tmp)['CALENDAR_TIME'],d=4)), TMRCA=as.numeric(round( -coef(tmp)['(Intercept)']/coef(tmp)['CALENDAR_TIME'], d=1 )) )
	ggplot(ph.info, aes(x=CALENDAR_TIME, y=ROOT2TIP)) + geom_point(alpha=0.5) + geom_line(aes(y=ROOT2TIP_LM)) +
			#scale_x_continuous(breaks=seq(1980,2020,2)) +						
			labs(x='Sequence sampling date', y='root-to-tip divergence\n(HIV-1 pol sequences)') +
			annotate("text", x=ph.info[, min(CALENDAR_TIME)], y=ph.info[, 0.9*max(ROOT2TIP)], label=paste("R2=", tmp2['R2'],'\nSlope=',tmp2['SLOPE'],'\nTMRCA=',tmp2['TMRCA'], sep=''), hjust = 0, size = 4) +
			theme(legend.position=c(0,1), legend.justification=c(0,1))		
	file			<- paste( outdir, '/', substr(infile,1,nchar(infile)-20),'INFO_simu_ExaMLpolR2.pdf', sep='' )
	ggsave(file=file, w=10, h=6)
	
	
}
##--------------------------------------------------------------------------------------------------------
##	check ancestral sequences from BEAST XML, create random draw to check
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.SIMU.SSAfg.checkancestralseq.createdataset<- function()
{
	require(phytools)
	
	tree.id.burnin			<- 2e7
	tree.id.labelsep		<- '|'
	tree.id.label.idx.ctime	<- 4
	dir.name				<- '/Users/Oliver/duke/2014_Gates'  	
	indir					<- paste(dir.name,'methods_comparison_rootseqsim/140813',sep='/')	
	infile					<- "PANGEA_SSAfgBwhRc-_140811_n390_AncSeq.R"	
	outdir					<- indir
	outfile.fg				<- infile
	outfile.partialpol		<- "PANGEA_SSApolBwhRc-_140811_n390_AncSeq.R"
	outsignat				<- "Mon_Aug_17_17:05:23_2014"	
	#	load ancestral sequences
	load( paste(indir, infile, sep='/') ) #	expect anc.seq.gag, anc.seq.pol, anc.seq.env, anc.seq.info
	#	load aligned HXB2 as outgroup
	load('/Users/Oliver/git/HPTN071sim/data_rootseq/PANGEA_SSAfg_140806_HXB2outgroup.R')	#expect "outgroup.seq.gag" "outgroup.seq.pol" "outgroup.seq.env"
	#
	#	basic checks
	#	
	if(0)
	{
		tmp					<- anc.seq.info[, 	sapply( strsplit(LABEL, tree.id.labelsep, fixed=TRUE),'[[', 1) ]
		anc.seq.info[, POOL:= regmatches(tmp, regexpr('pool[0-9]+', tmp))]
		set(anc.seq.info, NULL, 'POOL', anc.seq.info[, as.numeric(substr(POOL, 5, nchar(POOL)))])
		ggplot(anc.seq.info, aes(x=CALENDAR_TIME)) + geom_histogram(binwidth=2) + facet_grid(.~POOL, margins=0)
		file				<- paste(indir, '/', substr(infile,1,nchar(infile)-2), '_calendartime.pdf' , sep='')
		ggsave(file=file, w=8, h=6)
		#	most sequences between 1940 - 1980
		subset(anc.seq.info, floor(CALENDAR_TIME)==1940)		
	}
	#	sample data set to run ExaML
	anc.seq.info[, CALENDAR_YR:=anc.seq.info[, floor(CALENDAR_TIME)]]	
	#	~ 1400 sequences from 1940
	#	sample 10 times 1e3 sequences randomly from exp increasing prevalence and estimate tree to calculate root to tip divergence + clustering on fg and partial pol
	s.seed						<- 42
	s.PREV.MAX					<- 0.25
	s.PREV.MIN					<- 0.01
	s.RANGE						<- 40
	s.size						<- 1e3
	s.baseline.ancseq.time		<- 1940
	s.baseline.calendar.time	<- 1980
	s.LENGTH.PARTIAL.POL		<- 1500
	tree.id.labelidx.ctime		<- 4		
	tmp							<- log( 1+s.PREV.MAX-s.PREV.MIN ) / (s.RANGE-1)
	tmp							<- exp( tmp*seq.int(0,s.RANGE-1) ) - 1 + s.PREV.MIN
	seq.s						<- c(s.PREV.MIN, diff(tmp))
	seq.s						<- data.table( SEQ_N= round( seq.s/s.PREV.MAX*s.size ), ANCSEQ_YR=seq_along(seq.s)+s.baseline.ancseq.time-1 )
	anc.seq.info				<- subset(anc.seq.info,  CALENDAR_YR>=s.baseline.ancseq.time & CALENDAR_YR<=(s.baseline.ancseq.time+s.RANGE-1))
	
	#	draw partial pol genome sequences - length is first 1500 sites
	#	and draw full genome sequences
	set.seed(s.seed)	
	for(check.draw in 1:5)
	{
		cat(paste('\nprocess checkdraw',check.draw))
		#	draw a large enough number of ancseq labels from the 3 pools for each year to accommodate SEQ_N/3 anc seqs from each pool
		anc.seq.infodraw			<- anc.seq.info[, {
					tmp	<- seq.s$SEQ_N[ which(seq.s$ANCSEQ_YR==CALENDAR_YR) ]
					list(LABEL= sample(LABEL, 2*tmp, replace=FALSE), SEQ_N=tmp, SEQ_N_GRACE=2*tmp)
				}, by='CALENDAR_YR']	
		anc.seq.draw				<- anc.seq.pol[anc.seq.infodraw[, LABEL], seq_len(s.LENGTH.PARTIAL.POL)]
		#
		#	make sure the drawn sequences are unique on partial POL
		#
		anc.seq.draw				<- seq.unique(anc.seq.draw)
		anc.seq.infodraw			<- merge( data.table(LABEL=rownames(anc.seq.draw)), anc.seq.infodraw, by='LABEL' )	
		stopifnot( anc.seq.infodraw[, list(SEQ_N_GRACE=length(LABEL), SEQ_N=SEQ_N[1]), by='CALENDAR_YR'][, all(SEQ_N_GRACE>=SEQ_N)] )
		anc.seq.infodraw			<- anc.seq.infodraw[, list(LABEL= LABEL[seq_len(SEQ_N[1])]), by='CALENDAR_YR']
		#	set new calendar time for sequences
		set(anc.seq.infodraw, NULL, 'LABEL_NEW', anc.seq.infodraw[, as.numeric( sapply( strsplit(LABEL,tree.id.labelsep,fixed=TRUE), '[[', tree.id.labelidx.ctime) ) ])
		set(anc.seq.infodraw, NULL, 'LABEL_NEW', anc.seq.infodraw[, LABEL_NEW-s.baseline.ancseq.time+s.baseline.calendar.time])	
		anc.seq.infodraw			<- anc.seq.infodraw[,	{
					tmp							<- strsplit(LABEL,tree.id.labelsep,fixed=TRUE)[[1]]
					tmp[tree.id.labelidx.ctime]	<- LABEL_NEW
					list(LABEL_NEW=paste(tmp, collapse=tree.id.labelsep,sep=''))
				}, by='LABEL']
		setkey(anc.seq.infodraw, LABEL)
		#
		#	select partial POL seqs
		#
		anc.seq.draw				<- anc.seq.pol[anc.seq.infodraw[, LABEL], seq_len(s.LENGTH.PARTIAL.POL)]		
		rownames(anc.seq.draw)		<- anc.seq.infodraw[ rownames(anc.seq.draw), ][, LABEL_NEW]
		anc.seq.draw				<- rbind(anc.seq.draw, outgroup.seq.pol[, seq_len(s.LENGTH.PARTIAL.POL)])
		#	save to file, including outgroup
		file						<- paste( outdir, '/', substr(outfile.partialpol,1,nchar(outfile.partialpol)-2),'_checkdraw', check.draw,'_', insignat, '.R', sep='' )
		cat(paste('\nsave to file', file))
		save(anc.seq.draw, file=file)		
		#	get NJ tree	with HXB2 outgroup		
		tmp				<- dist.dna(anc.seq.draw)
		anc.seq.nj		<- nj(tmp)
		tmp				<- which(anc.seq.nj$tip.label=="HXB2")
		anc.seq.nj		<- reroot(anc.seq.nj, tmp, anc.seq.nj$edge.length[which(anc.seq.nj$edge[,2]==tmp)])
		anc.seq.nj		<- ladderize(anc.seq.nj)		
		file			<- paste( outdir, '/', substr(outfile.partialpol,1,nchar(outfile.partialpol)-2),'_checkdraw', check.draw,'_NJ.pdf', sep='' )	
		pdf(file=file, w=10, h=80)
		plot(anc.seq.nj, show.tip=TRUE, cex=0.5)
		dev.off()			
		#	get root to tip divergence
		anc.seq.nj		<- drop.tip(anc.seq.nj,'HXB2')		
		tmp				<- node.depth.edgelength(anc.seq.nj)
		anc.seq.nj.info	<- data.table(LABEL=anc.seq.nj$tip.label, ROOT2TIP=tmp[seq_len(Ntip(anc.seq.nj))] )
		#tmp<- distRoot(anc.seq.nj)
		#set(anc.seq.nj.info, NULL, 'distRoot', tmp)
		set(anc.seq.nj.info, NULL, 'CALENDAR_TIME', anc.seq.nj.info[, as.numeric(sapply(strsplit(LABEL, tree.id.labelsep, fixed=TRUE),'[[',tree.id.label.idx.ctime))] )
		tmp				<- lm(ROOT2TIP~CALENDAR_TIME, data=anc.seq.nj.info)		 
		set( anc.seq.nj.info, NULL, 'ROOT2TIP_LM', predict(tmp, type='response') ) 	
		tmp2			<- c( R2=round(summary(tmp)$r.squared,d=3), SLOPE= as.numeric(round(coef(tmp)['CALENDAR_TIME'],d=4)), TMRCA=as.numeric(round( -coef(tmp)['(Intercept)']/coef(tmp)['CALENDAR_TIME'], d=1 )) )
		ggplot(anc.seq.nj.info, aes(x=CALENDAR_TIME, y=ROOT2TIP)) + geom_point(alpha=0.5) + geom_line(aes(y=ROOT2TIP_LM)) +
				#scale_x_continuous(breaks=seq(1980,2020,2)) +						
				labs(x='Sequence sampling date', y='root-to-tip divergence') +
				annotate("text", x=anc.seq.nj.info[, min(CALENDAR_TIME)], y=anc.seq.nj.info[, 0.9*max(ROOT2TIP)], label=paste("R2=", tmp2['R2'],'\nSlope=',tmp2['SLOPE'],'\nTMRCA=',tmp2['TMRCA'], sep=''), hjust = 0, size = 4) +
				theme(legend.position=c(0,1), legend.justification=c(0,1))
		file			<- paste( outdir, '/', substr(outfile.partialpol,1,nchar(outfile.partialpol)-2),'_checkdraw', check.draw,'_Root2Tip.pdf', sep='' )
		ggsave(file=file, w=10, h=6)		
		#
		#	select the same full genome seqs plus aligned HXB2 as outgroup 
		#
		anc.seq.draw				<- do.call( 'cbind', list( anc.seq.gag[anc.seq.infodraw[, LABEL], ], anc.seq.pol[anc.seq.infodraw[, LABEL], ], anc.seq.env[anc.seq.infodraw[, LABEL], ] ) ) 
		rownames(anc.seq.draw)		<- anc.seq.infodraw[ rownames(anc.seq.draw), ][, LABEL_NEW]		
		anc.seq.draw				<- rbind(anc.seq.draw, cbind(outgroup.seq.gag[,seq_len(ncol(anc.seq.gag))], outgroup.seq.pol, outgroup.seq.env))
		#	get NJ tree	
		tmp				<- dist.dna(anc.seq.draw)
		anc.seq.nj		<- nj(tmp)
		tmp				<- which(anc.seq.nj$tip.label=="HXB2")
		anc.seq.nj		<- reroot(anc.seq.nj, tmp, anc.seq.nj$edge.length[which(anc.seq.nj$edge[,2]==tmp)])
		anc.seq.nj		<- ladderize(anc.seq.nj)		
		file			<- paste( outdir, '/', substr(outfile.fg,1,nchar(outfile.fg)-2),'_checkdraw', check.draw,'_NJ.pdf', sep='' )	
		pdf(file=file, w=10, h=80)
		plot(anc.seq.nj, show.tip=TRUE, cex=0.5)
		dev.off()			
		#	get root to tip divergence
		anc.seq.nj		<- drop.tip(anc.seq.nj,'HXB2')
		tmp				<- node.depth.edgelength(anc.seq.nj)
		anc.seq.nj.info	<- data.table(LABEL=anc.seq.nj$tip.label, ROOT2TIP=tmp[seq_len(Ntip(anc.seq.nj))] )
		#tmp<- distRoot(anc.seq.nj)
		#set(anc.seq.nj.info, NULL, 'distRoot', tmp)
		set(anc.seq.nj.info, NULL, 'CALENDAR_TIME', anc.seq.nj.info[, as.numeric(sapply(strsplit(LABEL, tree.id.labelsep, fixed=TRUE),'[[',tree.id.label.idx.ctime))] )
		tmp				<- lm(ROOT2TIP~CALENDAR_TIME, data=anc.seq.nj.info)		 
		set( anc.seq.nj.info, NULL, 'ROOT2TIP_LM', predict(tmp, type='response') ) 	
		tmp2			<- c( R2=round(summary(tmp)$r.squared,d=3), SLOPE= as.numeric(round(coef(tmp)['CALENDAR_TIME'],d=4)), TMRCA=as.numeric(round( -coef(tmp)['(Intercept)']/coef(tmp)['CALENDAR_TIME'], d=1 )) )
		ggplot(anc.seq.nj.info, aes(x=CALENDAR_TIME, y=ROOT2TIP)) + geom_point(alpha=0.5) + geom_line(aes(y=ROOT2TIP_LM)) +
				#scale_x_continuous(breaks=seq(1980,2020,2)) +						
				labs(x='Sequence sampling date', y='root-to-tip divergence') +
				annotate("text", x=anc.seq.nj.info[, min(CALENDAR_TIME)], y=anc.seq.nj.info[, 0.9*max(ROOT2TIP)], label=paste("R2=", tmp2['R2'],'\nSlope=',tmp2['SLOPE'],'\nTMRCA=',tmp2['TMRCA'], sep=''), hjust = 0, size = 4) +
				theme(legend.position=c(0,1), legend.justification=c(0,1))		
		file			<- paste( outdir, '/', substr(outfile.fg,1,nchar(outfile.fg)-2),'_checkdraw', check.draw, '_Root2Tip.pdf', sep='' )
		ggsave(file=file, w=10, h=6)		
		#	save to file including outgroup
		file						<- paste( outdir, '/', substr(outfile.fg,1,nchar(outfile.fg)-2),'_checkdraw', check.draw,'_',insignat,'.R', sep='' )
		cat(paste('\nsave to file', file))
		save(anc.seq.draw, file=file)
	}
}
##--------------------------------------------------------------------------------------------------------
##	check ancestral sequences from BEAST XML, run ExaML
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.SIMU.SSAfg.checkancestralseq.runExaML<- function()
{
	#DATA			<<- "/work/or105/Gates_2014"
	DATA				<<- '/Users/Oliver/duke/2014_Gates'
	dir.name			<- DATA  	
	indir				<- paste(dir.name,'methods_comparison_rootseqsim/140813',sep='/')
	#	ExaML bootstrap args
	bs.from		<- 0
	bs.to		<- 1
	bs.n		<- 100
	
	#	search for 'checkdraw' files
	infiles		<- list.files(indir)
	infiles		<- infiles[ sapply(infiles, function(x) grepl('.*checkdraw[0-9]+.*R$',x) ) ]	
	if(!length(infiles))	stop('cannot find files matching criteria')
	
	outdir		<- indir
	for(i in seq_along(infiles))
	{
		infile		<- infiles[i]
		infile		<- substr(infile, 1, nchar(infile)-2)
		insignat	<- regmatches(infile, regexpr('checkdraw[0-9]+_.*', infile))
		insignat	<- regmatches(insignat,regexpr('_.*',insignat))
		insignat	<- substr(insignat,2,nchar(insignat))
		infile		<- regmatches(infile, regexpr('.*checkdraw[0-9]+', infile))
		
		
		cmd			<- hivc.cmd.examl.bootstrap(indir, infile, insignat, insignat, bs.from=bs.from, bs.to=bs.to,bs.n=bs.n,outdir=outdir, resume=1, verbose=1)
		dummy		<- lapply(cmd, function(x)
				{				
					x		<- hivc.cmd.hpcwrapper(x, hpc.walltime=24, hpc.q= NA, hpc.mem="450mb", hpc.nproc=1)
					#x		<- hivc.cmd.hpcwrapper(x, hpc.walltime=24, hpc.q="pqeph", hpc.mem="3850mb", hpc.nproc=8)
					signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
					outfile	<- paste("exa",signat,sep='.')
					#cat(x)
					hivc.cmd.hpccaller(outdir, outfile, x)
					Sys.sleep(1)
				})
		stop()
	}
}
##--------------------------------------------------------------------------------------------------------
##	check ancestral sequences from BEAST XML, run ExaML
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.SIMU.SSAfg.checkancestralseq.evalExaML<- function()
{
	require(ape)
	require(phytools)
	require(ggplot2)
	#DATA			<<- "/work/or105/Gates_2014"
	DATA				<<- '/Users/Oliver/duke/2014_Gates'
	dir.name			<- DATA  	
	indir				<- paste(dir.name,'methods_comparison_rootseqsim/140813',sep='/')
	label.sep			<- '|'
	label.idx.ctime		<- 4
	clu.thresh.bs		<- 0.9
	clu.thresh.brl		<- 0.045	#0.035
	bs.n				<- 100
	#	search for 'checkdraw examl' files
	infiles				<- list.files(indir)
	infiles				<- infiles[ grepl('.*checkdraw[0-9]+_examl.*newick$',infiles)  ]	
	if(!length(infiles))	stop('cannot find files matching criteria')
	
	#	unexpected:		R2 differs largely from NJ tree to ExaML tree, especially for i=3
	for(i in seq_along(infiles))
	{
		#i		<- 3
		infile	<- infiles[i]
		file	<- paste(indir, infile, sep='/')
		ph		<- read.tree(file)
		tmp		<- which(ph$tip.label=="HXB2")
		ph		<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])
		ph		<- ladderize(ph)
		#
		file	<- paste( indir, '/', substr(infile, 1, nchar(infile)-7), '_Tree.pdf', sep='' )
		pdf(file=file, h=150, w=10)
		plot(ph, cex=0.7)
		dev.off()
		#
		#	check root to tip divergence
		#
		ph		<- drop.tip(ph,'HXB2')
		tmp		<- node.depth.edgelength(ph)
		ph.info	<- data.table(LABEL=ph$tip.label, ROOT2TIP=tmp[seq_len(Ntip(ph))] )
		set(ph.info, NULL, 'CALENDAR_TIME', ph.info[, as.numeric(sapply(strsplit(LABEL, label.sep, fixed=TRUE),'[[',label.idx.ctime))] )
		tmp		<- lm(ROOT2TIP~CALENDAR_TIME, data=ph.info)		 
		set( ph.info, NULL, 'ROOT2TIP_LM', predict(tmp, type='response') ) 	
		tmp2	<- c( R2=round(summary(tmp)$r.squared,d=3), SLOPE= as.numeric(round(coef(tmp)['CALENDAR_TIME'],d=4)), TMRCA=as.numeric(round( -coef(tmp)['(Intercept)']/coef(tmp)['CALENDAR_TIME'], d=1 )) )
		ggplot(ph.info, aes(x=CALENDAR_TIME, y=ROOT2TIP)) + geom_point(alpha=0.5) + geom_line(aes(y=ROOT2TIP_LM)) +
				#scale_x_continuous(breaks=seq(1980,2020,2)) +						
				labs(x='Sequence sampling date', y='root-to-tip divergence') +
				annotate("text", x=ph.info[, min(CALENDAR_TIME)], y=ph.info[, 0.9*max(ROOT2TIP)], label=paste("R2=", tmp2['R2'],'\nSlope=',tmp2['SLOPE'],'\nTMRCA=',tmp2['TMRCA'], sep=''), hjust = 0, size = 4) +
				theme(legend.position=c(0,1), legend.justification=c(0,1))
		file	<- paste( indir, '/', substr(infile, 1, nchar(infile)-7), '_Root2Tip.pdf', sep='' )
		ggsave(file=file, w=10, h=6)
		if(0)
		{
			#	check brl divergence
			#brl.tips	<- distTips(ph , method='patristic')
			#
			#	check clustering among root seqs - this should be minimal at cut-off 0.045
			#
			dist.brl		<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=hivc.clu.min.transmission.cascade )
			quantile(dist.brl,p=c(0.01,0.05,0.1,0.2,0.3,0.4,0.5))
			hist(dist.brl, breaks=30)
			ph.node.bs		<- as.numeric( ph$node.label )		
			ph.node.bs[is.na(ph.node.bs)]	<- 0
			ph.node.bs		<- ph.node.bs/bs.n
			ph$node.label	<- ph.node.bs	
			ph.clu			<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=clu.thresh.bs, thresh.brl=clu.thresh.brl, dist.brl=dist.brl, nodesupport=ph.node.bs, retval="all")
			#	plot clustering tree
			file	<- paste( indir, '/', substr(infile, 1, nchar(infile)-7), '_Clutree.pdf', sep='' )
			pdf(file=file, h=150, w=10)
			tmp		<- hivc.clu.plot(ph, ph.clu[["clu.mem"]], show.tip.label=TRUE, cex.edge.incluster=1.5)
			dev.off()	
		}		
	}
	#	consider all bootstrap trees	
	indirs	<- c(	'/Users/Oliver/duke/2014_Gates/methods_comparison_rootseqsim/140813/PANGEA_SSApolBwhRc-_140811_n390_AncSeq_checkdraw1_examlout_Mon_Aug_17_17:05:23_2014',
					'/Users/Oliver/duke/2014_Gates/methods_comparison_rootseqsim/140813/PANGEA_SSApolBwhRc-_140811_n390_AncSeq_checkdraw2_examlout_Mon_Aug_17_17:05:23_2014',
					'/Users/Oliver/duke/2014_Gates/methods_comparison_rootseqsim/140813/PANGEA_SSApolBwhRc-_140811_n390_AncSeq_checkdraw3_examlout_Mon_Aug_17_17:05:23_2014'	)
	df.r2t	<- lapply(seq_along(indirs), function(i)
		{
			indir				<- indirs[i]
			infiles				<- list.files(indir)
			infiles				<- infiles[ grepl('^ExaML_result.*finaltree*',infiles)  ]	
			if(!length(infiles))	stop('cannot find files matching criteria')
					
			df.r2t	<- lapply(seq_along(infiles), function(j)
					{					
						infile	<- infiles[j]
						cat(paste('\nprocess file', infile))
						file	<- paste(indir, infile, sep='/')
						ph		<- read.tree(file)
						tmp		<- which(ph$tip.label=="HXB2")
						ph		<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])
						ph		<- ladderize(ph)
						ph		<- drop.tip(ph,'HXB2')
						tmp		<- node.depth.edgelength(ph)
						ph.info	<- data.table(LABEL=ph$tip.label, ROOT2TIP=tmp[seq_len(Ntip(ph))] )
						set(ph.info, NULL, 'CALENDAR_TIME', ph.info[, as.numeric(sapply(strsplit(LABEL, label.sep, fixed=TRUE),'[[',label.idx.ctime))] )	
						set(ph.info, NULL, 'BS', as.numeric(substr(infile, nchar(infile)-2, nchar(infile))) )
						ph.info
					})
			df.r2t	<- do.call('rbind',df.r2t)
			set(df.r2t, NULL, 'FILE', i)
			df.r2t
		})	
	df.r2t	<- do.call('rbind',df.r2t)
	#	the EXAML root to tip divergence is highly variable!
	#ggplot( subset( df.r2t, LABEL=='PANGEA_SSAfgBwhRc-_140811_n390_pool1|STATE_20400000|250|1987.78' ), aes(x=ROOT2TIP)) + geom_histogram()
	#	get R2 across all bootstrap ExaML tree
	df.r2t.su<- df.r2t[,	{
					tmp		<- lm(ROOT2TIP ~ CALENDAR_TIME)		 					 	
					list( R2=round(summary(tmp)$r.squared,d=3) )				
				},by=c('BS','FILE')]
	set(df.r2t.su, NULL, 'FILE', df.r2t.su[, factor(FILE)])	
		
	ggplot( df.r2t.su, aes(x=R2, fill=FILE)) + geom_histogram(binwidth=0.05) + facet_grid(.~FILE) +
			labs(fill='draw of 1000\nancestral\nsequences', x='R2 of lm(ROOT2TIP ~ CALENDAR_TIME)')
	indir		<- paste(dir.name,'methods_comparison_rootseqsim/140813',sep='/')
	file		<- paste(indir, '/','PANGEA_SSApolBwhRc-_140811_n390_AncSeq_checkdraw1-3_examlout.pdf', sep='')
	ggsave(file=file, w=10, h=5)
}
##--------------------------------------------------------------------------------------------------------
##	get anecestral sequences from BEAST XML - developper version
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.BEAST.devel.getancestralseq.from.output<- function()
{
	#dir.name	<- "/work/or105/Gates_2014"
	dir.name	<- '/Users/Oliver/duke/2014_Gates'  
	
	if(1)	#	devel
	{		
		indir				<- paste(dir.name,'methods_comparison_rootseqsim/140730',sep='/')
		outdir				<- indir
		infile.xml			<- 'working.xml'
		infile.beastparsed	<- 'working.R'
		outfile				<- 'working_ancseq.R'
		
		#	load BEAST PARSER output
		file		<- paste(indir, '/',infile.beastparsed, sep='')
		load(file)	#	expect tree, node.stat		
		#	get original sequences
		file		<- paste(indir, '/',infile.xml, sep='')
		bxml		<- xmlTreeParse(file, useInternalNodes=TRUE, addFinalizer = TRUE)
		bseq		<- hivc.beast.get.sequences(bxml, verbose=1)	
		bseq		<- merge(bseq, data.table(ALIGNMENT_ID=paste('alignment',1:3,sep=''), GENE=c('env','gag','pol')), by='ALIGNMENT_ID')				
		#	compute gag pol env		
		tmp			<- PANGEA.RootSeqSim.get.ancestral.seq(tree, node.stat, bseq, tree.id.sep='_', tree.id.idx.mcmcit=2, tree.id.burnin=1, label.sep='|', label.idx.ctime=2)
		ancseq.gag	<- tmp$GAG
		ancseq.env	<- tmp$ENV
		ancseq.pol	<- tmp$POL
		#	save as R
		file		<- paste(outdir, outfile, sep='/')
		save(ancseq.gag, ancseq.env, ancseq.pol, file=file)			
		#	save as FASTA
		file		<- paste(outdir, paste(substr(outfile,1,nchar(outfile)-1),'fasta',sep=''), sep='/')		
		write.dna(cbind(ancseq.gag, ancseq.env, ancseq.pol), file, format = "fasta")		
		#
		#	sample ancestral sequences between 1980-2000 and reconstruct tree with RAxML
		#
		ancseq				<- cbind(ancseq.gag, ancseq.env, ancseq.pol)
		label.sep			<- '|'
		label.idx.tree.id	<- 1
		label.idx.node.id	<- 2
		label.idx.ctime		<- 3
		ancseq.label		<- data.table(	TREE_ID= sapply(strsplit(rownames(ancseq), label.sep, fixed=1),'[[',label.idx.tree.id),
				NODE_ID= as.numeric(sapply(strsplit(rownames(ancseq), label.sep, fixed=1),'[[',label.idx.node.id)),
				CALENDAR_TIME= as.numeric(sapply(strsplit(rownames(ancseq), label.sep, fixed=1),'[[',label.idx.ctime)))
		hist( ancseq.label[, CALENDAR_TIME], breaks=100 )
	}
	if(0)
	{
		
		tmp	<- subset(node.stat, select=c(TREE_ID, NODE_ID, CALENDAR_TIME))
		setkey(tmp, TREE_ID, NODE_ID)
		hist(tmp[, CALENDAR_TIME], breaks=seq(1930,2011,1))
		#	reconstruct gene sequences and store in ape format
		#	get calendar time for gene sequence
		set(node.stat, NULL, 'VALUE', node.stat[, gsub('\"','',VALUE)])
		#	check seq lengths
		tmp		<- node.stat[, list(NCHAR=nchar(VALUE)), by=c('STAT','NODE_ID','TREE_ID')]
		stopifnot( tmp[, list(CHECK=all(NCHAR[1]==NCHAR)), by='STAT'][, all(CHECK)] )
		
		tmp[, list(CHECK=unique(NCHAR)), by='STAT']
		
		ENV.CP1<- "AGA"
		ENV.CP2<- "XYZ"
		ENV.CP3<- "KLM"
		tmp		<- do.call('rbind',sapply(list(ENV.CP1,ENV.CP2,ENV.CP3), strsplit, ''))
		tmp		<- paste(as.vector(tmp), collapse='')
		
		subset(node.stat, TREE_ID=='STATE_0' & NODE_ID==which(btree[[1]]$tip.label=='C.BW.AF443074|1996'))
		
	}
	if(0)	#devel
	{
		dir.name			<- '/Users/Oliver/duke/2014_Gates'  
		indir				<- paste(dir.name,'methods_comparison_rootseqsim/140801',sep='/')
		infile				<- 'ALLv06.n97.rlx.gmrf' 		
		insignat			<- 'Sun_Jul_27_09-00-00_2014'	
		file				<- paste(indir, '/', infile, '_', insignat, '.timetrees', sep='')
		
		indir				<- paste(dir.name,'methods_comparison_rootseqsim/140730',sep='/')		
		infile.timetrees	<- 'working.timetrees'
		outdir				<- indir
		
		file				<- paste(indir, '/', infile.timetrees, sep='')
		tmp					<- hivc.beast2out.read.nexus.and.stats(file, tree.id=NA, method.node.stat='any.node')
		tree				<- tmp$tree
		node.stat			<- tmp$node.stat
		
		file				<- paste(indir,'/',paste(substr(infile.timetrees,1,nchar(infile.timetrees)-9),'R',sep=''),sep='')
		save(tree, node.stat, file=file)
	}	
}
##--------------------------------------------------------------------------------------------------------
##	prepare LANL download
##	downloaded 556 full genome sequences from SSA on 06-08-2014
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.DATA.SSAfg.process.LosAlamos<- function()
{
	DATA			<<- '/Users/Oliver/duke/2014_Gates'
	indir			<- paste(DATA,'methods_comparison_rootseqsim/140806',sep='/')
	infile.fasta	<- 'PANGEA_SSAfg_140806_n557.fasta'
	infile.fasta.gag<- 'PANGEA_SSAfg_gag_140806_n557.fasta'
	infile.fasta.pol<- 'PANGEA_SSAfg_pol_140806_n557.fasta'
	infile.fasta.env<- 'PANGEA_SSAfg_env_140806_n557.fasta'
	outfile.outgroup<- 'PANGEA_SSAfg_140806_HXB2outgroup.R'
	file	<- paste(indir, '/', infile.fasta,sep='')
	seq		<- read.dna(file, format='fasta')	
	
	label.sep				<- '|'
	label.idx.country.id	<- 2
	label.idx.label			<- 3
	label.idx.ctime			<- 4
	tmp						<- strsplit(names(seq), label.sep, fixed=1)
	
	seq.label				<- data.table(	LABEL= names(seq),
			COUNTRY_ID= sapply(tmp,'[[',label.idx.country.id),
			NAME= sapply(strsplit(names(seq), label.sep, fixed=1),'[[',label.idx.label),
			CALENDAR_TIME= sapply(strsplit(names(seq), label.sep, fixed=1),'[[',label.idx.ctime))
	#	remove sequences without a sampling date								
	seq.label	<- subset( seq.label, !is.na(as.numeric(CALENDAR_TIME)) )							
	set(seq.label, NULL, 'CALENDAR_TIME', seq.label[, as.numeric(CALENDAR_TIME)])
	#	histogram of sample dates and sample location	
	ggplot(seq.label, aes(x=CALENDAR_TIME)) + geom_histogram(binwidth=1)
	ggplot(seq.label, aes(x=COUNTRY_ID)) + geom_histogram()
	#
	#	some seqs are crap and have duplicate runs, need to edit manually. ISOLATING GAG POL ENV
	#
	data(refseq_hiv1_hxb2)
	seq.gag.c	<- as.character(read.dna(file='~/git/hivclust/pkg/data/LosAlamos_HIV1B_CONSENSUS_2004_gag_DNA.fasta', format='fasta'))
	seq.pol.c	<- as.character(read.dna(file='~/git/hivclust/pkg/data/LosAlamos_HIV1B_CONSENSUS_2004_pol_DNA.fasta', format='fasta'))
	seq.env.c	<- as.character(read.dna(file='~/git/hivclust/pkg/data/LosAlamos_HIV1B_CONSENSUS_2004_env_DNA.fasta', format='fasta'))
	#	located gag start at 1009 manually, rm everything before
	seq			<- as.character(seq)
	for(i in seq_along(seq))	seq[[i]][1:1008]<- "-"
	tmp		<- c(1041:1043, 1218:1220, 1321:1322, 1326, 1341:1352, 1366:1389, 1395:1397, 1402:1407, 1411:1413, 1417:1440, 1443:1445, 1447:1449, 1896:1898, 2271:2291, 2455:2466, 2507:2509, 
			2511:2513, 2521:2523, 2551:2577, 2596, 2640, 2644:2655, 2668:2670 )	
	for(i in seq_along(seq))	seq[[i]][tmp]	<- '-'	
	#	add '-' to get all seqs of same length					
	tmp			<- max(sapply(seq, length))
	seq			<- lapply(seq, function(x) c(x, rep('-',tmp-length(x))) 	)
	#	get into matrix DNAbin	
	seq			<- as.DNAbin(do.call('rbind', seq))
	#	rm '-' columns so far
	seq			<- seq.rmallchar(seq, rm.char='-', verbose=1)
	file		<- paste(indir, '/fixup1_', infile.fasta,sep='')
	write.dna(seq, file=file, format='fasta')
	#	located pol/prot start at 1657 manually, cut and rm gaps and re-align gag
	seq				<- as.character(seq)
	seq.gag			<- seq[, 1:1656]
	tmp				<- rownames(seq.gag)
	seq.gag			<- lapply(seq_len(nrow(seq.gag)), function(i){	seq.gag[i, seq.gag[i,]!="-" ]	})
	seq.gag[[length(seq.gag)+1]]	<- hxb2[, as.character(HXB2.K03455)]
	seq.gag[[length(seq.gag)+1]]	<- seq.gag.c['CONSENSUS_C',]
	names(seq.gag)	<- c(tmp,'HXB2','CONSENSUS_C')	
	seq.gag			<- as.DNAbin(seq.gag)
	file			<- paste(indir, '/', infile.fasta.gag,sep='')
	write.dna(seq.gag, file=file, format='fasta')
	#	align GAG
	tmp				<- hivc.cmd.clustalo(indir, infile.fasta.gag, signat='', outdir=indir, prog= PR.CLUSTALO, hmm=PR.CLUSTALO.HMM, nproc=1, verbose=1)
	#	continue fixup for POL	
	file			<- paste(indir, '/fixup1_', infile.fasta,sep='')
	seq				<- read.dna(file, format='fasta')
	seq				<- as.character(seq)
	seq[, 1:1656]	<- "-"
	#	located end of POL at 4597	- found length issues, align with own HXB2 
	seq.pol			<- seq[, 1:4700]	
	tmp				<- rownames(seq.pol)
	seq.pol			<- lapply(seq_len(nrow(seq.pol)), function(i){	seq.pol[i, seq.pol[i,]!="-" ]	})	
	seq.pol[[length(seq.pol)+1]]	<- hxb2[, as.character(HXB2.K03455)]
	seq.pol[[length(seq.pol)+1]]	<- seq.pol.c['CONSENSUS_C',]	
	names(seq.pol)	<- c(tmp,'HXB2','CONSENSUS_C')	
	seq.pol			<- as.DNAbin( seq.pol )
	file			<- paste(indir, '/', infile.fasta.pol,sep='')
	write.dna(seq.pol, file=file, format='fasta')
	#	align POL
	#	/Users/Oliver/git/hivclust/pkg/inst/mafft-mac/mafft.bat --op 1.8 --ep 0.4 --maxiterate 15 --thread 4 /Users/Oliver/duke/2014_Gates/methods_comparison_rootseqsim/140806/PANGEA_SSAfg_pol_140806_n557.fasta > /Users/Oliver/duke/2014_Gates/methods_comparison_rootseqsim/140806/PANGEA_SSAfg_pol_140806_n557.fasta.mafft
	#tmp				<- hivc.cmd.clustalo(indir, infile.fasta.pol, signat='', outdir=indir, prog= PR.CLUSTALO, hmm=PR.CLUSTALO.HMM, nproc=1, verbose=1)
	#
	#	located start of ENV around 5918			which( seq[nrow(seq), 4597:ncol(seq)]=='-' )		which( seq[nrow(seq), 5918:ncol(seq)]=='-' )	-> 660 -> 9058
	#	located end of ENV around 9176 
	file			<- paste(indir, '/fixup1_', infile.fasta,sep='')
	seq				<- read.dna(file, format='fasta')
	seq				<- as.character(seq)
	seq.env			<- seq[, 5850:9250]
	tmp				<- rownames(seq.env)
	seq.env			<- lapply(seq_len(nrow(seq.env)), function(i){	seq.env[i, seq.env[i,]!="-" ]	})	
	seq.env[[length(seq.env)+1]]	<- hxb2[, as.character(HXB2.K03455)]
	seq.env[[length(seq.env)+1]]	<- seq.env.c['CONSENSUS_C',]	
	names(seq.env)	<- c(tmp,'HXB2','CONSENSUS_C')
	seq.env			<- as.DNAbin( seq.env )
	file			<- paste(indir, '/', infile.fasta.env,sep='')
	write.dna(seq.env, file=file, format='fasta')
	#	align ENV
	tmp			<- hivc.cmd.clustalo(indir, infile.fasta.env, signat='', outdir=indir, prog= PR.CLUSTALO, hmm=PR.CLUSTALO.HMM, nproc=1, verbose=1)
	#
	#	process GAG
	#	
	file		<- paste(indir, '/', infile.fasta.gag,'.mafft',sep='')
	seq.gag		<- read.dna(file, format='fasta')
	seq.gag		<- as.character(seq.gag)
	seq.gag		<- seq.gag[,790:2698]
	tmp			<- rownames(seq.gag)
	seq.gag		<- lapply(seq_len(nrow(seq.gag)), function(i){	seq.gag[i, seq.gag[i,]!="-" ]	})
	names(seq.gag)	<- tmp
	seq.gag			<- as.DNAbin(seq.gag)
	infile.fasta.gag<- 'PANGEA_SSAfg_gag2_140806_n557.fasta'
	file			<- paste(indir, '/', infile.fasta.gag,sep='')
	write.dna(seq.gag, file=file, format='fasta')	
	#/Users/Oliver/git/hivclust/pkg/inst/mafft-mac/mafft.bat --op 3.0 --ep 1.5 --maxiterate 1000 --thread 4 /Users/Oliver/duke/2014_Gates/methods_comparison_rootseqsim/140806/PANGEA_SSAfg_gag2_140806_n557.fasta > /Users/Oliver/duke/2014_Gates/methods_comparison_rootseqsim/140806/PANGEA_SSAfg_gag2_140806_n557.fasta.mafft
	#	manual edits
	file			<- paste(indir, '/', infile.fasta.gag,'.mafft',sep='')
	seq.gag			<- read.dna(file, format='fasta')
	seq.gag			<- as.character(seq.gag)
	seq.gag[, c(1361:1404, 1406:1413, 1421:1458)]	<- '-'
	seq.gag			<- seq.rmgaps(seq.gag, rm.only.col.gaps=1, verbose=1)
	infile.fasta.gag<- 'PANGEA_SSAfg_gag3_140806_n557.fasta'
	file			<- paste(indir, '/', infile.fasta.gag,sep='')
	write.dna(seq.gag, file=file, format='fasta')
	#	final without HXB2 / CONSENSUS	
	seq.gag			<- seq.rmgaps(seq.gag, rm.only.col.gaps=1, verbose=1)	#len 1466 -- 3 more than expected; there is an AA insertion for the SA seqs at pos 367
	outgroup.seq.gag<- seq.gag[ 'HXB2'==rownames(seq.gag),  ]
	seq.gag			<- seq.gag[ -c( which(grepl('HXB2',rownames(seq.gag))), which(grepl('CONSENSUS',rownames(seq.gag))) ),  ]
	seq.gag			<- seq.rmgaps(seq.gag, rm.only.col.gaps=1, verbose=2)
	#	based on verbose output, cut the following sites from outgroup.seq.gag 
	#	373:381, 1123:1125, 1401:1403, 1408:1410, 1449:1454
	outgroup.seq.gag<- outgroup.seq.gag[, -c(373:381, 1123:1125, 1401:1403, 1408:1410, 1449:1454)]
	infile.fasta.gag<- 'PANGEA_SSAfg_gag_140806_n556_final.fasta'
	file			<- paste(indir, '/', infile.fasta.gag,sep='')
	write.dna(seq.gag, file=file, format='fasta')	
	#
	#	process POL
	#	
	file			<- paste(indir, '/', infile.fasta.pol,'.mafft',sep='')
	seq.pol			<- read.dna(file, format='fasta')
	#	start: 2253		end: 5096	
	seq.pol			<- seq.pol[,2253:5096]
	outgroup.seq.pol<- seq.pol[ 'HXB2'==rownames(seq.pol),  ]
	seq.pol			<- seq.pol[ -c( which(grepl('HXB2',rownames(seq.pol))), which(grepl('CONSENSUS',rownames(seq.pol))) ),  ]
	infile.fasta.pol<- 'PANGEA_SSAfg_pol_140806_n556_final.fasta'
	file			<- paste(indir, '/', infile.fasta.pol,sep='')
	write.dna(seq.pol, file=file, format='fasta')
	#
	#	process ENV
	#
	#	start 6216	end 9734
	file		<- paste(indir, '/', infile.fasta.env,'.mafft',sep='')
	seq.env		<- read.dna(file, format='fasta')
	seq.env		<- as.character(seq.env)
	tmp			<- names(seq.env)
	seq.env		<- t(sapply( seq_along(seq.env), function(i){	seq.env[[i]][ 6216:9734 ]	}))
	seq.env		<- lapply(seq_len(nrow(seq.env)), function(i){	seq.env[i, seq.env[i,]!="-" ]	})	
	names(seq.env)	<- tmp
	seq.env			<- as.DNAbin(seq.env)
	infile.fasta.env<- 'PANGEA_SSAfg_env2_140806_n557.fasta'
	file			<- paste(indir, '/', infile.fasta.env,sep='')
	write.dna(seq.env, file=file, format='fasta')	
	#	align on LosAlamos + manual edits
	file			<- paste(indir, '/', infile.fasta.env,'.hivalign',sep='')
	seq.env			<- read.dna(file, format='fasta')	
	seq.env			<- as.character(seq.env)
	seq.env[,c(1613:1645,2852:2932) ]	<- '-'  
	seq.env			<- seq.rmgaps(seq.env, rm.only.col.gaps=1, verbose=1)	
	infile.fasta.env<- 'PANGEA_SSAfg_env3_140806_n557.fasta'
	file			<- paste(indir, '/', infile.fasta.env,sep='')
	write.dna(seq.env, file=file, format='fasta', colsep='', nbcol=-1)
	#	align on LosAlamos + manual edits
	file			<- paste(indir, '/', infile.fasta.env,'.hivalignnt',sep='')
	seq.env			<- read.dna(file, format='fasta')	
	seq.env			<- as.character(seq.env)
	tmp				<- max(sapply(seq.env, length))
	seq.env			<- t(sapply(seq.env, function(x) c(x, rep('-',tmp-length(x))) 	))
	seq.env[,c(430:543) ]	<- '-'		#rm non-consensus jitter
	seq.env[,c(1294:1344) ]	<- '-'		#rm non-consensus jitter
	seq.env			<- seq.rmgaps(seq.env, rm.only.col.gaps=1, verbose=1)
	infile.fasta.env<- 'PANGEA_SSAfg_env4_140806_n557.fasta'
	file			<- paste(indir, '/', infile.fasta.env,sep='')
	write.dna(seq.env, file=file, format='fasta', colsep='', nbcol=-1)
	#seq.env[,c(1:9) ]	<- '-'			#C consensus start
	#!! There are two extra AA that are not in consensus C TAGTAG at pos 561
	#!! There are two extra AA that are not in consensus C at pos 2310
	#fixed some final issues manually
	seq.env			<- read.dna(file, format='fasta')
	outgroup.seq.env<- seq.env[ 'HXB2'==rownames(seq.env),  ]
	seq.env 		<- seq.env[ -c( which(grepl('HXB2',rownames(seq.env))), which(grepl('CONSENSUS',rownames(seq.env))) ),  ]
	seq.env			<- seq.rmgaps(seq.env, rm.only.col.gaps=1, verbose=2)
	#	based on verbose output, cut the following sites from outgroup.seq.env 
	#	55:66, 919:924, 1057:1059
	outgroup.seq.env<- outgroup.seq.env[, -c(55:66, 919:924, 1057:1059)]
	infile.fasta.env<- 'PANGEA_SSAfg_env_140806_n556_final.fasta'
	file			<- paste(indir, '/', infile.fasta.env,sep='')
	write.dna(seq.env, file=file, format='fasta', colsep='', nbcol=-1)
	#	write outgroups
	file			<- paste(indir, '/', outfile.outgroup, sep='')
	save(file=file, outgroup.seq.gag, outgroup.seq.pol, outgroup.seq.env)
}
##--------------------------------------------------------------------------------------------------------
##	get sequences from BEAST XML, create concatenated file
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.DATA.checkDrugResistanceMask<- function()
{
	require(XML)
	require(ape)
	DATA		<<- "/work/or105/Gates_2014"
	DATA		<<- '/Users/Oliver/duke/2014_Gates'
	indir		<- paste(DATA,'methods_comparison_rootseqsim/140727',sep='/')
	outdir		<- paste(DATA,'methods_comparison_rootseqsim/140728',sep='/')
	infile		<- 'ALLv02.n100.rlx.gmrf' 
	insignat	<- 'Sun_Jul_27_09-00-00_2014'
	
	file			<- '/Users/Oliver/duke/2014_Gates/methods_comparison_rootseqsim/140727/ALLv01.n100.rlx.gmrf_Sun_Jul_27_09-00-00_2014.xml'
	file			<- '/Users/Oliver/duke/2014_Gates/methods_comparison_rootseqsim/140801/ALLv06.n97.rlx.gmrf_Sun_Jul_27_09-00-00_2014.xml'
	bxml			<- xmlTreeParse(file, useInternalNodes=TRUE, addFinalizer = TRUE)
	bseq			<- hivc.beast.get.sequences(bxml, verbose=1)	
	bseq			<- merge(bseq, data.table(ALIGNMENT_ID=paste('alignment',1:3,sep=''), GENE=c('env','gag','pol')), by='ALIGNMENT_ID')
	#	check all seq of same length per gene
	bseq			<- merge(bseq, bseq[, list(SEQ_N=nchar(SEQ)), by=c('GENE','TAXON_ID')], by=c('GENE','TAXON_ID'))
	stopifnot( bseq[, all(SEQ_N==SEQ_N[1]), by='GENE'][, all(V1)] )
	#	check 3 genes per taxon
	stopifnot( bseq[, length(SEQ)==3, by='TAXON_ID'][, all(V1)] )
	#	check if indeed patterns are compressed
	if(0)
	{
		bseq	<- subset(bseq, GENE=='env')
		tmp		<- tolower(do.call('rbind',strsplit(bseq[, SEQ],'')))
		bseq.CP1<- tmp[, seq.int(1, ncol(tmp), by=3)]
		bseq.CP2<- tmp[, seq.int(2, ncol(tmp), by=3)]
		bseq.CP3<- tmp[, seq.int(3, ncol(tmp), by=3)]
		
		tmp		<- bseq.CP3
		tmp2	<- apply( tmp, 1, function(x) paste(x,sep='',collapse=''))	#identical sequences?	
		cat(paste('\nunique sequences, n=',length(unique(tmp2))))		
		tmp2	<- apply( tmp, 2, function(x) paste(x,sep='',collapse=''))	#identical patterns? 
		cat(paste('\nunique patterns, n=',length(unique(tmp2))))
		tmp2	<- apply( tmp, 2, function(x) all(x==x[1]))					#invariant sites?
		length(which(tmp2))
	}
	#	check for drug resistance mutations
	if(0)
	{
		tmp				<- subset(bseq, GENE=='pol')
		bseq.pol.m		<- do.call('rbind',strsplit(tmp[, SEQ],''))
		
		file			<- '/Users/Oliver/git/hivclust/pkg/data/IAS_primarydrugresistance_201303.rda'
		alignment.start	<- 2085
		load(file)
		IAS_primarydrugresistance_201303		<- as.data.table(IAS_primarydrugresistance_201303)
		set(IAS_primarydrugresistance_201303, NULL, "Alignment.nuc.pos", IAS_primarydrugresistance_201303[,Alignment.nuc.pos]-alignment.start+1)
		#pol not in HXB2 consensus coordinates - TODO would have to align against consensus
		z<- seq.rm.drugresistance(bseq.pol.m, IAS_primarydrugresistance_201303, verbose=1, rtn.DNAbin=0)		
	}
	#	concatenate into single DNAbin matrix and save
	tmp		<- dcast.data.table(bseq, TAXON_ID ~ GENE, value.var="SEQ")	
	tmp[, SEQ_ALL:=paste(gag, pol, env, sep='')]
	tmp2	<- tolower(do.call('rbind',strsplit(tmp[, SEQ_ALL],'')))
	rownames(tmp2)	<- tmp[, TAXON_ID]
	seq		<- as.DNAbin(tmp2)		
	outfile	<- paste(infile,'_conc_',insignat,'.R',sep='')
	save(seq, file=paste(outdir, outfile, sep='/'))
}
##--------------------------------------------------------------------------------------------------------
##	check for recombinants
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.DATA.checkRecombinants<- function()
{
	require(XML)
	require(ape)
	require(r3SEQ)
	DATA			<<- "/work/or105/Gates_2014"
	DATA			<<- '/Users/Oliver/duke/2014_Gates'
	indir			<- paste(DATA,'methods_comparison_rootseqsim/140728',sep='/')	
	infile			<- 'ALLv02.n100.rlx.gmrf_conc'
	outfile			<- 'ALLv03.n97.rlx.gmrf'
	insignat		<- 'Sun_Jul_27_09-00-00_2014'
	
	#file		<- paste(indir, '/', infile, '_', insignat, '.R', sep='')
	#load(file)
	#	run 3SEQ
	infile.3seq		<- paste(infile, '_', insignat, '.R', sep='')
	pipeline.recom.run.3seq(indir, infile.3seq, batch.n=1, hpc.walltime=1, hpc.q=NA, hpc.mem="500mb", hpc.nproc=1)
	#	parse 3SEQ output
	argv			<<-	cmd.recombination.process.3SEQ.output(indir, infile.3seq, '', resume=1, verbose=1) 
	argv			<<- unlist(strsplit(argv,' '))
	df.recomb		<- prog.recom.process.3SEQ.output()	
	#	select potential recombinants with p-value < 1e-3
	df.recomb		<- subset( df.recomb, adjp<1e-3 )
	cat(paste('\nfound potential recombinants, n=',nrow(df.recomb)))
	#
	#	remove potential recombinants from XML file
	#
	file			<- '/Users/Oliver/duke/2014_Gates/methods_comparison_rootseqsim/140727/ALLv02.n100.rlx.gmrf_Sun_Jul_27_09-00-00_2014.xml'
	bxml.template	<- xmlTreeParse(file, useInternalNodes=TRUE, addFinalizer = TRUE)
	bseq			<- hivc.beast.get.sequences(bxml.template, verbose=1)	
	bseq			<- merge(bseq, data.table(ALIGNMENT_ID=paste('alignment',1:3,sep=''), GENE=c('env','gag','pol')), by='ALIGNMENT_ID')
	bseq			<- subset(bseq, !TAXON_ID%in%df.recomb[, child])
	#
	#	write XML file with new sequences
	#		
	bxml			<- newXMLDoc(addFinalizer=T)
	bxml.beast		<- newXMLNode("beast", doc=bxml, addFinalizer=T)
	newXMLCommentNode(text=paste("Generated by HIVCLUST from template",file), parent=bxml.beast, doc=bxml, addFinalizer=T)
	#	add new set of ENV sequences into alignment ID 1	
	set( bseq, NULL, 'SEQ', bseq[, tolower(SEQ)] )
	tmp				<- subset(bseq, GENE=='env')
	tmp2			<- tmp[, do.call('rbind',strsplit(tmp[, SEQ],''))]
	rownames(tmp2)	<- tmp[, TAXON_ID]
	tmp2			<- as.DNAbin(tmp2)
	bxml			<- hivc.beast.add.seq(bxml, tmp2, df=NULL, beast.label.datepos= 2, beast.label.sep= '|', beast.date.direction= "forwards", beast.date.units= "years", beast.alignment.id="alignment1", beast.alignment.dataType= "nucleotide", verbose=1)
	#	add new set of GAG sequences into alignment ID 2
	tmp				<- subset(bseq, GENE=='gag')
	tmp2			<- tmp[, do.call('rbind',strsplit(tmp[, SEQ],''))]
	rownames(tmp2)	<- tmp[, TAXON_ID]
	tmp2			<- as.DNAbin(tmp2)
	bxml			<- hivc.beast.add.seq(bxml, tmp2, df=NULL, beast.label.datepos= 2, beast.label.sep= '|', beast.date.direction= "forwards", beast.date.units= "years", beast.alignment.id="alignment2", beast.alignment.dataType= "nucleotide", verbose=1)
	#	add new set of POL sequences into alignment ID 3
	tmp				<- subset(bseq, GENE=='pol')
	tmp2			<- tmp[, do.call('rbind',strsplit(tmp[, SEQ],''))]
	rownames(tmp2)	<- tmp[, TAXON_ID]
	tmp2			<- as.DNAbin(tmp2)
	bxml			<- hivc.beast.add.seq(bxml, tmp2, df=NULL, beast.label.datepos= 2, beast.label.sep= '|', beast.date.direction= "forwards", beast.date.units= "years", beast.alignment.id="alignment3", beast.alignment.dataType= "nucleotide", verbose=1)
	#	copy rest from template	
	bt.beast		<- getNodeSet(bxml.template, "//beast")[[1]]
	dummy			<- sapply(seq.int( max(which( xmlSApply(bt.beast, xmlName)=="alignment" ))+1, xmlSize(bt.beast) ), function(i)
			{
				if( class(bt.beast[[i]])[1]=="XMLInternalCommentNode" )
					dummy<- newXMLCommentNode(text=xmlValue(bt.beast[[i]]), parent=bxml.beast, doc=bxml, addFinalizer=T)
				else
					dummy<- addChildren( bxml.beast, xmlClone( bt.beast[[i]], addFinalizer=T, doc=bxml ) )
			})
	#	change outfile name 
	bxml.onodes	<- getNodeSet(bxml, "//*[@fileName]")
	tmp			<- sapply(bxml.onodes, function(x) xmlGetAttr(x,"fileName"))
	tmp			<- gsub("(time).","time",tmp,fixed=1)
	tmp			<- gsub("(subst).","subst",tmp,fixed=1)	
	tmp			<- sapply(strsplit(tmp,'.',fixed=1), function(x)	paste(outfile,'_',insignat, '.', tail(x,1), sep=''))		
	dummy		<- sapply(seq_along(bxml.onodes), function(i){		xmlAttrs(bxml.onodes[[i]])["fileName"]<- tmp[i]		})
	#	write to file
	file		<- paste(indir,'/',outfile,'_',insignat,".xml", sep='')
	if(verbose)	cat(paste("\nwrite xml file to",file))
	saveXML(bxml, file=file)

	if(0)
	{
		#	get RAXML trees
		triplets			<- 1
		#triplets			<- 147:nrow(df.recomb)
		dummy	<- lapply(triplets, function(i)
				{				
					i<- 1
					if(verbose)	cat(paste("\nprocess triplet number",i,"\n"))
					argv				<<- cmd.recombination.check.candidates(indir, infile.3seq, '', i, resume=resume, verbose=1, hpc.walltime=1, hpc.q=NA, hpc.mem='500mb', hpc.nproc=1)
					argv				<<- unlist(strsplit(argv,' '))
					prog.recom.get.incongruence()		#this starts ExaML for the ith triplet			
				})	
	}		
}
##--------------------------------------------------------------------------------------------------------
##	simulate sequence sampling and break up transmission chains by imports
##--------------------------------------------------------------------------------------------------------
project.PANGEA.SampleSeq.simulate<- function()
{
	indir			<- system.file(package="rPANGEAHIVsim", "misc")	
	tmpdir.HPTN071	<- '/Users/Oliver/git/HPTN071sim/tmp140908'
	infile.ind		<- '140716_RUN001_IND.csv'
	infile.trm		<- '140716_RUN001_TRM.csv'	
	outfile.ind		<- '140716_RUN001_IND.csv'
	outfile.trm		<- '140716_RUN001_TRM.csv'
		
	dir.create(tmpdir.HPTN071, showWarnings=FALSE)
	#get input into 'argv'. this is needed because the input parser is usually called from the command line, and 'argv' mimics the way input is provided when the parser is called from the command line
	cmd				<- cmd.HPTN071.input.parser(indir, infile.trm, infile.ind, tmpdir.HPTN071,  infile.trm, infile.ind)				 
	argv			<<- unlist(strsplit(cmd,' '))
	prog.HPTN071.input.parser.v1()
	
}
##--------------------------------------------------------------------------------------------------------
##	pipeline to generate sequence data sets from HPTN071 output
##--------------------------------------------------------------------------------------------------------
pipeline.HPTN071<- function()
{
	stop()
	#	example input files
	indir			<- system.file(package="rPANGEAHIVsim", "misc")
	indir			<- ifelse(indir=='','/Users/Oliver/git/HPTN071sim/raw_trchain',indir)
	infile.ind		<- '140716_RUN001_IND.csv'
	infile.trm		<- '140716_RUN001_TRM.csv'	
	#	re-name the following:
	tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140908'
	#
	#	pipeline start
	#	
	##	sample sequences and draw imports 
	cmd				<- paste('mkdir -p ', tmpdir,'\n',sep='')	
	tmpdir.HPTN071	<- paste(tmpdir,'/HPTN071parser',sep='')
	cmd				<- paste(cmd, 'mkdir -p ', tmpdir.HPTN071,'\n',sep='')
	cmd				<- paste(cmd, cmd.HPTN071.input.parser.v2(indir, infile.trm, infile.ind, tmpdir.HPTN071,  infile.trm, infile.ind), sep='\n')
	##	run virus tree simulator
	tmpdir.VTS		<- paste(tmpdir,'/VirusTreeSimulator',sep='')
	cmd				<- paste(cmd, 'mkdir -p ', tmpdir.VTS,'\n',sep='')
	outfile			<- substr(infile.ind, 1, nchar(infile.ind)-7)
	prog.args		<- '-demoModel Logistic -N0 100000 -growthRate 0.0001 -t50 -0.04'
	cmd				<- paste(cmd, cmd.VirusTreeSimulator(tmpdir.HPTN071, infile.trm, infile.ind, tmpdir.VTS, outfile, prog.args=prog.args), sep='\n')	
	##	create seq gen input files 
	tmpdir.SG		<- paste(tmpdir,'/SeqGen',sep='')
	cmd				<- paste(cmd, 'mkdir -p ', tmpdir.SG,'\n',sep='')
	infile.epi		<- paste( substr(infile.ind, 1, nchar(infile.ind)-7),'SAVE.R', sep='' )
	infile.vts		<- substr(infile.ind, 1, nchar(infile.ind)-7)
	cmd				<- paste(cmd, cmd.SeqGen.createInputFiles(indir, infile.epi, tmpdir.VTS, infile.vts, tmpdir.SG), sep='\n')
	
	stop()
	#	currently requires output from last step so this is a separate batch file
	cmd				<- ''
	infile.epi		<- paste( substr(infile.ind, 1, nchar(infile.ind)-7),'SAVE.R', sep='' )
	plot.file		<- paste(indir, '/', substr(infile.epi,1,nchar(infile.epi)-6),'INFO_sg_ERPOSTERIOR_by_MODEL.pdf', sep='')
	cmd				<- paste(cmd, cmd.PANGEA.SeqGen(tmpdir.SG, infile.vts, plot.file=plot.file), sep='\n')		
	tmpdir.sim		<- paste(tmpdir,'/Simu',sep='')
	cmd				<- paste(cmd, 'mkdir -p ', tmpdir.sim,'\n',sep='')		
	cmd				<- paste(cmd, cmd.SeqGen.readOutputFiles(indir, infile.epi, tmpdir.SG, tmpdir.sim), sep='')
	cat(cmd)
	
	outfile			<- paste("pngea",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
	cmd.hpccaller(indir, outfile, cmd)
	quit("no")
}
##--------------------------------------------------------------------------------------------------------
##	devel code to process SeqGen output
##--------------------------------------------------------------------------------------------------------
project.PANGEA.SeqGen.readOutput<- function(indir.sg, infile.prefix)
{
	label.idx.codonpos	<- 1
	label.idx.gene		<- 2
	label.idx.clu		<- 3
	treelabel.idx.idpop	<- 1
	treelabel.idx.sep	<- '|'

	indir.epi		<- '/Users/Oliver/git/HPTN071sim/data_HPTN071epimodel_output'
	infile.epi		<- '140716_RUN001_SAVE.R'	
	indir.sg		<- '/Users/Oliver/git/HPTN071sim/tmp/SeqGen'
	outdir			<- '/Users/Oliver/git/HPTN071sim/data_HPTN071epimodel_output'
	#	load simulated epi data
	file		<- paste(indir.epi, '/', infile.epi, sep='')
	load(file)	#expect "df.epi"    "df.trms"   "df.inds"   "df.sample"
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
	#	merge simulated data
	simu.df		<- merge(subset(df.seq, select=c(IDPOP, GAG, POL, ENV)), subset( df.inds, !is.na(TIME_SEQ) ), by='IDPOP', all.x=TRUE)
	#
	#	save simulated data -- internal
	#
	outfile.prefix	<- substr(infile.epi,1,nchar(infile.epi)-6)
	file			<- paste(outdir, '/', outfile.prefix, 'SIMULATED_INTERNAL.R', sep='')
	cat(paste('\nwrite to file', file))
	save(df.epi, df.trms, df.inds, df.sample, df.seq, file=file)
	#
	#	save simulated data -- to be shared
	#	
	tmp				<- subset( df.inds, !is.na(TIME_SEQ), select=c(IDPOP, GENDER, CIRCM, DOB, DOD, TIME_SEQ ) )
	file			<- paste(outdir, '/', outfile.prefix, 'SIMULATED_metadata.csv', sep='')
	cat(paste('\nwrite to file', file))
	write.csv(tmp, file)		
	tmp				<- tolower(do.call('rbind',strsplit(df.seq[, GAG],'')))
	rownames(tmp)	<- df.seq[, LABEL]
	tmp				<- as.DNAbin(tmp)
	file			<- paste(outdir, '/', outfile.prefix, 'SIMULATED_gag.fa', sep='')
	write.dna(tmp, file, format = "fasta")	
	tmp				<- tolower(do.call('rbind',strsplit(df.seq[, POL],'')))
	rownames(tmp)	<- df.seq[, LABEL]
	tmp				<- as.DNAbin(tmp)
	file			<- paste(outdir, '/', outfile.prefix, 'SIMULATED_pol.fa', sep='')
	write.dna(tmp, file, format = "fasta")	
	tmp				<- tolower(do.call('rbind',strsplit(df.seq[, ENV],'')))
	rownames(tmp)	<- df.seq[, LABEL]
	tmp				<- as.DNAbin(tmp)
	file			<- paste(outdir, '/', outfile.prefix, 'SIMULATED_env.fa', sep='')
	write.dna(tmp, file, format = "fasta")	
	#	zip simulated files
	tmp				<- c( paste(outdir, '/', outfile.prefix, 'SIMULATED_metadata.csv', sep=''), paste(outdir, '/', outfile.prefix, 'SIMULATED_env.fa', sep=''), paste(outdir, '/', outfile.prefix, 'SIMULATED_pol.fa', sep=''), paste(outdir, '/', outfile.prefix, 'SIMULATED_gag.fa', sep='') )
	zip( paste(outdir, '/', outfile.prefix, 'SIMULATED.zip', sep=''), tmp, flags = "-FSr9XTj")
	file.remove(tmp)
}
##--------------------------------------------------------------------------------------------------------
##	devel code to call VirusTreeSimulator
##--------------------------------------------------------------------------------------------------------
project.PANGEA.VirusTreeSimulator.v1<- function()	
{
	require(data.table)
	indir		<- '/Users/Oliver/git/HPTN071sim/tmp/HPTN071parser'
	outdir		<- '/Users/Oliver/duke/2014_Gates/methods_comparison_trchphylosim/140819'
	infile.ind	<- '140716_RUN001_IND.csv'
	infile.trm	<- '140716_RUN001_TRM.csv'
	outfile		<- '140716_RUN001_VirusTreeSim_'
			
	#	N0			effective pop size at time 0
	#	growthRate 	the effective population size growth rate
	#	t50 		the time point, relative to the time of infection in backwards time, at which the population is equal to half its final asymptotic value
	cmd			<- cmd.VirusTreeSimulator(indir, infile.trm, infile.ind, outdir, outfile, prog.args='-demoModel Logistic -N0 0.1 -growthRate 1.5 -t50 -4')
	cmd			<- cmd.VirusTreeSimulator(indir, infile.trm, infile.ind, outdir, outfile, prog.args='-demoModel Logistic -N0 100000 -growthRate 0.0001 -t50 -0.04')
			
	outdir		<- '/Users/Oliver/duke/2014_Gates/methods_comparison_trchphylosim/140827'
	cmd			<- cmd.VirusTreeSimulator(indir, infile.trm, infile.ind, outdir, outfile, prog.args='-demoModel Constant -N0 15000')
	cat(cmd)
}


