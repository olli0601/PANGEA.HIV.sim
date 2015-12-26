


prog.hello<- function()	
{
	print('hello')
	
	file	<- '/Users/Oliver/Library/R/2.15/library/rPANGEAHIVsim/libs/x86_64/rPANGEAHIVsim.so'
	dyn.load(file)
	is.loaded('seqgen')
	argv	<- "-n1 -k1 -on -z42 -mGTR -a1 -g4 -i0 -s1 -f0.25,0.25,0.25,0.25 -r1,1,1,1,1,1 </Users/Oliver/git/HPTN071sim/tmp140914-3/SeqGen/140716_RUN001-3_911_POL_CP3.phy> /Users/Oliver/git/HPTN071sim/tmp140914-3/SeqGen/140716_RUN001-3_911_POL_CP3.seqgen"
	z		<- strsplit(argv, ' ')[[1]]
	z2		<- as.integer(length(z))
	.C('seqgen', z2, z)
	
	indir	<- '/Users/Oliver/git/HPTN071sim/tmp140914-3/SeqGen'
	infile	<- '140716_RUN001-3_911_POL_CP3.phy'
	outdir	<- indir
	outfile	<- '140716_RUN001-3_911_POL_CP3.seqgen'
	cmd.SeqGen(indir, infile, outdir, outfile,  prog.args='-n1 -k1 -on -z42', alpha=1, gamma=4, invariable=0, scale=1, 
			freq.A=0.25, freq.C=0.25, freq.G=0.25, freq.T=0.25,
			rate.AC=1, rate.AG=1, rate.AT=1, rate.CG=1, rate.CT=1, rate.GT=1)

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
	#DATA		<<- '~/git/HPTN071sim/data_rootseq'
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
	#DATA			<<- '/Users/Oliver/duke/2014_Gates'
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
	#DATA		<<- '/Users/Oliver/duke/2014_Gates'	
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
	#DATA		<<- '/Users/Oliver/duke/2014_Gates'	
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
	#DATA		<<- '/Users/Oliver/duke/2014_Gates'	
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
	#DATA			<<- '/Users/Oliver/duke/2014_Gates'
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
		tmp				<- dist.dna( seq, model='raw' )
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
	#DATA			<<- '/Users/Oliver/duke/2014_Gates'
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
		tmp				<- dist.dna( seq, model='raw' )
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
	indir		<- '/Users/Oliver/duke/2014_Gates/methods_comparison_pipeline/140914'  
	outdir		<- '/Users/Oliver/duke/2014_Gates/methods_comparison_pipeline/140914'	
	
	
	infiles		<- list.files(indir, '.*INTERNAL.R$', full.names=FALSE)
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
		#	get R2 for df.seq.pol
		#
		seq				<- df.seq.pol
		seq				<- rbind(seq, outgroup.seq.pol[, seq_len(ncol(seq))])
		#	get NJ tree	
		tmp				<- dist.dna(seq, model='raw')
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
	}
	#
	#	get R2 for concatenated genome
	#
	seq				<- cbind(df.seq.gag,df.seq.pol,df.seq.env)
	tmp				<- cbind(outgroup.seq.gag[,1:ncol(df.seq.gag)], outgroup.seq.pol, outgroup.seq.env)
	seq				<- rbind(seq,tmp)
	#	get NJ tree	
	tmp				<- dist.dna(seq, model='raw')
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
##	olli 18.09.14
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.scenarios<- function()
{	
	if(0)
	{
		indir			<- system.file(package="rPANGEAHIVsim", "misc")		
		infile.ind		<- '170914_HPTN071_scA_rep1'
		infile.trm		<- '170914_HPTN071_scA_rep1'	
		label			<- ''
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
														s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.2, 
														epi.model='HPTN071', epi.dt=1/48, epi.import=0.1,
														v.N0tau=1, v.r=2.851904, v.T50=-2,
														wher.mu=log(0.006716145)-0.37^2/2, wher.sigma=0.37,
														bwerm.mu=log(0.002239075)-0.05^2/2, bwerm.sigma=0.05)	
		tmpdir			<- '/Users/Oliver/git/HPTN071sim/HscABase140918'
		tmpdir			<- paste(tmpdir,label,sep='')
		dir.create(tmpdir, showWarnings=FALSE)						
		#						
		file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
		file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
		file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
		#system(file)
	}
}
##--------------------------------------------------------------------------------------------------------
##	rename simulated sequences
##	olli 27.10.14
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.November2014.renameVillage<- function()
{
	indir	<- '/Users/Oliver/duke/2014_Gates/meeting_2014_09_DSPSupdate/CD4-toBlind'	
	infiles	<- data.table(FILE= list.files(path=indir))
	
	set(infiles, NULL, 'EPI', infiles[, sapply(strsplit(FILE, '_'),'[[',1)])
	set(infiles, NULL, 'EPI', infiles[, factor(EPI, levels=c('A','B','C'), labels=c(1,2,3))])
	set(infiles, NULL, 'TIMEPOINT', infiles[, sapply(strsplit(FILE, '_'),'[[',3)])
	set(infiles, NULL, 'TIMEPOINT', infiles[, substr(TIMEPOINT,1,1)])
	set(infiles, NULL, 'END', infiles[, sapply(strsplit(FILE, '.', fixed=TRUE),'[[',2)])
	
	tmp		<- unique(subset(infiles, select=c(EPI, TIMEPOINT)))	
	set(tmp, NULL, 'TIMEPOINT_BLINDED', c('B','A','C',  'D','F','E',   'I','H','G'))	
	infiles	<- merge(infiles, tmp, by=c('EPI','TIMEPOINT'))
	
	set(infiles, NULL, 'SAMPLE', 1)
	set(infiles, infiles[, which(EPI==3 & grepl('60', FILE))], 'SAMPLE', 2)
	
	set(infiles, NULL, 'TO', infiles[, paste('281014','_Village_sc',TIMEPOINT_BLINDED,'_sample',SAMPLE,'_epi',EPI,'.',END,sep='')])
	
	infiles[, file.rename(paste(indir,'/',FILE,sep=''), paste(indir,'/',TO,sep='')), by='FILE']
	
}
##--------------------------------------------------------------------------------------------------------
##	rename simulated sequences
##	olli 26.10.14
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.November2014.rename<- function()
{
	indir	<- '/Users/Oliver/git/Pangea__Oct27_2014_code'
	indir	<- '/Users/Oliver/git/Pangea__Oct27_2014_code'
	infiles	<- data.table(FILE= list.files(path=indir, pattern='*SCENARIO.*.csv'))
	
	set(infiles, NULL, 'DATE', infiles[, substr(FILE, 1, 6)])
	set(infiles, NULL, 'END', infiles[, substr(FILE, nchar(FILE)-6, nchar(FILE))])
	set(infiles, NULL, 'REP', infiles[, regmatches(FILE,regexpr('RUN[0-9]+',FILE))])
	set(infiles, NULL, 'REP', infiles[, as.numeric(substr(REP, 4, nchar(REP)))])
	set(infiles, NULL, 'SCENARIO', infiles[, regmatches(FILE,regexpr('SCENARIO_[0-3]',FILE))])
	set(infiles, NULL, 'SCENARIO', infiles[, substr(SCENARIO, nchar(SCENARIO), nchar(SCENARIO))])
	
	set(infiles, infiles[, which(SCENARIO=='0')], 'SCENARIO', 'A')
	set(infiles, infiles[, which(SCENARIO=='3')], 'SCENARIO', 'B')
	set(infiles, infiles[, which(SCENARIO=='2')], 'SCENARIO', 'C')
	
	set(infiles, NULL, 'TO', infiles[, paste(DATE,'_HPTN071_sc',SCENARIO,'_rep',REP,'_',END,sep='')])
	
	infiles[, file.rename(paste(indir,'/',FILE,sep=''), paste(indir,'/',TO,sep='')), by='FILE']
}
##--------------------------------------------------------------------------------------------------------
##	check simulated sequences: create ExaML tree and estimate R2
##	olli 03.11.14
##--------------------------------------------------------------------------------------------------------
project.PANGEA.logo<- function()
{	
	treelabel.idx.sep	<- '|'
	infile				<- "/Users/Oliver/Dropbox\ (Infectious Disease)/PANGEAHIVsim_internal/freeze_Nov14/freeze1028_HPTN071_seqsimA/281014_HPTN071_scA_rep1_INTERNAL/281014_HPTN071_scA_rep1_SIMULATED_INTERNAL.R"
	load(infile)
		
	
	#	concatenate sequences
	tmp				<- tolower(do.call('rbind',strsplit(df.seq[, paste(GAG,POL,ENV,sep='')],'')))
	rownames(tmp)	<- df.seq[, paste(IDCLU,treelabel.idx.sep,LABEL,sep='')]	
	seq				<- as.DNAbin(tmp)
	#tmp				<- cbind(outgroup.seq.gag[,1:ncol(df.seq.gag)], outgroup.seq.pol, outgroup.seq.env)
	#seq				<- rbind(seq,tmp)	
	seq.ph			<- nj(dist.dna(seq, model='raw'))		
	seq.ph			<- ladderize(seq.ph)
	#	plot
	file			<- paste('/Users/Oliver/duke/2014_Gates', '/', 'logo_ph.pdf', sep='')				
	pdf(file=file, w=10, h=10)
	plot(seq.ph, type='fan',show.tip.label=0)
	dev.off()	
	
}
##--------------------------------------------------------------------------------------------------------
##	check simulated sequences: create ExaML tree and estimate R2
##	olli 23.10.14
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.November2014.plotrepclicates<- function()
{		
	indir	<- '/Users/Oliver/git/HPTN071sim/freeze1028'
	infiles	<- list.files(path=indir, pattern='*INTERNAL.R', recursive=TRUE)
	infiles	<- data.table(FILE=infiles[grepl('281014',infiles)])
	infiles[, SCENARIO:=regmatches(FILE,regexpr('sc[A-C]',FILE))]
	infiles[, REP:=regmatches(FILE,regexpr('rep[0-9]+',FILE))]
	set(infiles, NULL, 'SCENARIO', infiles[, substr(SCENARIO, 3, 3)])
	set(infiles, NULL, 'REP', infiles[, as.numeric(substr(REP, 4, nchar(REP)))])
	set(infiles, NULL, 'REP2', 1:3)
	setkey(infiles, SCENARIO, REP)
	
	
	df.epi	<- infiles[, {
								 z<- load( paste(indir, '/', FILE, sep='') )
								 df.epi[, SCENARIO:=SCENARIO]
								 df.epi[, REP2:=REP2]
								 df.epi				 
							}, by=c('SCENARIO','REP')]
					
	ggplot( df.epi, aes(x=YR, y=INC, group=interaction(SCENARIO,REP), colour=SCENARIO) ) + geom_line() + theme_bw() + 
			scale_x_continuous(limits=c(1980,2020), breaks=seq(1980,2020,2)) +
			scale_y_continuous(breaks=seq(0,2000,100))
	
	
	ggplot( df.epi, aes(x=YR, y=INCp*100, group=interaction(SCENARIO,REP), colour=interaction(SCENARIO,REP)) ) + geom_line(size=1.5) + theme_bw() +
			#scale_fill_manual (values=colours, guide=FALSE) +
			scale_colour_hue(c=50, l=60, guide=FALSE) +
			#scale_colour_brewer(palette='Accent', guide=FALSE) +
			scale_x_continuous(limits=c(1980,2020), breaks=seq(1980,2020,10), expand=c(0,0)) +
			scale_y_continuous(breaks=seq(0,5,0.5)) +
			labs(x='', y='% Incidence') +
			theme(panel.grid.major=element_line(colour="grey30", size=0.2))
	
	file<- paste(indir,'/inc_scenarios.pdf',sep='')
	ggsave(file=file, width=6, height=4)
}
##--------------------------------------------------------------------------------------------------------
##	check simulated sequences: create ExaML tree and estimate R2
##	olli 26.01.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Feb2015.dev<- function()
{	
	if(0)
	{
		#check simu tree
		file	<- '~/git/HPTN071sim/tmp151027-o0PU/150127_RUN001-o0PU_SIMULATED_DATEDTREE/150127_RUN001-o0PU_DATEDTREE.newick'
		ph		<- read.tree(file)
		df		<- data.table( 	LABEL=ph$tip.label,
								TIME_SEQ= as.numeric(sapply( strsplit(ph$tip.label, '|' , fixed=TRUE), '[[', 4 )),					
								TIME2ROOT=node.depth.edgelength(ph)[seq_len(Ntip(ph))] )
		ggplot(df, aes(x=TIME_SEQ, y=TIME2ROOT)) + geom_point()		
		#straight line as expected
	}
	if(0)
	{
		indir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/PANGEAHIVsim_internal/freeze_Jan15/regional/150125'
		infile	<- 'Annual_outputs_CL19_SA_A_V0_Rand1_Run48.csv'
		
		file	<- paste(indir,infile,sep='/')
		df		<- as.data.table(read.csv(file))
		set(df, NULL, 'X', NULL)
		tmp		<- names(df)[ which(as.character(sapply(df, class))!='numeric') ]
		for(x in tmp)
			set(df, NULL, x, as.numeric(df[[x]]))
		df[, 'Percent\nINC':= df[, NewCasesThisYear/(TotalPopulation-Prevalence)]]
		
		df		<- melt(df, id.vars=c('Year'))
		set(df, NULL, 'variable', df[, gsub('_','\n',variable)])
		set(df, NULL, 'variable', df[, gsub('MeanN','MeanN\n',variable)])
		set(df, NULL, 'variable', df[, gsub('NewCases','NewCases\n',variable)])	
		set(df, NULL, 'variable', df[, gsub('NHIVTested','NHIVTested\n',variable)])	
		set(df, NULL, 'variable', df[, gsub('Cumulative','Cumulative\n',variable)])
		set(df, NULL, 'variable', df[, gsub('Number','Number\n',variable)])
		set(df, NULL, 'variable', df[, gsub('PropHIV','PropHIV\n',variable)])
		set(df, NULL, 'variable', df[, gsub('PropAnnual','PropAnnual\n',variable)])
		#	paste(df[, unique(variable)], collapse='\',\'')
		tmp		<- c('Prevalence','Number\nPositive','NewCases\nThisYear','Percent\nINC','PropAnnual\nAcute','PropHIV\nPosONART','NAnnual','TotalPopulation','Number\nPositiveM','PopulationM','Number\nPositiveF','PopulationF','Cumulative\nNonPopartHIVtests','Cumulative\nPopartHIVtests','Cumulative\nNonPopartCD4tests','Cumulative\nPopartCD4tests','NHIVTested\nThisYear','NOnARTM','NNeedARTM','NOnARTF','NNeedARTF','PropMenCirc','Prop\nriskLow','Prop\nriskMed','Prop\nriskHigh','Prevalence\nriskLow','Prevalence\nriskMed','Prevalence\nriskHigh','MeanN\ncurrentpart\nriskLow','MeanN\ncurrentpart\nriskMed','MeanN\ncurrentpart\nriskHigh','MeanN\ncurrentsdpart\nriskLow','MeanN\ncurrentsdpart\nriskMed','MeanN\ncurrentsdpart\nriskHigh','MeanN\newPartnersthisyear\nriskLow','MeanN\newPartnersthisyear\nriskMed','MeanN\newPartnersthisyear\nriskHigh','MeanN\nlifetimepart\nriskLow','MeanN\nlifetimepart\nriskMed','MeanN\nlifetimepart\nriskHigh')
		set(df, NULL, 'variable', df[, factor(variable, levels=tmp, labels=tmp)])
		
		ggplot(df, aes(x=Year, y=value, group=variable)) + geom_step(with.guide=FALSE) + 
				facet_grid(variable~., scales='free_y')  + 
				scale_x_continuous(breaks=seq(1950,2050,by=10)) + theme_bw() +
				theme(strip.text=element_text(size=7))
		file	<- paste(indir,gsub('\\.csv','\\.pdf',infile),sep='/')
		ggsave(file=file, w=8, h=40)	
	}
	#
	#	new epi model; new Sampling model based on diagnoses; prelim scenario Acute Low + ART Fast
	#
	if(0)
	{
		indir			<- '/Users/Oliver/git/HPTN071sim/source/rPANGEAHIVsim/inst/misc'
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42, s.MODEL='Prop2DiagB4I',
				s.PREV.max=0.08, s.INTERVENTION.start=2015, s.INTERVENTION.mul= 2, s.ARCHIVAL.n=50,
				epi.model='HPTN071', epi.dt=1/48, epi.import=0.05,
				v.N0tau=1, v.r=2.851904, v.T50=-2,
				wher.mu=log(0.003358613)-0.3^2/2, wher.sigma=0.3, bwerm.mu=log(0.002239075)-0.13^2/2, bwerm.sigma=0.13, er.gamma=NA,
				sp.prop.of.sexactive= 0.05, report.prop.recent=1.0, 
				dbg.GTRparam=0, dbg.rER=0, index.starttime.mode='fix1955', startseq.mode='one', seqtime.mode='AtDiag')								
		
		# proposed standard run and control simulation
		pipeline.vary	<- data.table(	seqtime.mode=	c('AtDiag','AtDiag','AtDiag','AtDiag','AtDiag','AtTrm'),	
										s.MODEL=		c('Prop2DiagB4I','Prop2DiagB4I','Prop2DiagB4I','Prop2DiagB4I','Prop2DiagB4I','Prop2Untreated'),
										startseq.mode=	c('one','one', 'one','many','one','one'),
										epi.import=		c(0.05, 0.05, 0.05, 0.05, 0, 0),
										wher.mu=		c(log(0.003358613), log(0.003358613)-0.3^2/2, log(0.003358613)-0.3^2/2, log(0.003358613)-0.3^2/2, log(0.002239075)-0.13^2/2, log(0.002239075)-0.13^2/2),
										wher.sigma=		c(0, 0.3, 0.3, 0.3,0.13,0.13),
										bwerm.mu=		c(log(0.002239075), log(0.002239075)-0.13^2/2, log(0.002239075)-0.13^2/2, log(0.002239075)-0.13^2/2, log(0.002239075)-0.13^2/2, log(0.002239075)-0.13^2/2),
										bwerm.sigma=	c(0, 0.13, 0.13, 0.13, 0.13, 0.13), 
										dbg.GTRparam= 	c(1, 1, 0, 0, 0, 0),
										dbg.rER=		c(1, 1, 0, 0, 0, 0),										
										label=			c('-o5111DI', '-o5011DI', '-o5DI','-m5DI','-o0DI','-o0PU'))						
		dummy			<- pipeline.vary[, {				
											set(pipeline.args, which( pipeline.args$stat=='seqtime.mode' ), 'v', as.character(seqtime.mode))
											set(pipeline.args, which( pipeline.args$stat=='s.MODEL' ), 'v', as.character(s.MODEL))											
											set(pipeline.args, which( pipeline.args$stat=='startseq.mode' ), 'v', as.character(startseq.mode))
											set(pipeline.args, which( pipeline.args$stat=='epi.import' ), 'v', as.character(epi.import))
											set(pipeline.args, which( pipeline.args$stat=='wher.mu' ), 'v', as.character(wher.mu))
											set(pipeline.args, which( pipeline.args$stat=='wher.sigma' ), 'v', as.character(wher.sigma))
											set(pipeline.args, which( pipeline.args$stat=='bwerm.mu' ), 'v', as.character(bwerm.mu))
											set(pipeline.args, which( pipeline.args$stat=='bwerm.sigma' ), 'v', as.character(bwerm.sigma))											
											set(pipeline.args, which( pipeline.args$stat=='dbg.GTRparam' ), 'v', as.character(dbg.GTRparam))
											set(pipeline.args, which( pipeline.args$stat=='dbg.rER' ), 'v', as.character(dbg.rER))
											
											print(pipeline.args)
											#	scenario Acute Low + ART Fast	
											infile.ind		<- '150127_RUN001'
											infile.trm		<- '150127_RUN001'
											tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp151030'
											tmpdir			<- paste(tmpdir,label,sep='')
											dir.create(tmpdir, showWarnings=FALSE)																													
											file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind, label,'_IND.csv',sep=''))
											file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm, label,'_TRM.csv',sep=''))
											file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
											system(file)
										}, by='label']
	}
	#
	#	check control: why so much variation in root to tip divergence?
	#
	if(0)
	{
		indir			<- '/Users/Oliver/git/HPTN071sim/source/rPANGEAHIVsim/inst/misc'
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42, s.MODEL='Prop2DiagB4I',
				s.PREV.max=0.08, s.INTERVENTION.start=2015, s.INTERVENTION.mul= 2, s.ARCHIVAL.n=50,
				epi.model='HPTN071', epi.dt=1/48, epi.import=0.05,
				v.N0tau=1, v.r=2.851904, v.T50=-2,
				wher.mu=log(0.003358613)-0.3^2/2, wher.sigma=0.3, bwerm.mu=log(0.002239075)-0.13^2/2, bwerm.sigma=0.13, er.gamma=NA,
				dbg.GTRparam=0, dbg.rER=0, index.starttime.mode='fix1955', startseq.mode='one', seqtime.mode='AtDiag')								
		
		# proposed standard run and control simulation
		pipeline.vary	<- data.table(	seqtime.mode=	c('AtTrm','AtTrm','AtTrm','AtTrm'),	
				s.MODEL=		c('Prop2Untreated','Prop2Untreated','Prop2Untreated','Prop2Untreated'),
				startseq.mode=	c('one','one','one','one'),
				epi.import=		c(0, 0, 0, 0),
				wher.mu=		c(log(0.002239075)-0.13^2/2, log(0.002239075)-0.13^2/2, log(0.002239075)-0.13^2/2, log(0.002239075)),
				wher.sigma=		c(0.13, 0.13, 0.13, 0),
				dbg.GTRparam= 	c(0, 1, 1, 1),
				dbg.rER=		c(0, 0, 1, 1),
				label=			c('-o0PU000','-o0PU010', '-o0PU011', '-o0PU111'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='seqtime.mode' ), 'v', as.character(seqtime.mode))
					set(pipeline.args, which( pipeline.args$stat=='s.MODEL' ), 'v', as.character(s.MODEL))											
					set(pipeline.args, which( pipeline.args$stat=='startseq.mode' ), 'v', as.character(startseq.mode))
					set(pipeline.args, which( pipeline.args$stat=='epi.import' ), 'v', as.character(epi.import))
					set(pipeline.args, which( pipeline.args$stat=='wher.mu' ), 'v', as.character(wher.mu))
					set(pipeline.args, which( pipeline.args$stat=='wher.sigma' ), 'v', as.character(wher.sigma))					
					set(pipeline.args, which( pipeline.args$stat=='dbg.GTRparam' ), 'v', as.character(dbg.GTRparam))
					set(pipeline.args, which( pipeline.args$stat=='dbg.rER' ), 'v', as.character(dbg.rER))
					
					print(pipeline.args)
					#	scenario Acute Low + ART Fast	
					infile.ind		<- '150127_RUN001'
					infile.trm		<- '150127_RUN001'
					tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp150130'
					tmpdir			<- paste(tmpdir,label,sep='')
					dir.create(tmpdir, showWarnings=FALSE)																													
					file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind, label,'_IND.csv',sep=''))
					file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm, label,'_TRM.csv',sep=''))
					file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
					#system(file)
				}, by='label']
	}
	#
	#	after bugfixes, R2 is 66%! Need to add more variability ..
	#
	if(1)
	{
		indir			<- '/Users/Oliver/git/HPTN071sim/source/rPANGEAHIVsim/inst/misc'
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42, s.MODEL='Prop2DiagB4I',
				s.PREV.max=0.08, s.INTERVENTION.start=2015, s.INTERVENTION.mul= 2, s.ARCHIVAL.n=50,
				epi.model='HPTN071', epi.dt=1/48, epi.import=0.05,
				v.N0tau=1, v.r=2.851904, v.T50=-2,
				wher.mu=log(0.003358613)-0.3^2/2, wher.sigma=0.3, bwerm.mu=log(0.002239075)-0.13^2/2, bwerm.sigma=0.13, er.gamma=NA,
				dbg.GTRparam=0, dbg.rER=0, index.starttime.mode='fix1955', startseq.mode='one', seqtime.mode='AtDiag')								
		
		# proposed standard run and control simulation
		pipeline.vary	<- data.table(	seqtime.mode=	c('AtDiag', 'AtDiag','AtDiag','AtDiag','Exp3','Exp3'),	
										s.MODEL=		c('Prop2DiagB4I','Prop2DiagB4I','Prop2DiagB4I','Prop2DiagB4I','Prop2DiagB4I','Prop2DiagB4I'),									
										wher.mu=		c(log(0.00447743)-0.3^2/2, log(0.00447743)-0.5^2/2, log(0.00447743)-0.3^2/2, log(0.00447743)-0.5^2/2, log(0.002239075)-0.13^2/2, log(0.002239075)-0.13^2/2),
										wher.sigma=		c(0.3, 0.5, 0.3, 0.5, 0.3, 0.5),
										er.gamma=		c(NA, NA, 4, 4, 4, 4),
										label=			c('-o30DI', '-o50DI', '-o34DI','-o54DI','-o34E3','-o54E3'))						
		dummy			<- pipeline.vary[, {				
						set(pipeline.args, which( pipeline.args$stat=='seqtime.mode' ), 'v', as.character(seqtime.mode))
						set(pipeline.args, which( pipeline.args$stat=='s.MODEL' ), 'v', as.character(s.MODEL))											
						set(pipeline.args, which( pipeline.args$stat=='wher.mu' ), 'v', as.character(wher.mu))
						set(pipeline.args, which( pipeline.args$stat=='wher.sigma' ), 'v', as.character(wher.sigma))
						set(pipeline.args, which( pipeline.args$stat=='er.gamma' ), 'v', as.character(er.gamma))
						
						print(pipeline.args)
						#	scenario Acute Low + ART Fast	
						infile.ind		<- '150127_RUN001'
						infile.trm		<- '150127_RUN001'
						tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp150130b'
						tmpdir			<- paste(tmpdir,label,sep='')
						dir.create(tmpdir, showWarnings=FALSE)																													
						file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind, label,'_IND.csv',sep=''))
						file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm, label,'_TRM.csv',sep=''))
						file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
						#system(file)
					}, by='label']
	}
	#
	#	-o54DI brings R2 down to 51%, need more
	#
	if(1)
	{
		indir			<- '/Users/Oliver/git/HPTN071sim/source/rPANGEAHIVsim/inst/misc'
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42, s.MODEL='Prop2DiagB4I',
				s.PREV.max=0.08, s.INTERVENTION.start=2015, s.INTERVENTION.mul= 2, s.ARCHIVAL.n=50,
				epi.model='HPTN071', epi.dt=1/48, epi.import=0.05,
				v.N0tau=1, v.r=2.851904, v.T50=-2,
				wher.mu=log(0.003358613)-0.3^2/2, wher.sigma=0.3, bwerm.mu=log(0.002239075)-0.13^2/2, bwerm.sigma=0.13, er.gamma=NA,
				dbg.GTRparam=0, dbg.rER=0, index.starttime.mode='fix1955', startseq.mode='one', seqtime.mode='AtDiag')								
		
		# proposed standard run and control simulation
		pipeline.vary	<- data.table(	seqtime.mode=	c('AtDiag','DUnif5', 'DUnif9','DUnif5','DUnif9','DUnif5','DUnif9'),	
				s.MODEL=		c('Prop2DiagB4I','Prop2DiagB4I','Prop2DiagB4I','Prop2DiagB4I','Prop2DiagB4I','Prop2DiagB4I','Prop2DiagB4I'),									
				wher.mu=		c(log(0.00447743)-0.5^2/2, log(0.00447743)-0.5^2/2, log(0.00447743)-0.5^2/2, log(0.00447743)-0.5^2/2, log(0.00447743)-0.5^2/2, log(0.00447743)-0.5^2/2, log(0.00447743)-0.5^2/2),
				wher.sigma=		c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
				er.gamma=		c(4, 4, 4, 4, 4, 4, 4),
				index.starttime.mode= c('fix1955','fix1955','fix1955','fix1955','fix1955','normal','normal'),
				startseq.mode=	c('one','one','one','many','many','many','many'),
				label=			c('-o54DI', '-o54D5', '-o54D9', '-m54D5', '-m54D9', '-f54D5', '-f54D9'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='seqtime.mode' ), 'v', as.character(seqtime.mode))
					set(pipeline.args, which( pipeline.args$stat=='s.MODEL' ), 'v', as.character(s.MODEL))											
					set(pipeline.args, which( pipeline.args$stat=='wher.mu' ), 'v', as.character(wher.mu))
					set(pipeline.args, which( pipeline.args$stat=='wher.sigma' ), 'v', as.character(wher.sigma))
					set(pipeline.args, which( pipeline.args$stat=='er.gamma' ), 'v', as.character(er.gamma))					
					set(pipeline.args, which( pipeline.args$stat=='index.starttime.mode' ), 'v', as.character(index.starttime.mode))
					set(pipeline.args, which( pipeline.args$stat=='startseq.mode' ), 'v', as.character(startseq.mode))
					
					print(pipeline.args)
					#	scenario Acute Low + ART Fast	
					infile.ind		<- '150127_RUN001'
					infile.trm		<- '150127_RUN001'
					tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp150131b'
					tmpdir			<- paste(tmpdir,label,sep='')
					dir.create(tmpdir, showWarnings=FALSE)																													
					file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind, label,'_IND.csv',sep=''))
					file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm, label,'_TRM.csv',sep=''))
					file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
					#system(file)
				}, by='label']
	}
	#
	#	-m54D5 brings R2 down to 42%. not too bad. try and get div overall a bit higher
	#
	if(0)
	{
		indir			<- '/Users/Oliver/git/HPTN071sim/source/rPANGEAHIVsim/inst/misc'
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42, s.MODEL='Prop2DiagB4I',
				s.PREV.max=0.08, s.INTERVENTION.start=2015, s.INTERVENTION.mul= 2, s.ARCHIVAL.n=50,
				epi.model='HPTN071', epi.dt=1/48, epi.import=0.05,
				v.N0tau=1, v.r=2.851904, v.T50=-2,
				wher.mu=log(0.00447743)-0.5^2/2, wher.sigma=0.5, bwerm.mu=log(0.002239075)-0.13^2/2, bwerm.sigma=0.13, er.gamma=4,
				dbg.GTRparam=0, dbg.rER=0, index.starttime.mode='fix1955', startseq.mode='many', seqtime.mode='AtDiag')								
		
		# proposed standard run and control simulation
		pipeline.vary	<- data.table(	seqtime.mode=	c('DUnif5','DUnif5','DUnif5','DUnif5','DUnif5','DUnif5','AtTrm'),	
				s.MODEL=				c('Prop2DiagB4I','Prop2DiagB4I','Prop2DiagB4I','Prop2DiagB4I','Prop2DiagB4I','Prop2DiagB4I','Prop2Untreated'),													
				index.starttime.mode= 	c('fix1955','fix1955','fix1955','fix1955','normal','normal','fix1955'),
				bwerm.mu=				c(log(0.002239075)-0.16^2/2, log(0.002239075)-0.2^2/2, log(0.002239075)-0.25^2/2, log(0.002239075)-0.3^2/2, log(0.002239075)-0.2^2/2, log(0.002239075)-0.3^2/2,log(0.002239075)-0.16^2/2),
				bwerm.sigma=			c(0.16, 0.2, 0.25, 0.3, 0.2, 0.3, 0.16),
				startseq.mode=			c('many','many','many','many','many','many','one'),
				label=					c('-m16', '-m20', '-m25', '-m30', '-n20', '-n30', '-o16'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='seqtime.mode' ), 'v', as.character(seqtime.mode))
					set(pipeline.args, which( pipeline.args$stat=='s.MODEL' ), 'v', as.character(s.MODEL))											
					set(pipeline.args, which( pipeline.args$stat=='bwerm.mu' ), 'v', as.character(bwerm.mu))
					set(pipeline.args, which( pipeline.args$stat=='bwerm.sigma' ), 'v', as.character(bwerm.sigma))					
					set(pipeline.args, which( pipeline.args$stat=='index.starttime.mode' ), 'v', as.character(index.starttime.mode))
					set(pipeline.args, which( pipeline.args$stat=='startseq.mode' ), 'v', as.character(startseq.mode))
					print(pipeline.args)
					#	scenario Acute Low + ART Fast	
					infile.ind		<- '150127_RUN001'
					infile.trm		<- '150127_RUN001'
					tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp150201'
					tmpdir			<- paste(tmpdir,label,sep='')
					dir.create(tmpdir, showWarnings=FALSE)																													
					file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind, label,'_IND.csv',sep=''))
					file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm, label,'_TRM.csv',sep=''))
					file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
					#system(file)
				}, by='label']
	}
	#
	#	-m30 brings R2 down to 32%. that s seem OK to use.
	#
	if(0)
	{
		indir			<- '/Users/Oliver/git/HPTN071sim/source/rPANGEAHIVsim/inst/misc'
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1985, yr.end=2020, seed=42, s.MODEL='Prop2DiagB4I',
				s.PREV.max=0.08, s.INTERVENTION.start=2015, s.INTERVENTION.mul= 2, s.ARCHIVAL.n=50,
				epi.model='HPTN071', epi.dt=1/48, epi.import=0.05,
				v.N0tau=1, v.r=2.851904, v.T50=-2,
				wher.mu=log(0.00447743)-0.5^2/2, wher.sigma=0.5, bwerm.mu=log(0.002239075)-0.3^2/2, bwerm.sigma=0.3, er.gamma=4,
				dbg.GTRparam=0, dbg.rER=0, index.starttime.mode='fix1955', startseq.mode='many', seqtime.mode='DUnif5')								
		
		# proposed standard run and control simulation
		pipeline.vary	<- data.table(	seqtime.mode=	c('DUnif5','AtDiag','AtDiag','AtDiag','AtTrm'),	
				s.MODEL=				c('Prop2DiagB4I','Prop2DiagB4I','Prop2DiagB4I','Prop2DiagB4I','Prop2Untreated'),													
				startseq.mode=			c('many','many','many','many','many'),
				s.multi=				c(1,1,2,4,1),
				label=					c('-mD5', '-mDI','-mDIs2x','-mDIs4x', '-mPU'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='seqtime.mode' ), 'v', as.character(seqtime.mode))
					set(pipeline.args, which( pipeline.args$stat=='s.MODEL' ), 'v', as.character(s.MODEL))											
					set(pipeline.args, which( pipeline.args$stat=='startseq.mode' ), 'v', as.character(startseq.mode))
					print(pipeline.args)						
					if(1)
					{
						#	scenario A					
						infile.ind		<- '150129_HPTN071_scA'
						infile.trm		<- '150129_HPTN071_scA'
						tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp150203-A'
						tmpdir			<- paste(tmpdir,label,sep='')
						dir.create(tmpdir, showWarnings=FALSE)																		
						set(pipeline.args, which( pipeline.args$stat=='s.PREV.max' ), 'v', as.character(s.multi*0.073))
						file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
						file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
						file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
						#system(file)													
					}
					if(1)
					{
						#	scenario B					
						infile.ind		<- '150129_HPTN071_scB'
						infile.trm		<- '150129_HPTN071_scB'
						tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp150203-B'
						tmpdir			<- paste(tmpdir,label,sep='')
						dir.create(tmpdir, showWarnings=FALSE)																		
						set(pipeline.args, which( pipeline.args$stat=='s.PREV.max' ), 'v', as.character(s.multi*0.084))
						file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
						file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
						file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
						#system(file)													
					}
					if(1)
					{
						#	scenario C					
						infile.ind		<- '150129_HPTN071_scC'
						infile.trm		<- '150129_HPTN071_scC'
						tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp150203-C'
						tmpdir			<- paste(tmpdir,label,sep='')
						dir.create(tmpdir, showWarnings=FALSE)																		
						set(pipeline.args, which( pipeline.args$stat=='s.PREV.max' ), 'v', as.character(s.multi*0.075))
						file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
						file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
						file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
						#system(file)													
					}
					if(1)
					{
						#	scenario D					
						infile.ind		<- '150129_HPTN071_scD'
						infile.trm		<- '150129_HPTN071_scD'
						tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp150203-D'
						tmpdir			<- paste(tmpdir,label,sep='')
						dir.create(tmpdir, showWarnings=FALSE)																		
						set(pipeline.args, which( pipeline.args$stat=='s.PREV.max' ), 'v', as.character(s.multi*0.08))
						file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
						file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
						file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
						#system(file)													
					}
					if(1)
					{
						#	scenario E					
						infile.ind		<- '150129_HPTN071_scE'
						infile.trm		<- '150129_HPTN071_scE'
						tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp150203-E'
						tmpdir			<- paste(tmpdir,label,sep='')
						dir.create(tmpdir, showWarnings=FALSE)																		
						set(pipeline.args, which( pipeline.args$stat=='s.PREV.max' ), 'v', as.character(s.multi*0.078))
						file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
						file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
						file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
						#system(file)													
					}
					if(1)
					{
						#	scenario F					
						infile.ind		<- '150129_HPTN071_scF'
						infile.trm		<- '150129_HPTN071_scF'
						tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp150203-F'
						tmpdir			<- paste(tmpdir,label,sep='')
						dir.create(tmpdir, showWarnings=FALSE)																		
						set(pipeline.args, which( pipeline.args$stat=='s.PREV.max' ), 'v', as.character(s.multi*0.083))
						file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
						file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
						file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
						#system(file)													
					}					
				}, by='label']
		}
		#
		#	fixed seq sampling fractions etc
		#
		if(1)
		{
			indir			<- '/Users/Oliver/git/HPTN071sim/source/rPANGEAHIVsim/inst/misc'
			pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1985, yr.end=2020, seed=42, s.MODEL='Prop2DiagB4I', report.prop.recent=1.0,
					s.PREV.max=NA, s.INTERVENTION.prop=NA, s.INTERVENTION.start=2015, s.INTERVENTION.mul= NA, s.ARCHIVAL.n=50,
					epi.model='HPTN071', epi.dt=1/48, epi.import=0.05,
					v.N0tau=1, v.r=2.851904, v.T50=-2,
					wher.mu=log(0.00447743)-0.5^2/2, wher.sigma=0.5, bwerm.mu=log(0.002239075)-0.3^2/2, bwerm.sigma=0.3, er.gamma=4,
					dbg.GTRparam=0, dbg.rER=0, index.starttime.mode='fix1955', startseq.mode='many', seqtime.mode='AtDiag')								
			
			# proposed standard run and control simulation
			pipeline.vary	<- data.table(	seqtime.mode=	c('AtDiag','AtDiag','AtDiag','AtDiag','AtDiag','AtDiag','AtDiag','AtDiag','AtDiag'),	
					s.MODEL=				c('Prop2DiagB4I','Prop2DiagB4I','Prop2Untreated','Fixed2Prop','Fixed2Prop','Fixed2Prop','Fixed2Prop','Fixed2Prop','Fixed2Prop'),													
					startseq.mode=			c('many','many','many','many','many','many','many','many','many'),
					s.multi=				c(1,2,1,1,1,1,2,1,1),
					yr.end=					c(2020,2020,2020,2020,2020,2020,2020,2020,2018),
					epi.import=				c(0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.2,0.05),
					s.PREV.max.n=			c(1600,1600,1600,1600,1600,1600,1600,1600,1280),
					s.INTERVENTION.prop=	c(NA,NA,NA,0.85,0.5,0.25,0.5,0.5,0.375),
					s.INTERVENTION.mul=		c(1,1,1,NA,NA,NA,NA,NA,NA),
					label=					c('-mDI','-mDIs2x','-mPU','-mFP85','-mFP50','-mFP25','-mFP50s2x','-mFP50c20','-mFP50e18'))						
			dummy			<- pipeline.vary[, {				
						set(pipeline.args, which( pipeline.args$stat=='seqtime.mode' ), 'v', as.character(seqtime.mode))
						set(pipeline.args, which( pipeline.args$stat=='s.MODEL' ), 'v', as.character(s.MODEL))											
						set(pipeline.args, which( pipeline.args$stat=='startseq.mode' ), 'v', as.character(startseq.mode))
						set(pipeline.args, which( pipeline.args$stat=='s.INTERVENTION.mul' ), 'v', as.character(s.INTERVENTION.mul))
						set(pipeline.args, which( pipeline.args$stat=='s.INTERVENTION.prop' ), 'v', as.character(s.INTERVENTION.prop))
						set(pipeline.args, which( pipeline.args$stat=='epi.import' ), 'v', as.character(epi.import))
						print(pipeline.args)						
						if(1)
						{
							#	scenario A					
							infile.ind		<- '150129_HPTN071_scA'
							infile.trm		<- '150129_HPTN071_scA'
							tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp150204-A'
							tmpdir			<- paste(tmpdir,label,sep='')
							dir.create(tmpdir, showWarnings=FALSE)																								
							set(pipeline.args, which( pipeline.args$stat=='s.PREV.max.n' ), 'v', as.character(s.multi*s.PREV.max.n))
							file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
							file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
							file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
							#system(file)													
						}
						if(1)
						{
							#	scenario B					
							infile.ind		<- '150129_HPTN071_scB'
							infile.trm		<- '150129_HPTN071_scB'
							tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp150204-B'
							tmpdir			<- paste(tmpdir,label,sep='')
							dir.create(tmpdir, showWarnings=FALSE)																		
							set(pipeline.args, which( pipeline.args$stat=='s.PREV.max.n' ), 'v', as.character(s.multi*s.PREV.max.n))
							file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
							file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
							file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
							#system(file)													
						}
						if(1)
						{
							#	scenario C					
							infile.ind		<- '150129_HPTN071_scC'
							infile.trm		<- '150129_HPTN071_scC'
							tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp150204-C'
							tmpdir			<- paste(tmpdir,label,sep='')
							dir.create(tmpdir, showWarnings=FALSE)																		
							set(pipeline.args, which( pipeline.args$stat=='s.PREV.max.n' ), 'v', as.character(s.multi*s.PREV.max.n))
							file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
							file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
							file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
							#system(file)													
						}
						if(1)
						{
							#	scenario D					
							infile.ind		<- '150129_HPTN071_scD'
							infile.trm		<- '150129_HPTN071_scD'
							tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp150204-D'
							tmpdir			<- paste(tmpdir,label,sep='')
							dir.create(tmpdir, showWarnings=FALSE)																		
							set(pipeline.args, which( pipeline.args$stat=='s.PREV.max.n' ), 'v', as.character(s.multi*s.PREV.max.n))
							file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
							file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
							file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
							#system(file)													
						}
						if(1)
						{
							#	scenario E					
							infile.ind		<- '150129_HPTN071_scE'
							infile.trm		<- '150129_HPTN071_scE'
							tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp150204-E'
							tmpdir			<- paste(tmpdir,label,sep='')
							dir.create(tmpdir, showWarnings=FALSE)																		
							set(pipeline.args, which( pipeline.args$stat=='s.PREV.max.n' ), 'v', as.character(s.multi*s.PREV.max.n))
							file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
							file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
							file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
							#system(file)													
						}
						if(1)
						{
							#	scenario F					
							infile.ind		<- '150129_HPTN071_scF'
							infile.trm		<- '150129_HPTN071_scF'
							tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp150204-F'
							tmpdir			<- paste(tmpdir,label,sep='')
							dir.create(tmpdir, showWarnings=FALSE)																									
							set(pipeline.args, which( pipeline.args$stat=='s.PREV.max.n' ), 'v', as.character(s.multi*s.PREV.max.n))
							file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
							file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
							file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
							#system(file)													
						}					
					}, by='label']
		}
		#
		#	double check R2 -- now very low!
		#
		if(1)
		{
			indir			<- '/Users/Oliver/git/HPTN071sim/source/rPANGEAHIVsim/inst/misc'
			pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1985, yr.end=2020, seed=42, s.MODEL='Fixed2Prop', report.prop.recent=1.0,
					s.PREV.max.n=1600, s.INTERVENTION.prop=0.5, s.INTERVENTION.start=2015, s.INTERVENTION.mul= NA, s.ARCHIVAL.n=50,
					epi.model='HPTN071', epi.dt=1/48, epi.import=0.05, root.edge.fixed=0,
					v.N0tau=1, v.r=2.851904, v.T50=-2,
					wher.mu=log(0.00447743)-0.5^2/2, wher.sigma=0.5, bwerm.mu=log(0.002239075)-0.3^2/2, bwerm.sigma=0.3, er.gamma=4,
					dbg.GTRparam=0, dbg.rER=0, index.starttime.mode='fix1955', startseq.mode='many', seqtime.mode='AtDiag')								
			
			# proposed standard run and control simulation
			pipeline.vary	<- data.table(	wher.mu=				c(log(0.00447743)-0.5^2/2,log(0.00447743)-0.5^2/2,log(0.00447743)-0.5^2/2,log(0.00447743)-0.3^2/2,log(0.00447743)-0.5^2/2,log(0.00447743)-0.5^2/2,log(0.00447743)-0.5^2/2,log(0.00447743)-0.3^2/2,log(0.002239075)-0.16^2/2,log(0.00447743)-0.5^2/2,log(0.00447743)-0.5^2/2,log(0.00447743)-0.5^2/2,log(0.00447743)-0.5^2/2),
											wher.sigma=				c(0.5,0.5,0.5,0.3,0.5,0.5,0.5,0.3, 0.16, 0.5, 0.5, 0.5, 0.5),
											bwerm.mu=				c(log(0.002239075)-0.3^2/2, log(0.002239075)-0.2^2/2, log(0.002239075)-0.16^2/2, log(0.002239075)-0.16^2/2,log(0.002239075)-0.3^2/2, log(0.002239075)-0.2^2/2, log(0.002239075)-0.16^2/2, log(0.002239075)-0.16^2/2, log(0.002239075)-0.16^2/2,log(0.002239075)-0.3^2/2, log(0.002239075)-0.2^2/2,log(0.002239075)-0.3^2/2, log(0.002239075)-0.2^2/2),
											bwerm.sigma=			c(0.3, 0.2, 0.16, 0.16,0.3, 0.2, 0.16, 0.16, 0.16, 0.3, 0.2, 0.3, 0.2),
											startseq.mode=			c('many','many','many','many','one','one','one','one','one','one','one','one','one'),
											index.starttime.mode=	c('fix1955','fix1955','fix1955','fix1955','fix1955','fix1955','fix1955','fix1955','fix1955','fix1965','fix1965','fix1972','fix1972'),
											label=					c('-bh30','-bh20','-bh16','-wh30','-obh30','-obh20','-obh16','-owh30','-owh16','-obh3065','-obh2065','-obh3072','-obh2072'))						
			dummy			<- pipeline.vary[, {				
						set(pipeline.args, which( pipeline.args$stat=='wher.mu' ), 'v', as.character(wher.mu))
						set(pipeline.args, which( pipeline.args$stat=='wher.sigma' ), 'v', as.character(wher.sigma))											
						set(pipeline.args, which( pipeline.args$stat=='bwerm.mu' ), 'v', as.character(bwerm.mu))
						set(pipeline.args, which( pipeline.args$stat=='bwerm.sigma' ), 'v', as.character(bwerm.sigma))
						set(pipeline.args, which( pipeline.args$stat=='startseq.mode' ), 'v', as.character(startseq.mode))
						set(pipeline.args, which( pipeline.args$stat=='index.starttime.mode' ), 'v', as.character(index.starttime.mode))
						
						print(pipeline.args)						
						if(1)
						{
							#	scenario A					
							infile.ind		<- '150129_HPTN071_scA'
							infile.trm		<- '150129_HPTN071_scA'
							tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp150206-A'
							tmpdir			<- paste(tmpdir,label,sep='')
							dir.create(tmpdir, showWarnings=FALSE)																															
							file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
							file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
							file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
							#system(file)													
						}					
						if(1)
						{
							#	scenario E					
							infile.ind		<- '150129_HPTN071_scE'
							infile.trm		<- '150129_HPTN071_scE'
							tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp150206-E'
							tmpdir			<- paste(tmpdir,label,sep='')
							dir.create(tmpdir, showWarnings=FALSE)																									
							file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
							file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
							file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
							#system(file)													
						}											
					}, by='label']
		}
		#
		#	high R2 runs -- these are used for the tree comparison
		#	for tree comparison TRAIN1: seed=42
		#	for tree comparison TRAIN2 and higher: seed=101 
		#
		if(1)
		{			
			indir			<- '/Users/Oliver/git/HPTN071sim/source/rPANGEAHIVsim/inst/misc'
			pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1985, yr.end=2020, seed=42, s.MODEL='Fixed2Prop', report.prop.recent=1.0,
					s.PREV.max.n=1600, s.INTERVENTION.prop=0.25, s.INTERVENTION.start=2015, s.INTERVENTION.mul= NA, s.ARCHIVAL.n=50,
					epi.model='HPTN071', epi.dt=1/48, epi.import=0.05, root.edge.fixed=1,
					v.N0tau=1, v.r=2.851904, v.T50=-2,
					wher.mu=log(0.002239075)-0.3^2/2, wher.sigma=0.3, 
					bwerm.mu=log(0.002239075)-0.3^2/2, bwerm.sigma=0.3, er.gamma=4,
					dbg.GTRparam=1, dbg.rER=1, index.starttime.mode='fix1955', startseq.mode='many', seqtime.mode='AtDiag')								
			
			# proposed standard run and control simulation
			pipeline.vary	<- data.table(	startseq.mode=			c('many','one'),					
					label=					c('-m111','-o111'))						
			dummy			<- pipeline.vary[, {				
						set(pipeline.args, which( pipeline.args$stat=='startseq.mode' ), 'v', as.character(startseq.mode))
						
						print(pipeline.args)						
						if(0)
						{
							#	scenario F					
							infile.ind		<- '150129_HPTN071_scF'
							infile.trm		<- '150129_HPTN071_scF'
							tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp150227-F'
							tmpdir			<- paste(tmpdir,label,sep='')
							dir.create(tmpdir, showWarnings=FALSE)																															
							file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
							file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
							file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
							#system(file)													
						}					
						if(1)
						{
							#	scenario E					
							infile.ind		<- '150129_HPTN071_scE'
							infile.trm		<- '150129_HPTN071_scE'
							tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp150227-E'
							tmpdir			<- paste(tmpdir,label,sep='')
							dir.create(tmpdir, showWarnings=FALSE)																									
							file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
							file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
							file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
							#system(file)													
						}											
					}, by='label']
		}
}
##--------------------------------------------------------------------------------------------------------
##	olli 22.06.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.treecomparison<- function()
{		
	if(1)
	{
		indir			<- '/Users/Oliver/git/HPTN071sim/source/rPANGEAHIVsim/inst/misc'
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1985, yr.end=2020, seed=42, s.MODEL='Fixed2Prop', report.prop.recent=1.0,
				s.PREV.max.n=1600, s.INTERVENTION.prop=0.25, s.INTERVENTION.start=2015, s.INTERVENTION.mul= NA, s.ARCHIVAL.n=50,
				epi.model='HPTN071', epi.dt=1/48, epi.import=0.05, root.edge.fixed=1,
				v.N0tau=1, v.r=2.851904, v.T50=-2,
				wher.mu=log(0.002239075)-0.3^2/2, wher.sigma=0.3, 
				bwerm.mu=log(0.002239075)-0.3^2/2, bwerm.sigma=0.3, er.gamma=4,
				dbg.GTRparam=1, dbg.rER=1, index.starttime.mode='fix1955', startseq.mode='many', seqtime.mode='AtDiag')								
		
		# proposed standard run and control simulation
		pipeline.vary	<- data.table(	startseq.mode=			c('many','one'),					
				label=					c('-m111','-o111'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='startseq.mode' ), 'v', as.character(startseq.mode))
					
					print(pipeline.args)						
					if(1)
					{
						#	scenario F					
						infile.ind		<- '150129_HPTN071_scF'
						infile.trm		<- '150129_HPTN071_scF'
						tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp150227-F'
						tmpdir			<- paste(tmpdir,label,sep='')
						dir.create(tmpdir, showWarnings=FALSE)																															
						file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
						file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
						file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
						#system(file)													
					}					
					if(1)
					{
						#	scenario E					
						infile.ind		<- '150129_HPTN071_scE'
						infile.trm		<- '150129_HPTN071_scE'
						tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp150227-E'
						tmpdir			<- paste(tmpdir,label,sep='')
						dir.create(tmpdir, showWarnings=FALSE)																									
						file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
						file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
						file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
						#system(file)													
					}											
				}, by='label']
	}
}
##--------------------------------------------------------------------------------------------------------
##	olli 22.06.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.treecomparison.gaps<- function()
{	
	#
	#	prepare gappy PANGEA files at different coverage
	#
	if(0)
	{
		callconsensus.cmd	<- 'python /Users/Oliver/Dropbox\\ \\(Infectious\\ Disease\\)/PANGEA-BEEHIVE-SHARED/ChrisCode/CallConsensus.py'
		callconsensus.minc	<- 10
		callconsensus.maxc	<- 10
		mergealignments.cmd	<- 'python /Users/Oliver/Dropbox\\ \\(Infectious\\ Disease\\)/PANGEA-BEEHIVE-SHARED/ChrisCode/MergeAlignments2.0.py'
		mergealignments.opt	<- '-d'
		
		indir.refalgn		<- '/Users/Oliver/Dropbox\\ \\(Infectious\\ Disease\\)/PANGEA_data'
		infile.refalgn		<- 'HIV1_COM_2012_genome_DNA_WithExtraA1UG.fasta'
		indir.basefreq		<- '/Users/Oliver/Dropbox (Infectious Disease)/PANGEA_data/mapping/2681/BaseFreqs'
		outdir				<- '/Users/Oliver/git/HPTN071sim/treecomparison/PANGEAcov'
		
		tmp					<- c(10,20,30,40,50)
		tmp					<- c(70,90,100)
		tmp					<- c(200)
		invisible(sapply(tmp, function(callconsensus.minc)
						{
							outfile				<- '150623_PANGEAGlobal2681'
							outfile				<- paste(outfile, '_C', callconsensus.minc, '.fa', sep='')
							seq					<- PANGEA.align.from.basefreq(indir.refalgn, infile.refalgn, indir.basefreq, callconsensus.minc, outdir, verbose=1)
							write.dna(seq, file=paste(outdir, '/', outfile, sep=''), format='fasta', colsep='', nbcol=-1)
							save(seq, file=paste(outdir, '/', gsub('.fa','.R',outfile,fixed=TRUE), sep=''))
							NULL
						}))
		
		#	proportion of '?' calls
		files	<- list.files(outdir, pattern='PANGEAGlobal.*fa$')
		ss		<- lapply(files, function(x) read.dna(paste(outdir,'/',x,sep=''), format='fasta'))	
		dfc		<- lapply(seq_along(ss), function(i)
				{
					tmp		<- as.character(ss[[i]])
					data.table(FILE=files[i], SEQ=rownames(tmp), NSEQ=apply(tmp!='-', 1, sum), NOCOV=apply(tmp=='?' | tmp=='n', 1, sum))				
				})
		dfc		<- do.call('rbind', dfc)
		dfc[, COVERAGE:= dfc[, regmatches(FILE,regexpr('C[0-9]+',FILE))]]
		set(dfc, NULL, 'COVERAGE', dfc[,as.numeric(substring(COVERAGE,2))])		
		dfc[, SITE:= dfc[, substr(SEQ, 6,7)]]
		dfci	<- dfc[, list(CALLED= median( (NSEQ-NOCOV)/NSEQ) ), by=c('SITE','COVERAGE')]
		setkey(dfci, SITE, COVERAGE)
		save(dfci, dfc, file='~/git/HPTN071sim/treecomparison/withgaps/150623_HPTN071_CoverageComparison.R')
		ggplot(dfci, aes(x=COVERAGE, y=CALLED*100, colour=SITE, group=SITE)) + 
				geom_line() + geom_point() +
				scale_x_continuous(expand=c(0,0)) +
				scale_y_continuous(breaks=seq(0,100,20),limits=c(0,100),expand=c(0,0)) +
				labs(x='majority coverage\n(minimum number of aligned short read calls that agree)', y='Non-?N calls out of ACGTRWS?N calls\n(median %)') + theme_bw()
		ggsave(file='~/git/HPTN071sim/treecomparison/withgaps/150623_HPTN071_CoverageComparison.pdf', w=6, h=6)
	}
	#
	#	select largest cluster & create partition files for regional 
	#
	if(0)
	{		
		indir.nogaps	<- '/Users/Oliver/git/HPTN071sim/treecomparison/nogaps'
		indir.wgaps		<- '/Users/Oliver/git/HPTN071sim/treecomparison/withgaps_150701'
		infiles.wgaps	<- list.files(indir.wgaps, pattern='^150623_HPTN071_TRAIN.*[0-9]\\.fa$')
		#	create tree for largest cluster		
		infiles.nogapsR	<- list.files(indir.nogaps, pattern='DATEDTREE\\.newick')
		for(infile.nogapsR in infiles.nogapsR)
		{
			#infile.nogapsR	<- infiles.nogapsR[1]
			ph				<- read.tree(paste(indir.nogaps,'/',infile.nogapsR,sep=''))
			load(paste(indir.nogaps, '/', gsub('DATEDTREE.*','SIMULATED_INTERNAL.R',infile.nogapsR), sep=''))
			df.clu			<- subset(df.inds, !is.na(TIME_SEQ))[, list(IDPOP=IDPOP, CLUN=length(IDPOP)), by='IDCLU']
			df.clu			<- subset(df.clu, max(CLUN)==CLUN)
			tmp				<- paste('IDPOP_',df.clu[1,IDPOP],'|',sep='')
			tmp				<- which(sapply(seq_along(ph), function(i)	any(grepl(tmp,ph[[i]]$tip.label,fixed=TRUE))	))
			stopifnot(tmp>0)		
			write.tree(ph[[tmp]],paste(indir.nogaps,'/',gsub('\\.newick','_lrgstclu\\.newick',infile.nogapsR),sep=''))
		}		
		#	create partition files 
		for(infile.wgaps in infiles.wgaps)
		{			
			infile.nogapsR	<- gsub('[^_]*\\.fa','INTERNAL.R',infile.wgaps)
			#	select large transmission cluster
			load(paste(indir.nogaps, '/', infile.nogapsR, sep=''))
			df.clu			<- subset(df.inds, !is.na(TIME_SEQ))[, list(IDPOP=IDPOP, CLUN=length(IDPOP)), by='IDCLU']
			df.clu			<- subset(df.clu, max(CLUN)==CLUN)
			tmp				<- read.dna(file=paste(indir.wgaps, '/', infile.wgaps, sep=''), format='fasta')
			dfs				<- data.table(LABEL= rownames(tmp))
			dfs[, IDPOP:= as.integer(substring(sapply(strsplit(LABEL,'|',fixed=TRUE), '[[', 1),7))]
			dfs				<- merge(df.clu, dfs, by='IDPOP')
			seq		<- tmp[dfs[,LABEL],]
			file	<- paste(indir.wgaps, '/', gsub('\\.','_lrgstclu\\.',infile.wgaps), sep='')		
			write.dna(seq, file=file, format='fasta', colsep='', nbcol=-1)
			save(seq, file=gsub('\\.fa.*','\\.R',file))
			#	determine partition lengths
			#	pol start
			infile.nogaps	<- list.files(indir.nogaps, pattern='pol\\.fa$')
			infile.nogaps	<- infile.nogaps[ grepl(gsub('_SIMULATED.*','',infile.wgaps),infile.nogaps) ]
			stopifnot(length(infile.nogaps)==1)				
			tmp				<- read.dna(paste(indir.nogaps, '/', infile.nogaps, sep=''), format='fasta')[rownames(seq)[1],]
			stopifnot(nrow(tmp)==1)
			key				<- paste(as.character(tmp)[1,1:8], collapse='')
			for(i in seq_len(nrow(seq)))
			{			
				pol.start	<- regexpr(key, paste(as.character(seq[i,]), collapse=''))
				if(pol.start>1)
					break
			}
			stopifnot(pol.start>1)
			#	env start
			infile.nogaps	<- list.files(indir.nogaps, pattern='env\\.fa$')
			infile.nogaps	<- infile.nogaps[ grepl(gsub('_SIMULATED.*','',infile.wgaps),infile.nogaps) ]
			stopifnot(length(infile.nogaps)==1)				
			tmp				<- read.dna(paste(indir.nogaps, '/', infile.nogaps, sep=''), format='fasta')[rownames(seq)[1],]
			stopifnot(nrow(tmp)==1)
			key				<- paste(as.character(tmp)[1,1:8], collapse='')
			for(i in seq_len(nrow(seq)))
			{			
				env.start	<- regexpr(key, paste(as.character(seq[i,]), collapse=''))
				if(env.start>1)
					break
			}
			stopifnot(env.start>1)
			#	write codon partition file
			tmp				<- paste(c(	'DNA, gagcodon1 = 1-',pol.start-1,'\\3\n',
							'DNA, gagcodon2 = 2-',pol.start-1,'\\3\n',
							'DNA, gagcodon3 = 3-',pol.start-1,'\\3\n',
							'DNA, polcodon1 = ',pol.start,'-',env.start-1,'\\3\n',
							'DNA, polcodon2 = ',pol.start+1,'-',env.start-1,'\\3\n',
							'DNA, polcodon3 = ',pol.start+2,'-',env.start-1,'\\3\n',
							'DNA, envcodon1 = ',env.start,'-',ncol(seq),'\\3\n',
							'DNA, envcodon2 = ',env.start+1,'-',ncol(seq),'\\3\n',
							'DNA, envcodon3 = ',env.start+2,'-',ncol(seq),'\\3\n'), collapse='' )
			infile.partition	<- gsub('\\.fa.*','_codon\\.txt',basename(file))
			cat(file=paste(indir.wgaps,'/',infile.partition,sep=''), tmp)
		}		
	}	
	#
	#	create RAxML tree with -fo and codon partitions on mock data
	#	
	if(0)
	{
		#
		#	run ExamML with partition
		#
		indir		<- indir.wgaps
		signat.in	<- "ZAC5_lrgstclu"
		signat.out	<- "ZAC5_lrgstclu"
		infile		<- gsub(paste('_',signat.in,sep=''),'',gsub('\\..*','',basename(file)))
		args.parser	<- paste("-m DNA -q",infile.partition)
		cmd			<- hivc.cmd.examl.bootstrap(indir, infile, signat.in, signat.out, bs.from=0, bs.to=0, prog.bscreate=PR.EXAML.BSCREATE, prog.parser= PR.EXAML.PARSER, prog.starttree= PR.EXAML.STARTTREE, prog.examl=PR.EXAML.EXAML, opt.bootstrap.by="codon", args.parser=args.parser, args.examl="-m GAMMA -D", prog.supportadder=PR.EXAML.BS, tmpdir.prefix="examl")
		
		#	run RAxML
		#	change seed several times to get replicate
		#	compare trees to true tree tmw
		#	set up big trees on cluster for best performing RF
		cmd	<- paste('raxmlHPC-AVX -s ',basename(file),'-m GTRGAMMAI -p42 -n ',gsub('fa','newick',basename(file)), sep='')
		cmd	<- paste('raxmlHPC-AVX -s ',basename(file),' -q ', paste('TMP_',gsub('fa','txt',basename(file)),sep=''),' -m GTRGAMMAI -p42 -n ',gsub('\\.fa','_codon\\.newick',basename(file)), sep='')
		cmd	<- paste('raxmlHPC-AVX -f o -s ',basename(file),' -q ', paste('TMP_',gsub('fa','txt',basename(file)),sep=''),' -m GTRGAMMAI -p42 -n ',gsub('\\.fa','_codon_hcold\\.newick',basename(file)), sep='')
		
		cat(cmd)		
	}
	#	duplicate ExaML files
	if(0)
	{
		indir.wgaps	<- '/Users/Oliver/git/HPTN071sim/treecomparison/withgapstrees'
		infiles		<- data.table(FILE=list.files(indir.wgaps, pattern='\\.R$'))		
		infiles[, {
					file.copy( paste(indir.wgaps,'/',FILE,sep=''), paste(indir.wgaps,'/',gsub('Apr15','Apr15_NOQ',gsub('SIMULATED','SIMULATED_NOQ',FILE)),sep=''))					
				}, by='FILE']
	}
	#	duplicate ExaML files: remove taxa with a high % of '?' 
	if(0)
	{
		indir.wgaps	<- '/Users/Oliver/git/HPTN071sim/treecomparison/withgapstrees'
		infiles		<- data.table(FILE=list.files(indir.wgaps, pattern='\\.R$'))		
		infiles		<- subset(infiles, !grepl('NOQ',FILE) & grepl('BWC5_|UGC5_',FILE))		
		invisible(lapply(infiles[, FILE],function(file)
				{
					load( paste(indir.wgaps,'/',file,sep='') )		
					dc		<- data.table(LABEL=rownames(seq), CVG=seq.length(seq, exclude=c('-','?')) / seq.length(seq, exclude=c('-')))
					#	keep taxa with ACGT >10%, >30%, >60%, >80%
					dc		<- do.call('rbind',lapply(c(0.1,0.3,0.6,0.8), function(cut)
									{
										tmp	<- subset(dc, CVG>=cut)
										tmp[, CUT:=cut]
										tmp
									}))
					#	for each cut, create extra data set
					seqc	<- copy(seq)
					invisible(dc[,{
										seq	<- seqc[LABEL,]	
										cat(paste('\nNumber of taxa in reduced alignment n=',nrow(seq)))
										save(seq, file=paste(indir.wgaps,'/',gsub('SIMULATED',paste('SIMULATED_CUT',CUT*100,sep=''),file),sep='') )
									}, by='CUT'])
					NULL
				}))
		infiles		<- data.table(FILE=list.files(indir.wgaps, pattern='\\.R$'))		
		infiles		<- subset(infiles, grepl('CUT',FILE))
		infiles[, PARTITION:= gsub('\\.R','_codon.txt',FILE)]
		infiles[, {
						file.copy( paste(indir.wgaps,'/',gsub('_CUT[0-9]+','',PARTITION),sep=''), paste(indir.wgaps,'/',PARTITION,sep=''))					
					}, by='FILE']
	}
	#
	#	create ExaML script files
	#		
	if(0)
	{
		indir.wgaps	<- '/Users/Oliver/git/HPTN071sim/treecomparison/withgapstrees'
		infiles		<- data.table(FILE=list.files(indir.wgaps, pattern='\\.R$'))		
		infiles		<- data.table(FILE=list.files(indir.wgaps, pattern='\\.R$'))
		infiles[,SIGNAT:=infiles[, regmatches(FILE,regexpr('BWC.*|UGC.*', FILE))]]
		tmp			<- infiles[, list(BASE=gsub(paste('_',SIGNAT,sep=''),'',FILE)), by='FILE']
		infiles		<- merge(infiles, tmp, by='FILE')
		set(infiles, NULL, 'SIGNAT', infiles[, gsub('\\.R','',SIGNAT)])
		infiles[, PARTITION:= gsub('\\.R','_codon.txt',FILE)]
		bs.from		<- 0
		bs.to		<- 0
		bs.n		<- 1
		#infiles		<- subset(infiles, grepl('NOQ',FILE) & grepl('BWC5_|UGC5_',FILE))
		infiles		<- subset(infiles, grepl('CUT',FILE))
		#infiles		<- subset(infiles, grepl('BWC0',SIGNAT))
		invisible(infiles[, {					
					infile		<- BASE
					signat.in	<- signat.out	<- SIGNAT
					args.parser	<- paste("-m DNA -q",PARTITION)
					args.parser	<- paste("-m DNA")
					args.parser	<- paste("-m DNA")
					args.examl	<- "-f o -m GAMMA"
					#args.examl	<- "-f k -m GAMMA -M"
					cmd			<- hivc.cmd.examl.bootstrap(indir.wgaps, infile, signat.in, signat.out, bs.from=bs.from, bs.to=bs.to, prog.bscreate=PR.EXAML.BSCREATE, prog.parser= PR.EXAML.PARSER, prog.starttree= PR.EXAML.STARTTREE, prog.examl=PR.EXAML.EXAML, opt.bootstrap.by="codon", args.parser=args.parser, args.examl=args.examl, prog.supportadder=PR.EXAML.BS, tmpdir.prefix="examl")					
					invisible(lapply(cmd, function(x)
							{												
								x		<- hivc.cmd.hpcwrapper(x, hpc.walltime=21, hpc.q="pqeelab", hpc.mem="950mb", hpc.nproc=1)
								signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
								outfile	<- paste("ex",signat,sep='.')
								#cat(x)
								hivc.cmd.hpccaller(outdir, outfile, x)
								Sys.sleep(1)
							}))
					NULL					
				}, by='FILE'])				
	}
}
##--------------------------------------------------------------------------------------------------------
##	olli 27.06.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.treecomparison.gaps.RobinsonFould.ByMajorityCoverageThreshold<- function()
{		
	#
	#	compare ExaML trees to true topology
	#
	if(1)
	{
		require(phangorn)
		indir.wgaps		<- '/Users/Oliver/git/HPTN071sim/treecomparison/withgapstrees_150630'
		indir.nogaps	<- '/Users/Oliver/git/HPTN071sim/treecomparison/nogaps_150630'
		files			<- data.table(SIMFILE=list.files(indir.wgaps, pattern='ExaML_result.*finaltree\\.[0-9]+', recursive=TRUE))
		files[, TRUEFILE:= gsub('SIMULATED.*','DATEDTREE_lrgstclu.newick',gsub('ExaML_result\\.','',basename(SIMFILE)))]
		set(files, NULL, 'TRUEFILE', files[, gsub('Vill_99_Apr15.*','Vill_99_Apr15.newick',TRUEFILE)])
		files[, CNFG:= regmatches(basename(SIMFILE), regexpr('BWC[0-9]+|UGC[0-9]+',basename(SIMFILE)))]
		files[, PARTITION:= 'None']
		set(files, files[, which(!grepl('NOQ', basename(SIMFILE)))], 'PARTITION', '9ByGeneCodon')
		files[, EXCLSEQ:= 'None']		
		tmp				<- files[, which(grepl('CUT[0-9]+', basename(SIMFILE)))]
		set(files, tmp, 'EXCLSEQ', files[tmp, regmatches(basename(SIMFILE), regexpr('CUT[0-9]+',basename(SIMFILE)))])
		files[, DT:= regmatches(basename(SIMFILE), regexpr('TRAIN[0-9]+|Vill_99',basename(SIMFILE)))]
		files[, REP:= as.numeric(regmatches(basename(SIMFILE), regexpr('[0-9]+$',basename(SIMFILE))))]
		files[, CVRG:= as.numeric(regmatches(CNFG, regexpr('[0-9]+',CNFG)))]
		files[, SITE:= regmatches(CNFG, regexpr('BW|UG',CNFG))]
		files			<- subset(files, PARTITION=='9ByGeneCodon' & EXCLSEQ=='None')
		#robinson fould distance
		#SIMFILE<- subset(files, SIMFILE=="Vill_99_Apr15_examlout_BWC0/ExaML_result.Vill_99_Apr15_BWC0.finaltree.000")[, SIMFILE]
		#TRUEFILE<- subset(files, SIMFILE=="Vill_99_Apr15_examlout_BWC0/ExaML_result.Vill_99_Apr15_BWC0.finaltree.000")[, TRUEFILE]
		drf				<- files[,{
					#print(SIMFILE)
					stree		<- unroot(read.tree(paste(indir.wgaps,'/',SIMFILE,sep='')))
					otree		<- unroot(read.tree(paste(indir.nogaps,'/',TRUEFILE,sep='')))					
					if(!is.binary.tree(stree))
						stree	<- multi2di(stree)
					#https://groups.google.com/forum/#!topic/raxml/JgvxgknTeqw
					#normalize with 2n-6		
					rf			<- RF.dist(otree, stree, check.labels = TRUE)
					list(RF=rf, NRF=rf/(2*Ntip(otree)-6))					
				}, by='SIMFILE']
		drf				<- merge(files, drf, by='SIMFILE')
		#plot trees in one pdf		
		invisible(subset(files, REP<15)[,{
					strees		<- lapply(SIMFILE, function(x)	read.tree(paste(indir.wgaps,'/',x,sep=''))	)					
					otree		<- read.tree(paste(indir.nogaps,'/',TRUEFILE[1],sep=''))					
					if(length(strees)+1==16)
					{
						cat(paste('\nplotting',DT,CNFG))
						pdf(file=paste(indir.wgaps,'/',gsub('finaltree.[0-9]+','CMPTREES.pdf',SIMFILE[1]),sep=''), w=8, h=8)					
						def.par <- par(no.readonly = TRUE)
						par(mar=c(0.5,0.5,0.5,0.5))					
						layout( matrix(c(1:16), 4, 4, byrow = TRUE) )					
						plot(otree, show.tip.label=FALSE,edge.color="red")
						invisible(sapply(strees, function(x) plot(x, show.tip.label=FALSE)))					
						par(def.par)
						dev.off()						
					}
					NULL
				}, by=c('DT','CNFG')])
		#plot RF distance
		ggplot(drf, aes(x=CVRG, y=NRF*100, colour=SITE)) + 
				geom_point() +
				scale_x_continuous(expand=c(0.01,0.01), breaks=seq(0,200,20)) +
				scale_y_continuous(breaks=seq(0,100,20),limits=c(0,100),expand=c(0,0)) +
				labs(x='majority coverage\n(minimum number of aligned short read calls that agree)', y='normalized RF distance btw est and true trees\n(RF/(2*taxa-6)') +
				facet_grid(~DT) + 
				theme_bw() + theme(panel.margin.x = unit(0.8, "lines"))
		ggsave(file= paste(indir.wgaps,'/150623_ExaML_Codon_RF_pts.pdf',sep=''), w=9, h=4)
		ggplot(drf, aes(x=factor(CVRG), y=NRF*100, colour=SITE)) + 
				geom_point(colour='grey50') +	geom_boxplot() +
				scale_x_discrete(expand=c(0.01,0.01)) +
				scale_y_continuous(breaks=seq(0,100,20),limits=c(0,100),expand=c(0,0)) +
				labs(x='majority coverage\n(minimum number of aligned short read calls that agree)', y='normalized RF distance btw est and true trees\n(RF/(2*taxa-6)') +
				facet_grid(~DT) + 
				theme_bw() + theme(panel.margin.x = unit(0.8, "lines"))
		ggsave(file= paste(indir.wgaps,'/150630_ExaML_Codon_RF_box.pdf',sep=''), w=12, h=4)
		
	}
}
##--------------------------------------------------------------------------------------------------------
##	olli 30.06.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.treecomparison.gaps.RobinsonFould.ByPartitionScheme<- function()
{		
	#
	#	compare ExaML trees to true topology
	#
	require(phangorn)
	indir.wgaps		<- '/Users/Oliver/git/HPTN071sim/treecomparison/withgapstrees_150630'
	indir.nogaps	<- '/Users/Oliver/git/HPTN071sim/treecomparison/nogaps'
	files			<- data.table(SIMFILE=list.files(indir.wgaps, pattern='ExaML_result.*finaltree\\.[0-9]+', recursive=TRUE))
	files[, TRUEFILE:= gsub('SIMULATED.*','DATEDTREE_lrgstclu.newick',gsub('ExaML_result\\.','',basename(SIMFILE)))]
	set(files, NULL, 'TRUEFILE', files[, gsub('Vill_99_Apr15.*','Vill_99_Apr15.newick',TRUEFILE)])
	files[, CNFG:= regmatches(basename(SIMFILE), regexpr('BWC[0-9]+|UGC[0-9]+',basename(SIMFILE)))]
	files[, PARTITION:= 'None']
	set(files, files[, which(!grepl('NOQ', basename(SIMFILE)))], 'PARTITION', '9ByGeneCodon')
	files[, EXCLSEQ:= 'None']		
	tmp				<- files[, which(grepl('CUT[0-9]+', basename(SIMFILE)))]
	set(files, tmp, 'EXCLSEQ', files[tmp, regmatches(basename(SIMFILE), regexpr('CUT[0-9]+',basename(SIMFILE)))])
	files[, DT:= regmatches(basename(SIMFILE), regexpr('TRAIN[0-9]+|Vill_99',basename(SIMFILE)))]
	files[, REP:= as.numeric(regmatches(basename(SIMFILE), regexpr('[0-9]+$',basename(SIMFILE))))]
	files[, CVRG:= as.numeric(regmatches(CNFG, regexpr('[0-9]+',CNFG)))]
	files[, SITE:= regmatches(CNFG, regexpr('BW|UG',CNFG))]
	filess			<- subset(files, (CNFG=='BWC5' | CNFG=='UGC5') & EXCLSEQ=='None' & DT=='TRAIN2') 						
	#robinson fould distance
	drf				<- filess[,{
				#print(SIMFILE)
				stree		<- unroot(read.tree(paste(indir.wgaps,'/',SIMFILE,sep='')))
				otree		<- unroot(read.tree(paste(indir.nogaps,'/',TRUEFILE,sep='')))					
				if(!is.binary.tree(stree))
					stree	<- multi2di(stree)
				#https://groups.google.com/forum/#!topic/raxml/JgvxgknTeqw
				#normalize with 2n-6		
				rf			<- RF.dist(otree, stree, check.labels = TRUE)
				list(RF=rf, NRF=rf/(2*Ntip(otree)-6))					
			}, by='SIMFILE']
	drf				<- merge(filess, drf, by='SIMFILE')	
	#plot RF distance	
	ggplot(drf, aes(x=factor(PARTITION), y=NRF*100, colour=SITE)) + 
			geom_point(colour='grey50') +	geom_boxplot() +
			scale_x_discrete(expand=c(0.01,0.01)) +
			scale_y_continuous(breaks=seq(0,100,20),limits=c(0,100),expand=c(0,0)) +
			labs(x='partition scheme used', y='normalized RF distance btw est and true trees\n(RF/(2*taxa-6)') +
			facet_grid(~DT) + 
			theme_bw() + theme(panel.margin.x = unit(0.8, "lines"))
	ggsave(file= paste(indir.wgaps,'/150630_ExaML_PartitionNoneOrCodon_RF_box.pdf',sep=''), w=5, h=4)			
}
##--------------------------------------------------------------------------------------------------------
##	olli 01.07.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.treecomparison.gaps.simulate<- function()
{		
	indir.simu		<- '/Users/Oliver/git/HPTN071sim/treecomparison/nogaps'
	indir.gap		<-	'~/git/HPTN071sim/treecomparison/PANGEAcov'
	infile.gap		<- '150623_PANGEAGlobal2681_C5.fa'
	outdir			<- '/Users/Oliver/git/HPTN071sim/treecomparison/withgaps_150701'					
	gap.symbol		<- '?'
	gap.seed		<- 42
	
	if(0)
	{
		gap.country		<- 'BW'
		infile.gaps		<- sapply(c(5),function(x) paste('150623_PANGEAGlobal2681_C',x,'.fa',sep=''))
		invisible(lapply(infile.gaps, function(infile.gap)
						{
							outfile.cov		<- regmatches(infile.gap,regexpr('C[0-9]+',basename(infile.gap)))
							infile.simu		<- '150701_Regional_TRAIN1_SIMULATED'
							outfile			<- paste(infile.simu, '_', gap.country, outfile.cov, '.fa', sep='')	
							ans				<- PANGEA.add.gaps.maintain.triplets(indir.simu, indir.gap, infile.simu, infile.gap, gap.country, gap.symbol, gap.seed, outfile=outfile, verbose=1)
							write.dna(ans, file=paste(outdir, outfile, sep='/'), format='fasta', colsep='', nbcol=-1)					
							infile.simu		<- '150701_Regional_TRAIN2_SIMULATED'
							outfile			<- paste(infile.simu, '_', gap.country, outfile.cov, '.fa', sep='')
							ans				<- PANGEA.add.gaps.maintain.triplets(indir.simu, indir.gap, infile.simu, infile.gap, gap.country, gap.symbol, gap.seed, outfile=outfile, verbose=1)
							write.dna(ans, file=paste(outdir, outfile, sep='/'), format='fasta', colsep='', nbcol=-1)
							infile.simu		<- '150701_Regional_TRAIN3_SIMULATED'
							outfile			<- paste(infile.simu, '_', gap.country, outfile.cov, '.fa', sep='')
							ans				<- PANGEA.add.gaps.maintain.triplets(indir.simu, indir.gap, infile.simu, infile.gap, gap.country, gap.symbol, gap.seed, outfile=outfile, verbose=1)
							write.dna(ans, file=paste(outdir, outfile, sep='/'), format='fasta', colsep='', nbcol=-1)					
						}))	
	}
	if(0)
	{
		gap.country		<- 'UG'
		infile.gaps		<- sapply(c(5),function(x) paste('150623_PANGEAGlobal2681_C',x,'.fa',sep=''))
		invisible(lapply(infile.gaps, function(infile.gap)
						{
							outfile.cov		<- regmatches(infile.gap,regexpr('C[0-9]+',basename(infile.gap)))
							infile.simu		<- '150701_Regional_TRAIN1_SIMULATED'
							outfile			<- paste(infile.simu, '_', gap.country, outfile.cov, '.fa', sep='')	
							ans				<- PANGEA.add.gaps.maintain.triplets(indir.simu, indir.gap, infile.simu, infile.gap, gap.country, gap.symbol, gap.seed, outfile=outfile, verbose=1)
							write.dna(ans, file=paste(outdir, outfile, sep='/'), format='fasta', colsep='', nbcol=-1)					
							infile.simu		<- '150701_Regional_TRAIN2_SIMULATED'
							outfile			<- paste(infile.simu, '_', gap.country, outfile.cov, '.fa', sep='')
							ans				<- PANGEA.add.gaps.maintain.triplets(indir.simu, indir.gap, infile.simu, infile.gap, gap.country, gap.symbol, gap.seed, outfile=outfile, verbose=1)
							write.dna(ans, file=paste(outdir, outfile, sep='/'), format='fasta', colsep='', nbcol=-1)
							infile.simu		<- '150701_Regional_TRAIN3_SIMULATED'
							outfile			<- paste(infile.simu, '_', gap.country, outfile.cov, '.fa', sep='')
							ans				<- PANGEA.add.gaps.maintain.triplets(indir.simu, indir.gap, infile.simu, infile.gap, gap.country, gap.symbol, gap.seed, outfile=outfile, verbose=1)
							write.dna(ans, file=paste(outdir, outfile, sep='/'), format='fasta', colsep='', nbcol=-1)
						}))	
	}
	if(0)
	{
		gap.country		<- 'BW'
		infile.gaps		<- sapply(c(5),function(x) paste('150623_PANGEAGlobal2681_C',x,'.fa',sep=''))
		invisible(lapply(infile.gaps, function(infile.gap)
						{
							outfile.cov		<- regmatches(infile.gap,regexpr('C[0-9]+',basename(infile.gap)))
							infile.simu		<- 'Vill_99_Apr15'
							outfile			<- paste(infile.simu, '_', gap.country, outfile.cov, '.fa', sep='')	
							ans				<- PANGEA.add.gaps.maintain.triplets(indir.simu, indir.gap, infile.simu, infile.gap, gap.country, gap.symbol, gap.seed, outfile=outfile, verbose=1)
							write.dna(ans, file=paste(outdir, outfile, sep='/'), format='fasta', colsep='', nbcol=-1)																		
						}))						
		gap.country		<- 'UG'
		infile.gaps		<- sapply(c(5),function(x) paste('150623_PANGEAGlobal2681_C',x,'.fa',sep=''))
		invisible(lapply(infile.gaps, function(infile.gap)
						{
							outfile.cov		<- regmatches(infile.gap,regexpr('C[0-9]+',basename(infile.gap)))
							infile.simu		<- 'Vill_99_Apr15'
							outfile			<- paste(infile.simu, '_', gap.country, outfile.cov, '.fa', sep='')	
							ans				<- PANGEA.add.gaps.maintain.triplets(indir.simu, indir.gap, infile.simu, infile.gap, gap.country, gap.symbol, gap.seed, outfile=outfile, verbose=1)
							write.dna(ans, file=paste(outdir, outfile, sep='/'), format='fasta', colsep='', nbcol=-1)																		
						}))		
	}
	if(1)
	{	
		gap.country		<- 'BW'
		infile.gaps		<- sapply(c(5),function(x) paste('150623_PANGEAGlobal2681_C',x,'.fa',sep=''))
		invisible(lapply(infile.gaps, function(infile.gap)
						{
							outfile.cov		<- regmatches(infile.gap,regexpr('C[0-9]+',basename(infile.gap)))
							infile.simu		<- 'Vill_98_Jul15'
							outfile			<- paste(infile.simu, '_', gap.country, outfile.cov, '.fa', sep='')	
							ans				<- PANGEA.add.gaps.maintain.triplets(indir.simu, indir.gap, infile.simu, infile.gap, gap.country, gap.symbol, gap.seed, outfile=outfile, verbose=1)
							write.dna(ans, file=paste(outdir, outfile, sep='/'), format='fasta', colsep='', nbcol=-1)																		
						}))						
		gap.country		<- 'UG'
		infile.gaps		<- sapply(c(5),function(x) paste('150623_PANGEAGlobal2681_C',x,'.fa',sep=''))
		invisible(lapply(infile.gaps, function(infile.gap)
						{
							outfile.cov		<- regmatches(infile.gap,regexpr('C[0-9]+',basename(infile.gap)))
							infile.simu		<- 'Vill_98_Jul15'
							outfile			<- paste(infile.simu, '_', gap.country, outfile.cov, '.fa', sep='')	
							ans				<- PANGEA.add.gaps.maintain.triplets(indir.simu, indir.gap, infile.simu, infile.gap, gap.country, gap.symbol, gap.seed, outfile=outfile, verbose=1)
							write.dna(ans, file=paste(outdir, outfile, sep='/'), format='fasta', colsep='', nbcol=-1)																		
						}))		
	}
	if(0)
	{
		#	store non-gappy alignments as BWC0 		
		files			<- data.table(FILE=list.files(indir.simu, pattern='TRAIN.*\\.fa.*$'))			
		files[, SIMU:=files[, regmatches(FILE,regexpr('TRAIN[0-9]', FILE))]]
		files[, GENE:=files[, sapply(strsplit(FILE,'_'), '[[', 5)]]
		set(files, NULL, 'GENE', toupper(files[, substr(GENE,1,nchar(GENE)-3)]))
		set(files, NULL, 'GENE', files[, factor(GENE, levels=c('GAG','POL','ENV'), labels=c('GAG','POL','ENV'))])
		setkey(files, SIMU, GENE)
		#	concatenate simu seqs
		files[, {
					gag	<- read.dna( paste(indir.simu, FILE[1], sep='/'), format='fasta' )
					pol	<- read.dna( paste(indir.simu, FILE[2], sep='/'), format='fasta' )
					env	<- read.dna( paste(indir.simu, FILE[3], sep='/'), format='fasta' )
					tmp	<- cbind(gag, pol, env)
					write.dna( tmp, file=paste(outdir, '/', gsub('gag','BWC0',FILE[1]), sep=''), format='fasta', colsep='', nbcol=-1)
					NULL
				}, by='SIMU']				
	}
}	
##--------------------------------------------------------------------------------------------------------
##	olli 03.07.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.treecomparison.gaps.Cut.150703<- function()
{
	#
	#	cut regional data sets to 950 sequences
	#
	indir.nogaps	<- '/Users/Oliver/git/HPTN071sim/treecomparison/nogaps'	
	indir.wgaps		<- '/Users/Oliver/git/HPTN071sim/treecomparison/withgaps_150701'			
	infiles.nogapsR	<- list.files(indir.nogaps, pattern='DATEDTREE|SUBSTTREE\\.newick')
	infiles.nogapsR	<- infiles.nogapsR[grepl('TRAIN',infiles.nogapsR)]
	for(infile.nogapsR in infiles.nogapsR)
	{
		#infile.nogapsR	<- infiles.nogapsR[1]		
		#	create tree for 980 taxa
		ph				<- read.tree(paste(indir.nogaps,'/',infile.nogapsR,sep=''))
		load(paste(indir.nogaps, '/', gsub('DATEDTREE.*|SUBSTTREE.*','SIMULATED_INTERNAL.R',infile.nogapsR), sep=''))
		df.clu			<- subset(df.inds, !is.na(TIME_SEQ))[, list(IDPOP=IDPOP, CLUN=length(IDPOP)), by='IDCLU']
		setkey(df.clu, IDCLU, CLUN)
		tmp				<- subset(unique(df.clu), cumsum(CLUN)<980)
		df.clu			<- merge(df.clu, subset(tmp, select=c(IDCLU)), by='IDCLU')
		dfl				<- data.table(LABEL= ph$tip.label)
		dfl[, ID:= seq_len(nrow(dfl))]
		dfl[, IDPOP:= dfl[, as.integer(substring(sapply(strsplit(LABEL, '|', fixed=1),'[[',1),7))]]
		df.clu			<- merge(dfl, df.clu, by='IDPOP')		
		ph				<- drop.tip(ph, setdiff(ph$tip.label,df.clu[, LABEL]))		
		write.tree(ph,paste(indir.nogaps,'/',gsub('\\.newick','_TAXA980\\.newick',infile.nogapsR),sep=''))
		#	create alignment of 980 taxa 		
		seq				<- read.dna(paste(indir.wgaps, '/', gsub('DATEDTREE.*|SUBSTTREE.*','SIMULATED_BWC5.fa',infile.nogapsR), sep=''), format='fasta')
		seq				<- seq[df.clu[,LABEL],]
		write.dna(seq, file=paste(indir.wgaps, '/', gsub('DATEDTREE.*|SUBSTTREE.*','SIMULATED_BWC5_TAXA980.fa',infile.nogapsR), sep=''), format='fasta', colsep='', nbcol=-1)
		seq				<- read.dna(paste(indir.wgaps, '/', gsub('DATEDTREE.*|SUBSTTREE.*','SIMULATED_UGC5.fa',infile.nogapsR), sep=''), format='fasta')
		seq				<- seq[df.clu[,LABEL],]
		write.dna(seq, file=paste(indir.wgaps, '/', gsub('DATEDTREE.*|SUBSTTREE.*','SIMULATED_UGC5_TAXA980.fa',infile.nogapsR), sep=''), format='fasta', colsep='', nbcol=-1)
		
		#	copy partition files
		file.copy( 	paste(indir.wgaps, '/', gsub('DATEDTREE.*|SUBSTTREE.*','SIMULATED_BWC5_gene.txt',infile.nogapsR), sep=''),
					paste(indir.wgaps, '/', gsub('DATEDTREE.*|SUBSTTREE.*','SIMULATED_BWC5_TAXA980_gene.txt',infile.nogapsR), sep='')	)
		file.copy( 	paste(indir.wgaps, '/', gsub('DATEDTREE.*|SUBSTTREE.*','SIMULATED_BWC5_codon.txt',infile.nogapsR), sep=''),
					paste(indir.wgaps, '/', gsub('DATEDTREE.*|SUBSTTREE.*','SIMULATED_BWC5_TAXA980_codon.txt',infile.nogapsR), sep='')	)
		file.copy( 	paste(indir.wgaps, '/', gsub('DATEDTREE.*|SUBSTTREE.*','SIMULATED_UGC5_gene.txt',infile.nogapsR), sep=''),
					paste(indir.wgaps, '/', gsub('DATEDTREE.*|SUBSTTREE.*','SIMULATED_UGC5_TAXA980_gene.txt',infile.nogapsR), sep='')	)
		file.copy( 	paste(indir.wgaps, '/', gsub('DATEDTREE.*|SUBSTTREE.*','SIMULATED_UGC5_codon.txt',infile.nogapsR), sep=''),
					paste(indir.wgaps, '/', gsub('DATEDTREE.*|SUBSTTREE.*','SIMULATED_UGC5_TAXA980_codon.txt',infile.nogapsR), sep='')	)
	}
	#
	#	create CUTs
	#
	indir.wgaps	<- '/Users/Oliver/git/HPTN071sim/treecomparison/withgapstrees_150703'
	infiles		<- data.table(FILE=list.files(indir.wgaps, pattern='\\.fa$'))		
	infiles		<- subset(infiles, grepl('BWC5|UGC5',FILE))		
	invisible(lapply(infiles[, FILE],function(file)
					{
						print(file)
						seq		<- read.dna( paste(indir.wgaps,'/',file,sep=''), format='fasta' )		
						dc		<- data.table(LABEL=rownames(seq), CVG=seq.length(seq, exclude=c('-','?')) / seq.length(seq, exclude=c('-')))
						#	keep taxa with ACGT >x%
						dc		<- do.call('rbind',lapply(c(0.05, 0.1, 0.2, 0.3, 0.6), function(cut)
										{
											tmp	<- subset(dc, CVG>=cut)
											tmp[, CUT:=cut]
											tmp
										}))
						#	for each cut, create extra data set
						seqc	<- copy(seq)
						invisible(dc[,{
											seq	<- seqc[LABEL,]	
											cat(paste('\nNumber of taxa in reduced alignment n=',nrow(seq)))
											save(seq, file=paste(indir.wgaps,'/',gsub('\\.fa',paste('_CUT',CUT*100,'\\.R',sep=''),file),sep='') )
										}, by='CUT'])
						NULL
					}))
	infiles		<- data.table(FILE=list.files(indir.wgaps, pattern='\\.R$'))		
	infiles		<- subset(infiles, grepl('CUT',FILE))
	infiles[, PARTITION:= gsub('\\.R','_gene.txt',FILE)]
	infiles[, {
				file.copy( paste(indir.wgaps,'/',gsub('_CUT[0-9]+','',PARTITION),sep=''), paste(indir.wgaps,'/',PARTITION,sep=''))					
			}, by='FILE']
}
##--------------------------------------------------------------------------------------------------------
##	olli 02.07.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.treecomparison.gaps.PartitionRegional<- function()
{
	indir.nogaps	<- '/Users/Oliver/git/HPTN071sim/treecomparison/nogaps'
	indir.wgaps		<- '/Users/Oliver/git/HPTN071sim/treecomparison/withgaps_150701'
	infiles.wgaps	<- list.files(indir.wgaps, pattern='^150701_Regional_TRAIN.*[0-9]\\.fa$')
		
	#	create partition files 
	for(infile.wgaps in infiles.wgaps)
	{			
		seq				<- read.dna(file=paste(indir.wgaps, '/', infile.wgaps, sep=''), format='fasta')
		#	determine partition lengths
		#	pol start
		infile.nogaps	<- list.files(indir.nogaps, pattern='pol\\.fa$')
		infile.nogaps	<- infile.nogaps[ grepl(gsub('_SIMULATED.*','',infile.wgaps),infile.nogaps) ]
		stopifnot(length(infile.nogaps)==1)				
		tmp				<- read.dna(paste(indir.nogaps, '/', infile.nogaps, sep=''), format='fasta')[rownames(seq)[1],]
		stopifnot(nrow(tmp)==1)
		key				<- paste(as.character(tmp)[1,1:15], collapse='')
		for(i in seq_len(nrow(seq)))
		{			
			pol.start	<- regexpr(key, paste(as.character(seq[i,]), collapse=''))
			if(pol.start>1)
				break
		}
		stopifnot(pol.start>1)
		#	env start
		infile.nogaps	<- list.files(indir.nogaps, pattern='env\\.fa$')
		infile.nogaps	<- infile.nogaps[ grepl(gsub('_SIMULATED.*','',infile.wgaps),infile.nogaps) ]
		stopifnot(length(infile.nogaps)==1)				
		tmp				<- read.dna(paste(indir.nogaps, '/', infile.nogaps, sep=''), format='fasta')[rownames(seq)[1],]
		stopifnot(nrow(tmp)==1)
		key				<- paste(as.character(tmp)[1,1:15], collapse='')
		for(i in seq_len(nrow(seq)))
		{			
			env.start	<- regexpr(key, paste(as.character(seq[i,]), collapse=''))
			if(env.start>1)
				break
		}
		stopifnot(env.start>1)
		print(c(pol.start, env.start))
		#	write gene partition file
		tmp				<- paste(c(	'DNA, gag = 1-',pol.start-1,'\n',
									'DNA, pol = ',pol.start,'-',env.start-1,'\n',
									'DNA, env = ',env.start,'-',ncol(seq),'\n'), collapse='' )
		infile.partition	<- gsub('\\.fa.*','_gene\\.txt',infile.wgaps)
		cat(file=paste(indir.wgaps,'/',infile.partition,sep=''), tmp)
		#	write codon partition file
		tmp				<- paste(c(	'DNA, gagcodon1 = 1-',pol.start-1,'\\3\n',
									'DNA, gagcodon2 = 2-',pol.start-1,'\\3\n',
									'DNA, gagcodon3 = 3-',pol.start-1,'\\3\n',
									'DNA, polcodon1 = ',pol.start,'-',env.start-1,'\\3\n',
									'DNA, polcodon2 = ',pol.start+1,'-',env.start-1,'\\3\n',
									'DNA, polcodon3 = ',pol.start+2,'-',env.start-1,'\\3\n',
									'DNA, envcodon1 = ',env.start,'-',ncol(seq),'\\3\n',
									'DNA, envcodon2 = ',env.start+1,'-',ncol(seq),'\\3\n',
									'DNA, envcodon3 = ',env.start+2,'-',ncol(seq),'\\3\n'), collapse='' )
		infile.partition	<- gsub('\\.fa.*','_codon\\.txt',infile.wgaps)
		cat(file=paste(indir.wgaps,'/',infile.partition,sep=''), tmp)
	}	
}
##--------------------------------------------------------------------------------------------------------
##	olli 01.07.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.treecomparison.gaps.PartitionVillage<- function()
{		
	indir.nogaps	<- '/Users/Oliver/git/HPTN071sim/treecomparison/nogaps'
	indir.wgaps		<- '/Users/Oliver/git/HPTN071sim/treecomparison/withgaps_150701'
	infiles.wgaps	<- list.files(indir.wgaps, pattern='^Vill.*[0-9]\\.fa$')
	
	for(infile.wgaps in infiles.wgaps)
	{				
		print(infile.wgaps)
		#infile.wgaps	<- infiles.wgaps[1]
		if(grepl('99_Apr15',infile.wgaps))
		{
			pol.key			<- "tttttcagggaaaattt"
			env.key			<- "atgagagtgagagggaa"			
		}
		if(grepl('98_Jul15',infile.wgaps))
		{
			pol.key			<- "ttttttagggaaaatttagcc"#
			env.key			<- "atgagtgtgagaaggat"			
		}		
		seq				<- read.dna(paste(indir.wgaps, '/', infile.wgaps, sep=''), format='fasta')
		save(seq, file=paste(indir.wgaps,'/',gsub('\\.fa.*','\\.R',infile.wgaps),sep=''))
		for(i in seq_len(nrow(seq)))
		{			
			#print(i)
			pol.start	<- regexpr(tolower(pol.key), paste(as.character(seq[i,]), collapse=''))
			if(pol.start>1)
				break
		}
		stopifnot(pol.start>1)
		for(i in seq_len(nrow(seq)))
		{			
			#print(i) 
			env.start		<- regexpr(tolower(env.key), paste(as.character(seq[i,]), collapse=''))
			if(env.start<4000)
				env.start	<- -1
			if(env.start>1)
				break
		}
		stopifnot(env.start>1)
		print(c(pol.start, env.start))
		#	write gene partition file
		tmp				<- paste(c(	'DNA, gag = 1-',pol.start-1,'\n',
									'DNA, pol = ',pol.start,'-',env.start-1,'\n',
									'DNA, env = ',env.start,'-',ncol(seq),'\n'), collapse='' )
		infile.partition	<- gsub('\\.fa.*','_gene\\.txt',infile.wgaps)
		cat(file=paste(indir.wgaps,'/',infile.partition,sep=''), tmp)
		#	write codon partition file
		tmp				<- paste(c(	'DNA, gagcodon1 = 1-',pol.start-1,'\\3\n',
						'DNA, gagcodon2 = 2-',pol.start-1,'\\3\n',
						'DNA, gagcodon3 = 3-',pol.start-1,'\\3\n',
						'DNA, polcodon1 = ',pol.start,'-',env.start-1,'\\3\n',
						'DNA, polcodon2 = ',pol.start+1,'-',env.start-1,'\\3\n',
						'DNA, polcodon3 = ',pol.start+2,'-',env.start-1,'\\3\n',
						'DNA, envcodon1 = ',env.start,'-',ncol(seq),'\\3\n',
						'DNA, envcodon2 = ',env.start+1,'-',ncol(seq),'\\3\n',
						'DNA, envcodon3 = ',env.start+2,'-',ncol(seq),'\\3\n'), collapse='' )
		infile.partition	<- gsub('\\.fa.*','_codon\\.txt',infile.wgaps)
		cat(file=paste(indir.wgaps,'/',infile.partition,sep=''), tmp)
	}
}
##--------------------------------------------------------------------------------------------------------
##	olli 03.07.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.treecomparison.gaps.RobinsonFould.ExaML3014.byMajorityCoverageThreshold<- function()
{		
	#
	#	compare ExaML trees to true topology
	#
	require(phangorn)
	indir.wgaps		<- '/Users/Oliver/git/HPTN071sim/treecomparison/withgapstrees_150701'
	indir.nogaps	<- '/Users/Oliver/git/HPTN071sim/treecomparison/nogaps_150630'
	files			<- data.table(SIMFILE=list.files(indir.wgaps, pattern='ExaML_result.*finaltree\\.[0-9]+', recursive=TRUE))
	files[, TRUEFILE:= gsub('SIMULATED.*','DATEDTREE_lrgstclu.newick',gsub('ExaML_result\\.','',basename(SIMFILE)))]
	set(files, NULL, 'TRUEFILE', files[, gsub('Vill_99_Apr15.*','Vill_99_Apr15.newick',TRUEFILE)])
	files[, CNFG:= regmatches(basename(SIMFILE), regexpr('BWC[0-9]+|UGC[0-9]+',basename(SIMFILE)))]
	files[, PARTITION:= 'None']
	set(files, files[, which(!grepl('NOQ', basename(SIMFILE)))], 'PARTITION', '9ByGeneCodon')
	files[, EXCLSEQ:= 'None']		
	tmp				<- files[, which(grepl('CUT[0-9]+', basename(SIMFILE)))]
	set(files, tmp, 'EXCLSEQ', files[tmp, regmatches(basename(SIMFILE), regexpr('CUT[0-9]+',basename(SIMFILE)))])
	set(files, tmp, 'EXCLSEQ', files[tmp, substr(EXCLSEQ,4,nchar(EXCLSEQ))])
	files[, DT:= regmatches(basename(SIMFILE), regexpr('TRAIN[0-9]+|Vill_99',basename(SIMFILE)))]
	files[, REP:= as.numeric(regmatches(basename(SIMFILE), regexpr('[0-9]+$',basename(SIMFILE))))]
	files[, CVRG:= as.numeric(regmatches(CNFG, regexpr('[0-9]+',CNFG)))]
	files[, SITE:= regmatches(CNFG, regexpr('BW|UG',CNFG))]
	#filess			<- subset(files, EXCLSEQ!='None')	 						
	#robinson fould distance
	#SIMFILE		<- subset(filess, 'BWC5'==CNFG & EXCLSEQ=='80')[, SIMFILE][1]
	#TRUEFILE	<- subset(filess, 'BWC5'==CNFG & EXCLSEQ=='80')[, TRUEFILE][1]
	drf				<- files[,{
				#print(SIMFILE)
				stree		<- unroot(read.tree(paste(indir.wgaps,'/',SIMFILE,sep='')))
				otree		<- unroot(read.tree(paste(indir.nogaps,'/',TRUEFILE,sep='')))					
				if(!is.binary.tree(stree))
					stree	<- multi2di(stree)
				z			<- setdiff(otree$tip.label, stree$tip.label)
				if(length(z))
					otree	<- unroot(drop.tip(otree, z))				
				rf			<- RF.dist(otree, stree, check.labels = TRUE)
				list(RF=rf, NRF=rf/(2*Ntip(otree)-6))					
			}, by='SIMFILE']
	drf				<- merge(files, drf, by='SIMFILE')	
	#plot RF distance	
	ggplot(drf, aes(x=factor(CVRG), y=NRF*100, colour=SITE)) + 
			geom_point(colour='grey50') +	geom_boxplot() +
			scale_x_discrete(expand=c(0.01,0.01)) +
			scale_y_continuous(breaks=seq(0,100,20),limits=c(0,100),expand=c(0,0)) +
			labs(x='majority coverage\n(minimum number of aligned short read calls that agree)', y='normalized RF distance btw est and true trees\n(RF/(2*taxa-6)') +			
			facet_grid(~DT) + 
			theme_bw() + theme(panel.margin.x = unit(0.8, "lines"))
	ggsave(file= paste(indir.wgaps,'/150630_ExaML_ExaML3014_MJCT_box.pdf',sep=''), w=5, h=4)			
}
##--------------------------------------------------------------------------------------------------------
##	olli 30.06.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.treecomparison.gaps.RobinsonFould.ByCut<- function()
{		
	#
	#	compare ExaML trees to true topology
	#
	require(phangorn)
	indir.wgaps		<- '/Users/Oliver/git/HPTN071sim/treecomparison/withgapstrees_150630'
	indir.nogaps	<- '/Users/Oliver/git/HPTN071sim/treecomparison/nogaps'
	files			<- data.table(SIMFILE=list.files(indir.wgaps, pattern='ExaML_result.*finaltree\\.[0-9]+', recursive=TRUE))
	files[, TRUEFILE:= gsub('SIMULATED.*','DATEDTREE_lrgstclu.newick',gsub('ExaML_result\\.','',basename(SIMFILE)))]
	set(files, NULL, 'TRUEFILE', files[, gsub('Vill_99_Apr15.*','Vill_99_Apr15.newick',TRUEFILE)])
	files[, CNFG:= regmatches(basename(SIMFILE), regexpr('BWC[0-9]+|UGC[0-9]+',basename(SIMFILE)))]
	files[, PARTITION:= 'None']
	set(files, files[, which(!grepl('NOQ', basename(SIMFILE)))], 'PARTITION', '9ByGeneCodon')
	files[, EXCLSEQ:= 'None']		
	tmp				<- files[, which(grepl('CUT[0-9]+', basename(SIMFILE)))]
	set(files, tmp, 'EXCLSEQ', files[tmp, regmatches(basename(SIMFILE), regexpr('CUT[0-9]+',basename(SIMFILE)))])
	set(files, tmp, 'EXCLSEQ', files[tmp, substr(EXCLSEQ,4,nchar(EXCLSEQ))])
	files[, DT:= regmatches(basename(SIMFILE), regexpr('TRAIN[0-9]+|Vill_99',basename(SIMFILE)))]
	files[, REP:= as.numeric(regmatches(basename(SIMFILE), regexpr('[0-9]+$',basename(SIMFILE))))]
	files[, CVRG:= as.numeric(regmatches(CNFG, regexpr('[0-9]+',CNFG)))]
	files[, SITE:= regmatches(CNFG, regexpr('BW|UG',CNFG))]
	filess			<- subset(files, EXCLSEQ!='None')	 						
	#robinson fould distance
	#SIMFILE		<- subset(filess, 'BWC5'==CNFG & EXCLSEQ=='80')[, SIMFILE][1]
	#TRUEFILE	<- subset(filess, 'BWC5'==CNFG & EXCLSEQ=='80')[, TRUEFILE][1]
	drf				<- filess[,{
				#print(SIMFILE)
				stree		<- unroot(read.tree(paste(indir.wgaps,'/',SIMFILE,sep='')))
				otree		<- unroot(read.tree(paste(indir.nogaps,'/',TRUEFILE,sep='')))					
				if(!is.binary.tree(stree))
					stree	<- multi2di(stree)
				z			<- setdiff(otree$tip.label, stree$tip.label)
				if(length(z))
					otree	<- unroot(drop.tip(otree, z))				
				#https://groups.google.com/forum/#!topic/raxml/JgvxgknTeqw
				#normalize with 2n-6		
				rf			<- RF.dist(otree, stree, check.labels = TRUE)
				list(RF=rf, NRF=rf/(2*Ntip(otree)-6))					
			}, by='SIMFILE']
	drf				<- merge(filess, drf, by='SIMFILE')	
	#plot RF distance	
	ggplot(drf, aes(x=factor(EXCLSEQ), y=NRF*100, colour=SITE)) + 
			geom_point(colour='grey50') +	geom_boxplot() +
			scale_x_discrete(expand=c(0.01,0.01)) +
			scale_y_continuous(breaks=seq(0,100,20),limits=c(0,100),expand=c(0,0)) +
			labs(x='taxa excluded if proportion of ? sites exceeds threshold\n(threshold value)', y='normalized RF distance btw est and true trees\n(RF/(2*taxa-6)') +
			facet_grid(~DT) + 
			theme_bw() + theme(panel.margin.x = unit(0.8, "lines"))
	ggsave(file= paste(indir.wgaps,'/150630_ExaML_TaxaExcluded_RF_box.pdf',sep=''), w=5, h=4)			
}
##--------------------------------------------------------------------------------------------------------
##	olli 12.02.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Feb2015.final<- function()
{	
	if(0)
	{
		indir			<- '/Users/Oliver/git/HPTN071sim/source/rPANGEAHIVsim/inst/misc'
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1985, yr.end=NA, seed=NA, s.MODEL='Fixed2Prop', report.prop.recent=1.0,
				s.PREV.max.n=NA, s.INTERVENTION.prop=NA, s.INTERVENTION.start=2015, s.INTERVENTION.mul= NA, s.ARCHIVAL.n=50,
				epi.model='HPTN071', epi.dt=1/48, epi.import=NA, root.edge.fixed=0,
				v.N0tau=1, v.r=2.851904, v.T50=-2,
				wher.mu=log(0.00447743)-0.5^2/2, wher.sigma=0.5, bwerm.mu=log(0.002239075)-0.3^2/2, bwerm.sigma=0.3, er.gamma=4,
				dbg.GTRparam=0, dbg.rER=0, index.starttime.mode='fix1970', startseq.mode='one', seqtime.mode='AtDiag')								
		
		# proposed standard run and control simulation
		pipeline.vary	<- data.table(	label=					c('-sq','-ph','-tr20','-s2x','-mFP85','-y3', '-mFP25'),										
										yr.end=					c(2020, 2020,  2020,  2020,   2020,    2018, 2020),
										epi.import=				c(0.05, 0.05,  0.2,   0.05,   0.05,    0.05, 0.05),
										s.PREV.max.n=			c(1600, 1600,  1600,  3200,   1600,    1280, 1600),
										s.INTERVENTION.prop=	c(0.5,  0.5,   0.5,   0.5,    0.85,    0.375, 0.25),
										seed=                   c(42,   5,     7,     11,     13,      17, 42)
										)		
		pipeline.vary	<- pipeline.vary[7,] 				
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='yr.end' ), 'v', as.character(yr.end))
					set(pipeline.args, which( pipeline.args$stat=='epi.import' ), 'v', as.character(epi.import))
					set(pipeline.args, which( pipeline.args$stat=='s.PREV.max.n' ), 'v', as.character(s.PREV.max.n))											
					set(pipeline.args, which( pipeline.args$stat=='s.INTERVENTION.prop' ), 'v', as.character(s.INTERVENTION.prop))
					set(pipeline.args, which( pipeline.args$stat=='s.seed' ), 'v', as.character(seed))
					print(pipeline.args)						
					if(1)
					{
						#	scenario A					
						infile.ind		<- '150129_HPTN071_scA'
						infile.trm		<- '150129_HPTN071_scA'
						tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp150208-A'
						tmpdir			<- paste(tmpdir,label,sep='')
						dir.create(tmpdir, showWarnings=FALSE)																														
						file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
						file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
						file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
						#system(file)	
						stop()
					}
					if(1)
					{
						#	scenario B					
						infile.ind		<- '150129_HPTN071_scB'
						infile.trm		<- '150129_HPTN071_scB'
						tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp150208-B'
						tmpdir			<- paste(tmpdir,label,sep='')
						dir.create(tmpdir, showWarnings=FALSE)																		
						file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
						file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
						file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
						#system(file)													
					}
					if(1)
					{
						#	scenario C					
						infile.ind		<- '150129_HPTN071_scC'
						infile.trm		<- '150129_HPTN071_scC'
						tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp150208-C'
						tmpdir			<- paste(tmpdir,label,sep='')
						dir.create(tmpdir, showWarnings=FALSE)																		
						file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
						file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
						file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
						#system(file)													
					}
					if(1)
					{
						#	scenario D					
						infile.ind		<- '150129_HPTN071_scD'
						infile.trm		<- '150129_HPTN071_scD'
						tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp150208-D'
						tmpdir			<- paste(tmpdir,label,sep='')
						dir.create(tmpdir, showWarnings=FALSE)																		
						file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
						file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
						file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
						#system(file)													
					}
					if(1)
					{
						#	scenario E					
						infile.ind		<- '150129_HPTN071_scE'
						infile.trm		<- '150129_HPTN071_scE'
						tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp150208-E'
						tmpdir			<- paste(tmpdir,label,sep='')
						dir.create(tmpdir, showWarnings=FALSE)																		
						file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
						file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
						file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
						#system(file)													
					}
					if(1)
					{
						#	scenario F					
						infile.ind		<- '150129_HPTN071_scF'
						infile.trm		<- '150129_HPTN071_scF'
						tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp150208-F'
						tmpdir			<- paste(tmpdir,label,sep='')
						dir.create(tmpdir, showWarnings=FALSE)																									
						file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
						file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
						file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
						#system(file)													
					}					
				}, by='label']
		
	}
}
##--------------------------------------------------------------------------------------------------------
##	runs for Manon's clustering analysis
##	olli 31.03.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Apr2015.Manon<- function()
{	
	if(0)
	{
		indir			<- '/Users/Oliver/git/HPTN071sim/source/rPANGEAHIVsim/inst/misc'
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1985, yr.end=2020, seed=42, s.MODEL='Prop2SuggestedSampling', report.prop.recent=1.0,
				s.PREV.max=NA, s.INTERVENTION.start=2015, 
				epi.model='HPTN071', epi.dt=1/48, epi.import=0.2, root.edge.fixed=0,
				v.N0tau=1, v.r=2.851904, v.T50=-2,
				wher.mu=log(0.00447743)-0.5^2/2, wher.sigma=0.5, bwerm.mu=log(0.002239075)-0.3^2/2, bwerm.sigma=0.3, er.gamma=4,
				dbg.GTRparam=0, dbg.rER=0, index.starttime.mode='fix1970', startseq.mode='one', seqtime.mode=NA)								
		
		# proposed standard run and control simulation
		pipeline.vary	<- data.table(	label=					c('-f70s3','-f70s6','-f80s3','-f80s6','-f80s3f80','-f80s6f80','-f80s3f60','-f80s6f60','-f80s3f40','-f80s6f40','-f80s3f20','-f80s6f20'),
										seqtime.mode=			c('AtYear3', 'AtYear6', 'AtYear3', 'AtYear6', 'AtYear3', 'AtYear6', 'AtYear3', 'AtYear6', 'AtYear3', 'AtYear6', 'AtYear3', 'AtYear6'),
										s.PREV.max= 			c(1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.6, 0.6, 0.4, 0.4, 0.2, 0.2),
										index.starttime.mode=	c('fix1970', 'fix1970', 'fix1980', 'fix1980', 'fix1980', 'fix1980', 'fix1980', 'fix1980', 'fix1980', 'fix1980', 'fix1980', 'fix1980')										
										)	
		pipeline.vary	<- pipeline.vary[5:12,]
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='index.starttime.mode' ), 'v', as.character(index.starttime.mode))
					set(pipeline.args, which( pipeline.args$stat=='seqtime.mode' ), 'v', as.character(seqtime.mode))
					set(pipeline.args, which( pipeline.args$stat=='s.PREV.max' ), 'v', as.character(s.PREV.max))
					print(pipeline.args)										
					if(1)
					{
						#	scenario E					
						infile.ind		<- '150129_HPTN071_scE'
						infile.trm		<- '150129_HPTN071_scE'
						outfile.ind		<- '150407_Regional4Manon'
						outfile.trm		<- '150407_Regional4Manon'
						tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp150401-E'
						tmpdir			<- paste(tmpdir,label,sep='')
						dir.create(tmpdir, showWarnings=FALSE)																		
						file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',outfile.ind,label,'_IND.csv',sep=''))
						file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',outfile.trm,label,'_TRM.csv',sep=''))
						file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(outfile.ind,label,'_IND.csv',sep=''), paste(outfile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
						#system(file)													
					}								
				}, by='label']
		
	}
}
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Apr2015.Manon.postprocess<- function()
{
	indir			<- '/Users/Oliver/git/HPTN071sim'		
	files			<- list.files(path=indir, pattern='*INTERNAL.R$', recursive='TRUE')
	files			<- files[ grepl('150401', files)]
	for(file in files)
	{
		load( paste(indir, file, sep='/') )
		df.trms		<- subset(df.trms, select=c(IDREC, IDTR, TIME_TR, SAMPLED_REC, SAMPLED_TR, IDCLU))
		df.inds		<- subset(df.inds, select=c(IDPOP, TIME_TR, GENDER, DOB, DOD, DIAG_T, TIME_SEQ, IDCLU))
		tmp			<- gsub('SIMULATED_INTERNAL','SIM',file)
		save(df.trms, df.inds, df.seq, file= paste(indir, tmp, sep='/'))
	}
}
##--------------------------------------------------------------------------------------------------------
##	check simulated sequences: create ExaML tree and estimate R2
##	olli 08.05.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Feb2015.plotsim<- function()
{
	#	get randomization names
	set.seed(42)
	indir		<- '/Users/Oliver/duke/2014_Gates/methods_comparison_pipeline/FINAL'	
	dfi			<- data.table(FILE=list.files(indir, '.*zip$', full.names=FALSE))
	
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
	
	dfi[, GSUB_FROM:= sapply(strsplit(FILE, '_'),'[[',3)]
	dfi[, GSUB_TO:= paste(OBJECTIVE,'_sc',SC_RND,sep='')]
	
	#	read df.epi for true simulations 
	indir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/PANGEAHIVsim_internal/freeze_Jan15/regional/150208'
	infiles	<- data.table(FILE=list.files(indir, '.*INTERNAL.R$', full.names=FALSE, recursive=TRUE))
	infiles[, BASENAME:= basename(FILE)]
	infiles[, SC:= sapply(strsplit(BASENAME, '_'),'[[',3)]
	infiles[, CONFIG:= sapply(strsplit(SC, '-'),'[[',2)]
	set(infiles, NULL, 'SC', infiles[, sapply(strsplit(SC, '-'),'[[',1)])	
	infiles	<- merge(infiles, subset(dfi, select=c(SC, CONFIG, DATAT, SC_RND, GSUB_FROM, GSUB_TO)), by=c('SC','CONFIG'))
	set(infiles, NULL, 'BASENAME', NULL)	
	infiles	<- subset(infiles, CONFIG=='sq' | SC%in%c('scE','scF'))
	df.epi	<- do.call('rbind',lapply(seq_len(nrow(infiles)), function(i)
					{
						file	<- paste(indir, '/', infiles[i, FILE], sep='')
						cat(paste('\nread file', infiles[i, FILE]))
						load(file)						
						epi.adult	<- 13
						tp			<- 1/12
						suppressWarnings( df.trms[, MNTH:= df.trms[, floor(TIME_TR) + floor( (TIME_TR%%1)*100 %/% (tp*100) ) * tp]] )						
						df.epi		<- df.trms[, list(INC=length(IDREC), INC_ACUTE=length(which(TR_ACUTE=='Yes')),IMPORT=length(which(IDTR<0))), by='MNTH']
						setkey(df.epi, MNTH)
						tmp			<- df.epi[, 	{
									sexactive		<- which( (df.inds[['DOB']]+epi.adult-tp/2)<=MNTH  &  (df.inds[['DOD']]+tp/2)>MNTH )
									infected.ever	<- which( (df.inds[['TIME_TR']]-tp/2)<=MNTH )
									infected		<- which( (df.inds[['DOB']]-tp/2)<=MNTH  &  (df.inds[['DOD']]-tp/2)>MNTH  &  (df.inds[['TIME_TR']]-tp/2)<=MNTH )
									diag			<- which( (df.inds[['DOB']]-tp/2)<=MNTH  &  (df.inds[['DOD']]-tp/2)>MNTH  &  (df.inds[['DIAG_T']]-tp/2)<=MNTH )
									diag.new		<- which( (df.inds[['DIAG_T']]-tp/2)==MNTH )
									treated			<- which( (df.inds[['DOB']]-tp/2)<=MNTH  &  (df.inds[['DOD']]-tp/2)>MNTH  &  (df.inds[['ART1_T']]-tp/2)<=MNTH & (is.na(df.inds[['VLS1_TE']]) | (df.inds[['VLS1_TE']]-tp/2)>MNTH) )
									infdead			<- which( (df.inds[['DOD']]-tp/2)==MNTH  &  (df.inds[['TIME_TR']]-tp/2)<=MNTH )									
									sampled			<- which( (df.inds[['DOB']]-tp/2)<=MNTH  &  (df.inds[['DOD']]-tp/2)>MNTH  &  (df.inds[['TIME_SEQ']]-tp/2)<=MNTH )
									list(POP=length(sexactive), PREV=length(infected), PREV_EVER=length(infected.ever), PREVDIED=length(infdead), DIAG=length(diag), NEW_DIAG=length(diag.new), TREATED=length(treated), SEQ=length(sampled))				
								},by='MNTH']
						df.epi		<- merge( tmp, df.epi, by='MNTH' )		
						set(df.epi, NULL, 'PREVp', df.epi[, PREV/(POP-PREV)])	
						set(df.epi, NULL, 'SEQp', df.epi[, SEQ/PREV])
						set(df.epi, NULL, 'ARTcov', df.epi[, TREATED/PREV])
						set(df.epi, NULL, 'UNDIAGp', df.epi[, (PREV-DIAG)/PREV])
						set(df.epi, NULL, 'GROWTHr', c(NA_real_, df.epi[, diff(log(PREV))]))
						df.epi[, DUMMY:= seq_len(nrow(df.epi))]
						tmp			<- df.epi[, {
									tmp						<- seq.int(DUMMY-6,DUMMY+5)
									tmp[tmp<1]				<- 1
									tmp[tmp>nrow(df.epi)]	<- nrow(df.epi) 
									list(	INCp= sum(df.epi[['INC']][tmp]) / mean(df.epi[['POP']][tmp]-df.epi[['PREV']][tmp]),
											IMPORTp= sum(df.epi[['IMPORT']][tmp]) / sum(df.epi[['INC']][tmp]),
											ACUTEp= sum(df.epi[['INC_ACUTE']][tmp]) / sum(df.epi[['INC']][tmp])
											)	
									}, by='DUMMY']
						df.epi	<- merge(df.epi, tmp, by='DUMMY')
						df.epi[, SC:= infiles[i, substr(SC,3,3)]] 
						df.epi[, DUMMY:=NULL]
						df.epi			
					}))
	#	plot prevalence				
	df.epi	<- melt(df.epi, id.vars=c('SC', 'MNTH'))	
	df		<- subset(df.epi, variable%in%c('POP','PREV','INC','PREVp','INCp','SEQp','ARTcov','UNDIAGp','ACUTEp','IMPORTp') & MNTH>=1985)
	tmp		<- data.table(	variable=c('POP','PREV','INC','PREVp','INCp','SEQp','ARTcov','UNDIAGp','ACUTEp','IMPORTp'),
							legend=c('adult\npopulation', 'HIV infected', 'new cases\n(month)','Prevalence\n(%)','%Incidence\n(year)','Sequence\ncoverage\n(%)','ART\ncoverage\n(%)','Undiagnosed\n(%)','% transm\nfrom < 3m\n(yr)','% transm\nfrom outside\n(yr)'))
	set(tmp, NULL, 'legend', tmp[, factor(legend, levels=tmp$legend, labels=tmp$legend)])									
	df		<- merge(df, tmp, by='variable')
	tmp		<- df[, which(!variable%in%c('POP','PREV','INC'))]
	set(df, tmp,'value', df[tmp, value*100])
	ggplot(df, aes(x=MNTH, y=value, group=SC, colour=SC)) + geom_line() +
			scale_colour_brewer(palette='Set1') +
			facet_wrap(~legend, scales='free',ncol=3) +
			theme_bw() + labs(x='', y='',colour='epidemic scenario', title='Regional\n(since 1985)\n') +
			theme(legend.position='bottom')
	outdir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim_internal/documents/external/2015_05_results/results'
	ggsave(file=paste(outdir,'/regional_sincestart.pdf',sep=''), w=10, h=10)
	
	df		<- subset(df.epi, variable%in%c('POP','PREV','INC','PREVp','INCp','SEQp','ARTcov','UNDIAGp','ACUTEp','IMPORTp') & MNTH>=2014.5)
	tmp		<- data.table(	variable=c('POP','PREV','INC','PREVp','INCp','SEQp','ARTcov','UNDIAGp','ACUTEp','IMPORTp'),
			legend=c('adult\npopulation', 'HIV infected', 'new cases\n(month)','Prevalence\n(%)','%Incidence\n(year)','Sequence\ncoverage\n(%)','ART\ncoverage\n(%)','Undiagnosed\n(%)','% transm\nfrom < 3m\n(yr)','% transm\nfrom outside\n(yr)'))
	set(tmp, NULL, 'legend', tmp[, factor(legend, levels=tmp$legend, labels=tmp$legend)])									
	df		<- merge(df, tmp, by='variable')
	tmp		<- df[, which(!variable%in%c('POP','PREV','INC'))]
	set(df, tmp,'value', df[tmp, value*100])
	ggplot(df, aes(group=SC, colour=SC)) +
			geom_rect(data=NULL, aes(xmin=2015, xmax=2018, ymin=-Inf, ymax=Inf), fill='grey80', colour='transparent', alpha=0.2) +
			geom_line(aes(x=MNTH, y=value)) +
			scale_colour_brewer(palette='Set1') +
			facet_wrap(~legend, scales='free',ncol=3) +
			theme_bw() + labs(x='', y='',colour='epidemic scenario', title='Regional\n(since start of intervention)\n') +
			theme(legend.position='bottom')
	outdir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim_internal/documents/external/2015_05_results/results'
	ggsave(file=paste(outdir,'/regional_sinceintervention.pdf',sep=''), w=10, h=10)
	
}	
##--------------------------------------------------------------------------------------------------------
##	check simulated sequences: create ExaML tree and estimate R2
##	olli 08.05.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Feb2015.resultsforsim<- function()
{
	#	get randomization names
	set.seed(42)
	indir		<- '/Users/Oliver/duke/2014_Gates/methods_comparison_pipeline/FINAL'	
	dfi			<- data.table(FILE=list.files(indir, '.*zip$', full.names=FALSE))
	
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
	
	dfi[, GSUB_FROM:= sapply(strsplit(FILE, '_'),'[[',3)]
	dfi[, GSUB_TO:= paste(OBJECTIVE,'_sc',SC_RND,sep='')]
	
	#	read df.epi for true simulations 
	indir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/PANGEAHIVsim_internal/freeze_Jan15/regional/150208'
	infiles	<- data.table(FILE=list.files(indir, '.*INTERNAL.R$', full.names=FALSE, recursive=TRUE))
	infiles[, BASENAME:= basename(FILE)]
	infiles[, SC:= sapply(strsplit(BASENAME, '_'),'[[',3)]
	infiles[, CONFIG:= sapply(strsplit(SC, '-'),'[[',2)]
	set(infiles, NULL, 'SC', infiles[, sapply(strsplit(SC, '-'),'[[',1)])	
	infiles	<- merge(infiles, subset(dfi, select=c(SC, CONFIG, DATAT, SC_RND, GSUB_FROM, GSUB_TO)), by=c('SC','CONFIG'))
	set(infiles, NULL, 'BASENAME', NULL)	
	df.epi	<- do.call('rbind',lapply(seq_len(nrow(infiles)), function(i)
			{
				file	<- paste(indir, '/', infiles[i, FILE], sep='')
				cat(paste('\nread file', infiles[i, FILE]))
				load(file)
				tmp		<- cbind( df.epi, infiles[i, -which(grepl('FILE', names(infiles))), with=0] )
				tmp			
			}))
	cat(paste('\nread data for scenarios, n=',df.epi[, length(unique(SC_RND))]))
	
	# 	get objectives for true simulations
	#	obj i: A-D decreasing, E-F stable
	#subset(df.epi, YR==2019)
	tmp		<- data.table(	SC=paste('sc',c('A','B','C','D','E','F'),sep=''), OBJ_i= c('decreasing','decreasing','decreasing','decreasing','stable','stable'))
	dfo		<- merge(df.epi, tmp, by='SC')
	#	obj ii: %incidence in last year of eval period
	tmp		<- unique(subset(df.epi, select=CONFIG)) 
	tmp[, YR:=2018]
	set(tmp, tmp[, which(CONFIG=='y3')], 'YR', 2016)
	tmp		<- subset(merge(df.epi, tmp, by=c('CONFIG','YR')), select=c(SC, CONFIG, INCp, YR))
	setkey(tmp, SC)
	setnames(tmp, c('INCp','YR'), c('OBJ_ii','OBJ_ii_te'))
	dfo		<- merge(dfo, tmp, by=c('SC','CONFIG'))
	#	obj iii: ratio incidence
	tmp		<- merge(tmp, subset(df.epi, YR==2014, select=c(SC, CONFIG, INCp)), by=c('SC','CONFIG'))
	set(tmp, NULL, 'OBJ_iii', tmp[, OBJ_ii/INCp] )
	dfo		<- merge(dfo, subset(tmp, select=c(SC, CONFIG, OBJ_iii)), by=c('SC','CONFIG'))
	#	obj iv: %acute 
	#subset(df.epi, YR==2014)
	tmp		<- data.table(	SC=paste('sc',c('A','B','C','D','E','F'),sep=''), OBJ_iv= c('<10%','>30%','<10%','>30%','<10%','>30%'))
	dfo		<- merge(dfo, tmp, by='SC')
	#	obj v: %acute at baseline
	tmp		<- subset(df.epi, YR==2014, select=c(SC, CONFIG, ACUTEp))
	setkey(tmp, SC)
	setnames(tmp, 'ACUTEp', 'OBJ_v')
	dfo		<- merge(dfo, tmp, by=c('SC','CONFIG'))
	#	obj vi: %acute at last yr of eval period
	tmp		<- unique(subset(df.epi, select=CONFIG)) 
	tmp[, YR:=2018]
	set(tmp, tmp[, which(CONFIG=='y3')], 'YR', 2016)
	tmp		<- subset(merge(df.epi, tmp, by=c('CONFIG','YR')), select=c(SC, CONFIG, ACUTEp, YR))
	setkey(tmp, SC)
	setnames(tmp, c('ACUTEp','YR'),c('OBJ_vi','OBJ_vi_te'))
	dfo		<- merge(dfo, tmp, by=c('SC','CONFIG'))
	
	#	write csv file
	setkey(dfo, SC, CONFIG)
	ans		<- unique(dfo)
	ans		<- subset(ans, select=c(SC, CONFIG, DATAT, SC_RND, GSUB_FROM, GSUB_TO, OBJ_i, OBJ_ii, OBJ_ii_te, OBJ_iii, OBJ_iv, OBJ_v, OBJ_vi, OBJ_vi_te))
	ans[, TEAM:='True']
	ans[, SUBMISSION_DATE:='08.05.2015']
	ans[, SIM_SCENARIO:= paste('150129_PANGEAsim_Regional_',GSUB_TO,'_SIMULATED_',DATAT,sep='')]
	ans[, USED_GENES:='all']
	ans[, ESTIMATE:='central']
	
	outdir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/PANGEAHIVsim_internal/documents/external/2015_05_results'
	file	<- paste(outdir, '/answers_Regional_Feb2015_rFormat.csv', sep='')
	write.csv(subset(ans, select=c(TEAM, SUBMISSION_DATE, SIM_SCENARIO, USED_GENES, OBJ_i, OBJ_ii, OBJ_iii, OBJ_iv, OBJ_v, OBJ_vi, ESTIMATE)), file=file, row.names=FALSE)
	
	setnames(ans, 'OBJ_ii_te','te')
	ans[, tb:=2014]
	file	<- paste(outdir, '/Regional_Feb2015_tb_te.csv', sep='')
	write.csv(subset(ans, select=c(SIM_SCENARIO, tb, te)), file=file, row.names=FALSE)
	
	
}
##--------------------------------------------------------------------------------------------------------
##	evaluate results
##	olli 08.05.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Feb2015.evaluate<- function()
{
	require(RColorBrewer)
	dfa		<- project.PANGEA.TEST.pipeline.Feb2015.evaluate.read()
	outdir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim_internal/documents/external/2015_05_results/results'
	save(dfa, file=paste(outdir,'/submissions.R',sep=''))
	load(paste(outdir,'/submissions.R',sep=''))
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
	set(dfa, NULL, 'IMPRT', dfa[, factor(IMPRT*100, levels=c(5,20,2,0), labels=paste(c(5,20,2,0),'%',sep=''))])
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
	#	add team color		
	TEAM_CL	<- brewer.pal(dfa[,length(unique(TEAM))], 'Paired')
	names(TEAM_CL)	<- dfa[, unique(TEAM)]
	TEAM_CL[7]		<- "#386CB0"
	TEAM_CL[3]		<- "#FF7F00"
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
	#
	#	show village / regional: Imports, Total Sequence, Sequence Coverage
	#
	
	
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
	
	#	for each objective
	#	compare results across teams
	require("grid")
	#	compare objectives with / without seq data, village + regional
	df	<- subset(dfa, DATA_T=='seq')	
	#	regional, trees corresponding to seq data sets
	tmp	<- subset(dfa, DATA_T=='phy' & SMPL_D=='5' & SMPL_M=='much' & SMPL_C=='8%' & SMPL_N==1600 & IMPRT=='5%')
	df	<- rbind(df, tmp)
	#	village, trees corresponding to seq data sets
	tmp	<- subset(dfa, SC_RND%in%c('11','09','12','10','00'))
	df	<- rbind(df, tmp)
	tmp	<- subset(df, OBJ%in%c('OBJ_ii','OBJ_iii','OBJ_v','OBJ_vi') & !grepl('(',TEAM,fixed=1) & DATAT_L=='Regional' & USED_GENES=='all')
	tmp2<- as.data.table(expand.grid(central=c(0.01,0.55), AC_L=c('%Acute\nlow','%Acute\nhigh'), INT_L=c('ART scale up\nfast','ART scale up\nslow','ART scale up\nnone'), TEAM='dummy', DATA_T2='using\ntrue tree', OBJ_L2=c('%Acute\n(baseline)','%Acute\n(endpoint)')))
	set(tmp2, NULL, 'INT_L', tmp2[, factor(INT_L, levels=c('ART scale up\nfast','ART scale up\nslow','ART scale up\nnone'), labels=c('ART scale up\nfast','ART scale up\nslow','ART scale up\nnone'))])
	setnames(tmp2, 'TEAM','team')
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')		
	ggplot(subset(tmp, TEAM!='True'), aes(y=INT_L, x=central)) +
			geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) + 
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM+DATA_T2~OBJ_L2+AC_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.5, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=2)) +
			labs(x= '\nEstimates', y='', title='Primary objectives: quantitative\n(Regional)\n')	
	ggsave(file=paste(outdir,'/res_acrossTEAM_primaryregional_nmbrs','.pdf',sep=''), w=13, h=10)
	#	qualitative
	tmp		<- subset(df, OBJ%in%c('OBJ_i','OBJ_iv') & !grepl('(',TEAM,fixed=1) & DATAT_L=='Regional' & USED_GENES=='all')
	set(tmp, NULL, 'central', tmp[, as.character(central)])
	set(tmp, tmp[, which(OBJ=='OBJ_i' & central=='-1')], 'central', 'declining')
	set(tmp, tmp[, which(OBJ=='OBJ_i' & central=='0')], 'central', 'stable')
	set(tmp, tmp[, which(OBJ=='OBJ_i' & central=='1')], 'central', 'increasing')
	set(tmp, tmp[, which(OBJ=='OBJ_iv' & central=='-1')], 'central', '<15%')
	set(tmp, tmp[, which(OBJ=='OBJ_iv' & central=='0')], 'central', '15%-30%')
	set(tmp, tmp[, which(OBJ=='OBJ_iv' & central=='1')], 'central', '>30%')
	tmp2	<- c('declining','stable','increasing','<15%','15%-30%','>30%')
	set(tmp, NULL, 'central', tmp[, factor(central, levels=tmp2, labels=tmp2)])	
	tmp2	<- as.data.table(expand.grid(central=tmp2, AC_T=c('low','high'), INT_L=c('ART scale up\nfast','ART scale up\nslow','ART scale up\nnone'), TEAM='dummy', DATA_T2='using\ntrue tree', OBJ_L2='Incidence\nTrend'))
	set(tmp2, tmp2[, which(grepl('%', central))], 'OBJ_L2', '%Acute Ctgr\n(baseline)')
	set(tmp2, NULL, 'INT_L', tmp2[, factor(INT_L, levels=c('ART scale up\nfast','ART scale up\nslow','ART scale up\nnone'), labels=c('ART scale up\nfast','ART scale up\nslow','ART scale up\nnone'))])	
	setkey(tmp2, central)
	ggplot(tmp2, aes(y=INT_L, x=central, colour=TEAM)) +
			#geom_errorbarh(height=0.3) + 
			geom_point(data=tmp2, size=1, colour='transparent') +
			geom_jitter(data=subset(tmp, TEAM!='True'), size=3, pch=13, position = position_jitter(width=0, height=.15)) +
			geom_point(data=subset(tmp, TEAM=='True'), size=3, colour='black', pch=18) +			
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(DATA_T2~OBJ_L2+AC_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.5, "lines"), legend.position='bottom', panel.grid.minor= element_blank(), panel.grid.major= element_blank()) +
			guides(colour=guide_legend(ncol=2)) +
			labs(x= '\nEstimates', y='', title='Primary objectives: qualitatitve\n(Regional)\n')	
	ggsave(file=paste(outdir,'/res_acrossTEAM_primaryregional_qultv','.pdf',sep=''), w=13, h=6)
	#invisible(lapply(c('OBJ_ii','OBJ_iii','OBJ_v','OBJ_vi'), function(x)
	#		{				
	#			#x	<- 'OBJ_ii'
	#			tmp	<- subset(df, OBJ==x & !grepl('(',TEAM,fixed=1) & DATAT_L=='Regional' & USED_GENES=='all')
	#			if(x=='OBJ_i')
	#				set(tmp, NULL, 'central', tmp[, factor(central, levels=c(1,0,-1), labels=c('increasing','stable','decreasing'))])
	#			if(x=='OBJ_iv')
	#				set(tmp, NULL, 'central', tmp[, factor(central, levels=c(1,0,-1), labels=c('>30%','15%-30%','<15%'))])	
	#			ggplot(subset(tmp, TEAM!='True'), aes(y=DATA_T, x=central, xmin=lower95, xmax=upper95, colour=TEAM)) +
	#					geom_errorbarh(height=0.3) + geom_point(size=3, pch=13) +
	#					geom_point(data=subset(tmp, TEAM=='True'), size=3, colour='black', pch=18) +
	#					scale_colour_manual(values=TEAM_CL) +
	#					facet_grid(AC_L~INT_L, scales='free', space='free_y') +
	#					theme_bw() + theme(panel.margin.x= unit(2, "lines")) + 
	#					labs(x= paste('\n',gsub('\n',' ',tmp[1, OBJ_L]),sep=''), y='')
	#			ggsave(file=paste(outdir,'/res_acrossTEAM_primaryregional_',x,'.pdf',sep=''), w=10, h=3)
	#		}))
	#	%INCIDENCE
	tmp	<- subset(df, OBJ%in%c('OBJ_ii') & !grepl('(',TEAM,fixed=1) & DATAT_L=='Village' & USED_GENES=='all')
	set(tmp, NULL, 'OBJ_L2', tmp[, factor(as.character(OBJ_L2))])
	tmp2<- as.data.table(expand.grid(central=c(0.01,0.03), AC_L=c('%Acute\nlow','%Acute\nhigh'), INT_L=c('ART scale up\nfast','ART scale up\nslow','ART scale up\nnone'), TEAM='dummy', DATA_T2='using\ntrue tree', OBJ_L2=c('%Incidence')))
	set(tmp2, NULL, 'INT_L', tmp2[, factor(INT_L, levels=c('ART scale up\nfast','ART scale up\nslow','ART scale up\nnone'), labels=c('ART scale up\nfast','ART scale up\nslow','ART scale up\nnone'))])
	setnames(tmp2, 'TEAM','team')
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	ggplot(subset(tmp, TEAM!='True'), aes(y=INT_L, x=central*100, colour=TEAM)) +
			geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95*100, xmax=upper95*100), height=0.3) + geom_point(size=3, pch=13) +
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(AC_L+DATA_T2~TEAM, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.5, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=2)) +
			labs(x= '%', y='', title='Primary objectives: %Incidence\n(Village)\n')	
	ggsave(file=paste(outdir,'/res_acrossTEAM_primaryvillage_pcinc','.pdf',sep=''), w=13, h=8)
	ggplot(subset(tmp, TEAM!='True'), aes(y=INT_L, x=central*100, colour=TEAM)) +
			geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95*100, xmax=upper95*100), height=0.3) + geom_point(size=3, pch=13) +
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(AC_L+DATA_T2~TEAM, scales='free_y', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.5, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=2)) +
			labs(x= '%', y='', title='Primary objectives: %Incidence\n(Village)\n')	
	ggsave(file=paste(outdir,'/res_acrossTEAM_primaryvillage_pcinc_ssc','.pdf',sep=''), w=13, h=8)	
	#	REDUCTION INCIDENCE
	tmp	<- subset(df, OBJ%in%c('OBJ_iii') & !grepl('(',TEAM,fixed=1) & DATAT_L=='Village' & USED_GENES=='all')
	set(tmp, NULL, 'OBJ_L2', tmp[, factor(as.character(OBJ_L2))])
	tmp2<- as.data.table(expand.grid(central=c(0.01,0.03), AC_L=c('%Acute\nlow','%Acute\nhigh'), INT_L=c('ART scale up\nfast','ART scale up\nslow','ART scale up\nnone'), TEAM='dummy', DATA_T2='using\ntrue tree', OBJ_L2=c('Incidence\nreduction')))
	set(tmp2, NULL, 'INT_L', tmp2[, factor(INT_L, levels=c('ART scale up\nfast','ART scale up\nslow','ART scale up\nnone'), labels=c('ART scale up\nfast','ART scale up\nslow','ART scale up\nnone'))])
	setnames(tmp2, 'TEAM','team')
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	ggplot(subset(tmp, TEAM!='True'), aes(y=INT_L, x=central, colour=TEAM)) +
			geom_point(data=tmp2, size=1, colour='transparent') +
			geom_vline(xintercept=1, colour='grey50', lwd=1) +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95), height=0.3) + geom_point(size=3, pch=13) +
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(AC_L+DATA_T2~TEAM, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.5, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=2)) +
			labs(x= '%', y='', title='Primary objectives: Incidence Reduction\n(Village)\n')	
	ggsave(file=paste(outdir,'/res_acrossTEAM_primaryvillage_redinc','.pdf',sep=''), w=13, h=8)
	ggplot(subset(tmp, TEAM!='True'), aes(y=INT_L, x=central, colour=TEAM)) +
			geom_point(data=tmp2, size=1, colour='transparent') +
			geom_vline(xintercept=1, colour='grey50', lwd=1) +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95), height=0.3) + geom_point(size=3, pch=13) +
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(AC_L+DATA_T2~TEAM, scales='free_y', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.5, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=2)) +
			labs(x= '%', y='', title='Primary objectives: Incidence Reduction\n(Village)\n')	
	ggsave(file=paste(outdir,'/res_acrossTEAM_primaryvillage_redinc_ssc','.pdf',sep=''), w=13, h=8)	
	#	%ACUTE BASELINE
	tmp	<- subset(df, OBJ%in%c('OBJ_v') & !grepl('(',TEAM,fixed=1) & DATAT_L=='Village' & USED_GENES=='all')
	set(tmp, NULL, 'OBJ_L2', tmp[, factor(as.character(OBJ_L2))])
	tmp2<- as.data.table(expand.grid(central=c(0.01,0.03), AC_L=c('%Acute\nlow','%Acute\nhigh'), INT_L=c('ART scale up\nfast','ART scale up\nslow','ART scale up\nnone'), TEAM='dummy', DATA_T2='using\ntrue tree', OBJ_L2=c('%Acute\n(baseline)')))
	set(tmp2, NULL, 'INT_L', tmp2[, factor(INT_L, levels=c('ART scale up\nfast','ART scale up\nslow','ART scale up\nnone'), labels=c('ART scale up\nfast','ART scale up\nslow','ART scale up\nnone'))])
	setnames(tmp2, 'TEAM','team')
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	ggplot(subset(tmp, TEAM!='True'), aes(y=INT_L, x=central*100, colour=TEAM)) +
			geom_point(data=tmp2, size=1, colour='transparent') +			
			geom_errorbarh(aes(xmin=lower95*100, xmax=upper95*100), height=0.3) + geom_point(size=3, pch=13) +
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(AC_L+DATA_T2~TEAM, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=2)) +
			labs(x= '%', y='', title='Primary objectives: %Acute at baseline\n(Village)\n')	
	ggsave(file=paste(outdir,'/res_acrossTEAM_primaryvillage_pcacutebaseline','.pdf',sep=''), w=8, h=7)
	ggplot(subset(tmp, TEAM!='True'), aes(y=INT_L, x=central*100, colour=TEAM)) +
			geom_point(data=tmp2, size=1, colour='transparent') +			
			geom_errorbarh(aes(xmin=lower95*100, xmax=upper95*100), height=0.3) + geom_point(size=3, pch=13) +
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			coord_cartesian(xlim=c(0,50)) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(AC_L+DATA_T2~TEAM, scales='free_y', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=2)) +
			labs(x= '%', y='', title='Primary objectives: %Acute at baseline\n(Village)\n')	
	ggsave(file=paste(outdir,'/res_acrossTEAM_primaryvillage_pcacutebaseline_ssc','.pdf',sep=''), w=10, h=7)
	#	%ACUTE ENDPOINT
	tmp	<- subset(df, OBJ%in%c('OBJ_vi') & !grepl('(',TEAM,fixed=1) & DATAT_L=='Village' & USED_GENES=='all')
	set(tmp, NULL, 'OBJ_L2', tmp[, factor(as.character(OBJ_L2))])
	tmp2<- as.data.table(expand.grid(central=c(0.01,0.03), AC_L=c('%Acute\nlow','%Acute\nhigh'), INT_L=c('ART scale up\nfast','ART scale up\nslow','ART scale up\nnone'), TEAM='dummy', DATA_T2='using\ntrue tree', OBJ_L2=c('%Acute\n(endpoint)')))
	set(tmp2, NULL, 'INT_L', tmp2[, factor(INT_L, levels=c('ART scale up\nfast','ART scale up\nslow','ART scale up\nnone'), labels=c('ART scale up\nfast','ART scale up\nslow','ART scale up\nnone'))])
	setnames(tmp2, 'TEAM','team')
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	ggplot(subset(tmp, TEAM!='True'), aes(y=INT_L, x=central, colour=TEAM)) +
			geom_point(data=tmp2, size=1, colour='transparent') +			
			geom_errorbarh(aes(xmin=lower95, xmax=upper95), height=0.3) + geom_point(size=3, pch=13) +
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(AC_L+DATA_T2~TEAM, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=2)) +
			labs(x= '%', y='', title='Primary objectives: %Acute at endpoint\n(Village)\n')	
	ggsave(file=paste(outdir,'/res_acrossTEAM_primaryvillage_pcacuteend','.pdf',sep=''), w=8, h=7)
	ggplot(subset(tmp, TEAM!='True'), aes(y=INT_L, x=central*100, colour=TEAM)) +
			geom_point(data=tmp2, size=1, colour='transparent') +			
			geom_errorbarh(aes(xmin=lower95*100, xmax=upper95*100), height=0.3) + geom_point(size=3, pch=13) +
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(AC_L+DATA_T2~TEAM, scales='free_y', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=2)) +
			labs(x= '%', y='', title='Primary objectives: %Acute at endpoint\n(Village)\n')	
	ggsave(file=paste(outdir,'/res_acrossTEAM_primaryvillage_pcacuteend_ssc','.pdf',sep=''), w=10, h=7)
	#	OVERALL	
	tmp	<- subset(df, OBJ%in%c('OBJ_ii','OBJ_iii','OBJ_v','OBJ_vi') & !grepl('(',TEAM,fixed=1) & DATAT_L=='Village' & USED_GENES=='all')
	tmp2<- as.data.table(expand.grid(central=c(0.01,0.55), AC_L=c('%Acute\nlow','%Acute\nhigh'), INT_L=c('ART scale up\nfast','ART scale up\nslow','ART scale up\nnone'), TEAM='dummy', DATA_T2='using\ntrue tree', OBJ_L2=c('%Acute\n(baseline)','%Acute\n(endpoint)')))
	set(tmp2, NULL, 'INT_L', tmp2[, factor(INT_L, levels=c('ART scale up\nfast','ART scale up\nslow','ART scale up\nnone'), labels=c('ART scale up\nfast','ART scale up\nslow','ART scale up\nnone'))])
	setnames(tmp2, 'TEAM','team')
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	ggplot(subset(tmp, TEAM!='True'), aes(y=INT_L, x=central, colour=TEAM)) +
			geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95), height=0.3) + geom_point(size=3, pch=13) +
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM+DATA_T2~OBJ_L2+AC_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.5, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=2)) +
			labs(x= 'Estimates', y='', title='Primary objectives: quantitative\n(Village)\n')	
	ggsave(file=paste(outdir,'/res_acrossTEAM_primaryvillage_nmbrs','.pdf',sep=''), w=13, h=13)
	#	qualitative
	tmp		<- subset(df, OBJ%in%c('OBJ_i','OBJ_iv') & !grepl('(',TEAM,fixed=1) & DATAT_L=='Village' & USED_GENES=='all')
	set(tmp, NULL, 'central', tmp[, as.character(central)])
	set(tmp, tmp[, which(OBJ=='OBJ_i' & central=='-1')], 'central', 'declining')
	set(tmp, tmp[, which(OBJ=='OBJ_i' & central=='0')], 'central', 'stable')
	set(tmp, tmp[, which(OBJ=='OBJ_i' & central=='1')], 'central', 'increasing')
	set(tmp, tmp[, which(OBJ=='OBJ_iv' & central=='-1')], 'central', '<15%')
	set(tmp, tmp[, which(OBJ=='OBJ_iv' & central=='0')], 'central', '15%-30%')
	set(tmp, tmp[, which(OBJ=='OBJ_iv' & central=='1')], 'central', '>30%')
	tmp2	<- c('declining','stable','increasing','<15%','15%-30%','>30%')
	set(tmp, NULL, 'central', tmp[, factor(central, levels=tmp2, labels=tmp2)])	
	tmp2	<- as.data.table(expand.grid(central=tmp2, AC_L=c('%Acute\nlow','%Acute\nhigh'), INT_L=c('ART scale up\nfast','ART scale up\nslow','ART scale up\nnone'), TEAM='dummy', DATA_T2='using\ntrue tree', OBJ_L2='Incidence\nTrend'))
	set(tmp2, tmp2[, which(grepl('%', central))], 'OBJ_L2', '%Acute Ctgr\n(baseline)')
	set(tmp2, NULL, 'INT_L', tmp2[, factor(INT_L, levels=c('ART scale up\nfast','ART scale up\nslow','ART scale up\nnone'), labels=c('ART scale up\nfast','ART scale up\nslow','ART scale up\nnone'))])	
	setkey(tmp2, central)
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp2, 'TEAM','team')
	setnames(tmp3, 'TEAM','team')
	ggplot(tmp2, aes(y=INT_L, x=central)) +
			#geom_errorbarh(height=0.3) + 
			geom_point(data=tmp2, size=1, colour='transparent') +
			geom_jitter(data=subset(tmp, TEAM!='True'), aes(colour=TEAM), size=3, pch=13, position = position_jitter(width=0, height=.15)) +
			geom_point(data=tmp3, size=3, colour='black', pch=18) +			
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM+DATA_T2~OBJ_L2+AC_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.5, "lines"), legend.position='bottom', panel.grid.minor= element_blank(), panel.grid.major= element_blank()) +
			guides(colour=guide_legend(ncol=2)) +
			labs(x= '\nEstimates', y='', title='Primary objectives: qualitative\n(Village)\n')	
	ggsave(file=paste(outdir,'/res_acrossTEAM_primaryvillage_qultv','.pdf',sep=''), w=13, h=13)	
	#invisible(lapply(c('OBJ_ii','OBJ_iii','OBJ_v','OBJ_vi'), function(x)
	#				{				
	#					tmp	<- subset(df, OBJ==x & !grepl('(',TEAM,fixed=1) & DATAT_L=='Village' & USED_GENES=='all')
	#					if(x=='OBJ_i')
	#						set(tmp, NULL, 'central', tmp[, factor(central, levels=c(1,0,-1), labels=c('increasing','stable','decreasing'))])
	#					if(x=='OBJ_iv')
	#						set(tmp, NULL, 'central', tmp[, factor(central, levels=c(1,0,-1), labels=c('>30%','15%-30%','<15%'))])	
	#					ggplot(subset(tmp, TEAM!='True'), aes(y=DATA_T, x=central, xmin=lower95, xmax=upper95, colour=TEAM)) +
	#							geom_errorbarh(height=0.3) + geom_point(size=3, pch=13) +
	#							geom_point(data=subset(tmp, TEAM=='True'), size=3, colour='black', pch=18) +
	#							scale_colour_manual(values=TEAM_CL) +
	#							facet_grid(AC_L~INT_L, scales='free', space='free_y') +
	#							theme_bw() + theme(panel.margin.x= unit(2, "lines")) + 
	#							labs(x= paste('\n',gsub('\n',' ',tmp[1, OBJ_L]),sep=''), y='')
	#					ggsave(file=paste(outdir,'/res_acrossTEAM_primaryvillage_',x,'.pdf',sep=''), w=10, h=3)
	#				}))
	

	#	SECONDARY: compare oversampling during intervention on regional
	df	<- subset(dfa, DATA_T=='phy' & SMPL_M=='extreme' & DATAT_L=='Regional')	
	tmp	<- subset(dfa, DATA_T=='phy' & SMPL_M=='much' & SMPL_C=='8%' & SMPL_N==1600 & IMPRT=='5%' & AC_T=='low' & INT_T!='none')
	df	<- rbind(tmp, df)
	tmp	<- subset(df, OBJ%in%c('OBJ_ii','OBJ_iii','OBJ_v','OBJ_vi') & !grepl('(',TEAM,fixed=1) & DATAT_L=='Regional' & USED_GENES=='all')
	tmp[, SMPL_L:= NA_character_]
	set(tmp, tmp[, which(SMPL_M=='much')], 'SMPL_L', 'sampling\n50% during intervention')
	set(tmp, tmp[, which(SMPL_M=='extreme')], 'SMPL_L', 'sampling\n85% during intervention')
	set(tmp, NULL, 'SMPL_L', tmp[, factor(SMPL_L, levels=c('sampling\n85% during intervention','sampling\n50% during intervention'), labels=c('sampling\n85% during intervention','sampling\n50% during intervention'))])
	ggplot(subset(tmp, TEAM!='True'), aes(y=SMPL_L, x=central, xmin=lower95, xmax=upper95, colour=TEAM)) +
			geom_errorbarh(height=0.3) + geom_point(size=3, pch=13) +
			geom_point(data=subset(tmp, TEAM=='True'), size=3, colour='black', pch=18) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(INT_L~OBJ_L2+AC_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.5, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=2)) +
			labs(x= 'Estimates', y='', title='Secondary objective: oversampling during intervention\n(Regional, using true tree)\n')	
	ggsave(file=paste(outdir,'/res_acrossTEAM_secondary_oversmplintrvntnregional','.pdf',sep=''), w=13, h=5.5)
	
		
	#	SECONDARY: compare imports high / low 
	df	<- subset(dfa, SC_RND%in%c('P','E','L','H','12','10'))
	tmp	<- subset(df, OBJ%in%c('OBJ_ii','OBJ_iii','OBJ_v','OBJ_vi') & !grepl('(',TEAM,fixed=1) & USED_GENES=='all')
	tmp[, IMPRT_L:= NA_character_]
	set(tmp, tmp[, which(IMPRT=='0%')], 'IMPRT_L', '0% trns/year from outside\n(Village)')
	set(tmp, tmp[, which(IMPRT=='5%')], 'IMPRT_L', '5% trns/year from outside\n(Regional)')
	set(tmp, tmp[, which(IMPRT=='20%')], 'IMPRT_L', '20% trns/year from outside\n(Regional)')		
	set(tmp, NULL, 'IMPRT_L', tmp[, factor(IMPRT_L, levels=rev(c('0% trns/year from outside\n(Village)','5% trns/year from outside\n(Regional)','20% trns/year from outside\n(Regional)')), labels=rev(c('0% trns/year from outside\n(Village)','5% trns/year from outside\n(Regional)','20% trns/year from outside\n(Regional)')))])
	set(tmp, NULL, 'INT_L', tmp[, factor(as.character(INT_L), levels=c('ART scale up\nfast','ART scale up\nslow'), labels=c('ART scale up\nfast','ART scale up\nslow'))])
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')
	tmp2	<- as.data.table(expand.grid(	central=0.1, IMPRT_L=rev(c('0% trns/year from outside\n(Village)','5% trns/year from outside\n(Regional)','20% trns/year from outside\n(Regional)')), 
											AC_L='%Acute\nhigh', INT_L=c('ART scale up\nfast','ART scale up\nslow'), TEAM=c('Imperial', 'Vancouver', 'Cambridge/Imperial', 'ETH Zurich'), OBJ_L2='%Incidence'))	
	
	ggplot(tmp2, aes(y=IMPRT_L, x=central, colour=TEAM)) +
			geom_point(size=1, colour='transparent') +
			geom_errorbarh(data=subset(tmp, TEAM!='True'), aes(xmin=lower95, xmax=upper95), height=0.3) + 			
			geom_point(data=subset(tmp, TEAM!='True'), size=3, pch=13) +
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM+INT_L~OBJ_L2+AC_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.5, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=2)) +
			labs(x= '\nEstimates', y='', title='Secondary objective: transmissions from outside\n(using true tree)\n')	
	ggsave(file=paste(outdir,'/res_acrossTEAM_secondary_cntm','.pdf',sep=''), w=13, h=12)
	
	
	# 	SECONDARY: compare duration sampling
	df	<- subset(dfa, SC_RND%in%c('O','F','T','L'))
	tmp	<- subset(df, OBJ%in%c('OBJ_ii','OBJ_iii','OBJ_v','OBJ_vi') & !grepl('(',TEAM,fixed=1) & USED_GENES=='all')
	tmp[, SMPL_L:= NA_character_]
	set(tmp, tmp[, which(SMPL_D=='3')], 'SMPL_L', '3 yr sampling duration\nafter intervention start')
	set(tmp, tmp[, which(SMPL_D=='5')], 'SMPL_L', '5 yr sampling duration\nafter intervention start')		
	set(tmp, NULL, 'SMPL_L', tmp[, factor(SMPL_L, levels=rev(c('3 yr sampling duration\nafter intervention start','5 yr sampling duration\nafter intervention start')), labels=rev(c('3 yr sampling duration\nafter intervention start','5 yr sampling duration\nafter intervention start')))])
	ggplot(subset(tmp, TEAM!='True'), aes(y=SMPL_L, x=central, xmin=lower95, xmax=upper95, colour=TEAM)) +
			geom_errorbarh(height=0.3) + geom_point(size=3, pch=13) +
			geom_point(data=subset(tmp, TEAM=='True'), size=3, colour='black', pch=18) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(INT_L~OBJ_L2+AC_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.5, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=2)) +
			labs(x= 'Estimates', y='', title='Secondary objective: sampling duration after intervention start\n(Regional, using true tree)\n')	
	ggsave(file=paste(outdir,'/res_acrossTEAM_secondary_sdurregional','.pdf',sep=''), w=13, h=4)
	
	
	# 	SECONDARY: compare seq coverage 
	df	<- subset(dfa, SC_RND%in%c('I','J','G','K','T','R','L','H','05','08','06','07','11','09','12','10'))
	tmp	<- subset(df, OBJ%in%c('OBJ_ii','OBJ_iii','OBJ_v','OBJ_vi') & !grepl('(',TEAM,fixed=1) & USED_GENES=='all')
	tmp[, SMPL_L:= NA_character_]
	set(tmp, tmp[, which(SMPL_C=='8%')], 'SMPL_L', '8% coverage (Regional)')
	set(tmp, tmp[, which(SMPL_C=='16%')], 'SMPL_L', '16% coverage (Regional)')
	set(tmp, tmp[, which(SMPL_C=='30%')], 'SMPL_L', '30% coverage (Village)')
	set(tmp, tmp[, which(SMPL_C=='60%')], 'SMPL_L', '60% coverage (Village)')		
	tmp2	<- c('8% coverage (Regional)','16% coverage (Regional)', '30% coverage (Village)', '60% coverage (Village)')	
	set(tmp, NULL, 'SMPL_L', tmp[, factor(SMPL_L, levels=rev(tmp2), labels=rev(tmp2))])
	set(tmp, NULL, 'INT_L', tmp[, factor(as.character(INT_L), levels=c('ART scale up\nfast','ART scale up\nslow'), labels=c('ART scale up\nfast','ART scale up\nslow'))])
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')			
	ggplot(subset(tmp, TEAM!='True'), aes(y=SMPL_L, x=central, xmin=lower95, xmax=upper95, colour=TEAM)) +
			geom_errorbarh(height=0.3) + geom_point(size=3, pch=13) +
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM+INT_L~OBJ_L2+AC_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.5, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=2)) +
			labs(x= 'Estimates', y='', title='Secondary objective: sampling coverage\n(using true tree)\n')	
	ggsave(file=paste(outdir,'/res_acrossTEAM_secondary_scvrg','.pdf',sep=''), w=13, h=10)
	tmp		<- subset(tmp, !(OBJ=='OBJ_ii' & central>0.1) & !(OBJ=='OBJ_v' & central>0.3 & AC_T=='low') & !(OBJ=='OBJ_vi' & central>0.3 & AC_T=='low'))
	tmp2	<- as.data.table(expand.grid(	central=0.1, SMPL_L=rev(c('8% coverage (Regional)','16% coverage (Regional)', '30% coverage (Village)', '60% coverage (Village)')), 
											AC_L='%Acute\nlow', INT_L=c('ART scale up\nfast','ART scale up\nslow'), TEAM=c('Imperial', 'Vancouver', 'Cambridge/Imperial', 'ETH Zurich'), OBJ_L2='%Incidence'))	
	ggplot(tmp2, aes(y=SMPL_L, x=central, colour=TEAM)) +
			geom_point(size=1, colour='transparent') +
			geom_point(data=subset(tmp, TEAM!='True'), size=3, pch=13) + geom_errorbarh(data=subset(tmp, TEAM!='True'), aes(xmin=lower95, xmax=upper95), height=0.3) + 
			geom_point(data=tmp3, size=3, colour='black', pch=18) +			
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM+INT_L~OBJ_L2+AC_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.5, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=2)) +
			labs(x= '\nEstimates', y='', title='Secondary objective: sampling coverage\n(using true tree, range cut)\n')	
	ggsave(file=paste(outdir,'/res_acrossTEAM_secondary_scvrgcut','.pdf',sep=''), w=13, h=11)
	
	
	
	dfi<- subset( dfa, select=c(SC_RND, DATAT_L, DATA_T, AC_T, INT_T, IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D) )
	setkey(dfi, DATAT_L, DATA_T, AC_T, INT_T, IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D)
	dfi	<- unique(dfi)
	file<- paste(outdir,'/../docs/SC_RND_info','.csv',sep='')
	write.csv(dfi, file=file, row.names=FALSE)
	
	tmp	<- subset(dfa, select=c(TEAM, DATA_T, DATAT_L, SIM_SCENARIO))
	setkey(tmp, TEAM, DATA_T, DATAT_L, SIM_SCENARIO)
	tmp	<- unique(tmp)
	tmp[, table(TEAM, DATA_T, DATAT_L )]
}
##--------------------------------------------------------------------------------------------------------
##	evaluate results
##	olli 08.05.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Feb2015.evaluate.read<- function()
{
	#	read truth for regional simus
	indir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/PANGEAHIVsim_internal/documents/external/2015_05_results'	
	file	<- paste(indir, '/answers_Regional_Feb2015_rFormat.csv', sep='')
	df		<- read.submission(file, verbose=0, reset.OBJiv.conservative=1)
	#	read truth for village simus
	file	<- paste(indir, '/answers_Village_Feb2015-yr43_rFormat.csv', sep='')
	tmp		<- read.submission(file, verbose=0, reset.OBJiv.conservative=1)
	set(tmp, NULL, 'TEAM', 'True')
	df		<- rbind(df, tmp)
	#	read submissions
	tmp		<- list.files(indir, pattern='csv$')
	tmp		<- tmp[!grepl('answers',tmp)]
	#	read Eriks multiple submissions
	tmp2	<- data.table(FILE=tmp[grepl('cambImp',tmp)])
	tmp2[, RUN:= tmp2[,  sapply( strsplit(FILE,'_'), function(x) rev(x)[1] )]]
	set(tmp2, NULL, 'RUN', tmp2[, substr(RUN, 1, nchar(RUN)-4)])
	set(tmp2, NULL, 'RUN', tmp2[, gsub('results0','',RUN)])
	dfs		<- do.call('rbind',lapply(seq_len(nrow(tmp2)), function(i)
				{
					z	<- read.submission( paste(indir, '/', tmp2[i, FILE], sep=''), verbose=0, reset.OBJiv.conservative=1 )
					set(z, NULL, 'TEAM', z[, paste(TEAM, ' (', tmp2[i, RUN], ')', sep='')])
					z
				}))
	tmp		<- tmp[!grepl('cambImp',tmp)]
	tmp		<- do.call('rbind',lapply(tmp, function(x) read.submission(paste(indir,'/',x,sep=''), verbose=0, reset.OBJiv.conservative=1)))
	dfs		<- rbind(dfs, tmp)
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
	dfi			<- data.table(FILE=list.files('/Users/Oliver/duke/2014_Gates/methods_comparison_pipeline/FINAL', '.*zip$', full.names=FALSE))	
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
	
	dfa
}
##--------------------------------------------------------------------------------------------------------
##	check simulated sequences: create ExaML tree and estimate R2
##	olli 09.02.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Feb2015.randomize<- function()
{	
	set.seed(42)
	indir		<- '/Users/Oliver/duke/2014_Gates/methods_comparison_pipeline/FINAL'
	outdir		<- '/Users/Oliver/duke/2014_Gates/methods_comparison_pipeline/FINAL_RND'
	dfi			<- data.table(FILE=list.files(indir, '.*zip$', full.names=FALSE))
	
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
	
	dfi[, GSUB_FROM:= sapply(strsplit(FILE, '_'),'[[',3)]
	dfi[, GSUB_TO:= paste(OBJECTIVE,'_sc',SC_RND,sep='')]
	
	tmp		<- dfi[, {
				cat(paste('process file',FILE))
				tmpdir			<- paste(indir,'/tmp',SC_RND,sep='')
				unzip( paste(indir,'/',FILE,sep=''), exdir=tmpdir)
				files.name.o	<- list.files( paste(indir,'/tmp',SC_RND,sep='') )
				files.name.n	<- gsub(GSUB_FROM,GSUB_TO,files.name.o)
				files.name.n	<- gsub('HPTN071','PANGEAsim_Regional',files.name.n)
				file.rename( paste(tmpdir, files.name.o, sep='/'), paste(tmpdir, files.name.n, sep='/'))
				newzip			<- gsub(GSUB_FROM,GSUB_TO,FILE)
				newzip			<- gsub('HPTN071','PANGEAsim_Regional',newzip)
				zip(  paste(outdir,'/',newzip,sep=''), paste(tmpdir, files.name.n, sep='/'), flags = "-FSr9XTj")
				file.remove( paste(tmpdir, files.name.n, sep='/'), showWarnings=FALSE)
				file.remove( tmpdir, showWarnings=FALSE )				
			}, by='DUMMY']
}
##--------------------------------------------------------------------------------------------------------
##	check simulated sequences: create ExaML tree and estimate R2
##	olli 23.10.14
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.November2014<- function()
{		
	#	with CD4
	if(0)
	{
		indir			<- system.file(package="rPANGEAHIVsim", "misc")
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
														s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0, 
														epi.model='HPTN071', epi.dt=1/48, epi.import=0.05,
														v.N0tau=1, v.r=2.851904, v.T50=-2,
														wher.mu=log(0.00447743)-0.3^2/2, wher.sigma=0.3, bwerm.mu=log(0.002239075)-0.13^2/2, bwerm.sigma=0.13, er.gamma=NA,
														dbg.GTRparam=0, dbg.rER=0, index.starttime.mode='normal', startseq.mode='sample', seqtime.mode=NA)						
		# standard run, fixed GTR param + fixed relative rates for each transmission chain
		# no WH Er to equal BH ER, and both being constant
		pipeline.vary	<- data.table(	seqtime.mode=c('Gamma9','Gamma3','Unif12'),				
										label=c('-n0005g9','-n0005g3','-n0005u'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='seqtime.mode' ), 'v', as.character(seqtime.mode))
					print(pipeline.args)
					#	scenario A		
					infile.ind		<- '211014_RUN123_SCENARIO_0'
					infile.trm		<- '211014_RUN123_SCENARIO_0'
					tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp141024'
					tmpdir			<- paste(tmpdir,label,sep='')
					dir.create(tmpdir, showWarnings=FALSE)																		
					set(pipeline.args, which( pipeline.args$stat=='s.PREV.max' ), 'v', '0.11')
					file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
					file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
					file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
					system(file)					
				}, by='label']
	}
	#	with CD4, lower WH rate
	if(1)
	{
		indir			<- system.file(package="rPANGEAHIVsim", "misc")
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
				s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0, 
				epi.model='HPTN071', epi.dt=1/48, epi.import=0.20,
				v.N0tau=1, v.r=2.851904, v.T50=-2,
				wher.mu=log(0.003358613)-0.3^2/2, wher.sigma=0.3, bwerm.mu=log(0.002239075)-0.13^2/2, bwerm.sigma=0.13, er.gamma=NA,
				dbg.GTRparam=0, dbg.rER=0, index.starttime.mode='fix', startseq.mode='sample', seqtime.mode=NA)						
		# standard run, fixed GTR param + fixed relative rates for each transmission chain
		# no WH Er to equal BH ER, and both being constant
		pipeline.vary	<- data.table(	seqtime.mode=c('Unif12'),				
										label=c('-n0005ul20'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='seqtime.mode' ), 'v', as.character(seqtime.mode))
					print(pipeline.args)
					if(1)
					{
						#	scenario A					
						infile.ind		<- '271014_HPTN071_scA_rep1'
						infile.trm		<- '271014_HPTN071_scA_rep1'
						tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp141025a'
						tmpdir			<- paste(tmpdir,label,sep='')
						dir.create(tmpdir, showWarnings=FALSE)																		
						set(pipeline.args, which( pipeline.args$stat=='s.PREV.max' ), 'v', '0.11')
						file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
						file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
						file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
						system(file)													
					}
				}, by='label']
	}
	#	with CD4, lower WH rate
	if(1)
	{
		indir			<- system.file(package="rPANGEAHIVsim", "misc")
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
														s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0, 
														epi.model='HPTN071', epi.dt=1/48, epi.import=0.05,
														v.N0tau=1, v.r=2.851904, v.T50=-2,
														wher.mu=log(0.003358613)-0.3^2/2, wher.sigma=0.3, bwerm.mu=log(0.002239075)-0.13^2/2, bwerm.sigma=0.13, er.gamma=NA,
														dbg.GTRparam=0, dbg.rER=0, index.starttime.mode='normal', startseq.mode='sample', seqtime.mode=NA)						
		# standard run, fixed GTR param + fixed relative rates for each transmission chain
		# no WH Er to equal BH ER, and both being constant
		pipeline.vary	<- data.table(	seqtime.mode=c('Gamma9','Gamma3','Unif12'),				
										label=c('-n0005g9l','-n0005g3l','-n0005ul'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='seqtime.mode' ), 'v', as.character(seqtime.mode))
					print(pipeline.args)
					if(1)
					{
						#	scenario A					
						infile.ind		<- '271014_HPTN071_scA_rep1'
						infile.trm		<- '271014_HPTN071_scA_rep1'
						tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp141025a'
						tmpdir			<- paste(tmpdir,label,sep='')
						dir.create(tmpdir, showWarnings=FALSE)																		
						set(pipeline.args, which( pipeline.args$stat=='s.PREV.max' ), 'v', '0.11')
						file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
						file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
						file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
						system(file)													
					}
					if(0)
					{
						#	scenario B						
						infile.ind		<- '271014_HPTN071_scB_rep1'
						infile.trm		<- '271014_HPTN071_scB_rep1'	
						tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp141025b'
						tmpdir			<- paste(tmpdir,label,sep='')
						dir.create(tmpdir, showWarnings=FALSE)																														
						set(pipeline.args, which( pipeline.args$stat=='s.PREV.max' ), 'v', '0.15')
						file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
						file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
						file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
						system(file)						
					}
					if(0)
					{
						#	scenario C						
						infile.ind		<- '271014_HPTN071_scC_rep1'
						infile.trm		<- '271014_HPTN071_scC_rep1'	
						tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp141025c'
						tmpdir			<- paste(tmpdir,label,sep='')
						dir.create(tmpdir, showWarnings=FALSE)																														
						set(pipeline.args, which( pipeline.args$stat=='s.PREV.max' ), 'v', '0.185')
						file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
						file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
						file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
						system(file)	
					}					
				}, by='label']
	}
	#	with CD4, no introductions after baseline for Leventhal et al
	if(1)
	{
		indir			<- system.file(package="rPANGEAHIVsim", "misc")
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
				s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0, 
				epi.model='HPTN071', epi.dt=1/48, epi.import=0,
				v.N0tau=1, v.r=2.851904, v.T50=-2,
				wher.mu=log(0.003358613)-0.3^2/2, wher.sigma=0.3, bwerm.mu=log(0.002239075)-0.13^2/2, bwerm.sigma=0.13, er.gamma=NA,
				dbg.GTRparam=0, dbg.rER=0, index.starttime.mode='normal', startseq.mode='sample', seqtime.mode=NA)						
		# standard run, fixed GTR param + fixed relative rates for each transmission chain
		# no WH Er to equal BH ER, and both being constant
		pipeline.vary	<- data.table(	seqtime.mode=c('Unif12'),				
										label=c('-Leventhal00ul'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='seqtime.mode' ), 'v', as.character(seqtime.mode))
					print(pipeline.args)
					if(1)
					{
						#	scenario A					
						infile.ind		<- '271014_HPTN071_scA_rep1'
						infile.trm		<- '271014_HPTN071_scA_rep1'
						tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp141114a'
						tmpdir			<- paste(tmpdir,label,sep='')
						dir.create(tmpdir, showWarnings=FALSE)																		
						set(pipeline.args, which( pipeline.args$stat=='s.PREV.max' ), 'v', '0.11')
						file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
						file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
						file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
						#system(file)													
					}
					if(1)
					{
						#	scenario B						
						infile.ind		<- '271014_HPTN071_scB_rep1'
						infile.trm		<- '271014_HPTN071_scB_rep1'	
						tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp141114b'
						tmpdir			<- paste(tmpdir,label,sep='')
						dir.create(tmpdir, showWarnings=FALSE)																														
						set(pipeline.args, which( pipeline.args$stat=='s.PREV.max' ), 'v', '0.15')
						file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
						file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
						file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
						#system(file)						
					}
					if(1)
					{
						#	scenario C						
						infile.ind		<- '271014_HPTN071_scC_rep1'
						infile.trm		<- '271014_HPTN071_scC_rep1'	
						tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp141114c'
						tmpdir			<- paste(tmpdir,label,sep='')
						dir.create(tmpdir, showWarnings=FALSE)																														
						set(pipeline.args, which( pipeline.args$stat=='s.PREV.max' ), 'v', '0.185')
						file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
						file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
						file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
						#system(file)	
					}					
				}, by='label']
	}
	#	with CD4, same WH rate for Suchard et al
	if(1)
	{
		indir			<- system.file(package="rPANGEAHIVsim", "misc")
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
				s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0, 
				epi.model='HPTN071', epi.dt=1/48, epi.import=0.05,
				v.N0tau=1, v.r=2.851904, v.T50=-2,
				wher.mu=log(0.002239075)-0.13^2/2, wher.sigma=0.13, bwerm.mu=log(0.002239075)-0.13^2/2, bwerm.sigma=0.13, er.gamma=NA,
				dbg.GTRparam=1, dbg.rER=1, index.starttime.mode='normal', startseq.mode='sample', seqtime.mode=NA)						
		# standard run, fixed GTR param + fixed relative rates for each transmission chain
		# no WH Er to equal BH ER, and both being constant
		pipeline.vary	<- data.table(	seqtime.mode=c('Unif12'),				
										label=c('-Suchard05ul'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='seqtime.mode' ), 'v', as.character(seqtime.mode))
					print(pipeline.args)
					if(1)
					{
						#	scenario A					
						infile.ind		<- '271014_HPTN071_scA_rep1'
						infile.trm		<- '271014_HPTN071_scA_rep1'
						tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp141114a'
						tmpdir			<- paste(tmpdir,label,sep='')
						dir.create(tmpdir, showWarnings=FALSE)																		
						set(pipeline.args, which( pipeline.args$stat=='s.PREV.max' ), 'v', '0.11')
						file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
						file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
						file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
						#system(file)													
					}
					if(1)
					{
						#	scenario B						
						infile.ind		<- '271014_HPTN071_scB_rep1'
						infile.trm		<- '271014_HPTN071_scB_rep1'	
						tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp141114b'
						tmpdir			<- paste(tmpdir,label,sep='')
						dir.create(tmpdir, showWarnings=FALSE)																														
						set(pipeline.args, which( pipeline.args$stat=='s.PREV.max' ), 'v', '0.15')
						file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
						file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
						file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
						#system(file)						
					}
					if(1)
					{
						#	scenario C						
						infile.ind		<- '271014_HPTN071_scC_rep1'
						infile.trm		<- '271014_HPTN071_scC_rep1'	
						tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp141114c'
						tmpdir			<- paste(tmpdir,label,sep='')
						dir.create(tmpdir, showWarnings=FALSE)																														
						set(pipeline.args, which( pipeline.args$stat=='s.PREV.max' ), 'v', '0.185')
						file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
						file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
						file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
						#system(file)	
					}					
				}, by='label']
	}
	#	all replicates
	if(1)
	{
		indir			<- '/Users/Oliver/duke/2014_Gates/methods_comparison_pipeline/Pangea_Oct28_2014_sim'
		#indir			<- '/Users/Oliver/git/Pangea__Oct27_2014_code'
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
				s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0, 
				epi.model='HPTN071', epi.dt=1/48, epi.import=0.05,
				v.N0tau=1, v.r=2.851904, v.T50=-2,
				wher.mu=log(0.003358613)-0.3^2/2, wher.sigma=0.3, bwerm.mu=log(0.002239075)-0.13^2/2, bwerm.sigma=0.13, er.gamma=0,
				dbg.GTRparam=0, dbg.rER=0, index.starttime.mode='normal', startseq.mode='sample', seqtime.mode='Unif12')						
		
		pipeline.vary	<- data.table(	REP=1:100, label=paste('-r',1:100,sep='')  )		
		dummy			<- pipeline.vary[, {									
					if(1)
					{
						#	scenario A					
						infile.ind		<- paste('281014_HPTN071_scA_rep',REP,sep='')
						infile.trm		<- paste('281014_HPTN071_scA_rep',REP,sep='')
						tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp141026a'
						tmpdir			<- paste(tmpdir,label,sep='')
						dir.create(tmpdir, showWarnings=FALSE)																		
						set(pipeline.args, which( pipeline.args$stat=='s.PREV.max' ), 'v', '0.11')
						file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,'_IND.csv',sep=''))
						file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,'_TRM.csv',sep=''))
						file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,'_IND.csv',sep=''), paste(infile.trm,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
						system(file)	
					}
					if(1)
					{
						#	scenario B						
						infile.ind		<- paste('281014_HPTN071_scB_rep',REP,sep='')
						infile.trm		<- paste('281014_HPTN071_scB_rep',REP,sep='')
						tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp141026b'
						tmpdir			<- paste(tmpdir,label,sep='')
						dir.create(tmpdir, showWarnings=FALSE)																														
						set(pipeline.args, which( pipeline.args$stat=='s.PREV.max' ), 'v', '0.15')
						file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,'_IND.csv',sep=''))
						file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,'_TRM.csv',sep=''))
						file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,'_IND.csv',sep=''), paste(infile.trm,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
						system(file)						
					}
					if(1)
					{
						#	scenario C						
						infile.ind		<- paste('281014_HPTN071_scC_rep',REP,sep='')
						infile.trm		<- paste('281014_HPTN071_scC_rep',REP,sep='')	
						tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp141026c'
						tmpdir			<- paste(tmpdir,label,sep='')
						dir.create(tmpdir, showWarnings=FALSE)																														
						set(pipeline.args, which( pipeline.args$stat=='s.PREV.max' ), 'v', '0.185')
						file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,'_IND.csv',sep=''))
						file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,'_TRM.csv',sep=''))
						file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,'_IND.csv',sep=''), paste(infile.trm,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
						system(file)	
					}					
				}, by='label']
	}
}
##--------------------------------------------------------------------------------------------------------
##	check simulated sequences: create ExaML tree and estimate R2
##	olli 14.09.14
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.October2014<- function()
{	
	if(0)
	{
		indir			<- system.file(package="rPANGEAHIVsim", "misc")
		infile.ind		<- '140716_RUN001'
		infile.trm		<- '140716_RUN001'	
		
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
				s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.2, 
				epi.model='HPTN071', epi.dt=1/48, epi.import=0.1,
				v.N0tau=3.58e4, v.r=2, v.T50=-1,
				wher.mu=NA, wher.sigma=NA,
				bwerm.mu=1.5, bwerm.sigma=0.12 )	
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
	if(0)
	{
		indir			<- system.file(package="rPANGEAHIVsim", "misc")
		infile.ind		<- '140716_RUN001'
		infile.trm		<- '140716_RUN001'	
		
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
														s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.2, 
														epi.model='HPTN071', epi.dt=1/48, epi.import=0.1,
														v.N0tau=3.58e4, v.r=2, v.T50=-1,
														wher.mu=0.004, wher.sigma=0.7,
														bwerm.mu=NA, bwerm.sigma=NA )	
		pipeline.vary	<- data.table(bwerm.mu=c(1.5, 1.75, 2.0), bwerm.sigma=c(0.12, 0.103, 0.09), label=c('-b15','-b17','-b20'))						
		dummy			<- pipeline.vary[, {				
												set(pipeline.args, which( pipeline.args$stat=='bwerm.mu' ), 'v', as.character(bwerm.mu))
												set(pipeline.args, which( pipeline.args$stat=='bwerm.sigma' ), 'v', as.character(bwerm.sigma))
												print(pipeline.args)
												#	re-name the following:
												tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140915'
												tmpdir			<- paste(tmpdir,label,sep='')
												dir.create(tmpdir, showWarnings=FALSE)						
												#						
												file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
												file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
												file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
												system(file)
											}, by='label']
	}
	if(0)
	{
		indir			<- system.file(package="rPANGEAHIVsim", "misc")
		infile.ind		<- '140716_RUN001'
		infile.trm		<- '140716_RUN001'	
		
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
														s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.2, 
														epi.model='HPTN071', epi.dt=1/48, epi.import=0.1,
														v.N0tau=1, v.r=NA, v.T50=-2,
														wher.mu=0.004, wher.sigma=0.7,
														bwerm.mu=1.75, bwerm.sigma=0.103 )	
		pipeline.vary	<- data.table(v.r=c(6.305779, 5.154461, 4.003191, 2.851904), label=c('-n5','-n4','-n3','-n2'))						
		dummy			<- pipeline.vary[, {				
												set(pipeline.args, which( pipeline.args$stat=='v.r' ), 'v', as.character(v.r))					
												print(pipeline.args)
												#	re-name the following:
												tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140915b'
												tmpdir			<- paste(tmpdir,label,sep='')
												dir.create(tmpdir, showWarnings=FALSE)						
												#						
												file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
												file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
												file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
												system(file)
											}, by='label']
	}							
	if(0)
	{
		indir			<- system.file(package="rPANGEAHIVsim", "misc")
		infile.ind		<- '140716_RUN001'
		infile.trm		<- '140716_RUN001'	
		
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
														s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.2, 
														epi.model='HPTN071', epi.dt=1/48, epi.import=0.1,
														v.N0tau=1, v.r=2.851904, v.T50=-2,
														wher.mu=NA, wher.sigma=NA,
														bwerm.mu=4, bwerm.sigma=0.06 )	
		pipeline.vary	<- data.table(wher.mu=exp(c(-5.477443+0.18, -5.071977+0.08, -4.784295+0.045)), wher.sigma=c(0.6, 0.4, 0.3), label=c('-2','-3','-4'))						
		dummy			<- pipeline.vary[, {				
												set(pipeline.args, which( pipeline.args$stat=='wher.mu' ), 'v', as.character(wher.mu))
												set(pipeline.args, which( pipeline.args$stat=='wher.sigma' ), 'v', as.character(wher.sigma))
												print(pipeline.args)
												#	re-name the following:
												tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140916'
												tmpdir			<- paste(tmpdir,label,sep='')
												dir.create(tmpdir, showWarnings=FALSE)						
												#						
												file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
												file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
												file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
												system(file)
											}, by='label']
	}
	#	changed: transmission edge model; fixed: transmission edges are within transmitter
	if(0)
	{
		indir			<- system.file(package="rPANGEAHIVsim", "misc")
		infile.ind		<- '140716_RUN001'
		infile.trm		<- '140716_RUN001'	
		
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
														s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.2, 
														epi.model='HPTN071', epi.dt=1/48, epi.import=0.1,
														v.N0tau=1, v.r=2.851904, v.T50=-2,
														wher.mu=exp(-5.071977+0.08), wher.sigma=0.4,
														bwerm.mu=NA, 
														bwerm.sigma=NA)	
		#	bwerm.mu is now shift										
		pipeline.vary	<- data.table(bwerm.mu=c(0, -2e-4, -4e-4, -6e-4, -8e-4, -1e-3), label=c('-s0','-s2','-s4','-s6','-s8','-s10'))						
		dummy			<- pipeline.vary[, {				
												set(pipeline.args, which( pipeline.args$stat=='bwerm.mu' ), 'v', as.character(bwerm.mu))												
												print(pipeline.args)
												#	re-name the following:
												tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140916b'
												tmpdir			<- paste(tmpdir,label,sep='')
												dir.create(tmpdir, showWarnings=FALSE)						
												#						
												file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
												file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
												file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
												system(file)
											}, by='label']
		#	rerun: bug?? NOPE!
		indir			<- system.file(package="rPANGEAHIVsim", "misc")
		indir			<- ifelse(indir=='','/Users/Oliver/git/HPTN071sim/raw_trchain',indir)
		infile.ind		<- '140716_RUN001'
		infile.trm		<- '140716_RUN001'	
		
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
				s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.2, 
				epi.model='HPTN071', epi.dt=1/48, epi.import=0.1,
				v.N0tau=1, v.r=2.851904, v.T50=-2,
				wher.mu=exp(-5.071977+0.08), wher.sigma=0.4,
				bwerm.mu=NA, 
				bwerm.sigma=NA)	
		#	bwerm.mu is now shift										
		pipeline.vary	<- data.table(bwerm.mu=c(0, -5e-4,  -1e-3), label=c('-s0','-s5','-s10'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='bwerm.mu' ), 'v', as.character(bwerm.mu))												
					print(pipeline.args)
					#	re-name the following:
					tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140917'
					tmpdir			<- paste(tmpdir,label,sep='')
					dir.create(tmpdir, showWarnings=FALSE)						
					#						
					file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
					file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
					file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
					system(file)
				}, by='label']
	}
	#	changed: only edges leading to tip are under within host evol (dead ends); added sdlog to transmission ER model; changed params to log scale
	if(0)
	{
		indir			<- system.file(package="rPANGEAHIVsim", "misc")
		infile.ind		<- '140716_RUN001'
		infile.trm		<- '140716_RUN001'	
		
		#	scenario 3* faster WH up to 0.02
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
														s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.2, 
														epi.model='HPTN071', epi.dt=1/48, epi.import=0.1,
														v.N0tau=1, v.r=2.851904, v.T50=-2,
														wher.mu=log(0.006716145)-0.37^2/2, wher.sigma=0.37,
														bwerm.mu=NA, 
														bwerm.sigma=NA)	
		#	bwerm.mu is now shift										
		pipeline.vary	<- data.table(bwerm.mu=c(log(0.002239075)-0.05^2/2, log(0.002239075-0.0005)-0.065^2/2,  log(0.002239075-0.001)-0.09^2/2), bwerm.sigma=c(0.05, 0.065, 0.09), label=c('-sh0','-sh5','-sh10'))						
		dummy			<- pipeline.vary[, {				
												set(pipeline.args, which( pipeline.args$stat=='bwerm.mu' ), 'v', as.character(bwerm.mu))
												set(pipeline.args, which( pipeline.args$stat=='bwerm.sigma' ), 'v', as.character(bwerm.sigma))
												print(pipeline.args)
												#	re-name the following:
												tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140918'
												tmpdir			<- paste(tmpdir,label,sep='')
												dir.create(tmpdir, showWarnings=FALSE)						
												#						
												file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
												file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
												file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
												system(file)
											}, by='label']
		#	scenario 2* faster WH limited to 0.01
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
				s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.2, 
				epi.model='HPTN071', epi.dt=1/48, epi.import=0.1,
				v.N0tau=1, v.r=2.851904, v.T50=-2,
				wher.mu=log(0.00447743)-0.3^2/2, wher.sigma=0.3,
				bwerm.mu=NA, bwerm.sigma=NA)	
		#	bwerm.mu is now shift										
		pipeline.vary	<- data.table(bwerm.mu=c(log(0.002239075)-0.05^2/2, log(0.002239075-0.0005)-0.065^2/2,  log(0.002239075-0.001)-0.09^2/2), bwerm.sigma=c(0.05, 0.065, 0.09), label=c('-sl0','-sl5','-sl10'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='bwerm.mu' ), 'v', as.character(bwerm.mu))
					set(pipeline.args, which( pipeline.args$stat=='bwerm.sigma' ), 'v', as.character(bwerm.sigma))
					print(pipeline.args)
					#	re-name the following:
					tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140918'
					tmpdir			<- paste(tmpdir,label,sep='')
					dir.create(tmpdir, showWarnings=FALSE)						
					#						
					file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
					file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
					file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
					system(file)
				}, by='label']
	}
	#	changed: time to seq limited to max 6 years instead of 35 ;-); sample fraction 10%
	if(0)
	{
		indir			<- system.file(package="rPANGEAHIVsim", "misc")
		infile.ind		<- '170914_HPTN071_scA_rep1'
		infile.trm		<- '170914_HPTN071_scA_rep1'	
		#	scenario 2* faster WH limited to 0.01 
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
														s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.1, 
														epi.model='HPTN071', epi.dt=1/48, epi.import=0.1,
														v.N0tau=1, v.r=2.851904, v.T50=-2,
														wher.mu=log(0.00447743)-0.3^2/2, wher.sigma=0.3,
														bwerm.mu=NA, bwerm.sigma=NA)	
		#	bwerm.mu is now shift										
		pipeline.vary	<- data.table(bwerm.mu=c(log(0.002239075)-0.05^2/2, log(0.002239075-0.0005)-0.065^2/2,  log(0.002239075-0.001)-0.09^2/2), bwerm.sigma=c(0.05, 0.065, 0.09), label=c('-sl0','-sl5','-sl10'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='bwerm.mu' ), 'v', as.character(bwerm.mu))
					set(pipeline.args, which( pipeline.args$stat=='bwerm.sigma' ), 'v', as.character(bwerm.sigma))
					print(pipeline.args)
					#	re-name the following:
					tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140918b'
					tmpdir			<- paste(tmpdir,label,sep='')
					dir.create(tmpdir, showWarnings=FALSE)						
					#						
					file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
					file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
					file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
					system(file)
				}, by='label']
	}
	if(0)
	{
		indir			<- system.file(package="rPANGEAHIVsim", "misc")
		infile.ind		<- '170914_HPTN071_scA_rep1'
		infile.trm		<- '170914_HPTN071_scA_rep1'	
		#	scenario 2* faster WH limited to 0.01 
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
														s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.1, 
														epi.model='HPTN071', epi.dt=1/48, epi.import=0.1,
														v.N0tau=1, v.r=2.851904, v.T50=-2,
														wher.mu=NA, wher.sigma=NA,
														bwerm.mu=log(0.002239075)-0.13^2/2, bwerm.sigma=0.13)	
		#	bwerm.mu is now shift										
		pipeline.vary	<- data.table(wher.mu=c(log(0.00447743)-0.3^2/2, log(0.00447743)-0.4^2/2,  log(0.00447743)-0.5^2/2), wher.sigma=c(0.3, 0.4, 0.5), label=c('-w3','-w4','-w5'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='wher.mu' ), 'v', as.character(wher.mu))
					set(pipeline.args, which( pipeline.args$stat=='wher.sigma' ), 'v', as.character(wher.sigma))
					print(pipeline.args)
					#	re-name the following:
					tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140918c'
					tmpdir			<- paste(tmpdir,label,sep='')
					dir.create(tmpdir, showWarnings=FALSE)						
					#						
					file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
					file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
					file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
					system(file)
				}, by='label']
	}
	#	sense check on CoV
	if(0)
	{
		indir			<- system.file(package="rPANGEAHIVsim", "misc")
		infile.ind		<- '180914_HPTN071_scA_rep1'
		infile.trm		<- '180914_HPTN071_scA_rep1'			 
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
														s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.1, 
														epi.model='HPTN071', epi.dt=1/48, epi.import=0.1,
														v.N0tau=1, v.r=2.851904, v.T50=-2,
														wher.mu=NA, wher.sigma=NA, bwerm.mu=NA, bwerm.sigma=NA)
		# no WH Er to equal BH ER, and both being constant										
		pipeline.vary	<- data.table(wher.mu=c(log(0.002239075)-0.13^2/2, log(0.002239075)-0.01^2/2), wher.sigma=c(0.13, 0.01), label=c('-c13','-c0'))						
		dummy			<- pipeline.vary[, {				
												set(pipeline.args, which( pipeline.args$stat=='wher.mu' ), 'v', as.character(wher.mu))
												set(pipeline.args, which( pipeline.args$stat=='wher.sigma' ), 'v', as.character(wher.sigma))
												set(pipeline.args, which( pipeline.args$stat=='bwerm.mu' ), 'v', as.character(wher.mu))
												set(pipeline.args, which( pipeline.args$stat=='bwerm.sigma' ), 'v', as.character(wher.sigma))
												
												print(pipeline.args)
												#	re-name the following:
												tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140919'
												tmpdir			<- paste(tmpdir,label,sep='')
												dir.create(tmpdir, showWarnings=FALSE)						
												#						
												file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
												file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
												file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
												system(file)
											}, by='label']
	}
	#	changed: ignore variation in meanRate by gene, added debug options
	if(0)
	{
		indir			<- system.file(package="rPANGEAHIVsim", "misc")
		infile.ind		<- '180914_HPTN071_scA_rep1'
		infile.trm		<- '180914_HPTN071_scA_rep1'			 
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
														s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.1, 
														epi.model='HPTN071', epi.dt=1/48, epi.import=0.1,
														v.N0tau=1, v.r=2.851904, v.T50=-2,
														wher.mu=NA, wher.sigma=NA, bwerm.mu=NA, bwerm.sigma=NA,
														dbg.GTRparam=1, dbg.rER=1)						
		# standard run, fixed GTR param + fixed relative rates for each transmission chain
		# no WH Er to equal BH ER, and both being constant
		pipeline.vary	<- data.table(	wher.mu=c(log(0.00447743)-0.3^2/2, log(0.002239075)-0.13^2/2, log(0.002239075)-0.01^2/2), wher.sigma=c(0.3, 0.13, 0.01), 
										bwerm.mu=c(log(0.002239075)-0.13^2/2, log(0.002239075)-0.13^2/2, log(0.002239075)-0.01^2/2), bwerm.sigma=c(0.13, 0.13, 0.01), label=c('-d1st','-d1s','-d1f'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='wher.mu' ), 'v', as.character(wher.mu))
					set(pipeline.args, which( pipeline.args$stat=='wher.sigma' ), 'v', as.character(wher.sigma))
					set(pipeline.args, which( pipeline.args$stat=='bwerm.mu' ), 'v', as.character(bwerm.mu))
					set(pipeline.args, which( pipeline.args$stat=='bwerm.sigma' ), 'v', as.character(bwerm.sigma))
					
					print(pipeline.args)
					#	re-name the following:
					tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140920a'
					tmpdir			<- paste(tmpdir,label,sep='')
					dir.create(tmpdir, showWarnings=FALSE)						
					#						
					file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
					file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
					file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
					system(file)
				}, by='label']
		#	debug only GTRparam
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
														s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.1, 
														epi.model='HPTN071', epi.dt=1/48, epi.import=0.1,
														v.N0tau=1, v.r=2.851904, v.T50=-2,
														wher.mu=NA, wher.sigma=NA, bwerm.mu=NA, bwerm.sigma=NA,
														dbg.GTRparam=1, dbg.rER=0)						
		# standard run, fixed GTR param + fixed relative rates for each transmission chain
		# no WH Er to equal BH ER, and both being constant
		pipeline.vary	<- data.table(	wher.mu=c(log(0.00447743)-0.3^2/2, log(0.002239075)-0.13^2/2, log(0.002239075)-0.01^2/2), wher.sigma=c(0.3, 0.13, 0.01), 
										bwerm.mu=c(log(0.002239075)-0.13^2/2, log(0.002239075)-0.13^2/2, log(0.002239075)-0.01^2/2), bwerm.sigma=c(0.13, 0.13, 0.01), label=c('-d10st','-d10s','-d10f'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='wher.mu' ), 'v', as.character(wher.mu))
					set(pipeline.args, which( pipeline.args$stat=='wher.sigma' ), 'v', as.character(wher.sigma))
					set(pipeline.args, which( pipeline.args$stat=='bwerm.mu' ), 'v', as.character(bwerm.mu))
					set(pipeline.args, which( pipeline.args$stat=='bwerm.sigma' ), 'v', as.character(bwerm.sigma))
					
					print(pipeline.args)
					#	re-name the following:
					tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140920a'
					tmpdir			<- paste(tmpdir,label,sep='')
					dir.create(tmpdir, showWarnings=FALSE)						
					#						
					file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
					file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
					file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
					system(file)
				}, by='label']
		#	no debug
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
														s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.1, 
														epi.model='HPTN071', epi.dt=1/48, epi.import=0.1,
														v.N0tau=1, v.r=2.851904, v.T50=-2,
														wher.mu=NA, wher.sigma=NA, bwerm.mu=NA, bwerm.sigma=NA,
														dbg.GTRparam=0, dbg.rER=0)						
		# standard run, fixed GTR param + fixed relative rates for each transmission chain
		# no WH Er to equal BH ER, and both being constant
		pipeline.vary	<- data.table(	wher.mu=c(log(0.00447743)-0.3^2/2, log(0.002239075)-0.13^2/2, log(0.002239075)-0.01^2/2), wher.sigma=c(0.3, 0.13, 0.01), 
										bwerm.mu=c(log(0.002239075)-0.13^2/2, log(0.002239075)-0.13^2/2, log(0.002239075)-0.01^2/2), bwerm.sigma=c(0.13, 0.13, 0.01), label=c('-d00st','-d00s','-d00f'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='wher.mu' ), 'v', as.character(wher.mu))
					set(pipeline.args, which( pipeline.args$stat=='wher.sigma' ), 'v', as.character(wher.sigma))
					set(pipeline.args, which( pipeline.args$stat=='bwerm.mu' ), 'v', as.character(bwerm.mu))
					set(pipeline.args, which( pipeline.args$stat=='bwerm.sigma' ), 'v', as.character(bwerm.sigma))
					
					print(pipeline.args)
					#	re-name the following:
					tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140920a'
					tmpdir			<- paste(tmpdir,label,sep='')
					dir.create(tmpdir, showWarnings=FALSE)						
					#						
					file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
					file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
					file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
					system(file)
				}, by='label']
		
	}
	#	debug 11, no import 
	if(0)
	{
		indir			<- system.file(package="rPANGEAHIVsim", "misc")
		infile.ind		<- '180914_HPTN071_scA_rep1'
		infile.trm		<- '180914_HPTN071_scA_rep1'			
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
														s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.1, 
														epi.model='HPTN071', epi.dt=1/48, epi.import=0,
														v.N0tau=1, v.r=2.851904, v.T50=-2,
														wher.mu=NA, wher.sigma=NA, bwerm.mu=NA, bwerm.sigma=NA,
														dbg.GTRparam=1, dbg.rER=1)						
		# standard run, fixed GTR param + fixed relative rates for each transmission chain
		# no WH Er to equal BH ER, and both being constant
		pipeline.vary	<- data.table(	wher.mu=c(log(0.00447743)-0.3^2/2, log(0.002239075)-0.13^2/2, log(0.002239075)-0.01^2/2), wher.sigma=c(0.3, 0.13, 0.01), 
										bwerm.mu=c(log(0.002239075)-0.13^2/2, log(0.002239075)-0.13^2/2, log(0.002239075)-0.01^2/2), bwerm.sigma=c(0.13, 0.13, 0.01), label=c('-i11st','-i11s','-i11f'))						
		dummy			<- pipeline.vary[, {				
												set(pipeline.args, which( pipeline.args$stat=='wher.mu' ), 'v', as.character(wher.mu))
												set(pipeline.args, which( pipeline.args$stat=='wher.sigma' ), 'v', as.character(wher.sigma))
												set(pipeline.args, which( pipeline.args$stat=='bwerm.mu' ), 'v', as.character(bwerm.mu))
												set(pipeline.args, which( pipeline.args$stat=='bwerm.sigma' ), 'v', as.character(bwerm.sigma))
												
												print(pipeline.args)
												#	re-name the following:
												tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140921b'
												tmpdir			<- paste(tmpdir,label,sep='')
												dir.create(tmpdir, showWarnings=FALSE)						
												#						
												file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
												file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
												file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
												system(file)
											}, by='label']
	}
	#	debug 11, no import, index.starttime.mode=fix
	if(0)
	{
		indir			<- system.file(package="rPANGEAHIVsim", "misc")
		infile.ind		<- '180914_HPTN071_scA_rep1'
		infile.trm		<- '180914_HPTN071_scA_rep1'			
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
														s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.1, 
														epi.model='HPTN071', epi.dt=1/48, epi.import=0,
														v.N0tau=1, v.r=2.851904, v.T50=-2,
														wher.mu=NA, wher.sigma=NA, bwerm.mu=NA, bwerm.sigma=NA,
														dbg.GTRparam=1, dbg.rER=1, index.starttime.mode='fix')						
		# standard run, fixed GTR param + fixed relative rates for each transmission chain
		# no WH Er to equal BH ER, and both being constant
		pipeline.vary	<- data.table(	wher.mu=c(log(0.00447743)-0.3^2/2, log(0.002239075)-0.13^2/2, log(0.002239075)-0.01^2/2), wher.sigma=c(0.3, 0.13, 0.01), 
				bwerm.mu=c(log(0.002239075)-0.13^2/2, log(0.002239075)-0.13^2/2, log(0.002239075)-0.01^2/2), bwerm.sigma=c(0.13, 0.13, 0.01), label=c('-S11st','-S11s','-S11f'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='wher.mu' ), 'v', as.character(wher.mu))
					set(pipeline.args, which( pipeline.args$stat=='wher.sigma' ), 'v', as.character(wher.sigma))
					set(pipeline.args, which( pipeline.args$stat=='bwerm.mu' ), 'v', as.character(bwerm.mu))
					set(pipeline.args, which( pipeline.args$stat=='bwerm.sigma' ), 'v', as.character(bwerm.sigma))
					
					print(pipeline.args)
					#	re-name the following:
					tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140921b'
					tmpdir			<- paste(tmpdir,label,sep='')
					dir.create(tmpdir, showWarnings=FALSE)						
					#						
					file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
					file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
					file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
					system(file)
				}, by='label']
	}
	#	debug 11, 10% import, index.starttime.mode=fix
	if(0)
	{
		indir			<- system.file(package="rPANGEAHIVsim", "misc")
		infile.ind		<- '180914_HPTN071_scA_rep1'
		infile.trm		<- '180914_HPTN071_scA_rep1'			
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
														s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.1, 
														epi.model='HPTN071', epi.dt=1/48, epi.import=0.1,
														v.N0tau=1, v.r=2.851904, v.T50=-2,
														wher.mu=NA, wher.sigma=NA, bwerm.mu=NA, bwerm.sigma=NA,
														dbg.GTRparam=1, dbg.rER=1, index.starttime.mode='fix')						
		# standard run, fixed GTR param + fixed relative rates for each transmission chain
		# no WH Er to equal BH ER, and both being constant
		pipeline.vary	<- data.table(	wher.mu=c(log(0.00447743)-0.3^2/2, log(0.002239075)-0.13^2/2, log(0.002239075)-0.01^2/2), wher.sigma=c(0.3, 0.13, 0.01), 
										bwerm.mu=c(log(0.002239075)-0.13^2/2, log(0.002239075)-0.13^2/2, log(0.002239075)-0.01^2/2), bwerm.sigma=c(0.13, 0.13, 0.01), label=c('-w11st','-w11s','-w11f'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='wher.mu' ), 'v', as.character(wher.mu))
					set(pipeline.args, which( pipeline.args$stat=='wher.sigma' ), 'v', as.character(wher.sigma))
					set(pipeline.args, which( pipeline.args$stat=='bwerm.mu' ), 'v', as.character(bwerm.mu))
					set(pipeline.args, which( pipeline.args$stat=='bwerm.sigma' ), 'v', as.character(bwerm.sigma))
					
					print(pipeline.args)
					#	re-name the following:
					tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140921b'
					tmpdir			<- paste(tmpdir,label,sep='')
					dir.create(tmpdir, showWarnings=FALSE)						
					#						
					file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
					file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
					file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
					system(file)
				}, by='label']
	}
	#	debug 11, no import, index.starttime.mode=fix, no rate variation
	if(0)
	{
		indir			<- system.file(package="rPANGEAHIVsim", "misc")
		infile.ind		<- '180914_HPTN071_scA_rep1'
		infile.trm		<- '180914_HPTN071_scA_rep1'			
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
														s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.1, 
														epi.model='HPTN071', epi.dt=1/48, epi.import=0,
														v.N0tau=1, v.r=2.851904, v.T50=-2,
														wher.mu=NA, wher.sigma=NA, bwerm.mu=NA, bwerm.sigma=NA, er.gamma=0,
														dbg.GTRparam=1, dbg.rER=1, index.starttime.mode='fix')						
		# standard run, fixed GTR param + fixed relative rates for each transmission chain
		# no WH Er to equal BH ER, and both being constant
		pipeline.vary	<- data.table(	wher.mu=c(log(0.00447743)-0.3^2/2, log(0.002239075)-0.13^2/2, log(0.002239075)-0.01^2/2), wher.sigma=c(0.3, 0.13, 0.01), 
										bwerm.mu=c(log(0.002239075)-0.13^2/2, log(0.002239075)-0.13^2/2, log(0.002239075)-0.01^2/2), bwerm.sigma=c(0.13, 0.13, 0.01), label=c('-R11st','-R11s','-R11f'))						
		dummy			<- pipeline.vary[, {				
												set(pipeline.args, which( pipeline.args$stat=='wher.mu' ), 'v', as.character(wher.mu))
												set(pipeline.args, which( pipeline.args$stat=='wher.sigma' ), 'v', as.character(wher.sigma))
												set(pipeline.args, which( pipeline.args$stat=='bwerm.mu' ), 'v', as.character(bwerm.mu))
												set(pipeline.args, which( pipeline.args$stat=='bwerm.sigma' ), 'v', as.character(bwerm.sigma))
												
												print(pipeline.args)
												#	re-name the following:
												tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140922'
												tmpdir			<- paste(tmpdir,label,sep='')
												dir.create(tmpdir, showWarnings=FALSE)						
												#						
												file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
												file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
												file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
												system(file)
											}, by='label']
	}
	#	debug 11, no import, index.starttime.mode=shift, no rate variation
	if(0)
	{
		indir			<- system.file(package="rPANGEAHIVsim", "misc")
		infile.ind		<- '180914_HPTN071_scA_rep1'
		infile.trm		<- '180914_HPTN071_scA_rep1'			
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
				s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.1, 
				epi.model='HPTN071', epi.dt=1/48, epi.import=0,
				v.N0tau=1, v.r=2.851904, v.T50=-2,
				wher.mu=NA, wher.sigma=NA, bwerm.mu=NA, bwerm.sigma=NA, er.gamma=0,
				dbg.GTRparam=1, dbg.rER=1, index.starttime.mode='shift')						
		# standard run, fixed GTR param + fixed relative rates for each transmission chain
		# no WH Er to equal BH ER, and both being constant
		pipeline.vary	<- data.table(	wher.mu=c(log(0.00447743)-0.3^2/2, log(0.002239075)-0.13^2/2, log(0.002239075)-0.01^2/2), wher.sigma=c(0.3, 0.13, 0.01), 
										bwerm.mu=c(log(0.002239075)-0.13^2/2, log(0.002239075)-0.13^2/2, log(0.002239075)-0.01^2/2), bwerm.sigma=c(0.13, 0.13, 0.01), label=c('-sh11st','-sh11s','-sh11f'))						
		dummy			<- pipeline.vary[, {				
												set(pipeline.args, which( pipeline.args$stat=='wher.mu' ), 'v', as.character(wher.mu))
												set(pipeline.args, which( pipeline.args$stat=='wher.sigma' ), 'v', as.character(wher.sigma))
												set(pipeline.args, which( pipeline.args$stat=='bwerm.mu' ), 'v', as.character(bwerm.mu))
												set(pipeline.args, which( pipeline.args$stat=='bwerm.sigma' ), 'v', as.character(bwerm.sigma))
												
												print(pipeline.args)
												#	re-name the following:
												tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140923'
												tmpdir			<- paste(tmpdir,label,sep='')
												dir.create(tmpdir, showWarnings=FALSE)						
												#						
												file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
												file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
												file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
												system(file)
											}, by='label']
	}
	#	debug 11, no import, same starting seq, index.starttime.mode=shift, no rate variation, 
	if(0)
	{
		indir			<- system.file(package="rPANGEAHIVsim", "misc")
		infile.ind		<- '180914_HPTN071_scA_rep1'
		infile.trm		<- '180914_HPTN071_scA_rep1'			
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
														s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.1, 
														epi.model='HPTN071', epi.dt=1/48, epi.import=0,
														v.N0tau=1, v.r=2.851904, v.T50=-2,
														wher.mu=NA, wher.sigma=NA, bwerm.mu=NA, bwerm.sigma=NA, er.gamma=0,
														dbg.GTRparam=1, dbg.rER=1, index.starttime.mode='shift', startseq.mode='fix')						
		# standard run, fixed GTR param + fixed relative rates for each transmission chain
		# no WH Er to equal BH ER, and both being constant
		pipeline.vary	<- data.table(	wher.mu=c(log(0.00447743)-0.3^2/2, log(0.002239075)-0.13^2/2, log(0.002239075)-0.01^2/2), wher.sigma=c(0.3, 0.13, 0.01), 
				bwerm.mu=c(log(0.002239075)-0.13^2/2, log(0.002239075)-0.13^2/2, log(0.002239075)-0.01^2/2), bwerm.sigma=c(0.13, 0.13, 0.01), label=c('-z11st','-z11s','-z11f'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='wher.mu' ), 'v', as.character(wher.mu))
					set(pipeline.args, which( pipeline.args$stat=='wher.sigma' ), 'v', as.character(wher.sigma))
					set(pipeline.args, which( pipeline.args$stat=='bwerm.mu' ), 'v', as.character(bwerm.mu))
					set(pipeline.args, which( pipeline.args$stat=='bwerm.sigma' ), 'v', as.character(bwerm.sigma))
					
					print(pipeline.args)
					#	re-name the following:
					tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140924'
					tmpdir			<- paste(tmpdir,label,sep='')
					dir.create(tmpdir, showWarnings=FALSE)						
					#						
					file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
					file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
					file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
					system(file)
				}, by='label']
		
		indir			<- system.file(package="rPANGEAHIVsim", "misc")
		infile.ind		<- '180914_HPTN071_scA_rep1'
		infile.trm		<- '180914_HPTN071_scA_rep1'			
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
														s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.1, 
														epi.model='HPTN071', epi.dt=1/48, epi.import=0,
														v.N0tau=1, v.r=2.851904, v.T50=-2,
														wher.mu=NA, wher.sigma=NA, bwerm.mu=NA, bwerm.sigma=NA, er.gamma=0,
														dbg.GTRparam=1, dbg.rER=1, index.starttime.mode='fix', startseq.mode='fix')						
		# standard run, fixed GTR param + fixed relative rates for each transmission chain
		# no WH Er to equal BH ER, and both being constant
		pipeline.vary	<- data.table(	wher.mu=c(log(0.00447743)-0.3^2/2, log(0.002239075)-0.13^2/2, log(0.002239075)-0.01^2/2), wher.sigma=c(0.3, 0.13, 0.01), 
										bwerm.mu=c(log(0.002239075)-0.13^2/2, log(0.002239075)-0.13^2/2, log(0.002239075)-0.01^2/2), bwerm.sigma=c(0.13, 0.13, 0.01), label=c('-f11st','-f11s','-f11f'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='wher.mu' ), 'v', as.character(wher.mu))
					set(pipeline.args, which( pipeline.args$stat=='wher.sigma' ), 'v', as.character(wher.sigma))
					set(pipeline.args, which( pipeline.args$stat=='bwerm.mu' ), 'v', as.character(bwerm.mu))
					set(pipeline.args, which( pipeline.args$stat=='bwerm.sigma' ), 'v', as.character(bwerm.sigma))
					
					print(pipeline.args)
					#	re-name the following:
					tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140924'
					tmpdir			<- paste(tmpdir,label,sep='')
					dir.create(tmpdir, showWarnings=FALSE)						
					#						
					file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
					file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
					file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
					system(file)
				}, by='label']
	}
	#	try and finalize:
	#	#startseq=1, starttime~1960, k=4 (did not make a big difference), standard run, import 2.5% 5% GTR 00/11
	if(0)
	{
		indir			<- system.file(package="rPANGEAHIVsim", "misc")
		infile.ind		<- '180914_HPTN071_scA_rep1'
		infile.trm		<- '180914_HPTN071_scA_rep1'			
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
														s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.1, 
														epi.model='HPTN071', epi.dt=1/48, epi.import=NA,
														v.N0tau=1, v.r=2.851904, v.T50=-2,
														wher.mu=log(0.00447743)-0.3^2/2, wher.sigma=0.3, bwerm.mu=log(0.002239075)-0.13^2/2, bwerm.sigma=0.13, er.gamma=4,
														dbg.GTRparam=1, dbg.rER=1, index.starttime.mode='fix', startseq.mode='fix')						
		# standard run, fixed GTR param + fixed relative rates for each transmission chain
		# no WH Er to equal BH ER, and both being constant
		pipeline.vary	<- data.table(	epi.import=c(0, 0.025, 0.05), label=c('-i110','-i112','-i115'))						
		dummy			<- pipeline.vary[, {				
												set(pipeline.args, which( pipeline.args$stat=='epi.import' ), 'v', as.character(epi.import))					
												print(pipeline.args)
												#	re-name the following:
												tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140925'
												tmpdir			<- paste(tmpdir,label,sep='')
												dir.create(tmpdir, showWarnings=FALSE)						
												#						
												file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
												file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
												file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
												system(file)
											}, by='label']
		#
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
														s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.1, 
														epi.model='HPTN071', epi.dt=1/48, epi.import=NA,
														v.N0tau=1, v.r=2.851904, v.T50=-2,
														wher.mu=log(0.00447743)-0.3^2/2, wher.sigma=0.3, bwerm.mu=log(0.002239075)-0.13^2/2, bwerm.sigma=0.13, er.gamma=4,
														dbg.GTRparam=0, dbg.rER=0, index.starttime.mode='fix', startseq.mode='fix')						
		# standard run, fixed GTR param + fixed relative rates for each transmission chain
		# no WH Er to equal BH ER, and both being constant
		pipeline.vary	<- data.table(	epi.import=c(0, 0.025, 0.05), label=c('-i000','-i002','-i005'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='epi.import' ), 'v', as.character(epi.import))					
					print(pipeline.args)
					#	re-name the following:
					tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140925'
					tmpdir			<- paste(tmpdir,label,sep='')
					dir.create(tmpdir, showWarnings=FALSE)						
					#						
					file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
					file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
					file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
					system(file)
				}, by='label']
		#	try improve TMRCA?
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
														s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.1, 
														epi.model='HPTN071', epi.dt=1/48, epi.import=NA,
														v.N0tau=1, v.r=2.851904, v.T50=-2,
														wher.mu=log(0.00447743)-0.3^2/2, wher.sigma=0.3, bwerm.mu=log(0.002239075)-0.13^2/2, bwerm.sigma=0.13, er.gamma=4,
														dbg.GTRparam=1, dbg.rER=1, index.starttime.mode='fix45', startseq.mode='fix')						
		# standard run, fixed GTR param + fixed relative rates for each transmission chain
		# no WH Er to equal BH ER, and both being constant
		pipeline.vary	<- data.table(	epi.import=c(0, 0.025, 0.05), label=c('-e110','-e112','-e115'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='epi.import' ), 'v', as.character(epi.import))					
					print(pipeline.args)
					#	re-name the following:
					tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140926'
					tmpdir			<- paste(tmpdir,label,sep='')
					dir.create(tmpdir, showWarnings=FALSE)						
					#						
					file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
					file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
					file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
					system(file)
				}, by='label']		
	}
	#	bugfix in SeqGen NEWICK file
	if(0)
	{
		indir			<- system.file(package="rPANGEAHIVsim", "misc")
		infile.ind		<- '180914_HPTN071_scA_rep1'
		infile.trm		<- '180914_HPTN071_scA_rep1'					
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
														s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.1, 
														epi.model='HPTN071', epi.dt=1/48, epi.import=NA,
														v.N0tau=1, v.r=2.851904, v.T50=-2,
														wher.mu=log(0.00447743)-0.3^2/2, wher.sigma=0.3, bwerm.mu=log(0.002239075)-0.13^2/2, bwerm.sigma=0.13, er.gamma=4,
														dbg.GTRparam=1, dbg.rER=1, index.starttime.mode='normal', startseq.mode='sample')						
		# standard run, fixed GTR param + fixed relative rates for each transmission chain
		# no WH Er to equal BH ER, and both being constant
		pipeline.vary	<- data.table(	epi.import=c(0, 0.05, 0.1), label=c('-n110','-n115','-n111'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='epi.import' ), 'v', as.character(epi.import))					
					print(pipeline.args)
					#	re-name the following:
					tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140926'
					tmpdir			<- paste(tmpdir,label,sep='')
					dir.create(tmpdir, showWarnings=FALSE)						
					#						
					file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
					file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
					file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
					system(file)
				}, by='label']
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
				s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.1, 
				epi.model='HPTN071', epi.dt=1/48, epi.import=NA,
				v.N0tau=1, v.r=2.851904, v.T50=-2,
				wher.mu=log(0.00447743)-0.3^2/2, wher.sigma=0.3, bwerm.mu=log(0.002239075)-0.13^2/2, bwerm.sigma=0.13, er.gamma=4,
				dbg.GTRparam=1, dbg.rER=1, index.starttime.mode='fix', startseq.mode='sample')						
		# standard run, fixed GTR param + fixed relative rates for each transmission chain
		# no WH Er to equal BH ER, and both being constant
		pipeline.vary	<- data.table(	epi.import=c(0, 0.05, 0.1), label=c('-f110','-f115','-f111'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='epi.import' ), 'v', as.character(epi.import))					
					print(pipeline.args)
					#	re-name the following:
					tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140926'
					tmpdir			<- paste(tmpdir,label,sep='')
					dir.create(tmpdir, showWarnings=FALSE)						
					#						
					file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
					file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
					file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
					system(file)
				}, by='label']
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
				s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.1, 
				epi.model='HPTN071', epi.dt=1/48, epi.import=NA,
				v.N0tau=1, v.r=2.851904, v.T50=-2,
				wher.mu=log(0.00447743)-0.3^2/2, wher.sigma=0.3, bwerm.mu=log(0.002239075)-0.13^2/2, bwerm.sigma=0.13, er.gamma=4,
				dbg.GTRparam=1, dbg.rER=1, index.starttime.mode='fix', startseq.mode='fix')						
		# standard run, fixed GTR param + fixed relative rates for each transmission chain
		# no WH Er to equal BH ER, and both being constant
		pipeline.vary	<- data.table(	epi.import=c(0, 0.05, 0.1), label=c('-s110','-s115','-s111'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='epi.import' ), 'v', as.character(epi.import))					
					print(pipeline.args)
					#	re-name the following:
					tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140926'
					tmpdir			<- paste(tmpdir,label,sep='')
					dir.create(tmpdir, showWarnings=FALSE)						
					#						
					file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
					file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
					file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
					system(file)
				}, by='label']
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
				s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.1, 
				epi.model='HPTN071', epi.dt=1/48, epi.import=NA,
				v.N0tau=1, v.r=2.851904, v.T50=-2,
				wher.mu=log(0.00447743)-0.3^2/2, wher.sigma=0.3, bwerm.mu=log(0.002239075)-0.13^2/2, bwerm.sigma=0.13, er.gamma=4,
				dbg.GTRparam=1, dbg.rER=1, index.starttime.mode='fix45', startseq.mode='fix')						
		# standard run, fixed GTR param + fixed relative rates for each transmission chain
		# no WH Er to equal BH ER, and both being constant
		pipeline.vary	<- data.table(	epi.import=c(0, 0.05, 0.1), label=c('-e110','-e115','-e111'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='epi.import' ), 'v', as.character(epi.import))					
					print(pipeline.args)
					#	re-name the following:
					tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140926'
					tmpdir			<- paste(tmpdir,label,sep='')
					dir.create(tmpdir, showWarnings=FALSE)						
					#						
					file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
					file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
					file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
					system(file)
				}, by='label']
	}
	#	try and finalize simulations
	if(0)
	{
		indir			<- system.file(package="rPANGEAHIVsim", "misc")
		infile.ind		<- '180914_HPTN071_scA_rep1'
		infile.trm		<- '180914_HPTN071_scA_rep1'	
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
				s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.1, 
				epi.model='HPTN071', epi.dt=1/48, epi.import=NA,
				v.N0tau=1, v.r=2.851904, v.T50=-2,
				wher.mu=log(0.00447743)-0.3^2/2, wher.sigma=0.3, bwerm.mu=log(0.002239075)-0.13^2/2, bwerm.sigma=0.13, er.gamma=0,
				dbg.GTRparam=1, dbg.rER=1, index.starttime.mode='fix', startseq.mode='sample')						
		# standard run, fixed GTR param + fixed relative rates for each transmission chain
		# no WH Er to equal BH ER, and both being constant
		pipeline.vary	<- data.table(	epi.import=c(0, 0.05, 0.1), label=c('-f0110','-f0115','-f0111'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='epi.import' ), 'v', as.character(epi.import))					
					print(pipeline.args)
					#	re-name the following:
					tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140927'
					tmpdir			<- paste(tmpdir,label,sep='')
					dir.create(tmpdir, showWarnings=FALSE)						
					#						
					file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
					file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
					file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
					system(file)
				}, by='label']
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
				s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.1, 
				epi.model='HPTN071', epi.dt=1/48, epi.import=NA,
				v.N0tau=1, v.r=2.851904, v.T50=-2,
				wher.mu=log(0.00447743)-0.3^2/2, wher.sigma=0.3, bwerm.mu=log(0.002239075)-0.13^2/2, bwerm.sigma=0.13, er.gamma=0,
				dbg.GTRparam=0, dbg.rER=0, index.starttime.mode='fix', startseq.mode='sample')						
		# standard run, fixed GTR param + fixed relative rates for each transmission chain
		# no WH Er to equal BH ER, and both being constant
		pipeline.vary	<- data.table(	epi.import=c(0, 0.05, 0.1), label=c('-f0000','-f0005','-f0001'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='epi.import' ), 'v', as.character(epi.import))					
					print(pipeline.args)
					#	re-name the following:
					tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140927'
					tmpdir			<- paste(tmpdir,label,sep='')
					dir.create(tmpdir, showWarnings=FALSE)						
					#						
					file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
					file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
					file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
					system(file)
				}, by='label']
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
				s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.1, 
				epi.model='HPTN071', epi.dt=1/48, epi.import=NA,
				v.N0tau=1, v.r=2.851904, v.T50=-2,
				wher.mu=log(0.00447743)-0.3^2/2, wher.sigma=0.3, bwerm.mu=log(0.002239075)-0.13^2/2, bwerm.sigma=0.13, er.gamma=4,
				dbg.GTRparam=0, dbg.rER=0, index.starttime.mode='fix', startseq.mode='sample')						
		# standard run, fixed GTR param + fixed relative rates for each transmission chain
		# no WH Er to equal BH ER, and both being constant
		pipeline.vary	<- data.table(	epi.import=c(0, 0.05, 0.1), label=c('-f4000','-f4005','-f4001'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='epi.import' ), 'v', as.character(epi.import))					
					print(pipeline.args)
					#	re-name the following:
					tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140927'
					tmpdir			<- paste(tmpdir,label,sep='')
					dir.create(tmpdir, showWarnings=FALSE)						
					#						
					file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
					file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
					file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
					system(file)
				}, by='label']
		
		
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
				s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.1, 
				epi.model='HPTN071', epi.dt=1/48, epi.import=NA,
				v.N0tau=1, v.r=2.851904, v.T50=-2,
				wher.mu=log(0.00447743)-0.3^2/2, wher.sigma=0.3, bwerm.mu=log(0.002239075)-0.13^2/2, bwerm.sigma=0.13, er.gamma=0,
				dbg.GTRparam=1, dbg.rER=1, index.starttime.mode='normal', startseq.mode='sample')						
		# standard run, fixed GTR param + fixed relative rates for each transmission chain
		# no WH Er to equal BH ER, and both being constant
		pipeline.vary	<- data.table(	epi.import=c(0, 0.05, 0.1), label=c('-n0110','-n0115','-n0111'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='epi.import' ), 'v', as.character(epi.import))					
					print(pipeline.args)
					#	re-name the following:
					tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140927'
					tmpdir			<- paste(tmpdir,label,sep='')
					dir.create(tmpdir, showWarnings=FALSE)						
					#						
					file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
					file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
					file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
					system(file)
				}, by='label']
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
				s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.1, 
				epi.model='HPTN071', epi.dt=1/48, epi.import=NA,
				v.N0tau=1, v.r=2.851904, v.T50=-2,
				wher.mu=log(0.00447743)-0.3^2/2, wher.sigma=0.3, bwerm.mu=log(0.002239075)-0.13^2/2, bwerm.sigma=0.13, er.gamma=0,
				dbg.GTRparam=0, dbg.rER=0, index.starttime.mode='normal', startseq.mode='sample')						
		# standard run, fixed GTR param + fixed relative rates for each transmission chain
		# no WH Er to equal BH ER, and both being constant
		pipeline.vary	<- data.table(	epi.import=c(0, 0.05, 0.1), label=c('-n0000','-n0005','-n0001'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='epi.import' ), 'v', as.character(epi.import))					
					print(pipeline.args)
					#	re-name the following:
					tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140927'
					tmpdir			<- paste(tmpdir,label,sep='')
					dir.create(tmpdir, showWarnings=FALSE)						
					#						
					file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
					file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
					file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
					system(file)
				}, by='label']
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
				s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.1, 
				epi.model='HPTN071', epi.dt=1/48, epi.import=NA,
				v.N0tau=1, v.r=2.851904, v.T50=-2,
				wher.mu=log(0.00447743)-0.3^2/2, wher.sigma=0.3, bwerm.mu=log(0.002239075)-0.13^2/2, bwerm.sigma=0.13, er.gamma=4,
				dbg.GTRparam=0, dbg.rER=0, index.starttime.mode='normal', startseq.mode='sample')						
		# standard run, fixed GTR param + fixed relative rates for each transmission chain
		# no WH Er to equal BH ER, and both being constant
		pipeline.vary	<- data.table(	epi.import=c(0, 0.05, 0.1), label=c('-n4000','-n4005','-n4001'))						
		dummy			<- pipeline.vary[, {				
					set(pipeline.args, which( pipeline.args$stat=='epi.import' ), 'v', as.character(epi.import))					
					print(pipeline.args)
					#	re-name the following:
					tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140927'
					tmpdir			<- paste(tmpdir,label,sep='')
					dir.create(tmpdir, showWarnings=FALSE)						
					#						
					file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
					file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
					file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
					system(file)
				}, by='label']
	}
	#	final three on all scenarios
	if(1)
	{
		indir			<- system.file(package="rPANGEAHIVsim", "misc")
		pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
														s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=NA, 
														epi.model='HPTN071', epi.dt=1/48, epi.import=0.05,
														v.N0tau=1, v.r=2.851904, v.T50=-2,
														wher.mu=log(0.00447743)-0.3^2/2, wher.sigma=0.3, bwerm.mu=log(0.002239075)-0.13^2/2, bwerm.sigma=0.13, er.gamma=NA,
														dbg.GTRparam=NA, dbg.rER=NA, index.starttime.mode=NA, startseq.mode='sample')						
		# standard run, fixed GTR param + fixed relative rates for each transmission chain
		# no WH Er to equal BH ER, and both being constant
		pipeline.vary	<- data.table(	dbg.GTRparam=c(1,1,0), 
										dbg.rER=c(1,1,0),
										index.starttime.mode=c('fix','normal','normal'),
										er.gamma=c(4,4,0),
										label=c('-f4115','-n4115','-n0005'))						
		dummy			<- pipeline.vary[, {				
												set(pipeline.args, which( pipeline.args$stat=='dbg.GTRparam' ), 'v', as.character(dbg.GTRparam))
												set(pipeline.args, which( pipeline.args$stat=='dbg.rER' ), 'v', as.character(dbg.rER))
												set(pipeline.args, which( pipeline.args$stat=='index.starttime.mode' ), 'v', as.character(index.starttime.mode))
												set(pipeline.args, which( pipeline.args$stat=='er.gamma' ), 'v', as.character(er.gamma))
												print(pipeline.args)
												#	scenario A						
												infile.ind		<- '180914_HPTN071_scA_rep1'
												infile.trm		<- '180914_HPTN071_scA_rep1'
												tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140928a'
												tmpdir			<- paste(tmpdir,label,sep='')
												dir.create(tmpdir, showWarnings=FALSE)																		
												set(pipeline.args, which( pipeline.args$stat=='s.PREV.max' ), 'v', '0.11')
												file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
												file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
												file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
												system(file)
												#	scenario B						
												infile.ind		<- '180914_HPTN071_scB_rep1'
												infile.trm		<- '180914_HPTN071_scB_rep1'	
												tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140928b'
												tmpdir			<- paste(tmpdir,label,sep='')
												dir.create(tmpdir, showWarnings=FALSE)																														
												set(pipeline.args, which( pipeline.args$stat=='s.PREV.max' ), 'v', '0.15')
												file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
												file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
												file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
												system(file)
												#	scenario C						
												infile.ind		<- '180914_HPTN071_scC_rep1'
												infile.trm		<- '180914_HPTN071_scC_rep1'	
												tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140928c'
												tmpdir			<- paste(tmpdir,label,sep='')
												dir.create(tmpdir, showWarnings=FALSE)																														
												set(pipeline.args, which( pipeline.args$stat=='s.PREV.max' ), 'v', '0.185')
												file.copy(paste(indir,'/',infile.ind,'_IND.csv',sep=''), paste(tmpdir,'/',infile.ind,label,'_IND.csv',sep=''))
												file.copy(paste(indir,'/',infile.trm,'_TRM.csv',sep=''), paste(tmpdir,'/',infile.trm,label,'_TRM.csv',sep=''))
												file			<- rPANGEAHIVsim.pipeline(tmpdir, paste(infile.ind,label,'_IND.csv',sep=''), paste(infile.trm,label,'_TRM.csv',sep=''), tmpdir, pipeline.args=pipeline.args)
												system(file)
											}, by='label']
	}
	
}
##--------------------------------------------------------------------------------------------------------
##	check simulated sequences: create ExaML tree and estimate R2
##	olli 02.02.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.SSApg.CLUSTERBEAST.skygrid.hky<- function()
{
	require(rBEAST)
	
	tree.id.labelsep		<- '|'
	tree.id.label.idx.ctime	<- 4 	
	#select		<- 'grid-mseq500'
	#select		<- 'grid-cseq3'
	select		<- 'grid-mseq400'
	indir		<- '/Users/Oliver/duke/2014_Gates/methods_comparison_pipeline/150414'
	#indir		<- '/Users/Oliver/git/HPTN071sim/tmp140914/140716_RUN001_INTERNAL'  
	outdir		<- indir
	infiles		<- list.files(indir, '.*INTERNAL.R$', full.names=FALSE)
	#stopifnot(length(infiles)==1)
	#
	#	run  
	#	
	selects		<- c( paste('grid-mseq',seq(600, 1200, 200), sep=''), paste('grid-clsmseq',seq(600, 1200, 200), sep=''), paste('grid-clrndseq',seq(600, 1200, 200), sep=''))
	for(select in selects)
	{
		for(i in seq_along(infiles))
		{
			infile			<- infiles[i]
			#	load simulated data
			file			<- paste(indir, '/', infile, sep='')
			file.name		<- paste(indir, gsub('_SIMULATED_INTERNAL.R',paste('_TEST_pol_HKY_fixedtree_',select,'.xml',sep=''),infile), sep='/')
			cat(paste('\nLoading file', file))
			load(file)		#expect "df.epi"    "df.trms"   "df.inds"   "df.sample" "df.seq"
			set( df.seq, NULL, 'IDCLU', df.seq[, as.integer(IDCLU)] )
			setnames(df.seq, c("LABEL", "IDCLU", "IDPOP"), c("TAXON_NAME", "CLU_ID", "TAXON_ID"))
			#	
			seq.select		<- beast.choose.seq.by.clusters(df.seq, select, verbose=1)
			#
			#	read NEWICK trees for each cluster phylogeny, if there
			#
			phd		<- NULL
			tmp		<- list.files(indir, '_DATEDTREE.newick$', full.names=FALSE)
			tmp		<- tmp[ grepl(substr(infile, 1, regexpr('_SIMULATED',infile)), tmp) ]
			if(length(tmp))
			{
				# select
				phd					<- read.tree(paste(indir, tmp, sep='/'))
				
				tmp2				<- data.table(TAXON_ID= sapply(phd, function(x) x$tip.label[1]), IDX=seq_along(phd))
				set( tmp2, NULL, 'TAXON_ID', tmp2[, as.integer(substring(sapply(strsplit(TAXON_ID, tree.id.labelsep, fixed=TRUE),'[[',1),7)) ] )
				tmp2				<- merge(subset(seq.select, select=c(TAXON_ID, CLU_ID)), tmp2, by='TAXON_ID')			
				phd					<- lapply(tmp2[,IDX], function(i) phd[[i]] )
				names(phd)			<- tmp2[, CLU_ID]
				# plot
				phd.plot			<- eval(parse(text=paste('phd[[',seq_along(phd),']]', sep='',collapse='+')))			
				#phd.plot			<- drop.tip(phd.plot, which(grepl('NOEXIST', phd.plot$tip.label)), root.edge=1)
				phd.plot			<- ladderize(phd.plot)
				tmp					<- paste(indir, '/', gsub('DATEDTREE','BEASTDATEDTREE',tmp), sep='')
				pdf(file=gsub('.xml','_BEASTDATEDTREE.pdf',file.name), w=10, h=Ntip(phd.plot)*0.1)
				plot(phd.plot, show.tip=TRUE, cex=0.5)
				dev.off()									
			}
			#
			#	create BEAST XML POL
			#
			seq.select.pol	<- subset(seq.select, select=c("CLU_ID", "TAXON_ID", "TAXON_NAME", "POL" ))
			setnames(seq.select.pol, 'POL', 'SEQ')									
			bxml			<- beastscript.multilocus.hky( file.name, seq.select.pol, phd, verbose=1 )
			cat(paste("\nwrite xml file to",file.name))
			saveXML(bxml, file=file.name)				
		}
	}
}
##--------------------------------------------------------------------------------------------------------
##	check simulated sequences: create ExaML tree and estimate R2
##	olli 02.02.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.SSApg.CLUSTERBEAST.skygrid.codon.gtr<- function()
{
	require(rBEAST)
	
	tree.id.labelsep		<- '|'
	tree.id.label.idx.ctime	<- 4 	
	#select		<- 'grid-mseq500'
	select		<- 'grid-cseq3'
	select		<- 'grid-mseq400'
	indir		<- '/Users/Oliver/duke/2014_Gates/methods_comparison_pipeline/150414'
	#indir		<- '/Users/Oliver/git/HPTN071sim/tmp140914/140716_RUN001_INTERNAL'  
	outdir		<- indir
	infiles		<- list.files(indir, '.*INTERNAL.R$', full.names=FALSE)
	#stopifnot(length(infiles)==1)
	#
	#	run  
	#	
	selects		<- c( paste('grid-mseq',seq(600, 1200, 200), sep=''), paste('grid-clsmseq',seq(600, 1200, 200), sep=''), paste('grid-clrndseq',seq(600, 1200, 200), sep=''))
	selects		<- paste('grid-mseq',seq(600, 1200, 200), sep='')
	for(select in selects)
	{		
		for(i in seq_along(infiles))
		{
			infile			<- infiles[i]
			#	load simulated data
			file			<- paste(indir, '/', infile, sep='')
			cat(paste('\nLoading file', file))
			load(file)		#expect "df.epi"    "df.trms"   "df.inds"   "df.sample" "df.seq"
			file.name		<- paste(indir, gsub('_SIMULATED_INTERNAL.R',paste('_TEST_pol_CODON-GTR_fixedtree_',select,'.xml',sep=''),infile), sep='/')
			set( df.seq, NULL, 'IDCLU', df.seq[, as.integer(IDCLU)] )
			setnames(df.seq, c("LABEL", "IDCLU", "IDPOP"), c("TAXON_NAME", "CLU_ID", "TAXON_ID"))
			#	
			seq.select		<- beast.choose.seq.by.clusters(df.seq, select, verbose=1)
			#
			#	read NEWICK trees for each cluster phylogeny, if there
			#
			phd		<- NULL
			tmp		<- list.files(indir, '_DATEDTREE.newick$', full.names=FALSE)
			tmp		<- tmp[ grepl(substr(infile, 1, regexpr('_SIMULATED',infile)), tmp) ]
			if(length(tmp))
			{
				# select
				phd					<- read.tree(paste(indir, tmp, sep='/'))
				
				tmp2				<- data.table(TAXON_ID= sapply(phd, function(x) x$tip.label[1]), IDX=seq_along(phd))
				set( tmp2, NULL, 'TAXON_ID', tmp2[, as.integer(substring(sapply(strsplit(TAXON_ID, tree.id.labelsep, fixed=TRUE),'[[',1),7)) ] )
				tmp2				<- merge(subset(seq.select, select=c(TAXON_ID, CLU_ID)), tmp2, by='TAXON_ID')			
				phd					<- lapply(tmp2[,IDX], function(i) phd[[i]] )
				names(phd)			<- tmp2[, CLU_ID]
				# plot
				phd.plot			<- eval(parse(text=paste('phd[[',seq_along(phd),']]', sep='',collapse='+')))			
				#phd.plot			<- drop.tip(phd.plot, which(grepl('NOEXIST', phd.plot$tip.label)), root.edge=1)
				phd.plot			<- ladderize(phd.plot)									
				pdf(file=gsub('.xml','_BEASTDATEDTREE.pdf',file.name), w=10, h=Ntip(phd.plot)*0.1)
				plot(phd.plot, show.tip=TRUE, cex=0.5)
				dev.off()									
			}
			#
			#	create BEAST XML POL
			#
			seq.select.pol	<- subset(seq.select, select=c("CLU_ID", "TAXON_ID", "TAXON_NAME", "POL" ))
			setnames(seq.select.pol, 'POL', 'SEQ')							
			bxml			<- beastscript.multilocus.codon.gtr( file.name, seq.select.pol, phd, verbose=1 )
			cat(paste("\nwrite xml file to",file.name))
			saveXML(bxml, file=file.name)		
		}
	}
}
##--------------------------------------------------------------------------------------------------------
##	check simulated sequences: create ExaML tree and estimate R2
##	olli 04.02.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.SeqNo<- function()
{
	indir			<- '/Users/Oliver/duke/2014_Gates/methods_comparison_pipeline/150205'
	dfi				<- data.table(FILE=list.files(indir, '.*SIMULATED_INTERNAL.R$', full.names=FALSE))
	dfi[, SC:= sapply(strsplit(FILE, '_'),'[[',3)]
	dfi[, CONFIG:= sapply(strsplit(SC, '-'),'[[',2)]
	set(dfi, NULL, 'SC', dfi[, sapply(strsplit(SC, '-'),'[[',1)])	
	dfs				<- dfi[,{
								file			<- paste(indir, FILE, sep='/' )
								cat(paste('\nprocess file=',file))
								load(file)
								df.inds
								tmp2				<- subset( df.inds, DOD>pipeline.args['yr.end',][, as.numeric(v)] & floor(TIME_TR)<pipeline.args['yr.end',][, as.numeric(v)])
								sc.alive20.infl20	<- tmp2[, length(which(!is.na(TIME_SEQ)))/length(TIME_SEQ) ]					
								tmp2				<- subset( df.inds, floor(TIME_TR)>=2000 & floor(TIME_TR)<pipeline.args['yr.end',][, as.numeric(v)])
								sc.infg99l20		<- tmp2[, length(which(!is.na(TIME_SEQ)))/length(TIME_SEQ) ]				
								sn.g14				<- subset( df.inds, TIME_SEQ>=pipeline.args['s.INTERVENTION.start',][, as.numeric(v)] & floor(TIME_TR)<pipeline.args['yr.end',][, as.numeric(v)])[, length(which(!is.na(TIME_SEQ))) ] 
								sn.total			<- subset(df.inds, floor(TIME_TR)<pipeline.args['yr.end',][, as.numeric(v)])[, length(which(!is.na(TIME_SEQ))) ]
								list(SN_TOTAL=sn.total, SN_G14=sn.g14, SC_20=sc.alive20.infl20, SC_9920=sc.infg99l20)								
							}, by=c('SC','CONFIG')]
	setkey(dfs, CONFIG, SC)
	dfs[, PC_G14:= SN_G14/SN_TOTAL]
}
##--------------------------------------------------------------------------------------------------------
##	check simulated sequences: create ExaML tree and estimate R2
##	olli 04.02.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.Surveys<- function()
{
	indir			<- '/Users/Oliver/duke/2014_Gates/methods_comparison_pipeline/150205'
	dfi				<- data.table(FILE=list.files(indir, '.*SIMULATED_INTERNAL.R$', full.names=FALSE))
	dfi[, SC:= sapply(strsplit(FILE, '_'),'[[',3)]
	dfi[, CONFIG:= sapply(strsplit(SC, '-'),'[[',2)]
	set(dfi, NULL, 'SC', dfi[, sapply(strsplit(SC, '-'),'[[',1)])
	#
	#	complete census
	#
	df.sp			<- dfi[,{
								file			<- paste(indir, FILE, sep='/' )
								cat(paste('\nprocess file=',file))
								load(file)
								df.sp	<- PANGEA.SeroPrevalenceSurvey(df.inds, epi.adult=13, s.INTERVENTION.start=2015, sp.prop.of.sexactive=1, sp.times=c(12, 8, 4, 0), test=1 )								
								df.sp[, list(ALIVE_N=sum(ALIVE_N), ALIVE_AND_HIV_N=sum(ALIVE_AND_HIV_N), ALIVE_AND_DIAG_N=sum(ALIVE_AND_DIAG_N), ALIVE_AND_ART_N=sum(ALIVE_AND_ART_N), ALIVE_AND_SEQ_N=sum(ALIVE_AND_SEQ_N)), by='YR']				
							}, by=c('SC','CONFIG')]
	df.sp[, DIAG_PC:= ALIVE_AND_DIAG_N/ALIVE_AND_HIV_N]
	df.sp[, ART_PC:= ALIVE_AND_ART_N/ALIVE_AND_HIV_N]
	df.sp[, SEQ_PC:= ALIVE_AND_SEQ_N/ALIVE_AND_HIV_N]
	df.sp			<- melt( df.sp, id.vars=c('SC','CONFIG','YR'), measure.vars=c('DIAG_PC','ART_PC','SEQ_PC') )	
	ggplot(df.sp, aes(x=YR, y=value, colour=SC, group=SC)) + geom_point() + geom_line() +
			scale_color_brewer(name='scenario', palette='Paired') + 
			scale_x_continuous(breaks=df.sp[, unique(YR)], expand=c(.2,.2)) +
			facet_wrap(CONFIG~variable, scales='free_y', ncol=3) +
			theme_bw() + labs(x='', y='proportion') 
	file	<- paste(indir, '150205_TEST_SurveyComplete.pdf')
	ggsave(file=file, w=12, h=6)
	#
	#	survey on 5%
	#
	df.sp			<- dfi[,{
				file			<- paste(indir, FILE, sep='/' )
				cat(paste('\nprocess file=',file))
				load(file)
				df.sp									
				df.sp[, list(ALIVE_N=sum(ALIVE_N), ALIVE_AND_DIAG_N=sum(ALIVE_AND_DIAG_N), ALIVE_AND_ART_N=sum(ALIVE_AND_ART_N), ALIVE_AND_SEQ_N=sum(ALIVE_AND_SEQ_N)), by='YR'] 
			}, by=c('SC','CONFIG')]
	df.sp			<- melt( df.sp, id.vars=c('SC','CONFIG','YR'), measure.vars=c('ALIVE_N','ALIVE_AND_DIAG_N','ALIVE_AND_ART_N','ALIVE_AND_SEQ_N') )
	ggplot(df.sp, aes(x=YR, y=value, colour=SC, group=SC)) + geom_point() + geom_line() +
			scale_color_brewer(name='scenario', palette='Paired') + 
			scale_x_continuous(breaks=df.sp[, unique(YR)], expand=c(.2,.2)) +
			facet_wrap(CONFIG~variable, scales='free_y', ncol=4) +
			theme_bw() + labs(x='', y='proportion') 
	file	<- paste(indir, '150205_TEST_SurveySampled.pdf')
	ggsave(file=file, w=14, h=7)
}
##--------------------------------------------------------------------------------------------------------
##	check simulated sequences: create ExaML tree and estimate R2
##	olli 03.02.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.PropCD4atDiag<- function()
{	
	#check time to coalescence in true dated tree
	label.sep		<- '|'
	label.idx.date	<- 4
	label.idx.idpop	<- 1
	indir			<- '/Users/Oliver/duke/2014_Gates/methods_comparison_pipeline/150205'
	dfi				<- data.table(FILE=list.files(indir, '.*SIMULATED_INTERNAL.R$', full.names=FALSE))
	dfi[, SC:= sapply(strsplit(FILE, '_'),'[[',3)]
	dfi[, CONFIG:= sapply(strsplit(SC, '-'),'[[',2)]
	set(dfi, NULL, 'SC', dfi[, sapply(strsplit(SC, '-'),'[[',1)])	
	
	df				<- dfi[, {
								file			<- paste(indir, FILE, sep='/' )
								cat(paste('\nprocess',file))
								load(file)
								ans				<- subset(df.inds, !is.na(DIAG_T) & DIAG_T<2020, c(TIME_SEQ, RECENT_TR, DIAG_T, DIAG_CD4 ))
								tmp				<- sample(seq_len(nrow(ans)), round(nrow(ans)*.5))
								ans[, RECENT_TRr:= NA_character_]
								set(ans, tmp, 'RECENT_TRr', ans[tmp, as.character(RECENT_TR)])
								ans
							}, by=c('SC', 'CONFIG')]
	df[, DIAG_CD4c:= df[, cut(DIAG_CD4, breaks=c(0,200,350,500,700,2000))]]
	set( df, NULL, 'RECENT_TRr', df[, factor(RECENT_TRr)] )
	#	plot CD4 counts at time of diagnosis
	df	<- subset(df, CONFIG%in%c('mFP25','mFP50','mFP85'))
	ggplot(df, aes(x=floor(DIAG_T), fill=DIAG_CD4c)) + geom_bar(binwidth=1, position='fill', alpha=0.5) + 
			labs(fill='CD4 category', y='percent of diagnosed', x='year diagnosed') + scale_y_continuous(breaks=seq(0,1,0.1)) + scale_fill_brewer(palette='Set1') +
			theme_bw() + theme(panel.grid.major.y=element_line(colour="black", size=0.2)) + facet_grid(SC~CONFIG)
	file	<- paste(indir, '150205_TEST_CD4atDiagnosis.pdf')
	ggsave(file=file, w=12, h=12)
	
	ggplot(subset(df, !is.na(TIME_SEQ)), aes(x=floor(DIAG_T), fill=DIAG_CD4c)) + geom_bar(binwidth=1, position='fill', alpha=0.5) + 
			labs(fill='CD4 category', y='percent of diagnosed\namong sequenced', x='year diagnosed') + scale_y_continuous(breaks=seq(0,1,0.1)) + scale_fill_brewer(palette='Set1') +
			theme_bw() + theme(panel.grid.major.y=element_line(colour="black", size=0.2)) + facet_grid(SC~CONFIG)
	file	<- paste(indir, '150205_TEST_CD4atDiagnosisAmongSequenced.pdf')
	ggsave(file=file, w=12, h=12)
	
	#	plot recent at time of diagnosis
	ggplot(df, aes(x=floor(DIAG_T), fill=RECENT_TR)) + geom_bar(binwidth=1, position='fill', alpha=0.5) + 
			labs(fill='Infected in\nlast 6 months', y='percent of diagnosed', x='year diagnosed') + scale_y_continuous(breaks=seq(0,1,0.1)) + scale_fill_brewer(palette='Set1') +
			theme_bw() + theme(panel.grid.major.y=element_line(colour="black", size=0.2)) + facet_grid(SC~CONFIG)
	file	<- paste(indir, '150205_TEST_RECENT.pdf')
	ggsave(file=file, w=12, h=12)
	
	ggplot(subset(df, !is.na(TIME_SEQ) & !is.na(RECENT_TRr)), aes(x=floor(DIAG_T), fill=RECENT_TRr)) + geom_bar(binwidth=1, position='fill', alpha=0.5) + 
			labs(fill='Infected in\nlast 6 months', y='percent of diagnosed\nfor whom recency known and sequenced', x='year diagnosed') + scale_y_continuous(breaks=seq(0,1,0.1)) + scale_fill_brewer(palette='Set1') +
			theme_bw() + theme(panel.grid.major.y=element_line(colour="black", size=0.2)) + facet_grid(SC~CONFIG)
	file	<- paste(indir, '150205_TEST_RECENTAmongSequenced50.pdf')
	ggsave(file=file, w=12, h=12)
	
}
##--------------------------------------------------------------------------------------------------------
##	check simulated sequences: create ExaML tree and estimate R2
##	olli 03.02.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.CoClusteringAcute<- function()
{
	label.sep		<- '|'
	label.idx.date	<- 4
	label.idx.idpop	<- 1
	indir			<- '/Users/Oliver/duke/2014_Gates/methods_comparison_pipeline/150205'
	dfi				<- data.table(FILE=list.files(indir, '.*SIMULATED_INTERNAL.R$', full.names=FALSE))
	dfi[, SC:= sapply(strsplit(FILE, '_'),'[[',3)]
	dfi[, CONFIG:= sapply(strsplit(SC, '-'),'[[',2)]
	set(dfi, NULL, 'SC', dfi[, sapply(strsplit(SC, '-'),'[[',1)])
	#
	#	true phylo
	#
	clu				<- dfi[, 	{
									#FILE<- dfi[2,FILE]
									file	<- paste(indir, FILE, sep='/' )
									cat(paste('\nprocess',file))										
									load( file )
									#
									clu		<- subset(df.inds, !is.na(TIME_SEQ))[, list(RECENT_Y=sum(length(which(RECENT_TR=='Y'))), SIZE=length(IDPOP) ), by='IDCLU']
									subset(clu, SIZE>1)
								}, by=c('SC','CONFIG')]
	clu[, SIZEc:= cut(SIZE, breaks=c(0,4,10,1e3), labels=c('<=4','5-10','>10'))]					
	clup			<- clu[, list(RECENT_P=mean(RECENT_Y/SIZE)), by=c('SC','CONFIG','SIZEc')]
	ggplot(clup, aes(x=SIZEc, y=RECENT_P, group=SC, colour=SC)) + geom_point() + geom_line() + facet_grid(.~CONFIG) + 
			theme_bw() + scale_colour_brewer(palette='Paired', name='scenario') +
			labs(x='Number of seqs in true dated cluster phylogeny', y='individuals within 6 mo of HIV infection\n(%)')
	file	<- paste(indir, '150205_TEST_ClustWithin6Mo_truephylo.pdf')
	ggsave(file=file, w=8, h=6)
	
	ggplot(subset(clup, grepl('s2x',CONFIG)), aes(x=SIZEc, y=RECENT_P, group=SC, colour=SC)) + geom_point() + geom_line() + facet_grid(.~CONFIG) + 
			theme_bw() + scale_colour_brewer(palette='Paired', name='scenario') +
			labs(x='Number of seqs in true dated cluster phylogeny', y='individuals within 6 mo of HIV infection\n(%)')
	
	
	#
	#	NJ tree
	#
	#	load outgroup sequences
	file			<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfg_HXB2outgroup.R')
	cat(paste('\nLoading outgroup seq from file', file))
	load(file)		#expect "outgroup.seq.gag" "outgroup.seq.pol" "outgroup.seq.env"
	
	clu				<- dfi[, 	{				
				file	<- paste(indir, FILE, sep='/' )
				cat(paste('\nprocess',file))										
				load( file )
				#	concatenate sequences	
				tmp				<- tolower(do.call('rbind',strsplit(df.seq[, paste(GAG,POL,ENV,sep='')],'')))		
				rownames(tmp)	<- df.seq[, paste(IDCLU,label.sep,LABEL,sep='')]	
				seq				<- as.DNAbin(tmp)
				tmp				<- cbind(outgroup.seq.gag[,seq_len(df.seq[1, nchar(GAG)])], outgroup.seq.pol, outgroup.seq.env)
				seq				<- rbind(seq,tmp)
				#	NJ tree
				seq.ph			<- nj(dist.dna(seq, model='raw'))		
				tmp				<- which(seq.ph$tip.label=="HXB2")
				seq.ph			<- reroot(seq.ph, tmp, seq.ph$edge.length[which(seq.ph$edge[,2]==tmp)])
				seq.ph			<- ladderize(seq.ph)
				#	clustering at 4%
				seq.brdist		<- hivc.clu.brdist.stats(seq.ph, eval.dist.btw="leaf", stat.fun=hivc.clu.min.transmission.cascade)
				seq.clu			<- lapply(c(0.04, 0.06, 0.08), function(x)
						{
							clustering		<- hivc.clu.clusterbythresh(seq.ph, thresh.brl=x, dist.brl=seq.brdist, retval="all")
							seq.clu		<- subset( data.table( PH_NODE_ID=seq_len(Ntip(seq.ph)), CLU_ID=clustering$clu.mem[ seq_len(Ntip(seq.ph))] ), !is.na(CLU_ID) )
							seq.clu[, LABEL:= seq.clu[, seq.ph$tip.label[PH_NODE_ID]] ]
							seq.clu[, PH_NODE_ID:=NULL]
							seq.clu[, BRL:=x]
							seq.clu
						})
				seq.clu			<- do.call('rbind',seq.clu)
				seq.clu[, IDPOP:=seq.clu[, as.integer(substring(sapply(strsplit(LABEL,label.sep,fixed=TRUE),'[[',2),7))]]
				
				seq.clu			<- merge(seq.clu, subset(df.inds, !is.na(TIME_SEQ), c(IDPOP, RECENT_TR)), by='IDPOP')				
				seq.clu			<- seq.clu[, list(RECENT_Y=sum(length(which(RECENT_TR=='Y'))), SIZE=length(IDPOP) ), by=c('CLU_ID','BRL')]
				setnames(seq.clu, 'CLU_ID','IDPHCLU')
				seq.clu
			}, by=c('SC','CONFIG')]
	file	<- paste(indir, '150205_TEST_ClustWithin6Mo_clusters.R')
	save(file=file, clu)
	clu[, SIZEc:= cut(SIZE, breaks=c(0,4,10,1e3), labels=c('<=4','5-10','>10'))]					
	clup			<- clu[, list(RECENT_P=mean(RECENT_Y/SIZE)), by=c('SC','CONFIG','SIZEc','BRL')]
	ggplot(clup, aes(x=SIZEc, y=RECENT_P, group=SC, colour=SC)) + geom_point() + geom_line() + facet_grid(BRL~CONFIG) + 
			theme_bw() + scale_colour_brewer(palette='Paired', name='scenario') +
			labs(x='Number of seqs in true dated cluster phylogeny', y='individuals within 6 mo of HIV infection\n(%)')
	file	<- paste(indir, '150205_TEST_ClustWithin6Mo_clusters.pdf')
	ggsave(file=file, w=8, h=10)
	
}
##--------------------------------------------------------------------------------------------------------
##	check reproducibility
##	olli 03.02.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.Reproducible<- function()
{	
	#require(devtools)
	#devtools::install("/Users/Oliver/git/HPTN071sim/source/rPANGEAHIVsim")
	#check time to coalescence in true dated tree
	label.sep		<- '|'
	label.idx.date	<- 4
	label.idx.idpop	<- 1
	
	file		<- "/Users/Oliver/git/HPTN071sim/tmp150227-E-o111/150129_HPTN071_scE-o111_INTERNAL/150129_HPTN071_scE-o111_SIMULATED_INTERNAL.R"
	load(file)
	df.epi2		<- copy(df.epi)    
	df.trms2	<- copy(df.trms)
	df.inds2	<- copy(df.inds)   
	df.sample2	<- copy(df.sample)
	df.seq2		<- copy(df.seq)
	df.sp2		<- copy(df.sp)	
	file		<-  "/Users/Oliver/git/HPTN071sim/tmp150227-E-o111-check/150129_HPTN071_scE-o111_INTERNAL/150129_HPTN071_scE-o111_SIMULATED_INTERNAL.R"
	load(file)
	
	all(df.epi[, PREV]==df.epi2[, PREV])
	all(df.epi[, IMPORT]==df.epi2[, IMPORT])
	all(df.sample[, s.nTOTAL]==df.sample2[, s.nTOTAL])
	all(subset(df.inds, !is.na(TIME_SEQ))[, TIME_SEQ]==subset(df.inds2, !is.na(TIME_SEQ))[, TIME_SEQ])
	all(subset(df.inds, !is.na(TIME_SEQ))[, IDCLU]==subset(df.inds2, !is.na(TIME_SEQ))[, IDCLU])	
	all(df.sp[, ALIVE_N]==df.sp2[, ALIVE_N])
	all(df.sp[, ALIVE_AND_SEQ_N]==df.sp2[, ALIVE_AND_SEQ_N])
	
	all(df.seq[, IDPOP]==df.seq2[, IDPOP])
	all(df.seq[, GAG]==df.seq2[, GAG])
	all(df.seq[, POL]==df.seq2[, POL])
	all(df.seq[, ENV]==df.seq2[, ENV])
	#
	#
	#
	file1	<- '/Users/Oliver/git/HPTN071sim/tmp150227-E-o111 copy/150129_HPTN071_scE-o111_INTERNAL/150129_HPTN071_scE-o111_SIMULATED_INTERNAL.R'
	file2	<- '/Users/Oliver/git/HPTN071sim/tmp150227-E-o111/150129_HPTN071_scE-o111_INTERNAL/150129_HPTN071_scE-o111_SIMULATED_INTERNAL.R'
	load(file2)
	df.epi2		<- copy(df.epi)    
	df.trms2	<- copy(df.trms)
	df.inds2	<- copy(df.inds)   
	df.sample2	<- copy(df.sample)
	df.seq2		<- copy(df.seq)
	df.sp2		<- copy(df.sp)		
	load(file1)
	all(df.seq[, IDPOP]==df.seq2[, IDPOP])
	all(df.seq[, GAG]==df.seq2[, GAG])
	all(df.seq[, POL]==df.seq2[, POL])
	all(df.seq[, ENV]==df.seq2[, ENV])
	
	ph1	<- read.tree('/Users/Oliver/git/HPTN071sim/tmp150227-E-o111-check/150129_HPTN071_scE-o111_SIMULATED_TREE/150129_HPTN071_scE-o111_DATEDTREE.newick')
	ph2	<- read.tree('/Users/Oliver/git/HPTN071sim/tmp150227-E-o111/150129_HPTN071_scE-o111_SIMULATED_TREE/150129_HPTN071_scE-o111_DATEDTREE.newick')
	setequal(ph1$tip.label,ph2$tip.label)
	tmp	<- merge( 	data.table(LABEL=ph1$tip.label, H1=node.depth.edgelength(ph1)[seq_len(Ntip(ph1))]), 
					data.table(LABEL=ph2$tip.label, H2=node.depth.edgelength(ph2)[seq_len(Ntip(ph2))]), by='LABEL')
	subset(tmp, round(H1,d=3)!=round(H2,d=3))		
	ph1	<- read.tree('/Users/Oliver/git/HPTN071sim/tmp150227-E-o111/150129_HPTN071_scE-o111_SIMULATED_TREE/150129_HPTN071_scE-o111_SUBSTTREE.newick')
	ph2	<- read.tree('/Users/Oliver/git/HPTN071sim/tmp150227-E-o111 copy/150129_HPTN071_scE-o111_SIMULATED_TREE/150129_HPTN071_scE-o111_SUBSTTREE.newick')
	setequal(ph1$tip.label,ph2$tip.label)
	tmp	<- merge( 	data.table(LABEL=ph1$tip.label, H1=node.depth.edgelength(ph1)[seq_len(Ntip(ph1))]), 
			data.table(LABEL=ph2$tip.label, H2=node.depth.edgelength(ph2)[seq_len(Ntip(ph2))]), by='LABEL')
	subset(tmp, round(H1,d=3)!=round(H2,d=3))		
	
	load('/Users/Oliver/git/HPTN071sim/tmp150227-E-o111 copy/SeqGen/150129_HPTN071_scE-o111_seqgen.R')
	df.seqgen2	<- copy(df.seqgen)
	load('/Users/Oliver/git/HPTN071sim/tmp150227-E-o111/SeqGen/150129_HPTN071_scE-o111_seqgen.R')
	setnames(df.seqgen2, 'NEWICK', 'NEWICK2')
	tmp	<- merge(subset(df.seqgen, GENE=='GAG' & CODON_POS=='CP1', select=c(IDCLU, NEWICK)), subset(df.seqgen2, GENE=='GAG' & CODON_POS=='CP1', select=c(IDCLU, NEWICK2)), by='IDCLU')
	for(clu in tmp[,IDCLU])
	{
		clu<-1
		ph1	<- subset(tmp, IDCLU==clu)[, read.tree(text=NEWICK)]
		ph2	<- subset(tmp, IDCLU==clu)[, read.tree(text=NEWICK2)]
		stopifnot( setequal(ph1$tip.label,ph2$tip.label) )
		xx	<- merge( 	data.table(LABEL=ph1$tip.label, H1=node.depth.edgelength(ph1)[seq_len(Ntip(ph1))]), 
						data.table(LABEL=ph2$tip.label, H2=node.depth.edgelength(ph2)[seq_len(Ntip(ph2))]), by='LABEL')
		subset(xx, round(H1,d=4)!=round(H2,d=4))		
	}	
	
	x1	<- data.table(FILE=list.files('/Users/Oliver/git/HPTN071sim/tmp150227-E-o111 copy 2/VirusTreeSimulator', full.names=1))
	x2	<- data.table(FILE=list.files('/Users/Oliver/git/HPTN071sim/tmp150227-E-o111 copy/VirusTreeSimulator', full.names=1))
	x1[, B:= basename(FILE)]
	x2[, B:= basename(FILE)]
	x1	<- merge(x1,x2,by='B')
	for(i in seq_len(nrow(x1)))
	{
		print(i)
		ph1	<- hivc.beast2out.read.nexus.and.stats(x1[i, FILE.x], method.node.stat='any.node')$tree
		ph2	<- hivc.beast2out.read.nexus.and.stats(x1[i, FILE.y], method.node.stat='any.node')$tree
		ph1	<- seq.collapse.singles(ph1)
		ph2	<- seq.collapse.singles(ph2)
		stopifnot( setequal(ph1$tip.label,ph2$tip.label) )
		if(Nnode(ph1))
		{
			xx	<- merge( 	data.table(LABEL=ph1$tip.label, H1=node.depth.edgelength(ph1)[seq_len(Ntip(ph1))]), 
					data.table(LABEL=ph2$tip.label, H2=node.depth.edgelength(ph2)[seq_len(Ntip(ph2))]), by='LABEL')
			stopifnot( nrow(subset(xx, round(H1,d=4)!=round(H2,d=4)))==0 )	
		}
		if(!Nnode(ph1))
			stopifnot( ph1$root.edge==ph2$root.edge )
		
	}
}
##--------------------------------------------------------------------------------------------------------
##	check simulated sequences: create ExaML tree and estimate R2
##	olli 03.02.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.PropAcute<- function()
{	
	#check time to coalescence in true dated tree
	label.sep		<- '|'
	label.idx.date	<- 4
	label.idx.idpop	<- 1
	indir			<- '/Users/Oliver/duke/2014_Gates/methods_comparison_pipeline/150205'
	dfi				<- data.table(FILE=list.files(indir, '.*DATEDTREE.newick$', full.names=FALSE))
	dfi[, SC:= sapply(strsplit(FILE, '_'),'[[',3)]
	dfi[, CONFIG:= sapply(strsplit(SC, '-'),'[[',2)]
	set(dfi, NULL, 'SC', dfi[, sapply(strsplit(SC, '-'),'[[',1)])	
	scl		<- data.table(SC=c('scA','scB','scC','scD','scE','scF'),SCL=c('Intervention=fast, Acute=low','Intervention=fast, Acute=high','Intervention=slow, Acute=low','Intervention=slow, Acute=high','Intervention=none, Acute=low','Intervention=none, Acute=high'))
	coal	<- dfi[, {
						file			<- paste(indir, FILE, sep='/' )
						cat(paste('\nprocess',file))
						phd				<- read.tree(file)						
						load( gsub('DATEDTREE.newick','SIMULATED_INTERNAL.R',file) )
						singletons.n	<- length(which(sapply(phd, Ntip)==1))
						phd				<- phd[ which(sapply(phd, Ntip)>1) ]
						#	sampling times of tips and node.depth of tips for every tree
						tmp				<- lapply(seq_along(phd), function(k)
								{
									ph	<- phd[[k]]
									tmp	<- data.table(	LABEL=ph$tip.label, CLU_IDX=k, 
														IDPOP= as.integer(substring(sapply(strsplit(ph$tip.label,label.sep,fixed=T),'[[',label.idx.idpop),7)),
														TIME_SEQ=as.numeric(sapply(strsplit(ph$tip.label,label.sep,fixed=T),'[[',label.idx.date)))
									tmp[, BRL:=ph$edge.length[ sort(ph$edge[,2],index.return=T)$ix ][ seq_len(Ntip(ph)) ] ]
									tmp
								})
						coal			<- do.call('rbind',tmp)
						coal			<- merge(coal, subset(df.inds, !is.na(DIAG_T), select=c(IDPOP, DIAG_T)), by='IDPOP')
						coal
					}, by=c('SC','CONFIG')]
	coal	<- merge(coal, scl, by='SC')
	coal[, PERIOD_SEQ:= coal[, cut(TIME_SEQ, breaks=c(1980,2004,2014,Inf), labels=c('<=2004','2005-2014','>=2015'))]]
	coal[, PERIOD_DIAG:= coal[, cut(DIAG_T, breaks=c(1980,2004,2014,Inf), labels=c('<=2004','2005-2014','>=2015'))]]
	#ggplot(coalb, aes(x=BRL)) + geom_histogram(binwidth=.5) + facet_grid(SC~CONFIG) + theme_bw()
	ggplot(coal, aes(x=BRL, colour=SCL)) + stat_ecdf() + facet_grid(CONFIG~PERIOD_DIAG) + theme_bw() + labs(x='time to coalescence of tips',y='empirical CDF') +
			scale_color_brewer(name='scenarios', palette='Paired')	
	file	<- paste(indir, '150205_TEST_Time2CoalescenceOfTips.pdf')
	ggsave(file=file, w=12, h=15)
	
	coal[, PERIOD_SEQ:= coal[, cut(TIME_SEQ, breaks=c(1980,2004,2012,2016, Inf), labels=c('<=2004','2005-2012','2013-2016','>=2017'))]]
	coal[, PERIOD_DIAG:= coal[, cut(DIAG_T, breaks=c(1980,2004,2012,2016, Inf), labels=c('<=2004','2005-2012','2013-2016','>=2017'))]]
	#ggplot(coalb, aes(x=BRL)) + geom_histogram(binwidth=.5) + facet_grid(SC~CONFIG) + theme_bw()
	ggplot(coal, aes(x=BRL, colour=SCL)) + stat_ecdf() + facet_grid(CONFIG~PERIOD_DIAG) + theme_bw() + labs(x='time to coalescence of tips',y='empirical CDF') +
			scale_color_brewer(name='scenarios', palette='Paired')	
	file	<- paste(indir, '150205_TEST_Time2CoalescenceOfTips2.pdf')
	ggsave(file=file, w=12, h=15)
	
	coal[, PERIOD_SEQ:= coal[, cut(TIME_SEQ, breaks=c(1980,2004,2009, Inf), labels=c('<=2004','2005-2009','>=2010'))]]
	coal[, PERIOD_DIAG:= coal[, cut(DIAG_T, breaks=c(1980,2004,2009, Inf), labels=c('<=2004','2005-2009','>=2010'))]]
	#ggplot(coalb, aes(x=BRL)) + geom_histogram(binwidth=.5) + facet_grid(SC~CONFIG) + theme_bw()
	ggplot(coal, aes(x=BRL, colour=SCL)) + stat_ecdf() + facet_grid(CONFIG~PERIOD_DIAG) + theme_bw() + labs(x='time to coalescence of tips',y='empirical CDF') +
			scale_color_brewer(name='scenarios', palette='Paired')	
	file	<- paste(indir, '150204_TEST_Time2CoalescenceOfTips3.pdf')
	ggsave(file=file, w=12, h=15)
}
##--------------------------------------------------------------------------------------------------------
##	check simulated sequences: create ExaML tree and estimate R2
##	olli 02.02.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.SSApg.CLUSTERBEAST.extendedskyline<- function()
{
	require(phytools)
	require(hivclust)
	require(XML)
	#	load outgroup sequences
	file			<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfg_HXB2outgroup.R')
	cat(paste('\nLoading outgroup seq from file', file))
	load(file)		#expect "outgroup.seq.gag" "outgroup.seq.pol" "outgroup.seq.env"	
	
	tree.id.labelsep		<- '|'
	tree.id.label.idx.ctime	<- 4 	
	select		<- 'mseq500'
	#select		<- 'mseq100'
	indir		<- '/Users/Oliver/duke/2014_Gates/methods_comparison_pipeline/150203'
	#indir		<- '/Users/Oliver/git/HPTN071sim/tmp140914/140716_RUN001_INTERNAL'  
	outdir		<- indir
	infiles		<- list.files(indir, '.*INTERNAL.R$', full.names=FALSE)
	#stopifnot(length(infiles)==1)
	#	read BEAST template files	
	infile.beast.pol	<- system.file(package="rPANGEAHIVsim", "misc",'BEAST_template_vTESTpol.xml')
	bxml.template.polut	<- xmlTreeParse(infile.beast.pol, useInternalNodes=TRUE, addFinalizer = TRUE)
	infile.beast.pol	<- system.file(package="rPANGEAHIVsim", "misc",'BEAST_template_vTESTpol_fixedtree.xml')
	#infile.beast.pol	<- '~/git/HPTN071sim/source/rPANGEAHIVsim/inst/misc/BEAST_template_vTESTpol_fixedtree.xml'
	bxml.template.polft	<- xmlTreeParse(infile.beast.pol, useInternalNodes=TRUE, addFinalizer = TRUE)
	infile.beast.gag	<- system.file(package="rPANGEAHIVsim", "misc",'BEAST_template_vTESTgag.xml')
	bxml.template.gag	<- xmlTreeParse(infile.beast.gag, useInternalNodes=TRUE, addFinalizer = TRUE)
	infile.beast.env	<- system.file(package="rPANGEAHIVsim", "misc",'BEAST_template_vTESTenv.xml')		
	bxml.template.env	<- xmlTreeParse(infile.beast.env, useInternalNodes=TRUE, addFinalizer = TRUE)	
	#
	#	run  
	#	
	for(i in seq_along(infiles))
	{
		infile			<- infiles[i]
		#	load simulated data
		file			<- paste(indir, '/', infile, sep='')
		cat(paste('\nLoading file', file))
		load(file)		#expect "df.epi"    "df.trms"   "df.inds"   "df.sample" "df.seq"
		set( df.seq, NULL, 'IDCLU', df.seq[, as.integer(IDCLU)] )
		#			
		if(grepl('nseq',select))
		{
			thresh.NSEQ		<- as.numeric(substring(select, 5)) 
			thresh.brl		<- c(seq(0.001, 0.05, 0.001), seq(0.06, 0.5, 0.1))
			thresh.nseq		<- sapply(thresh.brl, function(x)
					{
						clustering		<- hivc.clu.clusterbythresh(seq.ph, thresh.brl=x, dist.brl=tmp, retval="all")
						which(clustering$size.tips>10)
						length(which(!is.na(clustering$clu.mem[ seq_len(Ntip(seq.ph))] )))					
					})
			thresh.brl		<- thresh.brl[ which(thresh.nseq>=thresh.NSEQ)[1] ]
			clustering		<- hivc.clu.clusterbythresh(seq.ph, thresh.brl=thresh.brl, dist.brl=tmp, retval="all")		
			cat(paste('\nFound clusters, n=', length(clustering$clu.idx)))
			seq.select		<- subset( data.table( PH_NODE_ID=seq_len(Ntip(seq.ph)), CLU_ID=clustering$clu.mem[ seq_len(Ntip(seq.ph))] ), !is.na(CLU_ID) )
			seq.select[, LABEL:= seq.select[, seq.ph$tip.label[PH_NODE_ID]] ]
			seq.select		<- merge(df.seq, seq.select, by='LABEL')
		}
		if(grepl('mseq',select))
		{
			thresh.NSEQ		<- as.numeric(substring(select, 5))
			seq.select		<- subset(df.inds, !is.na(IDCLU))
			seq.select		<- subset(merge(seq.select, seq.select[, list(CLU_N=-length(which(!is.na(TIME_SEQ)))), by='IDCLU'], by='IDCLU'), CLU_N<0 & !is.na(TIME_SEQ))
			setkey(seq.select, CLU_N, IDCLU)
			tmp				<- unique(seq.select)
			tmp[, CLU_CN:= tmp[,cumsum(-CLU_N)]]
			tmp				<- tmp[seq_len( tmp[, which(CLU_CN>=thresh.NSEQ)[1]] ), ] 
			seq.select		<- merge( seq.select, subset(tmp, select=IDCLU), by='IDCLU' )
			seq.select		<- merge(df.seq, seq.select, by=c('IDCLU','IDPOP'))			
			cat(paste('\nFound clusters, n=', seq.select[, length(unique(IDCLU))])) 
			cat(paste('\nFound sequences, n=', seq.select[, length(unique(IDPOP))]))					
		}
		#
		#	read NEWICK trees for each cluster phylogeny, if there
		#
		tmp		<- list.files(indir, '_DATEDTREE.newick$', full.names=FALSE)
		tmp		<- tmp[ grepl(substr(infile, 1, regexpr('_SIMULATED',infile)), tmp) ]
		if(length(tmp))
		{
			# select
			phd					<- read.tree(paste(indir, tmp, sep='/'))
			
			tmp2				<- data.table(IDPOP= sapply(phd, function(x) x$tip.label[1]), IDX=seq_along(phd))
			set( tmp2, NULL, 'IDPOP', tmp2[, as.integer(substring(sapply(strsplit(IDPOP, tree.id.labelsep, fixed=TRUE),'[[',1),7)) ] )
			tmp2				<- merge(subset(seq.select, select=c(IDPOP, IDCLU)), tmp2, by='IDPOP')			
			phd					<- lapply(tmp2[,IDX], function(i) phd[[i]] )
			names(phd)			<- tmp2[, IDCLU]
			# plot
			phd.plot			<- eval(parse(text=paste('phd[[',seq_along(phd),']]', sep='',collapse='+')))			
			#phd.plot			<- drop.tip(phd.plot, which(grepl('NOEXIST', phd.plot$tip.label)), root.edge=1)
			phd.plot			<- ladderize(phd.plot)
			tmp					<- paste(indir, '/', gsub('DATEDTREE','BEASTDATEDTREE',tmp), sep='')						
			pdf(file=gsub('newick','pdf',tmp), w=10, h=Ntip(phd.plot)*0.1)
			plot(phd.plot, show.tip=TRUE, cex=0.5)
			dev.off()
			
			bxml.template.pol	<- bxml.template.polft						
		}
		else
		{
			phd					<- NULL
			bxml.template.pol	<- bxml.template.polut
		}
		#
		#	create BEAST XML
		#		
		#
		#	POL
		#
		if(1)
		{
			cat(paste('\ncreate POL BEAST XML file for seqs=',paste( seq.select[,LABEL], collapse=' ')))
			pool.infile		<- paste(  substr(infile,1,nchar(infile)-21),'_TEST_pol', sep='' )
			#	write XML file with new sequences
			bxml			<- newXMLDoc(addFinalizer=T)
			bxml.beast		<- newXMLNode("beast", doc=bxml, addFinalizer=T)
			tmp				<- newXMLCommentNode(text=paste("Generated by HIVCLUST from template",infile.beast.pol), parent=bxml.beast, doc=bxml, addFinalizer=T)
			#	add new taxa 
			tmp 			<- subset(seq.select, select=LABEL)
			tmp[, BEASTlabel:=LABEL]
			setnames(tmp, 'LABEL', 'FASTASampleCode')
			bxml			<- hivc.beast.add.taxa(bxml, tmp, beast.label.datepos= 4, beast.label.sep= '|', beast.date.direction= "forwards", beast.date.units= "years", verbose=1)
			#	add alignment for every transmission cluster
			for(clu in seq.select[, unique(IDCLU)])
			{
				tmp				<- tolower(do.call('rbind',strsplit(subset(seq.select, IDCLU==clu)[, POL],'')))
				rownames(tmp)	<- subset(seq.select, IDCLU==clu)[, LABEL]
				tmp				<- as.DNAbin(tmp)
				bxml			<- hivc.beast.add.alignment(bxml, tmp, beast.alignment.id=paste("POL.alignment_CLU",clu,sep=''), beast.alignment.dataType= "nucleotide", verbose=1)				
			}
			#	add starting tree for every transmission cluster			
			for(clu in names(phd))
			{
				bxml			<- hivc.beast.add.startingtree(bxml, phd[[clu]], beast.rootHeight= NA, beast.usingDates="true", beast.newickid=paste("startingTree_CLU",clu,sep=''), beast.brlunits="years", verbose=1)	
			}
			#	add CODON patterns for every alignment
			for(clu in seq.select[, unique(IDCLU)])
				for(k in 1:3)
				{					
					bxml		<- hivc.beast.add.patterns(bxml, paste("POL.n100.rlx.gmrf.CP",k,".patterns_CLU",clu,sep=''), paste("POL.alignment_CLU",clu,sep=''), k, beast.patterns.every=3, beast.patterns.strip='false', verbose=1)
				}
			#	add tree model for every alignment
			for(clu in seq.select[, unique(IDCLU)])
			{
				bxml			<- hivc.beast.add.treemodel(bxml, paste("treeModel_CLU",clu,sep=''), paste("treeModel.rootHeight_CLU",clu,sep=''), paste("treeModel.internalNodeHeights_CLU",clu,sep=''), paste("treeModel.allInternalNodeHeights_CLU",clu,sep=''), newick.id=paste("startingTree_CLU",clu,sep=''), internalNodes='true', rootNode='true', verbose=1)
			}
			#	add variableDemographic
			tmp					<- paste('treeModel_CLU',seq.select[, unique(IDCLU)], sep='')
			bxml				<- hivc.beast.add.variableDemographic(bxml, 'demographic', tmp, 'coalescent', type='stepwise', useMidpoints='true', 
					popSize.id='demographic.popSize', popSize.value='500.0', 
					demographic.indicators.id='demographic.indicators', demographic.indicators.value='0.0',
					demographic.popMeanDist.id='demographic.populationMeanDist', demographic.popMean.id='demographic.populationMean', demographic.popMean.value='500.0', 
					sumStatistic.id='demographic.populationSizeChanges', sumStatistic.elementwise='true', 
					verbose=1)
#	population mean value: fix?										
			#	add uncorrelated relaxed clock for every transmission cluster
			tmp					<- seq.select[, unique(IDCLU)]
			clu					<- tmp[1]
			bxml				<- hivc.beast.add.discretizedBranchRates(bxml, paste("branchRates_CLU",clu,sep=''), paste("treeModel_CLU",clu,sep=''), paste("branchRates.categories_CLU",clu,sep=''),
					mean.id="ucld.mean", mean.value="0.0025", mean.lower="0.0", sd.id="ucld.stdev", sd.value="0.3333333333333333", sd.lower="0.0", meanInRealSpace="true", 
					rateCategories.dimension= subset(seq.select, IDCLU==clu)[, 2*(length(unique(IDPOP))-1)], verbose=1)
			for(clu in tmp[-1])
			{
				bxml			<- hivc.beast.add.discretizedBranchRates(bxml, paste("branchRates_CLU",clu,sep=''), paste("treeModel_CLU",clu,sep=''), paste("branchRates.categories_CLU",clu,sep=''),
						mean.idref="ucld.mean", sd.idref="ucld.stdev", meanInRealSpace="true", 
						rateCategories.dimension= subset(seq.select, IDCLU==clu)[, 2*(length(unique(IDPOP))-1)], verbose=1)
			}		
			#	add rate statistic for first cluster
			bxml				<- hivc.beast.add.rateStatistics(bxml, 'branchRates', paste('treeModel_CLU',seq.select[1, IDCLU], sep=''), paste('branchRates_CLU',seq.select[1, IDCLU], sep=''))
			#	add GTR CP1 2 3  models
			for(k in 1:3)
			{
				clu			<- seq.select[, unique(IDCLU)][1]
				dummy		<- hivc.beast.add.gtrModel(bxml, paste("POL.n100.rlx.gmrf.CP",k,".gtr_CLU", clu,sep=""), paste("POL.n100.rlx.gmrf.CP",k,".patterns_CLU",clu,sep=''), paste("POL.n100.rlx.gmrf.CP",k,".frequencies_CLU",clu,sep=''),
						rateAC.id=paste("POL.n100.rlx.gmrf.CP",k,".ac",sep=""), rateAC.value="0.3", rateAC.lower="0.0",
						rateAG.id=paste("POL.n100.rlx.gmrf.CP",k,".ag",sep=""), rateAG.value="1.0", rateAG.lower="0.0",
						rateAT.id=paste("POL.n100.rlx.gmrf.CP",k,".at",sep=""), rateAT.value="0.2", rateAT.lower="0.0",
						rateCG.id=paste("POL.n100.rlx.gmrf.CP",k,".cg",sep=""), rateCG.value="0.2", rateCG.lower="0.0",
						rateGT.id=paste("POL.n100.rlx.gmrf.CP",k,".gt",sep=""), rateGT.value="0.1", rateGT.lower="0.0",
						frequencies.dimension='4', frequencyModel.dataType='nucleotide')
				for(clu in seq.select[, unique(IDCLU)][-1])
					dummy	<- hivc.beast.add.gtrModel(bxml, paste("POL.n100.rlx.gmrf.CP",k,".gtr_CLU", clu,sep=""), paste("POL.n100.rlx.gmrf.CP",k,".patterns_CLU",clu,sep=''), paste("POL.n100.rlx.gmrf.CP",k,".frequencies_CLU",clu,sep=''),
							rateAC.idref=paste("POL.n100.rlx.gmrf.CP",k,".ac",sep=""), 
							rateAG.idref=paste("POL.n100.rlx.gmrf.CP",k,".ag",sep=""), 
							rateAT.idref=paste("POL.n100.rlx.gmrf.CP",k,".at",sep=""), 
							rateCG.idref=paste("POL.n100.rlx.gmrf.CP",k,".cg",sep=""), 
							rateGT.idref=paste("POL.n100.rlx.gmrf.CP",k,".gt",sep=""), 
							frequencies.dimension='4', frequencyModel.dataType='nucleotide')				
			}	
			#	add CP1 2 3  site models 
			for(k in 1:3)
			{
				clu			<- seq.select[, unique(IDCLU)][1]
				dummy		<- hivc.beast.add.siteModel(bxml, paste("POL.n100.rlx.gmrf.CP",k,".siteModel_CLU", clu,sep=""), paste("POL.n100.rlx.gmrf.CP",k,".gtr_CLU", clu,sep=""), 
						relativeRate.id=paste("POL.n100.rlx.gmrf.CP",k,".mu",sep=""), relativeRate.value="1.0", relativeRate.lower="0.0",
						gamma.id=paste("POL.n100.rlx.gmrf.CP",k,".alpha",sep=""), gamma.value="0.5", gamma.lower="0.0",
						gammaCategories="4", verbose=1)
				for(clu in seq.select[, unique(IDCLU)][-1])
					dummy	<- hivc.beast.add.siteModel(bxml, paste("POL.n100.rlx.gmrf.CP",k,".siteModel_CLU", clu,sep=""), paste("POL.n100.rlx.gmrf.CP",k,".gtr_CLU", clu,sep=""), 
							relativeRate.idref=paste("POL.n100.rlx.gmrf.CP",k,".mu",sep=""), gamma.idref=paste("POL.n100.rlx.gmrf.CP",k,".alpha",sep=""), gammaCategories="4", verbose=1)				
			}
			#	copy compoundParameter
			dummy				<- addChildren( bxml.beast, xmlClone(getNodeSet(bxml.template.pol, paste("//*[@id='ALL.n100.rlx.gmrf.allMus']",sep=''))[[1]], addFinalizer=T, doc=bxml) )
			#	add treeLikelihood
			for(clu in seq.select[, unique(IDCLU)])
				for(k in 1:3)
				{
					dummy		<- hivc.beast.add.treeLikelihood(bxml, paste("POL.n100.rlx.gmrf.CP",k,".treeLikelihood_CLU",clu,sep=""), 
							paste("POL.n100.rlx.gmrf.CP",k,".patterns_CLU",clu,sep=''), 
							paste("treeModel_CLU",clu,sep=''), 
							paste('POL.n100.rlx.gmrf.CP',k,'.siteModel_CLU',clu,sep=''), 
							paste("branchRates_CLU",clu,sep=''), useAmbiguities='false', verbose=1)
				}
			#	add compound likelihood
			tmp					<- as.vector(sapply(seq.select[, unique(IDCLU)], function(clu) sapply(1:3,function(k){ paste("POL.n100.rlx.gmrf.CP",k,".treeLikelihood_CLU",clu,sep="") }) ))
			bxml				<- hivc.beast.add.likelihood(bxml, 'likelihood', tmp)
			#	copy operators
			dummy				<- addChildren( bxml.beast, xmlClone(getNodeSet(bxml.template.pol, paste("//operators",sep=''))[[1]], addFinalizer=T, doc=bxml) )
			#	add branch rate category operators
			for(clu in seq.select[, unique(IDCLU)])
			{
				bxml			<- hivc.beast.add.uniformIntegerOperator(bxml, paste("branchRates.categories_CLU",clu,sep=''), "10")
				bxml			<- hivc.beast.add.swapOperator(bxml, paste("branchRates.categories_CLU",clu,sep=''), "1", "10", "false")
			}
			#	copy mcmc
			dummy				<- addChildren( bxml.beast, xmlClone(getNodeSet(bxml.template.pol, paste("//mcmc",sep=''))[[1]], addFinalizer=T, doc=bxml) )
			#	add tree log for every transmission cluster
			for(clu in seq.select[, unique(IDCLU)])
			{				
				bxml			<- hivc.beast.add.logTree(bxml, paste("logTree_CLU",clu,sep=''), paste("treeModel_CLU",clu,sep=''), paste("branchRates_CLU",clu,sep=''), 'posterior', logEvery='3000', fileName=paste(pool.infile,'_CLU', clu, '_', select, '.trees', sep=''), nexusFormat='true', sortTranslationTable='true', verbose=1)	
			}
			#	change MCMC attributes
			bxml				<- hivc.beast.adjust.mcmc(bxml, beast.mcmc.chainLength=3000000, beast.mcmc.logEvery=3000, verbose=1)
			#	change outfile name of log 
			tmp							<- getNodeSet(bxml, "//*[@id='fileLog']")[[1]]
			xmlAttrs(tmp)["fileName"]	<- paste(pool.infile,'_', select, '.log', sep='')			
			#	add VDAnalysis
			tmp					<- sapply(seq.select[, unique(IDCLU)], function(clu) paste(pool.infile,'_CLU', clu, '_', select, '.trees', sep=''))			
			bxml				<- hivc.beast.add.VDAnalysis(bxml, 'demographic.analysis', paste(pool.infile,'_', select, '.log', sep=''), tmp, paste(pool.infile,'_', select, '_VDA.csv', sep=''), populationModelType='stepwise', populationFirstColumn='demographic.popSize1', indicatorsFirstColumn='demographic.indicators1', burnIn="0.1", useMidpoints="true", verbose=1)
			#	write to file
			file		<- paste(indir,'/',pool.infile,'_',select,".xml", sep='')
			cat(paste("\nwrite xml file to",file))
			saveXML(bxml, file=file)	
		}				
	}
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
	#	load outgroup sequences
	file			<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfg_HXB2outgroup.R')
	cat(paste('\nLoading outgroup seq from file', file))
	load(file)		#expect "outgroup.seq.gag" "outgroup.seq.pol" "outgroup.seq.env"	
	
	tree.id.labelsep		<- '|'
	tree.id.label.idx.ctime	<- 4 	
	select		<- 'nseq800'
	select		<- 'nseq500'
	select		<- 'mseq100'
	indir		<- '/Users/Oliver/duke/2014_Gates/methods_comparison_pipeline/150201b'
	#indir		<- '/Users/Oliver/git/HPTN071sim/tmp140914/140716_RUN001_INTERNAL'  
	outdir		<- indir
	infiles		<- list.files(indir, '.*INTERNAL.R$', full.names=FALSE)
	#stopifnot(length(infiles)==1)
	#	read BEAST template files	
	infile.beast.pol	<- system.file(package="rPANGEAHIVsim", "misc",'BEAST_template_vTESTpol.xml')
	bxml.template.polut	<- xmlTreeParse(infile.beast.pol, useInternalNodes=TRUE, addFinalizer = TRUE)
	infile.beast.pol	<- system.file(package="rPANGEAHIVsim", "misc",'BEAST_template_vTESTpol_fixedtree.xml')
	bxml.template.polft	<- xmlTreeParse(infile.beast.pol, useInternalNodes=TRUE, addFinalizer = TRUE)
	infile.beast.gag	<- system.file(package="rPANGEAHIVsim", "misc",'BEAST_template_vTESTgag.xml')
	bxml.template.gag	<- xmlTreeParse(infile.beast.gag, useInternalNodes=TRUE, addFinalizer = TRUE)
	infile.beast.env	<- system.file(package="rPANGEAHIVsim", "misc",'BEAST_template_vTESTenv.xml')		
	bxml.template.env	<- xmlTreeParse(infile.beast.env, useInternalNodes=TRUE, addFinalizer = TRUE)	
	#
	#	run  
	#	
	for(i in seq_along(infiles))
	{
		infile			<- infiles[i]
		#	load simulated data
		file			<- paste(indir, '/', infile, sep='')
		cat(paste('\nLoading file', file))
		load(file)		#expect "df.epi"    "df.trms"   "df.inds"   "df.sample" "df.seq"
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
		tmp				<- dist.dna( seq, model='raw' )
		seq.ph			<- nj(tmp)		
		tmp				<- which(seq.ph$tip.label=="HXB2")
		seq.ph			<- reroot(seq.ph, tmp, seq.ph$edge.length[which(seq.ph$edge[,2]==tmp)])
		#	
		tmp				<- hivc.clu.brdist.stats(seq.ph, eval.dist.btw="leaf", stat.fun=hivc.clu.min.transmission.cascade)
		if(grepl('nseq',select))
		{
			thresh.NSEQ		<- as.numeric(substring(select, 5)) 
			thresh.brl		<- c(seq(0.001, 0.05, 0.001), seq(0.06, 0.5, 0.1))
			thresh.nseq		<- sapply(thresh.brl, function(x)
					{
						clustering		<- hivc.clu.clusterbythresh(seq.ph, thresh.brl=x, dist.brl=tmp, retval="all")
						which(clustering$size.tips>10)
						length(which(!is.na(clustering$clu.mem[ seq_len(Ntip(seq.ph))] )))					
					})
			thresh.brl		<- thresh.brl[ which(thresh.nseq>=thresh.NSEQ)[1] ]
			clustering		<- hivc.clu.clusterbythresh(seq.ph, thresh.brl=thresh.brl, dist.brl=tmp, retval="all")		
			cat(paste('\nFound clusters, n=', length(clustering$clu.idx)))
			seq.select		<- subset( data.table( PH_NODE_ID=seq_len(Ntip(seq.ph)), CLU_ID=clustering$clu.mem[ seq_len(Ntip(seq.ph))] ), !is.na(CLU_ID) )
			seq.select[, LABEL:= seq.select[, seq.ph$tip.label[PH_NODE_ID]] ]
			seq.select		<- merge(df.seq, seq.select, by='LABEL')
		}
		if(grepl('mseq',select))
		{
			thresh.NSEQ		<- as.numeric(substring(select, 5))
			seq.select		<- subset(df.inds, !is.na(IDCLU))
			seq.select		<- subset(merge(seq.select, seq.select[, list(CLU_N=-length(which(!is.na(TIME_SEQ)))), by='IDCLU'], by='IDCLU'), CLU_N<0 & !is.na(TIME_SEQ))
			setkey(seq.select, CLU_N, IDCLU)
			tmp				<- unique(seq.select)
			tmp[, CLU_CN:= tmp[,cumsum(-CLU_N)]]
			tmp				<- tmp[seq_len( tmp[, which(CLU_CN>=thresh.NSEQ)[1]] ), ] 
			seq.select		<- merge( seq.select, subset(tmp, select=IDCLU), by='IDCLU' )
			cat(paste('\nFound clusters, n=', seq.select[, length(unique(IDCLU))])) 
			cat(paste('\nFound sequences, n=', seq.select[, length(unique(IDPOP))]))							
		}
		if(select=='same')
		{
			thresh.brl		<- seq(0.01, 0.05, 0.001)
			thresh.nclu		<- sapply(thresh.brl, function(x)
					{
						clustering		<- hivc.clu.clusterbythresh(seq.ph, thresh.brl=x, dist.brl=tmp, retval="all")
						max( clustering$size.tips )					
					})
			thresh.brl		<- thresh.brl[ which(thresh.nclu>=100)[1] ]
			clustering		<- hivc.clu.clusterbythresh(seq.ph, thresh.brl=thresh.brl, dist.brl=tmp, retval="all")		
			cat(paste('\nFound clusters, n=', length(clustering$clu.idx)))
			#	Take 1 sequence from each cluster
			tmp				<- which( clustering$clu.mem==which.max( clustering$size.tips ) )
			seq.select		<- data.table( PH_NODE_ID=tmp[ which(tmp<=Ntip(seq.ph)) ], CLU_ID=which.max(clustering$size.tips) )
			seq.select		<- seq.select[1:100,]
			seq.select[, LABEL:=seq.ph$tip.label[PH_NODE_ID] ]			
			seq.select		<- merge(df.seq, seq.select, by='LABEL')	
			cat(paste('\nSelected sequences, n=',nrow(seq.select)))
		}
		if(select=='divergent')
		{
			#	find thresh with ~100 clusters
			thresh.brl		<- c(seq(0.001, 0.05, 0.001), seq(0.06, 0.5, 0.1))
			thresh.nclu		<- sapply(thresh.brl, function(x)
					{
						clustering		<- hivc.clu.clusterbythresh(seq.ph, thresh.brl=x, dist.brl=tmp, retval="all")
						length(clustering$clu.idx)					
					})			
			tmp2			<- which(thresh.nclu>100)
			thresh.brl		<- ifelse(length(tmp2), thresh.brl[tmp2][1], thresh.brl[which.max(thresh.nclu)])
			clustering		<- hivc.clu.clusterbythresh(seq.ph, thresh.brl=thresh.brl, dist.brl=tmp, retval="all")		
			cat(paste('\nFound clusters, n=', length(clustering$clu.idx)))
			#	Take 1 sequence from each cluster
			seq.select		<- data.table( PH_NODE_ID=seq_len(Ntip(seq.ph)), CLU_ID=clustering$clu.mem[ seq_len(Ntip(seq.ph)) ] )
			seq.select		<- subset(seq.select, !is.na(CLU_ID))[, list(LABEL= seq.ph$tip.label[PH_NODE_ID[1]]), by='CLU_ID']
			seq.select		<- merge(df.seq, seq.select, by='LABEL')			
		}
		#
		#	read NEWICK trees for each cluster phylogeny, if there
		#
		tmp		<- list.files(indir, '_DATEDTREE.newick$', full.names=FALSE)
		tmp		<- tmp[ grepl(substr(infile, 1, regexpr('_SIMULATED',infile)), tmp) ]
		if(length(tmp))
		{
			# select
			phd					<- read.tree(paste(indir, tmp, sep='/'))
			tmp2				<- data.table(IDCLU= as.integer(substring(names(phd),7)), IDX=seq_along(phd))
			tmp2				<- merge(unique(seq.select), tmp2, by='IDCLU')[, unique(IDX)]
			phd					<- lapply(tmp2, function(i) phd[[i]] )
			# plot
			phd.plot			<- eval(parse(text=paste('phd[[',seq_along(phd),']]', sep='',collapse='+')))			
			#phd.plot			<- drop.tip(phd.plot, which(grepl('NOEXIST', phd.plot$tip.label)), root.edge=1)
			phd.plot			<- ladderize(phd.plot)
			tmp					<- paste(indir, '/', gsub('DATEDTREE','BEASTDATEDTREE',tmp), sep='')						
			pdf(file=gsub('newick','pdf',tmp), w=10, h=Ntip(phd.plot)*0.1)
			plot(phd.plot, show.tip=TRUE, cex=0.5)
			dev.off()
			
			bxml.template.pol	<- bxml.template.polft						
		}
		else
		{
			phd					<- NULL
			bxml.template.pol	<- bxml.template.polut
		}
		#
		#	create BEAST XML
		#
		if(0)
		{
			#
			#	GAG
			#
			cat(paste('\ncreate GAG BEAST XML file for seqs=',paste( seq.select[,LABEL], collapse=' ')))
			pool.infile		<- paste(  substr(infile,1,nchar(infile)-21),'_TEST_gag', sep='' )
			#	write XML file with new sequences
			bxml			<- newXMLDoc(addFinalizer=T)
			bxml.beast		<- newXMLNode("beast", doc=bxml, addFinalizer=T)
			tmp				<- newXMLCommentNode(text=paste("Generated by HIVCLUST from template",infile.beast.pol), parent=bxml.beast, doc=bxml, addFinalizer=T)
			#	add new set of GAG sequences into GAG alignment
			tmp				<- tolower(do.call('rbind',strsplit(seq.select[, GAG],'')))
			rownames(tmp)	<- seq.select[, LABEL]
			tmp				<- as.DNAbin(tmp)
			bxml			<- hivc.beast.add.seq(bxml, tmp, df=NULL, beast.label.datepos=4, beast.label.sep= '|', beast.date.direction= "forwards", beast.date.units= "years", beast.alignment.id="GAG.alignment", beast.alignment.dataType= "nucleotide", verbose=1)
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
			tmp			<- sapply(strsplit(tmp,'.',fixed=1), function(x)	paste(pool.infile,'_',select, '.', tail(x,1), sep=''))
			dummy		<- sapply(seq_along(bxml.onodes), function(i){		xmlAttrs(bxml.onodes[[i]])["fileName"]<- tmp[i]		})
			#	write to file
			file		<- paste(indir,'/',pool.infile,'_',select,".xml", sep='')
			cat(paste("\nwrite xml file to",file))
			saveXML(bxml, file=file)	
		}		
		#
		#	POL
		#
		if(1)
		{
			cat(paste('\ncreate POL BEAST XML file for seqs=',paste( seq.select[,LABEL], collapse=' ')))
			pool.infile		<- paste(  substr(infile,1,nchar(infile)-21),'_TEST_pol', sep='' )
			#	write XML file with new sequences
			bxml			<- newXMLDoc(addFinalizer=T)
			bxml.beast		<- newXMLNode("beast", doc=bxml, addFinalizer=T)
			tmp				<- newXMLCommentNode(text=paste("Generated by HIVCLUST from template",infile.beast.pol), parent=bxml.beast, doc=bxml, addFinalizer=T)
			#	add new set of POL sequences into POL alignment
			tmp				<- tolower(do.call('rbind',strsplit(seq.select[, POL],'')))
			rownames(tmp)	<- seq.select[, LABEL]
			tmp				<- as.DNAbin(tmp)
			bxml			<- hivc.beast.add.seq(bxml, tmp, df=NULL, beast.label.datepos=4, beast.label.sep= '|', beast.date.direction= "forwards", beast.date.units= "years", beast.alignment.id="POL.alignment", beast.alignment.dataType= "nucleotide", verbose=1)
			
			df.seq[, BEASTlabel:= LABEL]
			setnames(df.seq, 'LABEL', 'FASTASampleCode')			
			bxml			<- hivc.beast.add.startingtree(bxml, phd, df.seq, beast.rootHeight= NA, beast.usingDates="true", beast.newickid= "startingTree", beast.brlunits="years", verbose=1)				
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
			#	change MCMC attributes
			if(!is.null(phd))
			{
				bxml	<- hivc.beast.adjust.mcmc(bxml, beast.mcmc.chainLength=3000000, beast.mcmc.logEvery=3000, verbose=1)				
			}
			#	change outfile name 
			bxml.onodes	<- getNodeSet(bxml, "//*[@fileName]")
			tmp			<- sapply(bxml.onodes, function(x) xmlGetAttr(x,"fileName"))
			tmp			<- gsub("(time).","time",tmp,fixed=1)
			tmp			<- gsub("(subst).","subst",tmp,fixed=1)	
			tmp			<- sapply(strsplit(tmp,'.',fixed=1), function(x)	paste(pool.infile,'_',select, '.', tail(x,1), sep=''))
			dummy		<- sapply(seq_along(bxml.onodes), function(i){		xmlAttrs(bxml.onodes[[i]])["fileName"]<- tmp[i]		})
			#	write to file
			file		<- paste(indir,'/',pool.infile,'_',select,".xml", sep='')
			cat(paste("\nwrite xml file to",file))
			saveXML(bxml, file=file)	
		}		
		#
		#	ENV
		#
		if(0)
		{
			cat(paste('\ncreate ENV BEAST XML file for seqs=',paste( seq.select[,LABEL], collapse=' ')))
			pool.infile		<- paste(  substr(infile,1,nchar(infile)-21),'_TEST_env', sep='' )
			#	write XML file with new sequences
			bxml			<- newXMLDoc(addFinalizer=T)
			bxml.beast		<- newXMLNode("beast", doc=bxml, addFinalizer=T)
			tmp				<- newXMLCommentNode(text=paste("Generated by HIVCLUST from template",infile.beast.pol), parent=bxml.beast, doc=bxml, addFinalizer=T)
			#	add new set of GAG sequences into ENV alignment
			tmp				<- tolower(do.call('rbind',strsplit(seq.select[, ENV],'')))
			rownames(tmp)	<- seq.select[, LABEL]
			tmp				<- as.DNAbin(tmp)
			bxml			<- hivc.beast.add.seq(bxml, tmp, df=NULL, beast.label.datepos=4, beast.label.sep= '|', beast.date.direction= "forwards", beast.date.units= "years", beast.alignment.id="ENV.alignment", beast.alignment.dataType= "nucleotide", verbose=1)
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
			tmp			<- sapply(strsplit(tmp,'.',fixed=1), function(x)	paste(pool.infile,'_',select, '.', tail(x,1), sep=''))
			dummy		<- sapply(seq_along(bxml.onodes), function(i){		xmlAttrs(bxml.onodes[[i]])["fileName"]<- tmp[i]		})
			#	write to file
			file		<- paste(indir,'/',pool.infile,'_',select,".xml", sep='')
			cat(paste("\nwrite xml file to",file))
			saveXML(bxml, file=file)	
		}
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
	indir		<- '/Users/Oliver/duke/2014_Gates/methods_comparison_pipeline/150227'  
	outdir		<- '/Users/Oliver/duke/2014_Gates/methods_comparison_pipeline/150227'
	infiles		<- list.files(indir, '.*INTERNAL.R$', full.names=FALSE)
	#stopifnot(length(infiles)==1)
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
		#	run ExaML on gag
		#
		seq				<- df.seq.gag
		seq				<- rbind(seq, outgroup.seq.gag[, seq_len(ncol(seq))])
		infile.seq.sig	<- "Sun_Sep_14_12:59:06_2013"
		infile.seq		<- paste(substr(infile,1,nchar(infile)-20),'INFO_simu_gagseq',sep='')
		file			<- paste( outdir, '/', infile.seq,'_',gsub('/',':',infile.seq.sig),'.R', sep='' )
		save(seq, file=file)
		#	run ExaML
		cmd				<- hivc.cmd.examl.bootstrap.on.one.machine(indir, infile.seq, infile.seq.sig, infile.seq.sig, bs.from=0, bs.to=0, verbose=1)
		cmd				<- hivc.cmd.hpcwrapper(cmd, hpc.walltime=21, hpc.q= NA, hpc.mem="450mb", hpc.nproc=1)
		cmd.hpccaller(outdir, paste("exa",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.'), cmd)
		Sys.sleep(1)	
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
		#	run ExaML on env
		#
		seq				<- df.seq.env
		seq				<- rbind(seq, outgroup.seq.env[, seq_len(ncol(seq))])
		infile.seq.sig	<- "Sun_Sep_14_12:59:06_2013"
		infile.seq		<- paste(substr(infile,1,nchar(infile)-20),'INFO_simu_envseq',sep='')
		file			<- paste( outdir, '/', infile.seq,'_',gsub('/',':',infile.seq.sig),'.R', sep='' )
		save(seq, file=file)
		#	run ExaML
		cmd				<- hivc.cmd.examl.bootstrap.on.one.machine(indir, infile.seq, infile.seq.sig, infile.seq.sig, bs.from=0, bs.to=0, verbose=1)
		cmd				<- hivc.cmd.hpcwrapper(cmd, hpc.walltime=21, hpc.q= NA, hpc.mem="450mb", hpc.nproc=1)
		cmd.hpccaller(outdir, paste("exa",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.'), cmd)
		Sys.sleep(1)
		if(1)
		{
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
	}	
	#
	#	evaluate R2 for pol
	#
	gene			<- 'pol'
	infiles			<- list.files(indir, paste('^ExaML_result.*',gene,'seq.*finaltree.000$',sep=''), full.names=FALSE)
	for(i in seq_along(infiles))
	{
		#	 read files
		infile			<- infiles[i]
		file			<- paste(indir,'/',infile,sep='')
		ph				<- read.tree(file)
		tmp				<- regmatches(infile,regexpr('.*_INFO',infile))
		file			<- paste(indir,'/', substr(tmp, 14, nchar(tmp)-4),"SIMULATED_INTERNAL.R",sep='')
		#tmp				<- substring(infile, 14)		
		#tmp				<- list.files(indir, paste(substr(tmp, 1, regexpr('_INFO',tmp) ),'.*INTERNAL.R$', sep=''), full.names=FALSE)
		#cat(paste('\nLoad file=',tmp))		
		load(file)
		#	re-root at HXB2 (outgroup, because the simulated sequences are HIV-1C)
		tmp				<- which(ph$tip.label=="HXB2")
		ph				<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])
		ph				<- ladderize(ph)		
		file			<- paste( outdir, '/', substr(infile,1,nchar(infile)-20),'INFO_simu_ExaML',gene,'.pdf', sep='' )	
		pdf(file=file, w=10, h=150)
		plot(ph, show.tip=TRUE, cex=0.5)
		add.scale.bar()
		dev.off()	
		#	get root to tip divergence
		ph				<- drop.tip(ph,'HXB2')
		file			<- paste( outdir, '/', substr(infile,1,nchar(infile)-20),'INFO_simu_ExaML',gene,'.newick', sep='' )
		write.tree(file=file, ph)		
		tmp				<- node.depth.edgelength(ph)
		ph.info			<- data.table(LABEL=ph$tip.label, ROOT2TIP=tmp[seq_len(Ntip(ph))] )
		set(ph.info, NULL, 'IDPOP', ph.info[, sapply(strsplit(LABEL, tree.id.labelsep, fixed=TRUE),'[[',1) ])
		set(ph.info, NULL, 'IDPOP', ph.info[, as.integer(substr(IDPOP, 7, nchar(IDPOP)))])
		ph.info			<- merge(ph.info, subset(df.inds, select=c(IDPOP, IDCLU)), by='IDPOP')
		set(ph.info, NULL, 'IDCLU', ph.info[, factor(IDCLU)])
		set(ph.info, NULL, 'CALENDAR_TIME', ph.info[, as.numeric(sapply(strsplit(LABEL, tree.id.labelsep, fixed=TRUE),'[[',tree.id.label.idx.ctime))] )
		tmp				<- lm(ROOT2TIP~CALENDAR_TIME, data=ph.info)		 
		set( ph.info, NULL, 'ROOT2TIP_LM', predict(tmp, type='response') ) 	
		tmp2			<- c( 	R2=round(summary(tmp)$r.squared,d=3), 
								SLOPE= as.numeric(round(coef(tmp)['CALENDAR_TIME'],d=4)), 
								TMRCA=as.numeric(round( -coef(tmp)['(Intercept)']/coef(tmp)['CALENDAR_TIME'], d=1 )) )
		ggplot(ph.info, aes(x=CALENDAR_TIME, y=ROOT2TIP, colour=IDCLU)) + geom_point(alpha=0.75) + geom_line(alpha=0.1, aes(group=IDCLU)) + geom_line(aes(y=ROOT2TIP_LM)) +
				#scale_x_continuous(breaks=seq(1980,2020,2)) +						
				scale_colour_discrete(guide=FALSE) +
				labs(x='Sequence sampling date', y=paste('root-to-tip divergence\n(HIV-1',gene,'sequences)')) +
				annotate("text", x=ph.info[, min(CALENDAR_TIME)], y=ph.info[, 0.9*max(ROOT2TIP)], label=paste("R2=", tmp2['R2'],'\nSlope=',tmp2['SLOPE'],'\nTMRCA=',tmp2['TMRCA'], sep=''), hjust = 0, size = 4) +
				theme(legend.position=c(0,1), legend.justification=c(0,1))		
		file			<- paste( outdir, '/', substr(infile,1,nchar(infile)-20),'INFO_simu_ExaML',gene,'R2.pdf', sep='' )
		ggsave(file=file, w=10, h=6)	
		#
		if(0)
		{
			nclu			<- subset(df.inds, !is.na(IDCLU))[, length(unique(IDCLU))]
			cat(paste('\nNumber of transmission clusters, n=', nclu))		
			tmp				<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=hivc.clu.min.transmission.cascade)
			#	almost random choice: 0.05 
			clustering		<- hivc.clu.clusterbythresh(ph, thresh.brl=0.05, dist.brl=tmp, retval="all")
			tmp				<- data.table( LABEL=ph$tip.label, EXACLUID= clustering$clu.mem[ seq_len(Ntip(ph))] )
			tmp[, IDPOP:= tmp[, as.integer(substring(sapply(strsplit(LABEL,'|',fixed=TRUE),'[[',1),7))]]
			ph.inds			<- merge( df.inds, tmp, by='IDPOP' )		
			tmp				<- subset(ph.inds, !is.na(EXACLUID))[, {
						tmp				<- extract.clade(ph, hivc.clu.mrca(ph, LABEL)$mrca, root.edge=1)
						tmp$root.edge	<- 0					
						list( LABEL=tmp$tip.label, ROOT2TIP=node.depth.edgelength(tmp)[seq_len(Ntip(tmp))] )					
					}, by='EXACLUID']
			ph.inds			<- merge(ph.inds, tmp, by=c('EXACLUID','LABEL'), all.x=TRUE)
			ph.inds			<- merge(ph.inds, subset(ph.inds, !is.na(EXACLUID))[, list(EXACLUID_N=length(IDPOP)), by='EXACLUID'], by='EXACLUID', all.x=TRUE)
			clx.info		<- subset(ph.inds, EXACLUID_N==max(EXACLUID_N, na.rm=1))
			ggplot(clx.info, aes(x=TIME_SEQ, y=ROOT2TIP, colour=IDCLU)) + geom_point(alpha=0.75) 	
		}		
	}	
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
		tmp				<- dist.dna(anc.seq.draw, model='raw')
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
		tmp				<- dist.dna(anc.seq.draw, model='raw')
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
	#DATA				<<- '/Users/Oliver/duke/2014_Gates'
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
	#DATA				<<- '/Users/Oliver/duke/2014_Gates'
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
	#DATA			<<- '/Users/Oliver/duke/2014_Gates'
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
	#DATA		<<- "/work/or105/Gates_2014"
	#DATA		<<- '/Users/Oliver/duke/2014_Gates'
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
	#DATA			<<- "/work/or105/Gates_2014"
	#DATA			<<- '/Users/Oliver/duke/2014_Gates'
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


