PR.PACKAGE					<- "rPANGEAHIVsim"
#PR.STARTME					<- system.file(package=PR.PACKAGE, "misc", "rPANGEAHIV.startme.R")
PR.STARTME					<- '/Users/Oliver/git/HPTN071sim/source/rPANGEAHIVsim/misc/rPANGEAHIV.startme.R'
#PR.STARTME					<- '/work/or105/libs/HPTN071sim/source/rPANGEAHIVsim/misc/rPANGEAHIV.startme.R'
PR.HPTN071.INPUT.PARSER4	<- paste('Rscript',system.file(package=PR.PACKAGE, "HPTN071.input.parser.v4.Rscript"),sep=' ')
PR.SEQGEN.FILECREATOR		<- paste('Rscript',system.file(package=PR.PACKAGE, "PANGEA.SeqGen.createInputFile.Rscript"),sep=' ')
PR.SEQGEN.SIMULATOR			<- paste('Rscript',system.file(package=PR.PACKAGE, "PANGEA.SeqGen.run.v4.Rscript"),sep=' ')
PR.VIRUSTREESIMULATOR		<- system.file(package=PR.PACKAGE, "ext", "VirusTreeSimulator.jar")
PR.SEQGEN					<- system.file(package=PR.PACKAGE, "ext", "seq-gen")
PR.HPTN071.LOWACUTE			<- system.file(package=PR.PACKAGE, "ext", "popart-lowacute")
PR.HPTN071.LOWACUTE.PAR		<- system.file(package=PR.PACKAGE, "ext", "PangeaParamsLowAcute")
PR.HPTN071.HIGHACUTE		<- system.file(package=PR.PACKAGE, "ext", "popart-highacute")
PR.HPTN071.HIGHACUTE.PAR	<- system.file(package=PR.PACKAGE, "ext", "PangeaParamsHighAcute")
PR.VARIOUS					<- paste(PR.STARTME," -exe=VARIOUS",sep='')

HPC.MPIRUN					<- {tmp<- c("mpirun","mpiexec"); names(tmp)<- c("debug","cx1.hpc.ic.ac.uk"); tmp}
HPC.CX1.IMPERIAL			<- "cx1.hpc.ic.ac.uk"		#this is set to system('domainname',intern=T) for the hpc cluster of choice
HPC.MEM						<- "1750mb"
HPC.CX1.IMPERIAL.LOAD		<- "module load intel-suite mpi mafft/7 R/3.2.0 qdist/2.0"
######################################################################################
cmd.hpcsys<- function()
{
	tmp<- system('domainname',intern=T)
	if(!nchar(tmp))	tmp<- "debug"
	tmp
}
##--------------------------------------------------------------------------------------------------------
##	call to various 
##	olli originally written 06-08-2015
##--------------------------------------------------------------------------------------------------------
cmd.various<- function(prog= PR.VARIOUS)
{
	cmd		<- "#######################################################
# start: run VARIOUS
#######################################################"
	cmd		<- paste(cmd, '\n', prog, '\n', sep='')
	cmd		<- paste(cmd,"#######################################################
# end: run VARIOUS
#######################################################\n",sep='')
	cmd
}
##--------------------------------------------------------------------------------------------------------
##	command line generator for 'prog.HPTN071.input.parser.v4'
##	olli originally written 16-08-2015
##--------------------------------------------------------------------------------------------------------
cmd.haircut.call<- function(indir.st, indir.al, outdir, mfile, trainfile=NA, batch.n=NA, batch.id=NA, prog=PR.HAIRCUT.CALL )	
{
	cmd<- "#######################################################
# start: run haircutprog.get.call.for.PNG_ID
#######################################################"
	cmd		<- paste(cmd, paste("\necho \'run ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd, paste(prog,' -mfile=', mfile,' -indir.st=',indir.st,' -indir.al=',indir.al,' -outdir=',outdir, sep=''))
	if(!is.na(trainfile))
		cmd	<- paste(cmd, ' -trainfile=',trainfile, sep='')
	if(!is.na(batch.n) & !is.na(batch.id))
		cmd	<- paste(cmd, ' -batch.n=',batch.n, ' -batch.id=',batch.id, sep='')
	cmd		<- paste('\n',cmd,paste("echo \'end ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd,"#######################################################
# end: run haircutprog.get.call.for.PNG_ID
#######################################################\n",sep='')
	cmd
}
##--------------------------------------------------------------------------------------------------------
##	batch file wrapper
##	olli originally written 26-08-2014
##--------------------------------------------------------------------------------------------------------
cmd.hpcwrapper<- function(cmd, hpcsys= cmd.hpcsys(), hpc.walltime=24, hpc.mem="1750mb", hpc.nproc=1, hpc.q=NA)
{
	wrap<- "#!/bin/sh"
	#hpcsys<- HPC.CX1.IMPERIAL
	if(hpcsys%in%c(HPC.CX1.IMPERIAL,'(none)'))
	{				
		tmp	<- paste("#PBS -l walltime=",hpc.walltime,":59:59,pcput=",hpc.walltime,":45:00",sep='')
		wrap<- paste(wrap, tmp, sep='\n')		
		tmp	<- paste("#PBS -l select=1:ncpus=",hpc.nproc,":mem=",hpc.mem,sep='')
		wrap<- paste(wrap, tmp, sep='\n')
		wrap<- paste(wrap, "#PBS -j oe", sep='\n')
		if(!is.na(hpc.q))
			wrap<- paste(wrap, paste("#PBS -q",hpc.q), sep='\n\n')
		wrap<- paste(wrap, HPC.CX1.IMPERIAL.LOAD, sep='\n')
	}
	else if(hpcsys=='debug')
		cat(paste("\ndetected no HPC system and no hpcwrapper generated, domain name is",hpcsys))
	else
		stop(paste("unknown hpc system with domain name",hpcsys))
	
	cmd<- lapply(seq_along(cmd),function(i){	paste(wrap,cmd[[i]],sep='\n')	})
	if(length(cmd)==1)
		cmd<- unlist(cmd)
	cmd	
}
##--------------------------------------------------------------------------------------------------------
##	batch file caller
##	olli originally written 26-08-2014
##--------------------------------------------------------------------------------------------------------
cmd.hpccaller<- function(outdir, outfile, cmd)
{
	if( nchar( Sys.which("qsub") ) )
	{
		file	<- paste(outdir,'/',gsub(':','',outfile),'.qsub',sep='')
		cat(paste("\nwrite HPC script to",file,"\n"))
		cat(cmd,file=file)
		cmd		<- paste("qsub",file)
		cat( cmd )
		cat( system(cmd, intern=TRUE) )
		Sys.sleep(1)
	}
	else
	{
		file	<- paste(outdir,'/',gsub(':','',outfile),'.sh',sep='')
		cat(paste("\nwrite Shell script to\n",file,"\nStart this shell file manually\n"))
		cat(cmd,file=file)
		Sys.chmod(file, mode = "777")	
		Sys.sleep(1)
	}
	file
}
##--------------------------------------------------------------------------------------------------------
##	olli originally written 22-12-2015
##--------------------------------------------------------------------------------------------------------
#' @title Command line generator for \code{HPTN071.simulator}
#' @return command line string
#' @example example/ex.HPTN071.simulator.R
cmd.HPTN071.simulator<- function(outdir, seed=NA, opt.acute='high', opt.intervention='fast')	
{
	stopifnot(opt.acute%in%c('low','high'), opt.intervention%in%c('none','slow','fast'))
	if(opt.intervention=='fast')
		arg.intervention	<- 1
	if(opt.intervention=='slow')
		arg.intervention	<- 2
	if(opt.intervention=='none')
		arg.intervention	<- 3
	if(opt.acute=='low')
	{
		prog				<- PR.HPTN071.LOWACUTE
		pars				<- PR.HPTN071.LOWACUTE.PAR
		if(is.na(seed))
			seed			<- 2
	}
	if(opt.acute=='high')
	{
		prog				<- PR.HPTN071.HIGHACUTE
		pars				<- PR.HPTN071.HIGHACUTE.PAR
		if(is.na(seed))
			seed			<- 1
	}	 
	outfile.start			<- paste(outdir,'/150129_HPTN071_sc',toupper(substring(opt.acute,1,1)),toupper(substring(opt.intervention,1,1)),sep='')
	cmd<- "#######################################################
# start: run HPTN071.simulator
#######################################################"
	cmd		<- paste(cmd, paste("\necho \'run ",basename(prog),"\'\n",sep=''))	
	cmd		<- paste(cmd, 'cp -R ',pars,' ',basename(pars), ' \n',sep='')	
	cmd		<- paste(cmd, paste(prog,' ', basename(pars),' ',seed,' ', arg.intervention,' 1\n', sep=''))
	cmd		<- paste(cmd,'mv phylogenetic_individualdata* ',outfile.start,'_IND.csv\n',sep='')
	cmd		<- paste(cmd,'mv phylogenetic_transmission* ',outfile.start,'_TRM.csv\n',sep='')
	cmd		<- paste(cmd,'mv Annual_outputs* ',outfile.start,'_EPI.csv\n',sep='')
	cmd		<- paste(cmd,'rm -rf ',basename(pars),' \n',sep='')
	cmd		<- paste(cmd,paste("echo \'end ",basename(prog),"\'\n",sep=''))
	cmd		<- paste(cmd,"#######################################################
# end: run HPTN071.simulator
#######################################################\n",sep='')
	cmd
}

##--------------------------------------------------------------------------------------------------------
##	command line generator for 'prog.HPTN071.input.parser.v4'
##	olli originally written 26-01-2015
##--------------------------------------------------------------------------------------------------------
#' @title Command line generator for \code{HPTN071.input.parser.v2}
#' @return command line string
#' @example example/ex.seq.sampler.v4.R
cmd.HPTN071.input.parser.v4<- function(indir, infile.trm, infile.ind, infile.args, outdir, outfile.trm, outfile.ind, prog=PR.HPTN071.INPUT.PARSER4 )	
{
	cmd<- "#######################################################
# start: run HPTN071.input.parser.v4
#######################################################"
	cmd		<- paste(cmd, paste("\necho \'run ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd, paste(prog,' -indir=', indir,' -infile.args=',infile.args,' -outdir=',outdir, sep=''))
	if(!is.na(infile.ind))
		cmd	<- paste(cmd, ' -infile.ind=',infile.ind, sep='')
	if(!is.na(infile.trm))
		cmd	<- paste(cmd, ' -infile.trm=',infile.trm, sep='')	
	if(!is.na(outfile.ind))
		cmd	<- paste(cmd, ' -outfile.ind=',outfile.ind, sep='')
	if(!is.na(outfile.trm))
		cmd	<- paste(cmd, ' -outfile.trm=',outfile.trm, sep='')	
	cmd		<- paste(cmd,paste("\necho \'end ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd,"#######################################################
# end: run HPTN071.input.parser.v4
#######################################################\n",sep='')
	cmd
}
######################################################################################
#	return command line to run the VirusTreeSimulator	
#	olli originally written 19-08-2014
#	return 		character string 
#' @title Command line generator for the \code{VirusTreeSimulator}
#' @description The \code{VirusTreeSimulator} reads files from the sequence sampler in directory \code{indir} and writes detailed nexus files
#' in directory \code{outdir} for the virus tree simulator. The program generates within-host phylogenies
#' for sampled and unsampled individuals in a transmission chain along the specified within host coalescent
#' model. Within-host phylogenies are then concatenated into a between-host phylogeny for each transmission chain.
#' @return command line string
#' @example example/ex.virus.tree.simulator.R
cmd.VirusTreeSimulator<- function(indir, infile.trm, infile.ind, outdir, outfile, prog=PR.VIRUSTREESIMULATOR, prog.args='-demoModel Logistic -N0 0.1 -growthRate 1.5 -t50 -4')
{
	cmd<- "#######################################################
# start: run VirusTreeSimulator
#######################################################"
	cmd		<- paste(cmd, paste("\necho \'run ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd, paste('java -Xms64m -Xmx400m -jar ',prog,' ',prog.args,'  ', indir,'/',infile.trm,' ',indir,'/',infile.ind,' ',outdir,'/',outfile, '\n', sep=''))
	cmd		<- paste(cmd, 'find ',outdir,' -name "*simple*" -delete\n', sep='')
	cmd		<- paste(cmd,paste("echo \'end ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd,"#######################################################
# end: run VirusTreeSimulator
#######################################################\n",sep='')
	cmd
}
######################################################################################
#	return command line to run SeqGen Input File Creator	
#	olli originally written 26-08-2014
#	return 		character string 
#' @title Command line generator for \code{prog.PANGEA.SeqGen.createInputFile}
#' @description The \code{prog.PANGEA.SeqGen.createInputFile} reads files from the virus tree simulator in directory \code{indir.vts} and writes input files for \code{SeqGen}
#' to directory \code{outdir}. The program reads simulated transmission chain phylogenies with branches in units of calendar time
#' for sampled and unsampled individuals in a transmission chain. Within host evolutionary rates are drawn from a distribution, and
#' within host branch lengths are converted into the expected number of substitutions along the branch. Transmission branches are
#' multiplied with a multiplier to allow for slower evolution between hosts. The multiplier is drawn from a distribution. Starting sequences
#' are drawn from a pool of precomputed sequences. GTR parameters are drawn from a distribution. This is all that s needed to specify 
#' the SeqGen input files for each transmission chain.
#' @return command line string
#' @example example/ex.seqgen.inputfilecreator.R
cmd.SeqGen.createInputFiles<- function(indir.epi, infile.epi, indir.vts, infile.vts, infile.args, outdir, prog=PR.SEQGEN.FILECREATOR)
{
	cmd<- "#######################################################
# start: run SeqGen.createInputFile 
#######################################################"
	cmd		<- paste(cmd, paste("\necho \'run ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd, paste(prog,' -indir.epi=', indir.epi,' -infile.epi=',infile.epi,' -indir.vts=', indir.vts,' -infile.vts=',infile.vts,' -infile.args=',infile.args,' -outdir=',outdir,' \n', sep=''))
	cmd		<- paste(cmd,paste("echo \'end ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd,"#######################################################
# end: run SeqGen.createInputFile
#######################################################\n",sep='')
	cmd
}
######################################################################################
#	return command line to run SeqGen and process SeqGen output 	
#	olli originally written 09-09-2014
#	return 		character string 
#' @title Command line generator for \code{prog.PANGEA.SeqGen.run}
#' @description \code{prog.PANGEA.SeqGen.run} reads file \code{infile.sg} in directory \code{indir.sg} that was
#' created with the \code{SeqGen} input file creator. The simulated partial sequences are collected, coerced back
#' into Gag, Pol, Env genes, and written in fasta format to directory \code{outdir}. Patient Metavariables are 
#' stored in the same directory, and zip files are created.
#' @return command line string
#' @example example/ex.seqgen.run.R
cmd.SeqGen.run<- function(indir.epi, infile.epi, indir.sg, infile.sg, infile.args, outdir, prog=PR.SEQGEN.SIMULATOR)
{
	cmd<- "#######################################################
# start: run SeqGen.run 
#######################################################"
	cmd		<- paste(cmd, paste("\necho \'run ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd, paste(prog,' -indir.epi=', indir.epi,' -infile.epi=',infile.epi, ' -indir.sg=', indir.sg,' -infile.sg=',infile.sg,' -infile.args=',infile.args,' -outdir=',outdir,' \n', sep=''))
	cmd		<- paste(cmd,paste("echo \'end ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd,"#######################################################
# end: run SeqGen.run
#######################################################\n",sep='')
	cmd
}
######################################################################################
#	return command line to run Seq-Gen-1.3.3	
#	olli originally written 26-08-2014
#	return 		character string 
cmd.SeqGen<- function(indir, infile, outdir, outfile, prog=PR.SEQGEN, prog.args='-n1 -k1 -on -z42', alpha=1, gamma=4, invariable=0, scale=1, 
						freq.A=0.25, freq.C=0.25, freq.G=0.25, freq.T=0.25,
						rate.AC=1, rate.AG=1, rate.AT=1, rate.CG=1, rate.CT=1, rate.GT=1)
{
	cmd<- "#######################################################
# start: run SeqGen
#######################################################"
	cmd		<- paste(cmd, paste("\necho \'run ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd, paste(prog,' ',prog.args,' ',sep=''),sep='')
	#	add substitution model
	#cmd		<- paste(cmd, paste("-mHKY -t3.0 -f0.3,0.2,0.2,0.3",sep=''),sep='')
	tmp		<- ifelse(!is.na(gamma) & gamma>0, paste(' -g',gamma,' -a',alpha,' -i',invariable,sep=''),'')
	cmd		<- paste(cmd, paste('-mGTR',tmp,' -s', scale,
									' -f',freq.A,',',freq.C,',',freq.G,',',freq.T,
									' -r',rate.AC, ',', rate.AG, ',', rate.AT, ',', rate.CG, ',', rate.CT, ',', rate.GT, sep=''),sep='')
	#	add I/O
	cmd		<- paste(cmd, paste(' < ', indir,'/',infile,' > ', outdir,'/',outfile, '\n', sep=''))
	#cmd		<- paste(cmd, 'rm ',indir,'/',infile,'\n', sep='')
	cmd		<- paste(cmd,paste("echo \'end ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd,"#######################################################
# end: run SeqGen
#######################################################\n",sep='')
	cmd
}


