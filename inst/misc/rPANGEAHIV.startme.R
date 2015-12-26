#! /Library/Frameworks/R.framework/Versions/2.15/Resources/bin/Rscript
##
##	first line in shell script starts with #! and points to absolute path to Rscript
##	CHANGE  as needed
##
##! /apps/R/2.15/lib64/R/bin/Rscript
###############################################################################
#
#	project scripts that can be run from command line, without re-building the package all the time,
# 	because the R files are re-loaded below
#
# usage from R:
#> setwd("/Users/Oliver/git/HPTN071sim/source/rPANGEAHIVsim"); source("misc/rPANGEAHIV.startme.R")
#
###############################################################################
args <- commandArgs()
if(!any(args=='--args'))
	args<- vector("numeric",0)
if(any(args=='--args'))
	args<- args[-(1:match("--args", args)) ]
#the package directory (local working copy of the code, not the installed package directory within the R directory 
CODE.HOME	<<- "/Users/Oliver/git/HPTN071sim/source/rPANGEAHIVsim"

#the home directory of all projects
HOME		<<- "/Users/Oliver/git/HPTN071sim/"
#HOME		<<- "/work/or105/UKCA_1309"
#HOME		<<- "/work/or105/ATHENA_2013"
DATA		<<- paste(HOME,"data",sep='/')

DEBUG		<<- 0		#If 1, a dump file is created that can be loaded and computations can be inspected at the point of error.
LIB.LOC		<<- NULL
#LIB.LOC	<<- paste(CODE.HOME,"../",sep='')
EPS			<<- 1e-12	#Machine precision	

#the default script to be called if -exe is not specified on the command line	
default.fun		<- 'pipeline.HPTN071'
###############################################################################
#	select script specified with -exe on the command line. If missing, start default script 'default.fun'.
argv<- list()
if(length(args))
{
	tmp<- na.omit(sapply(args,function(arg)
					{
						switch(substr(arg,2,4),
								exe= return(substr(arg,6,nchar(arg))),
								NA)
					}))
	if(length(tmp)!=0)
	{
		if(length(tmp)>1) stop("hivclu.startme.R: duplicate -exe")
		else default.fun<- switch(tmp[1],
					ROXYGENIZE				= "package.roxygenize",
					MAKE.RDATA				= "package.generate.rdafiles",					
					HPTN071.INPUT.PARSER1	= "prog.HPTN071.input.parser.v1",
					HPTN071.INPUT.PARSER2	= "prog.HPTN071.input.parser.v2",
					DSPS.INPUT.PARSER2		= "prog.DSPS.input.parser.v2",
					PR.SEQGEN.FILECREATOR	= "prog.PANGEA.SeqGen.createInputFile",
					PR.SEQGEN.SIMULATOR		= "prog.PANGEA.SeqGen.run"
					)
	}
	tmp<- na.omit(sapply(args,function(arg)
					{
						switch(substr(arg,2,10),
								code.home= return(substr(arg,12,nchar(arg))),
								NA)
					}))	
	if(length(tmp)!=0)	CODE.HOME<<- tmp[1]
	tmp<- na.omit(sapply(args,function(arg)
					{
						switch(substr(arg,2,6),
								debug= 1,
								NA)
					}))		
	if(length(tmp)!=0)	DEBUG<<- tmp[1]	
	argv<<- args
}
###############################################################################
.ls.objects <- function (pos = 1, pattern, order.by,
		decreasing=FALSE, head=FALSE, n=5) {
	napply <- function(names, fn) sapply(names, function(x)
					fn(get(x, pos = pos)))
	names <- ls(pos = pos, pattern = pattern)
	obj.class <- napply(names, function(x) as.character(class(x))[1])
	obj.mode <- napply(names, mode)
	obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
	obj.prettysize <- napply(names, function(x) {
				capture.output(print(object.size(x), units = "auto")) })
	obj.size <- napply(names, object.size)
	obj.dim <- t(napply(names, function(x)
						as.numeric(dim(x))[1:2]))
	vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
	obj.dim[vec, 1] <- napply(names, length)[vec]
	out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
	names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
	if (!missing(order.by))
		out <- out[order(out[[order.by]], decreasing=decreasing), ]
	if (head)
		out <- head(out, n)
	out
}
# from: http://stackoverflow.com/questions/1358003/tricks-to-manage-the-available-memory-in-an-r-session
lsos <- function(..., n=10) {
	.ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}
my.mkdir<-function(root,data.name)
{
	if(length(dir(root,pattern=paste('^',data.name,'$',sep='')))==0)
		system(paste("mkdir ",paste(root,data.name,sep='/'),sep=''))
}
my.dumpframes<- function()
{
	geterrmessage()
	dump.frames()
	cat(paste("\nmy.dumpframes dump 'last.dump' to file",paste(DATA,paste("debug_",paste(strsplit(date(),' ')[[1]],collapse='_'),".rda\n",sep=''),sep='/')))
	save(last.dump, file=paste(DATA,paste("debug_",paste(strsplit(date(),' ')[[1]],collapse='_'),".rda",sep=''),sep='/'))
	q()
}
###############################################################################
#	re-load all R files
require(rPANGEAHIVsim)
print(CODE.HOME)
function.list<-c(list.files(path= paste(CODE.HOME,"R",sep='/'), pattern = ".R$", all.files = FALSE,
				full.names = TRUE, recursive = FALSE),paste(CODE.HOME,"misc","rPANGEAHIV.prjcts.R",sep='/'))
sapply(function.list,function(x){ source(x,echo=FALSE,print.eval=FALSE, verbose=FALSE) })
###############################################################################
#	run script
#stop()
if(DEBUG)	options(error= my.dumpframes)	
cat(paste("\nrPANGEAHIV: ",ifelse(DEBUG,"debug",""),"call",default.fun,"\n"))
do.call(default.fun,list()) 	
cat("\nrPANGEAHIV: ",ifelse(DEBUG,"debug","")," end\n")
