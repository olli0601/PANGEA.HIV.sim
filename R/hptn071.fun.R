##--------------------------------------------------------------------------------------------------------
#	olli copied from hivclust
##--------------------------------------------------------------------------------------------------------
seq.unique<- function(seq.DNAbin.matrix)
{
	x<- as.character(seq.DNAbin.matrix)
	x<- apply(x, 1, function(z) paste(z,collapse=''))
	seq.DNAbin.matrix[!duplicated(x),]			
}
##--------------------------------------------------------------------------------------------------------
#	olli copied from hivclust
##--------------------------------------------------------------------------------------------------------
seq.rmgaps<- function(seq.DNAbin.matrix, rm.only.col.gaps=1, rm.char='-', verbose=0)
{
	if(class(seq.DNAbin.matrix)=='DNAbin')
		seq.DNAbin.matrix		<- as.character(seq.DNAbin.matrix)		
	if(!rm.only.col.gaps)
	{	
		if(is.matrix(seq.DNAbin.matrix))
		{
			tmp					<- lapply(seq_len(nrow(seq.DNAbin.matrix)), function(i){	seq.DNAbin.matrix[i, !seq.DNAbin.matrix[i,]%in%rm.char]	})
			names(tmp)			<- rownames(seq.DNAbin.matrix)
		}
		else
		{
			tmp					<- lapply(seq_along(seq.DNAbin.matrix), function(i){	seq.DNAbin.matrix[[i]][ !seq.DNAbin.matrix[[i]]%in%rm.char]	})
			names(tmp)			<- names(seq.DNAbin.matrix)
		}		
		seq.DNAbin.matrix	<- tmp
	}
	else
	{		
		gap					<- apply(seq.DNAbin.matrix,2,function(x) all(x%in%rm.char)) 
		if(verbose)		cat(paste("\nremove gaps, n=",length(which(gap))))
		if(verbose>1)	cat(paste("\nremove gaps, at pos=",which(gap)))
		seq.DNAbin.matrix	<- seq.DNAbin.matrix[,!gap]	
	}
	as.DNAbin( seq.DNAbin.matrix )
}
##--------------------------------------------------------------------------------------------------------
#	olli copied from hivclust
##--------------------------------------------------------------------------------------------------------
seq.singleton2bifurcatingtree<- function(ph.s, dummy.label=NA)
{	
	if(!Nnode(ph.s))
	{
		stopifnot(!nrow(ph.s$edge))
		if(is.na(dummy.label))
			dummy.label		<- paste('DUMMY',ph.s$tip.label, sep='_')
		ph.s$edge			<- matrix(c(3,1,3,2), nrow=2, ncol=2, byrow=TRUE) 	
		ph.s$edge.length	<- c(ph.s$root.edge,0)
		ph.s$tip.label		<- c(ph.s$tip.label, dummy.label)
		ph.s$root.edge		<- 0
		ph.s$Nnode			<- 1
	}
	ph.s
}
##--------------------------------------------------------------------------------------------------------
#	olli copied from hivclust
##--------------------------------------------------------------------------------------------------------
seq.addrootnode<- function(ph, dummy.label)
{	
	ph$tip.label	<- c(ph$tip.label, dummy.label)
	tmp				<- which( ph$edge[,1]>=Ntip(ph) )
	ph$edge[tmp,1]	<- ph$edge[tmp,1] + 2 
	tmp				<- which( ph$edge[,2]>=Ntip(ph) )
	ph$edge[tmp,2]	<- ph$edge[tmp,2] + 2 
	if(Nnode(ph)>0)
		ph$edge		<- rbind( matrix(c(Ntip(ph)+1,Ntip(ph)+1,Ntip(ph),Ntip(ph)+2), nrow=2, ncol=2), ph$edge)
	if(Nnode(ph)==0)
		ph$edge		<- rbind( matrix(c(Ntip(ph)+1,Ntip(ph)+1,Ntip(ph),Ntip(ph)-1), nrow=2, ncol=2), ph$edge)
	ph$edge.length	<- c(0,ph$root.edge, ph$edge.length)
	ph$root.edge	<- 0
	ph$Nnode		<- ph$Nnode+1
	ph
}
##--------------------------------------------------------------------------------------------------------
#	olli copied from hivclust
##--------------------------------------------------------------------------------------------------------
seq.collapse.singles<- function (tree) 
{
	elen 		<- tree$edge.length
	xmat 		<- tree$edge
	node.lab 	<- tree$node.label
	nnode 		<- tree$Nnode
	ntip 		<- length(tree$tip.label)
	root		<- 0
	singles 	<- NA
	while (length(singles) > 0) 
	{
		tx <- tabulate(xmat[, 1])
		singles <- which(tx == 1)
		if (length(singles) > 0) 
		{
			i 					<- singles[1]
			prev.node 			<- which(xmat[, 2] == i)
			next.node 			<- which(xmat[, 1] == i)
			xmat[prev.node, 2] 	<- xmat[next.node, 2]
			xmat 				<- xmat[xmat[, 1] != i, , drop=0]
			xmat[xmat > i] 		<- xmat[xmat > i] - 1L
			if(!length(prev.node))
				root			<- root + elen[next.node]
			if(length(prev.node))
				elen[prev.node] <- elen[prev.node] + elen[next.node]
			if (!is.null(node.lab)) 
				node.lab <- node.lab[-c(i - ntip)]
			nnode <- nnode - 1L
			elen <- elen[-next.node]
		}
	}
	tree$edge 			<- xmat
	tree$edge.length 	<- elen
	tree$node.label 	<- node.lab
	tree$Nnode 			<- nnode
	tree$root.edge		<- root
	tree
}
##--------------------------------------------------------------------------------------------------------
#	olli copied from hivclust
##--------------------------------------------------------------------------------------------------------
seq.read.newick<- function (file = "", text) 
{
	if (file != "") 
		text <- scan(file, sep = "\n", what = "character")	
	Nnode		<- length(gregexpr("\\(", text)[[1]])
	Ntip		<- 1+length(gregexpr(",", text)[[1]])	
	tree 		<- unlist(strsplit(text, NULL))
	tip.label 	<- vector(mode = "character")
	node.label	<- vector(mode = 'character')
	edge 		<- matrix(data = 0, Nnode + Ntip - 1, 2)
	edge.length <- rep(0, Nnode + Ntip - 1)
	ntip 		<- vector(mode = "numeric")
	currnode 	<- Ntip + 1
	nodecount 	<- currnode
	i 	<- 1
	j 	<- 1
	k 	<- 1
	while(tree[i] != ";") 
	{
		if(tree[i] == "(") 
		{
			edge[j, 1] <- currnode
			i <- i + 1
			if(is.na(match(tree[i], c("(", ")", ",", ":", ";")))) 
			{
				l				<- gregexpr(",|:|\\)", substr(text, i, nchar(text)))[[1]][1]
				stopifnot(l>0)
				tip.label[k] 	<- substr(text, i, i+l-2)	
				i				<- i+l-1
				edge[j, 2] 		<- k
				k 				<- k + 1
				ntip[j] 		<- 1
				if (tree[i] == ":") 
				{
					i				<- i + 1
					l				<- gregexpr(",|\\)", substr(text, i, nchar(text)))[[1]][1]
					stopifnot(l>0)
					edge.length[j]	<- as.numeric(substr(text, i, i+l-2))
					i				<- i+l-1					
				}
			}
			else if(tree[i] == "(") 
			{
				nodecount 	<- nodecount + 1
				currnode 	<- nodecount
				edge[j, 2] 	<- currnode
			}
			j <- j + 1
		}
		else if(tree[i] == ")") 
		{
			i <- i + 1			
			if(is.na(match(tree[i], c(":", ")")))) 
			{
				l	<- gregexpr(":|;|\\)", substr(text, i, nchar(text)))[[1]][1]
				stopifnot(l>0)
				node.label[currnode-Ntip]	<- substr(text, i, i+l-2)
				i	<- i+l-1				
			}
			if(tree[i] == ":") 
			{
				i 	<- i + 1
				l	<- gregexpr(",|\\)", substr(text, i, nchar(text)))[[1]][1]
				stopifnot(l>0)
				edge.length[match(currnode, edge[, 2])] <- as.numeric(substr(text, i, i+l-2))
				i	<- i+l-1	
			}
			ntip[match(currnode, edge[, 2])] 	<- sum(ntip[which(edge[, 1] == currnode)])
			currnode 							<- edge[match(currnode, edge[, 2]), 1]
		}
		else if(tree[i]==",") 
		{
			edge[j, 1] 	<- currnode
			i 			<- i + 1
			if(is.na(match(tree[i], c("(", ")", ",", ":", ";")))) 
			{
				l				<- gregexpr(",|:|\\)", substr(text, i, nchar(text)))[[1]][1]
				stopifnot(l>0)
				tip.label[k] 	<- substr(text, i, i+l-2)	
				i				<- i+l-1
				edge[j, 2] 		<- k
				k 				<- k + 1
				ntip[j] 		<- 1
				if (tree[i] == ":") 
				{
					i				<- i + 1
					l				<- gregexpr(",|\\)", substr(text, i, nchar(text)))[[1]][1]
					stopifnot(l>0)
					edge.length[j]	<- as.numeric(substr(text, i, i+l-2))
					i				<- i+l-1
				}
			}
			else if (tree[i] == "(") 
			{
				nodecount 	<- nodecount + 1
				currnode 	<- nodecount
				edge[j, 2] 	<- currnode
			}
			j <- j + 1
		}
	}
	tmp	<- which(edge[,1]==0)
	if(length(tmp))
	{
		edge		<- edge[-tmp,]
		edge.length	<- edge.length[-tmp]		
		tmp			<- sort( unique( as.numeric( edge ) ) )
		tmp			<- rbind(tmp, seq_along(tmp))		
		tmp			<- sapply( as.numeric( edge ), function(j)	tmp[2, match(j, tmp[1,])] )
		edge		<- matrix(tmp, ncol=2)
	}
	phy <- list(edge = edge, Nnode = as.integer(Nnode), tip.label = tip.label, Ndesc = ntip)
	if(sum(edge.length) > 1e-08) 
		phy$edge.length	<- edge.length
	if(length(node.label))
		phy$node.label	<- node.label
	class(phy) <- "phylo"
	return(phy)
}
##--------------------------------------------------------------------------------------------------------
#	olli copied from hivclust
##--------------------------------------------------------------------------------------------------------
hivc.beast2out.read.nodestats <- function(bstr) 
{
	#	remove anything before first '('
	bstr	<- regmatches(bstr, regexpr('\\(.*',bstr))
	# 	store meta info for inner nodes that is given in [], and not in :[] which is meta info for edges	
	tmp		<- unlist(regmatches(bstr,gregexpr('[^:]\\[[^]]+',bstr)))
	stopifnot(length(tmp)>0)
	tmp		<- sapply( tmp, function(x) substr(x, 4, nchar(x)) ) 
	#	for each inner node, extract stats
	tmp		<- strsplit(tmp, ',')
	tmp		<- lapply(seq_along(tmp), function(i)
			{
				z<- strsplit(tmp[[i]],'=')				
				data.table(NODE_PARSE_ID=i, STAT=sapply(z,'[',1), VALUE=sapply(z,'[',2))
			})
	node.stat	<- do.call('rbind', tmp)
	tmp			<- node.stat[, unique(STAT)]
	cat(paste('\nFound node statistics=',paste(tmp,collapse=' ')))
	tmp			<- node.stat[, list(has.all.stats= !length(setdiff(tmp, STAT))  ) , by='NODE_PARSE_ID']
	tmp			<- subset(tmp, !has.all.stats)[, NODE_PARSE_ID]
	cat(paste('\nSome statistics missing for nodes=',paste(tmp,collapse=' ')))
	node.stat 
}
##--------------------------------------------------------------------------------------------------------
#	olli copied from hivclust
##--------------------------------------------------------------------------------------------------------
hivc.beast2out.read.nodeidtree <- function(bstr, method.node.stat='any.node') 
{
	# strip all meta variables and ; at end
	bstr		<- gsub("\\[[^]]*\\]", "", bstr)
	bstr		<- gsub(';','',bstr)
	# for each node, add a dummy node label NODE_PARSE_IDxx	
	dummy.tree	<- unlist(strsplit(bstr, ":"))
	if(method.node.stat=='inner.node')
	{
		#	interior branch length: 	previous index ends in ). so tmp is the index of the dummy.tree chunks that gives the start of a branch length of an inner node
		tmp			<- which( c(FALSE, grepl(')$',dummy.tree)[-length(dummy.tree)]) )
		#	prepend NODE_PARSE_IDxx before the branch length of an inner node
		tmp			<- tmp-1			
	}
	if(method.node.stat=='any.node')
		tmp			<- seq_along(dummy.tree)
	dummy.tree	<- sapply(seq_along(dummy.tree), function(i)
			{
				z<- which(i==tmp)
				ifelse(length(z),	paste(dummy.tree[i],'NODE_PARSE_ID',z,sep=''),	dummy.tree[i] )
			}) 			
	dummy.tree	<- paste(dummy.tree, collapse=':',sep='')
	dummy.tree	<- regmatches(dummy.tree, regexpr('\\(.*',dummy.tree))
	dummy.tree	<- paste(dummy.tree, ';', sep='')	
	ph			<-  seq.read.newick(text=dummy.tree)
	ph
}
##--------------------------------------------------------------------------------------------------------
#	olli copied from hivclust
##--------------------------------------------------------------------------------------------------------
hivc.beast2out.read.nexus.and.stats<- function(file, tree.id=NA, method.node.stat='any.node') 
{	
	stopifnot(method.node.stat%in%c('any.node','inner.node'))
	
	X				<- scan(file = file, what = "", sep = "\n", quiet = TRUE)	
	#	read TRANSLATE chunk
	X.endblock		<- grep("END;|ENDBLOCK;|End;", X, ignore.case = TRUE)
	X.semico 		<- grep(";", X)
	X.i1 			<- grep("BEGIN TREES;|Begin trees;", X, ignore.case = TRUE)
	X.i2 			<- grep("TRANSLATE|Translate", X, ignore.case = TRUE)	
	tmp 			<- X.semico[X.semico > X.i2][1]
	tmp 			<- X[(X.i2 + 1):tmp]
	tmp				<- gsub('[,;]$','',gsub('^\\s+','',tmp))
	tmp				<- tmp[nzchar(tmp)]
	tmp				<- strsplit(tmp, ' ')
	df.translate	<- data.table(NEXUS_ID= sapply(tmp, '[[', 1), NEXUS_LABEL=sapply(tmp, '[[', 2) )
	set(df.translate, NULL, 'NEXUS_LABEL', df.translate[, gsub("\'","",NEXUS_LABEL)])
	cat(paste('\nFound taxa, n=', nrow(df.translate)))
	
	if(!is.na(tree.id))
	{
		#	read one newick tree with id 'tree.id'
		bstr		<- X[grep(paste(tree.id,"[[:space:]]+",sep=''), X)]
		node.stat	<- hivc.beast2out.read.nodestats(bstr)
		cat(paste('\nFound node statistics, n=', nrow(node.stat)))
		set(node.stat, NULL, 'tree.id', tree.id[i] )		
		btree		<- hivc.beast2out.read.nodeidtree(bstr, method.node.stat=method.node.stat) 
		#
		# link node.stats with tree nodes (tip + inner node)
		# NODE_ID is index of node in 'btree' phylo object
		#
		tmp			<- strsplit( btree$tip.label, 'NODE_PARSE_ID' )
		df.link		<- data.table(NODE_ID=seq_along(btree$tip.label), NEXUS_ID=sapply(tmp,'[[',1), NODE_PARSE_ID=sapply(tmp,'[[',2))
		df.link		<- merge(df.link, df.translate, by='NEXUS_ID')
		cat(paste('\nFound tree tips with taxon name, n=', nrow(df.link)))
		tmp			<- strsplit( btree$node.label, 'NODE_PARSE_ID' )
		tmp			<- data.table(NODE_ID=Ntip(btree)+seq_along(btree$node.label), NODE_PARSE_ID=sapply(tmp,'[[',2), NEXUS_LABEL=NA_character_)
		df.link		<- rbind(subset(df.link,select=c(NODE_ID, NODE_PARSE_ID, NEXUS_LABEL)), tmp)
		set(df.link,NULL,'NODE_PARSE_ID',df.link[, as.integer(NODE_PARSE_ID)])
		set(df.link,NULL,'NODE_ID',df.link[, as.integer(NODE_ID)])
		set(df.link,NULL,'TREE_ID',tree.id)
		node.stat	<- merge( node.stat, subset(df.link, select=c(NODE_PARSE_ID, NODE_ID, TREE_ID)), by='NODE_PARSE_ID' )
		set(node.stat,NULL,'NODE_PARSE_ID',NULL)
		cat(paste('\nLinked node statistics to tree nodes, n=', nrow(node.stat)))
		#
		# set tip.labels and rm node.labels
		#
		setkey(df.link, NODE_ID)
		btree$tip.label		<- df.link[seq_len(Ntip(btree)),][,NEXUS_LABEL]
		btree$node.label	<- NULL		
	}
	if(is.na(tree.id))
	{
		#	read all newick trees in nexus file
		tmp			<- regexpr('^tree\\s\\S+',X)
		tree.id		<- sapply( regmatches(X,tmp), function(x) substr(x, 5, nchar(x)))
		tree.id		<- gsub('\\s','',tree.id)		
		cat(paste('\nFound tree id=', paste(tree.id, collapse=' ')))
		X			<- X[ which(tmp>0) ]
		cat(paste('\nFound trees, n=',length(tree.id)))
		node.stat	<- lapply(seq_along(tree.id), function(i)
				{
					
					bstr	<- X[grep(paste(tree.id[i],"[[:space:]]+",sep=''), X)]
					cat(paste('\nGet node statistics for tree id=',tree.id[i]))
					tmp		<- hivc.beast2out.read.nodestats(bstr)
					set(tmp, NULL, 'TREE_ID', tree.id[i] )
					tmp
				})
		if(length(node.stat)>1)
			node.stat	<- do.call('rbind',node.stat)
		if(length(node.stat)==1)
			node.stat	<- node.stat[[1]]
		suppressWarnings({ node.stat[, NODE_ID:=NA_integer_] })				
		setkey(node.stat, TREE_ID, NODE_PARSE_ID)
		btree		<- vector('list',length(tree.id))
		for(i in seq_along(tree.id))
		{
			bstr		<- X[grep(paste(tree.id[i],"[[:space:]]+",sep=''), X)]
			cat(paste('\nRead tree for tree id=',tree.id[i]))
			btree.i		<- hivc.beast2out.read.nodeidtree(bstr, method.node.stat=method.node.stat)
			#
			# link node.stats with tree nodes (tip + inner node)
			# NODE_ID is index of node in 'btree.i' phylo object
			#
			tmp			<- strsplit( btree.i$tip.label, 'NODE_PARSE_ID' )
			df.link		<- data.table(NODE_ID=seq_along(btree.i$tip.label), NEXUS_ID=sapply(tmp,'[[',1), NODE_PARSE_ID=sapply(tmp,'[[',2))
			df.link		<- merge(df.link, df.translate, by='NEXUS_ID')
			cat(paste('\nFound tree tips with taxon name, n=', nrow(df.link)))
			tmp			<- strsplit( btree.i$node.label, 'NODE_PARSE_ID' )
			tmp			<- data.table(NODE_ID=Ntip(btree.i)+seq_along(btree.i$node.label), NODE_PARSE_ID=sapply(tmp,'[[',2), NEXUS_LABEL=NA_character_)
			df.link		<- rbind(subset(df.link,select=c(NODE_ID, NODE_PARSE_ID, NEXUS_LABEL)), tmp)
			set(df.link,NULL,'NODE_PARSE_ID',df.link[, as.integer(NODE_PARSE_ID)])
			set(df.link,NULL,'NODE_ID',df.link[, as.integer(NODE_ID)])		
			for(j in seq_len(nrow(df.link)))
				set(node.stat, node.stat[, which(TREE_ID==tree.id[i] & NODE_PARSE_ID==df.link[j,NODE_PARSE_ID])], 'NODE_ID', df.link[j,NODE_ID])		
			tmp			<- node.stat[, length(which(!is.na(NODE_ID)))]
			cat(paste('\nTotal linked node statistics to tree nodes, n=', tmp  ))
			#
			# set tip.labels and rm node.labels
			#
			setkey(df.link, NODE_ID)
			btree.i$tip.label	<- df.link[seq_len(Ntip(btree.i)),][,NEXUS_LABEL]
			btree.i$node.label	<- NULL	
			btree[[i]]			<- btree.i
		}
		if(length(btree)>=2)
		{
			names(btree)	<- tree.id
			class(btree)	<- "multiPhylo"			
		}
		if(length(btree)<2)
			btree	<- btree[[1]]
		tmp				<- node.stat[, length(which(is.na(NODE_ID)))]
		cat(paste('\nTotal unlinked node statistics [should be zero], n=', tmp  ))
		set(node.stat,NULL,'NODE_PARSE_ID',NULL)
	}
	list(tree=btree, node.stat=node.stat)	 
}
##--------------------------------------------------------------------------------------------------------
#	return distribution of GTR parameters	
#	olli originally written 20-09-2014
##--------------------------------------------------------------------------------------------------------
#' @title Create data.table of GTR parameters
#' @description Returns a data.table of GTR parameters. 
#' @return data.table
PANGEA.GTR.params<- function()
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
	#	check that relative rates have mean 1
	stopifnot( log.df[, list(CHECK=sum(mu)), by=c('FILE','GENE','state')][, !any(abs(CHECK-3)>2*1e-12)] )
	#	get mean rate. need to be a bit careful here: rate is per site, so length of gene does not matter
	#	but we may have different numbers of samples for each gene
	tmp		<- log.df[, list(meanRate=mean(meanRate)), by='GENE'][, mean(meanRate)]	
	#	ignore variation in meanRate by gene, keep variation across codon_pos for each state
	log.df	<- merge(log.df, log.df[, list(mu.gene= mean(meanRate)/tmp), by='GENE'], by='GENE')	
	set(log.df, NULL, 'mu', log.df[, mu*mu.gene])
	set(log.df, NULL, 'meanRate', tmp)
	set(log.df, NULL, 'mu.gene', NULL)
	log.df
}
##--------------------------------------------------------------------------------------------------------
#	simulate imports during the epidemic	
#	olli originally written 27-01-2015
##--------------------------------------------------------------------------------------------------------
PANGEA.ImportSimulator.SimulateStartingTimeOfIndexCase.v2<- function(df.ind, df.trm, index.starttime.mode='normal')
{
	stopifnot(grepl('normal|fix|shift', index.starttime.mode))
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
	if( grepl('fix',index.starttime.mode) ) 
	{
		tmp2		<- as.numeric(substring(index.starttime.mode, 4))
		stopifnot(tmp2<=1980)
		cat(paste('\nUsing index.starttime.mode rep=', tmp2))
		set(tmp, NULL, 'IDTR_TIME_INFECTED.new', rep(tmp2, nrow(tmp)) )
	}	
	if( grepl('shift',index.starttime.mode))
	{
		tmp2		<- as.numeric(substring(index.starttime.mode, 6))
		cat(paste('\nUsing index.starttime.mode rep=', tmp2))
		set(tmp, NULL, 'IDTR_TIME_INFECTED.new', rep(tmp2, nrow(tmp)) )		
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
#	simulates seroprevalence survey of a proportion of a population	
#	olli originally written 29-01-2015
##--------------------------------------------------------------------------------------------------------
PANGEA.SeroPrevalenceSurvey<- function(df.inds, epi.adult=13, s.INTERVENTION.start=2015, sp.prop.of.sexactive=0.05, sp.times=c(12, 8, 4, 0), verbose=1, test=0)
{		
	df.sp	<- lapply( s.INTERVENTION.start-sp.times, function(yr)
			{
				yr			<- yr-1/365
				sxon.all	<- which( (df.inds[['DOB']]+epi.adult)<=yr  &  df.inds[['DOD']]>yr )
				sxon.sp		<- sort(sample(sxon.all, round(length(sxon.all)*sp.prop.of.sexactive) ))
				if(verbose)
				{
					cat(paste('\nSero prevalence survey in year', yr))
					cat(paste('\nTotal number of individuals=', length(sxon.all)))
					cat(paste('\nProp to include in survey=', sp.prop.of.sexactive))
					cat(paste('\nTotal number of individuals in survey=', length(sxon.sp)))
				}
				df.sp		<- df.inds[sxon.sp, ]				
				df.sp[, YR:=yr]
				df.sp[, AGE:= df.sp[, cut(YR-DOB, breaks=c(epi.adult-1, 16, 20, 25, 30, 35, 40, 50, 60, Inf))]]
				df.sp
			})
	df.sp	<- do.call('rbind', df.sp)	
	if(test)
		df.sp	<- df.sp[, list( ALIVE_N=length(IDPOP), ALIVE_AND_HIV_N=length(which(TIME_TR<=YR)), ALIVE_AND_DIAG_N=length(which(DIAG_T<=YR)), ALIVE_AND_ART_N=length(which(ART1_T<=YR)), ALIVE_AND_SEQ_N=length(which(TIME_SEQ<=YR)) ), by=c('YR', 'GENDER', 'AGE')]
	if(!test)
		df.sp	<- df.sp[, list( ALIVE_N=length(IDPOP), ALIVE_AND_DIAG_N=length(which(DIAG_T<=YR)), ALIVE_AND_ART_N=length(which(ART1_T<=YR)), ALIVE_AND_SEQ_N=length(which(TIME_SEQ<=YR)) ), by=c('YR', 'GENDER', 'AGE')]	
	setkey(df.sp, YR, GENDER, AGE)
	df.sp
}
##--------------------------------------------------------------------------------------------------------
#	simulate imports during the epidemic	
#	olli originally written 11-09-2014
##--------------------------------------------------------------------------------------------------------
PANGEA.ImportSimulator.SimulateIndexCase<- function(df.ind, df.trm, epi.import)
{	
	#	model imports by re-assigning a fraction of infecteds as 'index cases'
	#	assume there are no imports at the moment, and that we know the time of infection of the imports
	tmp		<- df.trm[, which(IDTR>0)]	  
	cat(paste('\nFound transmissions within population, n=', length(tmp)))
	tmp2	<- round(nrow(df.trm)*epi.import) - ( nrow(df.trm)-length(tmp)-1 )
	tmp2	<- max(tmp2, 0)
	cat(paste('\nRe-setting infecteds as index cases after imports, n=', tmp2))
	stopifnot(length(tmp)>1)
	stopifnot(length(tmp2)>=0)
	tmp2	<- as.integer(sample( tmp, tmp2, replace=FALSE ))	
	#	update df.trm
	setkey(df.trm, TIME_TR)
	set(df.trm, tmp2, 'IDTR', df.trm[, min(IDTR)]-rev(seq_along(tmp2)) )
	if('TR_ACUTE'%in%colnames(df.trm))
		set(df.trm, df.trm[, which(IDTR<0)], 'TR_ACUTE', NA_character_)
	#	update df.ind
	tmp		<- subset(df.trm, select=c(IDTR, IDTR_TIME_INFECTED))
	setnames(tmp, c('IDTR', 'IDTR_TIME_INFECTED'), c('IDPOP', 'TIME_TR'))
	tmp2	<- subset(df.trm, select=c(IDREC, TIME_TR))
	setnames(tmp2, 'IDREC', 'IDPOP')
	tmp		<- rbind( tmp, tmp2 )
	setkey(tmp, IDPOP)
	df.ind	<- merge( df.ind, unique(tmp), by=c('IDPOP','TIME_TR'), all.x=TRUE, all.y=TRUE )
	#	
	cat(paste('\nProportion of imported transmissions, p=', (nrow(subset(df.trm, IDTR<0))-1)/nrow(df.trm) ))
	stopifnot( length(setdiff(df.trm[, IDTR], df.ind[, IDPOP]))==0 )
	stopifnot( length(setdiff(df.trm[, IDREC], df.ind[, IDPOP]))==0 )
	list(df.ind=df.ind, df.trm=df.trm)
}
##--------------------------------------------------------------------------------------------------------
#	simulate guide to sequence sampling times. if not NA, then in every year, individuals in +-6mo to the guide are sampled	
#	olli originally written 26-01-2015
##--------------------------------------------------------------------------------------------------------
PANGEA.Seqsampler.SimulateGuideToSamplingTimes.v2<- function(df.ind, seqtime.mode)
{
	stopifnot(grepl('DUnif|Exp|AtDiag|AtART|AtTrm|AtYear',seqtime.mode))
	cat(paste('\nUsing seqtime.mode=', seqtime.mode ))
	if(grepl('Exp',seqtime.mode))
	{
		tmp		<- as.numeric(substring(seqtime.mode,4))
		stopifnot(is.finite(tmp))
		cat(paste('\nPANGEA.Seqsampler.SimulateGuideToSamplingTimes.v2: avg time to sequencing since diagnosis=', tmp))		
		df.ind[, T1_SEQ:= rexp(nrow(df.ind), rate=1/tmp) + df.ind[,DIAG_T]]
		#	need T1_SEQ even when no diagnosis for archival samples
		tmp		<- df.ind[, which(is.na(T1_SEQ) & DOD-TIME_TR>0.5)]
		set( df.ind, tmp, 'T1_SEQ', df.ind[tmp, runif(length(tmp), TIME_TR+0.5, DOD)] )		
		tmp		<- df.ind[, which( is.na(DIAG_T) & T1_SEQ>min(DIAG_T, na.rm=1) )]
		set( df.ind, tmp, 'T1_SEQ', NA_real_)
		tmp	<- df.ind[, which(T1_SEQ>ART1_T)]
		set(df.ind, tmp, 'T1_SEQ', df.ind[tmp,ART1_T])		
	}
	if(grepl('AtYear',seqtime.mode))
	{
		tmp		<- as.numeric(substring(seqtime.mode,7))
		stopifnot(is.finite(tmp))
		cat(paste('\nPANGEA.Seqsampler.SimulateGuideToSamplingTimes.v2: sequencing years after infection=', tmp))		
		df.ind[, T1_SEQ:= tmp + df.ind[,TIME_TR]]
		tmp		<- df.ind[, which(!is.na(T1_SEQ) & DOD<T1_SEQ)]
		set(df.ind, tmp, 'T1_SEQ', df.ind[tmp, (DOD-TIME_TR)*.9 + TIME_TR])
	}
	if(grepl('DUnif',seqtime.mode))
	{
		tmp		<- as.numeric(substring(seqtime.mode,6))
		stopifnot(is.finite(tmp))
		cat(paste('\nPANGEA.Seqsampler.SimulateGuideToSamplingTimes.v2: max time to sequencing since diagnosis=', tmp))		
		df.ind[, T1_SEQ:= runif(nrow(df.ind), min=0, max=tmp) + df.ind[,DIAG_T]]
		#	need T1_SEQ even when no diagnosis for archival samples
		tmp		<- df.ind[, which(is.na(T1_SEQ) & DOD-TIME_TR>0.5)]
		set( df.ind, tmp, 'T1_SEQ', df.ind[tmp, runif(length(tmp), TIME_TR+0.5, DOD)] )		
		tmp		<- df.ind[, which( is.na(DIAG_T) & T1_SEQ>min(DIAG_T, na.rm=1) )]
		set( df.ind, tmp, 'T1_SEQ', NA_real_)
		tmp	<- df.ind[, which(T1_SEQ>ART1_T)]
		set(df.ind, tmp, 'T1_SEQ', df.ind[tmp,ART1_T])		
	}
	if(seqtime.mode=='AtTrm')
	{
		df.ind[, T1_SEQ:= df.ind[, TIME_TR+0.1]]	
		tmp		<- df.ind[, which(is.na(T1_SEQ) & DOD<T1_SEQ)]
		set( df.ind, tmp, 'T1_SEQ', df.ind[tmp, TIME_TR+(TIME_TR+DOD)/2] )
		tmp	<- df.ind[, which(T1_SEQ>ART1_T)]
		set(df.ind, tmp, 'T1_SEQ', df.ind[tmp,ART1_T])		
	}	
	if(seqtime.mode=='AtDiag')
	{
		df.ind[, T1_SEQ:= df.ind[, DIAG_T]]	
		tmp		<- df.ind[, which(is.na(T1_SEQ) & DOD-TIME_TR>0.5)]
		set( df.ind, tmp, 'T1_SEQ', df.ind[tmp, runif(length(tmp), TIME_TR+0.5, DOD)] )
		tmp		<- df.ind[, which( is.na(DIAG_T) & T1_SEQ>min(DIAG_T, na.rm=1) )]
		set( df.ind, tmp, 'T1_SEQ', NA_real_)
		tmp	<- df.ind[, which(T1_SEQ>ART1_T)]
		set(df.ind, tmp, 'T1_SEQ', df.ind[tmp,ART1_T])		
	}
	if(seqtime.mode=='AtART')
	{
		df.ind[, T1_SEQ:= df.ind[, ART1_T]]		
		tmp		<- df.ind[, which(is.na(T1_SEQ) & DOD-TIME_TR>0.5)]
		set( df.ind, tmp, 'T1_SEQ', df.ind[tmp, runif(length(tmp), TIME_TR+0.5, DOD)] )
		tmp		<- df.ind[, which( is.na(DIAG_T) & T1_SEQ>min(DIAG_T, na.rm=1) )]
		set( df.ind, tmp, 'T1_SEQ', NA_real_)
		tmp	<- df.ind[, which(T1_SEQ>ART1_T)]
		set(df.ind, tmp, 'T1_SEQ', df.ind[tmp,ART1_T])		
	}
	df.ind
}
##--------------------------------------------------------------------------------------------------------
#	sample proportional to untreated
#
#	olli originally written 27-01-2015
##--------------------------------------------------------------------------------------------------------
PANGEA.Seqsampler.sample.prop.to.untreated<- function(df.ind, df.epi, pipeline.args, sp, verbose=1)
{	
	s.total					<- round( df.epi[nrow(df.epi), PREV_EVER] * sp )
	if(verbose)
	{
		cat(paste('\nSample proportional to untreated population'))
		cat(paste('\nSampling sequences, target is n=', s.total))		
	}
	#	setup number of sequences to be sampled for each year
	df.sample	<- subset( df.epi, YR>= pipeline.args['yr.start',][, as.numeric(v)] & YR<pipeline.args['yr.end',][, as.numeric(v)] )	
	set(df.sample, NULL, 's.nTOTAL', as.numeric(rmultinom(1, s.total, df.sample[, PREV-TREATED]/ df.epi[, sum(PREV-TREATED)])) )
	if(verbose)
		cat(paste('\nSampling sequences, scheduled number is n=', df.sample[, sum(s.nTOTAL)]))
	stopifnot(df.sample[, all(s.nTOTAL>=0)])
	if(verbose)
		cat(paste('\n prop of sequences sampled among HIV+=', df.sample[, sum( s.nTOTAL )] / df.sample[, rev(PREV_EVER)[1]]))
	#
	#	SAMPLE INFECTED INDIVIDUALS BASED ON NUMBERS PER YEAR
	#
	#	sample non-incident cases by year
	df.inds		<- copy(df.ind)
	df.inds[, TIME_SEQ:= NA_real_]
	for(yr in df.sample[, YR])		#TODO took out [-1] because there are s.n.notINC for DSPS in 1980
	{
		#	of all infected and not incident and not yet sampled, sample
		if(verbose)
			cat(paste('\nadd samples in year',yr,', required=', subset( df.sample, YR==yr )[, s.nTOTAL]))
		tmp		<- df.inds[, which(is.na(TIME_SEQ) & floor(DOB)<=yr & ceiling(DOD)>yr+1 & yr-floor(TIME_TR)>=0)]
		if(verbose)
			cat(paste('\navailable non-sampled HIV+ individuals in year=', length(tmp)))
		tmp		<- df.inds[, which(is.na(TIME_SEQ) & floor(DOB)<=yr & ceiling(DOD)>yr+1 & yr-floor(T1_SEQ)==0)]
		if(verbose)
			cat(paste('\navailable non-sampled HIV+ individuals with suggested sampling date in year=', length(tmp)))
		stopifnot(length(tmp)>0)
		tmp		<- sample(tmp, subset( df.sample, YR==yr )[, s.nTOTAL])
		set( df.inds, tmp, 'TIME_SEQ', df.inds[tmp, T1_SEQ] )				
	}
	stopifnot( df.inds[, !any(TIME_SEQ>DOD, na.rm=TRUE)] )
	stopifnot( df.inds[, !any(TIME_SEQ<TIME_TR, na.rm=TRUE)] )
	if(verbose)
	{
		cat(paste('\n total number of HIV+ in df.inds=', nrow(subset(df.inds, !is.na(TIME_TR)))))
		cat(paste('\n total number of sampled HIV+ in df.inds=', nrow(subset(df.inds, !is.na(TIME_TR) & !is.na(TIME_SEQ)))))		
	}
	list(df.inds=df.inds, df.sample=df.sample)
}
##--------------------------------------------------------------------------------------------------------
#	sample proportional to infected
#	olli originally written 03-04-2015
##--------------------------------------------------------------------------------------------------------
PANGEA.Seqsampler.sample.prop.to.T1SEQ<- function(df.ind, df.epi, pipeline.args, verbose=0)
{	
	sp						<- pipeline.args['s.PREV.max',][, as.numeric(v)]
	if(verbose)
	{
		cat(paste('\nSample proportional to T1SEQ'))
		cat(paste('\nSampling fraction per year is n=', sp))		
	}
	#	setup number of sequences to be sampled for each year
	df.sample	<- subset( df.epi, YR>= pipeline.args['yr.start',][, as.numeric(v)] & YR<pipeline.args['yr.end',][, as.numeric(v)] )
	set(df.sample, NULL, 's.nTOTAL', 0)	
	set(df.sample, NULL, 's.nTOTAL', rbinom(nrow(df.sample), df.sample[, T1_SEQ], sp) )
	if(verbose)
		cat(paste('\nSampling sequences, scheduled number is n=', df.sample[, sum(s.nTOTAL)]))
	stopifnot(df.sample[, all(s.nTOTAL>=0)])
	#
	#	SAMPLE INFECTED INDIVIDUALS BASED ON NUMBERS PER YEAR
	#
	#	sample non-incident cases by year
	df.inds		<- copy(df.ind)
	df.inds[, TIME_SEQ:= NA_real_]
	for(yr in df.sample[, YR])		#TODO took out [-1] because there are s.n.notINC for DSPS in 1980
	{
		#	of all infected and not incident and not yet sampled, sample
		if(verbose)
			cat(paste('\nadd samples in year',yr,', required=', subset( df.sample, YR==yr )[, s.nTOTAL]))
		tmp		<- df.inds[, which(is.na(TIME_SEQ) & floor(T1_SEQ)==yr & TIME_TR>=pipeline.args['yr.start',][, as.numeric(v)]) ]
		if(verbose)
			cat(paste('\navailable non-sampled HIV+ individuals in year=', length(tmp)))
		stopifnot(length(tmp)>0)
		tmp		<- sample(tmp, subset( df.sample, YR==yr )[, s.nTOTAL])
		if(verbose)
			cat(paste('\nfound samples in year',yr,', required=', length(tmp)))		
		set( df.inds, tmp, 'TIME_SEQ', df.inds[tmp, T1_SEQ] )				
	}
	stopifnot( df.inds[, !any(TIME_SEQ>DOD, na.rm=TRUE)] )
	stopifnot( df.inds[, !any(TIME_SEQ<TIME_TR, na.rm=TRUE)] )
	if(verbose)
	{
		cat(paste('\n total number of HIV+ after start in df.inds=', nrow(subset(df.inds, IDPOP>=0 & !is.na(TIME_TR) & floor(TIME_TR)>=pipeline.args['yr.start',][, as.numeric(v)] ))))
		cat(paste('\n total number of sampled HIV+ after start in df.inds=', nrow(subset(df.inds, IDPOP>=0 & !is.na(TIME_TR) & floor(TIME_TR)>=pipeline.args['yr.start',][, as.numeric(v)] & !is.na(TIME_SEQ)))))		
		cat(paste('\n total number of non-sampled HIV+ after start in df.inds=', nrow(subset(df.inds, !is.na(TIME_TR) & is.na(TIME_SEQ) & IDPOP>=0 & floor(TIME_TR)>=pipeline.args['yr.start',][, as.numeric(v)] ))))
		cat(paste('\n total number of non-sampled HIV+ after start with TIME_TR before end in df.inds=', nrow(subset(df.inds, !is.na(TIME_TR) & is.na(TIME_SEQ) & IDPOP>=0 & TIME_TR<pipeline.args['yr.end',][, as.numeric(v)] & floor(TIME_TR)>=pipeline.args['yr.start',][, as.numeric(v)] ))))
		cat(paste('\n total number of non-sampled HIV+ after start with T1_SEQ before end in df.inds=', nrow(subset(df.inds, !is.na(TIME_TR) & is.na(TIME_SEQ) & IDPOP>=0 & T1_SEQ<pipeline.args['yr.end',][, as.numeric(v)] & floor(TIME_TR)>=pipeline.args['yr.start',][, as.numeric(v)] ))))
	}
	list(df.inds=df.inds, df.sample=df.sample)
}
##--------------------------------------------------------------------------------------------------------
#	sample proportional to diagnoses before and after interventions
#	s% of those newly diagnosed per year until 2015
#	2*s% of those newly diagnosed per year after 2015
#	from before then: 50 (uniform)	
#
#	olli originally written 27-01-2015
##--------------------------------------------------------------------------------------------------------
PANGEA.Seqsampler.sample.prop.to.diagnosis<- function(df.ind, df.epi, pipeline.args, sp)
{
	s.total					<- round( df.epi[nrow(df.epi), PREV] * sp )
	s.archival.yr			<- subset(df.epi, DIAG==0)[, tail(YR,1)]
	s.diagb4intervention.n	<- subset( df.epi, YR<pipeline.args['s.INTERVENTION.start',][, as.numeric(v)] )[, tail(DIAG,1)]
	s.diagb4intervention	<- (s.total-pipeline.args['s.ARCHIVAL.n',][, as.numeric(v)]) / (s.diagb4intervention.n+pipeline.args['s.INTERVENTION.mul'][, as.numeric(v)]*(df.epi[nrow(df.epi), DIAG]-s.diagb4intervention.n))
	cat(paste('\nSample proportional to diagnoses up to intervention, and proportional to diagnoses after intervention start'))
	cat(paste('\nSampling sequences, target is n=', s.total))
	cat(paste('\nNo diagnoses up to yr, sampling archival sequences till then. yr=', s.archival.yr))
	cat(paste('\nSampling archival sequences before any diagnoses, n=', pipeline.args['s.ARCHIVAL.n',][, as.numeric(v)]))	
	cat(paste('\nSampling sequences before intervention start from diagnosed, p=', s.diagb4intervention))
	cat(paste('\nSampling sequences after intervention start from diagnosed, p=', 2*s.diagb4intervention))
	#	setup number of sequences to be sampled for each year
	df.sample	<- subset( df.epi, YR>= pipeline.args['yr.start',][, as.numeric(v)] & YR<pipeline.args['yr.end',][, as.numeric(v)] )
	set(df.sample, NULL, 's.nTOTAL', 0)	
	tmp			<- df.sample[, which(YR<=s.archival.yr)]
	set(df.sample, tmp, 's.nTOTAL', as.numeric(rmultinom(1, pipeline.args['s.ARCHIVAL.n',][, as.numeric(v)], rep(1/length(tmp), length(tmp)))))
	tmp			<- df.sample[, which( YR>s.archival.yr & YR<pipeline.args['s.INTERVENTION.start',][, as.numeric(v)] )]	
	set(df.sample, tmp, 's.nTOTAL', as.numeric(rmultinom(1, round( s.diagb4intervention.n*s.diagb4intervention ), df.sample[tmp, NEW_DIAG]/df.sample[tmp, sum(NEW_DIAG)])) )
	tmp			<- df.sample[, which(YR>=pipeline.args['s.INTERVENTION.start',][, as.numeric(v)]) ]
	set(df.sample, tmp, 's.nTOTAL', as.numeric(rmultinom(1, round( (df.sample[nrow(df.sample), DIAG]-s.diagb4intervention.n)*pipeline.args['s.INTERVENTION.mul'][, as.numeric(v)]*s.diagb4intervention ), df.sample[tmp, NEW_DIAG]/df.sample[tmp, sum(NEW_DIAG)])) )
	cat(paste('\nSampling sequences, scheduled number is n=', df.sample[, sum(s.nTOTAL)]))
	stopifnot(df.sample[, all(s.nTOTAL>=0)])
	cat(paste('\n prop of sequences sampled among HIV+=', df.sample[, sum( s.nTOTAL )] / df.sample[, rev(PREV)[1]]))
	stopifnot( abs(sp-df.sample[, sum( s.nTOTAL )] / df.sample[, rev(PREV)[1]])<=sp*0.1 )
	#
	#	SAMPLE INFECTED INDIVIDUALS BASED ON NUMBERS PER YEAR
	#
	#	sample non-incident cases by year
	df.inds		<- copy(df.ind)
	df.inds[, TIME_SEQ:= NA_real_]
	for(yr in df.sample[, YR])		#TODO took out [-1] because there are s.n.notINC for DSPS in 1980
	{
		#	of all infected and not incident and not yet sampled, sample
		cat(paste('\nadd samples in year',yr,', required=', subset( df.sample, YR==yr )[, s.nTOTAL]))
		tmp		<- df.inds[, which(is.na(TIME_SEQ) & floor(DOB)<=yr & ceiling(DOD)>yr+1 & yr-floor(TIME_TR)>=0)]
		cat(paste('\navailable non-sampled HIV+ individuals in year=', length(tmp)))
		tmp		<- df.inds[, which(is.na(TIME_SEQ) & floor(DOB)<=yr & ceiling(DOD)>yr+1 & yr-floor(T1_SEQ)==0)]
		cat(paste('\navailable non-sampled HIV+ individuals with suggested sampling date in year=', length(tmp)))
		stopifnot(length(tmp)>0)
		tmp		<- sample(tmp, subset( df.sample, YR==yr )[, s.nTOTAL])
		set( df.inds, tmp, 'TIME_SEQ', df.inds[tmp, T1_SEQ] )				
	}
	stopifnot( df.inds[, !any(TIME_SEQ>DOD, na.rm=TRUE)] )
	stopifnot( df.inds[, !any(TIME_SEQ<TIME_TR, na.rm=TRUE)] )	
	cat(paste('\n total number of HIV+ in df.inds=', nrow(subset(df.inds, !is.na(TIME_TR)))))
	cat(paste('\n total number of sampled HIV+ in df.inds=', nrow(subset(df.inds, !is.na(TIME_TR) & !is.na(TIME_SEQ)))))
	list(df.inds=df.inds, df.sample=df.sample)
}
##--------------------------------------------------------------------------------------------------------
#	sample a fixed proportion after the start of the intervention (2015)
#	after 2015, a fixed amount per year
#	olli originally written 27-01-2015
##--------------------------------------------------------------------------------------------------------
PANGEA.Seqsampler.sample.fixed.to.prop<- function(df.ind, df.epi, pipeline.args, sp, verbose=1)
{
	# 	SAMPLING PROBABILITIES and TOTALS PER YEAR
	s.total					<- round( df.epi[nrow(df.epi), PREV_EVER] * sp )
	s.archival.yr			<- subset(df.epi, DIAG==0)[, tail(YR,1)]
	s.b4intervention.n		<- round(s.total*(1-pipeline.args['s.INTERVENTION.prop',][, as.numeric(v)]))
	stopifnot(s.b4intervention.n>=pipeline.args['s.ARCHIVAL.n',][, as.numeric(v)])
	if(verbose)
	{
		cat(paste('\nSample a fixed proportion after the start of the intervention (2015), and of those a fixed amount per year'))
		cat(paste('\nSampling sequences, target is n=', s.total))
		cat(paste('\nNo diagnoses up to yr, sampling archival sequences till then. yr=', s.archival.yr))
		cat(paste('\nSampling archival sequences before any diagnoses, n=', pipeline.args['s.ARCHIVAL.n',][, as.numeric(v)]))
		cat(paste('\nSampling sequences before intervention start from diagnosed, n=', s.b4intervention.n))
		cat(paste('\nSampling sequences after intervention start from diagnosed, n=', s.total	- s.b4intervention.n))		
	}
	#	setup number of sequences to be sampled for each year
	df.sample	<- subset( df.epi, YR>= pipeline.args['yr.start',][, as.numeric(v)] & YR<pipeline.args['yr.end',][, as.numeric(v)] )
	set(df.sample, NULL, 's.nTOTAL', 0)	
	tmp			<- df.sample[, which(YR<=s.archival.yr)]
	set(df.sample, tmp, 's.nTOTAL', as.numeric(rmultinom(1, pipeline.args['s.ARCHIVAL.n',][, as.numeric(v)], rep(1/length(tmp), length(tmp)))))
	if(verbose)
		cat(paste('\nSampling sequences from archival, scheduled number is n=', df.sample[tmp, sum(s.nTOTAL)]))
	tmp			<- df.sample[, which( YR>s.archival.yr & YR<pipeline.args['s.INTERVENTION.start',][, as.numeric(v)] )]	
	set(df.sample, tmp, 's.nTOTAL', as.numeric(rmultinom(1, s.b4intervention.n-pipeline.args['s.ARCHIVAL.n',][, as.numeric(v)], df.sample[tmp, NEW_DIAG]/df.sample[tmp, sum(NEW_DIAG)])) )
	if(verbose)
		cat(paste('\nSampling sequences before intervention, scheduled number is n=', df.sample[tmp, sum(s.nTOTAL)]))
	tmp			<- df.sample[, which(YR>=pipeline.args['s.INTERVENTION.start',][, as.numeric(v)]) ]
	set(df.sample, tmp, 's.nTOTAL', as.numeric(rmultinom(1, s.total	- s.b4intervention.n, rep(1/length(tmp),length(tmp))) ) )
	if(verbose)
		cat(paste('\nSampling sequences after intervention, scheduled number is n=', df.sample[tmp, sum(s.nTOTAL)]))
	stopifnot(df.sample[, all(s.nTOTAL>=0)])
	#stopifnot( abs(sp-df.sample[, sum( s.nTOTAL )] / df.sample[, rev(PREV)[1]])<=sp*0.1 )
	#
	#	SAMPLE INFECTED INDIVIDUALS BASED ON NUMBERS PER YEAR
	# 
	#	sample non-incident cases by year
	df.inds		<- copy(df.ind)
	df.inds[, TIME_SEQ:= NA_real_]
	for(yr in df.sample[, YR])		#TODO took out [-1] because there are s.n.notINC for DSPS in 1980
	{
		#	of all infected and not incident and not yet sampled, sample
		if(verbose)
			cat(paste('\nadd samples in year',yr,', required=', subset( df.sample, YR==yr )[, s.nTOTAL]))
		tmp		<- df.inds[, which(is.na(TIME_SEQ) & floor(DOB)<=yr & ceiling(DOD)>yr+1 & yr-floor(TIME_TR)>=0)]
		if(verbose)
			cat(paste('\navailable non-sampled HIV+ individuals in year=', length(tmp)))
		tmp		<- df.inds[, which(is.na(TIME_SEQ) & floor(DOB)<=yr & ceiling(DOD)>yr+1 & yr-floor(T1_SEQ)==0)]
		if(verbose)
			cat(paste('\navailable non-sampled HIV+ individuals with suggested sampling date in year=', length(tmp)))				
		k		<- 1
		while(length(tmp)<subset( df.sample, YR==yr )[, s.nTOTAL])
		{	
			if(verbose)
				cat(paste('\nCannot find samples, fall back to ',k,'th previous year; n=', subset( df.sample, YR==yr )[, s.nTOTAL]-length(tmp) ))
			tmp2	<- df.inds[, which(is.na(TIME_SEQ) & floor(DOB)<=yr & ceiling(DOD)>yr+1 & yr-floor(T1_SEQ)==k & (is.na(ART1_T) | ART1_T-yr<=1)) ]
			if(verbose)
				cat(paste('\navailable non-sampled HIV+ individuals with suggested sampling date one year before=', length(tmp2)))
			tmp2	<- sample(tmp2,  min(subset( df.sample, YR==yr )[, s.nTOTAL]-length(tmp), length(tmp2)) )
			tmp		<- c(tmp, tmp2)
			k		<- k+1
			stopifnot(k<=10)
		}	
		stopifnot(length(tmp)>=subset( df.sample, YR==yr )[, s.nTOTAL])
		tmp		<- sample(tmp, subset( df.sample, YR==yr )[, s.nTOTAL])
		set( df.inds, tmp, 'TIME_SEQ', df.inds[tmp, T1_SEQ+yr-floor(T1_SEQ)] )			
	}
	stopifnot( df.inds[, !any(TIME_SEQ>DOD, na.rm=TRUE)] )
	stopifnot( df.inds[, !any(TIME_SEQ<TIME_TR, na.rm=TRUE)] )
	if(verbose)
	{
		cat(paste('\n total number of HIV+ in df.inds=', nrow(subset(df.inds, !is.na(TIME_TR)))))
		cat(paste('\n total number of sampled HIV+ in df.inds=', nrow(subset(df.inds, !is.na(TIME_TR) & !is.na(TIME_SEQ)))))
	}
	list(df.inds=df.inds, df.sample=df.sample)
}
##--------------------------------------------------------------------------------------------------------
#	sample proportional to diagnoses before and after interventions
#	s% of those newly diagnosed per year until 2015
#	from before any diagnosed: 50 (uniform)
#	after 2015, a fixed amount per year
#	olli originally written 27-01-2015
##--------------------------------------------------------------------------------------------------------
PANGEA.Seqsampler.sample.prop.to.diagnosis.b4intervention<- function(df.ind, df.epi, pipeline.args, sp, verbose=1)
{
	# 	SAMPLING PROBABILITIES and TOTALS PER YEAR
	s.total					<- round( df.epi[nrow(df.epi), PREV_EVER] * sp )
	s.archival.yr			<- subset(df.epi, DIAG==0)[, tail(YR,1)]
	s.diagb4intervention.n	<- subset( df.epi, YR<pipeline.args['s.INTERVENTION.start',][, as.numeric(v)] )[, tail(DIAG,1)]
	s.diagb4intervention	<- (s.total-pipeline.args['s.ARCHIVAL.n',][, as.numeric(v)]) / (s.diagb4intervention.n+pipeline.args['s.INTERVENTION.mul'][, as.numeric(v)]*(df.epi[nrow(df.epi), DIAG]-s.diagb4intervention.n))
	if(verbose)
	{
		cat(paste('\nSample proportional to diagnoses up to intervention, and then a fixed amount per year'))
		cat(paste('\nSampling sequences, target is n=', s.total))
		cat(paste('\nNo diagnoses up to yr, sampling archival sequences till then. yr=', s.archival.yr))
		cat(paste('\nSampling archival sequences before any diagnoses, n=', pipeline.args['s.ARCHIVAL.n',][, as.numeric(v)]))
		cat(paste('\nNumber of diagnoses before intervention, n=', s.diagb4intervention.n))
		cat(paste('\nNumber of diagnoses after intervention start, n=', df.epi[nrow(df.epi), DIAG]-s.diagb4intervention.n))
		cat(paste('\nSampling sequences before intervention start from diagnosed, p=', s.diagb4intervention))
		cat(paste('\nSampling sequences after intervention start from diagnosed, n=', round(s.diagb4intervention*pipeline.args['s.INTERVENTION.mul'][, as.numeric(v)]*(df.epi[nrow(df.epi), DIAG]-s.diagb4intervention.n))))		
	}
	#	setup number of sequences to be sampled for each year
	df.sample	<- subset( df.epi, YR>= pipeline.args['yr.start',][, as.numeric(v)] & YR<pipeline.args['yr.end',][, as.numeric(v)] )
	set(df.sample, NULL, 's.nTOTAL', 0)	
	tmp			<- df.sample[, which(YR<=s.archival.yr)]
	set(df.sample, tmp, 's.nTOTAL', as.numeric(rmultinom(1, pipeline.args['s.ARCHIVAL.n',][, as.numeric(v)], rep(1/length(tmp), length(tmp)))))
	if(verbose)
		cat(paste('\nSampling sequences from archival, scheduled number is n=', df.sample[tmp, sum(s.nTOTAL)]))
	tmp			<- df.sample[, which( YR>s.archival.yr & YR<pipeline.args['s.INTERVENTION.start',][, as.numeric(v)] )]	
	set(df.sample, tmp, 's.nTOTAL', as.numeric(rmultinom(1, round( s.diagb4intervention.n*s.diagb4intervention ), df.sample[tmp, NEW_DIAG]/df.sample[tmp, sum(NEW_DIAG)])) )
	if(verbose)
		cat(paste('\nSampling sequences before intervention, scheduled number is n=', df.sample[tmp, sum(s.nTOTAL)]))
	tmp			<- df.sample[, which(YR>=pipeline.args['s.INTERVENTION.start',][, as.numeric(v)]) ]
	set(df.sample, tmp, 's.nTOTAL', as.numeric(rmultinom(1, round( (df.sample[nrow(df.sample), DIAG]-s.diagb4intervention.n)*pipeline.args['s.INTERVENTION.mul'][, as.numeric(v)]*s.diagb4intervention ), rep(1/length(tmp),length(tmp))) ) )
	if(verbose)
		cat(paste('\nSampling sequences after intervention, scheduled number is n=', df.sample[tmp, sum(s.nTOTAL)]))
	stopifnot(df.sample[, all(s.nTOTAL>=0)])
	#stopifnot( abs(sp-df.sample[, sum( s.nTOTAL )] / df.sample[, rev(PREV)[1]])<=sp*0.1 )
	#
	#	SAMPLE INFECTED INDIVIDUALS BASED ON NUMBERS PER YEAR
	# 
	#	sample non-incident cases by year
	df.inds		<- copy(df.ind)
	df.inds[, TIME_SEQ:= NA_real_]
	for(yr in df.sample[, YR])		#TODO took out [-1] because there are s.n.notINC for DSPS in 1980
	{
		#	of all infected and not incident and not yet sampled, sample
		if(verbose)
			cat(paste('\nadd samples in year',yr,', required=', subset( df.sample, YR==yr )[, s.nTOTAL]))
		tmp		<- df.inds[, which(is.na(TIME_SEQ) & floor(DOB)<=yr & ceiling(DOD)>yr+1 & yr-floor(TIME_TR)>=0)]
		if(verbose)
			cat(paste('\navailable non-sampled HIV+ individuals in year=', length(tmp)))
		tmp		<- df.inds[, which(is.na(TIME_SEQ) & floor(DOB)<=yr & ceiling(DOD)>yr+1 & yr-floor(T1_SEQ)==0)]
		if(verbose)
			cat(paste('\navailable non-sampled HIV+ individuals with suggested sampling date in year=', length(tmp)))				
		k		<- 1
		while(length(tmp)<subset( df.sample, YR==yr )[, s.nTOTAL])
		{	
			if(verbose)
				cat(paste('\nCannot find samples, fall back to ',k,'th previous year; n=', subset( df.sample, YR==yr )[, s.nTOTAL]-length(tmp) ))
			tmp2	<- df.inds[, which(is.na(TIME_SEQ) & floor(DOB)<=yr & ceiling(DOD)>yr+1 & yr-floor(T1_SEQ)==k & (is.na(ART1_T) | ART1_T-yr<=1)) ]
			if(verbose)
				cat(paste('\navailable non-sampled HIV+ individuals with suggested sampling date one year before=', length(tmp2)))
			tmp2	<- sample(tmp2,  min(subset( df.sample, YR==yr )[, s.nTOTAL]-length(tmp), length(tmp2)) )
			tmp		<- c(tmp, tmp2)
			k		<- k+1
			stopifnot(k<=10)
		}	
		stopifnot(length(tmp)>=subset( df.sample, YR==yr )[, s.nTOTAL])
		tmp		<- sample(tmp, subset( df.sample, YR==yr )[, s.nTOTAL])
		set( df.inds, tmp, 'TIME_SEQ', df.inds[tmp, T1_SEQ+yr-floor(T1_SEQ)] )			
	}
	stopifnot( df.inds[, !any(TIME_SEQ>DOD, na.rm=TRUE)] )
	stopifnot( df.inds[, !any(TIME_SEQ<TIME_TR, na.rm=TRUE)] )
	if(verbose)
	{
		cat(paste('\n total number of HIV+ in df.inds=', nrow(subset(df.inds, !is.na(TIME_TR)))))
		cat(paste('\n total number of sampled HIV+ in df.inds=', nrow(subset(df.inds, !is.na(TIME_TR) & !is.na(TIME_SEQ)))))
	}
	list(df.inds=df.inds, df.sample=df.sample)
}
##--------------------------------------------------------------------------------------------------------
#	return ancestral sequence sampler	
#	olli originally written 27-01-2015
##--------------------------------------------------------------------------------------------------------
PANGEA.Seqsampler.v4<- function(df.ind, df.trm, pipeline.args, outfile.ind, outfile.trm, with.plot=1)
{	
	stopifnot(grepl('Fixed2Prop|Prop2DiagB4I|Prop2Diag|Prop2Untreated|Prop2SuggestedSampling',pipeline.args['s.MODEL',][, v]))
	# 	compute prevalence and incidence by year	
	epi.adult	<- 13
	suppressWarnings( df.trm[, YR:= df.trm[, floor(TIME_TR)]] )
	df.epi		<- df.trm[, list(INC=length(IDREC), INC_ACUTE=length(which(TR_ACUTE=='Yes')),IMPORT=length(which(IDTR<0))), by='YR']
	tmp			<- df.epi[, 	{
				sexactive		<- which( floor(df.ind[['DOB']]+epi.adult)<=YR  &  ceiling(df.ind[['DOD']])>YR )
				infected.ever	<- which( floor(df.ind[['TIME_TR']])<=YR )
				infected		<- which( floor(df.ind[['DOB']])<=YR  &  floor(df.ind[['DOD']])>YR  &  floor(df.ind[['TIME_TR']])<=YR )
				diag			<- which( floor(df.ind[['DOB']])<=YR  &  floor(df.ind[['DOD']])>YR  &  floor(df.ind[['DIAG_T']])<=YR )
				diag.new		<- which( floor(df.ind[['DIAG_T']])==YR )
				treated			<- which( floor(df.ind[['DOB']])<=YR  &  floor(df.ind[['DOD']])>YR  &  floor(df.ind[['ART1_T']])<=YR & (is.na(df.ind[['VLS1_TE']]) | floor(df.ind[['VLS1_TE']])>YR) )
				infdead			<- which( floor(df.ind[['DOD']])==YR  &  floor(df.ind[['TIME_TR']])<=YR )
				suggsampling	<- which( floor(df.ind[['T1_SEQ']])==YR & floor(df.ind[['TIME_TR']])>=pipeline.args['yr.start',][, as.numeric(v)])
				list(POP=length(sexactive), PREV=length(infected), PREV_EVER=length(infected.ever), PREVDIED=length(infdead), DIAG=length(diag), NEW_DIAG=length(diag.new), TREATED=length(treated), T1_SEQ=length(suggsampling))				
			},by='YR']
	df.epi		<- merge( tmp, df.epi, by='YR' )		
	set(df.epi, NULL, 'PREVp', df.epi[, PREV/(POP-PREV)])	
	set(df.epi, NULL, 'INCp', df.epi[, INC/(POP-PREV)])
	set(df.epi, NULL, 'IMPORTp', df.epi[, IMPORT/INC])
	set(df.epi, NULL, 'ACUTEp', df.epi[, INC_ACUTE/INC])
	set(df.epi, NULL, 'ARTcov', df.epi[, TREATED/PREV])
	set(df.epi, NULL, 'UNDIAGp', df.epi[, (PREV-DIAG)/PREV])
	set(df.epi, NULL, 'GROWTHr', c(NA_real_, df.epi[, diff(log(PREV))]))
	
	stopifnot( !is.na(pipeline.args['s.PREV.max.n',][, as.numeric(v)]) | !is.na(pipeline.args['s.PREV.max',][, as.numeric(v)]) )
	s.PREV.max.guess	<- pipeline.args['s.PREV.max.n',][, as.numeric(v)] / df.epi[nrow(df.epi), PREV_EVER]	
	if(is.na(s.PREV.max.guess) & pipeline.args['s.MODEL',][, v]=='Prop2DiagB4I')
	{
		#	calibration runs to determine sampling coverage
		#	set(pipeline.args, pipeline.args[, which(stat=='s.PREV.max')], 'v', '0.08')
		s.PREV.max.guess	<- seq( pipeline.args['s.PREV.max',][, as.numeric(v)]*0.4, pipeline.args['s.PREV.max',][, as.numeric(v)]*0.6, length.out=40) 
		tmp					<- sapply(s.PREV.max.guess, function(x)
				{					
					cat(paste('\ntry s.PREV.max.guess=',x))
					tmp					<- PANGEA.Seqsampler.sample.prop.to.diagnosis.b4intervention(df.ind, df.epi, pipeline.args, x, verbose=0)
					
					tmp2				<- subset( tmp$df.inds, DOD>pipeline.args['s.INTERVENTION.start',][, as.numeric(v)] & floor(TIME_TR)<pipeline.args['s.INTERVENTION.start',][, as.numeric(v)])
					sc.alive15.infl15	<- tmp2[, length(which(!is.na(TIME_SEQ)))/length(TIME_SEQ) ]					
					tmp2				<- subset( tmp$df.inds, DOD>pipeline.args['s.INTERVENTION.start',][, as.numeric(v)] & floor(TIME_TR)<pipeline.args['yr.end',][, as.numeric(v)])
					sc.alive15.infl20	<- tmp2[, length(which(!is.na(TIME_SEQ)))/length(TIME_SEQ) ]
					tmp2				<- subset( tmp$df.inds, DOD>pipeline.args['yr.end',][, as.numeric(v)] & floor(TIME_TR)<pipeline.args['yr.end',][, as.numeric(v)])
					sc.alive20.infl20	<- tmp2[, length(which(!is.na(TIME_SEQ)))/length(TIME_SEQ) ]					
					tmp2				<- subset( tmp$df.inds, floor(TIME_TR)>=2000 & floor(TIME_TR)<pipeline.args['yr.end',][, as.numeric(v)])
					sc.infg99l20		<- tmp2[, length(which(!is.na(TIME_SEQ)))/length(TIME_SEQ) ]
					tmp2				<- subset( tmp$df.inds, floor(TIME_TR)>=pipeline.args['s.INTERVENTION.start',][, as.numeric(v)] & floor(TIME_TR)<pipeline.args['yr.end',][, as.numeric(v)])
					sc.infg14l20		<- tmp2[, length(which(!is.na(TIME_SEQ)))/length(TIME_SEQ) ]
					sn.g14				<- subset( tmp$df.inds, TIME_SEQ>=pipeline.args['s.INTERVENTION.start',][, as.numeric(v)] & floor(TIME_TR)<pipeline.args['yr.end',][, as.numeric(v)])[, length(which(!is.na(TIME_SEQ))) ] 
					sn.total			<- subset(tmp$df.inds, floor(TIME_TR)<pipeline.args['yr.end',][, as.numeric(v)])[, length(which(!is.na(TIME_SEQ))) ]
					sp.g14				<- sn.g14 / sn.total
					c(sc.alive15.infl15=sc.alive15.infl15, sc.alive15.infl20=sc.alive15.infl20, sc.alive20.infl20=sc.alive20.infl20, sc.infg99l20=sc.infg99l20, sc.infg14l20=sc.infg14l20, sn.g14=sn.g14, sn.total=sn.total, sp.g14=sp.g14)					
				})
		tmp					<- as.data.table(t(tmp))
		tmp[, s.PREV.max.guess:=s.PREV.max.guess]
		#	use best guess
		s.PREV.max.guess	<- tmp[which.min(abs(sc.alive20.infl20-pipeline.args['s.PREV.max',][, as.numeric(v)])), ][, s.PREV.max.guess]		
	}	
	if(!is.na(s.PREV.max.guess))
		cat(paste('\nFound best s.PREV.max.guess=',s.PREV.max.guess))
	if(pipeline.args['s.MODEL',][, v]=='Prop2SuggestedSampling')
		tmp				<- PANGEA.Seqsampler.sample.prop.to.T1SEQ(df.ind, df.epi, pipeline.args, verbose=1)	
	if(pipeline.args['s.MODEL',][, v]=='Fixed2Prop')
		tmp				<- PANGEA.Seqsampler.sample.fixed.to.prop(df.ind, df.epi, pipeline.args, s.PREV.max.guess, verbose=1)
	if(pipeline.args['s.MODEL',][, v]=='Prop2DiagB4I')
		tmp				<- PANGEA.Seqsampler.sample.prop.to.diagnosis.b4intervention(df.ind, df.epi, pipeline.args, s.PREV.max.guess, verbose=1)
	if(pipeline.args['s.MODEL',][, v]=='Prop2Diag')
		tmp				<- PANGEA.Seqsampler.sample.prop.to.diagnosis(df.ind, df.epi, pipeline.args, s.PREV.max.guess, verbose=1)
	if(pipeline.args['s.MODEL',][, v]=='Prop2Untreated')
		tmp				<- PANGEA.Seqsampler.sample.prop.to.untreated(df.ind, df.epi, pipeline.args, s.PREV.max.guess, verbose=1)	
	df.inds				<- copy(tmp$df.inds)
	df.sample			<- copy(tmp$df.sample)
	#
	tmp2				<- subset( df.inds, DOD>pipeline.args['yr.end',][, as.numeric(v)] & floor(TIME_TR)<pipeline.args['yr.end',][, as.numeric(v)])
	sc.alive20.infl20	<- tmp2[, length(which(!is.na(TIME_SEQ)))/length(TIME_SEQ) ]					
	tmp2				<- subset( df.inds, floor(TIME_TR)>=2000 & floor(TIME_TR)<pipeline.args['yr.end',][, as.numeric(v)])
	sc.infg99l20		<- tmp2[, length(which(!is.na(TIME_SEQ)))/length(TIME_SEQ) ]
	sn.g14				<- subset( df.inds, TIME_SEQ>=pipeline.args['s.INTERVENTION.start',][, as.numeric(v)] & floor(TIME_TR)<pipeline.args['yr.end',][, as.numeric(v)])[, length(which(!is.na(TIME_SEQ))) ] 
	sn.total			<- subset(df.inds, floor(TIME_TR)<pipeline.args['yr.end',][, as.numeric(v)])[, length(which(!is.na(TIME_SEQ))) ]
	cat(paste('\ncoverage alive20.infl20=', sc.alive20.infl20))
	cat(paste('\ncoverage infg99l20=', sc.infg99l20))
	cat(paste('\nseqs sn.g14=', sn.g14))
	cat(paste('\nseqs all=', sn.total))
	
	#	set sampling in df.trm
	tmp		<- subset( df.inds, !is.na(TIME_SEQ), select=c(IDPOP, TIME_SEQ) )
	setnames(tmp, c('IDPOP','TIME_SEQ'), c('IDREC','SAMPLED_REC'))
	df.trms	<- merge(df.trm, tmp, by='IDREC', all.x=TRUE)
	setnames(tmp, c('IDREC','SAMPLED_REC'), c('IDTR','SAMPLED_TR'))
	df.trms	<- merge(df.trms, tmp, by='IDTR', all.x=TRUE)	
	#
	#	seroprevalence survey
	#
	df.sp	<- PANGEA.SeroPrevalenceSurvey(df.inds, epi.adult=epi.adult, s.INTERVENTION.start=pipeline.args['s.INTERVENTION.start', ][, as.numeric(v)], sp.prop.of.sexactive=pipeline.args['sp.prop.of.sexactive', ][, as.numeric(v)], sp.times=c(12, 8, 4, 0) )
	file	<- gsub('IND.csv','CROSS_SECTIONAL_SURVEY.csv', outfile.ind)
	cat(paste('\nwrite to file', file))
	write.csv(df.sp, file=file)
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
		tmp		<- names(df.sample)[ which(as.character(sapply(df.sample, class))!='numeric') ]
		for(x in tmp)
			set(df.sample, NULL, x, as.numeric(df.sample[[x]]))
		tmp		<- melt(df.sample, id.vars=c('YR'))
		set(tmp, NULL, 'variable', tmp[, factor(variable, levels=tmp[, unique(variable)], labels=tmp[, unique(variable)])])
		
		ggplot(tmp, aes(x=YR, y=value, group=variable)) + geom_step(with.guide=FALSE) + 
				facet_grid(variable~., scales='free_y')  + 
				theme_bw() + scale_x_continuous(name='year', breaks=seq(1980,pipeline.args['yr.end',][, as.numeric(v)],2)) +
				theme(strip.text=element_text(size=7))
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_Totals.pdf',sep='')
		cat(paste('\nPlotting to file',file))
		ggsave(file=file, w=16, h=16)	
		#	plot CD4 counts at time of diagnosis
		tmp	<- subset(df.inds, !is.na(DIAG_T), select=c(DIAG_T, DIAG_CD4))
		tmp[, YR:= floor(DIAG_T)]
		tmp	<- merge(tmp, tmp[, list(DIAG_CD4_ME=mean(DIAG_CD4)), by='YR'], by='YR')
		setkey(tmp, YR)
		ggplot(tmp) + geom_point(aes(x=DIAG_T, y=DIAG_CD4)) + geom_step(data=unique(tmp), aes(x=YR, y=DIAG_CD4_ME), col='red', size=1.5) + scale_y_continuous(breaks=seq(100,1000,100))		
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_CD4.pdf',sep='')
		cat(paste('\nPlotting to file',file))
		ggsave(file=file, w=8, h=8)					
		ggplot(subset(df.inds, !is.na(DIAG_T)), aes(x=floor(DIAG_T), fill=cut(DIAG_CD4, breaks=c(0,200,350,500,700,2000)))) + geom_bar(binwidth=1, position='fill') + 
				labs(fill='CD4 category', y='percent of diagnosed', x='year diagnosed') + scale_y_continuous(breaks=seq(0,1,0.1)) + scale_fill_brewer(palette='Set1')		
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_CD4PROP_AMONGDIAG.pdf',sep='')
		cat(paste('\nPlotting to file',file))
		ggsave(file=file, w=8, h=8)		
		ggplot(subset(df.inds, !is.na(DIAG_T) & !is.na(TIME_SEQ)), aes(x=floor(DIAG_T), fill=cut(DIAG_CD4, breaks=c(0,200,350,500,700,2000)))) + geom_bar(binwidth=1, position='fill') + 
				labs(fill='CD4 category', y='percent of sequenced', x='year diagnosed') + scale_y_continuous(breaks=seq(0,1,0.1)) + scale_fill_brewer(palette='Set1')		
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_CD4PROP_AMONGSEQ.pdf',sep='')
		cat(paste('\nPlotting to file',file))
		ggsave(file=file, w=8, h=8)		
		#	plot distribution between transmission time and sequencing time		
		ggplot(subset(df.inds, !is.na(TIME_SEQ)), aes(x=TIME_SEQ-TIME_TR)) + geom_histogram(binwidth=1) + 
				scale_x_continuous(name='time from transmission to sequence sampling\n(years)', breaks=seq(0,100,2))
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_Time2Seq.pdf',sep='')
		ggsave(file=file, w=8, h=8)
		#	plot distribution between transmission time and diagnosis
		ggplot(subset(df.inds, !is.na(DIAG_T)), aes(x=DIAG_T-TIME_TR)) + geom_histogram(binwidth=1) + 
				scale_x_continuous(name='time from transmission to diagnosis\n(years)', breaks=seq(0,100,2))
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_Time2Diag.pdf',sep='')
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
		if(0)
		{
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
		}
		#ggplot(df.trms, aes(x=IDTR, y=TIME_TR)) + geom_point()
		#ggplot(df.trms, aes(x=IDTR, y=IDCLU)) + geom_point()
	}
	#
	#	SAVE SAMPLED RECIPIENTS AND TRANSMISSIONS TO SAMPLED RECIPIENTS
	#
	#	save for us
	file		<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'SAVE.R',sep='')
	save(file=file, df.epi, df.trms, df.inds, df.sample, df.sp)
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
#	olli originally written 20-09-2014
##--------------------------------------------------------------------------------------------------------
#' @title Create starting sequence sampler
#' @description Returns a function and function arguments to draw ancestral sequences. 
#' @param root.ctime.grace	Sample a starting sequence with time that matches a query times +- this grace
#' @param sample.grace		Internal parameter to make sure the requested number of samples is obtained. Internally oversample by this multiplier to the sample size, and then check if sequences are unique.
#' @return list of the sampler \code{rANCSEQ} and its arguments \code{rANCSEQ.args}
PANGEA.RootSeq.create.sampler<- function(root.ctime.grace= 0.5, sample.grace= 3)
{	
	file			<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfgBwhRc-_140907_n390_AncSeq.R')
	cat(paste('\nLoading starting sequences from file', file))
	load(file)		#expect "anc.seq.gag"  "anc.seq.pol"  "anc.seq.env"  "anc.seq.info"
	setkey(anc.seq.info, CALENDAR_TIME)
	rANCSEQ.args	<- list(	root.ctime.grace=root.ctime.grace, sample.grace=sample.grace, anc.seq.info=anc.seq.info, anc.seq.gag=anc.seq.gag, anc.seq.pol=anc.seq.pol, anc.seq.env=anc.seq.env)	
	
	rANCSEQ<- function(root.ctime, rANCSEQ.args)
	{		
		tmp				<- lapply(seq_along(root.ctime), function(i)
				{
					tmp	<- subset(rANCSEQ.args$anc.seq.info, 	CALENDAR_TIME>root.ctime[i]-rANCSEQ.args$root.ctime.grace &  CALENDAR_TIME<=root.ctime[i]+rANCSEQ.args$root.ctime.grace 	)
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
#	return evolutionary rate sampler for transmission edges	
#	olli originally written 16-09-2014
##--------------------------------------------------------------------------------------------------------
#' @title Create sampler of Evolutionary Rates for transmission edges
#' @description Returns a function to draw evolutionary rates for transmission edges. Currently modelled with a log normal density.
#' @param er.shift		shift to mean of the log normal density
#' @return R function
PANGEA.TransmissionEdgeEvolutionaryRate.create.sampler<- function(er.meanlog, er.sdlog)
{
	if(0)
	{
		x		<- seq(0.0005, 0.003, 0.00001)		
		tmp		<- data.table(x=x, y1=dLOGNO(x, mu=-6.714086, sigma=0.15), y2=dLOGNO(x, mu=-6.714086-0.15^2/2+0.075^2/2, sigma=0.075), y4=dLOGNO(x, mu=-6.714086-0.15^2/2+0.0375^2/2, sigma=0.0375))
		tmp		<- melt(tmp, id.var='x')
		ggplot(tmp, aes(x=x, y=value, group=variable, colour=variable)) + geom_line() + scale_x_continuous(breaks=seq(0,0.003,0.0005))
		x		<- seq(0.0005, 0.003, 0.00001)
		tmp		<- data.table(x=x, y1=dLOGNO(x, mu=log(0.002239075)-0.05^2/2, sigma=0.05), y5=dLOGNO(x, mu=log(0.002239075-0.0005)-0.065^2/2, sigma=0.065), y10=dLOGNO(x, mu=log(0.002239075-0.001)-0.09^2/2, sigma=0.09))
		
		x		<- seq(0.0005, 0.01, 0.00001)
		tmp		<- data.table(	x=x, 
								y5=dLOGNO(x, mu=log(0.002239075)-0.05^2/2, sigma=0.05), y7=dLOGNO(x, mu=log(0.002239075)-0.07^2/2, sigma=0.07), y10=dLOGNO(x, mu=log(0.002239075)-0.1^2/2, sigma=0.1), y13=dLOGNO(x, mu=log(0.002239075)-0.13^2/2, sigma=0.13), 
								W441=dLOGNO(x, mu=log(0.00447743)-0.3^2/2, sigma=0.3), W442=dLOGNO(x, mu=log(0.00447743)-0.4^2/2, sigma=0.4), W443=dLOGNO(x, mu=log(0.00447743)-0.5^2/2, sigma=0.5))
		#	times 2				
		x		<- seq(0.0005, 0.01, 0.00001)
		tmp		<- data.table(	x=x, 
								TransmissionLineage=dLOGNO(x, mu=log(0.002239075)-0.13^2/2, sigma=0.13),
								TransmissionLineage2=dLOGNO(x, mu=log(0.002239075)-0.3^2/2, sigma=0.3),
								TransmissionLineage3=dLOGNO(x, mu=log(0.002239075)-0.2^2/2, sigma=0.2),
								TipBranch=dLOGNO(x, mu=log(0.00447743)-0.5^2/2, sigma=0.5)
								)
		#	times 1.5				
		#tmp		<- data.table(	x=x, 
		#						TransmissionLineage=dLOGNO(x, mu=log(0.002239075)-0.13^2/2, sigma=0.13), 
		#						TipBranch=dLOGNO(x, mu=log(0.003358613)-0.3^2/2, sigma=0.3))						
		#tmp		<- data.table(	x=x, 
		#						TransmissionLineage=dLOGNO(x, mu=log(0.002239075)-0.01^2/2, sigma=0.01), 
		#						TipBranch=dLOGNO(x, mu=log(0.003)-0.2^2/2, sigma=0.2))
						
		tmp		<- melt(tmp, id.var='x')
		ggplot(tmp, aes(x=x, y=value, group=variable, colour=variable)) + geom_line() + scale_x_continuous(breaks=seq(0,0.02,0.001))		
	}
	if(0)
	{
		require(gamlss)		
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
		#	get posterior distribution of overall meanRate across genes
		tmp		<- log.df[, list(meanRate=mean(meanRate)), by='GENE'][, mean(meanRate)]	
		log.df	<- log.df[, list(meanRate=meanRate/mean(meanRate)*tmp), by='GENE']
		#ggplot(log.df, aes(x=meanRate)) + geom_histogram() + facet_grid(.~GENE, scales='free')
		#	get lognormal density parameters
		sd.log	<- log.df[, sd(log(meanRate))]			# 0.1499262
		mean.log<- log.df[, mean(log(meanRate))]		# -6.113092; 0.002239075
		#	mean is exp( mean.log+sd.log*sd.log/2 ); median is exp(mean.log) 
	}
	rERbw.args	<- list(meanlog=er.meanlog, sdlog=er.sdlog)
	if(er.sdlog>0)
	{
		rERbw		<- function(n, rERbw.args)
				{	
					rlnorm(n, meanlog=rERbw.args$meanlog, sdlog=rERbw.args$sdlog)		
				}		
	}
	if(er.sdlog==0)
	{
		rERbw		<- function(n, rERbw.args)
				{	
					rep( exp(rERbw.args$meanlog), n)		
				}		
	}	
	list(rERbw=rERbw, rERbw.args=rERbw.args)	
}
##--------------------------------------------------------------------------------------------------------
#	Function to add gaps into sequences  	
#	olli originally written 23-06-2015
##--------------------------------------------------------------------------------------------------------
PANGEA.align.from.basefreq<- function(indir.refalgn, infile.refalgn, indir.basefreq, callconsensus.minc, outdir, verbose=1)
{		
	callconsensus.cmd	<- 'python /Users/Oliver/Dropbox\\ \\(Infectious\\ Disease\\)/PANGEA-BEEHIVE-SHARED/ChrisCode/CallConsensus.py'
	callconsensus.maxc	<- callconsensus.minc
	mergealignments.cmd	<- 'python /Users/Oliver/Dropbox\\ \\(Infectious\\ Disease\\)/PANGEA-BEEHIVE-SHARED/ChrisCode/MergeAlignments2.0.py'
	mergealignments.opt	<- '-d'
				
	#	create consensus at coverage cut off		
	files				<- data.table(FILE=list.files(indir.basefreq, pattern='dat$'))	
	invisible(files[,{								
						tmp					<- gsub(' ','\\ ',gsub(')','\\)',gsub('(','\\(',paste(indir.basefreq, '/', FILE, sep=''),fixed=1),fixed=1),fixed=1)				
						cmd					<- paste(callconsensus.cmd,' ',tmp,' ', callconsensus.minc, ' ', callconsensus.maxc, sep='')
						#cat(cmd)				
						ans					<- system(cmd, intern=TRUE)	
						file				<- paste('TMP_',gsub('BaseFreqs.dat',paste('consensus_C',callconsensus.minc,'.fa',sep=''), FILE),sep='')
						cat(paste(ans, collapse='\n'), file=paste(outdir, '/', file, sep=''))
						NULL
					}, by='FILE'])
	#	add to reference alignment
	files				<- data.table(FILE=list.files(outdir, pattern=paste('consensus_C',callconsensus.minc,'.fa',sep='')))
	ps					<- files[,{
				cmd	<- paste(mergealignments.cmd,' ',mergealignments.opt,' ',indir.refalgn,'/',infile.refalgn,' ',outdir,'/',FILE,sep='')
				#cat(cmd)				
				ans	<- system(cmd, intern=TRUE)	
				list(SEQ=paste(ans[-1], collapse=''))				
			}, by='FILE']
	#	convert consensus alignment to DNAbin
	set(ps, NULL, 'FILE', ps[, gsub('TMP_','',FILE)])
	set(ps, NULL, 'FILE', ps[, gsub('_consensus.*','',FILE)])	
	tmp					<- tolower(do.call('rbind',strsplit(ps[, SEQ],'')))
	rownames(tmp)		<- ps[, FILE]
	ps					<- as.DNAbin(tmp)
	#	clean up
	tmp			<- list.files(outdir, pattern='^TMP_', full.names=TRUE)
	invisible(file.remove(tmp))	
	
	ps
}
##--------------------------------------------------------------------------------------------------------
#	Function to add gaps into sequences  	
#	olli originally written 23-06-2015
##--------------------------------------------------------------------------------------------------------
PANGEA.add.gaps<- function(indir.simu, indir.gap, infile.simu, infile.gap, gap.country, gap.symbol, gap.seed, verbose=1)
{		
	#	find files 		
	files			<- data.table(FILE=list.files(indir.simu, pattern='\\.fa.*$'))
	files			<- subset(files, !grepl('^TMP',FILE) & grepl(infile.simu, FILE, fixed=1))
	if(files[, any(grepl('gag', FILE))])
	{
		files[, SIMU:=files[, regmatches(FILE,regexpr('TRAIN[0-9]', FILE))]]
		files[, GENE:=files[, sapply(strsplit(FILE,'_'), '[[', 5)]]
		set(files, NULL, 'GENE', toupper(files[, substr(GENE,1,nchar(GENE)-3)]))
		set(files, NULL, 'GENE', files[, factor(GENE, levels=c('GAG','POL','ENV'), labels=c('GAG','POL','ENV'))])
		setkey(files, SIMU, GENE)
		files			<- subset(files, grepl(infile.simu, FILE, fixed=1) & !is.na(GENE))
		#	concatenate simu seqs
		files[, {
					gag	<- read.dna( paste(indir.simu, FILE[1], sep='/'), format='fasta' )
					pol	<- read.dna( paste(indir.simu, FILE[2], sep='/'), format='fasta' )
					env	<- read.dna( paste(indir.simu, FILE[3], sep='/'), format='fasta' )
					tmp	<- cbind(gag, pol, env)
					write.dna( tmp, file=paste(indir.simu, '/', 'TMP', gsub('gag','conc',FILE[1]), sep=''), format='fasta', colsep='', nbcol=-1)
					NULL
				}, by='SIMU']
		infile.simu	<- paste('TMP', gsub('gag','conc',files[1,FILE]), sep='')
	}
	else
		infile.simu	<- files[1,FILE]
	#
	#	add to fixed PANGEA alignment
	#
	ps				<- read.dna(paste(indir.gap,'/',infile.gap,sep=''), format='fasta')
	tmp				<- ps[(nrow(ps)-20):nrow(ps),]
	tmp				<- as.character(tmp)
	tmp[tmp=='?']	<- 'n'	#need this for MAFFT
	tmp				<- as.DNAbin(tmp)
	infile.old		<- paste('TMP',gsub('\\.','_last20\\.',infile.gap),sep='')
	write.dna(tmp, file=paste(indir.simu, infile.old, sep='/'), format='fasta', colsep='', nbcol=-1)	
	infile.old	<- paste(indir.simu, infile.old, sep='/')	
	#	mafft --thread 1 --treeout --reorder --keeplength --mapout --ep 0.0 --add new_sequences input > output
	infile.new	<- paste(indir.simu,'/',infile.simu, sep='')	
	file		<- paste(indir.simu,'/','TMP', gsub('\\.','mergedpartial\\.',infile.simu), sep='')
	cmd			<- paste('mafft --thread 4 --treeout --reorder --keeplength --mapout --ep 0.0 --add ',infile.new,' ',infile.old,' > ',file, sep='')
	cat('\ncalling')
	cat(cmd)
	system(cmd)	
	#	remove missing taxa from PANGEA alignment
	ps			<- as.character(ps)
	tmp			<- apply(ps, 1, function(x) !all(x%in%c('-','?','n'))	)
	ps			<- as.DNAbin( ps[tmp, ] )
	#	merge all
	file		<- list.files(indir.simu, pattern='mergedpartial')
	stopifnot(length(file)==1)
	tmp	<- read.dna(paste(indir.simu, file, sep='/'), format='fasta')
	tmp	<- tmp[grepl('IDPOP|HOUSE', rownames(tmp)),]
	psm	<- rbind(ps, tmp)
	#
	#	randomly select gaps from sequences in infile.gap from non-missing sequences
	#
	set.seed(gap.seed)
	x			<- as.character(psm)	
	tmp			<- which(grepl(gap.country,rownames(x)))
	gp.df		<- data.table(IDXSIM=which(grepl('IDPOP|HOUSE',rownames(x))))
	gp.df[, IDXGAP:=sample(tmp, nrow(gp.df), replace=TRUE)]
	for(i in seq_len(nrow(gp.df)))
	{
		z						<- which(x[ gp.df[i, IDXGAP], ]==gap.symbol)
		x[gp.df[i, IDXSIM],z]	<- gap.symbol
	}
	x			<- x[ grepl('IDPOP|HOUSE', rownames(x)), ]
	#	remove leading and trailing gap rows
	tmp			<- which(apply(x, 2, function(z) !all(z%in%c('-','?','n'))))
	sgp			<- as.DNAbin(x[, seq.int(tmp[1], tmp[length(tmp)])])
	cat(paste('\nSeq length without leading/trailing gaps and with seq coverage (no first/last common "?","-"), n=', ncol(sgp)))
	#	clean up
	tmp			<- list.files(indir.simu, pattern='^TMP', full.names=TRUE)
	file.remove(tmp)
	
	sgp	
}
##--------------------------------------------------------------------------------------------------------
#	Function to add gaps into sequences  	
#	olli originally written 01-07-2015
##--------------------------------------------------------------------------------------------------------
PANGEA.add.gaps.maintain.triplets<- function(indir.simu, indir.gap, infile.simu, infile.gap, gap.country, gap.symbol, gap.seed, outfile=NA, verbose=1)
{			
	#	find files 		
	files			<- data.table(FILE=list.files(indir.simu, pattern='\\.fa.*$'))
	files			<- subset(files, !grepl('RMGPS',FILE) & !grepl('map$',FILE) & !grepl('^TMP',FILE) & grepl(infile.simu, FILE, fixed=1))
	if(files[, any(grepl('gag', FILE))])
	{
		files[, SIMU:=files[, regmatches(FILE,regexpr('TRAIN[0-9]', FILE))]]
		files[, GENE:=files[, sapply(strsplit(FILE,'_'), '[[', 5)]]
		set(files, NULL, 'GENE', toupper(files[, substr(GENE,1,nchar(GENE)-3)]))
		set(files, NULL, 'GENE', files[, factor(GENE, levels=c('GAG','POL','ENV'), labels=c('GAG','POL','ENV'))])
		setkey(files, SIMU, GENE)
		files			<- subset(files, grepl(infile.simu, FILE, fixed=1) & !is.na(GENE))
		#	concatenate simu seqs
		files[, {
					gag	<- read.dna( paste(indir.simu, FILE[1], sep='/'), format='fasta' )
					pol	<- read.dna( paste(indir.simu, FILE[2], sep='/'), format='fasta' )
					env	<- read.dna( paste(indir.simu, FILE[3], sep='/'), format='fasta' )
					tmp	<- cbind(gag, pol, env)
					write.dna( tmp, file=paste(indir.simu, '/', 'TMP', gsub('gag','conc',FILE[1]), sep=''), format='fasta', colsep='', nbcol=-1)
					NULL
				}, by='SIMU']
		infile.simu	<- paste('TMP', gsub('gag','conc',files[1,FILE]), sep='')
	}
	else
		infile.simu	<- files[1,FILE]
	#
	#	add to fixed PANGEA alignment
	#
	ps				<- read.dna(paste(indir.gap,'/',infile.gap,sep=''), format='fasta')
	tmp				<- ps[(nrow(ps)-20):nrow(ps),]
	tmp				<- as.character(tmp)
	tmp[tmp=='?']	<- 'n'	#need this for MAFFT
	tmp				<- as.DNAbin(tmp)
	infile.old		<- paste('TMP',gsub('\\.','_last20\\.',infile.gap),sep='')
	write.dna(tmp, file=paste(indir.simu, infile.old, sep='/'), format='fasta', colsep='', nbcol=-1)	
	infile.old	<- paste(indir.simu, infile.old, sep='/')	
	#	mafft --thread 1 --treeout --reorder --keeplength --mapout --ep 0.0 --add new_sequences input > output
	infile.new	<- paste(indir.simu,'/',infile.simu, sep='')	
	file		<- paste(indir.simu,'/','TMP', gsub('\\.','mergedpartial\\.',infile.simu), sep='')
	#	do not force length: mafft --reorder --anysymbol --add new_sequences --auto input
	cmd			<- paste('mafft --thread 4 --treeout --reorder --keeplength --mapout --ep 0.0 --add ',infile.new,' ',infile.old,' > ',file, sep='')
	cat('\ncalling')
	cat(cmd)
	system(cmd)	
	#	remove missing taxa from PANGEA alignment
	ps			<- as.character(ps)
	tmp			<- apply(ps, 1, function(x) !all(x%in%c('-','?','n'))	)
	ps			<- as.DNAbin( ps[tmp, ] )
	#	merge all
	file		<- list.files(indir.simu, pattern='mergedpartial')
	stopifnot(length(file)==1)
	tmp	<- read.dna(paste(indir.simu, file, sep='/'), format='fasta')
	tmp	<- tmp[grepl('IDPOP|HOUSE', rownames(tmp)),]
	psm	<- rbind(ps, tmp)
	x	<- as.character(psm)
	#	remove leading and trailing gap rows
	tmp			<- which(apply(x[ grepl('IDPOP|HOUSE', rownames(x)), ], 2, function(z) !all(z%in%c('-','?','n'))))
	cat(paste('\nfirst non-gap',tmp[1],'last non-gap',tmp[length(tmp)]))
	x			<- x[, seq.int(tmp[1], tmp[length(tmp)])]
	#	rm gap columns that are in the PANGEA alignment and in the simulations
	#	(should not be in PANGEA alignment)
	tmp			<- which(apply(x, 2, function(z) all(z%in%c('-','?','n'))))
	cat(paste('\nrm common gap columns, n=',length(tmp)))
	x			<- x[,-tmp]		 
	#	check that gaps are only added in triples of 3
	#	not true. curating alignment manually: can explain all 0110 discrepancies. 
	#	remove columns that break the codon structure
	#exp.nontriplegaps<- c(388, 1476)
	#tmp			<- paste(as.numeric(x[nrow(x),]=='-'),collapse='')
	#tmp			<- gregexpr('0110',tmp)[[1]] + 1
	#cat(paste('\nstarting positions of non-triplet gaps',paste(tmp, collapse=', ')))
	#cat(paste('\nexpect starting positions of non-triplet gaps',paste(exp.nontriplegaps, collapse=', ')))
	#stopifnot( all(tmp==exp.nontriplegaps) )
	#x			<- x[,-c(tmp, tmp+1)]
	#	check that gaps are only added in triples of 3 
	xx			<- paste(as.numeric(x[nrow(x),]=='-'),collapse='')
	dfg			<- data.table(STRT=gregexpr('01',xx)[[1]] + 1)	#gap start positions
	tmp			<- strsplit(xx,'0')[[1]]
	dfg[, LEN:=sapply( tmp[tmp!=""], nchar )]
	dfg[, CUT:= LEN%%3]
	dfg			<- subset(dfg, CUT>0)
	if(nrow(dfg))
	{
		tmp		<- unlist(lapply( seq_len(nrow(dfg)), function(i)	seq.int(from=dfg[i,STRT], length.out=dfg[i,CUT])))
		cat(paste('\npositions of large non-triplet gaps that are removed',paste(tmp, collapse=', ')))
		x		<- x[,-tmp]
	}
	if(!is.na(outfile))		
	{
		#	keep merged alignment for back-compatibility
		write.dna(as.DNAbin(x),file=paste(indir.simu,'/',gsub('\\.fa','_RMGPS\\.fa',outfile),sep=''),format='fasta', colsep='', nbcol=-1)	
	}
	#	check that total no of gaps is now divisible by 3
	tmp			<- apply( x[ grepl('IDPOP|HOUSE', rownames(x)), ], 1, function(z)	sum(z=='-')	)	
	#stopifnot( all(tmp%%3==0) )
	#
	#	randomly select gaps from sequences in infile.gap from non-missing sequences
	#
	set.seed(gap.seed)
	tmp			<- which(grepl(gap.country,rownames(x)))
	dfg		<- data.table(IDXSIM=which(grepl('IDPOP|HOUSE',rownames(x))))
	dfg[, IDXGAP:=sample(tmp, nrow(dfg), replace=TRUE)]
	for(i in seq_len(nrow(dfg)))
	{
		z						<- which(x[ dfg[i, IDXGAP], ]==gap.symbol)
		x[dfg[i, IDXSIM],z]	<- gap.symbol
	}
	x			<- x[ grepl('IDPOP|HOUSE', rownames(x)), ]
	sgp			<- as.DNAbin(x)	
	#	clean up
	tmp			<- list.files(indir.simu, pattern='^TMP', full.names=TRUE)
	file.remove(tmp)	
	sgp		
}
##--------------------------------------------------------------------------------------------------------
#	return evolutionary rate modifier sampler for transmission edges	
#	olli originally written 21-08-2014
##--------------------------------------------------------------------------------------------------------
PANGEA.BetweenHostEvolutionaryRateModifier.create.sampler.v1<- function(bwerm.mu=1.5, bwerm.sigma=0.12)
{
	require(gamlss)
	if(0)
	{
		#from Vrancken et al:
		#c(3.5/2.5, 7.0/4.2) #1.4, 1.67	draw this from lognormal with mean 1.5 and variance so that tail captures 1.8		
		x		<- seq(1.01, 15, 0.01)
		tmp		<- data.table(x=x, y5=dLOGNO(x, mu=log(1.5), sigma=0.12), y8=dLOGNO(x, mu=log(1.75), sigma=0.103), y2=dLOGNO(x, mu=log(2), sigma=0.09))
		tmp		<- data.table(x=x, y4=dLOGNO(x, mu=log(4), sigma=0.045), y3=dLOGNO(x, mu=log(3), sigma=0.06), y2=dLOGNO(x, mu=log(2), sigma=0.09))
		tmp		<- data.table(x=x, y4=dLOGNO(x, mu=log(4), sigma=0.06), y3=dLOGNO(x, mu=log(3), sigma=0.06), y2=dLOGNO(x, mu=log(2), sigma=0.06))		
		tmp		<- data.table(x=x, y6=dLOGNO(x, mu=log(6)+0.4^2/2, sigma=0.4), y5=dLOGNO(x, mu=log(5)+0.4^2/2, sigma=0.4), y4=dLOGNO(x, mu=log(4)+0.4^2/2, sigma=0.4), y3=dLOGNO(x, mu=log(3)+0.4^2/2, sigma=0.4))
		tmp		<- melt(tmp, id.var='x')
		ggplot(tmp, aes(x=x, y=value, group=variable, colour=variable)) + geom_line() + scale_x_log10(breaks=seq(1,20,1))
	}
	rER.bwm<- function(n)
	{
		rLOGNO(n, mu=log(bwerm.mu), sigma=bwerm.sigma)
	}
	rER.bwm
}
##--------------------------------------------------------------------------------------------------------
#	return within host evolutionary rate sampler	
#	olli originally written 21-08-2014
##--------------------------------------------------------------------------------------------------------
#' @title Create sampler of Within Host Evolutionary Rates 
#' @description Returns a function to draw within host evolutionary rates. Currently modelled with a log normal density.
#' @param wher.mu		mean of the log normal density
#' @param wher.sigma 	standard deviation of the log normal density
#' @return R function
PANGEA.WithinHostEvolutionaryRate.create.sampler.v1<- function(wher.mu=log(0.005), wher.sigma=0.8)
{
	require(gamlss)
	if(0)
	{	
		#extremely basic model of within host evolutionary rate from HIV-1B pol estimates in the literature
		#median log10 ER of pol in Alizon / Fraser
		df.er	<- data.table(ER= 10^c(-1.85, -2.2, -2.5, -2.7, -2.72, -3.2), GENE='POL')		
		tmp		<- gamlss(ER~1, data=df.er, family=LOGNO)
		x		<- seq(0.0005, 0.02, 0.0001)
		tmp		<- data.table(x=x, y5=dLOGNO(x, mu=log(0.005), sigma=0.8), y4=dLOGNO(x, mu=log(0.004), sigma=0.7), y3=dLOGNO(x, mu=log(0.003), sigma=0.6))
		tmp		<- data.table(x=x, y83=dLOGNO(x, mu=-4.784295+0.045, sigma=0.3), y62=dLOGNO(x, mu=-5.071977+0.08, sigma=0.4), y41=dLOGNO(x, mu=-5.477443+0.18, sigma=0.6))
		#	should have been '-' all along:
		#	mean(rLOGNO(1e4, mu=-5.071977+0.4^2/2, sigma=0.4))
		tmp		<- data.table(x=x, y62=dLOGNO(x, mu=-5.071977+0.08, sigma=0.4), y62b=dLOGNO(x, mu=-5.071977-0.4^2/2, sigma=0.4) )
		tmp		<- data.table(x=x, y62=dLOGNO(x, mu=log(0.006716145)-0.37^2/2, sigma=0.37), y621=dLOGNO(x, mu=-5.071977-0.2^2/2, sigma=0.2), y441=dLOGNO(x, mu=log(0.00447743)-0.3^2/2, sigma=0.3) )		
		tmp		<- melt(tmp, id.var='x')
		ggplot(tmp, aes(x=x, y=value, group=variable, colour=variable)) + geom_line() + scale_x_continuous(breaks=seq(0,0.02,0.002))
		#correlation model: need high multiplier if ER high
		require(compositions)		
		y		<- rlnorm.rplus(1e4, meanlog=c(-5.071977+0.08, log(3)+0.4^2/2),  varlog=diag(c(0.6,0.4))%*%matrix(c(1,0.95,0.95,1),2,2)%*%diag(c(0.6,0.4)) )		
		tmp		<- data.table(wh=y[,1], bm=y[,2])
		ggplot(tmp, aes(x= wh, y=bm)) + geom_point()
		#explore marginal ER along branches
		tmp		<- do.call('rbind',lapply(c(3,4,5,6), function(bm)
				{
					tmp	<- rlnorm.rplus(1e4, meanlog=c(-5.071977+0.08, log(bm)+0.4^2/2),  varlog=diag(c(0.6,0.4))%*%matrix(c(1,0.95,0.95,1),2,2)%*%diag(c(0.6,0.4)) )
					data.table(value=tmp[,1]/tmp[,2], variable=paste('y',bm,sep=''))
				}))
		ggplot(tmp, aes(x=value, fill=variable)) + geom_histogram(binwidth=0.00025) + facet_grid(.~variable, scales='free')
	}	
	if(wher.sigma>0)
	{
		rER.pol<- function(n)
			{		
				ans	<- numeric(0)
				while(length(ans)<n)
				{
					tmp		<- rLOGNO(2*n, mu=wher.mu, sigma=wher.sigma)
					tmp[which(tmp>0.02)]	<- NA
					ans		<- c(ans, na.omit(tmp))						
				}			
				ans[seq_len(n)]
			}
	}
	if(wher.sigma==0)
	{
		rER.pol<- function(n)
			{		
				rep(exp(wher.mu), n)			
			}
	}
	rER.pol
}
######################################################################################
#	return GAG POL ENV ancestral sequences from BEAST PARSER output	
#	olli originally written 06-08-2014
#	tree 		beast trees in ape format, needed to compute calendar time for each ancestral sequence
#	node.stat	data.table containing meta information in nexus file for nodes
#	bseq		data.table containing original sequences. only needed for BEAST decompression.
#	return 		list of GAG POL ENV sequences in ape format 
PANGEA.RootSeqSim.get.ancestral.seq.withDecompression<- function(tree, node.stat, bseq, tree.id.sep='_', tree.id.idx.mcmcit=2, tree.id.burnin=1, label.sep='|', label.idx.ctime=2)
{
	require(data.table)
	require(ape)
	
	tree.id				<- names(tree)
	#	add calendar time for inner node at NODE_ID to node.stat
	node.stat[, CALENDAR_TIME:=NA_real_]		
	setkey(node.stat, TREE_ID, NODE_ID)
	for(i in seq_along(tree.id))
	{
		label.ctime			<- sapply( strsplit(tree[[i]]$tip.label, label.sep, fixed=TRUE), '[[', label.idx.ctime ) 
		label.ctime			<- as.numeric(label.ctime)			
		depth				<- node.depth.edgelength( tree[[ i ]] )
		tmp					<- which.max(depth)
		depth				<- depth-depth[tmp]+label.ctime[tmp]
		for(j in seq_along(depth))
			set(node.stat, node.stat[, which(TREE_ID==tree.id[i] & NODE_ID==j)], 'CALENDAR_TIME', depth[j])					
	}
	tmp			<- node.stat[, length(which(is.na(CALENDAR_TIME)))]
	cat(paste('\nTotal node statistics with no CALENDAR_TIME [should be zero], n=', tmp  ))
	#	keep only inner nodes
	tmp			<- lapply(seq_along(tree.id), function(i)
			{
				subset(node.stat, TREE_ID==tree.id[i] & NODE_ID>Ntip(tree[[i]]))
			})
	node.stat	<- do.call('rbind',tmp)
	set(node.stat, NULL, 'VALUE', node.stat[, gsub('\"','',VALUE)])
	#
	#	reconstruct ancestral sequences, need to decompress patterns that were compressed with BEAST
	#	TODO ** this results in duplicate columns and should be removed at a later point **
	#
	bseq			<- merge(bseq, bseq[, list(SEQ_N=nchar(SEQ)), by=c('GENE','TAXON_ID')], by=c('GENE','TAXON_ID'))
	bseq			<- bseq[, {
				tmp<- unlist(strsplit(SEQ,''))
				list(	CP1= paste(tmp[seq.int(1,length(tmp),by=3)], collapse='',sep=''), 
						CP2= paste(tmp[seq.int(2,length(tmp),by=3)], collapse='',sep=''), 
						CP3= paste(tmp[seq.int(3,length(tmp),by=3)], collapse='',sep='') 	)
			}, by=c('TAXON_ID','GENE')]		
	bseq			<- melt(bseq, measure.var=c('CP1','CP2','CP3'), variable.name='CODON_POS', value.name='SEQ')
	#	get index of orginal patterns and duplicate patterns
	bseq.decompress	<- bseq[, {
				#print(GENE)
				#print(CODON_POS)
				tmp		<- do.call('rbind',strsplit(SEQ,''))
				tmp2	<- apply( tmp, 2, function(x) paste(x,sep='',collapse=''))	#identical patterns?
				tmp2	<- which(duplicated(tmp2))
				#for each duplicate, work out index of original
				tmp3	<- sapply(tmp2, function(j1) which(apply( tmp[,seq.int(1,j1-1,1), drop=FALSE] == tmp[,j1], 2, all))[1]	 )
				list(NSEQ=ncol(tmp), DUPLICATE_PATTERN=tmp2, MOTHER_PATTERN=tmp3)
			}, by=c('GENE','CODON_POS')]		
	set(bseq.decompress, bseq.decompress[, which(GENE=='env')], 'GENE', 'ENV')
	set(bseq.decompress, bseq.decompress[, which(GENE=='gag')], 'GENE', 'GAG')
	set(bseq.decompress, bseq.decompress[, which(GENE=='pol')], 'GENE', 'POL')
	set(bseq.decompress, NULL, 'xSTAT', bseq.decompress[, paste(GENE,CODON_POS,sep='.')])		
	#	reconstruct ancestral genes sequences - decompress patterns		
	ancseq	<- node.stat[,  {													
				#print(STAT)
				#TREE_ID<- 'STATE_0'
				#STAT<- 'GAG.CP3'
				tmp								<- subset(bseq.decompress, xSTAT==STAT)											
				seq								<- matrix(data=NA_character_, nrow=length(VALUE), ncol=tmp[,NSEQ])
				seq.compressed					<- setdiff( seq_len(ncol(seq)), tmp[, DUPLICATE_PATTERN])	
				#print(dim(seq))										
				tmp2							<- do.call('rbind',strsplit(VALUE,''))
				#print(dim(tmp2))
				stopifnot(length(seq.compressed)==ncol(tmp2))
				seq[, seq.compressed]			<- tmp2
				#print(seq[1,])
				seq[, tmp[, DUPLICATE_PATTERN]] <- seq[, tmp[, MOTHER_PATTERN]]
				#print(seq[1,])
				#stop()
				seq								<- apply(seq, 1, function(x) paste(x, sep='',collapse=''))
				list(TREE_ID=TREE_ID, NODE_ID=NODE_ID, CALENDAR_TIME=CALENDAR_TIME, SEQ=seq) 
			}, by=c('STAT')]
	#	checks of ancseq before we proceed
	tmp		<- ancseq[, list(NSEQ= nchar(SEQ)), by=c('TREE_ID', 'NODE_ID', 'STAT')]		
	stopifnot( tmp[, list(CHECK= all(NSEQ==NSEQ[1])), by='STAT'][, all(CHECK)] )
	set(tmp, NULL, 'GENE', tmp[, sapply(strsplit(STAT,'\\.'),'[[',1)])
	set(tmp, NULL, 'CODON_POS', tmp[, sapply(strsplit(STAT,'\\.'),'[[',2)])
	stopifnot( tmp[, list(CHECK=all(NSEQ==NSEQ[1])), by='GENE'][, all(CHECK)] )
	#	reconstruct genes from codon positions
	ancseq		<- dcast.data.table(ancseq, TREE_ID + NODE_ID + CALENDAR_TIME ~ STAT, value.var="SEQ")
	ancseq		<- ancseq[, {
				tmp		<- do.call('rbind',sapply(list(ENV.CP1,ENV.CP2,ENV.CP3), strsplit, ''))
				env		<- paste(as.vector(tmp), collapse='')
				tmp		<- do.call('rbind',sapply(list(GAG.CP1,GAG.CP2,GAG.CP3), strsplit, ''))
				gag		<- paste(as.vector(tmp), collapse='')
				tmp		<- do.call('rbind',sapply(list(POL.CP1,POL.CP2,POL.CP3), strsplit, ''))
				pol		<- paste(as.vector(tmp), collapse='')
				list(GAG=gag, POL=pol, ENV=env, CALENDAR_TIME=CALENDAR_TIME)
			}, by=c('TREE_ID','NODE_ID')]
	set(ancseq, NULL, 'LABEL', ancseq[, paste(TREE_ID, NODE_ID, round(CALENDAR_TIME,d=3), sep='|')])		
	#	remove tree id STATE_xx where xx is smaller than burn-in
	set(ancseq, NULL, 'BEAST_MCMC_IT', ancseq[, as.integer(sapply(strsplit(TREE_ID,tree.id.sep),'[[',tree.id.idx.mcmcit))])
	ancseq		<- subset(ancseq, BEAST_MCMC_IT>tree.id.burnin)
	#
	#	return DNAbin
	#
	ancseq.gag				<- tolower(do.call('rbind',strsplit(ancseq[, GAG],'')))
	rownames(ancseq.gag)	<- ancseq[, LABEL]
	ancseq.gag				<- as.DNAbin(ancseq.gag)		
	ancseq.pol				<- tolower(do.call('rbind',strsplit(ancseq[, POL],'')))
	rownames(ancseq.pol)	<- ancseq[, LABEL]
	ancseq.pol				<- as.DNAbin(ancseq.pol)		
	ancseq.env				<- tolower(do.call('rbind',strsplit(ancseq[, ENV],'')))
	rownames(ancseq.env)	<- ancseq[, LABEL]
	ancseq.env				<- as.DNAbin(ancseq.env)				
	#ancseq					<- cbind(ancseq.gag, ancseq.pol, ancseq.env)
	#
	list(GAG=ancseq.gag, POL=ancseq.pol, ENV=ancseq.env)
}
######################################################################################
#	return GAG POL ENV ancestral sequences from BEAST PARSER output	
#	olli originally written 09-09-2014
#	tree 		beast trees in ape format, needed to compute calendar time for each ancestral sequence
#	node.stat	data.table containing meta information in nexus file for nodes
#	return 		list of GAG POL ENV sequences in ape format 
PANGEA.RootSeqSim.get.ancestral.seq.pg<- function(tree, node.stat, tree.id.sep='_', tree.id.idx.mcmcit=2, tree.id.burnin=1, label.sep='|', label.idx.ctime=5)
{
	require(data.table)
	require(ape)
	
	tree.id				<- names(tree)
	#	add calendar time for inner node at NODE_ID to node.stat
	node.stat[, CALENDAR_TIME:=NA_real_]		
	setkey(node.stat, TREE_ID, NODE_ID)
	for(i in seq_along(tree.id))
	{
		cat(paste('\nProcess CALENDAR_TIME for tree id', tree.id[i]  ))
		label.ctime			<- sapply( strsplit(tree[[i]]$tip.label, label.sep, fixed=TRUE), '[[', label.idx.ctime ) 
		label.ctime			<- as.numeric(label.ctime)			
		depth				<- node.depth.edgelength( tree[[ i ]] )
		tmp					<- which.max(depth)
		depth				<- depth-depth[tmp]+label.ctime[tmp]
		tmp					<- node.stat[, which(TREE_ID==tree.id[i])]
		for(j in seq_along(depth))
		{
			set(node.stat, tmp[ node.stat[tmp, which(NODE_ID==j)] ], 'CALENDAR_TIME', depth[j])
		}								
	}
	tmp			<- node.stat[, length(which(is.na(CALENDAR_TIME)))]
	cat(paste('\nTotal node statistics with no CALENDAR_TIME [should be zero], n=', tmp  ))
	#	keep only inner nodes
	tmp			<- sapply(tree, Ntip)
	stopifnot(all(tmp==tmp[1]))
	node.stat	<- subset(node.stat, NODE_ID>tmp[1])	
	#
	set(node.stat, NULL, 'VALUE', node.stat[, gsub('\"','',VALUE)])
	#	checks of ancseq before we proceed
	tmp			<- node.stat[, list(NSEQ= nchar(VALUE)), by=c('TREE_ID', 'NODE_ID', 'STAT')]		
	stopifnot( tmp[, list(CHECK= all(NSEQ==NSEQ[1])), by='STAT'][, all(CHECK)] )
	set(tmp, NULL, 'GENE', tmp[, sapply(strsplit(STAT,'\\.'),'[[',1)])
	set(tmp, NULL, 'CODON_POS', tmp[, sapply(strsplit(STAT,'\\.'),'[[',2)])	
	tmp			<- dcast.data.table(tmp, TREE_ID + NODE_ID + GENE ~ CODON_POS, value.var='NSEQ')
	tmp			<- tmp[, list(CPM=min(CP1, CP2, CP3)), by=c('TREE_ID','NODE_ID','GENE')]
	stopifnot( tmp[, list(CHECK=all(CPM==CPM[1])), by='GENE'][, all(CHECK)] )
	setkey(tmp, GENE)
	#	truncate to following size of coding regions (if necessary)
	tmp			<- unique(tmp)[, list(STAT=paste(GENE,'.CP',1:3,sep=''), CPM=CPM), by='GENE']
	node.stat	<- merge(node.stat, subset(tmp, select=c(STAT, CPM)), by='STAT')
	set(node.stat, NULL, 'VALUE', node.stat[, substr(VALUE,1,CPM)])
	set(node.stat, NULL, 'CPM', NULL)
	set(node.stat, NULL, 'GENE', node.stat[, substr(STAT,1,nchar(STAT)-4)])
	set(node.stat, NULL, 'STAT', node.stat[, substr(STAT,nchar(STAT)-2,nchar(STAT))])
	#	reconstruct genes from codon positions
	node.stat	<- dcast.data.table(node.stat, TREE_ID + NODE_ID + GENE + CALENDAR_TIME ~ STAT, value.var="VALUE")
	node.stat	<- node.stat[, {
									tmp		<- do.call('rbind',sapply(list(CP1,CP2,CP3), strsplit, ''))
									tmp		<- paste(as.vector(tmp), collapse='')
									list(SEQ=tmp, GENE=GENE, CALENDAR_TIME=CALENDAR_TIME)
								}, by=c('TREE_ID','NODE_ID')]	
	set(node.stat, NULL, 'LABEL', node.stat[, paste(TREE_ID, NODE_ID, round(CALENDAR_TIME,d=3), sep='|')])		
	#	remove tree id STATE_xx where xx is smaller than burn-in
	set(node.stat, NULL, 'BEAST_MCMC_IT', node.stat[, as.integer(sapply(strsplit(TREE_ID,tree.id.sep),'[[',tree.id.idx.mcmcit))])
	node.stat			<- subset(node.stat, BEAST_MCMC_IT>tree.id.burnin)
	cat(paste('\nFound ancestral sequences, n=', nrow(node.stat)  ))
	tmp			<- node.stat[, unique(GENE)]
	stopifnot(length(tmp)==1)
	node.stat
}
##--------------------------------------------------------------------------------------------------------
##	program to generate files for Seq Gen from output of Matt s VirusTreeSimulator
##	olli originally written 18-09-2014
##	modified 17-01-2015
##--------------------------------------------------------------------------------------------------------
#' @title Program to generate \code{SeqGen} input files
#' @description The \code{prog.PANGEA.SeqGen.createInputFile} reads files from the virus tree simulator in directory \code{indir.vts} and writes input files for \code{SeqGen}
#' to directory \code{outdir}. The program reads simulated transmission chain phylogenies with branches in units of calendar time
#' for sampled and unsampled individuals in a transmission chain. Within host evolutionary rates are drawn from a distribution, and
#' within host branch lengths are converted into the expected number of substitutions along the branch. Transmission branches are
#' multiplied with a multiplier to allow for slower evolution between hosts. The multiplier is drawn from a distribution. Starting sequences
#' are drawn from a pool of precomputed sequences. GTR parameters are drawn from a distribution. This is all that s needed to specify 
#' the SeqGen input files for each transmission chain.
#' @return NULL. Input to call SeqGen is stored in an R file.
#' @example example/ex.seqgen.inputfilecreator.R
#' @export
PANGEA.SeqGen.createInputFile<- function(indir.epi, infile.epi, indir.vts, infile.prefix, outdir.sg, pipeline.args, verbose=1, with.plot=1, label.sep='|')
{
	EPS			<- 1e-12
	stopifnot( all( c('s.seed','wher.mu','wher.sigma','bwerm.mu','bwerm.sigma')%in%pipeline.args[, stat] ) )
	#
	set.seed( pipeline.args['s.seed',][, as.numeric(v)] )	
	#
	#	setup samplers
	#
	file				<- paste(indir.epi, '/', infile.epi, sep='')
	load(file)	#expect "df.epi"    "df.trms"   "df.inds"   "df.sample"
	#
	cat(paste('\ncreate sampler of evolutionary rates'))
	#	create sampler of within host evolutionary rates
	rER.pol				<- PANGEA.WithinHostEvolutionaryRate.create.sampler.v1(wher.mu=pipeline.args['wher.mu',][, as.numeric(v)], wher.sigma=pipeline.args['wher.sigma',][, as.numeric(v)])
	#	create sampler of between host evolutionary rate dampener
	tmp					<- PANGEA.TransmissionEdgeEvolutionaryRate.create.sampler(er.meanlog=pipeline.args['bwerm.mu',][, as.numeric(v)], er.sdlog=pipeline.args['bwerm.sigma',][, as.numeric(v)])
	rERbw				<- tmp$rERbw
	rERbw.args			<- tmp$rERbw.args
	#	create sampler of ancestral sequences
	cat(paste('\ncreate sampler of ancestral sequences'))
	tmp					<- subset( df.inds, !is.na(TIME_SEQ) )[, length(unique(IDCLU))]
	cat(paste('\nnumber of clusters for which a sequence is required, n=', tmp))
	tmp					<- ifelse(tmp<400, 0.5, ifelse(tmp<500, 1, 2))	
	cat(paste('\nUsing root.ctime.grace=', tmp))
	tmp					<- PANGEA.RootSeq.create.sampler(root.ctime.grace= tmp, sample.grace= 3 )	
	rANCSEQ				<- tmp$rANCSEQ
	rANCSEQ.args		<- tmp$rANCSEQ.args 	
	#	read GTR parameters
	log.df				<- PANGEA.GTR.params()
	if(pipeline.args['dbg.rER',][, as.numeric(v)]==1 )
	{
		cat(paste('\nSetting mus to mean per gene and codon_pos'))
		tmp				<- log.df[, list(mu= mean(mu)), by=c('GENE','CODON_POS')]
		#tmp[, ER:=mu*log.df$meanRate[1]]
		log.df			<- merge( subset(log.df, select=which(colnames(log.df)!='mu')), tmp, by=c('GENE','CODON_POS'))		
	}		
	if( pipeline.args['dbg.GTRparam',][, as.numeric(v)]==1 )
	{
		cat(paste('\nSetting GTR parameters to MEAN (except mu)'))
		tmp				<- mean	
		log.df			<- log.df[, list(state=state, mu=mu, alpha=tmp(alpha), at=tmp(at), ac=tmp(ac), cg=tmp(cg), ag=tmp(ag), gt=tmp(gt), meanRate=tmp(meanRate), a=tmp(a), c=tmp(c),  g=tmp(g), t=tmp(t) ), by=c('GENE','CODON_POS')]		
	}		
	log.df[, IDX:= seq_len(nrow(log.df))]
	log.df[, FILE:=NULL]	
	#
	#
	#	
	infiles				<- list.files(indir.vts)
	tmp					<- paste('^',infile.prefix,'.*nex$',sep='')
	infiles				<- sort(infiles[ grepl(tmp, infiles)  ])	
	if(!length(infiles))	stop('cannot find files matching criteria')
	#
	#	read from VirusTreeSimulator and convert branch lengths in time to branch lengths in subst/site
	#
	df.ph				<- vector('list', length(infiles))
	phd					<- vector('list', length(infiles))
	phs					<- vector('list', length(infiles))
	phd.plot			<- vector('list', length(infiles))
	df.nodestat			<- vector('list', length(infiles))
	cat(paste('\nUsing StartTimeMode',pipeline.args['index.starttime.mode',][,v]))
	root.edge.rate		<- NA
	if( pipeline.args['index.starttime.mode',][,v]=='shift' )
	{		
		root.edge.rate	<- 1e-6
		cat(paste('\nFix root edge rate to =',root.edge.rate))
	}		
	if( as.numeric(pipeline.args['root.edge.fixed',][,v] )) 
	{
		root.edge.rate	<- log.df[1,meanRate]
		cat(paste('\nFix root edge rate to =',root.edge.rate))		
	}	
	for(i in seq_along(infiles))
	{				
		#i	<- 40
		#i<- 47; i<- 9		
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
		#	produce collapsed tree with branch length in units of calendar time
		phd[[i]]		<- seq.collapse.singles(ph)
		#
		#	create collapsed Newick tree with expected substitutions / site for each branch 
		#
		#	draw within host evolutionary rates for every individual in the transmission chain, and	smaller ERs along the transmission lineages
		node.stat		<- merge(node.stat, data.table( IDPOP=node.stat[, unique(IDPOP)], ER= rER.pol(node.stat[, length(unique(IDPOP))]), BWM= rERbw(node.stat[, length(unique(IDPOP))], rERbw.args) ), by='IDPOP')
		#	re-set to previous notation in terms of BWM ( between host multiplier to ER, ie ER= within-host ER / BWM )
		set(node.stat, NULL, 'BWM', node.stat[, ER/BWM])	
		#	set BWM to 1 for all edges that are NOT leading to a transmission. 
		#	Because only one seq is sampled per patient, these are only edges that end in a tip.
		set(node.stat, node.stat[, which( NODE_ID%in%seq.int(1,Ntip(ph)) )], 'BWM', 1.)
		#	no ER possible for root node - there s no edge leading to it
		set(node.stat, node.stat[, which(NODE_ID==Ntip(ph)+1)], c('ER','BWM'), NA_real_)		
		#	set root edge evolutionary rate to overall mean between-host rate
		#	get NODE_ID of edge from root
		if(!is.na(root.edge.rate))
		{
			tmp				<- ph$edge[match(Ntip(ph)+1, ph$edge[1, ]), 2]
			tmp				<- node.stat[, which(NODE_ID==tmp)]		
			set(node.stat, tmp, 'ER', root.edge.rate )		
			set(node.stat, tmp, 'BWM', 1. )		# no need to further slow down root edge			
		}
		#	check root edge length
		if( pipeline.args['index.starttime.mode',][,v]=='fix45' )
		{
			stopifnot( ph$edge.length[ which( ph$edge[, 1] == Ntip(ph)+1 ) ]>=29.5 )			
		}
		#	check calendar time of root in simulated phylogeny for consistency
		tmp				<- seq.collapse.singles(ph)
		tmp2			<- regmatches(tmp$tip.label[1], regexpr('ID_[0-9]+',tmp$tip.label[1]))
		tmp2			<- as.numeric(substr(tmp2, 4, nchar(tmp2)))
		tmp2			<- subset(node.stat, IDPOP==tmp2)[1, TIME_SEQ]
		root.ctime		<- ifelse(Nnode(tmp), tmp2 - (node.depth.edgelength(tmp)[1] + tmp$root.edge), tmp2-tmp$root.edge)		
		tmp				<- subset(node.stat, IDPOP<0)[, unique(IDPOP)]
		stopifnot(length(tmp)==1)
		stopifnot(subset(df.trms, IDTR==tmp)[, round(IDTR_TIME_INFECTED, d=1)]==round(root.ctime, d=1))
		#	check if all sampling times are consistent with node height
		tmp				<- subset( node.stat, NODE_ID<=Ntip(ph) )
		setkey(tmp, NODE_ID)
		tmp2			<- seq.collapse.singles(ph) 
		if( Nnode(tmp2) )
			tmp[, NODE_DEPTH:=root.ctime + tmp2$root.edge + node.depth.edgelength(tmp2)[ seq_len(Ntip(tmp2)) ] ]
		if( Nnode(tmp2)==0 )
			tmp[, NODE_DEPTH:=root.ctime + tmp2$root.edge ]
		stopifnot( tmp[, max(abs(NODE_DEPTH-TIME_SEQ))<=1e-6 ] )
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
		ph$tip.label		<- node.stat[seq_len(Ntip(ph)), ][, LABEL]
		phd[[i]]$tip.label	<- node.stat[seq_len(Ntip(phd[[i]])), ][, LABEL]
		phd.plot[[i]]		<- seq.singleton2bifurcatingtree( phd[[i]] )
		phs[[i]]			<- seq.singleton2bifurcatingtree( ph )
		#phd[[i]]			<- seq.addrootnode( phd[[i]], dummy.label=paste('NOEXIST_IDCLU',node.stat[, unique(IDCLU)],'|NA|DOB_NA|',root.ctime,sep='') )
		#
		df.nodestat[[i]]	<- node.stat
		if(Nnode(ph))
		{
			tmp			<- write.tree(ph, digits = 10)			
			tmp			<- paste( '(',substr(tmp,1,nchar(tmp)-1),',NOEXIST_NA|NA|DOB_NA|',root.ctime,':0):0;', sep='')
			phd[[i]]	<- write.tree(phd[[i]], digits = 10)
		}
		if(!Nnode(ph))
		{
			tmp			<- paste( '(',ph$tip.label,':',ph$root.edge,',NOEXIST_NA|NA|DOB_NA|',root.ctime,':0):0;', sep='')
			phd[[i]]	<- paste( '(',phd[[i]]$tip.label,':',phd[[i]]$root.edge, ');', sep='' )
		}
		df.ph[[i]]		<- data.table(ROOT_CALENDAR_TIME= root.ctime, IDCLU=node.stat[, unique(IDCLU)], NEWICK=tmp)
		#readline()
	}
	#	write cluster trees of each transmission network into a single newick file
	phd			<- sapply(phd,'[[',1)
	file	<- paste(indir.epi, '/', substr(infile.epi,1,nchar(infile.epi)-6),'DATEDCLUTREES.newick', sep='')
	cat(paste('\nWrite to file=',file))
	writeLines(phd, con=file)
	#	get multifurcating tree with brl in units of calendar time
	options(expressions=1e4)
	phd			<- eval(parse(text=paste('phd.plot[[',seq_along(phd.plot),']]', sep='',collapse='+')))
	options(expressions=5e3)
	phd			<- drop.tip(phd, which(grepl('DUMMY', phd$tip.label)), root.edge=1)
	phd			<- ladderize(phd)
	#	write multifurcating tree with brl in units of calendar time
	file	<- paste(indir.epi, '/', substr(infile.epi,1,nchar(infile.epi)-6),'DATEDTREE.newick', sep='')
	cat(paste('\nWrite to file=',file))
	write.tree(phd, file=file)
	#	write multifurcating tree with brl in units of subst/site
	options(expressions=1e4)
	phs			<- eval(parse(text=paste('phs[[',seq_along(phs),']]', sep='',collapse='+')))
	options(expressions=5e3)
	phs			<- drop.tip(phs, which(grepl('DUMMY', phs$tip.label)), root.edge=1)
	phs			<- ladderize(phs)
	file	<- paste(indir.epi, '/', substr(infile.epi,1,nchar(infile.epi)-6),'SUBSTTREE.newick', sep='')
	cat(paste('\nWrite to file=',file))
	write.tree(phs, file=file)
	#
	df.ph		<- do.call('rbind',df.ph)
	df.nodestat	<- do.call('rbind',df.nodestat)
	#	
	#	
	#	check that we have exactly one root edge with overall mean between host rate per cluster
	if(!is.na(root.edge.rate))
		stopifnot( df.nodestat[, length(unique(IDCLU))]==nrow(subset(df.nodestat, ER==root.edge.rate)) )	
	
	#tmp	<- unique(subset(df.inds, IDPOP>=-110 & IDPOP<0, IDCLU))
	#tmp	<- merge(df.inds, tmp, by='IDCLU')[, length(which(!is.na(TIME_SEQ)))] / df.inds[, length(which(!is.na(TIME_SEQ)))]
	#cat(paste('\nProportion of sequences descending from no import after baseline=', tmp))
	
	#df.nodestat[, length(which(!is.na(TIME_SEQ))), by='IDCLU']
	#subset(df.nodestat, IDCLU==72)
	#subset(df.ph, IDCLU==72)
	#subset(df.ph, IDCLU==47)	
	if(with.plot)
	{
		#subset(df.nodestat, ER!=root.edge.rate)[, table(BWM==1.)]
		#ggplot(subset(df.nodestat, ER!=root.edge.rate) , aes(x=ER/BWM)) + geom_histogram(binwidth=0.001)
		#ggplot(subset(df.nodestat, ER!=root.edge.rate & BWM!=1.) , aes(x=ER/BWM)) + geom_histogram(binwidth=0.0001)
		#ggplot(subset(df.nodestat, ER!=root.edge.rate), aes(x=ER, y=BWM)) + geom_point()	
		#	plot used within-host ERs
		tmp		<- root.edge.rate
		if(is.na(tmp))
			tmp	<- Inf
		ggplot(subset(df.nodestat, ER!=tmp), aes(x=ER/BWM)) + geom_histogram(binwidth=0.001)	+ labs(x='simulated within-host evolutionary rate') +
				scale_x_continuous(breaks= seq(0, 0.02, 0.002))
		file	<- paste(indir.epi, '/', substr(infile.epi,1,nchar(infile.epi)-6),'INFO_sg_ER.pdf', sep='')
		cat(paste('\nWrite to file=',file))
		ggsave(file, w=6, h=6)
		#	plot used between host modifiers
		ggplot(subset(df.nodestat, ER!=tmp & BWM==1) , aes(x=ER/BWM)) + geom_histogram(binwidth=0.001) + labs(x='simulated within-host rate evolutionary rate\nwithout transmission edges') +
				scale_x_continuous(breaks= seq(0, 0.02, 0.002))
		file	<- paste(indir.epi, '/', substr(infile.epi,1,nchar(infile.epi)-6),'INFO_sg_BWM.pdf', sep='')
		cat(paste('\nWrite to file=',file))
		ggsave(file, w=6, h=6)
		#	plot used ERs along transmission edges
		ggplot(subset(df.nodestat, ER!=tmp & BWM!=1) , aes(x=ER/BWM)) + geom_histogram(binwidth=0.0001) + labs(x='simulated evolutionary rates along transmission edges') +
				scale_x_continuous(breaks= seq(0, 0.02, 0.0005))
		file	<- paste(indir.epi, '/', substr(infile.epi,1,nchar(infile.epi)-6),'INFO_sg_BWER.pdf', sep='')
		cat(paste('\nWrite to file=',file))
		ggsave(file, w=6, h=6)	
		
		if(grepl('fix',pipeline.args['index.starttime.mode',][,v]))
		{
			setkey(df.nodestat, IDPOP)
			tmp				<- merge(subset(unique(df.nodestat), select=c(IDCLU, LABEL)), unique(df.nodestat)[, list(IDCLU_N=length(which(!is.na(TIME_SEQ)))), by='IDCLU'], by='IDCLU')
			tmp[, LABEL_NEW:= tmp[, paste(IDCLU,'_',IDCLU_N,label.sep,LABEL,sep='')]]
			setkey(tmp, LABEL)			
			phd$tip.label	<- tmp[ phd$tip.label, ][, LABEL_NEW]
			file			<- paste(indir.epi, '/', substr(infile.epi,1,nchar(infile.epi)-6),'DATEDTREE.pdf', sep='')
			cat(paste('\nWrite to file=',file))
			pdf(file=file, w=10, h=Ntip(phd)*0.1)
			plot(phd, show.tip=TRUE, cex=0.5)			
			dev.off()
		}
	}
	#
	#	draw ancestral sequences and add to df.ph
	#
	if(pipeline.args['startseq.mode',][,v]=='many')
	{
		root.ctime		<- df.ph[, ROOT_CALENDAR_TIME]
		ancseq			<- rANCSEQ(root.ctime, rANCSEQ.args)
		ancseq			<- data.table(ANCSEQ= apply(as.character(ancseq),1,function(x) paste(x, collapse='')) )		
	}		
	if(pipeline.args['startseq.mode',][,v]=='one')
	{
		cat(paste('\nStartSeqModel=',pipeline.args['startseq.mode',][,v],'use first sampled starting sequence for all' ))
		stopifnot(max(abs(df.ph[1, ROOT_CALENDAR_TIME]-df.ph[, ROOT_CALENDAR_TIME]))<100*EPS)
		root.ctime		<- df.ph[1, ROOT_CALENDAR_TIME]
		tmp				<- rANCSEQ(root.ctime, rANCSEQ.args)
		tmp				<- data.table(ANCSEQ= apply(as.character(tmp),1,function(x) paste(x, collapse='')) )
		ancseq			<- data.table(ANCSEQ=rep(NA_character_, nrow(df.ph)))
		set(ancseq, NULL, 'ANCSEQ', tmp[1,ANCSEQ])
	}
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
#	Program to simulate sequences with Seq-Gen-1.3.3 	
#	olli originally written 27-01-2015
######################################################################################
#' @title Program to simulate gene sequences
#' @description \code{prog.PANGEA.SeqGen.run} reads file \code{infile.sg} in directory \code{indir.sg} that was
#' created with the \code{SeqGen} input file creator. The simulated partial sequences are collected, coerced back
#' into Gag, Pol, Env genes, and written in fasta format to directory \code{outdir}. Patient Metavariables are 
#' stored in the same directory, and zip files are created.
#' @return NULL. Saves zip files with simulations.
#' @example example/ex.seqgen.run.R
#' @export
PANGEA.SeqGen.run.v4<- function(indir.epi, infile.epi, indir.sg, infile.prefix, outdir, pipeline.args)
{	
	set.seed(pipeline.args['s.seed',][, as.numeric(v)])
	#
	file		<- paste(indir.epi, '/', infile.epi, sep='')
	load(file)	#expect "df.epi"    "df.trms"   "df.inds"   "df.sample"   "df.sp"
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
	stopifnot( nrow(df.ph.out)<65535 )	# running out of seeds?
	df.ph.out[, SEED:=sample.int(65535, nrow(df.ph.out))]
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
				cmd	<- cmd.SeqGen(indir.sg, FILE, indir.sg, gsub('seqgen','phy',FILE), prog=PR.SEQGEN, prog.args=paste('-n',1,' -k1 -or -z',SEED,sep=''), 
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
	save(df.epi, df.trms, df.inds, df.sample, df.seq, df.sp, file=file)
	#
	#	save simulated data -- to be shared
	#	
	if(pipeline.args['epi.model'][,v]=='HPTN071')
	{
		tmp			<- subset( df.inds, !is.na(TIME_SEQ), select=c(IDPOP, GENDER, DOB, DOD, DIAG_T, DIAG_CD4, ART1_T, ART1_CD4, TIME_SEQ, RECENT_TR ) )
		tmp2		<- tmp[, which(!is.na(RECENT_TR))]
		cat(paste('\nSet RECENT_TR to missing for p=',1-pipeline.args['report.prop.recent',][,as.numeric(v)]))
		tmp2		<- sample(tmp2, round(length(tmp2)*(1-pipeline.args['report.prop.recent',][,as.numeric(v)])))
		set(tmp, NULL, 'RECENT_TR', tmp[, as.character(RECENT_TR)])
		set(tmp, tmp2, 'RECENT_TR', NA_character_)
		set(tmp, NULL, 'RECENT_TR', tmp[, factor(RECENT_TR)])	
		
		set(tmp, NULL, 'GENDER', tmp[,as.character(GENDER)])
		tmp2		<- tmp[, which(is.na(DIAG_T) & TIME_SEQ<2000)]
		cat(paste('\nSet patient variables to NA for archival seq, n=',length(tmp2)))		
		set(tmp, tmp2, c('DOB','DOD'), NA_real_) 
		set(tmp, tmp2, 'GENDER', NA_character_)
		tmp2		<- tmp[, which(is.na(DIAG_T) & TIME_SEQ>=2000)]
		cat(paste('\nSet patient variables to NA after 200, n=',length(tmp2)))
		print(tmp[tmp2,])
		set(tmp, tmp2, c('DOB','DOD'), NA_real_) 
		set(tmp, tmp2, 'GENDER', NA_character_)
		set(tmp, NULL, 'GENDER', tmp[,factor(GENDER)])
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
	
	df.seq			<- merge(df.seq, df.seq[, list(IDCLU_N=length(IDPOP)), by='IDCLU'], by='IDCLU')
	if(with.NJ)
	{
		#
		#	create and plot NJ tree on conc seq
		#			
		#	load outgroup sequences
		file			<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfg_HXB2outgroup.R')
		cat(paste('\nLoading outgroup seq from file', file))
		load(file)		#expect "outgroup.seq.gag" "outgroup.seq.pol" "outgroup.seq.env"
		if(nrow(df.seq)<2000)
			df.seq.nj	<- copy(df.seq)
		if(nrow(df.seq)>=2000)
		{
			cat(paste('\nToo many seqs for quick NJ tree calculation, selecting first 2000'))
			df.seq.nj	<- unique(subset(df.seq, select=c(IDCLU, IDCLU_N)))			
			df.seq.nj[, IDCLU_CN:=df.seq.nj[, cumsum(IDCLU_N)]]
			df.seq.nj	<- df.seq.nj[seq_len(max(1, which(IDCLU_CN>2000)[1]-1)), ]
			df.seq.nj	<- merge(subset(df.seq.nj, select=IDCLU), df.seq, by='IDCLU')
		}
		#	concatenate sequences	
		df.seq.nj[, TIPCOLOR:='black']
		set(df.seq.nj, df.seq.nj[,which(IDCLU_N<3)],'TIPCOLOR','red')
		tmp				<- tolower(do.call('rbind',strsplit(df.seq.nj[, paste(GAG,POL,ENV,sep='')],'')))		
		rownames(tmp)	<- df.seq.nj[, paste(IDCLU,'_',IDCLU_N,treelabel.idx.sep,LABEL,sep='')]	
		seq				<- as.DNAbin(tmp)
		tmp				<- cbind(outgroup.seq.gag[,1:ncol(df.seq.gag)], outgroup.seq.pol, outgroup.seq.env)
		seq				<- rbind(seq,tmp)	
		seq.ph			<- nj(dist.dna(seq, model='raw'))		
		tmp				<- which(seq.ph$tip.label=="HXB2")
		seq.ph			<- reroot(seq.ph, tmp, seq.ph$edge.length[which(seq.ph$edge[,2]==tmp)])
		seq.ph			<- ladderize(seq.ph)
		file			<- paste(indir.sg, '/', infile.prefix, 'INFO_NJconc.pdf', sep='')	
		cat(paste('\nwrite to file',file))
		pdf(file=file, w=10, h=0.1*Ntip(seq.ph))
		plot(seq.ph, show.tip=TRUE, cex=0.5, tip.color=df.seq.nj[,TIPCOLOR])
		dev.off()	
		#
		#	create and plot NJ tree on pol seq
		#			
		tmp				<- tolower(do.call('rbind',strsplit(df.seq.nj[, POL],'')))		
		rownames(tmp)	<- df.seq.nj[, paste(IDCLU,'_',IDCLU_N,treelabel.idx.sep,LABEL,sep='')]	
		seq				<- as.DNAbin(tmp)		
		seq				<- rbind(seq,outgroup.seq.pol)	
		seq.ph			<- nj(dist.dna(seq, model='raw'))		
		tmp				<- which(seq.ph$tip.label=="HXB2")
		seq.ph			<- reroot(seq.ph, tmp, seq.ph$edge.length[which(seq.ph$edge[,2]==tmp)])
		seq.ph			<- ladderize(seq.ph)
		file			<- paste(indir.sg, '/', infile.prefix, 'INFO_NJpol.pdf', sep='')	
		cat(paste('\nwrite to file',file))
		pdf(file=file, w=10, h=0.1*Ntip(seq.ph))
		plot(seq.ph, show.tip=TRUE, cex=0.5, tip.color=df.seq.nj[,TIPCOLOR])
		dev.off()			
	}
	#
	#	zip simulated sequence files
	#
	tmp				<- c( paste(outdir, '/', infile.prefix, 'SIMULATED_metadata.csv', sep=''), paste(indir.epi, '/', gsub('SAVE.R','CROSS_SECTIONAL_SURVEY.csv',infile.epi), sep=''), paste(outdir, '/', infile.prefix, 'SIMULATED_env.fa', sep=''), paste(outdir, '/', infile.prefix, 'SIMULATED_pol.fa', sep=''), paste(outdir, '/', infile.prefix, 'SIMULATED_gag.fa', sep='') )
	zip( paste(outdir, '/', infile.prefix, 'SIMULATED_SEQ.zip', sep=''), tmp, flags = "-FSr9XTj")
	#
	#	zip simulated tree files
	#
	tmp2			<- c( paste(outdir, '/', infile.prefix, 'SIMULATED_metadata.csv', sep=''), paste(indir.epi, '/', gsub('SAVE.R','CROSS_SECTIONAL_SURVEY.csv',infile.epi), sep=''), 
			paste(indir.epi, '/', infile.prefix, 'DATEDTREE.newick', sep=''),
			paste(indir.epi, '/', infile.prefix, 'DATEDCLUTREES.newick', sep=''),
			paste(indir.epi, '/', infile.prefix, 'SUBSTTREE.newick', sep=''))
	zip( paste(outdir, '/', infile.prefix, 'SIMULATED_TREE.zip', sep=''), tmp2, flags = "-FSr9XTj")
	#
	tmp				<- unique(c(tmp, tmp2))
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
##--------------------------------------------------------------------------------------------------------
##	olli originally written 26-01-2015
##--------------------------------------------------------------------------------------------------------
#' @title HPTN071 parser (version 4, uses date of diagnosis)
#' @description Reads files from the epi simulator in directory \code{indir} and writes csv files
#' in directory \code{outdir} for the virus tree simulator. The program samples sequences according to
#' an exponentially increasing sampling fraction in the same way as \code{prog.HPTN071.input.parser.v1}.
#' In addition, transmissions are broken and treated as imported from outside the simulated population.
#' The infected of a broken transmission chain is considered a new index case of a transmission chain within the 
#' simulated population. All input arguments are specified via the \code{argv} 
#' string, see the Examples.
#' @return NULL. Saves simulations to file.
#' @example example/ex.seq.sampler.v4.R
#' @export
PANGEA.HPTN071.input.parser.v4<- function(indir, infile.ind, infile.trm, outdir, outfile.ind, outfile.trm, pipeline.args, verbose=1, with.plot=1)	
{
	require(data.table)			
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
	df.trm	<- as.data.table(read.csv(infile.trm, stringsAsFactors=FALSE, sep=',', dec='.'))
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
	setnames(df.ind, c(	"Id","Sex","DoB","DoD","RiskGp"), c('IDPOP','GENDER','DOB','DOD','RISK'))
	setnames(df.ind, 	c(	"HIV_pos", "t_diagnosed","cd4_diagnosis","cd4atfirstART","t_1stARTstart","t_1stVLsupp_start","t_1stVLsupp_stop"), 
			c( 'HIV', 'DIAG_T','DIAG_CD4','ART1_CD4','ART1_T',"VLS1_TS","VLS1_TE"))
	set(df.ind, NULL, 'GENDER', df.ind[, factor(GENDER)])
	set(df.ind, NULL, 'RISK', df.ind[, factor(RISK)])
	stopifnot( df.ind[, !any(is.na(DOB))] )
	set(df.ind, df.ind[, which(DOD==-1)], 'DOD', ceiling(max(df.trm$TIME_TR))+1.)
	set(df.ind, NULL, 'HIV', df.ind[, factor(HIV, levels=c(0,1), labels=c('N','Y'))])
	cat(paste('\nFound HIV+, n=', df.ind[, length(which(HIV=='Y'))]))
	set(df.ind, df.ind[, which(DIAG_T=='ND')], 'DIAG_T', NA_character_)
	set(df.ind, NULL, 'DIAG_T', df.ind[, as.numeric(DIAG_T)])
	set(df.ind, df.ind[, which(DIAG_CD4<0)], 'DIAG_CD4', NA_real_)
	stopifnot( df.ind[, !any(HIV=='N' & !is.na(DIAG_T)) ] )
	stopifnot( df.ind[, !any(!is.na(DIAG_T) & is.na(DIAG_CD4))] )
	cat(paste('\nFound not % undiagnosed , diagnosed=', paste( subset(df.ind, HIV=='Y')[, round(table(!is.na(DIAG_T)) / length(DIAG_T), d=2)], collapse=', ' )))	
	set(df.ind, df.ind[, which(ART1_CD4<0)], 'ART1_CD4', NA_real_)
	set(df.ind, df.ind[, which(ART1_T<0)], 'ART1_T', NA_real_)
	stopifnot( df.ind[, all(DIAG_T<ART1_T, na.rm=TRUE)] )
	cat(paste('\nFound not % on ART , not on ART among diagnosed=', paste(subset(df.ind, !is.na(DIAG_T))[, round( table(is.na(ART1_T))/length(ART1_T), d=2)], collapse=', ') ))
	set(df.ind, df.ind[, which(VLS1_TS<0)], 'VLS1_TS', NA_real_)
	set(df.ind, df.ind[, which(VLS1_TE<0)], 'VLS1_TE', NA_real_)
	tmp			<- df.ind[, which(ART1_T>=VLS1_TS)]
	cat(paste('\nFound ART1_T<VLS1_TS, setting VLS1_TS and VLS1_TE to NA, n=', length(tmp)))
	set(df.ind, tmp, 'VLS1_TS', NA_real_)
	set(df.ind, tmp, 'VLS1_TE', NA_real_)
	stopifnot( df.ind[, all(ART1_T<VLS1_TS, na.rm=TRUE)] )
	stopifnot( df.ind[, all(VLS1_TS<VLS1_TE, na.rm=TRUE)] )
	cat(paste('\nFound not % reached viral suppression , did not reach viral suppression among treated=',  paste(subset(df.ind, !is.na(ART1_T))[, round( table(is.na(VLS1_TS))/length(VLS1_TS), d=2)], collapse=', ') ))
	tmp			<- df.ind[, which(DIAG_CD4<ART1_CD4-DIAG_CD4*0.5 & DIAG_CD4>250)]
	cat(paste('\nFound individuals whose CD4 at ART start is much higher than at diagnosis, n=', length(tmp)))
	#stopifnot( df.ind[, all(DIAG_CD4>ART1_CD4, na.rm=TRUE)] )
	stopifnot( df.ind[, all(DIAG_CD4>0, na.rm=TRUE)] )
	stopifnot( df.ind[, all(ART1_CD4>0, na.rm=TRUE)] )	
	#	add transmission time
	tmp			<- subset(df.trm, select=c(IDREC, TIME_TR))	
	setnames(tmp, 'IDREC','IDPOP')
	df.ind		<- merge(df.ind, tmp, by='IDPOP', all.x=TRUE)
	stopifnot( df.ind[, all(TIME_TR<DIAG_T, na.rm=TRUE)] )
	stopifnot( df.ind[, !any(is.na(TIME_TR) & !is.na(DIAG_T))] )
	cat(paste('\nFound individuals with a valid record, n=', nrow(df.ind)))
	cat(paste('\nFound individuals with an infection event, n=', nrow(subset(df.ind,!is.na(TIME_TR)))))
	#	reset times if needed, because TIME_TR got randomized by a small bit above
	tmp			<- df.ind[, which(DIAG_T<TIME_TR)]
	set(df.ind, tmp, 'DIAG_T', df.ind[tmp, TIME_TR+(TIME_TR-DIAG_T)])
	tmp			<- df.ind[, which(ART1_T<DIAG_T)]
	set(df.ind, tmp, 'ART1_T', df.ind[tmp, DIAG_T+(DIAG_T-ART1_T)])
	stopifnot( df.ind[, all(ART1_T<VLS1_TS, na.rm=TRUE)] )	
	tmp			<- df.ind[, which(DOD<TIME_TR)]
	set(df.ind, tmp, 'DOD', df.ind[tmp, TIME_TR+(TIME_TR-DOD)])
	#	add if transmission in the last 6 mo	
	df.ind[, RECENT_TR:=df.ind[, factor(as.numeric((DIAG_T-TIME_TR)<.5), levels=c(0,1), labels=c('N','Y'))]]	
	#	add time of infection of transmitter to df.trm	
	tmp		<- subset(df.ind, select=c(IDPOP, TIME_TR))
	setnames(tmp, c('IDPOP','TIME_TR'), c('IDTR','IDTR_TIME_INFECTED') )
	setkey(tmp, IDTR)
	df.trm	<- merge(df.trm, unique(tmp), by='IDTR', all.x=TRUE)
	stopifnot( df.trm[, !any(TIME_TR<=IDTR_TIME_INFECTED, na.rm=TRUE)] )
	#	simulate time individual ready for sequencing
	df.ind	<- PANGEA.Seqsampler.SimulateGuideToSamplingTimes.v2(df.ind, seqtime.mode= pipeline.args['seqtime.mode',][,v])	
	cat(paste('\nFound % sampled at or after ART start=', subset(df.ind, !is.na(DIAG_T))[, mean(!is.na(ART1_T) & T1_SEQ>=ART1_T, na.rm=TRUE)] ))
	cat(paste('\nFound % sampled after end of first viral suppression=', subset(df.ind, !is.na(DIAG_T))[, mean(!is.na(VLS1_TE) & T1_SEQ>=VLS1_TE, na.rm=TRUE)] ))
	# 
	#
	#	reduce to time frame of interest
	#
	#tmp		<- subset( df.trm, TIME_TR>=as.numeric( pipeline.args['yr.end',][, as.numeric(v)] ) )[, IDREC]
	df.trm	<- subset( df.trm, TIME_TR<as.numeric( pipeline.args['yr.end',][, as.numeric(v)] ) )
	#df.ind	<- subset(df.ind, !IDPOP%in%tmp)
	df.ind	<- subset(df.ind, DOB<pipeline.args['yr.end',][, as.numeric(v)] )
	df.ind	<- subset(df.ind, is.na(DOD) | DOD >= floor(min(df.trm$TIME_TR)) )
	cat(paste('\nFound individuals born before',pipeline.args['yr.end',][, as.numeric(v)],', n=', nrow(df.ind)))
	cat(paste('\nFound transmissions before',pipeline.args['yr.end',][, as.numeric(v)],', n=', nrow(df.trm)))
	cat(paste('\nTotal transmitters in sampling frame, n=', df.trm[, length(unique(IDTR))]))		
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
	#tmp		<- df.trm[, which(is.na(IDTR_TIME_INFECTED))]	
	#set( df.trm, tmp, 'IDTR_TIME_INFECTED', df.trm[tmp, runif(length(tmp), TIME_TR-5, TIME_TR)] )	
	tmp		<- PANGEA.ImportSimulator.SimulateStartingTimeOfIndexCase.v2(df.ind, df.trm, index.starttime.mode= pipeline.args['index.starttime.mode',][,v])
	df.trm	<- tmp$df.trm
	df.ind	<- tmp$df.ind	
	#
	#	sample sequences 
	#
	PANGEA.Seqsampler.v4(df.ind, df.trm, pipeline.args, outfile.ind, outfile.trm, with.plot=with.plot)
	#
	return(1)
}
######################################################################################
#	return GAG POL ENV ancestral sequences from BEAST PARSER output	
#	olli originally written 13-08-2014
#	tree 		beast trees in ape format, needed to compute calendar time for each ancestral sequence
#	node.stat	data.table containing meta information in nexus file for nodes
#	return 		list of GAG POL ENV sequences in ape format 
PANGEA.RootSeqSim.get.ancestral.seq<- function(tree, node.stat, tree.id.sep='_', tree.id.idx.mcmcit=2, tree.id.burnin=1, label.sep='|', label.idx.ctime=5)
{
	require(data.table)
	require(ape)
	
	tree.id				<- names(tree)
	#	add calendar time for inner node at NODE_ID to node.stat
	node.stat[, CALENDAR_TIME:=NA_real_]		
	setkey(node.stat, TREE_ID, NODE_ID)
	for(i in seq_along(tree.id))
	{
		cat(paste('\nProcess CALENDAR_TIME for tree id', tree.id[i]  ))
		label.ctime			<- sapply( strsplit(tree[[i]]$tip.label, label.sep, fixed=TRUE), '[[', label.idx.ctime ) 
		label.ctime			<- as.numeric(label.ctime)			
		depth				<- node.depth.edgelength( tree[[ i ]] )
		tmp					<- which.max(depth)
		depth				<- depth-depth[tmp]+label.ctime[tmp]
		tmp					<- node.stat[, which(TREE_ID==tree.id[i])]
		for(j in seq_along(depth))
		{
			set(node.stat, tmp[ node.stat[tmp, which(NODE_ID==j)] ], 'CALENDAR_TIME', depth[j])
		}								
	}
	tmp			<- node.stat[, length(which(is.na(CALENDAR_TIME)))]
	cat(paste('\nTotal node statistics with no CALENDAR_TIME [should be zero], n=', tmp  ))
	#	keep only inner nodes
	tmp			<- sapply(tree, Ntip)
	stopifnot(all(tmp==tmp[1]))
	node.stat	<- subset(node.stat, NODE_ID>tmp[1])	
	#
	set(node.stat, NULL, 'VALUE', node.stat[, gsub('\"','',VALUE)])
	#	checks of ancseq before we proceed
	tmp			<- node.stat[, list(NSEQ= nchar(VALUE)), by=c('TREE_ID', 'NODE_ID', 'STAT')]		
	stopifnot( tmp[, list(CHECK= all(NSEQ==NSEQ[1])), by='STAT'][, all(CHECK)] )
	set(tmp, NULL, 'GENE', tmp[, sapply(strsplit(STAT,'\\.'),'[[',1)])
	set(tmp, NULL, 'CODON_POS', tmp[, sapply(strsplit(STAT,'\\.'),'[[',2)])	
	tmp			<- dcast.data.table(tmp, TREE_ID + NODE_ID + GENE ~ CODON_POS, value.var='NSEQ')
	tmp			<- tmp[, list(CPM=min(CP1, CP2, CP3)), by=c('TREE_ID','NODE_ID','GENE')]
	stopifnot( tmp[, list(CHECK=all(CPM==CPM[1])), by='GENE'][, all(CHECK)] )
	setkey(tmp, GENE)
	#	truncate to following size of coding regions (if necessary)
	tmp			<- unique(tmp)[, list(STAT=paste(GENE,'.CP',1:3,sep=''), CPM=CPM), by='GENE']
	node.stat	<- merge(node.stat, subset(tmp, select=c(STAT, CPM)), by='STAT')
	set(node.stat, NULL, 'VALUE', node.stat[, substr(VALUE,1,CPM)])
	set(node.stat, NULL, 'CPM', NULL)
	#	reconstruct genes from codon positions
	node.stat	<- dcast.data.table(node.stat, TREE_ID + NODE_ID + CALENDAR_TIME ~ STAT, value.var="VALUE")
	node.stat	<- node.stat[, {
				tmp		<- do.call('rbind',sapply(list(ENV.CP1,ENV.CP2,ENV.CP3), strsplit, ''))
				env		<- paste(as.vector(tmp), collapse='')
				tmp		<- do.call('rbind',sapply(list(GAG.CP1,GAG.CP2,GAG.CP3), strsplit, ''))
				gag		<- paste(as.vector(tmp), collapse='')
				tmp		<- do.call('rbind',sapply(list(POL.CP1,POL.CP2,POL.CP3), strsplit, ''))
				pol		<- paste(as.vector(tmp), collapse='')
				list(GAG=gag, POL=pol, ENV=env, CALENDAR_TIME=CALENDAR_TIME)
			}, by=c('TREE_ID','NODE_ID')]
	
	set(node.stat, NULL, 'LABEL', node.stat[, paste(TREE_ID, NODE_ID, round(CALENDAR_TIME,d=3), sep='|')])		
	#	remove tree id STATE_xx where xx is smaller than burn-in
	set(node.stat, NULL, 'BEAST_MCMC_IT', node.stat[, as.integer(sapply(strsplit(TREE_ID,tree.id.sep),'[[',tree.id.idx.mcmcit))])
	node.stat			<- subset(node.stat, BEAST_MCMC_IT>tree.id.burnin)
	cat(paste('\nFound ancestral sequences, n=', nrow(node.stat)  ))
	node.stat
}