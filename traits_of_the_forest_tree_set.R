
library(ape)

RcodePath <- "D:/PhD/GitHub/Rshiny_webapp/Rshiny_webapp/"
source( paste(RcodePath,"getEl.R",sep="") )
source( paste(RcodePath,"quick_time_scale_tree.R",sep="") )
source( paste(RcodePath,"HPD_stats.R",sep="") )


select_n	<- function( selstate, traits=traits, n=n ) {
		jj <- which(traits==selstate)
		if (length(jj) < n) {
			return( jj )
		} else {
			return( sample(jj, n) )
		}
	}

select_random <- function( utraits=utraits, traits=traits, n=n) {
		res <- lapply(utraits, select_n, traits=traits, n=n)
		res <- sort(unlist(res))
		return( res )
	}


# traits are factors
select_in_category <- function( traits=traits, nreps=10, numPerCat=10, evens=FALSE ) {


		tbl	  	<- table(traits)
		utraits 	<- rownames(tbl)
		ntraits 	<- length(utraits)

		minPerCat 	<- min(tbl)
		if (minPerCat < numPerCat) {
			if (evens) {
				print( paste("Warning cannot do proper evens sampling because minPerCat=",minPerCat,"and numPerCat=",numPerCat) )
			}
		}

		if (evens) {
			# the same number of samples per category
			# split into non-overlapping sets
			# max of nreps (could be less)
			assign_to_set <- array(0, length(traits))
			for (j in 1:ntraits) {
				jj <- which(traits==utraits[j])

				if (length(jj) > 1) {
					# re-order
					jj <- sample(jj, length(jj))
						
					# first numPerCat to first bin etc
					k <- ceiling( (1:length(jj)) / numPerCat )
					assign_to_set[jj] <- k
				}
				if (length(jj) == 1) {
					assign_to_set[jj] <- 1
				}
			}
			temp 		<- table(assign_to_set)
			ok_sets 	<- which(temp == numPerCat*ntraits)
			n_ok_sets	<- length(ok_sets)

			if (n_ok_sets >= 1) {
				selinds	<- vector("list",n_ok_sets)
				for (i in 1:n_ok_sets) {
					selinds[[i]] <- which(assign_to_set == ok_sets[i] )
				}
				if (n_ok_sets < nreps) {
					print( paste("Warning could only make",n_ok_sets,"non-overlapping sets") )
				}	
	
			} else {
				print( paste("There are no ok non-overlapping sets, so just returning the first one" ) )
				selinds <- vector("list",1)
				selinds[[1]] <- which(assign_to_set == 1)
			}

		} else {
			# a maximum of numPerCat samples per category (could be less if there are less)
			# potentially overlapping sets
			selinds <- vector("list", nreps)
			for (i in 1:nreps) {
				selinds[[i]] <- select_random(utraits=utraits, traits=traits, n=numPerCat)
			}
		}

		return( selinds )

	}

get_sampled_sequence_list <- function( traits=traits, sequenceNames=sequenceNames, nreps=10, numPerCat=10, evens=FALSE, 
								includeOutgroup=FALSE, outgroupName="Outgroup" ) {

		if (includeOutgroup) {
			if (length(outgroupName) > 1) {
				og_ind <- unlist(lapply(outgroupName, grep, sequenceNames))
			} else {
				og_ind <- grep(outgroupName,sequenceNames)
			}

			if (length(og_ind) >= 1) {
				incs		 <- setdiff(1:length(traits), og_ind)
				no_og_traits <- as.matrix(traits[ incs ])
				no_og_seqnames<- as.matrix(sequenceNames[ incs ])
			} else {
				no_og_traits <- traits
				no_og_seqnames<- sequenceNames
			}
			selinds		<- select_in_category(traits=no_og_traits, nreps=nreps, numPerCat=numPerCat, evens=evens)

			selseqnames		<- vector("list",length(selinds))
			for (j in 1:length(selinds)) {
				selseqnames[[j]] <- c(as.matrix(sequenceNames[ og_ind ]), no_og_seqnames[ selinds[[j]] ] )
			}


		} else {
			selinds		<- select_in_category(traits=traits, nreps=nreps, numPerCat=numPerCat, evens=evens)
			selseqnames		<- vector("list",length(selinds))
			for (j in 1:length(selinds)) {
				selseqnames[[j]] <- sequenceNames[ selinds[[j]] ]
			}
		}

		return( selseqnames )
	}


deduplicate_sequences <- function( seqs, decDates=c(-1), includeDate=TRUE, maxsamples=-1, gtol=0, dtol=0.5/365 ) {

		ntaxa <- length(seqs)
		dd <- dist.dna( seqs, model="TN93", gamma=1, pairwise.deletion=TRUE, as.matrix=TRUE)
		GG <- (dd <= gtol)*1
		M  <- GG

		if ( includeDate & (decDates[1] == -1)) {
			taxa 		<- attributes(seqs)$names
			decDates 	<- as.numeric(apply(as.matrix(taxa), 1, getEl, ind=1, fromEnd=TRUE, sep="\\|"))
		}
		
		if (includeDate) {
			TT		<- matrix( rep(decDates, length(decDates)), length(decDates), length(decDates) )
			TT		<- (abs( TT - t(TT) ) <= dtol)*1
			M		<- ((GG + TT)==2)*1
		}

		g 		<- graph.adjacency(M)
		clusts 	<- clusters(g)
		nclusts	<- clusts$no
		maxcsize	<- max(clusts$csize)
		dups		<- which(clusts$csize > 1)
		ndups		<- length(dups)

		if (maxcsize==1) {
			selinds 	 <- vector("list",1)
			selinds[[1]] <- 1:ntaxa
		} else {
			if (maxsamples < 0) {
				maxsamples <- maxcsize
			} else {
				if (maxcsize <= maxsamples) {
					maxsamples <- maxcsize
				}
			}
			assign_to_set <- matrix(0, nclusts, maxsamples)
			inds		  <- which(clusts$csize==1)
			singles	  <- match(inds,clusts$membership)
			assign_to_set[inds,] <- singles

			# check
			# sum(M[singles,singles])==length(singles)

			inds		<- which(clusts$csize > 1)
			for (j in 1:length(inds)) {
				multiples <- which(clusts$membership==inds[j])
				if (length(multiples) >= maxsamples) {
					assign_to_set[inds[j],] <- sample(multiples,maxsamples)
				} else {
					assign_to_set[inds[j],1:length(multiples)] <- multiples
					nremain <- maxsamples-length(multiples)
					if (nremain > 1) {
						assign_to_set[inds[j],(length(multiples)+1):maxsamples] <- sample(multiples,nremain,replace=TRUE)
					} else {
						assign_to_set[inds[j],(length(multiples)+1)] <- sample(multiples,1)
					}
				}
				
			}

			selinds <- vector("list",maxsamples)
			for (k in 1:maxsamples) {
				#print( sum(M[assign_to_set[,k],assign_to_set[,k]]) == length(assign_to_set[,k]) )
				selinds[[k]] <- assign_to_set[,k]
			}

		}

		return( selinds )
		
		
	}

tree_of_the_forest <- function(path,seqs,name,traitFileName,jointTrait1,sampleType,numPerCat,repsPerCat,jn,jr,outgroups = "Outgroup",removeOutgroup = TRUE ,lowerW = 5e-4,upperW = 1.5e-2,use.expect = TRUE,expectW = 2e-3)
{
  incProgress(0.1, detail = paste("Doing part", 2))
  numPerCat = numPerCat
 # print(numPerCat)
  print("avant_read_dna")
	#seqs 		<- read.dna( paste(path,paste(name,".fasta",sep=""),sep="/"), format="fasta", as.matrix=FALSE)
	taxa 		<- as.matrix(attributes(seqs)$names)

# LOAD TRAITS
# logging
	treeFileName<- paste(paste(path,name,sep="/"),"_nj_tree_set.trees.txt",sep="")
	logName	<- paste(paste(path,name,sep="/"),"_traits_of_the_forest_tree_set_log.txt",sep="")
	write( paste("Path =",path), file=logName, append=FALSE)
	write( paste("Name =",name), file=logName, append=TRUE)
	write( paste("Newick trees = ",name,"nj_tree_set.txt",sep=""), file=logName, append=TRUE)
	write( paste("Time start=",Sys.time()), file=logName, append=TRUE)
	write( "Getting traits", file=logName, append=TRUE)

	traitsTbl	<- read.table( paste(paste(path,traitFileName,sep="/"),".txt",sep=""), header=TRUE, sep="\t")
	cn		<- colnames(traitsTbl)
	if (jointTrait1 != "all") {
		c1		<- which(cn==jointTrait1)
#		c2		<- which(cn==jointTrait2)
	}

	# check
	all( taxa == traitsTbl[,1] )
	
# FULL JOINT SAMPLING FOR DIVERSITY
# logging

	write( "-------------------", file=logName, append=TRUE)
	write( paste("Time =",Sys.time()), file=logName, append=TRUE)
	write( paste("** Full Joint Sampling for Diversity with ",jointTrait1,"+",jointTrait2), file=logName, append=TRUE)
	write( paste("numPerCat=",jn), file=logName, append=TRUE)
	write( paste("nreps=",jr), file=logName, append=TRUE)	

	incProgress(0.1, detail = paste("Doing part", 3))
	if (jointTrait1 != "all") {
		joint_traits		<- paste( traitsTbl[,c1] )# , traitsTbl[,c2] )
	} else {
		joint_traits <- array("",length(traitsTbl[,1]))
		for (j in 2:length(cn)) {
#		  j=2
			joint_traits <- paste(joint_traits, traitsTbl[,j])
		}
	}
	jointList			<- get_sampled_sequence_list(traits=joint_traits, sequenceNames=taxa, numPerCat=jn, nreps=jr, 
						includeOutgroup=TRUE, outgroupName=outgroups)
	incProgress(0.1, detail = paste("Doing part", 4))
	print("After jointlist creation")
	print(jointList)
	#print( length(traitsTbl[,1]) )
	#print( table(unlist(lapply(jointList,length))) )

# SAMPLING BY A TRAIT

# logging
	write( "-------------------", file=logName, append=TRUE)
	write( paste("Time =",Sys.time()), file=logName, append=TRUE)
	write( paste("** Sampling by trait=",sampleType), file=logName, append=TRUE)

	print("before table print")
	print(cn)
	print(sampleType)
	print(traitsTbl[, which(cn==sampleType) ])
	print( table( traitsTbl[, which(cn==sampleType) ] ) )
	seqList 	<- c()

	for (kk in 1:length(jointList)) {
		jinds 		<- match(jointList[[kk]], traitsTbl[,1])
		jinds
		traits  		<- traitsTbl[jinds,which(cn==sampleType)]
		sequenceNames	<- traitsTbl[jinds,1]

		seqList_kk 		<- get_sampled_sequence_list(traits=traits, sequenceNames=sequenceNames, 
										nreps=repsPerCat, evens=FALSE, numPerCat=numPerCat,
										includeOutgroup=FALSE, outgroupName="Outgroup")
		if (kk==1) {
			seqList <- seqList_kk
		} else {
			for (k in 1:length(seqList_kk)) {
				seqList[[ length(seqList)+1 ]] <- seqList_kk[[k]]
			}
		}
	}

	# save list of chosen sequences in R format
	save(seqList, file=paste(paste(path,name,sep="/"),"_",sampleType,"_seqList.Rdata",sep=""))

	print( table(unlist(lapply(seqList, length))) )

# TIME SCALED TREES ONLY
	nreps <- length(seqList)
	taxa  <- traitsTbl[,1]

# logging
	write( "-------------------", file=logName, append=TRUE)
	write( paste("Time =",Sys.time()), file=logName, append=TRUE)
	write( "** Time Scaled Trees Only", file=logName, append=TRUE)
	write( paste("nreps=",nreps), file=logName, append=TRUE)
	incProgress(0.1, detail = paste("Doing part", 5))
	traits 	<- traitsTbl[,which(cn==sampleType)]
	full_tbl	<- table(traits)
	ustate	<- sort(unique( traits ))
	nstates	<- length(ustate)
	j		<- 1
	sinds 	<- match(seqList[[j]], taxa)
	all( seqList[[j]] == traitsTbl[sinds,1] )
	sub_tbl <-  table( traits[sinds] )
	print(sub_tbl)
	write("Trait distribution in full data", file=logName, append=TRUE)
	write.table( full_tbl, file=logName, col.names=FALSE, row.names=TRUE, sep="\t", append=TRUE)
	write("Trait distribution per subsample", file=logName, append=TRUE)
	#print("ici")
	write.table( sub_tbl, file=logName, col.names=FALSE, row.names=TRUE, sep="\t", append=TRUE)

	trStats <- matrix(0, nreps, 3)
	colnames(trStats) <- c("w","TreeLength","RootHeight")

# logging
	write( "-------------------", file=logName, append=TRUE)
	trs    <- vector("list",nreps)
	rn_trs <- vector("list",nreps)

	for (j in 1:nreps) {
#	  j=2
		# select sequences
		sinds 	<- match(seqList[[j]], taxa)
		selseqs	<- seqs[sinds]
	  
	  #selseqs	<- seqs[as.numeric(seqList[[j]])]
#print("ici")
		tr 		<- quick_time_scale_tree_from_sequences( selseqs, lowerW=lowerW, upperW=upperW, use.expect=use.expect, expectW=expectW,
											guessRoot=FALSE, useOutgroup=FALSE, outgroupName=outgroups, removeOutgroup=removeOutgroup )
		ntips 	<- length(tr$tip.label)

		# record tree results
		trStats[j,]			<- c(tr$w, sum(tr$edge.length), tr$nodeTimes[ntips+1])

		# match traits to tree tips
		tinds		<- match(tr$tip.label, taxa)

		tr$tip.state		<- traits[tinds]

		trs[[j]] <- tr

		# for tree set; want all the taxa names to be the same between subsamples
		# rename for traits
		rn_tr <- tr
		
		for (k in 1:nstates) {
			kinds 			<- which(rn_tr$tip.state==ustate[k])
			newNames 			<- paste(1:length(kinds),ustate[k],sep="_")
			rn_tr$tip.label[kinds] 	<- newNames
		}

		rn_trs[[j]] <- rn_tr
#print("ici")
# logging
		if ((j %% 100)==0) {
		  #print("ici")
		  #print(j)
			write( paste("Time =",Sys.time(),"completed",j,"trees and models"), file=logName, append=TRUE)
		}
	}
	incProgress(0.1, detail = paste("Doing part", 6))
# logging
	write( paste("Time =",Sys.time(),"completed",j,"trees and models"), file=logName, append=TRUE)
		
	outname <- paste(paste(path,name,sep="/"),"_",sampleType,"_trStats.txt",sep="")
	write.table( trStats, file=outname, col.names=TRUE, row.names=FALSE, sep="\t", append=FALSE )

	# empirical trees file
	rn_trs[2]
	write.nexus(rn_trs, translate=TRUE, file=treeFileName)
	temp <- readLines( treeFileName )
	tinds<- grep("TREE * UNTITLED",temp,fixed=TRUE)
	for (j in 1:length(tinds)) {
		temp[tinds[j]] <- gsub("TREE * UNTITLED",paste("tree STATE_",j,sep=""),temp[tinds[j]],fixed=TRUE)
	}
	write(temp, file=treeFileName)
  list_trees <- temp
	save(trs, file=paste(paste(path,name,sep="/"),"_orignames.trs.Rdata",sep=""))
	save(rn_trs, file=paste(paste(path,name,sep="/"),"_newnames.rn_trs.Rdata",sep="") )

	# use the first replicate as a template for BEAUTI upload
	# traits file
	#	  <- 1
	tempTbl <- cbind(rn_trs[[j]]$tip.label, paste(rn_trs[[j]]$tip.state))
	reinds  <- sort(paste(tempTbl[,1]), index.return=TRUE)$ix
	tempTbl <- tempTbl[reinds,]
	colnames(tempTbl) <- c("traits",sampleType)
	write.table(tempTbl, file=paste(paste(path,name,sep="/"),"_nj_tree_set_traits.txt",sep=""), sep="\t", col.names=TRUE, row.names=FALSE,quote = FALSE)
	# remove the "
	temp <- readLines( paste(paste(path,name,sep="/"),"_nj_tree_set_traits.txt",sep="") )
	temp <- gsub("\"", "", temp)
	traittaible2 <- matrix(gsub(".*\t", " ", temp))
	traittaible1 <- matrix(gsub("\t.*", " ", temp))
	traittaible <- cbind(traittaible1,traittaible2)
	traittaible <- traittaible[-1,]
	colnames(traittaible) <- c("taxa","trait")
	print(traittaible)
	write(temp, file=paste(paste(path,name,sep="/"),"_nj_tree_set_traits.txt",sep=""))

	# sequences file
	sinds 	<- match(trs[[j]]$tip.label, taxa)
	selseqs	<- seqs[sinds]
	attributes(selseqs)$names <- rn_trs[[j]]$tip.label
	selseqs	<- selseqs[reinds]
	write.dna( selseqs, file=paste(paste(path,name,sep="/"),"_nj_tree_set.fas",sep=""), format="fasta", nbcol=-1, colsep="")
  #return(selseqs)
	print("Avant ecrire table ")
	print(traittaible)
	return(list(first=selseqs, second=traittaible,third = list_trees))
	

# logging
	write( "-------------------", file=logName, append=TRUE)
	write( paste("End Time =",Sys.time()), file=logName, append=TRUE)

}


