# functions to get a quick approx time scaled tree
# S J Lycett
# 12 April 2016

# load ape
library(ape)

#################################################################################################
# INTERNAL FUNCTIONS
# do not modify
#################################################################################################


#################################################################################################
# from getEl.R
# useful function to split a sequence name
# S. J. Lycett
# 19 May 2011
# 6  Oct 2011

getEl	<- function( line, sep=",", ind=-1, final=FALSE, reconstruct=FALSE, ex=-1, fromEnd=FALSE ) {
	els	<- strsplit(line, sep)[[1]]

	if (ind[1] != -1) {
		if (fromEnd) {
			ind <- length(els)-(ind-1)
		}
	}

	if (final) {
		return( els[length(els)] )
	} else {

		if (reconstruct) {
			if (ex[1] > 0) {
				if (fromEnd) {
					ex <- length(els)-(ex-1)
				}
				ind <- setdiff((1:length(els)),ex)
			}

			newLine <- els[ind[1]]
			if (length(ind) > 1) {
				for (i in 2:length(ind)) {
					newLine <- paste(newLine, els[ind[i]], sep=sep)
				}
			}
			return ( newLine )
		} else {
			if ( ind[1] == -1 ) {
				return( els )
			} else {
				return( els[ind] )
			}
		}
	}
}

########################################################################################################
# from nodeDists.R (12 Sept 2014)

# function to calculate the distance from root in a tree
# (not nec time scaled tree)

distFromRoot	<- function( tr ) {

	rootNode	<- length(tr$tip.label)+1
	nodeDists	<- array(0, max(tr$edge))
	
	toProcess	<- c(rootNode)

	while ( length(toProcess) > 0 ) {
		einds 			<- which(tr$edge[,1]==toProcess[1])

		if (length(einds) > 0) {
			children 			<- tr$edge[einds,2]
			nodeDists[children] 	<- nodeDists[toProcess[1]] + tr$edge.length[einds]

			toProcess			<- c(toProcess, children)
		}
		toProcess			<- setdiff(toProcess, toProcess[1])
	}

	tr$nodeDists <- nodeDists
	tr$rootHeight<- max(tr$nodeDists) 	# 12 Sept 2014

	return ( tr )
}


# use the distance from root to add times to all internal nodes
# note assumes a time scaled tree

nodeTimes	<- function(tr, youngestTip=2017.164) {

	if ( !any(attributes(tr)$names == "nodeDists") ) {
		tr <- distFromRoot(tr)
	}

	 	# 12 Sept 2014 - included in distFromRoot
	 	# tr$rootHeight<- max(tr$nodeDists)

	tr$nodeTimes <- youngestTip - tr$rootHeight + tr$nodeDists

	return ( tr )
}


########################################################################################################

# tip.dates:	the tip dates for the tr$tip.label (in order) as a numeric array
# if not set, then tip dates are taken from the tr$tip.label assuming that they are in numeric format in the last element of the name
# e.g. 10664|H5N1|China|A/little_egret/Hong_Kong/718/2006|2006.500
# 
# sep:  		optional - the separator for the processing the date from the taxa name (see tip.dates)
#
# guessRoot:	if true then root of tree is guessed as the ancestor of the oldest set of sequences
# oldest_tol:	used in guess root
# 
quick_time_scale_tree	<- function( tr=tr, tip.dates=c(-1), sep="/", 
							guessRoot=TRUE, oldest_tol=0, 
							useOutgroup=FALSE, outgroupName="Outgroup", removeOutgroup=FALSE,
							doPlot=FALSE, lowerW=0, upperW=1, use.expect=FALSE, expectW=5e-3 ) {

	if (tip.dates[1] == -1) {
		taxa 		<- tr$tip.label
		tr$tip.date <- as.numeric(apply(as.matrix(taxa), 1, getEl, ind=1, fromEnd=TRUE, sep=sep))
	} else {
		tr$tip.date <- tip.dates
	}

	if (useOutgroup) {
		if (length(outgroupName) > 1) {
			og <- unlist(lapply(outgroupName, grep, sequenceNames))
		} else {
			og <- grep(outgroupName,sequenceNames)
		}

		if (length(og) >= 1) {
			
			if (length(og) > 1) {
				anc_of_og	<- getMRCA(tr, og)
				tr		<- root(tr, node=anc_of_og)
			} else {
				tr	<- root(tr, og)
			}	
		

			if (removeOutgroup) {
				tr	<- drop.tip(tr, tr$tip.label[og])
			}

			if (tip.dates[1] == -1) {
				taxa 		<- tr$tip.label
				tr$tip.date <- as.numeric(apply(as.matrix(taxa), 1, getEl, ind=1, fromEnd=TRUE, sep=sep))
			} else {
				tr$tip.date <- tip.dates
			}

		} else {
			print("Warning rooting by outgroup has not worked because cant identify outgroup from tip label")
		}

	} else {
		if (guessRoot) {
			oldest	<- min(tr$tip.date)
			oldest_tips <- which( abs(tr$tip.date-oldest) <= oldest_tol)
			if (length(oldest_tips)==1) {	
				anc_of_oldest	<- tr$edge[tr$edge[,2]==oldest_tips,1]
			} else {
				anc_of_oldest	<- getMRCA(tr, oldest_tips)
			}
			tr			<- root(tr, node=anc_of_oldest)
		}
	}


	tr 	<- distFromRoot(tr)
	ntips <- length(tr$tip.label)
	nnodes<- max(tr$edge)

	# date ~ w*dist
	x	<- tr$tip.date
	y	<- tr$nodeDists[1:ntips]
	fit 	<- lm(y ~ x)

	w	<- fit$coefficients[2]
	if (w < 0) {
		print( paste("Warning w=",w,"but setting to +ve") )
		w <- abs(w)
	}
	if (w < lowerW) {
		print( paste("Warning w=",w," which is < lower bounds",lowerW) )
		if (use.expect) {
			w <- expectW
		} else {
			w <- lowerW
		}
	}
	if (w > upperW) {
		print( paste("Warning w=",w," which is > upper bounds",upperW) )
		if (use.expect) {
			w <- expectW
		} else {
			w <- upperW
		}
	}
	newTr <- tr
	newTr$edge.length <- newTr$edge.length/w
	newTr <- distFromRoot(newTr)
	newTr <- nodeTimes(newTr, youngestTip=max(tr$tip.date))
	newTr$w <- w

	if (doPlot) {
		op <- par(mfrow=c(1,2))
		plot( tr$tip.date, tr$nodeDists[1:ntips] )
		abline(fit)

		plot(newTr, show.tip=FALSE)
		add.scale.bar(1)
		nodelabels(format(newTr$nodeTimes[(ntips+1):nnodes],digits=6))
		par(op)
	}

	return( newTr )
}

quick_time_scale_tree_from_sequences <- function( seqs=seqs, model="TN93", gamma=1, pairwise.deletion=TRUE, etol=1e-10,
									tip.dates=c(-1), sep="/", guessRoot=TRUE, oldest_tol=0, 
									useOutgroup=FALSE, outgroupName="Outgroup", removeOutgroup=FALSE,
									doPlot=FALSE,
									lowerW=0, upperW=1, use.expect=FALSE, expectW=5e-3 ) {

	dd <- dist.dna( seqs, model=model, gamma=gamma, pairwise.deletion=pairwise.deletion, as.matrix=TRUE )
	tr <- njs(dd)
	tr <- multi2di(tr)
	tr <- ladderize(tr)
	einds		<- which(tr$edge.length <= etol)
	if (length(einds) > 0) {
		tr$edge.length[einds] <- etol
	}

	tr <- quick_time_scale_tree( tr=tr, tip.dates=tip.dates, sep=sep, guessRoot=guessRoot, oldest_tol=oldest_tol, 
							useOutgroup=useOutgroup, outgroupName=outgroupName, removeOutgroup=removeOutgroup,
							doPlot=doPlot, lowerW=lowerW, upperW=upperW, use.expect=use.expect, expectW=expectW)
	

	tr <- multi2di(tr)
	tr <- ladderize(tr)
	einds		<- which(tr$edge.length <= etol)
	if (length(einds) > 0) {
		tr$edge.length[einds] <- etol
	}

	return( tr )
}





