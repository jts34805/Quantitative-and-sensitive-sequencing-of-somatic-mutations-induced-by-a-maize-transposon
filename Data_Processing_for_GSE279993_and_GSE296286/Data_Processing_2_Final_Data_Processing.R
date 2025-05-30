MergeFn = function(A, MoleculeCounting = T) {  # Function to merge CleanData list into a matrix of Mu TIR counts
	if (MoleculeCounting) {
		MergingData = lapply(A, function(A2) { table(paste(A2$Chromosome, A2$Position, A2$Direction)) })
	} else {
		MergingData = lapply(A, function(A2) { by(A2$Reads, paste(A2$Chromosome, A2$Position, A2$Direction), sum) })
	}
	
	allLocations = unique(unlist(lapply(MergingData, names)))
	XX = matrix(unlist(lapply(MergingData, function(x) { x[allLocations] })), ncol = length(MergingData))
	rownames(XX) = allLocations
	colnames(XX) = names(MergingData)
	XX[is.na(XX)] = 0
	return(XX)
}

MuCountFn = function(A) {  # Function to connect reads out of both transposon direction, converting Mu TIR counts to Mu insertion counts
	gLocs = t(matrix(unlist(strsplit(rownames(A), ' ')), nrow=3))
	gLocs = data.frame(Chr = gLocs[,1], Pos = as.numeric(gLocs[,2]), Dir = gLocs[,3])
	gLocs$MuCenter = paste(gLocs[,1], gLocs[,2] + 4*as.numeric(gLocs[,3]), sep=',')
	
	MuCountsRB = matrix(0, nrow = length(unique(gLocs$MuCenter)), ncol = ncol(A))
	rownames(MuCountsRB) = unique(gLocs$MuCenter)
	colnames(MuCountsRB) = colnames(A)
	MuCountsLB = MuCountsRB
	
	MuCountsRB[gLocs$MuCenter[gLocs$Dir == '1'],] = A[gLocs$Dir == '1',]
	MuCountsLB[gLocs$MuCenter[gLocs$Dir == '-1'],] = A[gLocs$Dir == '-1',]


	# By default, TIR borders on both sides of a Mu element are connected assuming they will map 9 bp apart -- this is true for non-reference elements given the 9 bp target site duplication; the following code connects both sides of an element for insertions that are present in the reference genome
	MU = read.table('blastW22.out', sep = '\t', comment ='#')
	MU = MU[MU[,7] == 1,]
	MU$ID = paste(MU[,2], MU[,9] + sign(MU[,9] - MU[,10]), sign(MU[,9] - MU[,10]), sep = ',')
	MU = MU[order(MU[,2], MU[,9]),-1]
	MU$ID2 = paste(MU[,1], MU[,8] + 5*sign(MU[,8] - MU[,9]), sep = ',')
	
	MU2 = MU[MU[,7] >= 65,]
	MU2$Pair = NA
	PR = which(abs(diff(c(0, MU2[, 8]))) < 10000)
	MU2$Pair[PR] = MU2$ID2[PR - 1]
	
	MU2 = rbind(MU2[1,],MU2)
	MU2$ID2[1] = 'chr3,122342449'
	MU2$Pair[1] = 'chr3,122342376'

	swaps = !is.na(MU2$Pair) & (MU2$ID2 %in% rownames(MuCountsRB)) & (MU2$Pair %in% rownames(MuCountsRB))
	MuCountsRB[MU2[swaps,]$Pair,] = MuCountsRB[MU2[swaps,]$ID2,]
	MuCountsRB = MuCountsRB[!(rownames(MuCountsRB) %in% MU2[swaps,]$ID2),]
	rownames(MuCountsRB)[match(MU2[swaps,]$Pair, rownames(MuCountsRB))] = MU2[swaps,]$ID2
	MuCountsLB = MuCountsLB[!(rownames(MuCountsLB) %in% MU2[swaps,]$ID2),]

	gLocs2 = data.frame(Chr = sub(',.+', '', rownames(MuCountsRB)), PosL = as.numeric(sub('.+,', '', rownames(MuCountsLB))) + 4, PosR = as.numeric(sub('.+,', '', rownames(MuCountsRB))) - 4)

	MuCountsAll = MuCountsLB + MuCountsRB
	rownames(MuCountsAll) = gsub(' ', '', paste(apply(gLocs2[,1:2],1,paste,collapse=':'),gLocs2[,3],sep=','))
	rownames(MuCountsLB) = gsub(' ', '', apply(gLocs2[,1:2],1,paste,collapse=':'))
	rownames(MuCountsRB) = gsub(' ', '', apply(gLocs2[,c(1,3)],1,paste,collapse=':'))

	return(list(MuCountsAll = MuCountsAll, MuCountsRB = MuCountsRB, MuCountsLB = MuCountsLB))
}


MuStagger = function(X, thresh = 50, distthresh = 15) {
	Locs = t(matrix(unlist(strsplit(sub(':',',',rownames(X[[1]])),',')), nrow = 3))
	Locs = data.frame(chr = Locs[,1], position = as.numeric(Locs[,2]), position2 = as.numeric(Locs[,3]))
	ords = order(Locs[,1], Locs[,2], Locs[,3])
	X[[1]] = X[[1]][ords,]
	X[[2]] = X[[2]][ords,]
	X[[3]] = X[[3]][ords,]
	Locs = Locs[ords,]

	mergedists = list()
	merges = NULL
	for (i in 1:ncol(X[[1]])) {
		flagRBi = which((X[[2]][,i] > X[[3]][,i]*thresh) & (X[[2]][,i] >= thresh))
		flagLBi = which((X[[2]][,i]*thresh < X[[3]][,i]) & (X[[3]][,i] >= thresh))
		xi = Locs[flagLBi,2]
		checks = sapply(Locs[flagRBi,2], function(xx) { xx - xi })
		xi = X[[3]][,i][flagLBi]
		sums = sapply(X[[2]][,i][flagRBi], function(xx) { xx + xi })
		
		
		if (is.null(dim(checks))) { mergedists[[colnames(X[[1]])[i]]] = integer(0) } else {
		sums[abs(checks) > distthresh] = NA
		checks[abs(checks) > distthresh] = NA

		if (all(is.na(checks))) { mergedists[[colnames(X[[1]])[i]]] = integer(0) } else {

		# when there are multiple LR-border pairs that qualify, match the two with the highest Mu counts in a sample
		sumsMax = sweep(!is.na(sums), 1, apply(sums,1,max,na.rm=T), '*')
		sumsMax = sweep(sumsMax > 0, 2, apply(sumsMax,2,max,na.rm=T), '*')
		sumsMax = sumsMax*(sums > 0)
		checks[sums < sumsMax] = NA
		
		mergedists[[colnames(X[[1]])[i]]] = table(checks)
		
		flagRBi = flagRBi[colSums(!is.na(checks)) > 0]
		checks = checks[,colSums(!is.na(checks)) > 0]
		if (is.null(dim(checks))) {
			flagLRi = cbind(flagRBi, flagLBi[which(!is.na(checks))],i)
		} else {
			flagLRi = cbind(flagRBi, flagLBi[apply(!is.na(checks), 2, which)],i)
		}

		merges = rbind(merges,flagLRi)
	}}}
	merges = as.matrix(merges)
	
	merges2 = unique(merges[,1:2])
	merges2 = cbind(merges2, as.numeric(table(merges[,2])[as.character(merges2[,2])]))
	
	mergesName = cbind(rownames(X[[2]][merges2[,1],]), rownames(X[[3]][merges2[,2],]))
	for (dups in unique(mergesName[duplicated(mergesName[,1]),1])) {
		dups = which(mergesName[,1] == dups)
		swaps = which(merges2[dups,3] == max(merges2[dups,3]))[1]
		mergesName[dups[-swaps],2] = mergesName[dups[swaps],2]
	}
	mergesName = cbind(paste(mergesName[,2], sub('ch.+:','',mergesName[,1]),sep=','), mergesName)
	
	for (i in nrow(merges2)) {
		X[[1]][merges2[i,1],] = X[[1]][merges2[i,1],] + X[[1]][merges2[i,2],]
		X[[2]][merges2[i,1],] = X[[2]][merges2[i,1],] + X[[2]][merges2[i,2],]
		X[[3]][merges2[i,1],] = X[[3]][merges2[i,1],] + X[[3]][merges2[i,2],]
	}
	rownames(X[[1]])[merges2[!duplicated(merges2[,1]),1]] = mergesName[!duplicated(merges2[,1]),1]
	rownames(X[[2]])[merges2[!duplicated(merges2[,1]),1]] = mergesName[!duplicated(merges2[,1]),2]
	rownames(X[[3]])[merges2[!duplicated(merges2[,1]),1]] = mergesName[!duplicated(merges2[,1]),3]
	X[[1]] = X[[1]][-merges2[,2],]
	X[[2]] = X[[2]][-merges2[,2],]
	X[[3]] = X[[3]][-merges2[,2],]

	print(table(unlist(lapply(mergedists, length))))
	bp = barplot(rowSums(matrix(unlist(lapply(mergedists, function(xx) { xx[as.character(-distthresh:distthresh)] })), nrow=distthresh*2+1), na.rm=T), border = NA, col = 'dark gray')
	axis(1, at = bp[seq(1,distthresh*2+1,5),1], labels = seq(-distthresh,distthresh,5))
	return(X)
}


adjCounts = function(X, thresh = 2, pseudo = 20) {
	filter = abs(log(rowSums(X[[2]]) + pseudo, 2) - log(rowSums(X[[3]]) + pseudo, 2)) > thresh
	cat(sum(filter), 'insertions estimated from only one border', sum(filter)/length(filter), '% of total\n')
	X[[1]][filter,] = abs(X[[2]][filter,] - X[[3]][filter,]) + (X[[2]][filter,] + X[[3]][filter,])
	return(X)
}


######################

load('CleanData.rda')
ReadsPerUMI = unlist(lapply(CleanData, function(x) { mean(x[,4]) }))
summary(ReadsPerUMI)

CleanData2 = lapply(CleanData, function(xx) { xx[xx$Reads >= mean(xx$Reads)/5,] })

MoleculeCounts = MergeFn(CleanData2)  # Merge CleanData list into a matrix of Mu TIR counts

MuCounts2 = MuCountFn(MoleculeCounts)  # Connect both TIRs and convert to matrix of Mu Insertion Counts
MuCounts3 = MuStagger(MuCounts2)  # Connect TIRs with shift from 9 bp TSD... this primarily affects inherited counts for small number of paternal insertions
nrow(MuCounts3[[1]]) - nrow(MuCounts2[[1]])  # 273 insertion sites were affected by the 'MuStagger' function (out of ~3.8 million)
MuCounts = adjCounts(MuCounts3, thresh = 1, pseudo = 500)  # For Mu insertions where there was a large discrepancy between the left and right borders, estimate MuCounts based on only the better-sampled border. This affects 422 insertions (0.00013% of total)
rm(MoleculeCounts, MuCounts2, MuCounts3)

keeps = rowSums(MuCounts[[1]]) > 0
MuCounts = lapply(MuCounts, function(xx) { xx[keeps,] })



## For GSE279993 data, output the following:
write.csv(MuCounts[[1]], file = 'MuSeq2 data for maize tissues.csv')
write.csv(MuCounts[[2]], file = 'MuSeq2 data for maize tissues _ Right Border.csv')
write.csv(MuCounts[[3]], file = 'MuSeq2 data for maize tissues _ Left Border.csv')



## For GSE296286 data, output the following:
