#### Analysis code for Scherer, et al. (2025)  ####
#                                                 #
# This code will reproduce all main text figures  #
# from the manuscript, starting from the          #
# processed data tables available from GEO        #
# accessions GSE279993 and GSE296286.             #
#                                                 #
###################################################
#                                                 #
# Figures created using:                          #
# R version 4.1.2                                 #
# ggplot2 version 3.5.1                           #
# ggpubr version 0.6.0                            #
# ComplexHeatmap version 2.10.0                   #
#                                                 #
###################################################

######################

genomesize = 2133880228  # W22 mappable genome size, calculated using genmap

require(ggplot2)
require(ggpubr)
require(ComplexHeatmap)
MuCounts = read.csv('MuSeq2 data for maize tissues.csv', row.names=1)  # .csv file available at NCBI GEO as a supplementary file from accession GSE279993


CPM2 = MuCounts[,c('L74','L74T','E74')]


poolTechReps = function(X) {
		samps = sub('T$', '', colnames(X))
		repsamps = unique(samps[duplicated(samps)])
		X[,repsamps] = (X[,repsamps] + X[,paste(repsamps,'T',sep='')])
		X = X[,!(colnames(X) %in% paste(repsamps,'T',sep=''))]
		return(X)
	}

MuCounts = poolTechReps(MuCounts)
CPM = sweep(MuCounts, 2, colSums(MuCounts), '/')*10^6


## Remove genome locations that were identified at moderate level in >= 50% of Mu inactive controls or >= 50% of all samples; 18 genome locations were masked by this criteria
blacklist = which((rowMeans((CPM >= 10) & (CPM < 1000)) >= .5) | (rowMeans((CPM[,grepl('MI',colnames(CPM))] >= 10) & (CPM[,grepl('MI',colnames(CPM))] < 1000)) >= .5))
CPM = CPM[-blacklist,]
CPM2 = CPM2[-blacklist,]
CPM = sweep(CPM, 2, colSums(CPM), '/')*10^6  # renormalize after removing blacklisted sites
CPM2 = sweep(CPM2, 2, colSums(CPM2), '/')*10^6  # renormalize after removing blacklisted sites
MuCounts = MuCounts[-blacklist,]
length(blacklist)  # 18 blacklisted sites

historicals = rowMeans(CPM[,grepl('MI', colnames(CPM))] >= 1000) == 1
sum(historicals)  # 29 historical insertions


################# Fig 1 #################

VAFnormalize = function(X, thresh = 1000) {
	plant = sub('_.+', '', gsub('[a-z]', '', colnames(X), ignore.case=T))
	samps = unique(plant[grepl('E', colnames(X))])
	X = X[,plant %in% samps]  # select only sample with endosperm when calculating VAF
	plant = plant[plant %in% samps]
	paternals = list()
	sems = NULL
	for (s in samps) {
		X2 = X[,plant == s] >= thresh
		if (ncol(X2) == 2) {
			paternals[[s]] = (rowMeans(X2) == 1) & !historicals
			} else {
			paternals[[s]] = X2[,grepl('E', colnames(X2))] & (rowSums(X2[,!grepl('E', colnames(X2))]) >= 1) & !historicals
		}
		X[,plant == s] = sweep(X[,plant == s], 2, colMeans(X[paternals[[s]],plant == s]), '/')
		sems1 = apply(X[paternals[[s]],plant == s], 2, function(xx) { sd(xx)/sqrt(length(xx)) }) / colMeans(X[paternals[[s]],plant == s])
		names(sems1) = colnames(X)[plant == s]
		sems = c(sems, sems1)
	}
	return(list(VAF = X, paternals = paternals, normErr = sems))
}
	
plotPairs = function(i, j, mn = .0001) {
	x = VAF[[1]][!historicals,i]/2
	y = VAF[[1]][!historicals,j]/2
	plants = gsub('[a-z]','',colnames(VAF[[1]])[c(i,j)],ignore.case=T)
	highlights = which((VAF[[2]][[plants[1]]] & VAF[[2]][[plants[2]]])[!historicals])
	plot(unique(cbind(x, y)[-highlights,])+mn, pch = 19, log = 'xy', cex = .3, xaxt = 'n', yaxt = 'n', xlim = c(mn, 1+mn), ylim = c(mn, 1+mn), las = 1, xlab = colnames(VAF[[1]])[i], ylab = colnames(VAF[[1]])[j])
	points(cbind(x, y)[highlights,]+mn, col='red', pch = 19, cex = .3)
	
	axis(1, at = c(0, 10^(log(mn,10):0)) + mn, labels = c(0, 10^(log(mn,10):0)), lwd = 0, lwd.tick = 2, tck = -.03)
	axis(1, at = rep(2:9, -log(mn,10)+1)*10^rep(log(mn,10):0, each = 8) + mn, labels = F, lwd = 0, lwd.tick = 1, tck = -.02)
	axis(2, at = c(0, 10^(log(mn,10):0)) + mn, labels = c(0, 10^(log(mn,10):0)), lwd = 0, lwd.tick = 2, tck = -.03)
	axis(2, at = rep(2:9, -log(mn,10)+1)*10^rep(log(mn,10):0, each = 8) + mn, labels = F, lwd = 0, lwd.tick = 1, tck = -.02)
	
	cat('R2 = ', cor(x,y)^2, '\n')
	cat('R2 (paternal) = ', cor(x[highlights],y[highlights])^2, '\n')
	cat('R2 (de novo) = ', cor(x[-highlights],y[-highlights])^2, '\n')
	cat('N = ', length(highlights), 'zygotic insertion sites\n')
	cat('N = ', sum(rowSums(cbind(x,y)[-highlights,]) > 0), 'de novo insertion sites\n')
	cat(sum((x[-highlights] > 0) & (y[-highlights] > 0)), ' insertions (', round(100*sum((x > 0) & (y > 0)) / sum((x > 0) | (y > 0)),2), '%) present in both\n', sep = '')
		cat(sum(x[-highlights] > 0) * mean(y[-highlights] > 0), ' insertions (', round(100*sum((x > 0) & (y > 0)) / sum((x > 0) | (y > 0)),2), '%) present in both\n', sep = '')
}


MuSub = CPM[,grepl('74',colnames(CPM)) | grepl('75',colnames(CPM)) | grepl('MI',colnames(CPM))]
MuSub[, grepl('MI',colnames(MuSub))] = MuSub[,grepl('MI',colnames(MuSub))]/3
MuSub = MuSub[rowSums(MuSub >= 2500) >= 1,]  # 2500 is the 1% quantile for zygotic insertions
MuSub[MuSub >= 15000] = 15000
MuSub[MuSub <= 500] = 500

Fig1B = draw(Heatmap(log(MuSub[,c(paste(rep(c('L','R','P','E'), 2), rep(75:74, each = 4), sep = ''), 'MI1','MI2','MI3')], 10), show_row_names = F, col = colorRampPalette(c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))(100), cluster_columns = F))


VAF = VAFnormalize(CPM2)

plotPairs(1,2)  # Fig 1C; R2 = .997 for all insertions, = 0.955 for de novo insertions
plotPairs(1,3)  # Fig 1D; 6.7 * 10^-5 for de novo insertions



################# Normalize data and calculate quantifications for results section #################

VAFdenovo = function(X) {
	plant = sub('_.+', '', gsub('[a-z]', '', colnames(X[[1]]), ignore.case=T))
	for (pl in unique(plant)) {
		X[[1]][X[[2]][[pl]], plant == pl] = NA
	}
	return(X[[1]][!historicals,])
}

summary(colSums(MuCounts[,grepl('L',colnames(MuCounts)) & !grepl('_',colnames(MuCounts))]))  # 1,698,397 UMIs per Mu-active leaf sample
summary(colSums(MuCounts[,grepl('MI',colnames(MuCounts))]))  # 2,782,187 UMIs per Mu-inactive leaf sample

leafSamps = which(grepl('L',colnames(MuCounts)) & !grepl('_',colnames(MuCounts)))
summary(colSums(MuCounts[,leafSamps] > 0))  # 184,432 insertion sites detected per Mu-active leaf sample


summary(colSums(CPM[historicals, grepl('MI',colnames(CPM))])/10^4)  # 99.8% UMIs map to historical insertions
summary(2 * colSums(CPM[!historicals, grepl('MI',colnames(CPM))]) / colMeans(CPM[historicals, grepl('MI',colnames(CPM))]))  # 0.1 false positive insertions per diploid cell
summary(2 * colSums(CPM[!historicals, grepl('MI',colnames(CPM))]) / colMeans(CPM[historicals, grepl('MI',colnames(CPM))])) / (2*genomesize)  # factor of 2 is because of ploidy. 2.7 x 10^-11 false positive insertions per bp
sd(2 * colSums(CPM[!historicals, grepl('MI',colnames(CPM))]) / colMeans(CPM[historicals, grepl('MI',colnames(CPM))])) / sqrt(3)  # SEM for the false positive rate is 0.04


VAF = VAFnormalize(CPM)  # This function calculates VAF per cell, and should be multiplied by 1/2 for diploid tissues (leaf, root) and pollen and multiplied by 1/3 for triploid endosperm
VAF[[1]] = VAF[[1]][, (apply(VAF[[1]],2,function(xx) { xx=xx[xx>0]; 1/min(xx)  }) >= 1000/2) | grepl('Dilution',colnames(VAF[[1]]))]  # Remove samples that do not reach a VAF of at least 10^-3. This removes 5 early pollen samples that had lower data quality metrics -- adding a second DNA purification of pollen DNA resulted in higher quality libraries with later pollen samples. This filter also removes 1 endosperm sample that reduced library complexity compared to the others
VAF[[1]] = VAF[[1]][,!(colnames(VAF[[1]]) %in% c('E61','E63','E65','E69','E74','E75','E76','E78','E79','E83','E85'))]  # Remove first 11 endosperm samples as they had signs of cross contamination. These samples were used only for determining inherited insertions. The endosperm cross-contamination issue was solved for later samples by using a stringent wash protocol with bleach on all items used for tissue grinding.

summary(apply(VAF$VAF[,!grepl('_',colnames(VAF$VAF))], 2, function(xx) { xx=xx[xx>0]; min(xx)  })*c(1/2,1/3)[1+grepl('E', colnames(VAF$VAF)[!grepl('_',colnames(VAF$VAF))])])  # median detection limit is a VAF of 6.0 x 10^-5; 1 part in 16,569

Vdn = VAFdenovo(VAF)
mutrate = colSums(Vdn[,grepl('P',colnames(Vdn))]/2, na.rm=T)  # Estimate per generation mutation rate from pollen
mean(mutrate)  # mean = 24.0 de novo insertions per generation
sd(mutrate)/sqrt(length(mutrate))  # SEM = 3.9 de novo insertions per generation
mean(mutrate)/genomesize  # 1.1 x 10^-8 de novo insertions per bp per generation
sum(Vdn > 0, na.rm=T)  # 2,577,130 de novo insertions in total



################# Fig 2 #################

aa = read.csv('Transmitted Mu insertions from bulk pollen outcross experiments.csv', row.names=1)  # .csv file available at NCBI GEO as a supplementary file from accession GSE296286
aa = aa[rowSums(aa[,1:5]) == 0, -c(1:5)]


patIns_emp = by(colSums(aa,na.rm=T),sub('.[0-9]$','',colnames(aa)),mean)
summary(patIns_emp)  # median of 27.8 paternal insertions / plant; empirical
boots_emp = sapply(1:2000, function(ii) {
	aa2 = aa[,sample(1:ncol(aa),replace=T)]
	ii = by(colSums(aa2,na.rm=T),sub('.[0-9]$','',colnames(aa2)),mean)
	median(ii)
	})
sd(boots_emp) # std error = 6.7

pollenSamps = which(grepl('^P', colnames(Vdn)))
patIns_VAF = colSums(Vdn[, pollenSamps]/2,na.rm=T)
summary(patIns_VAF)  # median of 22.5 paternal insertions / plant; predicted from VAF
boots_VAF = sapply(1:2000, function(ii) { median(sample(patIns_VAF, replace=T)) })
sd(boots_VAF) # std error = 3.8
wilcox.test(patIns_emp, patIns_VAF)  # p = 0.711

df = data.frame(type = c(rep('emp', length(patIns_emp)), rep('vaf', length(patIns_VAF))), Insertions = c(patIns_emp, patIns_VAF))
Fig2B = ggboxplot(df, 'type', 'Insertions', color = 'type') + scale_y_continuous(limits = c(0,60), expand = c(0,0)) + theme(legend.position="none")


aaM = sapply(1:ncol(aa), function(i) { sapply(1:ncol(aa), function(j) {
	if (any(duplicated(sub('.[0-9]$','',colnames(aa)[c(i,j)])))) {
		return(sum((aa[,i] + aa[,j]) == 2))
	} else { return(NA) }
	}) })
diag(aaM) = NA
colnames(aaM) = colnames(aa)
rownames(aaM) = colnames(aa)

shareIns_emp = by(rowMeans(aaM,na.rm=T),sub('.[0-9]$','',rownames(aaM)),mean)
summary(shareIns_emp)  # median of 0.18 insertions shared between siblings; empirical

shareIns_VAF = colSums((Vdn[, pollenSamps]/2)^2,na.rm=T)
summary(shareIns_VAF)  # median of 0.19 insertions shared between siblings; predicted from VAF
wilcox.test(shareIns_emp, shareIns_VAF)  # p = .938

df = data.frame(type = c(rep('emp', length(shareIns_emp)), rep('vaf', length(shareIns_VAF))), Shared = c(shareIns_emp, shareIns_VAF))
Fig2C = ggboxplot(df, 'type', 'Shared', color = 'type') + scale_y_continuous(limits = c(0,2.5), expand = c(0,0)) + theme(legend.position="none")

Fig2 = ggarrange(Fig2B, Fig2C, labels = c('B','C'))



################# Fig 3 #################

plothist = function(i, combineCategories = 0, lims = c(3*10^-5,1.9), mx = 2, cols = c('#aeaeae', '#1e78b3', '#92dd85')) {
	i = which(colnames(VAF$VAF) == i)
	X = VAF$VAF[,i]/mx
	pats = VAF$paternals[[sub('_.+','',gsub('[a-z]', '', sub('-.+', '', colnames(VAF$VAF)[i]), ignore.case=T))]]
	cat = factor(c('de novo', 'paternal', 'historical')[1 + pats + 2*historicals], levels = c('historical', 'paternal', 'de novo'))
	df = data.frame(VAF = X, cat = cat)
		
	if (combineCategories == 1) { df$cat[df$cat == 'paternal'] = 'de novo' }  # the combineCategories parameter was used to correct the ggplot error when plotting stacked histograms with a log-scale y-axis
	if (combineCategories == 2) { df$cat = 'de novo' }
		
	ggplot(df, aes(x=VAF, fill=cat, color=cat)) +
  geom_histogram(position="stack", aes(y = ..count.. * 2),binwidth=.225) + scale_x_log10(limits = lims, expand = c(0,0)) + scale_y_log10(labels = function(x) x/2, limits = c(.5,200000)*2, expand = c(0,0)) + theme_classic() + annotation_logticks(outside = T) + coord_cartesian(clip = 'off') + ylab('counts') + scale_fill_manual(values=cols) + scale_color_manual(values=cols)
	}

plotPMF = function(i, truncate = 0) {
	X = VAF$VAF[!historicals & !VAF$paternals[[sub('_.+','',gsub('[A-Z]', '', colnames(VAF$VAF)[i], ignore.case=T))]],i]/2
	X = table(X)
	X = X[length(X):(2+truncate)]
	cbind(as.numeric(names(X)),as.numeric(cumsum(X)))
}

interpolateVAF = function(i, truncate = 0) { approx(log(rbind(c(1,0),plotPMF(i, truncate)),10), xout = seq(-5,0,length.out = 200))$y }

tissueVAF = function(tissue, IDs, thresh = .75, truncate = 0) {
		out = 10^sapply(IDs, interpolateVAF, truncate = truncate)
		out[rowMeans(!is.na(out)) < thresh, ] = NA

		boots = sapply(1:2000, function(k) { rowMeans(out[,sample(1:ncol(out),replace=T)], na.rm=T) })
		
		data.frame(tissue = tissue, VAF = 10^seq(-5,0,length.out = 200), counts = rowMeans(boots, na.rm=T), lci = apply(boots,1,quantile,p=.025, na.rm=T), uci = apply(boots,1,quantile,p=.975, na.rm=T))	
	}

plotCDF2 = function(i = NULL, truncate = 0) {
	X2 = tissueVAF('', i, truncate = truncate)[,c(2,3,3)]
	colnames(X2) = c('f','Cf','Nf')
	X2 = X2[!is.na(X2[,2]) & (X2[,2] > 0),]
	fm1 = lm(Cf ~ f, log10(X2))
	coef = c(coefficients(fm1), summary(fm1)[['r.squared']])

	X2 = lapply(i, function(ii) {
		xx = plotPMF(ii, truncate)
		return(data.frame(lib = colnames(VAF[[1]])[ii], f = xx[,1], Cf = xx[,2]))
		})
	X2 = do.call(rbind, X2)
		
	cat('intercept = ', coef[1], '; slope = ', coef[2], '; R2 = ', coef[3], '\n', sep = '')
	
	ggplot(X2, aes(x=f, y = Cf, color = lib)) + geom_line(size = .8) + theme_classic() + scale_y_log10(limits = c(.5,250000)) + scale_x_log10(limits = c(10^-5,1)) + annotation_logticks(outside = T) + coord_cartesian(clip = 'off') + xlab('VAF') + ylab('C(f)') + scale_color_manual(values = rep('#00983d',length(i))) + geom_abline(slope = coef[2], intercept = coef[1], linewidth = 1, col = 'black', linetype = 'dashed') + geom_abline(slope = -1, intercept = log10(.5)*(coef[2]+1)+coef[1], linewidth = 1, col = '#666666', linetype = 'dotted')
	}

leaf = which(grepl('L',colnames(VAF[[1]])) & !grepl('_',colnames(VAF[[1]])))
Fig3A = plothist('L79')  # ggplot does not handle a log scale in a histogram well. Specifically, when there are multiple bars stacked on top of each other (e.g. 5 de novo insertions and 6 historical insertions), it will overestimate the total bar height substantially by treating the individual bars linearly (rather than log-scale). For the stacked portions of the histogram, the plots were edited after the fact to reflect the quantitative number accurately.
Fig3B = plotCDF2(leaf)  # intercept = -1.36; slope = -1.49;R2 = 0.995



################# Fig 4 #################

calcCDF = function(i, mx = 2) {
		X = VAF$VAF[,i]/mx
		pats = VAF$paternals[[sub('_.+','',gsub('[a-z]', '', colnames(VAF$VAF)[i], ignore.case=T))]]
		X = X[!historicals & !pats]
		X = X[X > 0]
		X2 = table(X)
		X2 = data.frame(f = as.numeric(names(X2)), Cf = as.numeric(X2), Nf = as.numeric(X2))
		X2 = X2[order(-X2$f),]
		X2$Cf = cumsum(X2$Cf)
		return(X2)
	}

plotCDF = function(i = NULL, mx = 2, plotit = T) {
	X2 = tissueVAF('', i)[,c(2,3,3)]
	colnames(X2) = c('f','Cf','Nf')
	X2 = X2[!is.na(X2[,2]) & (X2[,2] > 0),]
	fm1 = lm(Cf ~ f, log(X2,10))
	coef = c(coefficients(fm1), summary(fm1)[['r.squared']])

	if (plotit) {		
		alts = lapply(i, calcCDF, mx = mx)
		for (j in 1:length(alts)) { alts[[j]]$lib = j }
		alts = do.call(rbind, alts)
		alts$lib = as.factor(alts$lib)
		cat('intercept = ', coef[1], '; slope = ', coef[2], '; R2 = ', coef[3], '\n', sep = '')
		ggplot(alts, aes(x=f, y = Cf, color = lib)) + geom_line(size = .8) + theme_classic() + scale_y_log10(limits = c(.5,250000)) + scale_x_log10(limits = c(10^-5,1)) + annotation_logticks(outside = T) + coord_cartesian(clip = 'off') + xlab('VAF') + ylab('C(f)') + scale_color_manual(values = rep('gray', length(unique(alts$lib)))) + geom_abline(slope = coef[2], intercept = coef[1], linewidth = 1, col = '#ff0000')
	} else {
		data.frame(intercept = coef[1], slope = coef[2], R2 = coef[3])
	}
	}


set.seed(1)
dat = rbind(tissueVAF('Endosperm', which(grepl('^E', colnames(VAF$VAF)))),
tissueVAF('Root', which(grepl('^R', colnames(VAF$VAF)))),
tissueVAF('Leaf', which(grepl('^L', colnames(VAF$VAF)) & !grepl('_', colnames(VAF$VAF)))),
tissueVAF('Pollen', which(grepl('^P', colnames(VAF$VAF)))))  # 5 E, 6 R, 6 L, 9 P

levyObs = c(rep(0,8), 10, 9, 3, 17, 63, 90, 114, 85, 64, 61)  # Data from Levy and Walbot (1990)
levyKernels = c(rep(2000,8), 100, 50, rep(84900/2^18, 8))
levyX = c(levyObs[levyKernels > 1], cumsum(levyObs[levyKernels <= 1]))
dat = rbind(dat, data.frame(tissue = 'Levy', VAF = 2^-(1:18)/2, counts = 10*cumsum(levyObs/levyKernels), lci = 10*qchisq(0.025, 2*levyX)/2/levyKernels, uci = 10*qchisq(0.975, 2*(levyX+1))/2/levyKernels))


dat = dat[!is.na(dat$lci),]
dat$lci[dat$lci < 0.5] = 0

cols = c('#c076e9','#7cae00','#aaaaaa','#f0a519','#a35800')

Fig4A = ggplot(dat, aes(x = VAF, y = counts, ymin = lci, ymax = uci, group = tissue)) + geom_line(aes(color = tissue)) + geom_ribbon(aes(fill = tissue), alpha = 0.2) + scale_x_log10(limits = c(.5*10^-4, .55), expand = c(0,0)) + scale_y_log10(limits = c(.5,200000), expand = c(0,0)) + theme_classic() + annotation_logticks(outside = T) + coord_cartesian(clip = 'off') + scale_color_manual(values = cols) + scale_fill_manual(values = cols) + theme(legend.position="none")



zzBin = apply(Vdn[,!grepl('x',colnames(Vdn))], 2, function(xx) { 
	bin = ceiling(-log10(xx/2))
	by(xx, bin, sum)[as.character(1:5)] })
zzBin[,grepl('E',colnames(zzBin))] = apply(Vdn[,grepl('E',colnames(Vdn))], 2, function(xx) { 
	bin = ceiling(-log10(xx/3))
	by(xx, bin, sum)[as.character(1:5)] })
	
zzBin[is.na(zzBin)] = 0
zzBin[4,] = colSums(zzBin[4:5,])
zzBin = zzBin[-5,]
zzBin2 = t(apply(zzBin, 1, function(xx) { by(xx, gsub('[0-9]','',colnames(zzBin)), mean) }))
zzBin2 = zzBin2[,c('L','P','R','E')]
zzBin2[,2] = zzBin2[,2]/2  # factor of 2 scaling accounts for meiotic reduction in haploid pollen (if there are 50 insertions in a diploid pollen mother cell, there will be 25 on average per pollen grain)


Fig4B = barplot(zzBin2, border = NA, las = 1, ylab = 'de novo insertions per cell', col = c('#333333', '#666666', '#a8a8a8', '#cccccc'), ylim = c(0,35))

Fig4C = barplot(100*sweep(zzBin2, 2, colSums(zzBin2), '/'), border = NA, las = 1, ylab = 'de novo insertions per cell (%)', col = c('#333333', '#666666', '#a8a8a8', '#cccccc'), ylim = c(0,100))





################# Fig 5 #################

#####  Simulate Robertson (1980) using the empirical pollen allele frequency data
RobertsonOutcross = function(j = 1) {
	VAF1 = na.omit(VAFp[,sample(1:ncol(VAFp),1)])
	mutFreqs = table(unlist(sapply(1:50, function(xx) { which(VAF1 >= runif(length(VAF1))) })))
	mutF = VAF1[as.numeric(names(mutFreqs))]
	out = cbind(mutF, as.numeric(mutFreqs))
	colnames(out) = c('VAF','NumF1')
	return(out)
	}

simRobertson = function(j = 1) { do.call(rbind, sapply(1:30, RobertsonOutcross)) }  # Simulate 30 outcrosses with 50 F1 plante each, for a total of 1500 F1 plants

calcRobertsonResults = function(S) {
	out = do.call(rbind, lapply(S, function(xx) {
	xx = xx[,2]
	c(sum(xx[xx == 1])/sum(xx), sum(xx[xx == 2])/sum(xx), sum(xx[xx >= 3])/sum(xx))
	}))
	colnames(out) = c('0','1','2+')
	return(out)
}
	
calcRobertsonResults2 = function(S) {
	out = do.call(rbind, lapply(S, function(xx) {
	xx = xx[,2]
	c(mean(xx == 1), mean(xx == 2), mean(xx >= 3))
	}))
	colnames(out) = c('0','1','2+')
	return(out)
}
	
	
VAFp = Vdn[, grepl('^P', colnames(Vdn))]/2
VAFp[VAFp == 0] = NA
VAFp = VAFp[rowSums(!is.na(VAFp)) > 0,]  # de novo allelel frequencies for the pollen samples

set.seed(1)
sims = sapply(1:1000, simRobertson)  # Simulate Robertson's study 1000 times; this will take a while and the number of simulations could be reduced for speed by changing the '1000' to '10'


#####  Downsample mutations from simulations to match the numbers from Robertson (1980)
set.seed(1)
simsDownsample = lapply(sims, function(xx) {
		xx[sample(1:nrow(xx), round(nrow(xx)*.116/mean(mutrate))),]
	})  # Robertson found an estimated 154 mutations. This line downsamples the simulated mutations to match Robertson (1980).
mean(unlist(lapply(sims, nrow)))  # 31,954 mutations found prior to downsampling, on average
mean(unlist(lapply(simsDownsample, nrow)))  # 154.2 mutations found after downsampling, on average


simResultsAll = calcRobertsonResults(sims)
simResultsDown = calcRobertsonResults(simsDownsample)
# Both the full simulations and downsampled simulations yield results within error of each other:
colMeans(simResultsAll) # 83.1% F1 have unique mutations, 7.8% shared with 1 sibling, 9.0% shared with 2+ siblings
colMeans(simResultsDown) # 83.4% F1 have unique mutations, 7.8% shared with 1 sibling, 8.8% shared with 2+ siblings
# However, the downsampled results better reflect the variation expected under the conditions of Robertson (1980). Essentially, because Robertson could only track mutations with visible seedling phenotypes, there was greater counting error as only a subset of mutations could be assessed.

SEM = apply(simResultsDown, 2, sd)

r = c(64+107*.7246, 107*(1-.7246)*c(5*2,3*3)/19)  # Robertson's results, calculated from Table 2 and the top paragraph of pg. 975. In the paragraph, you can see he lists that 142 mutations are unique and there were 154 total. These values are equavalent but recalculated with less rounding. Notice that round(r[1]) = 142 and round(sum(r)) = 154
r = r/sum(r)


## Fig. 5B
df = data.frame(source = rep(c('Robertson (1980)', 'Simulation'), each = 3), cat = rep(c('0','1','2+'), 2), val = 100*c(r, colMeans(simResultsDown)), lci = 100*c(NA,NA,NA, colMeans(simResultsDown) - SEM), uci = 100*c(NA,NA,NA, colMeans(simResultsDown) + SEM))
Fig5B = ggbarplot(df, x = 'cat', y = 'val', fill = 'source', position = position_dodge()) + geom_errorbar(aes(x=cat, ymin=lci, ymax=uci)) + scale_y_continuous(limits = c(0,100), expand = c(0,0))

pvals = sweep(simResultsDown, 2, r, '-')
pvals = 2*apply(cbind(colMeans(pvals >= 0), colMeans(pvals <= 0)), 1, min)  # two-tailed p-value
pvals # p = 0.890, 0.632, 0.950




############### Fig. S15 ###############

#####  Simulate Robertson (1980) assuming a constant mutation rate and exponential cell divison (a Luria-Delbrück Process)

mutate = function(xx, mutrate) {  c(xx, sample(1:10^7, rpois(1,mutrate)))  } # 10^7 was arbitrarily chosen as the effective number of genome sites that can be mutated

expon = function(cells, mutrate, divs) {
	for (i in 1:divs) {
		cells = lapply(cells, mutate, mutrate = mutrate) # S-phase
		cells = c(cells,cells) # mitosis
	}
	return(cells)
}

simOutcross = function(i=0, mutrate = 2*171/1541/19) {  # mutation rate is the expected fraction of mutant plants (11%, or 171 / 1541 from Robertson 1980) * 2 (to account for meiotic reduction) / 15 (there are 15 cell divisions in this simulation)
	cells = list(A = NULL) # Initialize cells
	cells = expon(cells, mutrate = mutrate, divs = 16)  # 16 mitotic divisions
	cells = lapply(cells, mutate, mutrate = mutrate) # S-phase for meiosis
	cells = cells[sample(1:length(cells), 50)] # for computational speed, subsample to the 50 offspring now. The chance of sampling the same precursor 2x after 16 divisions is essentially 0
	#cells = c(cells,cells) # meiosis I; unnessecary because chance of sampling multiple products of the same meiosis is essentially 0
	#cells = c(cells,cells) # meiosis II; unnessecary because chance of sampling multiple products of the same meiosis is essentially 0
	cells = lapply(cells, function(xx) { xx[runif(length(xx)) <= .5] })  # Meiotic reductions
	cells = lapply(cells, mutate, mutrate = mutrate) # S-phase for pollen mitosis I
	cells = lapply(cells, mutate, mutrate = mutrate) # S-phase for pollen mitosis II

	F1 = cells[unlist(lapply(cells, length))>0]
	
	F1shared = unlist(lapply(F1, function(xx) { sum(unlist(lapply(F1, function(yy) { any(xx %in% yy) }))) - 1 }))
	F1shared[F1shared > 2] = 2
	out = table(F1shared)[as.character(0:2)]
	out[is.na(out)] = 0
	
	names(out) = c('0','1','2+')
	return(out)
}
	
simRobertsonLD = function(i = 0) {
	out2 = sapply(1:30, simOutcross)
	cat(rowSums(out2)/sum(out2)*100,'\n')
	rowSums(out2)
}


set.seed(1)
zz = t(sapply(1:1000, simRobertsonLD))  # Simulate Robertson's study 1000 times; this will take a VERY long time and the number of simulations could be reduced for speed by changing the '1000' to a lower number

summary(rowSums(zz))
zz2 = sweep(zz, 1, rowSums(zz), '/')  # Values for Fig. S15


## Fig. S15
df = data.frame(source = rep(c('Robertson (1980)', 'Simulation'), each = 3), cat = rep(c('0','1','2+'), 2), val = 100*c(r, colMeans(zz2)), lci = 100*c(NA,NA,NA, colMeans(zz2) - apply(zz2,2,sd)), uci = 100*c(NA,NA,NA, colMeans(zz2) + apply(zz2,2,sd)))
FigS15 = ggbarplot(df, x = 'cat', y = 'val', fill = 'source', position = position_dodge()) + geom_errorbar(aes(x=cat, ymin=lci, ymax=uci)) + scale_y_continuous(limits = c(0,100), expand = c(0,0))


plower = colMeans(sweep(zz2, 2, r, '-') <= 0)
phigher = colMeans(sweep(zz2, 2, r, '-') >= 0)
ptwotail = apply(cbind(plower,phigher), 1, min)*2
ptwotail # p = 0.438, 0.738, 0.358


