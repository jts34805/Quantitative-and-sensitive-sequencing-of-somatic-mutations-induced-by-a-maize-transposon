MuCounts = MuCounts = read.csv('MuSeq2 count data for bulk pollen outcross experiments.csv', row.names=1)
blacklist = read.table('blacklist.txt')[,1]
historicals = read.table('historicals.txt')[,1]
MuCounts = MuCounts[!(rownames(MuCounts) %in% blacklist),(colSums(MuCounts) >= 5000) | grepl('^D',colnames(MuCounts))]
CPM = sweep(MuCounts, 2, colSums(MuCounts), '/')*10^6



#Inherited calls
thresh = 500

CPM2 = 2*(CPM >= thresh) + (CPM >= thresh/2.5)
CPM2 = CPM2[rowSums(CPM2 > 0) > 1,]
CPM2 = CPM2[,(sub('_.','',colnames(CPM2)) %in% names(which(table(sub('_.','',colnames(CPM2))) == 2))) | grepl('D',colnames(CPM2))]

Inh = t(apply(CPM2, 1, function(xx) { by(xx, sub('_.', '', colnames(CPM2)), sum) }))
Inh[,grepl('^D',colnames(Inh))] = Inh[,grepl('^D',colnames(Inh))]*20
Inh = (Inh >= 4) + 0
Inh = Inh[rowSums(Inh[,!grepl('^D',colnames(Inh))]) > 0,]
Inh = cbind(0,Inh)
colnames(Inh)[1] = 'Mu-inactive parent'
Inh[historicals,1] = 1

write.csv(Inh, file = 'Transmitted Mu insertions from bulk pollen outcross experiments.csv')
