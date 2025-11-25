#------------------
# PARAMETERS
#------------------

# File containing the potential SINEUP lncRNAs in the whole transcriptome
sineup.file = '/home/remo/Dropbox/WORK/PROJECTS/overlappeR/SINEUPPOME/MUS/longer_noncoding_overlapping_SINEs.xls'

# Pattern present in the name of the expression table (must not be present in other files in the folder)
pattern = 'n2a'

# Cutoff on the expression level in RPM that should be satisfied by minexp.cutoff
exp.cutoff = 5
minexp.cutoff = 3

# Cutoff on transcript support level
tsl.cutoff = 1

#------------------
# FUNCTIONS
#------------------

# Function to read the counts from STAR and taking the expression values of the reverse
# Make also the calculation of the RPM and returns 3 tables
read.star.counts.reverse = function(filename) {
  t = read.table(filename)
  n = t[1:4,c(1,4)]
  t = t[5:nrow(t),c(1,4)]
  filename = gsub('_ReadsPerGene.out.tab','',filename)
  colnames(t) = c('gene_id',filename)
  colnames(n) = c('type',filename)
  tot.t = sum(t[,2])
  tot.n = sum(n[,2])
  tot = tot.t + tot.n
  nt = t
  nt[,2] = nt[,2]/tot*1000000
  n[5,] = c('counted',tot.t)
  n[6,] = c('uncounted',tot.n)
  n[7,] = c('tot',tot)
  return(list(t,n,nt))
}

#------------------
# LET'S THE CODE BEGIN
#------------------

# Read the count tables containing the pattern in their names
files = dir(pattern=pattern)

# Read the first sample and create the data frames
t = read.star.counts.reverse(files[1])
data = t[[1]]
summ = t[[2]]
norm = t[[3]]

# Add the other samples to the data frames
for(i in 2:length(files)) {
  t = read.star.counts.reverse(files[i])
  data = merge(data,t[[1]],by='gene_id',sort=F)
  summ = merge(summ,t[[2]],by='type',sort=F)
  norm = merge(norm,t[[3]],by='gene_id',sort=F)
}

# Write the tables containg info on all the experiments with pattern in their name
write.table(data,file=paste(pattern,'counts.xls',sep='_'),col.names=T,sep='\t',row.names=F,quote=F)
write.table(summ,file=paste(pattern,'summary.xls',sep='_'),col.names=T,sep='\t',row.names=F,quote=F)
write.table(data,file=paste(pattern,'normalized.xls',sep='_'),col.names=T,sep='\t',row.names=F,quote=F)

# Keep only the genes showing an expression level >= exp.cutoff in >= minexp.cutoff experiments
norm.sel = norm[rowSums(norm[,-1] >= exp.cutoff) >= minexp.cutoff,]

# Read the table with the potential SINEUPS
sineup = read.table(sineup.file,head=T,sep='\t')

# Merge the SINEUPS table with the filtered expression table
res = merge(norm.sel,sineup,by='gene_id',sort=F)
write.table(res,file=paste(pattern,exp.cutoff,'rpm_notslcutoff.xls',sep='_'),col.names=T,sep='\t',row.names=F,quote=F)

# Take out the NA changing them with 0
res[is.na(res)] = 0

# The analysis is a the gene level and we are interested in repeats overlapping them.
# Here we take out the info about transcript to eliminate redundancy.
# Consider to develop the pipeline to work at the transcript level
# If something is changed in the SINEUPS table you need to change numbering here!!!
uni.res = unique(res[,c(1:ncol(data),(ncol(data)+1):(ncol(data)+4),(ncol(data)+9):ncol(res))])

# Info and write table with the expression filter
dim(uni.res)
length(unique(uni.res$gene_id))
write.table(uni.res,file=paste(pattern,exp.cutoff,'rpm_notslcutoff_unique.xls',sep='_'),col.names=T,sep='\t',row.names=F,quote=F)

# Filter the table also on the TSL cutoff and write it
sel = subset(res,transcript_support<=tsl.cutoff)
write.table(sel,file=paste(pattern,exp.cutoff,'rpm',tsl.cutoff,'maxtsl.xls',sep='_'),col.names=T,sep='\t',row.names=F,quote=F)

# Take out transcript info to avoid redundancy
# If something is changed in the SINEUPS table you need to change numbering here!!!
uni.sel = unique(sel[,c(1:ncol(data),(ncol(data)+1):(ncol(data)+4),(ncol(data)+9):ncol(sel))])

# Info and write table with the expression and the TSL filter
dim(uni.sel)
length(unique(uni.sel$gene_id))
write.table(uni.sel,file=paste(pattern,exp.cutoff,'rpm',tsl.cutoff,'maxtsl_unique.xls',sep='_'),col.names=T,sep='\t',row.names=F,quote=F)
