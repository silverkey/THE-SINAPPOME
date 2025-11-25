# Code to create a GRanges annotation of the repeatmasker annotation from UCSC
# We made a little conversion of coordinates on the consensus to allow a more
# standardized annotation in which start is always smaller than end.

library(GenomicRanges)

# The file downloaded from UCSC database
rmskfile = 'rmsk.txt'
# rmt ---> RepeatMasker Table
rmt = read.table(rmskfile,head=F,comment.char='',sep='\t')

# This information are taken from the 'rmsk.sql' file from UCSC downloaded
# togheter with the 'rmsk.txt'
colnames(rmt) = c('bin','swScore','milliDiv','milliDel','milliIns',
                  'genoName','genoStart','genoEnd','genoLeft','strand',
                  'repName','repClass','repFamily','repStart','repEnd',
                  'repLeft','id')

# Workaround to represent the coordinates on the repeat consensus in a more
# standardized way with start always smalle than end
rmt$consStart = ifelse(rmt$strand=='+',rmt$repStart,rmt$repLeft)
rmt$consEnd = rmt$repEnd
rmt$consLeft = - ifelse(rmt$strand=='+',rmt$repLeft,rmt$repStart)

# rmr ---> RepeatMasker Ranges
rmr = GRanges(seqnames=rmt$genoName,
              ranges=IRanges(start=rmt$genoStart,end=rmt$genoEnd),
              strand=rmt$strand)

# Other important columns regarding the annotation of repeats are stored
# as metadata of the rmr object
mcols(rmr) = rmt[,c('swScore','milliDiv','milliDel','milliIns',
                    'repName','repClass','repFamily',
                    'consStart','consEnd','consLeft')]

# We save the Granges object (rmr). It can be reloaded using readRDS()
saveRDS(rmr,file='rmr.RDS')

# To generate a data.frame from the GRanges object
rmdf = as.data.frame(rmr)

