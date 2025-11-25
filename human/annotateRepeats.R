# Genome packages to eventually install and use
# BSgenome.Mmusculus.UCSC.mm10
# BSgenome.Mmusculus.UCSC.mm39
# BSgenome.Hsapiens.UCSC.hg19
# BSgenome.Hsapiens.UCSC.hg38

library(GenomicFeatures)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
source('utils.R')

genome = BSgenome.Hsapiens.UCSC.hg38

# The files with repeats and transcriptome db saved previously
rmrf = 'rmr.RDS'
txdbf = 'Homo_sapiens.GRCh38.108.sqlite'
gtfrds = 'Homo_sapiens.GRCh38.108.RDS'
# TE = c('DNA','LTR','Retroposon','SINE','LINE','RC')
TE = c('SINE')

# Load the databases
rmr = readRDS(rmrf)
txdb = loadDb(txdbf)
gtf = readRDS(gtfrds)

# Look at the repeats classes
table(as.data.frame(rmr)$repClass)
# Filter repeats keeping only specific classes defined in TE variable
rmrs = subset(rmr, repClass %in% TE)
# Remove the 'chr' string from the chromosome names
rmrs = removeChrFromGranges(rmrs)
# Keep only chromosomes, no aplotype, no unplaced
rmrs = keepOnlyChromosomesHg(rmrs)
# Look at the filtered repeats chromosomes
table(as.data.frame(rmrs)$seqnames)[table(as.data.frame(rmrs)$seqnames) > 0]

# Generate the ranges to work on
genes = genes(txdb)
transcripts = transcripts(txdb)
exons = exons(txdb)
introns = intronsByTranscript(txdb)
promoters = promoters(txdb)
utr5 = fiveUTRsByTranscript(txdb)
utr3 = threeUTRsByTranscript(txdb)
transgenes = transcriptsBy(txdb, 'gene')
exgenes = exonsBy(txdb, 'gene')
extranscripts = exonsBy(txdb, 'tx')

# Make dataframes from ranges
gdf = as.data.frame(genes)
tdf = as.data.frame(transcripts)
edf = as.data.frame(exons)
idf = as.data.frame(introns)
pdf = as.data.frame(promoters)
u5df = as.data.frame(utr5)
u3df = as.data.frame(utr3)
tbgdf = as.data.frame(transgenes)
ebgdf = as.data.frame(exgenes)
ebtdf = as.data.frame(extranscripts)

egtfdf = subset(as.data.frame(gtf), type=='exon')
table(egtfdf$source)
ggtfdf = subset(as.data.frame(gtf), type=='gene')
table(ggtfdf$source)

# Calculate overlaps
gov = findOverlaps(rmrs, genes, ignore.strand=T)
tov = findOverlaps(rmrs, transcripts, ignore.strand=T)
eov = findOverlaps(rmrs, exons, ignore.strand=T)
pov = findOverlaps(rmrs, promoters, ignore.strand=T)
u5ov = findOverlaps(rmrs, utr5, ignore.strand=T)
u3ov = findOverlaps(rmrs, utr3, ignore.strand=T)

# Create dataframe containing info of the exons and the overlapping feature
eovdf = cbind(as.data.frame(rmrs[queryHits(eov),]),
              as.data.frame(exons[subjectHits(eov),]))

cn = c('rchr','rstart','rend','rlength','rstrand','rswscore',
       'rmillidiv','rmillidel','rmillins','rname','rclass',
       'rfamily','rconstart','rconsend','rconsleft','echr',
       'estart','eend','elength','estrand','exon_id')
colnames(eovdf) = cn

# Calculate the overlap coordinates and its length
eovdf$ovstart = ifelse(eovdf$rstart > eovdf$estart, eovdf$rstart, eovdf$estart)
eovdf$ovend = ifelse(eovdf$rend < eovdf$eend, eovdf$rend, eovdf$eend)
eovdf$ovlength = eovdf$ovend - eovdf$ovstart + 1
eovdf$ovstrand = ifelse(eovdf$rstrand == eovdf$estrand, 'dir', 'inv')

# Prepare the merged table with all the info
eovdf = merge(eovdf, ebtdf[,c('exon_id','exon_name')], by.x='exon_id', by.y='exon_id')
eovdf = merge(eovdf, egtfdf, by.x='exon_name', by.y='exon_id')
#eovdf = eovdf[,c(1:26,32,33,36,38,40,41,43,45,48,49,51)]

# Calculate number of nucleotides upstream the subject repeat for potential BD
eovdf$bpforbd = ifelse(eovdf$estrand=='+',eovdf$ovstart-eovdf$estart,eovdf$eend-eovdf$ovend)

# Filter for repeat in the first exon where there is no space left for potential BD
# with filter on at least 30 bp left upstream
eovdf$spaceforbd = ifelse((eovdf$exon_number == 1 & eovdf$bpforbd >= 30) | eovdf$exon_number > 1, 1, 0)

# Start to aggregate. Here select the unique repeat bystart/end...
# Then we will need to aggregate also for start/end/length of the overlap because the same
# repeat might have different overlaps on different exons and/or be of different lengths...
runiq = unique(eovdf[,c('rchr','rstart','rend','rlength','rstrand','rswscore','rmillidiv',
                        'rmillidel','rmillins','rname','rclass','rfamily','rconstart',
                        'rconsend','rconsleft','gene_id','gene_name','gene_biotype','estrand',
                        'transcript_id','transcript_support_level','tag',
                        'ovstart','ovend','ovlength','ovstrand','exon_number','bpforbd','spaceforbd')])

# Filter for overlap length >= 50 and inverted orientation and biotype not protein coding
# and space left for BD
filtered = subset(runiq, ovlength>=50 & ovstrand=='inv' & gene_biotype!='protein_coding' & spaceforbd==1)

# Create the variable that will be used to aggregate
filtered$rid = apply(filtered,1,function(x) paste(x[1],as.numeric(x[2]),as.numeric(x[3]),sep='_'))

# Select only the overlapping region and the longest overlap for a given repeat
agg = aggregate(filtered$ovlength,list(filtered$rid),max)
longer = merge(filtered, agg, by.x=c('rid','ovlength'), by.y=c('Group.1','x'))

# Extract the sequences. We give one nucleotide more just in case...
seq = getSeq(genome, names=paste('chr',longer$rchr,sep=''), start=longer$ovstart-1, end=longer$ovend+1, strand=longer$rstrand)
longer$seq = as.character(seq)

# Take only unique sequences. Create a unique id for every unique sequence and add it to the df.
# This takes into account only perfect identity over the entire sequence. Does not check for one sequence contained in another...
useq = data.frame(seq=unique(longer$seq))
useq$seqid = 1:nrow(useq)
longer = merge(longer,useq,by='seq')

# Considerare anche filtro su milli...
longer$rchanges = (longer$rmillidiv+longer$rmillidel+longer$rmillins)/1000*longer$rlength
longer$rmillisum = (longer$rmillidiv+longer$rmillidel+longer$rmillins)
longer$transcript_support = gsub(' (.+)','',longer$transcript_support_level)
columns = c('seqid','rid','gene_id','gene_name','gene_biotype','transcript_id','transcript_support','tag',
            'bpforbd','rname','rclass','rfamily','rchanges','rmillisum','ovlength','ovstrand','seq')

# Write output tables
write.table(unique(eovdf),file='exon_SINE_overlap.xls',col.names=T,sep='\t',row.names=F,quote=F)
write.table(unique(runiq),file='overlapping_SINEs.xls',col.names=T,sep='\t',row.names=F,quote=F)
write.table(unique(longer[,columns]),file='longer_noncoding_overlapping_SINEs.xls',col.names=T,sep='\t',row.names=F,quote=F)

# A bit of plotting.... just testing....
sel = subset(eovdf, ovlength>=80 & ovlength<=150)
table(sel$rfamily,sel$ovstrand,sel$gene_biotype)

sel1 = subset(eovdf, ovlength<=20)
table(sel1$rfamily,sel1$ovstrand,sel1$gene_biotype)

sel2 = subset(eovdf, ovlength>=50)
table(sel2$rfamily,sel2$ovstrand,sel2$gene_biotype)

# Plot overlaps length divided by the strand
par(mfrow=c(2,1))
hist(subset(eovdf,ovstrand=='dir')$ovlength,breaks=1:550,xlim=c(1,350))
abline(v=80,col='red')
abline(v=130,col='red')
hist(subset(eovdf,ovstrand=='inv')$ovlength,breaks=1:550,xlim=c(1,350))
abline(v=80,col='red')
abline(v=130,col='red')

par(mfrow=c(2,1))
hist(subset(runiq,ovstrand=='dir')$ovlength,breaks=1:550,xlim=c(1,350))
abline(v=80,col='red')
abline(v=130,col='red')
abline(h=130,col='red')
hist(subset(runiq,ovstrand=='inv')$ovlength,breaks=1:550,xlim=c(1,350))
abline(v=80,col='red')
abline(v=130,col='red')
abline(h=130,col='red')

