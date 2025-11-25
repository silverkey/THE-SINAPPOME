# Code to build a TxDb using ensembl GTF annotation for the human genome
# The GTF gives less warnings than GFF3 at least in hg19
# We need to use the GenomicFeatures package
# https://www.rdocumentation.org/packages/GenomicFeatures/versions/1.24.4/topics/TxDb-class
# https://bioconductor.org/packages/devel/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.pdf

# Start R and build the db from Ensembl GTF...
library(rtracklayer)
library(GenomicFeatures)

gtf    = 'Mus_musculus.GRCm39.108.gtf'
sqlite = 'Mus_musculus.GRCm39.108.sqlite'
rds    = 'Mus_musculus.GRCm39.108.RDS'
source = 'ensembl'
org    = 'Mus musculus'

# The database is directly loaded by using the downloaded GTF. It requires a bit of time.
# http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz
txdb = makeTxDbFromGFF(file = gtf, dataSource = source, organism = org)

# Read also the GTF in order to have access to metadata (last column of the GTF)
# rtracklayer will make a GRanges object of all the info in the GTF
gtf = rtracklayer::import(gtf)

# We save tha database in an sqlite file. It can be reloaded using loadDb()
saveDb(txdb, file = sqlite)

# And let's save also the GTF object
saveRDS(gtf,file = rds)

# We make a data.frame of the gtf
gtfdf = as.data.frame(gtf)

# And we can use it to do all the subsetting we need to operate on the txdb
# or we can attach some info as mcols to the GRanges obtained from the txdb
table(subset(gtfdf,type=='gene')$gene_biotype)
table(gtfdf$seqnames)
