removeChrFromGranges = function(granges) {
  seqlevels(granges) = gsub('chr', '', seqlevels(granges))
  granges
}

keepOnlyChromosomesHg = function(granges) {
  granges = subset(granges, seqnames %in% c(1:22,'X','Y','MT','M'))
  granges
}

# Builds GRanges object from a data frame with only 3 columns containing
# in the following order the coordinate information: chrom, start, end.
# The strand is not used (we set it to '*') because we use the generated
# GRanges to calculate overlaps independently from the strand
make.granges = function(df) {
  rd = GRanges(seqnames=df[,1],
               ranges=IRanges(start=df[,2],end=df[,3]),
               strand='*')
  rd
}
