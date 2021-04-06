usage = "
--------------------------------------------------
Annotate GRIDSS VCF, then filtered PASS and Simple 
(INS, INV, DEL, DUP) records. 
--------------------------------------------------
Usage: 
    <Program> <input_gridss.vcf> <out_prefix>

Will write:
  anntated vcf to \"out_prefix_annotated.vcf\"
  PASS-simple-records to \"out_prefix_simple.bed\"
  PASS-itx-records to \"out_prefix_itx.bed\"
  unPASS records to \"out_prefix_unpass.bed\"
--------------------------------------------------
gst 2018-12-03
"
argv <- commandArgs(TRUE)
if(length(argv) != 2){stop(paste0("\nWrong num of input files\n",usage))}
input=argv[1]
pre=argv[2]
#source("https://bioconductor.org/biocLite.R")
#biocLite("VariantAnnotation")
#install_github("PapenfussLab/StructuralVariantAnnotation")
#install.packages("stringr")
library(VariantAnnotation)
library(StructuralVariantAnnotation)
library(stringr)

#' Simple SV type classifier
simpleEventType <- function(gr) {
  return(ifelse(seqnames(gr) != seqnames(partner(gr)), "ITX", # inter-chromosomosal
          ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS", # TODO: improve classification of complex events
           ifelse(strand(gr) == strand(partner(gr)), "INV",
            ifelse(xor(start(gr) < start(partner(gr)), strand(gr) == "-"), "DEL",
             "DUP")))))
}
# using the example in the GRIDSS /example directory
vcf <- readVcf(input)
gr <- breakpointRanges(vcf)
svtype <- simpleEventType(gr)
info(vcf)$SIMPLE_TYPE <- NA_character_
info(vcf[gr$vcfId])$SIMPLE_TYPE <- svtype
info(vcf[gr$vcfId])$SVLEN <- gr$svLen
writeVcf(vcf, paste0(pre,"_annotated.vcf"))

# TODO: perform event filtering here
# By default, GRIDSS is very sensitive but this comes at the cost of a high false discovery rate
# unpassgr <- gr[ ! gr$FILTER == "PASS" & partner(gr)$FILTER == "PASS"] # low confidence calls
gr <- gr[gr$FILTER == "PASS" & partner(gr)$FILTER == "PASS"] # Remove low confidence calls

simplegr <- gr[simpleEventType(gr) %in% c("INS", "INV", "DEL", "DUP")]
simplebed <- data.frame(
    chrom=seqnames(simplegr),
	# call the centre of the homology/inexact interval
    start=as.integer((start(simplegr) + end(simplegr)) / 2),
    end=as.integer((start(partner(simplegr)) + end(partner(simplegr))) / 2),
    name=simpleEventType(simplegr),
    score=simplegr$QUAL,
    strand=".",
    svID=paste0(pre,"_",simplegr$vcfId),
    svLen=simplegr$svLen,
    startf=abs(end(simplegr) - start(simplegr)),
    endf=abs(end(partner(simplegr)) - start(partner(simplegr)))
    )
# Just the lower of the two breakends so we don't output everything twice
simplebed <- simplebed[simplebed$start < simplebed$end,]
write.table(simplebed, paste0(pre,"_simple.bed"), quote=FALSE, sep='\t', row.names=FALSE, col.names=T)

itxgr <- gr[ ! simpleEventType(gr) %in% c("INS", "INV", "DEL", "DUP")]
itxbed <- data.frame(
    chrom=seqnames(itxgr),
  # call the centre of the homology/inexact interval
    start=as.integer((start(itxgr) + end(itxgr)) / 2),
    end=as.integer((start(partner(itxgr)) + end(partner(itxgr))) / 2),
    name=simpleEventType(itxgr),
    score=itxgr$QUAL,
    strand=".",
    svID=paste0(pre,"_",itxgr$vcfId),
    svLen=itxgr$svLen,
    startf=abs(end(itxgr) - start(itxgr)),
    endf=abs(end(partner(itxgr)) - start(partner(itxgr)))
    )
# Just the lower of the two breakends so we don't output everything twice
itxbed <- itxbed[itxbed$start < itxbed$end,]
write.table(itxbed, paste0(pre,"_itx.bed"), quote=FALSE, sep='\t', row.names=FALSE, col.names=T) 

# unpassbed <- data.frame(
#     chrom=seqnames(unpassgr),
#   # call the centre of the homology/inexact interval
#     start=as.integer((start(unpassgr) + end(unpassgr)) / 2),
#     end=as.integer((start(partner(unpassgr)) + end(partner(unpassgr))) / 2),
#     name=simpleEventType(unpassgr),
#     score=unpassgr$QUAL,
#     strand=".",
#     svID=unpassgr$vcfId,
#     svLen=unpassgr$svLen,
#     startf=abs(end(unpassgr) - start(unpassgr)),
#     endf=abs(end(partner(unpassgr)) - start(partner(unpassgr)))
#     )
# # Just the lower of the two breakends so we don't output everything twice
# unpassbed <- unpassbed[unpassbed$start < unpassbed$end,]
# write.table(unpassbed, paste0(pre,"_unpass.bed"), quote=FALSE, sep='\t', row.names=FALSE, col.names=T) 
