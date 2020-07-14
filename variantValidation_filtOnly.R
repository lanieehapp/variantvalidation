library(parallel)
library(stringr)
library(plyr)

##tmp

args = commandArgs(trailingOnly=TRUE)

vcf_path<-args[[1]]
rna_bam_path<-args[[2]]

ref<-args[[3]]

load(vcf_path)

filt.fix<-single.sample.merged[[4]]
filt.gt<-single.sample.merged[[6]]

var.length<-apply(filt.fix[,4:5], 1, function(x){max(nchar(x))})
end.pos<-as.numeric(filt.fix$POS)+var.length-1
filt.bed<-cbind(filt.fix$CHROM, filt.fix$POS, end.pos)
filt.bed<-unique(filt.bed)
write.table(filt.bed, sep="\t", file="filt_bed.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)



##################### dna filtered variants in rna ####################

file<-rna_bam_path
#comm<-paste0("bcftools mpileup -R whitelist_bed.txt -d 1000 --fasta-ref GRCh38.p12.genome.plus.ERCC.fa ", file, " > rna_vars.txt" )
comm<-paste0('bcftools mpileup -R filt_bed.txt -d 1000 --ignore-RG --fasta-ref ', ref, " ", file, ' | bcftools call -m -A -M | grep -v "##" > /data/rna_filt.txt' )
system(comm)

rna_filt<-read.csv(file="/data/rna_filt.txt", sep="\t", stringsAsFactors = FALSE)
rna_filt$CHROM_POS_REF_ALT<-paste(rna_filt$X.CHROM, rna_filt$POS, rna_filt$REF, rna_filt$ALT, sep="-")

all.info<-as.character(rna_filt$INFO)

print("Parsing INFO column")
tmp<-mclapply(all.info, function(x){
  
  num.values<-str_count(x, ";")
  parsed<-str_split_fixed(x, ";", n=num.values+1)
  
  colnames(parsed)<-unlist(lapply(parsed, function(x){strsplit(x, "=")[[1]][1]}))
  parsed<-gsub(".*=", "", parsed)
  parsed<-data.frame(parsed, stringsAsFactors = FALSE)
  return(parsed)
})
all.info.parsed<-rbind.fill(tmp)
dp4<-data.frame(str_split_fixed(all.info.parsed$DP4, ",", 4), stringsAsFactors = FALSE)
dp4$X1<-as.numeric(dp4$X1)
dp4$X2<-as.numeric(dp4$X2)
dp4$X3<-as.numeric(dp4$X3)
dp4$X4<-as.numeric(dp4$X4)

filt_rna<-NULL
filt_rna$CHROM_POS_REF_ALT<-rna_filt$CHROM_POS_REF_ALT
filt_rna$RNA_REF<-dp4$X1 + dp4$X2
filt_rna$RNA_ALT<-dp4$X3 + dp4$X4
filt_rna$RNA_DEPTH<-filt_rna$RNA_ALT+filt_rna$RNA_REF
filt_rna$RNA_AF<-filt_rna$RNA_ALT/filt_rna$RNA_DEPTH
filt_rna<-data.frame(filt_rna)
filt_rna$RNA_EVIDENCE<-FALSE
filt_rna$RNA_EVIDENCE[filt_rna$RNA_DEPTH<30]<-"LOW_COVERAGE"
filt_rna$RNA_EVIDENCE[dp4$X3>1 & dp4$X4>1]<-TRUE
filt_rna$RNA_EVIDENCE[str_count(rna_filt$ALT, ",")>0]<-"MULTIALLELIC"

filt_rna$CHROM_POS<-paste0(rna_filt$X.CHROM, "-", rna_filt$POS)
filt_rna$ALT<-rna_filt$ALT

fix_pos<-paste0(filt.fix$CHROM, "-", filt.fix$POS)

filt.gt$RNA_EVIDENCE<-NA
filt.gt$RNA_AF<-NA
filt.gt$RNA_DEPTH<-NA
for(i in 1:nrow(filt_rna)){
  if(filt_rna$CHROM_POS[i] %in% fix_pos){
    idx<-which(fix_pos == filt_rna$CHROM_POS[i])
    filt.gt$RNA_EVIDENCE[idx]<-filt_rna$RNA_EVIDENCE[i]
    filt.gt$RNA_AF[idx]<-filt_rna$RNA_AF[i]
    filt.gt$RNA_DEPTH[idx]<-filt_rna$RNA_DEPTH[i]
  }
}

samp.id<-unlist(str_split(colnames(filt.gt)[2], "[.]"))[1]

colnames(filt.gt)[29:31]<-c(paste0(samp.id, ".RNA_EVIDENCE"), paste0(samp.id, ".RNA_AF"), paste0(samp.id, ".RNA_DEPTH"))

single.sample.filt<-list(single.sample.merged[[4]], single.sample.merged[[5]], filt.gt)

save(single.sample.filt, file="/data/single_sample_filtered_withRNA.RData")

simple.single.sample<-cbind(filt.fix, filt.gt[,grep("nCallers|afMax|dpMax|RNA_EVIDENCE|RNA_AF|RNA_DEPTH", colnames(filt.gt))])
write.table(simple.single.sample, file="/data/simple_single_sample.tsv", quote=FALSE, row.names = FALSE, sep="\t")

