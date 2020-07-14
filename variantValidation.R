library(parallel)
library(stringr)
library(plyr)

##tmp

args = commandArgs(trailingOnly=TRUE)

vcf_path<-args[[1]]
dna_bam_path<-args[[2]]
rna_bam_path<-args[[3]]

ref<-args[[4]]

#source("multisample_merger_functions.R")
#get_vcf(file=vf_id)
#get_dna_bam(file=align_id)
#get_rna_bam(file=align_id)

#get variant location regions in bed format - from whitelist (RNA + DNA) and filtered variants (RNA only)

#whitelist - only need to run when whitelist is updated
#load("~/Dropbox (DaveLab)/Lanie_resources/variantFilter/whitelist.RData")
#var.length<-apply(whitelist[,4:5], 1, function(x) {max(nchar(x))})
#end.pos<-(whitelist$Start)+var.length-1
#whitelist.bed<- cbind(as.character(whitelist$Chr), whitelist$Start, end.pos )
#whitelist.bed<-unique(whitelist.bed)
#write.table(whitelist.bed, sep="\t", file="whitelist_bed.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
load(vcf_path)

filt.fix<-single.sample.merged[[4]]
filt.gt<-single.sample.merged[[6]]

var.length<-apply(filt.fix[,4:5], 1, function(x){max(nchar(x))})
end.pos<-as.numeric(filt.fix$POS)+var.length-1
filt.bed<-cbind(filt.fix$CHROM, filt.fix$POS, end.pos)
filt.bed<-unique(filt.bed)
write.table(filt.bed, sep="\t", file="filt_bed.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

################## dna whitelist ##################
file<-dna_bam_path
#file<-paste0("/data/", file)

#comm<-paste0("bcftools mpileup -R whitelist_bed.txt -d 1000 --fasta-ref GRCh38.p12.genome.plus.ERCC.fa ", file, " > rna_vars.txt" )
comm<-paste0('bcftools mpileup -R whitelist_bed.txt -d 1000 --ignore-RG --fasta-ref ',ref, ' ', file, ' | bcftools call -m -A -M | grep -v "##" > dna_whitelist.txt' )
system(comm)

dna_whitelist<-read.csv(file="dna_whitelist.txt", sep="\t", stringsAsFactors = FALSE)
dna_whitelist$CHROM_POS_REF_ALT<-paste(dna_whitelist$X.CHROM, dna_whitelist$POS, dna_whitelist$REF, dna_whitelist$ALT, sep="-")

all.info<-as.character(dna_whitelist$INFO)

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

whitelist<-NULL
whitelist$CHROM_POS_REF_ALT<-dna_whitelist$CHROM_POS_REF_ALT
whitelist$DNA_REF<-dp4$X1 + dp4$X2
whitelist$DNA_ALT<-dp4$X3 + dp4$X4
whitelist$DNA_DEPTH<-whitelist$DNA_ALT+whitelist$DNA_REF
whitelist$DNA_AF<-whitelist$DNA_ALT/whitelist$DNA_DEPTH

whitelist<-data.frame(whitelist)

whitelist$DNA_EVIDENCE<-FALSE
whitelist$DNA_EVIDENCE[whitelist$DNA_DEPTH<30]<-"LOW_COVERAGE"
whitelist$DNA_EVIDENCE[dp4$X3>0 & dp4$X4>0]<-TRUE
whitelist$DNA_EVIDENCE[str_count(dna_whitelist$ALT, ",")>0]<-"MULTIALLELIC"


##################### rna whitelist ####################
file<-rna_bam_path
#file<-paste0("/", file)

#comm<-paste0("bcftools mpileup -R whitelist_bed.txt -d 1000 --fasta-ref GRCh38.p12.genome.plus.ERCC.fa ", file, " > rna_vars.txt" )
comm<-paste0('bcftools mpileup -R whitelist_bed.txt -d 1000 --ignore-RG --fasta-ref ', ref, ' ', file, ' | bcftools call -m -A -M | grep -v "##" > rna_whitelist.txt' )
system(comm)

rna_whitelist<-read.csv(file="rna_whitelist.txt", sep="\t", stringsAsFactors = FALSE)
rna_whitelist$CHROM_POS_REF_ALT<-paste(rna_whitelist$X.CHROM, rna_whitelist$POS, rna_whitelist$REF, rna_whitelist$ALT, sep="-")

all.info<-as.character(rna_whitelist$INFO)

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

whitelist_rna<-NULL
whitelist_rna$CHROM_POS_REF_ALT<-rna_whitelist$CHROM_POS_REF_ALT
whitelist_rna$RNA_REF<-dp4$X1 + dp4$X2
whitelist_rna$RNA_ALT<-dp4$X3 + dp4$X4
whitelist_rna$RNA_DEPTH<-whitelist_rna$RNA_ALT+whitelist_rna$RNA_REF
whitelist_rna$RNA_AF<-whitelist_rna$RNA_ALT/whitelist_rna$RNA_DEPTH
whitelist_rna<-data.frame(whitelist_rna)
whitelist_rna$RNA_EVIDENCE<-FALSE
whitelist_rna$RNA_EVIDENCE[whitelist_rna$RNA_DEPTH<30]<-"LOW_COVERAGE"
whitelist_rna$RNA_EVIDENCE[dp4$X3>0 & dp4$X4>0]<-TRUE
whitelist_rna$RNA_EVIDENCE[str_count(rna_whitelist$ALT, ",")>0]<-"MULTIALLELIC"



all_whitelist<-merge(x=whitelist, y=whitelist_rna, by="CHROM_POS_REF_ALT", all.x=T, all.y=T)


##################### dna filtered variants in rna ####################

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
for(i in 1:nrow(filt_rna)){
  if(filt_rna$CHROM_POS[i] %in% fix_pos){
    idx<-which(fix_pos == filt_rna$CHROM_POS[i])
    filt.gt$RNA_EVIDENCE[idx]<-filt_rna$RNA_EVIDENCE[i]
  }
}

samp.id<-unlist(str_split(colnames(filt.gt)[2], "[.]"))[1]

colnames(filt.gt)[36]<-paste0(samp.id, ".RNA_EVIDENCE")

#all_whitelist<-merge(x=whitelist, y=whitelist_rna, by="CHROM_POS_REF_ALT", all.x=T, all.y=T)

single.sample.merged[[6]]<-filt.gt

single.sample.filt<-list(single.sample.merged[[4]], single.sample.merged[[5]], filt.gt)

save(single.sample.filt, file="/data/single_sample_filtered_withRNA.RData")
#write.table(filt_rna, file="/data/filtered_vars_in_RNA.tsv", quote=FALSE, sep="\t", row.names = FALSE)
write.table(all_whitelist, file="/data/all_whitelist_vars.tsv", quote=FALSE, sep="\t", row.names=FALSE)


simple.single.sample<-cbind(filt.fix, filt.gt[,grep("nCallers|afMax|dpMax|RNA_EVIDENCE", colnames(filt.gt))])
write.table(simple.single.sample, file="/data/simple_single_sample.tsv", quotes=FALSE, row.names = FALSE, sep="\t")

