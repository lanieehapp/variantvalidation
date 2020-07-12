get_vcf = function(file=file) {
    comm = paste("/Users/laniehapp/google-cloud-sdk/bin/gsutil cp gs://analysis_results/",file, "/variantfiltering/*.RData ~/tmp" ,sep="")
    
    print(comm)
    system(comm)
    
    comm = paste('/usr/local/bin/aws s3 cp s3://davelab-analysis-results/', file, '/variantfiltering/ ~/tmp --recursive --exclude "*" --include "*.RData"' , sep="")
    print(comm)
    system(comm)
}

get_dna_bam = function(file=file) {
  comm = paste("/Users/laniehapp/google-cloud-sdk/bin/gsutil cp gs://analysis_results/",file, "/gatk_filtersamreads__nonspliced__dna/*.bam ~/tmp" ,sep="")
  
  print(comm)
  system(comm)
  
  comm = paste("/Users/laniehapp/google-cloud-sdk/bin/gsutil cp gs://analysis_results/",file, "/samtools_index__gatk_filtersamreads__nonspliced__dna/*.bam.bai ~/tmp" ,sep="")
  
  print(comm)
  system(comm)
  
  comm = paste('/usr/local/bin/aws s3 cp s3://davelab-analysis-results/', file, '/gatk_filtersamreads__nonspliced__dna/ ~/tmp --recursive --exclude "*" --include "*.bam"' , sep="")
  print(comm)
  system(comm)
  
  comm = paste('/usr/local/bin/aws s3 cp s3://davelab-analysis-results/', file, '/samtools_index__gatk_filtersamreads__nonspliced__dna/ ~/tmp --recursive --exclude "*" --include "*.bam.bai"' , sep="")
  print(comm)
  system(comm)
}

get_rna_bam = function(file=file) {
  comm = paste("/Users/laniehapp/google-cloud-sdk/bin/gsutil cp gs://analysis_results/",file, "/picard_markduplicates__rna/*.bam ~/tmp" ,sep="")
  
  print(comm)
  system(comm)
  
  comm = paste("/Users/laniehapp/google-cloud-sdk/bin/gsutil cp gs://analysis_results/",file, "/samtools_index__picard_markduplicates__rna/*.bam.bai ~/tmp" ,sep="")
  
  print(comm)
  system(comm)
  
  comm = paste('/usr/local/bin/aws s3 cp s3://davelab-analysis-results/', file, '/picard_markduplicates__rna/ ~/tmp --recursive --exclude "*" --include "*.bam"' , sep="")
  print(comm)
  system(comm)
  
  comm = paste('/usr/local/bin/aws s3 cp s3://davelab-analysis-results/', file, '/samtools_index__picard_markduplicates__rna/ ~/tmp --recursive --exclude "*" --include "*.bam.bai"' , sep="")
  print(comm)
  system(comm)
}
