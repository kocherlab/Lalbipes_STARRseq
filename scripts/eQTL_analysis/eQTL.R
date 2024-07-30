library(MatrixEQTL); library(Biobase); library(devtools); library(broom); library(dplyr); library(reshape2)
library(tidyr); library(snpsettest); library(stringr); library(GenomicRanges); library(ggplot2); library(caret)
library(ggstatsplot); library(magrittr); library(gridExtra)

expr<-read.table("FCs_all_peaks.txt",header=T)
conv<-read.table("Filtered_Enhancers_36216.gff",header=F)
names(conv)<-c("chr","source","type","start","end","strand","score","something","name")
conv$enhancer<-gsub("ID=","",conv$name)
conv$region<-paste(conv$chr,":",conv$start,"-",conv$end,sep="")
dim(conv)
# 36216     11

# Need to reorder columns to match allele data
expr<-expr[,c(1,2,3,10,11,12,4,5,6,13,14,15,16,17,18,7,8,9)]

alleles<-read.table(gzfile("LALB_v3_STAR.ad.txt.gz"),header=T)
# dim(alleles) 2,756,082      19

snp_info<-data.frame(alleles$allele)
snp_info[c('chr', 'pos', 'alleles')] <- str_split_fixed(snp_info$alleles.allele, ':', 3)
names(snp_info)<-c("id","chr","pos","refalt")

enhancer_info<-data.frame(conv$enhancer,conv$chr,conv$start,conv$end)
names(enhancer_info)<-c("gene.id","chr","start","end")

snp_sets <- map_snp_to_gene(snp_info, enhancer_info,extend_start=0,extend_end=0)
map2<-snp_sets$map
map2$locus<-map2$id
map2$peak<-map2$gene.id

dim(map2[which(is.na(map2$gene.start)),])
# 2,411,181 are not in enhancers
# 344901, 12.5% of snps are in enhancers

snps_enh_tested<-unique(map2$peak)
# length 31825; so there are ~32k enhancers with SNPs 
# (of 36216, 87.9%)
snps_in_enh<-map2[which(map2$gene.start>=0),]$id

chrom_sizes<-read.table("LALB_genome_v3_chrom.sizes.txt",sep="\t",header=F)

snps_in_enh_df<-map2[which(map2$gene.start>=0),]
alleles_man_enh<-alleles[which((alleles$allele) %in% snps_in_enh_df$id),]
alleles_man_enh$locus<-alleles_man_enh$allele
alleles_conv<-merge(alleles_man_enh,map2,by="locus",all.x=TRUE,all.y=FALSE)[,c(1,3:20,22,23,24)]

# removing SNPs with no variation
freqs<-alleles_conv[,2:19]
keep<-apply(freqs, 1, function(x) (length(unique(x[!is.na(x)])) != 1))
# length(keep[keep==FALSE]) is 2072, so that's how many will be removed
alleles_conv_filt<-alleles_conv[keep,]
dim(alleles_conv_filt)
# [1] 342829    22
alleles_conv_filt<-alleles_conv_filt[,c(1:19,22)]

### MAKE FUNCTION
find_corr_snps<-function(enhancer) {
  enh_exp<-expr[enhancer,]
  enh_snps<-alleles_conv_filt[alleles_conv_filt$gene.id %in% rownames(enh_exp),]
  enh_snps$gene.id<-NULL
  rownames(enh_snps)<-enh_snps$locus
  enh_snps$locus<-NULL
  
  if(nrow(enh_snps) == 0) return(NULL)
  
  useModel=modelLINEAR
  
  snps = SlicedData$new()
  snps$fileDelimiter = "\t"     # the TAB character
  snps$fileOmitCharacters = "NA" # denote missing values;
  snps$fileSkipRows = 1          # one row of column labels
  snps$fileSkipColumns = 1       # one column of row labels
  snps$fileSliceSize = 2000     # read file in pieces of 2,000 rows
  snps$CreateFromMatrix(as.matrix(enh_snps))
  
  gene = SlicedData$new()
  gene$fileDelimiter = "\t"      # the TAB character
  gene$fileOmitCharacters = "NA" # denote missing values;
  gene$fileSkipRows = 1          # one row of column labels
  gene$fileSkipColumns = 1      # one column of row labels
  gene$fileSliceSize = 2000      # read file in pieces of 2,000 rows
  gene$CreateFromMatrix(as.matrix(enh_exp))
  
  covariates_file_name=character()
  cvrt = SlicedData$new()
  cvrt$fileDelimiter = "\t"      # the TAB character
  cvrt$fileOmitCharacters = "NA" # denote missing values;
  cvrt$fileSkipRows = 1          # one row of column labels
  cvrt$fileSkipColumns = 1      # one column of row labels
  cvrt$fileSliceSize = 2000      # read file in pieces of 2,000 rows
  if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
  }
  
  # Output file name
  output_file_name = tempfile();
  errorCovariance = numeric();
  ###
  
  me = Matrix_eQTL_engine(
    snps=snps,
    gene=gene,
    cvrt=cvrt,
    pvOutputThreshold=1,
    output_file_name="eQTL.out",
    useModel=useModel,
    errorCovariance=errorCovariance,
    verbose=TRUE,
    pvalue.hist=TRUE,
    min.pv.by.genesnp=TRUE,
    noFDRsaveMemory = FALSE)
  
  unlink(output_file_name);
  results<-me$all$eqtls
  write.table(results, "eQTL_output.txt" , append = TRUE, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  
}

# run this once:
lapply(rownames(expr),find_corr_snps)

snps<-read.table("eQTL_output.txt",header=F,sep="\t")
names(snps)<-c("snp","peak","statistic","pvalue","FDR","beta")
snps$uni_fdr<-p.adjust(snps$pvalue,method="fdr")
snps$FDR<-NULL

write.table(snps,"eQTL_output_wFDR.txt",sep="\t",quote=F,row.names=F)

sig_snps<-snps[snps$uni_fdr<0.1,]
dim(sig_snps)
# 7608 sig snps in 4071 enhancers with fdr<0.1
