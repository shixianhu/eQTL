setwd("Desktop")
##Load the epigenetic features information from Ensembl-derived annotation (recovered by BioMart, Ensembl 92 Regulation Database)
Epi_annotation <- read.table(file="test_run_QTLs/Epigenetic_annotation_Ensembl_hg19.txt", stringsAsFactors=F, sep="\t", header=T)

##Load the Ensembl annotation for the target eQTL variants, downloaded with VEP (GRCh37/hg19 or GRCh38 as needed)
eQTL_VEP <- read.table(file="eQTL_variants_VEP_report.txt", stringsAsFactors=F, sep="\t", header=T)

##Remove possible duplicate reads, due to VEP inherent redundancy
eQTL_VEP <- unique(eQTL_VEP)


##Main code; Check for overlap with regulatory elements
final_report <- NULL

for (x in 1:dim(eQTL_VEP)[1]){
  SNP <- eQTL_VEP$Uploaded_variation[x]
  varChrom <- unlist(strsplit(eQTL_VEP$Location[x],split=":"))[1]
  varDist <- unlist(strsplit(eQTL_VEP$Location[x],split=":"))[2]
  varStart <- as.numeric(unlist(strsplit(varDist,split="-"))[1])
  varFinish <- as.numeric(unlist(strsplit(varDist,split="-"))[2])
  chrom_ann <- Epi_annotation[Epi_annotation$Chromosome.scaffold.name == varChrom,]
  for (y in 1:dim(chrom_ann)[1]){
    if(varStart==varFinish){
      if(varStart >= chrom_ann$Start.bp.[y] && varStart <= chrom_ann$End.bp.[y]){
        record <- cbind(SNP,chrom_ann[y,])
        final_report <- rbind(final_report,record) ##This retrieves the regulatory element info and adds the SNP rsID for recognition
      }
    } else {if((varStart >= chrom_ann$Start.bp.[y] && varStart <= chrom_ann$End.bp.[y]) || 
               (varFinish >= chrom_ann$Start.bp.[y] && varFinish <= chrom_ann$End.bp.[y])){
                 record <- cbind(SNP,chrom_ann[y,])
                 final_report <- rbind(final_report,record) ##This retrieves the regulatory element info and adds the SNP rsID for recognition
        }
      }
  }  
}


##Save the complete report in a text file
write.table(final_report, file="test_run_QTLs/eQTL_variants_epi_overlap.txt",sep="\t", row.names=F, quote=F)
