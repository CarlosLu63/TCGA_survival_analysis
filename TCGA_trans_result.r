rm(list = ls())
library(tidyr)
library(biomaRt)
library(xlsx)

path = "/home/carloslu/Survival"
raw_path = "/home/carloslu/Survival/Rawdata"
result_path = "/home/carloslu/Survival/Cancer_type_result"

cancer_type = clinical$type[!duplicated(clinical$type)]

###Target kinases list###
kinase <- read.table(paste0(raw_path, "/target_kinase.txt"))
colnames(kinase) <- c("hgnc_symbol")

#####Prepare UCSC ID conversion list and corresponding ensembl ID#####
convert_list <- read.table(paste0(raw_path, "/knownToEnsembl.txt"))
colnames(convert_list) <- c("UCSC_transcript_id", "ensembl_transcript_id")

mart <- useMart("ensembl", "hsapiens_gene_ensembl")
ensembl2gene <- getBM(attributes=c("ensembl_transcript_id", "hgnc_symbol", "ensembl_gene_id"),
                       filters = "ensembl_transcript_id",
                       values = convert_list[,"ensembl_transcript_id"], 
                       mart = mart)
convert_list2 <- merge(convert_list, ensembl2gene, by = "ensembl_transcript_id")

#####Read Z-score files separately#####
i="ACC"
setwd(result_path)
trans_result <- read.csv(paste0(i, "_uni_coxph.csv"))
colnames(trans_result) <- c("UCSC_transcript_id", "Hazard_ratio", "HR_lower", "HR_upper", "Z.score", "p.value")

###filter for |Z| >= 5###
trans_sig <-trans_result %>% filter(Z.score > 5 | Z.score < (-5))

###Convert significant transcripts UCSC ID to Ensembl ID###
merged_gene <- merge(convert_list2, trans_sig, by = "UCSC_transcript_id")
merged_gene <- merged_gene[order(merged_gene$p.value, decreasing=FALSE),]
#write.csv(merged_gene, paste0(path, "/Significant_transcripts/", i, "_sig_transcripts.csv"), row.names=FALSE)
write.xlsx2(merged_gene, file=paste(path, "/Significant_transcripts/", i, "_sig_transcripts.xlsx", sep=""), sheetName = "Sig_trans_list", append=TRUE, row.names = FALSE)  

###Merge NHRI Target kinases list###
kinase_merge <- merge(kinase, merged_gene, by = "hgnc_symbol")
write.xlsx2(kinase_merge, file=paste(path, "/Significant_transcripts/", i, "_sig_transcripts.xlsx", sep=""), sheetName = "overlap_target_list", append=TRUE, row.names = FALSE)  







###################################################################################################################

convert_list <- read.table(paste0(raw_path, "/knownToEnsembl.txt"))
colnames(convert_list) <- c("Transcripts", "Ensembl_ID")
setwd(result_path)

#####Read Z-score files separately#####
for (i in cancer_type) {
    trans_result <- read.csv(paste0(i, "_uni_coxph.csv"))

    ###filter for |Z| >= 5###
    trans_sig <-trans_result %>% filter(Z.score > 5 | Z.score < (-5))
    
    ###Convert significant transcripts UCSC ID to Ensembl ID###
    merged_gene <- merge(convert_list, trans_sig, by = "Transcripts")
    merged_gene <- merged_gene[order(merged_gene$p.value, decreasing=FALSE),]
    write.csv(merged_gene, paste0(path, "/Significant_transcripts/", i, "_sig_transcripts.csv"), row.names=FALSE)
}
