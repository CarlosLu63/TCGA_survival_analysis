library(survival)
library(dplyr)
library(readxl)
library(data.table)
library(xlsx)

path = "/home/carloslu/Survival"
raw_path = "/home/carloslu/Survival/Rawdata"
result_path = "/home/carloslu/Survival/Cancer_type_result"

#####Prepare TCGA Clinical Data#####
setwd(path)
print("Preparing TCGA Clinical Data...")
TCGA_clinical <- as.data.frame(read_excel(paste0(raw_path, "/TCGA-CDR-SupplementalTableS1.xlsx"), sheet = "TCGA-CDR"))
clinical = TCGA_clinical[c("bcr_patient_barcode","type", "vital_status", "OS.time", "gender", "ajcc_pathologic_tumor_stage")]
colnames(clinical) <- c("sample_ID", "type", "status", "os_time", "gender", "stage")
clinical <- clinical[which(is.na(clinical$status)==FALSE),]
clinical$status[clinical$status=="Dead"] <- 2
clinical$status[clinical$status=="Alive"] <- 1
clinical$gender[clinical$gender=="FEMALE"] <- 2
clinical$gender[clinical$gender=="MALE"] <- 1
clinical$gender = as.numeric(clinical$gender)
clinical$status = as.numeric(clinical$status)
cancer_type = clinical$type[!duplicated(clinical$type)]

#####Prepare cancer-type clinical data#####
all_info <- data.frame()
for (i in cancer_type) {
  print(paste0("Collecting ", i, " clinical data and transcript information..."))
  cancer_samples <- clinical[clinical$type==i,]
  
  #####Prepare cancer-type transcript data#####
  start_time <- proc.time()
  transcript <- as.data.frame(fread(paste0(raw_path, "/TCGA-", i, "/TCGA_", i, "_transcript_exp.csv")))
  colnames(transcript)[1] <- "sample_ID"
  test <- as.data.frame(lapply(transcript[2:ncol(transcript)], as.numeric))
  colsum_0 <- which(colSums(test)==0 | colSums(test==0)>10  | colMeans(test)<1)
  test2 <- test[, -colsum_0]
  remaining_trans <- colnames(test2)
  test3 <- cbind(transcript["sample_ID"], test2)

  #####Merge clinical and transcript data#####
  test4 <- merge(cancer_samples, test3, by = "sample_ID")

  #####Construct Univariate Cox Harzard Model#####
  print(paste0("Constructing Univariate Cox Hazard Model for ", i, " with ", length(remaining_trans), " Transcripts" ))
  #hum.cox <- coxph(Surv(os_time, status) ~ uc010unu.1, data = test4)

  #test_trans = remaining_trans[1:100]
  univ_formulas <- sapply(remaining_trans, function(x) as.formula(paste('Surv(os_time, status)~', x)))

  result_all <- data.frame()
  p=1
  pb <- txtProgressBar(style = 3)
  for (j in univ_formulas) {
    #print(j)
    cox_model <- try(coxph(j, data = test4), silent = TRUE)
    if ('try-error' %in% class(cox_model)) {
      invisible(cox_model)
      next
    } else {
      mus.cox = cox_model
      x <- summary(mus.cox)
      trans_name <- as.character(row.names(x$coefficients))
      p.value <-signif(x$coef[1,5], digits = 6)
      HR <-signif(x$coef[1,2], digits = 6)
      HR.confint.lower <- signif(x$conf.int[1,3], 6)
      HR.confint.upper <- signif(x$conf.int[1,4], 6)
      z_score <- signif(x$coef[1,4], 6)
      res <- as.data.frame(c(trans_name, HR, HR.confint.lower, HR.confint.upper, z_score, p.value))
      colnames(res) <- NULL
      res <- as.data.frame(t(res))
    }
    result_all <- rbind(result_all, res)
    p=p+1
    setTxtProgressBar(pb, p/length(univ_formulas))
  }
  close(pb)
  
  ###Export Z-score file###
  colnames(result_all) <- c("Transcripts", "Hazard_ratio", "HR_lower", "HR_upper", "Z score", "p-value")
  write.csv(result_all, paste0(result_path, "/", i, "_uni_coxph.csv"), row.names=FALSE)

  ###Summary data usage###
  data_info <- as.data.frame(t(c(i, length(cancer_samples$sample_ID), length(remaining_trans))))
  all_info <- rbind(all_info, data_info)

  end_time <- proc.time()
  run_time <- (end_time - start_time)/60
  time <- round(run_time[3][[1]], digits = 2)
  print(paste0("RUN TIME: ", time, " mins"))
}
###Export data usage###
colnames(all_info) <- c("Cancer_type", "Sample_number", "Transcript_number")
write.csv(all_info, paste0(result_path, "/Summary_of_data_usage.csv"), row.names=FALSE)


##########Pickup significant transcripts in each cancer type##########
print("Pickup significant transcripts in each cancer type...")
library(biomaRt)

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
setwd(result_path)
all_overlap <- data.frame()
q=1
bar <- txtProgressBar(style = 3)
for (i in cancer_type) {
  trans_result <- read.csv(paste0(i, "_uni_coxph.csv"))
  colnames(trans_result) <- c("UCSC_transcript_id", "Hazard_ratio", "HR_lower", "HR_upper", "Z.score", "p.value")

  ###filter for |Z| >= 1.96###
  trans_sig <-trans_result %>% filter(Z.score >= 1.96 | Z.score <= (-1.96))
    
  ###Convert significant transcripts UCSC ID to Ensembl ID###
  merged_gene <- merge(convert_list2, trans_sig, by = "UCSC_transcript_id")
  merged_gene <- merged_gene[order(merged_gene$p.value, decreasing=FALSE),]
  write.xlsx2(merged_gene, file=paste(path, "/Significant_transcripts/", i, "_sig_transcripts.xlsx", sep=""), sheetName = "Sig_trans_list", append=TRUE, row.names = FALSE)
  
  ###Merge NHRI Target kinases list###
  kinase_merge <- merge(kinase, merged_gene, by = "hgnc_symbol")
  write.xlsx2(kinase_merge, file=paste(path, "/Significant_transcripts/", i, "_sig_transcripts.xlsx", sep=""), sheetName = "overlap_target_list", append=TRUE, row.names = FALSE)  

  catch_result = kinase_merge
  catch_result$cancer_type = rep(i, length(catch_result[,"hgnc_symbol"]))
  all_overlap <- rbind(all_overlap, catch_result)

  q=q+1
  setTxtProgressBar(bar, q/length(cancer_type))
}
write.xlsx2(all_overlap, paste(path, "/Significant_transcripts/Overlap_targets.xlsx", sep=""), sheetName = "overlap_list (|Z|>=1.96)", append=TRUE, row.names = FALSE)
close(bar)
print("All the process is finished!!")
