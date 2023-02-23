rm(list = ls())
library(TCGAbiolinks)
library(readxl)
library(tidyr)

raw_path = "/home/carloslu/Survival/Rawdata"

TCGA_samples <- as.data.frame(read_excel("/home/carloslu/Survival/Rawdata/TCGA-CDR-SupplementalTableS1.xlsx", sheet="TCGA-CDR", col_names = TRUE))
TCGA_samples = TCGA_samples[c("bcr_patient_barcode","type")]
cancer_type = TCGA_samples$type[!duplicated(TCGA_samples$type)]
test_type = cancer_type[-14] #LAMLL

setwd(raw_path)
for (i in test_type){
    print(paste("Collecting samples for", i))
    cancer_samples <- TCGA_samples[TCGA_samples$type==i,]
    cancer_samples = cancer_samples$bcr_patient_barcode
    
    query <- GDCquery(project = paste0("TCGA-", i),
                data.category = "Gene expression",
                data.type = "Isoform expression quantification",
                experimental.strategy = "RNA-Seq",
                sample.type = "Primary Tumor",
                platform = "Illumina HiSeq",
                barcode = cancer_samples,
                file.type  = "normalized_results",
                legacy = TRUE)
                
    print(paste("Downloading transcript data for", i, "samples"))
    GDCdownload(query, files.per.chunk = 100, directory = raw_path)

    #####changing filename#####
    ###Convert Nested List to a DataFrame###
    print("Changing file name and move them to cancer-type file...")
    query_results =  query$results
    query_results <- as.data.frame(do.call(cbind, query_results))
    dir_name = query_results$file_id
    file_name = query_results$file_name
    id = query_results$cases.submitter_id
    
    rename_files <- list(paste0(raw_path, "/TCGA-", i, "/legacy/Gene_expression/Isoform_expression_quantification/"
    , dir_name, "/", file_name))
    rename_files <- as.vector(unlist(rename_files))

    after_named <- list(paste0(raw_path, "/TCGA-", i, "/legacy/Gene_expression/Isoform_expression_quantification/",
    dir_name, "/", id, "_isoforms_quantification.txt"))
    after_named <- as.vector(unlist(after_named))
    file.rename(rename_files, after_named)

    ###move transcript mdata to cancer-type file###
    move_use <- as.data.frame(after_named)
    for (j in move_use$after_named){
        move_data = paste0("mv ", j, " -t ", raw_path, "/TCGA-", i)
        system(move_data)
    }
    ###merge all sample transcripts together###
    cancer_path = paste0("/home/carloslu/Survival/Rawdata/", "/TCGA-", i)
    setwd(cancer_path)
    unlink("legacy",recursive = TRUE)
    unlink("MANIFEST.txt")


    print(paste("Merging sample transcripts together for", i))
    isoform_all = NULL
    for (k in list.files(cancer_path)){
        sample_iso = read.table(k, header = TRUE)
        trans_name = sample_iso["isoform_id"]
        ids = substr(k, 1, 12)
        isoform_count = sample_iso["normalized_count"] 
        colnames(isoform_count) <- ids
        isoform_all = as.data.frame(cbind(isoform_all, as.matrix(isoform_count)))
    }
    row.names(isoform_all) <- as.matrix(trans_name)

    t_isoform_all <- as.data.frame(t(isoform_all))
    write.csv(t_isoform_all, paste0("TCGA_",i, "_transcript_exp.csv"), row.names=TRUE)
    print(paste(i, "is Done!!!"))
}
##############################################################################################
cancer_samples <- TCGA_samples[TCGA_samples$type=="LAML",]
cancer_samples = cancer_samples$bcr_patient_barcode

#####collecting data by giving sample IDs#####
query <- GDCquery(project = "TCGA-LAML",
                data.category = "Gene expression",
                data.type = "Isoform expression quantification",
                experimental.strategy = "RNA-Seq",
                platform = "Illumina HiSeq",
                sample.type = "Primary Blood Derived Cancer - Peripheral Blood",
                barcode = cancer_samples,
                file.type  = "normalized_results",
                legacy = TRUE)
#print(paste("Downloading", cancer_type, "sample transcript data..."))
GDCdownload(query, files.per.chunk = 100, directory = raw_path)

#####changing filename#####
###Convert Nested List to a DataFrame###
query_results =  query$results
query_results <- as.data.frame(do.call(cbind, query_results))
dir_name = query_results$file_id
file_name = query_results$file_name
id = query_results$cases.submitter_id

rename_files <- list(paste0(raw_path, "/TCGA-", i, "/legacy/Gene_expression/Isoform_expression_quantification/"
, dir_name, "/", file_name))
rename_files <- as.vector(unlist(rename_files))

after_named <- list(paste0(raw_path, "/TCGA-", i, "/legacy/Gene_expression/Isoform_expression_quantification/",
dir_name, "/", id, "_isoforms_quantification.txt"))
after_named <- as.vector(unlist(after_named))
file.rename(rename_files, after_named)

###move isoformdata to cancer-type file###
move_use <- as.data.frame(after_named)
for (j in move_use$after_named){
    move_data = paste0("mv ", j, " -t ", raw_path, "/TCGA-", i)
    system(move_data)
}

###merge all sample transcripts together###
cancer_path = paste0("/home/carloslu/Survival/Rawdata/", "/TCGA-", i)
setwd(cancer_path)
unlink("legacy",recursive = TRUE)

print(paste("Merging sample transcripts together for", cancer_type))
isoform_all = NULL
for (i in list.files(cancer_path)){
    sample_iso = read.table(i, header = TRUE)
    trans_name = sample_iso["isoform_id"]
    ids = substr(i, 1, 12)
    isoform_count = sample_iso["normalized_count"] 
    colnames(isoform_count) <- ids
    isoform_all = as.data.frame(cbind(isoform_all, as.matrix(isoform_count)))
}
row.names(isoform_all) <- as.matrix(trans_name)

t_isoform_all <- as.data.frame(t(isoform_all))
write.csv(t_isoform_all, paste0("TCGA_",cancer_type, "_transcript_exp.csv"), row.names=TRUE)
#isoform_all <- cbind(trans_name, isoform_all)
#isoform_all_id_converted <- merge(isoform_all, ID_convert, by = "isoform_id")
