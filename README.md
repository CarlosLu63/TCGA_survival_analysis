# TCGA_survival_analysis
TCGA survival analysis in transcript level  
Created: 2023.02.17  

## Environment and Version
Packages: TCGAbiolinks 3.16, survival 2.25.3, biomaRt 2.54.0  

## TCGAbiolinks docker
先下載 image，再進入 TCGAbiolinks container. 裡面已經有 R 以及 TCGAbiolinks function，需再額外安裝 survival package  
docker pull tiagochst/tcgabiolinksgui  
image id: tiagochst/tcgabiolinksgui  
References:  
> TCGAbiolinks: https://github.com/BioinformaticsFMRP/TCGAbiolinks  
> Survival R package: https://rdrr.io/cran/survival/

## 分析資料來源
### Clinical data  
Source: PanCanAtlas Publications  
Only primary tumor specimens (TCGA barcode: 01) were analyzed.  
For the analysis of LAML (Acute Myeloid Leukemia, 急性骨髓性白血病) data, blood cancer specimens were given the TCGA barcode ‘‘03’’ and cancers with this code were included.  
### Transcript data  
Source: Legacy Archive  
Category: Gene expression, Isoform expression quantification.  
Sample selection: primary tumor specimens (TCGA barcode: 01) and TCGA barcode ‘‘03” for LAML.  
File content: transcript_id (encoded by UCSC ID), normalized_count (RSEM normalization method, normalized by transcript length)  

## Notion for more information  
> https://www.notion.so/TCGA-Survival-Analysis-in-Transcript-Level-Feb-2023-20f8a793d25d4a09bdb3e938db22c067?pvs=4  
