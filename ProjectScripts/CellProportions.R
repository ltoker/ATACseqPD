source("/data/Rprojects/GeneralScripts/generalFunc.R")
source("ProjectScripts/ProjectFunctions.R")
plotMA = DESeq2::plotMA

Count2CPM <- function(countData){
  apply(countData, 2, function(smp){
    TotalCount = sum(smp)
    (10^6)*smp/TotalCount
  })
}

annoFileCollapsed <- GetGenomeAnno(genome = "hg19")
NarrowPeak <-  read.table("CellTypesATAC/data/dorsolateral_prefrontal_cortex.narrowPeak.gz", header = F, sep = "\t")
names(NarrowPeak) <- c("CHR", "START", "END", "PeakName", "Score", "Srand", "Enrichment", "pPvalue", "pQvalue", "OffsetFromStart")

NarrowPeak$CHR <- sapply(as.character(NarrowPeak$CHR), function(x){
  paste0("chr", x) 
})

NarrowPeak %<>% as(., "GRanges")

NarrowPeak <-  mergeByOverlaps(annoFileCollapsed, NarrowPeak, maxgap = 0, type = "any", select = "all") %>% data.frame()
NarrowPeak %<>% filter(pQvalue > 7) %>% .[!grepl("MT|random|GL|hs", .$NarrowPeak.seqnames), ] %>% droplevels()


HTseqCountsCellTypes <- read.table("CellTypesATAC/data/dorsolateral_prefrontal_cortex_counts.tsv.gz", header = T, sep = "\t")
HTseqCountsCellTypes %<>% filter(Geneid %in% NarrowPeak$PeakName)

names(HTseqCountsCellTypes) <- sapply(names(HTseqCountsCellTypes) , function(x){
  x <- gsub("X.data3.public_samples.fullard_2018.results.bamfiles.", "", x)
  strsplit(x, "_")[[1]][1]
})

MetaCellTypes <- read.table("CellTypesATAC/meta/atac_Cortical_sample_ids.csv", header = T, sep = "\t")

#Filter the relevant samples from the count matrix

HTseqCountsCellTypes %<>% select(c("Geneid", as.character(MetaCellTypes$SampleID)))
countMatrix <- as.matrix(HTseqCountsCellTypes[,-1])
rownames(countMatrix) <- as.character(HTseqCountsCellTypes$Geneid)

MetaCellTypes$RiP <- apply(HTseqCountsCellTypes[,-1], 2, sum)
MetaCellTypes$CellType <- relevel(MetaCellTypes$CellType, ref = "Non_Neuronal")

Model = as.formula("~CellType")

DESeqDS <- DESeqDataSetFromMatrix(countData = countMatrix, colData = MetaCellTypes, design = Model)
DESeqOut <- DESeq(DESeqDS)

DESeqOutCellTypeResults <- results(DESeqOut, name = "CellType_Neuronal_vs_Non_Neuronal", independentFiltering = T, alpha = 0.05) %>% data.frame %>% mutate(PeakName = rownames(.))
DESeqOutCellTypeResultsAnno <- merge(DESeqOutCellTypeResults, NarrowPeak, by = "PeakName")

CellTypePeaks <- DESeqOutCellTypeResultsAnno %>% filter(abs(log2FoldChange) > 2, padj < 0.05,  baseMean > 100)
write.table(CellTypePeaks, "CellTypesATAC/CellTypesATAC.tsv", sep = "\t", row.names = F, col.names = T)



