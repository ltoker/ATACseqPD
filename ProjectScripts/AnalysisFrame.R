source("generalFunc.R")
source("ProjectScripts/ProjectFunctions.R")
packageF("ggExtra")

Cohort = "PW"
cohort = Cohort
ResultsPath = paste0(Cohort, "_Results")

PeakLoc = "data/peaks/"
CountMatrixLoc = paste0("data/", Cohort, "_counts.tsv.gz")
CellTypePeakCountLoc = "data/NeuN_narrow_peak_counts.tsv"


if(!ResultsPath %in% list.dirs(full.names = F, recursive = F)){
  dir.create(ResultsPath)
}
ResultsPath = paste0(ResultsPath, "/")

plotMA = DESeq2::plotMA

annoFileCollapsed <- GetGenomeAnno(genome = "hg19")

source("ProjectScripts/PrepareMeta.R")

#Keep only the relevant samples

Metadata %<>% filter(Cohort == cohort) %>% droplevels



ContNarrowPeak <- read.table(paste0("data/peaks/Control_", Cohort, ".narrowPeak.gz"), header = F, sep = "\t") 
names(ContNarrowPeak) <-  c("CHR", "START", "END", "PeakName", "Score", "Srand", "Enrichment", "pPvalue", "pQvalue", "OffsetFromStart")
ContNarrowPeak %<>% filter(pPvalue > 7)  %>% mutate(Width = END-START) %>% .[!grepl("GL|hs|MT", .$CHR),]

PDNarrowPeak <- read.table(paste0("data/peaks/Case_", Cohort, ".narrowPeak.gz"), header = F, sep = "\t") 
names(PDNarrowPeak) <-  c("CHR", "START", "END", "PeakName", "Score", "Srand", "Enrichment", "pPvalue", "pQvalue", "OffsetFromStart")
PDNarrowPeak %<>% filter(pPvalue > 7) %>% mutate(Width = END-START) %>% .[!grepl("GL|hs|MT", .$CHR),]


########## Repeat for all samples #############

InputPeakAllCalled <- list()
InputPeakAllCalled$PeakData <- read.table(paste0(PeakLoc, Cohort, ".narrowPeak.gz"), header = F, sep = "\t")
names(InputPeakAllCalled$PeakData) <- c("CHR", "START", "END", "PeakName", "Score", "Srand", "Enrichment", "pPvalue", "pQvalue", "OffsetFromStart")[1:ncol(InputPeakAllCalled$PeakData)]
InputPeakAllCalled$PeakData %<>% mutate(pValue = 10^(-pPvalue), width = END-START) %>% filter(pValue < 10^(-7))
InputPeakAllCalled$PeakData <- InputPeakAllCalled$PeakData[!grepl("GL|hs", InputPeakAllCalled$PeakData$CHR),] %>% droplevels()

#Remove blacklisted peaks 
BlackListed <- read.table("data/H3K27Ac_black_list.bed", header = F, sep = "\t")
BlackListed %<>% mutate(Peak.Location = paste0("chr", V1, ":", V2, "-", V3)) 
names(BlackListed)[1:3] <- c("CHR", "START", "END")


#Filter the blacklisted peaks from the peak file
blackListedPeaks <- findOverlaps(query = InputPeakAllCalled$PeakData %>% as(., "GRanges"), subject = BlackListed %>% as(., "GRanges"), maxgap = 0, type = "any", select = "all")

InputPeakAllCalled$PeakData <- InputPeakAllCalled$PeakData[-queryHits(blackListedPeaks),]



InputPeakAllCalled$Summary <- data.frame(TotalPeaks = InputPeakAllCalled$PeakData %>% data.frame %>% nrow,
                                         TotalCoverage = InputPeakAllCalled$PeakData %>% data.frame %>% .$width %>% sum)

InputPeakAllCalled$Summary %<>% mutate(TotalCovPercent = signif(100*TotalCoverage/(2.7*10^9), digits = 3))
InputPeakAllCalled$PeakData %>% data.frame() %>% .$width %>% summary()                      

##################### HTseq counts ##########
HTseqCounts <- read.table(CountMatrixLoc, header = T, sep = "\t")
names(HTseqCounts) <- sapply(names(HTseqCounts), function(x){
  if(grepl("_", x)){
    strsplit(x, "_")[[1]][2]
  } else {
    x
  }
})

HTseqCounts %<>% mutate(Peak.Location = paste0("chr", Chr, ":", Start, "-", End))
names(HTseqCounts)[2:4] <- c("CHR", "START", "END")

#Arrange the samples to match Metadata order
HTseqCounts %<>% select(c("Geneid", "CHR", "START", "END", "Strand",  "Length", as.character(Metadata$SampleID)))



# Keep only filtered peaks 
HTseqCounts <- HTseqCounts %>% filter(Geneid %in% InputPeakAllCalled$PeakData$PeakName) %>% droplevels()
AllReads <- apply(HTseqCounts %>% select(matches("Sample")), 2, sum)

#Remove MT reads

MTpeaks <- HTseqCounts %>% filter(CHR == "MT")
MTpeakRead <- apply(MTpeaks %>% select(matches("Sample")), 2, sum)

HTseqCounts %<>% filter(CHR != "MT")
HTseqCounts
  
#Merge the replicates
HTseqCountsOrg <- HTseqCounts
MetadataOrg <- Metadata

Replicates <- sapply(unique(Metadata$activemotif_id), function(x){
  temp <- c(x, Metadata %>% filter(activemotif_id == x) %>% .$SampleID %>% as.character()) %>% t %>% data.frame()
  names(temp) <- c("activemotif_id", "Rep1", "Rep2")
  temp
}, simplify = F) %>% rbindlist

Metadata %<>% filter(!duplicated(activemotif_id))
Metadata$SampleID <- sapply(as.character(Metadata$SampleID), function(x){
  gsub("..$", "", x)
}) 

Metadata <- merge(Metadata, Replicates, by = "activemotif_id")

HTseqCounts <- sapply(Metadata$SampleID, function(subj){
  Rep <- Metadata %>% filter(SampleID == subj) %>% select(Rep1, Rep2) %>% unlist %>% as.character()
  HTseqCounts %>% select(Rep) %>% apply(1, sum)
})

HTseqCounts <- cbind(HTseqCountsOrg %>% select(-matches("Sample")), HTseqCounts)

HTseqCounts$MaxCount = apply(HTseqCounts %>% select(matches("Sample")), 1, max)
HTseqCounts %<>% mutate(NormalizedMaxCount = MaxCount/Length)  


AllCalledData <- GetCountMatrixHTseq(HTseqCounts, OtherNormRegEx = "^C1orf43_|^CHMP2A_|^EMC7_|^GPI_|^PSMB2_|^PSMB4_|^RAB7A_|^REEP5_|^SNRPD3_|^VCP_|^VPS29",
                                     MetaSamleCol = "SampleID", countSampleRegEx = "Sample", MetaCol = c("SampleID", "Sex", "Age", "Age", "PMI", "RIN", "Cohort", "Group", "Rep1", "Rep2"))

closeDev()

#Get cellular estimates 
AllCalledData$SampleInfo <- GetCellularProportions(AllCalledData$SampleInfo, MetaSamplCol = "SampleID")

pdf(paste0(ResultsPath, "SampleCorAllPeaks", Cohort, ".pdf"), useDingbats = F, width = 10, height = 8)
countMatrixFullAllCalled <- GetCollapsedMatrix(countsMatrixAnnot = AllCalledData$countsMatrixAnnot %>% filter(!duplicated(.$PeakName)), collapseBy = "PeakName",CorMethod = "pearson",countSampleRegEx = "Samp",MetaSamleCol = "SampleID", MetaSamleIDCol = "SampleID",
                                               FilterBy = "", meta = AllCalledData$SampleInfo, title = paste0("Sample correlation, ", Cohort))
closeDev()


countMatrixDF <- AllCalledData$countsMatrixAnnot %>% filter(!duplicated(.$PeakName)) %>% data.frame %>% select(matches("Peak|Sample"))
countMatrixDF$MedianCount <- apply(countMatrixDF %>% select(matches("Sample")), 1, mean)
countMatrixDF %<>% mutate(NormCount = 200*MedianCount/Peak.width)
countMatrixDF$baseMean <- apply(countMatrixDF %>% select(matches("Sample")), 1, mean)

countMatrix_filtered <- countMatrixDF %>% filter(NormCount > 5) %>% select(matches("Sample")) %>% as.matrix()
rownames(countMatrix_filtered) <- as.character(countMatrixDF %>% filter(NormCount > 5) %>% .$PeakName)

Model = as.formula(~Group + Sex + Age + Oligo_MSP + NeuNall_MSP + MTprop)
DESeqOutAll_Full <- RunDESeq(data = countMatrix_filtered, UseModelMatrix = T, MetaSamleCol = "SampleID",SampleNameCol = "SampleID",
                             meta = countMatrixFullAllCalled$Metadata, normFactor = "MeanRatioAll",
                             FullModel = Model, test = "Wald", FitType = "local")

DESegResultsSex_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "SexM") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsAge <- GetDESeqResults(DESeqOutAll_Full, coef = "Age") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsGroup <- GetDESeqResults(DESeqOutAll_Full, coef = "GroupPD") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")

ggMarginal(data = DESegResultsGroup %>% filter(!duplicated(PeakName)) %>% mutate(LogBaseMean = log(baseMean)), x = "LogBaseMean", y ="log2FoldChange", type = "densigram", margins = "both")

#Repeat with RLE normalization
DESeqOutAll_FullRLE <- RunDESeq(data = countMatrix_filtered, UseModelMatrix = T, MetaSamleCol = "SampleID",SampleNameCol = "SampleID",
                                meta = countMatrixFullAllCalled$Metadata,
                                FullModel = Model, test = "Wald", FitType = "local")
                            


DESegResultsSex_FullAllRLE <- GetDESeqResults(DESeqOutAll_FullRLE, coef = "SexM") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsAgeRLE <- GetDESeqResults(DESeqOutAll_FullRLE, coef = "Age") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsGroupRLE <- GetDESeqResults(DESeqOutAll_Full, coef = "GroupPD") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")

ggMarginal(data = DESegResultsGroupRLE %>% filter(!duplicated(PeakName)), x = "baseMean", y = "log2FoldChange", type = "densigram", margins = "both")

save.image(paste0(ResultsPath, Cohort, ".Rdata"))
saveRDS(DESeqOutAll_Full, file = paste0(ResultsPath, Cohort, "DEoutput.Rds"))
saveRDS(DESeqOutAll_FullRLE, file = paste0(ResultsPath, Cohort, "DEoutputRLE.Rds"))
saveRDS(DESegResultsGroup, file = paste0(ResultsPath, Cohort, "DEresults.Rds"))
saveRDS(DESegResultsGroupRLE, file = paste0(ResultsPath, Cohort, "DEresultsRLE.Rds"))