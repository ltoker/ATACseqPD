SampleData <- sapply(list.files("ReadCounts/"), function(x){
  temp <- read.table(paste0("ReadCounts/", x), header = F, quote = "'", comment.char = "!", sep = "\t")
  temp2 <- temp$V1[grepl("read1", temp$V1)] %>% sapply(function(Line){
    strsplit(as.character(Line), " ")[[1]][1]
  }) %>% as.numeric() %>% t %>% data.frame() 
  names(temp2) <- c("PairsAll", "PairsNoMM", "PairsNoMM_NoDuplicate")
  temp2
}, simplify = F) %>% rbindlist()

SampleData$SampleID <- sapply(list.files("ReadCounts/"), function(x){
  strsplit(x, "_")[[1]][3] %>% gsub(".txt", "", .) %>% make.names()
})

Metadata2 <- merge(AllCalledData$SampleInfo, SampleData, by = "SampleID")
Metadata2 %<>% mutate(RiPraw = c(apply(HTseqCountsRaw %>% select(matches("Sample")), 2, sum)), FRiP = RiPraw/PairsNoMM_NoDuplicate, FRip2 = TotalCount/PairsNoMM_NoDuplicate) %>% arrange(FRiP)




MetadataBoth$SubjectID <- sapply(as.character(MetadataBoth$SampleID), function(x){
  gsub("..$", "", x)
})
MetadataBothMerged <- MetadataBoth %>% group_by(SubjectID) %>% summarise(PairsNoMM_NoDuplicate = sum(PairsNoMM_NoDuplicate),
                                                                          AllRiP = sum(AllRiP), MTRiP = sum(MTRiP)) %>% mutate(FRiP = AllRiP/PairsNoMM_NoDuplicate, MTprop = MTRiP/AllRiP)
MetadataBothMerged$AMid = sapply(as.character(MetadataBothMerged$SubjectID), function(x){
  strsplit(x, "\\.")[[1]][[2]] 
}) %>% as.numeric()

MetadataBothMerged %<>% arrange(AMid)
