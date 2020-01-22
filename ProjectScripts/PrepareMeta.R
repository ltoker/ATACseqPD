Metadata <- read.table("meta/sample_list.csv", header = T, sep = "\t")
Meta2 <- read.table("meta/atac_sample_ids.csv", header = T, sep = "\t")

Meta2$activemotif_id <- sapply(as.character(Meta2$sample_id), function(x){
  strsplit(x, "-")[[1]][2]
})

Metadata <- merge(Metadata %>% select(-sample_id), Meta2 %>% select(-matches("cohort")), by = "activemotif_id")
Metadata$SampleID <- sapply(make.names(Metadata$sample_id), function(x){
  strsplit(x, "_")[[1]][2]
})

Metadata$Cohort <- sapply(Metadata$origin, function(x){
  if(x == "PV"){
    "PW"
  } else {
    as.character(x)
  }
})


names(Metadata) <- sapply(names(Metadata), function(x){
  if(grepl("age", x)){
    "Age"
  } else if(grepl("sex", x)){
    "Sex"
  } else if(grepl("pm_time", x)){
    "PMI"
  } else if (grepl("rin", x)){
    "RIN"
  } else {
    x
  }
})

Metadata$Group <- sapply(Metadata$condition, function(x){
  if(x == "Case"){
    "PD"
  } else if(x == "Control") {
    "Cont"
  }
}) %>% factor(levels  = c("Cont", "PD"))