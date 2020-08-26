# pre-process data #############################################################
read_and_select_lumi <- function(microarray_file, metadata_file) {
  lumibatch <- lumiR(fileName = microarray_file,
    detectionTh = 0.05,
    annotationColumn = c('ENTREZ_GENE_ID','ACCESSION', 'SYMBOL', 
      'PROBE_SEQUENCE', 'PROBE_START', 'CHROMOSOME', 'PROBE_CHR_ORIENTATION', 
      'PROBE_COORDINATES', 'DEFINITION'))

  pData <- readr::read_csv(metadata_file, col_types = cols(
      MicroarrayId = col_character(),
      StudyId = col_double(),
      Date = col_date(format = ""),
      Index = col_character(),
      RIN = col_double(),
      Enrollment = col_date(format = ""),
      Status = col_character(),
      Hsv2 = col_character(),
      Progesterone = col_double(),
      BV = col_character(),
      Visit = col_character(),
      Drop = col_character()
  ))
  pData <- as.data.frame(pData)
  rownames(pData) <- pData$MicroarrayId
  pData(lumibatch) <- pData
    
  # remove unwanted samples    
  remove <- str_detect(pData(lumibatch)$Drop, "Drop for microarray|all")
  
  lumibatch[ , !remove]
}

batch1 <- read_and_select_lumi(microarray_batch_1_file, 
    microarray_batch_1_metadata)

batch2 <- read_and_select_lumi(microarray_batch_2_file,
    microarray_batch_2_metadata)

# Create column indicating which batch samples came from
pData(batch1)$Batch <- "1"
pData(batch2)$Batch <- "2"

# convert Index column to character so it matches in both sets
pData(batch1)$Index <- as.character(pData(batch1)$Index)

# combine and normalize
combined <- lumi::combine(batch1, batch2) %>% 
    lumiT(method = "log2") %>% 
    lumiN(method = "quantile")

# keep probes detected in at least 10 samples
filtered <- combined[rowSums(detection(combined) < 0.05) > 10 , ]

fix_probe <- function(filtered, old_target_id, new_target_id, new_entrez_id) {
    fData(filtered)$ENTREZ_GENE_ID[fData(filtered)$TargetID == old_target_id] <- new_entrez_id
    fData(filtered)$SYMBOL[fData(filtered)$TargetID == old_target_id] <- new_target_id
    fData(filtered)$TargetID[fData(filtered)$TargetID == old_target_id] <- new_target_id
    filtered
}

# Manual correction for some probes that are annotated incorrectly
# the probe for HS.483183 maps to TNFAIP8, entrez 25816
filtered <- fix_probe(filtered, "HS.483183", "TNFAIP8", 25816)
# the probe for LOC388312 maps to LOC100132062, entrez 100132062
filtered <- fix_probe(filtered, "LOC388312", "LOC100132062", 100132062)
# the probe for LOC653887 maps to ZSWIM6, entrez 57688
filtered <- fix_probe(filtered, "LOC653887", "ZSWIM6", 57688)

save(filtered, file = "data/clean/CombinedBatches_lumibatch.Rda")

# save expression levels as a data frame ####################################################
# file is too big with all the fData so subset to key columns
genes <- fData(filtered) %>% 
    select(ProbeId = ProbeID, Gene = TargetID, EntrezId = ENTREZ_GENE_ID) %>% 
    canonicalize_gene_symbols()

# create participant level expression df
participant_level_probes <- exprs(filtered) %>% 
    reshape2::melt(value.name = "Expression", varnames = c("ProbeId", "MicroarrayId")) %>% 
    mutate(ProbeId = as.character(ProbeId), MicroarrayId = as.character(MicroarrayId)) %>% 
    left_join(pData(filtered), by = "MicroarrayId") %>%
    left_join(genes, by = "ProbeId")

save(participant_level_probes, file = "data/clean/participant_level_probes.Rda")
