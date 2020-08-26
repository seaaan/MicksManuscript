# metadata for microarray
# files from MicksSampleLog/data/clean
get_metadata <- function(batch) {
  if (is.numeric(batch)) {
      f <- paste0("data/required_files/Batch", batch, "_metadata.csv")
      read_csv(f, col_types = 
              cols(
                  MicroarrayId = col_character(),
                  StudyId = col_integer(),
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
  } else {
      f <- paste0("data/required_files/ddpcr_proteomics_metadata.csv")
      read_csv(f, col_types = 
          cols(
              StudyId = col_double(),
              Enrollment = col_date(format = ""),
              Date = col_date(format = ""),
              Index = col_character(),
              Visit = col_character(),
              Progesterone = col_double(),
              Hsv2 = col_character(),
              BV = col_character(),
              RIN = col_double(),
              Drop = col_character()
          ))
  }
  
}

get_filtered_expression_set <- function() {
    load("data/clean/CombinedBatches_lumibatch.Rda")
    filtered
}

get_participant_level_probe_expression <- function() {
    load("data/clean/participant_level_probes.Rda")
    participant_level_probes
}

# from http://bioinf.wehi.edu.au/software/MSigDB/
get_hallmark_gene_sets <- function() {
    load("data/clean/human_H_v5p2.rdata")
    Hs.H
}

get_hallmark_gene_set_names_and_categories <- function() {
    read_csv("data/clean/hallmark_gene_sets_renamed.csv", col_types = 
            cols(
                Name = col_character(),
                Category = col_character(),
                ShortName = col_character()
            ))
}

# look up what the file name should be based on the coefficient of interest
file_lookup <- function(coef) {
    original_coef <- coef
    coef <- str_to_lower(coef)
    
    result <- case_when(
        coef %in% c("fml", "rmf", "rml") ~ "menstrual_phase",
        coef == "bv" ~ "bv",
        coef %in% c("hsv", "hsv2") ~ "hsv", 
        TRUE ~ NA_character_)
    if (is.na(result)) stop(paste("invalid coef:", original_coef))
    
    result
}

get_probe_fold_changes <- function(coef) {
    file <- paste0("data/clean/", file_lookup(coef), "_probe_fold_changes.csv")
    
    read_csv(file,
        col_types = cols(
            ProbeId = col_integer(),
            Gene = col_character(),
            EntrezId = col_integer(),
            Log2FoldChange = col_double(),
            PValue = col_double(),
            FDR = col_double(),
            Coefficient = col_character()))
}

get_gene_fold_changes <- function(coef) {
    file <- paste0("data/clean/", file_lookup(coef), "_gene_fold_changes.csv")
    
    read_csv(file,
        col_types = cols(
            Coefficient = col_character(),
            Gene = col_character(),
            EntrezId = col_integer(),
            Log2FoldChange = col_double(),
            PValue = col_double(), 
            FDR = col_double()))
}

get_hallmark <- function(coef) {
    file <- paste0("data/clean/", file_lookup(coef), "_hallmark.csv")
    
    read_csv(file, 
        col_types = cols(
            NGenes = col_integer(),
            Direction = col_character(),
            PValue = col_double(),
            FDR = col_double(),
            GeneSet = col_character(),
            Coefficient = col_character()))    
}

get_other_studies_gene_fold_change <- function() {
    them <- .get_other_studies("gene_fold_changes", 
                cols(
                    Coefficient = col_character(),
                    Gene = col_character(),
                    EntrezId = col_character(),
                    Log2FoldChange = col_double(),
                    PValue = col_double(),
                    FDR = col_double()
                ))
    
    us <- get_gene_fold_changes("FmL") %>% 
        filter(Coefficient == "Follicular vs. luteal") %>% 
        mutate(EntrezId = as.character(EntrezId))
        
    add_us(them, us)    
}

get_other_studies_hallmark <- function() {
    them <- .get_other_studies("hallmark",  
                cols(
                    NGenes = col_integer(),
                    Direction = col_character(),
                    PValue = col_double(),
                    FDR = col_double(),
                    GeneSet = col_character(),
                    Coefficient = col_character()
                ))
    
    us <- get_hallmark("FmL") %>% 
        filter(Coefficient == "Follicular vs. luteal")
    
    add_us(them, us)    
}

get_other_studies_key <- function() {
    key <- read_csv("data/required_files/other_studies_key.csv", 
        col_types = cols(
            Study = col_character(),
            Paired = col_character(),
            Organ = col_character(),
            Sample = col_character(),
            Phasing = col_character(),
            Notes = col_character()
        )) %>% 
        factor_organ()    
    
    files <- list.files("data/clean/other_studies", full.names = TRUE) %>% 
        .[str_detect(., "pData")]
    
    read_pData <- function(f) {
        study <- str_remove(f, "_pData.csv")
        study <- str_remove(study, "data\\/clean\\/other_studies\\/")
        
        read.csv(f, stringsAsFactors = FALSE) %>% 
            mutate(Study = study)
    }
    
    ns <- map_df(files, read_pData) %>% count(Study, Phase) %>% spread(Phase, n)
    
    inner_join(key, ns, by = "Study") %>% 
        select(Study, Follicular, Luteal, everything())
}

factor_organ <- function(d) {
    d %>% 
        mutate(Organ = factor(Organ, 
            levels = c("Fallopian tube", "Endometrium", "Endocervix", "Ectocervix", "Vagina")))
}

.get_other_studies <- function(type, cols) {
    files <- list.files("data/clean/other_studies", full.names = TRUE)
    files <- files[str_detect(files, type)]    
    
    read_file <- function(file) {
        d <- read_csv(file, col_types = cols)
        d$Study <- str_extract(file, "[GSE|PMID][^_]*")
        d
    }
    
    key <- get_other_studies_key()
    
    map_df(files, read_file) %>% 
        inner_join(key, by = "Study")
}

add_us <- function(them, us) {
    us <- us %>% 
        mutate(Study = "Us", Organ = "Endocervix", Sample = "Cytobrush") %>% 
        factor_organ()
        
    bind_rows(them, us)    
}

get_ddpcr <- function() {
    read_csv("data/clean/ddpcr_granulysin.csv", col_types = cols(
            Plate = col_character(),
            Index = col_character(),
            StudyId = col_character(),
            Visit = col_factor(levels = c("Ring", "Luteal", "Follicular")),
            BV = col_character(),
            Hsv2 = col_character(),
            Drop = col_character(),
            Batch = col_character(),
            GNLY = col_double(),
            YWHAZ = col_double(),
            Log2GranulysinPerYWHAZ = col_double()))
}

get_ddpcr_vs_microarray <- function() {
   read_csv("data/clean/ddpcr_granulysin_vs_microarray.csv", col_types = cols(
       Plate = col_character(),
       Index = col_character(),
       StudyId = col_double(),
       Visit = col_character(),
       BV = col_character(),
       Hsv2 = col_character(),
       Drop = col_character(),
       GNLY = col_double(),
       YWHAZ = col_double(),
       Log2GranulysinPerYWHAZ = col_double(),
       ProbeId = col_double(),
       MicroarrayId = col_character(),
       Expression = col_double(),
       Date = col_date(format = ""),
       RIN = col_double(),
       Enrollment = col_date(format = ""),
       Status = col_character(),
       Progesterone = col_double(),
       EntrezId = col_double(),
       Gene = col_character())
   )
}

get_ddpcr_stats <- function() {
    read_csv("data/clean/ddpcr_granulysin_stats.csv", col_types = cols(
        Value = col_double(),
        Std.Error = col_double(),
        DF = col_double(),
        `t-value` = col_double(),
        `p-value` = col_double(),
        Coefficient = col_character()
    ))    
}

get_msd <- function() {
    read_csv("data/clean/msd_conc.csv", col_types = cols(
        Assay = col_character(),
        OfficialSymbol = col_character(),
        EntrezId = col_double(),
        Index = col_character(),
        StudyId = col_double(),
        Visit = col_character(),
        BV = col_character(),
        Hsv2 = col_character(),
        Drop = col_character(),
        Plate = col_character(),
        LOD = col_double(),
        Dilution = col_double(),
        Concentration = col_double(),
        Detectable = col_logical(),
        Log2Concentration = col_double(),
        Log2LOD = col_double()
    ))
}

get_elisa_granulysin <- function() {
    read_csv("data/clean/elisa_granulysin_conc.csv", col_types = cols(
        Assay = col_character(),
        EntrezId = col_double(),
        OfficialSymbol = col_character(),
        Plate = col_character(),
        StudyId = col_double(),
        Visit = col_character(),
        BV = col_character(),
        Hsv2 = col_character(),
        Drop = col_character(),
        Log2Concentration = col_double()
    ))
}

get_protein <- function() {
    bind_rows(
        get_elisa_granulysin(), 
        get_msd())
}

get_protein_stats <- function() {
    read_csv("data/clean/protein_stats.csv", col_types = cols(
            Assay = col_character(),
            OfficialSymbol = col_character(),
            EntrezId = col_double(),
            Value = col_double(),
            Std.Error = col_double(),
            DF = col_double(),
            `t-value` = col_double(),
            `p-value` = col_double(),
            Coefficient = col_character(),
            adjusted = col_double()
    ))
}
