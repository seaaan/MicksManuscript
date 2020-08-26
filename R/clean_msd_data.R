# Data import and analysis ----------------------------------------------------
# there is data from two sets of plates, one 5-plex and one 9-plex
read_plex <- function(n) {
    quiet_read_csv(paste0(msd_data_folder, "MSD_", n, "Plex_raw.csv"), skip = 1) %>% 
        mutate(Plex = n)
}

d <- map_df(c(5, 9), read_plex) %>% 
    mutate(Assay = str_replace(Assay, "TNF.*", "TNF-a"),
        Assay = str_replace(Assay, "MIP-3.*", "MIP-3a"),
        Assay = str_replace(Assay, "\U03B2", "b"), 
        Assay = str_replace(Assay, "IFN.*", "IFN-y"), 
        # this one does IL-1α
        Assay = str_replace(Assay, "IL-1[^b]$", "IL-1a"))

# add entrez IDs and official symbols
d <- left_join(d, read_csv("data/required_files/cytokine_lookup.csv", col_types = cols(
        Gene = col_character(),
        EntrezId = col_double(),
        OfficialSymbol = col_character()
    )), 
    by = c("Assay" = "Gene"))

# none of the samples start with the letter S
# and some unused wells were labeled as "X"
samples <- d %>% 
    filter(!str_detect(Sample, "^S")) %>% 
    filter(Sample != "X") 

# add metadata
key <- get_metadata("ddpcr_proteomics")

samples <- left_join(samples, key, by = "Sample")

# adjust for dilution
# 5 plex were 1:300
# 9 plex were 1:30
samples <- samples %>% 
    mutate(Dilution = ifelse(str_detect(`Plate Name`, "5Plex"), 300, 30)) 

valid_samples <- samples %>% 
    filter(Drop == "Don't drop") %>% 
    mutate(StudyId = factor(StudyId), Visit = factor(Visit, levels = c("Ring", "Luteal", "Follicular"))) %>% 
    rename(Plate = `Plate Name`) %>% 
    arrange(BV)

valid_samples_conc <- valid_samples %>% 
    group_by(Assay, Plate) %>% 
    # limit of detection is the lower of the calculated detection limit or the lowest observed
    mutate(LOD = min(c(`Calc. Concentration`, `Detection Limits: Calc. Low`), na.rm = TRUE)) %>% 
    group_by(Assay, Index, StudyId, Visit, BV, Hsv2, Drop, Plate, LOD, Dilution, OfficialSymbol, EntrezId) %>% 
    summarise(Concentration = mean(`Calc. Concentration`, na.rm = TRUE)) %>% 
    ungroup() %>% 
    mutate(Detectable = !is.na(Concentration)) %>% 
    # adjust for dilution
    mutate(Concentration = Concentration * Dilution, LOD = LOD * Dilution) %>% 
    # assign missing values to 1/2 the LOD
    mutate(Concentration = ifelse(Detectable, Concentration, LOD / 2)) %>% 
    mutate(Log2Concentration = log2(Concentration), Log2LOD = log2(LOD)) %>% 
    select(Assay, OfficialSymbol, EntrezId, everything())

write_csv(valid_samples_conc, "data/clean/msd_conc.csv")
