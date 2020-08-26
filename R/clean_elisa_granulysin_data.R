d <- plater::read_plates(paste0(elisa_data_folder, "ELISA_Granulysin_Plate", 1:2, ".csv"), c("Plate1", "Plate2"))

d <- d %>% 
    # subtract well background at 540 nm
    mutate(Signal = `450nm` - `540nm`) %>% 
    filter(Sample != "Empty") %>% 
    group_by(Plate, Sample) %>% 
    summarise(Signal = mean(Signal))

# subtract signal background from 0 samples
zeros <- d %>% 
    filter(Sample == "S0") %>% 
    select(Plate, Background = Signal)

d <- d %>% 
    inner_join(zeros, by = "Plate") %>% 
    mutate(BgSubtractedSignal = Signal - Background)

# calculate concentrations using standard curves
calculate_concentrations <- function(d) {
    # calculate standard curves
    standards <- d %>% 
        # Samples start with "S"
        filter(str_detect(Sample, "^S"), Sample != "S0") %>% 
        mutate(Concentration = as.numeric(str_remove(Sample, "S")))
    
    standards_model <- drm(BgSubtractedSignal ~ Concentration, fct = LL.4(names=c("Slope", "Lower", "Upper", "ED50")), data = standards)
    d$Concentration <- ED(standards_model, d$BgSubtractedSignal, type = "absolute", display = FALSE)[,"Estimate"]
    d    
}

# warnings okay, many 0s
concentrations <- d %>% 
    group_by(Plate) %>% 
    nest() %>% 
    # lots of 0s produce NaNs with warning, suppress those warnings (will later be replaced)
    mutate(result = list(map_df(data, ~suppressWarnings(calculate_concentrations(.x))))) %>% 
    unnest(result) %>% 
    select(-data) %>% 
    # remove standards and controls
    filter(!str_detect(Sample, "^S"), Sample != "Control")

concentrations <- concentrations %>% 
    mutate(Detectable = !is.na(Concentration)) %>% 
    ungroup() %>% 
    # set below LOD to LOD/3 
    mutate(Concentration = ifelse(!Detectable, min(Concentration, na.rm = TRUE) / 3, Concentration)) %>% 
    # adjust for sample dilution
    mutate(Concentration = 50 * Concentration)

key <- get_metadata("ddpcr_proteomics")

concentrations <- left_join(concentrations, key, by = "Sample")

# drop bad samples
valid_samples <- concentrations %>% 
    filter(Drop == "Don't drop") 

# clean up data frame
valid_samples_conc <- valid_samples %>% 
    mutate(Log2Concentration = log2(Concentration)) %>%
    mutate(Assay = "GNLY", EntrezId = 10578, OfficialSymbol = "GNLY") %>% 
    select(Assay, EntrezId, OfficialSymbol, Plate, StudyId, Visit, BV, Hsv2, Drop, Log2Concentration, Detectable) %>% 
    mutate(StudyId = factor(StudyId), Visit = factor(Visit, levels = c("Ring", "Luteal", "Follicular"))) %>% 
    arrange(BV)
    
write_csv(valid_samples_conc, "data/clean/elisa_granulysin_conc.csv")