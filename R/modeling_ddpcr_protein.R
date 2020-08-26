# ddpcr

# run mixed model function
# mixed model function is in "analysis_helpers.R" 
# tests RvF, RvL, FvL, BV, HSV
stats <- get_ddpcr() %>% 
    mutate(Readout = Log2GranulysinPerYWHAZ) %>% 
    lme_function()

write_csv(stats, "data/clean/ddpcr_granulysin_stats.csv")

# protein

# run mixed model function on each cytokine separately
# mixed model function is in "analysis_helpers.R" 
# tests RvF, RvL, FvL, BV, HSV
protein_stats <- get_protein() %>% 
    rename(Readout = Log2Concentration) %>% 
    group_by(Assay, EntrezId, OfficialSymbol) %>%
    nest() %>% 
    # in a couple of cases, StudyId has no effect and a singular
    # warning is emitted, okay to suppress
    mutate(stats = map(data, ~suppressMessages(lme_function(.x)))) %>% 
    select(-data) %>% 
    unnest(cols = c(stats)) %>% 
    arrange(`p-value`) %>% 
    # adjust separately for each coefficient
    group_by(Coefficient) %>% 
    mutate(adjusted = p.adjust(`p-value`, "fdr")) %>% 
    arrange(`p-value`) 

write_csv(protein_stats, "data/clean/protein_stats.csv")
