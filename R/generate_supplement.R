# MICROARRAY ---------------------------------------------------------------------------
us_gene_fold_changes <- bind_rows(
        get_gene_fold_changes("FmL"),
        get_gene_fold_changes("BV"),
        get_gene_fold_changes("HSV")) %>% 
    group_by(Coefficient) %>% 
    nest()

us_gene_fold_changes_list <- us_gene_fold_changes$data
names(us_gene_fold_changes_list) <- us_gene_fold_changes$Coefficient
writexl::write_xlsx(us_gene_fold_changes_list, "supplement/Microarray_fold_changes.xlsx")

us_hallmark <- bind_rows(
    get_hallmark("FmL"),
    get_hallmark("BV"),
    get_hallmark("HSV")) %>% 
    group_by(Coefficient) %>% 
    nest()

us_hallmark_list <- us_hallmark$data
names(us_hallmark_list) <- us_hallmark$Coefficient
writexl::write_xlsx(us_hallmark_list, "supplement/Microarray_gene_set_testing.xlsx")

us_ora <- bind_rows(read_csv("data/clean/hsv_ora_GO.csv"), 
    read_csv("data/clean/bv_ora_GO.csv")) %>% 
    group_by(Direction) %>% 
    nest()

us_ora_list <- us_ora$data
names(us_ora_list) <- us_ora$Direction 
writexl::write_xlsx(us_ora_list, "supplement/Microarray_overrepresentation_analysis.xlsx")

# PROTEIN AND DDPCR ---------------------------------------------------------------------
us_protein <- get_protein_stats() %>% 
    rename(Log2FoldChange = Value) %>% 
    group_by(Coefficient) %>% 
    nest()

us_protein_list <- us_protein$data
names(us_protein_list) <- us_protein$Coefficient
writexl::write_xlsx(us_protein_list, "supplement/MSD-ELISA_fold_changes.xlsx")

get_ddpcr_stats() %>% 
    rename(Log2FoldChange = Value) %>% 
    writexl::write_xlsx("supplement/ddPCR_fold_changes.xlsx")

# OTHER STUDIES -----------------------------------------------------------------------------------
other_studies_fold_changes <- get_other_studies_gene_fold_change() %>% 
    mutate(Name = paste(Organ, str_remove(Study, "-.*"), sep = "_")) %>% 
    group_by(Name) %>% 
    nest()

other_studies_fold_changes_list <- other_studies_fold_changes$data
names(other_studies_fold_changes_list) <- other_studies_fold_changes$Name
writexl::write_xlsx(other_studies_fold_changes_list, "supplement/Other_studies_fold_changes.xlsx")

other_studies_hallmark <- get_other_studies_hallmark() %>% 
    mutate(Name = paste(Organ, str_remove(Study, "-.*"), sep = "_")) %>% 
    group_by(Name) %>% 
    nest()

other_studies_hallmark_list <- other_studies_hallmark$data
names(other_studies_hallmark_list) <- other_studies_hallmark$Name
writexl::write_xlsx(other_studies_hallmark_list, "supplement/Other_studies_hallmark.xlsx")