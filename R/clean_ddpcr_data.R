p1 <- quiet_read_csv(paste0(ddpcr_data_folder, "ddpcr_plate_1.csv")) %>% 
    add_plate(paste0(ddpcr_data_folder, "ddpcr_plate_1_layout.csv"), "Well") %>% 
    rename(Index = ncol(.)) %>% 
    mutate(Plate = "Plate1")

p2 <- quiet_read_csv(paste0(ddpcr_data_folder, "ddpcr_plate_2.csv")) %>% 
    plater::add_plate(paste0(ddpcr_data_folder, "ddpcr_plate_2_layout.csv"), "Well") %>% 
    rename(Index = ncol(.)) %>% 
    mutate(Plate = "Plate2")

redos <- quiet_read_csv(paste0(ddpcr_data_folder, "ddpcr_plate_3.csv")) %>% 
    plater::add_plate(paste0(ddpcr_data_folder, "ddpcr_plate_3_layout.csv"), "Well") %>% 
    mutate(Plate = "Plate3") %>% 
    mutate(Target = ifelse(str_detect(TargetType, "Ch1"), "GNLY", "YWHAZ"))

d <- bind_rows(p1, p2, redos) %>% 
    # remove NA wells from redo plate
    filter(!is.na(Index)) %>% 
    mutate(Batch = ifelse(Plate == "Plate3", "Redo", "Original"))

controls <- filter(d, Index %in% c("NRT", "NTC"))

d <- filter(d, !Index %in% c("NRT", "NTC"))

meta <- get_metadata("ddpcr_proteomics")

d <- left_join(d, meta, by = "Index")

# remove wells with fewer than 10000 droplets
d <- filter(d, AcceptedDroplets > 1E4)

# remove donors with samples that fell below limit ofdetection
d <- d %>% filter(!StudyId %in% c(6036, 6046))

# batch 2 only for redos
redos <- d %>% filter(Batch == "Redo") %>% .$Index %>% unique()
d <- d %>% 
    filter(Batch == "Redo" | (!Index %in% redos))

averaged <- d %>% 
    select(Plate, Well, Target, Concentration, AcceptedDroplets, Index:Drop) %>% 
    filter(!str_detect(Drop, "Drop for all")) %>% 
    spread(Target, Concentration) %>% 
    # multiple ddpcr replicates per visit
    group_by(Plate, Index, StudyId, Visit, BV, Hsv2, Drop, Batch) %>% 
    summarise(GNLY = mean(GNLY), YWHAZ = mean(YWHAZ), 
        Log2GranulysinPerYWHAZ = log2(GNLY / YWHAZ)) 

write_csv(averaged, "data/clean/ddpcr_granulysin.csv")

microarray <- get_participant_level_probe_expression() %>% 
    filter(Gene == "GNLY")

ddpcr_vs_microarray <- inner_join(averaged, microarray, by = c("StudyId", "Index", "Hsv2", "BV", "Visit", "Drop")) %>% 
    select(-contains("Batch"))

write_csv(ddpcr_vs_microarray, "data/clean/ddpcr_granulysin_vs_microarray.csv")
