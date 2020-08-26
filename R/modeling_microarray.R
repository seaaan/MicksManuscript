# menstrual phase and ring comparison ------------------------------------------ 
to_model <- get_filtered_expression_set()

phase <- factor(pData(to_model)$Visit)
bv_status <- factor(pData(to_model)$BV)
hsv <- factor(pData(to_model)$Hsv)
batch <- factor(pData(to_model)$Batch)

design <- model.matrix(~ 0 + phase + bv_status + hsv + batch)

dup_cors <- duplicateCorrelation(to_model, design = design, 
    block = pData(to_model)$StudyId)

fit <- lmFit(to_model, design, block = pData(to_model)$StudyId, 
    correlation = dup_cors$consensus.correlation)

contrasts <- makeContrasts(FmL = phaseFollicular - phaseLuteal, 
    RmL = phaseRing - phaseLuteal,
    RmF = phaseRing - phaseFollicular,
    levels = design)

fit <- eBayes(contrasts.fit(fit, contrasts = contrasts))

menstrual_tt <- bind_rows(
        format_tt(topTable(fit, coef = "FmL", n = Inf), "Follicular vs. luteal"),
        format_tt(topTable(fit, coef = "RmL", n = Inf), "Ring vs. luteal"),
        format_tt(topTable(fit, coef = "RmF", n = Inf), "Ring vs. follicular")) %>% 
    arrange(Coefficient, FDR, PValue)

write_csv(menstrual_tt, path = "data/clean/menstrual_phase_probe_fold_changes.csv")

gene_fold_changes <- collapse_probes(menstrual_tt)

write_csv(gene_fold_changes, path = "data/clean/menstrual_phase_gene_fold_changes.csv")

identifiers <- as.character(fData(to_model)[,"ENTREZ_GENE_ID"])
index <- ids2indices(gene.sets = get_hallmark_gene_sets(), identifiers = identifiers)
menstrual_hallmark <- bind_rows(
    format_camera(cameraPR(fit$t[,"FmL"], index), "Follicular vs. luteal"),
    format_camera(cameraPR(fit$t[,"RmL"], index), "Ring vs. luteal"),
    format_camera(cameraPR(fit$t[,"RmF"], index), "Ring vs. follicular"))
    
write_csv(menstrual_hallmark, path = "data/clean/menstrual_phase_hallmark.csv")

## BV and HSV2-----------------------------------------------------------------
# HSV -------------------------------------------------------------------------
design <- model.matrix(~ 0 + hsv + bv_status + phase + batch)

dup_cors <- duplicateCorrelation(to_model, design = design, 
    block = pData(to_model)$StudyId)

fit <- lmFit(to_model, design, block = pData(to_model)$StudyId, 
    correlation = dup_cors$consensus.correlation)

contrasts <- makeContrasts(HSV2 = hsvPositive - hsvNegative,
    levels = design)

fit <- eBayes(contrasts.fit(fit, contrasts = contrasts))

hsv_tt <- format_tt(topTable(fit, n = Inf, coef = "HSV2"), "HSV-2 infection") %>% 
    arrange(FDR, PValue)

write_csv(hsv_tt, path = "data/clean/hsv_probe_fold_changes.csv")

hsv_gene_fold_changes <- collapse_probes(hsv_tt)
write_csv(hsv_gene_fold_changes, 
    path = "data/clean/hsv_gene_fold_changes.csv")

hsv_hallmark <- format_camera(cameraPR(fit$t[,"HSV2"], index), "HSV2")
write_csv(hsv_hallmark, path = "data/clean/hsv_hallmark.csv")

# BV --------------------------------------------------------------------------
design <- model.matrix(~ 0 + bv_status + hsv + phase + batch)

dup_cors <- duplicateCorrelation(to_model, design = design, 
    block = pData(to_model)$StudyId)

fit <- lmFit(to_model, design, block = pData(to_model)$StudyId, 
    correlation = dup_cors$consensus.correlation)

contrasts <- makeContrasts(BV = bv_statusY - bv_statusN,
    levels = design)

fit <- eBayes(contrasts.fit(fit, contrasts = contrasts))

bv_tt <- format_tt(topTable(fit, n = Inf, coef = "BV"), "Bacterial vaginosis") %>% 
    arrange(FDR, PValue)

write_csv(bv_tt, path = "data/clean/bv_probe_fold_changes.csv")

bv_gene_fold_changes <- collapse_probes(bv_tt)
write_csv(bv_gene_fold_changes, 
    path = "data/clean/bv_gene_fold_changes.csv")

bv_hallmark <- format_camera(cameraPR(fit$t[,"BV"], index), "Bacterial vaginosis")
write_csv(bv_hallmark, path = "data/clean/bv_hallmark.csv")

## Over representation analysis of BV and HSV-2 ------------------------------

# function to perform overrepresentation analysis

# d is vector of entrez IDs of the DEGs of interest
# f is either limma::goana or limma::kegga, which do ORA on GO and KEGG databases
# passing in f is just a cute way to cut down on code duplication
ora <- function(d, f, universe) {
    f(d, universe = universe) %>% 
        mutate(Id = rownames(.), FDR = p.adjust(P.DE, method = "fdr")) %>% 
        arrange(P.DE)
}

# all genes detected in BV data set
bv_universe <- bv_gene_fold_changes %>% 
    filter(!is.na(EntrezId))

# entrez IDs of up DEGs from BV
bv_up <- bv_universe %>% 
    filter(FDR < 0.05, Log2FoldChange > 0) %>% 
    pull(EntrezId) %>% unique()

# entrez IDs of down DEGs from BV
bv_down <- bv_universe %>% 
    filter(FDR < 0.05, Log2FoldChange < 0) %>% 
    pull(EntrezId) %>% unique()

# entrez IDs of universe
bv_universe <- bv_universe %>% pull(EntrezId) %>% unique()

bv_up_results <- ora(bv_up, goana, bv_universe)
bv_down_results <- ora(bv_down, goana, bv_universe)

bv_ora_results <- bind_rows(
    mutate(bv_up_results, Direction = "Up during BV", Which = "GO"),
    mutate(bv_down_results, Direction = "Down during BV", Which = "GO")
)

write_csv(bv_ora_results, "data/clean/bv_ora_GO.csv")

# HSV-2 ------------------------------------------------------------------

# all genes detected in hsv data set
hsv_universe <- hsv_tt %>% 
    filter(!is.na(EntrezId))

# entrez IDs of up DEGs from hsv
hsv_up <- hsv_universe %>% 
    filter(FDR < 0.05, Log2FoldChange > 0) %>% 
    pull(EntrezId) %>% unique()

# entrez IDs of down DEGs from hsv
hsv_down <- hsv_universe %>% 
    filter(FDR < 0.05, Log2FoldChange < 0) %>% 
    pull(EntrezId) %>% unique()

# entrez IDs of universe
hsv_universe <- hsv_universe %>% pull(EntrezId) %>% unique()

hsv_up_results <- ora(hsv_up, goana, hsv_universe)
hsv_down_results <- ora(hsv_down, goana, hsv_universe)

hsv_ora_results <- bind_rows(
    mutate(hsv_up_results, Direction = "Up during HSV-2", Which = "GO"),
    mutate(hsv_down_results, Direction = "Down during HSV-2", Which = "GO")
)

write_csv(hsv_ora_results, "data/clean/hsv_ora_GO.csv")
