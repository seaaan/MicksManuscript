# subset to useful columns, improve column names, add identifying coef
format_tt <- function(tt, coef) {
    tt %>% 
        dplyr::transmute(ProbeId = ProbeID, Gene = TargetID, EntrezId = ENTREZ_GENE_ID,
            Log2FoldChange = logFC, PValue = P.Value, FDR = adj.P.Val) %>% 
        canonicalize_gene_symbols() %>% 
        mutate(Coefficient = coef) %>% 
        arrange(FDR)
} 

# convert probe fold changes to gene fold changes
collapse_probes <- function(tt) {
    # count differentially expressed genes
    # group into unique study arms (Coefficient)
    # and into unique genes (Gene, EntrezId)
    # and into whether that probe is significant or not
    # drop if multiple significant in opposite directions
    # otherwise take probe with lowest p-value per gene
    d <- tt %>%
        group_by(Coefficient, Gene, EntrezId, Sig = FDR < 0.05) %>%
        # remove cases where multiple probes are significant in different directions
        mutate(Direction = ifelse(all(Log2FoldChange < 0) | all(Log2FoldChange > 0) | !Sig, "okay", "problem")) 
        
        if (any(d$Direction == "problem")) {
            warning(paste(sum(d$Direction == "problem"), "probes are significant in opposite directions for the same gene"))
        }
    
    d %>%
        dplyr::filter(Direction != "problem") %>%
        # collapse multiple probes into a single gene by taking lowest p-value
        group_by(Coefficient, Gene, EntrezId) %>%
        dplyr::arrange(FDR, PValue) %>%
        dplyr::slice(1) %>%
        ungroup() %>% 
        dplyr::arrange(Coefficient, PValue, FDR) %>% 
        select(Coefficient, Gene, EntrezId, Log2FoldChange, PValue, FDR)
}

# save row names as column, add identifying coef
format_camera <- function(camera, coef) {
    camera %>% 
        mutate(GeneSet = rownames(.), Coefficient = coef)
}

# update gene names
# requires d to have a column containing entrez IDs
# adds a column named "Gene" containing the corresponding gene symbols
canonicalize_gene_symbols <- function(d, entrez_id_column = "EntrezId", 
    current_symbols_column = "Gene") {
    
    # it breaks below if the input has a column called "SYMBOL"
    stopifnot(all(colnames(d) != "SYMBOL"))
    
    # replace Illumina-supplied gene symbols with up-to-date gene symbols
    lookup <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
        keys = as.character(unique(d[[entrez_id_column]])),
        columns = c("ENTREZID", "SYMBOL"), keytype = "ENTREZID"))
    
    combined <- merge(d, lookup, by.x = entrez_id_column, by.y = "ENTREZID", all.x = TRUE) 
    
    # some probes don't have an EntrezId associated because they measure 
    # something other than a gene, such as an miRNA or an RNA of another type
    # in these cases, they do have a symbol, so get those symbols to use
    # in cases where there is no entrez ID and so AnnotationDbi::select returns
    # NA for symbol
    gene_names <- ifelse(is.na(combined[[entrez_id_column]]), combined[[current_symbols_column]], combined$SYMBOL)
    
    combined %>%
        mutate(Gene = gene_names) %>%
        dplyr::select(-SYMBOL)
}

# x = vector to replace
# d = data frame containing columns for lookup
# desired = column name in d of what you want it to say
# current = column name in d of what it currently says
lookup <- function(x, d, desired, current) {
    lookup <- d[[desired]]
    names(lookup) <- d[[current]]
    
    unname(lookup[x])
}


clean_hallmark_gene_sets <- function(d, n_lump = 3) {
    hallmarks <- get_hallmark_gene_set_names_and_categories()
    
    d$Category <- lookup(d$GeneSet, hallmarks, "Category", "Name")
    d$GeneSet <- lookup(d$GeneSet, hallmarks, "ShortName", "Name")
    
    d %>% 
        ungroup() %>% # ungroup avoids warning about unequal Category levels
        mutate(Category = fct_lump(Category, n_lump)) %>% 
        mutate(Category = str_to_title(Category))
}

create_outputs <- function(gse, filtered, fit, file_prefix, coef = "FmL") {
    write_csv(pData(filtered), path = paste0(file_prefix, "_pData.csv"))
     
    tt <- format_tt(topTable(fit, coef = "FmL", n = Inf), "Follicular vs. luteal") %>% 
        arrange(Coefficient, FDR, PValue)
    
    write_csv(tt, path = paste0(file_prefix, "_probe_fold_changes.csv"))
    
    gene_fold_changes <- collapse_probes(tt)
    
    write_csv(gene_fold_changes, path = paste0(file_prefix, "_gene_fold_changes.csv"))
    
    # camera
    identifiers <- as.character(fData(filtered)[,"ENTREZ_GENE_ID"])
    index <- ids2indices(gene.sets = get_hallmark_gene_sets(), identifiers = identifiers)
    hallmark <- format_camera(cameraPR(fit$t[ , coef], index), "Follicular minus luteal")
    write_csv(hallmark, path = paste0(file_prefix, "_hallmark.csv"))
    
    invisible(NULL)
}


# runs mixed models on the data, assumes measurement variable is called Readout (e.g. concentration)
# Needs to have columns Visit, BV, Hsv2, StudyId in addition
# Can be used for ddpcr, msd, elisa data
# fits two separate models based on advice from stats consultation
# tests RvF, RvL, FvL, BV, HSV
lme_function <- function(d, model = gaussian_model) {
    d <- d %>% 
        ungroup() %>% 
        mutate(Visit = factor(Visit, levels = c("Luteal", "Follicular", "Ring")))
    
    first <- model(d) %>% 
        mutate(Coef = recode(Coef, 
            "VisitFollicular" = "Follicular vs. luteal", 
            "VisitRing" = "Ring vs. luteal"))
    
    d <- d %>% mutate(Visit = factor(Visit, levels = c("Follicular", "Ring", "Luteal")))
    
    second <- model(d) %>% 
        filter(Coef == "VisitRing") %>% 
        mutate(Coef = recode(Coef, "VisitRing" = "Ring vs. follicular")) 
    
    bind_rows(first, second) %>% 
        rename(Coefficient = Coef) %>% 
        mutate(Coefficient = recode(Coefficient,
            "Hsv2Positive" = "HSV-2 infection", 
            "BVY" = "Bacterial vaginosis"
        ))
}

gaussian_model <- function(d) {
    m <- lmerTest::lmer(Readout ~ Visit + BV + Hsv2 + (1|StudyId), data = d) %>% summary() %>% .$coefficients %>% 
        as.data.frame()  %>% mutate(Coef = rownames(.)) %>% 
        filter(!str_detect(Coef, "Intercept"))
    
    colnames(m) <- c("Value", "Std.Error", "DF", "t-value", "p-value", "Coef")
    m
}

logistic_model <- function(d) {
    m <- lme4::glmer(Detectable ~ Visit + BV + Hsv2 + (1|StudyId), data = d, family = binomial) %>% 
        summary() %>% .$coefficients %>% 
        as.data.frame()  %>% mutate(Coef = rownames(.)) %>% 
        filter(!str_detect(Coef, "Intercept")) 
    
    colnames(m) <- c("Value", "Std.Error", "z-value", "p-value", "Coef")
    m
}

# verbose output from read_csv sometimes needs to be suppressed
quiet_read_csv <- function(...) {
    suppressMessages(suppressWarnings(read_csv(...)))
}
