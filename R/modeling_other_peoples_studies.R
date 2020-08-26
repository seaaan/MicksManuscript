# helpers ----------------------------------------------------------------------
get_x_from_probe_id <- function(x, data, rownames = TRUE) {
   lookup <- x
   # Get the probe identifiers that are mapped to a gene symbol
   mapped_probes <- mappedkeys(lookup)
   # Convert to a list
   lookup <- as.list(lookup[mapped_probes])

   if (rownames) {
      probes <- rownames(data)
   } else {
      probes <- data
   }
   
   # nulls get silently dropped, so need to replace them with NA
   unname(unlist(lapply(lookup[probes], 
      function(x) ifelse(is.null(x), NA, x))))
}

# works for microarrays 
get_fData_probe <- function(data, symbol, entrez) {
   fData <- data.frame(ProbeID = rownames(data), 
      Symbol = get_x_from_probe_id(symbol, data), 
      EntrezId = get_x_from_probe_id(entrez, data))
   rownames(fData) <- fData$ProbeID   
   fData
}

# common code to preprocess cel files from affymetrix microarrays
# gse = path to cel files
# sample_threshold = minimum number of samples in which probe must be expressed
preprocess_cel_files <- function(gse, sample_threshold, pData) {
    x <- ReadAffy(filenames = paste0(gse, list.celfiles(gse)), compress = TRUE)
    
    x_detection <- se.exprs(mas5calls(x))
    
    x2 <- rma(x)
    
    if (!all(rownames(x2) == rownames(x_detection))) {
        stop("some rows in the detection calls differ from in the expression data")
    }
    
    filtered <- x2[rowSums(x_detection < 0.05) >= sample_threshold, ]
    
    fData(filtered) <- get_fData_probe(filtered, hgu133plus2SYMBOL, hgu133plus2ENTREZID)
    fData(filtered)$ProbeID <- rownames(fData(filtered))
    fData(filtered)$TargetID <- fData(filtered)$Symbol
    fData(filtered)$ENTREZ_GENE_ID <- fData(filtered)$EntrezId
    
    # remove probes with NA entrez IDs
    filtered <- filtered[!is.na(fData(filtered)$EntrezId), ]
    
    if (!all(rownames(pData) == colnames(filtered))) {
        stop("some samples in the pData don't match samples in the expression data")
    }
    
    pData(filtered) <- pData
    
    filtered
}

# GSE80455 ---------------------------------------------------------------------
# endocervical tissue from surgical discards

# initial file
# need to get detection p-values out of here
# file structure is probe IDs, then sets of two columns per PTID: 
#   col 1 = expression
#   col 2 = detection p-value
# need to reformat for lumi, create detection matrix
GSE80455_raw <- read.table(GSE80455_raw_file, sep = "\t", header = TRUE)
GSE80455_detection <- as.matrix(GSE80455_raw[ , seq(from = 3, to = 25, by = 2)])
rownames(GSE80455_detection) <- GSE80455_raw$ID_REF
colnames(GSE80455_detection) <- colnames(GSE80455_raw)[seq(from = 2, to = 24, by = 2)]

# need to get normalized data from GEO because raw data doesn't include se.exprs
GSE80455_normalized <- getGEO("GSE80455", destdir = "data/raw/other_studies/GSE80455/")[[1]]

# put matrices in same order
GSE80455_detection <- GSE80455_detection[rownames(GSE80455_normalized), ]

# filtering
GSE80455_filtered <- GSE80455_normalized[rowSums(GSE80455_detection < 0.05) > 4, ]

# more convenient names
pData(GSE80455_filtered)$Phase <- pData(GSE80455_filtered)$`phases of the menstrual cycle:ch1`
fData(GSE80455_filtered)$ENTREZ_GENE_ID <- fData(GSE80455_filtered)$Entrez_Gene_ID
fData(GSE80455_filtered)$ProbeID <- fData(GSE80455_filtered)$Probe_Id
fData(GSE80455_filtered)$TargetID <- fData(GSE80455_filtered)$Symbol

#  fit
phase <- factor(pData(GSE80455_filtered)$Phase)
GSE80455_design <- model.matrix(~ 0 + phase) 
GSE80455_fit <- lmFit(GSE80455_filtered, GSE80455_design)
GSE80455_fit <- eBayes(contrasts.fit(GSE80455_fit, 
   contrasts = makeContrasts(FmL = phaseFollicular - phaseLuteal, 
      levels = GSE80455_design)))

create_outputs("GSE80455", GSE80455_filtered, GSE80455_fit, "data/clean/other_studies/GSE80455")

# GSE 29981 --------------------------------------------------------------------
# endometrial biopsies then laser captured microdissection to glandular epithelium
GSE29981_pData <- read.csv(GSE29981_metadata_file)
rownames(GSE29981_pData) <- paste0(GSE29981_pData$MicroarrayId, ".CEL.gz")
GSE29981_pData$Phase <- ifelse(GSE29981_pData$Progesterone > 1, "Luteal", "Follicular")

GSE29981_filtered <- preprocess_cel_files(GSE29981_cel_directory, 6, GSE29981_pData)

#  fit
phase <- factor(pData(GSE29981_filtered)$Phase)
GSE29981_design <- model.matrix(~ 0 + phase) 
GSE29981_fit <- lmFit(GSE29981_filtered, GSE29981_design)
GSE29981_fit <- eBayes(contrasts.fit(GSE29981_fit, 
    contrasts = makeContrasts(FmL = phaseFollicular - phaseLuteal, 
        levels = GSE29981_design)))

# A few probes conflict, they are genuine conflicts (not eg NA Entrez IDs)
create_outputs("GSE29981", GSE29981_filtered, GSE29981_fit, "data/clean/other_studies/GSE29981")

# GSE28044 ---------------------------------------------------------------------
# fallopian tube epithelium, laser microcapture, from BRCA+ and BRCA-
# surgical discards
# epithelium proximal to fimbria
GSE28044_pData <- read.csv(GSE28044_metadata_file)
rownames(GSE28044_pData) <- paste0(GSE28044_pData$MicroarrayId, ".cel.gz")

GSE28044_filtered <- preprocess_cel_files(GSE28044_cel_directory, 9, GSE28044_pData)

#  fit
phase <- factor(pData(GSE28044_filtered)$Phase)
brca <- factor(str_extract(pData(GSE28044_filtered)$Brca, "yes|no"))
GSE28044_design <- model.matrix(~ 0 + phase + brca) 
GSE28044_fit <- lmFit(GSE28044_filtered, GSE28044_design)
GSE28044_fit <- eBayes(contrasts.fit(GSE28044_fit, 
    contrasts = makeContrasts(FmL = phaseFollicular - phaseLuteal, 
        levels = GSE28044_design)))

create_outputs("GSE28044", GSE28044_filtered, GSE28044_fit, "data/clean/other_studies/GSE28044")

# GSE86491-------------------------------------------------------------------
# Paired suction curette endometrial tissue 
# 6-8 days post start of menses ("proliferative" = follicular)
# 7-9 days after positive ovulation test ("post-ovulation" = luteal)
GSE86491_raw <- read.table(GSE86491_raw_file, sep = "\t", header = TRUE)

GSE86491_counts <- GSE86491_raw[,2:15]
rownames(GSE86491_counts) <- GSE86491_raw[,1]

GSE86491_fData <- tibble::tibble(EnsemblId = rownames(GSE86491_counts), 
    EntrezId = get_x_from_probe_id(org.Hs.egENSEMBL2EG, GSE86491_counts), 
    ENTREZ_GENE_ID = EntrezId,
    TargetID = get_x_from_probe_id(org.Hs.egSYMBOL, EntrezId, FALSE))
GSE86491_fData <- as.data.frame(GSE86491_fData)
rownames(GSE86491_fData) <- GSE86491_fData$EnsemblId
GSE86491_fData$ProbeID <- rownames(GSE86491_fData)

GSE86491_data <- DGEList(counts = GSE86491_counts, genes = GSE86491_fData)
# remove NA genes
GSE86491_data <- GSE86491_data[!is.na(GSE86491_data$genes$EntrezId), ]

GSE86491_pData <- data.frame(
    # phase is from GEO website
    Phase = factor(rep(c("Follicular", "Luteal"), 7)),
    # donor is based on an email from Hanna Åmark, who stated that samples from the 
    # same donor were sequential (14 Jan 2020)
    Donor = factor(rep(LETTERS[1:7], each = 2))
)
rownames(GSE86491_pData) <- colnames(GSE86491_data)

phase <- GSE86491_pData$Phase
donor <- GSE86491_pData$Donor

GSE86491_design <- model.matrix(~ 0 + phase + donor)
keep <- edgeR::filterByExpr(GSE86491_data, GSE86491_design)
GSE86491_data <- GSE86491_data[keep, , keep.lib.sizes = FALSE]

GSE86491_data <- calcNormFactors(GSE86491_data)
GSE86491_data <- voom(GSE86491_data, GSE86491_design, plot = FALSE)

GSE86491_fit <- lmFit(GSE86491_data, GSE86491_design)
GSE86491_fit <- eBayes(contrasts.fit(GSE86491_fit, 
    contrasts = makeContrasts(FmL = phaseFollicular - phaseLuteal, 
        levels = GSE86491_design)), trend = TRUE)

# create an expression set so that the create_outputs function will work
GSE86491_fake_filtered <- ExpressionSet(assayData = GSE86491_data$E)
fData(GSE86491_fake_filtered) <- GSE86491_data$genes
pData(GSE86491_fake_filtered) <- GSE86491_pData

create_outputs("GSE86491", GSE86491_fake_filtered, GSE86491_fit, "data/clean/other_studies/GSE86491")

# PMID26750085 ------------------------------------------------------------
# paired vaginal biopsies taken from healthy premenopausal women in the
# follicular and luteal phases based on average menstrual period length
# and serum hormone levels
# these data are not available in a public repository. They were received
# on 15 Mar 2018 from Irina Zalenskaya <ZalensIA@EVMS.EDU>, to whom I had
# been directed after emailing first author Andrea Thurman <ThurmaAR@EVMS.EDU>
# I found out that there was data from an additional study, which I received
# from Irina on 5 Dec 2018
# V2 = FOL
# V3 = LUT
# Studies = CONRAD119 and CONRAD124

if (ANALYZE_PMID26750085) {
    conrad124 <- ReadAffy(filenames = paste0(PMID26750085_124_study_cel_directory, 
        list.celfiles(PMID26750085_124_study_cel_directory)), compress = FALSE)
    
    # detection p-values
    conrad124_detection <- se.exprs(mas5calls(conrad124))
    
    conrad124 <- rma(conrad124)
    
    conrad119 <- ReadAffy(filenames = paste0(PMID26750085_119_study_cel_directory, 
        list.celfiles(PMID26750085_119_study_cel_directory)), compress = FALSE)
    
    # detection p-values
    conrad119_detection <- se.exprs(mas5calls(conrad119))
    
    conrad119 <- rma(conrad119)
    
    both_studies <- BiocGenerics::combine(conrad124, conrad119)
    
    # remove not detected in at least 10 samples
    PMID26750085_filtered <- both_studies[(rowSums(conrad119_detection < 0.05) + rowSums(conrad124_detection < 0.05)) > 10, ]
    fData(PMID26750085_filtered) <- get_fData_probe(PMID26750085_filtered, hgu133plus2SYMBOL, hgu133plus2ENTREZID)
    fData(PMID26750085_filtered)$TargetID <- fData(PMID26750085_filtered)$Symbol
    fData(PMID26750085_filtered)$ENTREZ_GENE_ID <- fData(PMID26750085_filtered)$EntrezId
    
    # remove rows with NA annotations
    PMID26750085_filtered <- PMID26750085_filtered[!is.na(fData(PMID26750085_filtered)$Symbol), ]
    
    conrad124_pData <- readxl::read_excel(PMID26750085_124_study_metadata_file)
    colnames(conrad124_pData) <- c("MicroarrayId", "Sample", "Key1", "Phase")
    conrad124_pData$Ptid <- paste0("124_", stringr::str_replace(conrad124_pData$Sample, ",.*", ""))
    conrad124_pData$Phase <- ifelse(conrad124_pData$Phase == "FOL", "Follicular", "Luteal")
    
    conrad119_pData <- readxl::read_excel(PMID26750085_119_study_metadata_file)
    colnames(conrad119_pData) <- c("MicroarrayId", "Sample", "Key1")
    conrad119_pData$Ptid <- paste0("119_", stringr::str_replace(conrad119_pData$Sample, ",.*", "")) 
    conrad119_pData$Phase <- ifelse(conrad119_pData$Key1 == "V2", "Follicular", "Luteal")
    
    PMID26750085_pData <- as.data.frame(bind_rows(conrad124_pData, conrad119_pData))
    PMID26750085_pData$Study <- str_sub(PMID26750085_pData$Ptid, 1, 3)
    rownames(PMID26750085_pData) <- paste0(PMID26750085_pData$MicroarrayId, ".CEL")
    
    # put pData in order of filtered
    PMID26750085_pData <- PMID26750085_pData[colnames(PMID26750085_filtered), ]
    pData(PMID26750085_filtered) <- PMID26750085_pData
    
    phase <- factor(pData(PMID26750085_filtered)$Phase)
    ptid <- factor(pData(PMID26750085_filtered)$Ptid)
    
    PMID26750085_design <- model.matrix(~ 0 + phase + ptid)
    
    PMID26750085_fit <- eBayes(contrasts.fit(lmFit(PMID26750085_filtered, PMID26750085_design), 
        contrasts = makeContrasts(FmL = phaseFollicular - phaseLuteal, 
            levels = PMID26750085_design)))
    
    create_outputs("PMID26750085", PMID26750085_filtered, PMID26750085_fit, "data/clean/other_studies/PMID26750085")
}

# GSE6364 ----------------------------------------------------------------------
# endometrial biopsies
# D = endometriosis; N = healthy control
# PE = proliferative; ESE = early secretory; MSE = mid secretory
GSE6364_pData <- read.csv(GSE6364_metadata_file)
rownames(GSE6364_pData) <- paste0(GSE6364_pData$MicroarrayId, ".CEL.gz")

GSE6364_filtered <- preprocess_cel_files(GSE6364_cel_directory, 9, GSE6364_pData)

# only want to compare proliferative vs mid secretory in healthy controls
# PE = proliferative
# MSE = mid secretory
# ESE = early secretory
# N = healthy control, D = disease
GSE6364_filtered <- GSE6364_filtered[ , !str_detect(pData(GSE6364_filtered)$Code, "ESE") & str_detect(pData(GSE6364_filtered)$Code, "N")]
pData(GSE6364_filtered)$Phase <- ifelse(str_detect(pData(GSE6364_filtered)$Code, "PE"), 
    "Follicular", "Luteal")

phase <- factor(pData(GSE6364_filtered)$Phase)
GSE6364_design <- model.matrix(~ 0 + phase)

GSE6364_fit <- eBayes(contrasts.fit(lmFit(GSE6364_filtered, GSE6364_design), 
    contrasts = makeContrasts(FmL = phaseFollicular - phaseLuteal, 
        levels = GSE6364_design)))

# a few genes come up in opposite directions, they're filtered out
create_outputs("GSE6364", GSE6364_filtered, GSE6364_fit, "data/clean/other_studies/GSE6364")

# GSE122248-------------------------------------------------------------------
# endo and ecto tissue from hysterectomies

GSE122248_raw_data <- read.table(GSE122248_raw_file, sep = "\t", header = TRUE)

GSE122248_counts <- GSE122248_raw_data[,7:31] %>% as.matrix(header = TRUE)
rownames(GSE122248_counts) <- GSE122248_raw_data$Geneid

GSE122248_fData <- tibble::tibble(EnsemblId = rownames(GSE122248_counts), 
    EntrezId = get_x_from_probe_id(org.Hs.egENSEMBL2EG, GSE122248_counts), 
    ENTREZ_GENE_ID = EntrezId,
    TargetID = get_x_from_probe_id(org.Hs.egSYMBOL, EntrezId, FALSE))
GSE122248_fData <- as.data.frame(GSE122248_fData)
rownames(GSE122248_fData) <- GSE122248_fData$EnsemblId
GSE122248_fData$ProbeID <- rownames(GSE122248_fData)

GSE122248_pData <- read_csv(GSE122248_pData_file, col_types = "cccc") %>% 
    rename(Donor = Ptid) %>% 
    mutate(Phase = factor(Phase, levels = c("Luteal", "Follicular"))) %>% 
    as.data.frame()

rownames(GSE122248_pData) <- GSE122248_pData$SampleName

GSE122248_data <- DGEList(counts = GSE122248_counts, genes = GSE122248_fData, samples = GSE122248_pData)
# remove NA genes
GSE122248_data <- GSE122248_data[!is.na(GSE122248_data$genes$EntrezId), ]

normalize_filter_model_and_export_GSE122248 <- function(tissue) {
    data <- GSE122248_data[ , GSE122248_data$samples$Tissue == tissue]
    
    pData <- data$samples
    
    phase <- data$samples$Phase
    
    design <- model.matrix(~ 0 + phase)
    keep <- edgeR::filterByExpr(data, design)
    data <- data[keep, , keep.lib.sizes = FALSE]
    
    data <- calcNormFactors(data)
    data <- voom(data, design, plot = FALSE)
    
    fit <- lmFit(data, design)
    fit <- eBayes(contrasts.fit(fit, 
        contrasts = makeContrasts(FmL = phaseFollicular - phaseLuteal, 
            levels = design)), trend = TRUE)
    
    # create an expression set so that the create_outputs function will work
    fake_filtered <- ExpressionSet(assayData = data$E)
    fData(fake_filtered) <- data$genes
    pData(fake_filtered) <- pData
    
    name <- paste0("GSE122248-", tissue)    
    
    create_outputs(name, fake_filtered, fit, paste0("data/clean/other_studies/", name))
}
normalize_filter_model_and_export_GSE122248("ENDO")
normalize_filter_model_and_export_GSE122248("ECTO")
