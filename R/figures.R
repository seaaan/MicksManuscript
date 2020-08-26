###################################################################################################
# Figure 1 ------------------------------------------------------------------------
###################################################################################################

# volcano plots -------------------------------------------------------------
to_plot <- get_gene_fold_changes("FmL") %>% 
    # add blank line after coefficient
    mutate(Coefficient = paste0(Coefficient, "\n")) %>% 
    # desired order: RvL, RvF, FvL
    mutate(Coefficient = factor(Coefficient, levels = 
            c("Ring vs. luteal\n", "Ring vs. follicular\n", "Follicular vs. luteal\n")))

volcano_plot(to_plot, 1.6, 5, ncol = 3) + 
    ylab("False discovery rate\n(-log10)") + 
    y_axis_all_facets() + 
    repellent_label(data = filter(to_plot, FDR < 0.05), aes(label = Gene), nudge_x = -Inf, direction = "y")

ggsave(file = "figures/images/Figure1A_volcano.png", dpi = 600, 
    width = 4.25, height = 2)

# significant probes
significant_probes <- get_probe_fold_changes("FmL") %>% filter(FDR < 0.05) %>% 
    .$ProbeId

significant_expression <- get_participant_level_probe_expression() %>% 
    filter(ProbeId %in% significant_probes) %>% 
    # simplify luteal/follicular/ring labels
    mutate(Visit = str_sub(Visit, 1, 1)) %>% 
    mutate(Visit = ifelse(Visit == "R", "C", Visit)) %>% 
    mutate(Visit = factor(Visit, levels = c("L", "F", "C")))

ggplot(significant_expression, aes(x = Visit, y = Expression, group = Visit)) + 
    geom_point(color = "grey75", size = pointSize()/2) + 
    geom_line(color = "grey75", aes(group = StudyId)) + 
    stat_summary(geom = "pointrange", fun.data = mean_cl_normal) + 
    theme_pub() + 
    ylab("Expression (log2)") + 
    facet_wrap(~Gene, scales = "free_y") + 
    x_axis_all_facets()

ggsave(file = "figures/images/Figure1B_expression.png", dpi = 600, 
    width = 4.25, height = 2.5)

# ddpcr -------------------------------------------------------------------
ggplot(get_ddpcr_vs_microarray(), aes(Expression, Log2GranulysinPerYWHAZ)) + 
    geom_point() + stat_smooth(method = "lm", color = highlight_line_color()) + 
    theme_pub() + 
    ylab("ddPCR (log2)") + 
    xlab("Microarray (log2)") +
    theme(aspect.ratio = 1) +
    facet_wrap(~ "Granulysin")

ggsave(file = "figures/images/Figure1C_ddpcr_granulysin_vs_microarray.png", dpi = 600, 
    width = 1.75, height = 1.75)

to_plot <- get_ddpcr() %>%     
    mutate(Visit = ifelse(Visit == "Ring", "CVR", as.character(Visit))) %>% 
    mutate(Visit = factor(Visit, levels = c("Luteal", "Follicular", "CVR"))) %>% 
    arrange(BV)

ggplot(to_plot, aes(x = Visit, y = Log2GranulysinPerYWHAZ, group = Visit)) + 
    geom_point(color = "grey75", size = pointSize()/2) + 
    geom_line(color = "grey75", aes(group = StudyId)) + 
    stat_summary(geom = "pointrange", fun.data = mean_cl_normal) + 
    xlab(NULL) +
    ylab("Granulysin ddPCR (log2)") + 
    scale_color_manual(values = c("grey50", "#e41a1c")) + 
    theme_pub()

ggsave(file = "figures/images/Figure1D_ddpcr_granulysin.png", dpi = 600, 
    width = 1.75, height = 2)

# hallmark ------------------------------------------------------------------
p_cutoff <- 0.01

to_plot <- get_hallmark("FmL") %>% 
    group_by(GeneSet) %>% 
    filter(sum(FDR < p_cutoff) > 1) %>% 
    mutate(Coefficient = str_replace(Coefficient, "Ring", "CVR")) %>% 
    # RvL, RvF, FvL
    mutate(Coefficient = factor(Coefficient, levels = 
            c("CVR vs. luteal", "CVR vs. follicular", "Follicular vs. luteal")))

to_plot %>% 
    clean_hallmark_gene_sets(1) %>%  
    hallmark_plot(p_cutoff) + guides(color = "none", fill = "none")

ggsave(file = "figures/images/Figure1E_hallmark.png", dpi = 600, 
    width = 5.65, height = 3)

###################################################################################################
# Figure 2 ------------------------------------------------------------------------
###################################################################################################

# protein ------------------------------------------------------------------------------------------------------
to_plot <- get_protein() %>%     
    mutate(Visit = ifelse(Visit == "Ring", "CVR", Visit)) %>% 
    mutate(Visit = factor(Visit, levels = c("Luteal", "Follicular", "CVR"))) %>% 
    arrange(BV, Hsv2)

ggplot(to_plot, aes(x = Visit, y = Log2Concentration)) + 
    geom_line(color = "grey75", size = 0.25, aes(group = StudyId)) + 
    geom_point(color = "grey75", size = pointSize()/2, aes(shape = Detectable)) + 
    scale_shape_manual(values = c(1, 16), guide = FALSE) + 
    stat_summary(geom = "pointrange", fun.data = mean_cl_normal) + 
    xlab(NULL) +
    ylab("Concentration (log2)") + 
    theme_pub() + facet_wrap(~ OfficialSymbol, scales = "free_y", ncol = 5) + 
    x_axis_all_facets() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave(file = "figures/images/Figure2a_protein_concentrations.png", dpi = 600, 
    width = 5, height = 4)

stats <- get_protein_stats() %>% 
    # rename to work with volcano_plot() function
    dplyr::rename(FDR = adjusted, Log2FoldChange = Value) %>% 
    # only the menstrual phase and ring visits have "vs" 
    filter(str_detect(Coefficient, "vs")) %>% 
    mutate(Coefficient = factor(Coefficient, levels = 
            c("Ring vs. luteal", "Ring vs. follicular", "Follicular vs. luteal")))

stats %>% 
    mutate(Significant = case_when(FDR < 0.05 ~ "Significant", TRUE ~ "Not")) %>% 
    ggplot(aes(x = Log2FoldChange, xmin = Log2FoldChange - 1.96*Std.Error, 
        xmax = Log2FoldChange + 1.96*Std.Error, y = reorder(OfficialSymbol, Log2FoldChange),
        color = Significant,
        fill = Significant)) + 
    geom_pointrange(pch = 21) + 
    facet_wrap(~Coefficient) + 
    ylab(NULL) + 
    xlab("Fold change (log2)") + 
    theme_pub() + 
    y_axis_all_facets() + 
    vline(0) + 
    scale_color_manual(values = c("Significant" = "black", "Not" = "grey50"), guide = FALSE) + 
    scale_fill_manual(values = c("Significant" = "black", "Not" = "white"), guide = FALSE) + 
    coord_cartesian(xlim = c(-6.5, 6.5))

ggsave(file = "figures/images/Figure2b_protein_stats.png", dpi = 600, 
    width = 4.25, height = 2.5)

###################################################################################################
# Figure 3 -----------------------------------------------------------------------------------
###################################################################################################
# other studies hallmark ------------------------------------------------------- 
to_plot <- get_other_studies_hallmark() %>% 
    group_by(GeneSet) %>% 
    mutate(Sig = sum(FDR < 0.05)) %>% 
    filter(Sig > 4) 

to_plot %>% 
    clean_hallmark_gene_sets(2) %>% 
    mutate(Category = factor(Category, levels = c("Immune", "Proliferation", "Other"))) %>% 
    mutate(Study = ifelse(Study == "Us", "Our study", Study)) %>% 
    mutate(StudyId = paste0(Organ, "\n", str_remove(Study, "\\-.*"))) %>% 
    mutate(StudyId = factor(StudyId, levels = c("Fallopian tube\nGSE28044", "Endometrium\nGSE29981",
        "Endometrium\nGSE6364", "Endometrium\nGSE86491", "Endocervix\nGSE80455", "Endocervix\nOur study", 
        "Endocervix\nGSE122248", "Ectocervix\nGSE122248", "Vagina\nPMID26750085"
    ))) %>% 
    hallmark_plot(transpose = FALSE) + 
        facet_grid(StudyId ~ Category, scales = "free_x", space = "free_x") +
        guides(color = "none", fill = "none") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1), 
            strip.text.y = element_text(angle = 0))

ggsave(file = "figures/images/Figure3A_hallmark.png", dpi = 600, 
    width = 5, height = 8)

###################################################################################################
# Figure 4 ------------------------------------------------------------------------
###################################################################################################
# volcano plots bv/hsv2 -------------------------------------------------------
# This plot is complicated... INPP5K has SUCH a different p-value than all the others that we need a broken
# axis, but ggplot doesn't support that well. I got around that problem by creating facets with separate axes,
# one facet for INPP5k and another for the rest. Then I facet without x-axes and it gives a broken axis
# Lots of weird code is needed for customizing the labels and stuff
to_plot <- bind_rows(
    get_gene_fold_changes("BV"),
    get_gene_fold_changes("Hsv2"))

to_plot <- to_plot %>% 
    # create separate facets for INPP5K and everything else
    mutate(facet = ifelse(FDR < 1E-5, "A", "B")) %>% 
    # pick out genes to highlight
    mutate(highlight = FDR < 1E-3 | 
            (Coefficient == "Bacterial vaginosis" & Log2FoldChange < 0 & FDR < 0.05) |
            (Coefficient == "HSV-2 infection" & Log2FoldChange > 1.46 & FDR < 0.03) | 
            (Coefficient == "HSV-2 infection" & Log2FoldChange < -0.59 & FDR < 0.05)) %>% 
    # the labels are going to be nudged and need to be in a specific order for that   
    mutate(NudgeOrder = case_when(
        Gene %in% c("BTF3L4", "CERNA2") ~ 1, 
        Gene %in% c("KRT6A", "LOC100132062") ~ 2,
        Gene == "INPP5K" ~ 3,
        Gene == "ARHGAP20" ~ 4, 
        Gene == "ZSWIM6" ~ 5,
        TRUE ~ as.numeric(NA)
        )) %>% 
    group_by(Coefficient) %>% 
    arrange(dplyr::desc(FDR))

# create dummy data frame that contains the axis bounds we want for each facet. 
# FDR of 1E-5 to 1 for lower facet, then 1E-12 to 1E-13 for upper
# will also use this data frame to put horizontal threshold line in lower but
# not upper facet, that's why the NAs
dummy <- data.frame(FDR = c(1, 1E-5, 1E-12, 1E-13), 
    facet = c("B", "B", "A", "A"), Log2FoldChange = 0, 
    Threshold = c(-log10(0.05), -log10(0.05), NA, NA))

# this function defines what gridlines/breaks should be used
# just want 12 and 13 for INPP5k, want 0 through 5 for the rest
custom_breaks <- function(x) {
    if (min(x) > 11) {
        c(12, 13)
    } else {
        0:5
    }
}

ggplot(to_plot, aes(x = Log2FoldChange, y = -log10(FDR), shape = FDR < 0.05, color = FDR < 0.05)) + 
    geom_point() +
    vline(0) + 
    xlim(c(-2, 2)) + 
    xlab("Fold change (log2)") + 
    ylab("False discovery\nrate (-log10)") + 
    # space = free sizes the facets based on the range of their y-axis
    facet_grid(facet ~ Coefficient, scales = "free_y", space = "free_y") + 
    scale_shape_manual(values = c(1, 16), guide = FALSE) + 
    scale_y_continuous(breaks = custom_breaks) + 
    scale_color_manual(values = c("grey50", "black"), guide = FALSE) + 
    y_axis_all_facets() +
    # hide the labels for the A/B facets
    theme_pub(strip.text.y = element_blank()) +
    # geom_blank sets the ylimits for each facet
    geom_blank(data = dummy) + 
    # this next line gives warnings about removing missing values, it's okay, just because of not having
    # lines for top facet
    geom_hline(data = dummy, aes(yintercept = Threshold)) +
    # add labels to select genes
    # nudge some to be at high end of x-axis, some at low end
    # nudge some up on y-axis
    # the nudging happens in order of the previously determined "NudgeOrder" which I put in to
    # make it possible to nudge the labels nicely
    repellent_label(data = filter(to_plot, highlight) %>% arrange(NudgeOrder), aes(label = Gene), 
        nudge_x = c(-Inf, Inf), nudge_y = c(0, 0.5), direction = "y", 
        min.segment.length = 0)
    
ggsave(file = "figures/images/Figure4A_bv_hsv2_volcano.png", dpi = 600, 
    width = 3.75, height = 2.25)

# plot of expression levels of selected DEGs ---------------------------------------------------------------
probe_expression <- get_participant_level_probe_expression() 

inpp5k <- probe_expression %>% filter(ProbeId == 1710438) %>% 
    mutate(Label = ifelse(BV == "Y", "BV+", "BV-"))

krt6a <- probe_expression %>% filter(Gene == "KRT6A") %>% 
    mutate(Label = ifelse(Hsv2 == "Positive", "HSV2+", "HSV2-"))

bind_rows(inpp5k, krt6a) %>% 
    ggplot(aes(x = Label, y = Expression, group = Label)) + 
        geom_point(color = "grey75", size = pointSize()/2) + 
        # only want connecting lines with INPP5K
        geom_line(data = inpp5k, color = "grey75", aes(group = StudyId)) + 
        stat_summary(geom = "pointrange", fun.data = mean_cl_normal) + 
        xlab(NULL) +
        ylab("Expression (log2)") + 
        scale_color_manual(values = c("grey50", "#e41a1c")) + 
        theme_pub() + 
        facet_wrap(~ Gene, scales = "free") + y_axis_all_facets()

ggsave(file = "figures/images/Figure4B_probe_expression.png", dpi = 600, 
    width = 2.5, height = 2)

# hallmark bv--------------------------------------------------
to_plot <- get_hallmark("BV") %>% 
    filter(FDR < 0.05) 

to_plot %>% 
    clean_hallmark_gene_sets(1) %>%  
    hallmark_plot() +
    guides(fill = "none", color = "none")

ggsave(file = "figures/images/Figure4C_left_bv_hallmark.png", dpi = 600, 
    width = 3, height = 3)

# hallmark hsv--------------------------------------------------
to_plot <- get_hallmark("HSV") %>% 
    filter(FDR < 0.05) 

to_plot %>% 
    mutate(Coefficient = "HSV-2 infection") %>% 
    clean_hallmark_gene_sets(2) %>%  
    mutate(Category = factor(Category, levels = c("Immune", "Signaling", "Other"))) %>% 
    hallmark_plot() +
    guides(fill = "none", color = "none") 

ggsave(file = "figures/images/Figure4C_right_hsv_hallmark.png", dpi = 600, 
    width = 3.5, height = 3)

# protein expression 
stats <- get_protein_stats() %>% 
    # rename to work with volcano_plot() function
    dplyr::rename(FDR = adjusted, Log2FoldChange = Value) %>% 
    # only the menstrual phase and ring visits have "vs" 
    filter(!str_detect(Coefficient, "vs")) 

bv_order <- stats %>% filter(Coefficient == "Bacterial vaginosis") %>% 
    arrange(Log2FoldChange) %>% pull(OfficialSymbol)

stats %>% 
    mutate(Significant = case_when(FDR < 0.05 ~ "Significant", TRUE ~ "Not")) %>% 
    mutate(OfficialSymbol = factor(OfficialSymbol, levels = bv_order)) %>% 
    ggplot(aes(x = Log2FoldChange, xmin = Log2FoldChange - 1.96*Std.Error, 
        xmax = Log2FoldChange + 1.96*Std.Error, y = OfficialSymbol,
        color = Significant,
        fill = Significant)) + 
    geom_pointrange(pch = 21) + 
    facet_wrap(~Coefficient) + 
    ylab(NULL) + 
    xlab("Fold change (log2)") + 
    theme_pub() + 
    y_axis_all_facets() + 
    vline(0) + 
    scale_color_manual(values = c("Significant" = "black", "Not" = "grey50"), guide = FALSE) + 
    scale_fill_manual(values = c("Significant" = "black", "Not" = "white"), guide = FALSE) + 
    coord_cartesian(xlim = c(-6.5, 6.5))

ggsave(file = "figures/images/Figure4D_protein_stats.png", dpi = 600, 
    width = 3.5, height = 2.5)

###################################################################################################
# Supplementary  Figure 1 ------------------------------------------------------------------------
###################################################################################################
# msd -----------------------------------------------------------------------------------
to_plot <- get_protein() %>% 
    mutate(BV = ifelse(BV == "N", "BV-", "BV+"), 
        Hsv2 = ifelse(Hsv2 == "Negative", "HSV-", "HSV+"))

ggplot(to_plot, aes(x = BV, y = Log2Concentration)) + 
    geom_line(color = "grey75", size = 0.25, aes(group = StudyId)) + 
    geom_point(color = "grey75", size = pointSize()/2) + 
    stat_summary(geom = "pointrange", fun.data = mean_cl_normal) + 
    xlab(NULL) +
    ylab("Cytokine concentration (log2)") + 
    theme_pub() + facet_wrap(~ OfficialSymbol, scales = "free_y", ncol = 3) + 
    x_axis_all_facets()

ggsave(file = "figures/images/Supplementary_Figure_1A_protein_concentrations_bv.png", dpi = 600, 
    width = 3, height = 5)

ggplot(to_plot, aes(x = Hsv2, y = Log2Concentration)) + 
    geom_point(color = "grey75", size = pointSize()/2) + 
    stat_summary(geom = "pointrange", fun.data = mean_cl_normal) + 
    xlab(NULL) +
    ylab("Cytokine concentration (log2)") + 
    theme_pub() + facet_wrap(~ OfficialSymbol, scales = "free_y", ncol = 3) + 
    x_axis_all_facets()

ggsave(file = "figures/images/Supplementary_Figure_1B_protein_concentrations_hsv.png", dpi = 600, 
    width = 3, height = 5)