# shows a volcano plot with color and shape determined by FDR < 0.05
volcano_plot <- function(d, xmax, ymax, ncol = 1) {
    d <- d %>% 
        group_by(Coefficient) %>% 
        arrange(dplyr::desc(FDR))
    
    p <- ggplot(d, aes(x = Log2FoldChange, y = -log10(FDR), shape = FDR < 0.05, color = FDR < 0.05))
    
    volcano_plot_helper(p, xmax, ymax, ncol) + 
        scale_color_manual(values = c("grey50", "black"), guide = FALSE)
}

# shows a volcano plot, using color to highlight the genes in the 
# provided hallmark gene set
# highlights should be of the form HALLMARK_INTERFERON_GAMMA_RESPONSE
highlight_volcano_plot <- function(d, xmax, ymax, highlights) {
    members <- get_hallmark_gene_sets()[[highlights]]
    
    to_plot <- d %>% 
        mutate(GeneSet = ifelse(EntrezId %in% members, "In", "Out")) %>% 
        arrange(dplyr::desc(GeneSet))
    
    p <- ggplot(to_plot, aes(x = Log2FoldChange, y = -log10(FDR), shape = FDR < 0.05, color = GeneSet))
    
    volcano_plot_helper(p, xmax, ymax) + 
        scale_color_manual(values = c("black", "grey50"))
}

# helper function called by the other volcano plot functions
volcano_plot_helper <- function(p, xmax, ymax, ncol = 1) {
    p + geom_point() +
        hline(-log10(0.05)) + 
        vline(0) + 
        ylim(c(0, ymax)) + 
        xlim(c(-xmax, xmax)) + 
        xlab("Fold change (log2)") + 
        ylab("False discovery\nrate (-log10)") + 
        facet_wrap(~ Coefficient, ncol = ncol) + 
        scale_shape_manual(values = c(1, 16), guide = FALSE) + 
        theme_pub() + 
        x_axis_all_facets()
}

hallmark_plot <- function(d, threshold = 0.05, transpose = TRUE) {
    # x and y are reversed with coord_flip so xlab, ylab, etc are all switched
    p <- ggplot(d, aes(y = -log10(FDR), x = reorder(GeneSet, FDR), fill = Direction, 
            color = Direction, alpha = FDR < threshold)) + 
        geom_col(width = 0.8) + 
        scale_fill_manual(values = c("Down" = "#0072b2", "Up" = "#d95f02")) + 
        scale_color_manual(values = c("Down" = "#0072b2", "Up" = "#d95f02")) + 
        scale_alpha_manual(values = c("FALSE" = 0, "TRUE" = 1), guide = FALSE) + 
        hline(-log10(threshold)) + 
        facet_grid(Category ~ Coefficient, scales = "free_y", space = "free_y") + 
        xlab(NULL) + ylab("False discovery rate (-log10)") + 
        theme_pub() + y_axis_all_facets()
    
    if (transpose) {
        p + coord_flip()
    } else {
        p
    }
}

repellent_label <- function(...) {
    ggrepel::geom_label_repel(..., size = 2.5, label.padding = 0.1, point.padding = 0.1, box.padding = 0.1, 
        segment.size = 0.2, segment.color = "black")
}
