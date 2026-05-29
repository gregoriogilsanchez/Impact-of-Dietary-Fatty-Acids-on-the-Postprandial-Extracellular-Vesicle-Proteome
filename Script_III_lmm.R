# Script III: Differential abundance analysis
# LIMMA package - Application of linear mixed models (LMM).

rm(list = ls())
options(stringsAsFactors = FALSE)

library(limma)
library(calibrate)
library(ggplot2)
library(ggrepel)
library(scales)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(svglite)
library(qvalue)
library(stringr)


## Data loading and formatting.

data_limma = read.csv("../datos_generados/clean_data.csv", header = TRUE)
PG_names = data_limma[,1]
rownames(data_limma) = PG_names
data_limma = data_limma[,-1]

## Experimental design matrix.

## Load the design matrix generated during preprocessing.  
design = read.csv("../datos_generados/design.csv", header = TRUE)
design_names = design[,1]
rownames(design) = design_names
design = design[,-1]
design$sample <- paste0("X", design$sample) #Format adaptation. 

## Adaptation to LIMMA format. 
design$group <- factor(design$group)
design_limma <- model.matrix(~ 0 + design$group)
colnames(design_limma) <- levels(design$group)
rownames(design_limma) <- design$sample
head(design_limma, 15)


## Limma expects a logarithmic transformation of the data.

datos_log = log2(data_limma)


corfit <- duplicateCorrelation(datos_log, design_limma, block = design$replicas) # intra-subject correlation as a random effect.
linear.fit <- lmFit(datos_log, design_limma, block = design$replicas, correlation = corfit$consensus.correlation)

head(linear.fit$coefficients) # Mean expression matrix


# Define the contrast matrix. B as the reference condition, NF as the reference diet.

contrast.matrix = makeContrasts(
  
  # --- COMPARISONS AT TIME POINT P ---
  
  SFA_NF_P = (SFAP - SFAB) - (NFP - NFB),
  MUFA_NF_P = (MUFAP - MUFAB) - (NFP - NFB),
  W3_NF_P = (W3P - W3B) - (NFP - NFB),
  
  # --- COMPARISONS AT TIME POINT PP ---
  SFA_NF_PP = (SFAPP - SFAB) - (NFPP - NFB),
  MUFA_NF_PP = (MUFAPP - MUFAB) - (NFPP - NFB),
  W3_NF_PP = (W3PP - W3B) - (NFPP - NFB),
  
  levels = c("MUFAB","MUFAP","MUFAPP","NFB","NFP","NFPP","SFAB","SFAP","SFAPP",
             "W3B","W3P","W3PP")
)


contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
contrast.results <- eBayes(contrast.linear.fit) 

# Results extraction.

# Due to the nature of the contrasts:
# Activated: expression increases (or decreases less)
# in SFAP compared to SFAB more than in NFP compared to NFB.

# Repressed: expression decreases (or increases less)
# in SFAP compared to SFAB more than in NFP compared to NFB.

contrasts = list("SFA_NF_P" = 1,
                 "MUFA_NF_P" = 2,
                 "ω3_NF_P" = 3,
                 "SFA_NF_PP" = 4,
                 "MUFA_NF_PP" = 5,
                 "ω3_NF_PP"  = 6)

for (contrast_name in names(contrasts)) {
  coef_idx <- contrasts[[contrast_name]]
  
  result <- topTable(contrast.results,
                     number = nrow(datos_log),
                     coef   = coef_idx,
                     sort.by = "logFC")
  
  limma.log2fc   <- result$logFC
  limma.pval     <- result$P.Value # nominal p value
  limma.prots.ids <- rownames(result)
  
  names(limma.log2fc) <- limma.prots.ids
  names(limma.pval) <- limma.prots.ids
  
  activated <- limma.prots.ids[limma.log2fc >  0.5849625 & limma.pval < 0.05] # Fold Change 1.5
  repressed  <- limma.prots.ids[limma.log2fc < -0.5849625 & limma.pval < 0.05]
  
  activated <- unique((na.omit(unlist(strsplit(activated, ";")))))
  repressed  <- unique((na.omit(unlist(strsplit(repressed,  ";")))))
  
  cat(contrast_name, "-> Activated:", length(activated),
      "| Repressed:", length(repressed), "\n")
  
  #write.table(activated, file = paste0("./results/", contrast_name, "_activated.txt"),
  #            quote = FALSE, row.names = FALSE, col.names = FALSE)
  #write.table(repressed,  file = paste0("./results/", contrast_name, "_repressed.txt"),
  #            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Volcano plot 

  #png(paste0("results/", contrast_name, "_volcano.png"), width = 800, height = 600)
  
  plot(limma.log2fc, -log10(limma.pval),
       main = contrast_name,
       xlab = "log2 Fold Change",
       ylab = "-log10(p-value)",
       pch  = 20, col = "grey")
  
  points(limma.log2fc[limma.log2fc >  0.5849625 & limma.pval < 0.05],
         -log10(limma.pval[limma.log2fc >  0.5849625 & limma.pval < 0.05]),
         pch = 20, col = "red")
  
  points(limma.log2fc[limma.log2fc < -0.5849625 & limma.pval < 0.05],
         -log10(limma.pval[limma.log2fc < -0.5849625 & limma.pval < 0.05]),
         pch = 20, col = "blue")
  
  abline(v = c(-0.5849625, 0.5849625), h = -log10(0.05),
         lty = 2, col = "black")
  
  dev.off()
  

}


#ID conversion -> Volcano plot with protein names

for (contrast_name in names(contrasts)) {
  coef_idx <- contrasts[[contrast_name]]
  
  result <- topTable(contrast.results,
                     number = nrow(datos_log),
                     coef   = coef_idx,
                     sort.by = "logFC")
  
  # data preparation
  plot_data <- data.frame(
    Protein = rownames(result),
    logFC   = as.numeric(result$logFC),
    pval    = as.numeric(result$P.Value),
    stringsAsFactors = FALSE
  )
  
  plot_data$logP <- -log10(plot_data$pval)
  plot_data <- plot_data[!is.na(plot_data$logFC) & !is.na(plot_data$logP), ]
  
  # Significance classification
  plot_data$Significance <- "Not Significant"
  plot_data$Significance[plot_data$logFC > 0.5849625 & plot_data$pval < 0.05] <- "Up"
  plot_data$Significance[plot_data$logFC < -0.5849625 & plot_data$pval < 0.05] <- "Down"
  
  num_up <- sum(plot_data$Significance == "Up")
  num_down <- sum(plot_data$Significance == "Down")
  
  cat(contrast_name, "-> Up:", num_up, "| Down:", num_down, "\n")
  
  # Selecting top 10 for labelling
  up_prots <- plot_data[plot_data$Significance == "Up", ]
  top10_up <- up_prots[order(up_prots$pval), "Protein"][1:10]
  
  #down_prots <- plot_data[plot_data$Significance == "Down", ]
  #top10_down <- down_prots[order(down_prots$pval), "Protein"][1:10]
  
  #plot_data$Label <- ""
  #proteins_to_label <- na.omit(c(top10_up, top10_down))
  #plot_data$Label[plot_data$Protein %in% proteins_to_label] <- plot_data$Protein[plot_data$Protein %in% proteins_to_label]
  
  # Mathematical transformation to shrink the lower quadrant
  trans_comprimida <- trans_new(
    name = "comprimir_bajos",
    transform = function(x) ifelse(x <= 1.30103, x * 0.2, 1.30103 * 0.2 + (x - 1.30103) * 1.8),
    inverse = function(x) ifelse(x <= 1.30103 * 0.2, x / 0.2, 1.30103 + (x - 1.30103 * 0.2) / 1.8)
  )
  

    # plotting
  p <- ggplot(plot_data, aes(x = logFC, y = logP, color = Significance)) +
    
    # Dot layer
    geom_point(alpha = 0.7, size = 2) +
    
    # threshold lines
    geom_vline(xintercept = c(-0.5849625, 0.5849625), linetype = "dashed", color = "black", alpha = 0.5) +
    geom_hline(yintercept = 1.30103, linetype = "dashed", color = "black", alpha = 0.5) +
    
    # pallette
    scale_color_manual(values = c("Up" = "darkred", "Down" = "darkblue", "Not Significant" = "grey75")) +
    
    # Up and Down labels
    annotate("label", x = min(plot_data$logFC) * 0.8, y = max(plot_data$logP) * 0.95, 
             label = paste0("Down: ", num_down), fill = "darkblue", color = "white", 
             fontface = "bold", size = 3.5, label.padding = unit(0.3, "lines")) +
    
    annotate("label", x = max(plot_data$logFC) * 0.8, y = max(plot_data$logP) * 0.95, 
             label = paste0("Up: ", num_up), fill = "darkred", color = "white", 
             fontface = "bold", size = 3.5, label.padding = unit(0.3, "lines")) +
    
    # Shrinked and continous scale for axis Y Y
    scale_y_continuous(trans = trans_comprimida, breaks = pretty_breaks(n = 6)) +
    
    # Format and labels
    labs(title = contrast_name, x = "log2 Fold Change", y = "-log10(p-value)") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      axis.title = element_text(face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      legend.position = "none" 
    )
  
  # Saving final plots in SVG format
  ggsave(paste0("results/", contrast_name, "_volcano.svg"), plot = p, width = 7, height = 6, dpi = 300)
}

###########################################################################

## Functional enrichment after LMM analysis -> Plotting Metascape Results 

SFA_NF_P_rep = read.csv("../../LMM/results/cytoscape/SFA_NF_P_rep/Enrichment_GO/_FINAL_GO.csv", check.names = F)

graph_list <- list()

colors_low = c(rep("salmon",5))
colors_high = c(rep("red",5))

for (i in 1) {
  
  col_logp <- paste0("_LogP_hits") #Metascape Nomenclature
  
  # Data preparation
  df_cluster <- SFA_NF_P_rep %>%
    # filtering only those with enrichment
    filter(.data[[col_logp]] < 0) %>%
    
    # Metascape provides log10(P).  Obtaining real p value
    mutate(P_value = 10^(.data[[col_logp]])) %>%
    
    # q value
    mutate(Q_value = qvalue(P_value, pi0 = 1)$qvalue) %>% # Since this is a Metascape filtered list,
    # the q value distribution is not normal
    
    # Converting to -log10(Q) 
    mutate(LogQ_pos = -log10(Q_value)) %>%
    
    #  Top 15
    mutate(Description = str_wrap(Description, width = 40)) %>%
    arrange(desc(LogQ_pos)) %>%
    head(10)
  
  # Plotting
  p <- ggplot(df_cluster, aes(x = reorder(Description, LogQ_pos), y = -log10(Q_value), fill = LogQ_pos)) +
    geom_col(color = "black", width = 0.7, linewidth = 0.3) +
    coord_flip() +
    # Colors
    scale_fill_gradient(low = colors_low[i], high = colors_high[i], name = "-log10(q)", limits = c(0, NA)) + 
    labs(
      x = "", 
      y = "-log10(q-value)",
      title = paste("SFA_NF_P_rep"),
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.y = element_text(color = "black", size = 9, lineheight = 0.7), 
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  # Saving the results
  graph_list[[i]] <- p
  print(p)
  
  nombre_archivo <- paste0("SFA_NF_P_rep", ".svg")
  ggsave(filename = nombre_archivo, plot = p, device = "svg", height = 5, width = 7, path = "../../LMM/results/cytoscape/SFA_NF_P_rep/")    
}

## SFA_NF_PP_act

SFA_NF_PP_act = read.csv("../../LMM/results/cytoscape/SFA_NF_PP_act/Enrichment_GO/_FINAL_GO.csv", check.names = F)

graph_list <- list()


colors_low = c(rep("salmon",5))
colors_high = c(rep("red",5))

for (i in 1) {
  
  col_logp <- paste0("_LogP_hits")
  
  df_cluster <- SFA_NF_PP_act %>%
    filter(.data[[col_logp]] < 0) %>%
    
    mutate(P_value = 10^(.data[[col_logp]])) %>%
    
    mutate(Q_value = qvalue(P_value, pi0 = 1)$qvalue) %>% 
    
    mutate(LogQ_pos = -log10(Q_value)) %>%
    
    mutate(Description = str_wrap(Description, width = 40)) %>%
    arrange(desc(LogQ_pos)) %>%
    head(10)
  
  p <- ggplot(df_cluster, aes(x = reorder(Description, LogQ_pos), y = -log10(Q_value), fill = LogQ_pos)) +
    geom_col(color = "black", width = 0.7, linewidth = 0.3) +
    coord_flip() +
    
    scale_fill_gradient(low = colors_low[i], high = colors_high[i], name = "-log10(q)", limits = c(0, NA)) + 
    labs(
      x = "", 
      y = "-log10(q-value)",
      title = paste("SFA_NF_PP_act"),
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.y = element_text(color = "black", size = 9, lineheight = 0.7), 
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  
  graph_list[[i]] <- p
  print(p)
  
  nombre_archivo <- paste0("SFA_NF_PP_act", ".svg")
  ggsave(filename = nombre_archivo, plot = p, device = "svg", height = 5, width = 7, path = "../../LMM/results/cytoscape/SFA_NF_PP_act/")    
}


## MUFA_NF_PP_rep

MUFA_NF_PP_rep = read.csv("../../LMM/results/cytoscape/MUFA_NF_PP_rep/Enrichment_GO/_FINAL_GO.csv", check.names = F)

graph_list <- list()


colors_low = c(rep("#E1E5D5",5))
colors_high = c(rep("#556B2F",5))

for (i in 1) {
  
  col_logp <- paste0("_LogP_hits")
  
  df_cluster <- MUFA_NF_PP_rep %>%
    filter(.data[[col_logp]] < 0) %>%
    
    mutate(P_value = 10^(.data[[col_logp]])) %>%
    
    mutate(Q_value = qvalue(P_value, pi0 = 1)$qvalue) %>% 
    
    mutate(LogQ_pos = -log10(Q_value)) %>%
    
    mutate(Description = str_wrap(Description, width = 40)) %>%
    arrange(desc(LogQ_pos)) %>%
    head(10)
  
  p <- ggplot(df_cluster, aes(x = reorder(Description, LogQ_pos), y = -log10(Q_value), fill = LogQ_pos)) +
    geom_col(color = "black", width = 0.7, linewidth = 0.3) +
    coord_flip() +
    scale_fill_gradient(low = colors_low[i], high = colors_high[i], name = "-log10(q)", limits = c(0, NA)) + 
    labs(
      x = "", 
      y = "-log10(q-value)",
      title = paste("MUFA_NF_PP_rep"),
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.y = element_text(color = "black", size = 9, lineheight = 0.7), 
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  graph_list[[i]] <- p
  print(p)
  
  nombre_archivo <- paste0("MUFA_NF_PP_rep", ".svg")
  ggsave(filename = nombre_archivo, plot = p, device = "svg", height = 5, width = 7, path = "../../LMM/results/cytoscape/MUFA_NF_PP_rep/")    
}

## W3_NF_PP_act

W3_NF_PP_act = read.csv("../../LMM/results/cytoscape/W3_NF_PP_act/Enrichment_GO/_FINAL_GO.csv", check.names = F)

graph_list <- list()


colors_low = c(rep("#7FFFD4",5))
colors_high = c(rep("darkblue",5))

for (i in 1) {
  
  col_logp <- paste0("_LogP_hits")
  
  df_cluster <- W3_NF_PP_act %>%
    filter(.data[[col_logp]] < 0) %>%
    
    mutate(P_value = 10^(.data[[col_logp]])) %>%
    
    mutate(Q_value = qvalue(P_value, pi0 = 1)$qvalue) %>% 
    
    mutate(LogQ_pos = -log10(Q_value)) %>%
    
    mutate(Description = str_wrap(Description, width = 40)) %>%
    arrange(desc(LogQ_pos)) %>%
    head(10)
  
  p <- ggplot(df_cluster, aes(x = reorder(Description, LogQ_pos), y = -log10(Q_value), fill = LogQ_pos)) +
    geom_col(color = "black", width = 0.7, linewidth = 0.3) +
    coord_flip() +
    scale_fill_gradient(low = colors_low[i], high = colors_high[i], name = "-log10(q)", limits = c(0, NA)) + 
    labs(
      x = "", 
      y = "-log10(q-value)",
      title = paste("W3_NF_PP_act"),
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.y = element_text(color = "black", size = 9, lineheight = 0.7), 
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  graph_list[[i]] <- p
  print(p)
  
  nombre_archivo <- paste0("W3_NF_PP_act", ".svg")
  ggsave(filename = nombre_archivo, plot = p, device = "svg", height = 5, width = 7, path = "../../LMM/results/cytoscape/W3_NF_PP_act/")    
}