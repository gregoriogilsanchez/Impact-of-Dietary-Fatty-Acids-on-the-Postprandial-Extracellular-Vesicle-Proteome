
# Script II: Clustering and Enrichment Preparation.
# Data preparation for VSClust and Metascape (online softwares).

rm(list = ls())
options(stringsAsFactors = FALSE)

library(ggplot2)
library(dplyr)
library(svglite)
library(qvalue)
library(stringr)

## Selecting samples from same dietary interventions.  
clean_data <- read.table("../datos_generados/clean_data.csv", header = TRUE, sep = ",")

PG_names = data [,1]
rownames(clean_data) = PG_names

NF2Vsclust = clean_data[, grep("NF", colnames(clean_data))]
write.csv(NF2Vsclust, file = "./datos_generados/NF2Vsclust.csv", row.names = T)

MUFA2Vsclust = clean_data[, grep("MUFA", colnames(clean_data))]
write.csv(MUFA2Vsclust, file = "./datos_generados/MUFA2Vsclust.csv", row.names = T)

SFA2Vsclust = clean_data[, grep("SFA", colnames(clean_data))]
write.csv(SFA2Vsclust, file = "./datos_generados/SFA2Vsclust.csv", row.names = T)

w32Vsclust = clean_data[, grep("W3", colnames(clean_data))]
write.csv(w32Vsclust, file = "./datos_generados/VSClust/w32Vsclust.csv", row.names = T)




# Selection of protein groups of each cluster -> Metascape



## NF

nf_vsclust = read.table(file = "./datos_generados/VSClust/VSClust_NF_result.tsv") 
## Data acquired via VSClust (online software) 


nf_clust_TRUE = nf_vsclust[nf_vsclust$isClusterMember == TRUE, ]

nf_cluster_1 = rownames(nf_clust_TRUE[nf_clust_TRUE$cluster == 1,]) 
nf_cluster_2 = rownames(nf_clust_TRUE[nf_clust_TRUE$cluster == 2,])
nf_cluster_3 = rownames(nf_clust_TRUE[nf_clust_TRUE$cluster == 3,])
nf_cluster_4 = rownames(nf_clust_TRUE[nf_clust_TRUE$cluster == 4,])
nf_cluster_5 = rownames(nf_clust_TRUE[nf_clust_TRUE$cluster == 5,])


## All Protein groups lists are required to have the same length to generate a data frame.
## Introducing NA values and replacing them with empty spaces

max_len = max(length(nf_cluster_1),length(nf_cluster_2),length(nf_cluster_3), 
              length(nf_cluster_4),length(nf_cluster_5))

length(nf_cluster_1) <- max_len
length(nf_cluster_2) <- max_len
length(nf_cluster_3) <- max_len
length(nf_cluster_4) <- max_len
length(nf_cluster_5) <- max_len

nf_vsc2metascape <- data.frame(
  cluster_1 = nf_cluster_1,
  cluster_2 = nf_cluster_2,
  cluster_3 = nf_cluster_3,
  cluster_4 = nf_cluster_4,
  cluster_5 = nf_cluster_5, stringsAsFactors = FALSE
)

nf_vsc2metascape[is.na(nf_vsc2metascape)] <- "" # NA values elimination


write.csv(nf_vsc2metascape, "./datos_generados/VSClust/nf_vsc2metascape.csv", row.names = F, quote = FALSE)


## SFA

sfa_vsclust = read.table(file = "./datos_generados/VSClust/SFA/VSClust_SFA_result.tsv")


sfa_clust_TRUE = sfa_vsclust[sfa_vsclust$isClusterMember == TRUE, ]

sfa_cluster_1 = rownames(sfa_clust_TRUE[sfa_clust_TRUE$cluster == 1,]) 
sfa_cluster_2 = rownames(sfa_clust_TRUE[sfa_clust_TRUE$cluster == 2,])
sfa_cluster_3 = rownames(sfa_clust_TRUE[sfa_clust_TRUE$cluster == 3,])
sfa_cluster_4 = rownames(sfa_clust_TRUE[sfa_clust_TRUE$cluster == 4,])


max_len = max(length(sfa_cluster_1),length(sfa_cluster_2),length(sfa_cluster_3), 
              length(sfa_cluster_4))

length(sfa_cluster_1) <- max_len
length(sfa_cluster_2) <- max_len
length(sfa_cluster_3) <- max_len
length(sfa_cluster_4) <- max_len

sfa_vsc2metascape <- data.frame(
  cluster_1 = sfa_cluster_1,
  cluster_2 = sfa_cluster_2,
  cluster_3 = sfa_cluster_3,
  cluster_4 = sfa_cluster_4,
   stringsAsFactors = FALSE
)

sfa_vsc2metascape[is.na(sfa_vsc2metascape)] <- "" 


write.csv(sfa_vsc2metascape, "./datos_generados/Metascape/sfa_vsc2metascape.csv", row.names = F, quote = FALSE)


## MUFA

mufa_vsclust = read.table(file = "./datos_generados/VSClust/MUFA/VSClust_MUFA_results.tsv")


mufa_clust_TRUE = mufa_vsclust[mufa_vsclust$isClusterMember == TRUE, ]

mufa_cluster_1 = rownames(mufa_clust_TRUE[mufa_clust_TRUE$cluster == 1,]) 
mufa_cluster_2 = rownames(mufa_clust_TRUE[mufa_clust_TRUE$cluster == 2,])
mufa_cluster_3 = rownames(mufa_clust_TRUE[mufa_clust_TRUE$cluster == 3,])
mufa_cluster_4 = rownames(mufa_clust_TRUE[mufa_clust_TRUE$cluster == 4,])
mufa_cluster_5 = rownames(mufa_clust_TRUE[mufa_clust_TRUE$cluster == 5,])


max_len = max(length(mufa_cluster_1),length(mufa_cluster_2),length(mufa_cluster_3), 
              length(mufa_cluster_4),length(mufa_cluster_5))

length(mufa_cluster_1) <- max_len
length(mufa_cluster_2) <- max_len
length(mufa_cluster_3) <- max_len
length(mufa_cluster_4) <- max_len
length(mufa_cluster_5) <- max_len

mufa_vsc2metascape <- data.frame(
  cluster_1 = mufa_cluster_1,
  cluster_2 = mufa_cluster_2,
  cluster_3 = mufa_cluster_3,
  cluster_4 = mufa_cluster_4,
  cluster_5 = mufa_cluster_5, stringsAsFactors = FALSE
)

mufa_vsc2metascape[is.na(mufa_vsc2metascape)] <- "" 


write.csv(mufa_vsc2metascape, "./datos_generados/Metascape/mufa_vsc2metascape.csv", row.names = F, quote = FALSE)


## w3

w3_vsclust = read.table(file = "./datos_generados/VSClust/W3/VSClust_w3_results.tsv")


w3_clust_TRUE = w3_vsclust[w3_vsclust$isClusterMember == TRUE, ]

w3_cluster_1 = rownames(w3_clust_TRUE[w3_clust_TRUE$cluster == 1,]) 
w3_cluster_2 = rownames(w3_clust_TRUE[w3_clust_TRUE$cluster == 2,])
w3_cluster_3 = rownames(w3_clust_TRUE[w3_clust_TRUE$cluster == 3,])
w3_cluster_4 = rownames(w3_clust_TRUE[w3_clust_TRUE$cluster == 4,])
w3_cluster_5 = rownames(w3_clust_TRUE[w3_clust_TRUE$cluster == 5,])


max_len = max(length(w3_cluster_1),length(w3_cluster_2),length(w3_cluster_3), 
              length(w3_cluster_4),length(w3_cluster_5))

length(w3_cluster_1) <- max_len
length(w3_cluster_2) <- max_len
length(w3_cluster_3) <- max_len
length(w3_cluster_4) <- max_len
length(w3_cluster_5) <- max_len

w3_vsc2metascape <- data.frame(
  cluster_1 = w3_cluster_1,
  cluster_2 = w3_cluster_2,
  cluster_3 = w3_cluster_3,
  cluster_4 = w3_cluster_4,
  cluster_5 = w3_cluster_5, stringsAsFactors = FALSE
)

w3_vsc2metascape[is.na(w3_vsc2metascape)] <- ""


write.csv(w3_vsc2metascape, "./datos_generados/Metascape/w3_vsc2metascape.csv", row.names = F, quote = FALSE)


#####################################################################################

# Metascape plots generation 

# MUFA
MUFAraw = read.csv("./MUFA_metascape_result/Enrichment_GO/_FINAL_GO.csv", check.names = F)

graph_list <- list()

colors_low = c(rep("#E1E5D5",5))
colors_high = c(rep("#556B2F",5))

for (i in 1:5) {
  
  # Generating column names: "_LogP_cluster_1", "_LogP_cluster_2", etc.
  col_logp <- paste0("_LogP_cluster_", i)
  
  # Data preparation
  df_cluster <- MUFAraw %>%
    # filtering only those with enrichment in this cluster.
    filter(.data[[col_logp]] < 0) %>%
    
    # Metascape provides log10(P). Obtaining real p value
    mutate(P_value = 10^(.data[[col_logp]])) %>%
    
    # Q value calculation
    mutate(Q_value = qvalue(P_value, pi0 = 1)$qvalue) %>% # Since this is a Metascape filtered list,
    # the q value distribution is not normal
    
    # Converting to  -log10(Q) 
    mutate(LogQ_pos = -log10(Q_value)) %>%
    
    # Top 15
    mutate(Description = str_wrap(Description, width = 40)) %>%
    arrange(desc(LogQ_pos)) %>%
    head(10)
  
  # Plotting
  p <- ggplot(df_cluster, aes(x = reorder(Description, LogQ_pos), y = -log10(Q_value), fill = LogQ_pos)) +
    geom_col(color = "black", width = 0.7, linewidth = 0.3) +
    coord_flip() +
    # Applying colors
    scale_fill_gradient(low = colors_low[i], high = colors_high[i], name = "-log10(q)", limits = c(0, NA)) + 
    labs(
      x = "", 
      y = "-log10(q-value)",
      title = paste("Top Enriched Terms: Cluster", i),
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.y = element_text(color = "black", size = 9, lineheight = 0.7), 
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  # Saving the results
  graph_list[[i]] <- p
  print(p)
  
  nombre_archivo <- paste0("MUFA_Enrichment_Cluster_", i, ".svg")
  ggsave(filename = nombre_archivo, plot = p, device = "svg", height = 5, width = 7, path = "./graficos_go/")    
}


## NF

NFraw = read.csv("./NF_metascape_result/Enrichment_GO/_FINAL_GO.csv", check.names = F)

graph_list <- list()

colors_low = c(rep("beige",5))
colors_high = c(rep("darkgrey",5))

for (i in 1:5) {
  
  col_logp <- paste0("_LogP_cluster_", i)
  
  df_cluster <- NFraw %>%
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
      title = paste("Top Enriched Terms: Cluster", i),
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.y = element_text(color = "black", size = 9, lineheight = 0.7), 
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  graph_list[[i]] <- p
  print(p)
  
  nombre_archivo <- paste0("NF_Enrichment_Cluster_", i, ".svg")
  ggsave(filename = nombre_archivo, plot = p, device = "svg", height = 5, width = 7, path = "./graficos_go/")    
}


## SFA

SFAraw = read.csv("./SFA_metascape_result/Enrichment_GO/_FINAL_GO.csv", check.names = F)

graph_list <- list()

colors_low = c(rep("salmon",5))
colors_high = c(rep("red",5))

for (i in 1:4) {
  
  col_logp <- paste0("_LogP_cluster_", i)
  
  df_cluster <- SFAraw %>%
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
      title = paste("Top Enriched Terms: Cluster", i),
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.y = element_text(color = "black", size = 9, lineheight = 0.7), 
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  graph_list[[i]] <- p
  print(p)
  
  nombre_archivo <- paste0("SFA_Enrichment_Cluster_", i, ".svg")
  ggsave(filename = nombre_archivo, plot = p, device = "svg", height = 5, width = 7, path = "./graficos_go/")    
}


## W3

W3raw = read.csv("./w3_metascape_result/Enrichment_GO/_FINAL_GO.csv", check.names = F)

graph_list <- list()

colors_low = c(rep("#7FFFD4",5))
colors_high = c(rep("darkblue",5))

for (i in 1:5) {
  
  col_logp <- paste0("_LogP_cluster_", i)
  
  df_cluster <- W3raw %>%
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
      title = paste("Top Enriched Terms: Cluster", i),
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.y = element_text(color = "black", size = 9, lineheight = 0.7), 
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  graph_list[[i]] <- p
  print(p)
  
  nombre_archivo <- paste0("W3_Enrichment_Cluster_", i, ".svg")
  ggsave(filename = nombre_archivo, plot = p, device = "svg", height = 5, width = 7, path = "./graficos_go/")    
}
