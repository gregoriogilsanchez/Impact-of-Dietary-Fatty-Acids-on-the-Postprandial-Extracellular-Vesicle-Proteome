
# Script I: Data preprocessing and preliminary analysis (PCA, correlation analysis,
# Z-score calculation, etc). 

rm(list = ls())
options(stringsAsFactors = FALSE)

library(dplyr)
library(stringr)
library(ggplot2)
library(svglite)
library("FactoMineR")
library("factoextra")
library(corrplot)
library(ComplexHeatmap)
library(circlize)
library(UniProt.ws)
library(enrichR)
library(qvalue)

## Data loading and formatting.
data <- read.delim("../datos_brutos/2024037-MagNetOT2_EV-enrichment_Report.tsv", header = TRUE)

PG_names = data [,1]
rownames(data) = PG_names

data = data[,-c(1:4)]


samplenames = c("41MUFAB", "41MUFAP", "41MUFAPP", "41NFB", "41NFP", 
                     "41NFPP", "41SFAB", "41SFAP", "41SFAPP", "41W3B", "41W3P", 
                     "41W3PP", "42MUFAB", "42MUFAP", "42MUFAPP", "42NFB", "42NFP", 
                     "42NFPP", "42SFAB", "42SFAP", "42SFAPP", "42W3B", "42W3P", 
                     "42W3PP", "43MUFAB", "43MUFAP", "43MUFAPP", "43NFB", "43NFP", 
                     "43NFPP", "43SFAB", "43SFAP", "43SFAPP", "43W3B", "43W3P", 
                     "43W3PP", "44MUFAB", "44MUFAP", "44MUFAPP", "44NFB", "44NFP", 
                     "44NFPP", "44SFAB", "44SFAP", "44SFAPP", "44W3B", "44W3P", 
                     "44W3PP", "45MUFAB", "45MUFAP", "45MUFAPP", "45NFB", "45NFP", 
                     "45NFPP", "45SFAB", "45SFAP", "45SFAPP", "45W3B", "45W3P", 
                     "45W3PP", "46MUFAB", "46MUFAP", "46MUFAPP", "46NFB", "46NFP", 
                     "46NFPP", "46SFAB", "46SFAP", "46SFAPP", "46W3B", "46W3P", 
                     "46W3PP", "47MUFAB", "47MUFAP", "47MUFAPP", "47NFB", "47NFP", 
                     "47NFPP", "47SFAB", "47SFAP", "47SFAPP", "47W3B", "47W3P", 
                     "47W3PP", "48MUFAB", "48MUFAP", "48MUFAPP", "48NFB", "48NFP", 
                     "48NFPP", "48SFAB", "48SFAP", "48SFAPP", "48W3B", "48W3P", 
                     "48W3PP", "49MUFAB", "49MUFAP", "49MUFAPP", "49NFB", "49NFP", 
                     "49NFPP", "49SFAB", "49SFAP", "49SFAPP", "49W3B", "49W3P", 
                     "49W3PP", "50MUFAB", "50MUFAP", "50MUFAPP", "50NFB", "50NFP", 
                     "50NFPP", "50SFAB", "50SFAP", "50SFAPP", "50W3B", "50W3P", "50W3PP")

colnames(data) = samplenames

data[1:6,1:6]

## Adjusting individuals names: 41:50 -> 1:10
samplenames <- colnames(data)
for (i in 41:50) {
  samplenames <- gsub(pattern = as.character(i), 
                  replacement = as.character(i - 40), 
                  x = samplenames)
}
colnames(data) <- samplenames

## Defining the experimental design. 

experimentalgroups <- gsub("^[0-9]{2}", "", samplenames)

replicas <- sub("^([0-9]{2}).*", "\\1", samplenames)

design <- data.frame(sample=colnames(data),
                     group=experimentalgroups,
                     replicas)
design

#write.csv(design, file = "./datos_generados/design.csv", row.names = T)


## Filtering proteins with a high proportion of missing values

unique_groups <- unique(design$group)

valid_proteins <- rep(FALSE, nrow(data)) # Indicates TRUE for proteins that pass the filter


for (g in unique_groups) { #Checking all 12 experimental conditions
  
  grouop_samples <- design$sample[design$group == g]
  
  data_sub <- data[, grouop_samples]
  
  # Counting the number of not missing values in the 12 samples. 
  detected_values <- rowSums(!is.na(data_sub))
  
  # Passes the filter if it has at least 6 non-missing values in at least one experimental condition
  follows_rule <- detected_values >= 6
  
  # Update the vector of proteins that pass the filter
  valid_proteins <- valid_proteins | follows_rule
}


clean_data <- data[valid_proteins, ]

nrow(clean_data)  

#write.csv(clean_data, file = "./datos_generados/clean_data.csv", row.names = T)

## Data characterization. 

## Histogram showing the number of identified proteins (PG) in each sample
num_PG =  colSums(!is.na(clean_data))


df_PG <- data.frame(
  Sample = names(num_PG),
  PGs = as.numeric(num_PG)
)

df_PG$Diet <- ifelse(grepl("NF", df_PG$Sample), "NF",
                     ifelse(grepl("SFA", df_PG$Sample), "SFA",
                            ifelse(grepl("MUFA", df_PG$Sample), "MUFA",
                                   ifelse(grepl("W3", df_PG$Sample), "W3","other"))))

df_PG$Diet <- factor(df_PG$Diet, levels = c("NF", "SFA", "MUFA", "W3")) #Order by diet
df_PG <- df_PG[order(df_PG$Diet), ]
df_PG$Sample <- factor(df_PG$Sample, levels = df_PG$Sample)

my_colors <- c(
  "NF" = "grey70",       
  "SFA" = "salmon",      
  "MUFA" = "olivedrab",  
  "W3" = "lightblue")    


PGs_plot <- ggplot(df_PG, aes(x = Sample, y = PGs, fill = Diet)) + #Plot details
  geom_col(color = "black", width = 0.75, linewidth = 0.3) + 
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(limits = c(0, 2400), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Detected protein groups per sample",
    x = "Sample",
    y = "Number of Protein Groups",
    fill = "Diet Intervention" 
  ) +
  theme_classic(base_size = 12) + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8), 
    axis.title = element_text(face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    legend.position = "top" 
  )

print(PGs_plot)

ggsave("./PGperSample.svg", plot = PGs_plot, width = 12, height = 6)

## Histogram showing the number of missing values per sample

num_na =  colSums(is.na(clean_data))


df_na <- data.frame(
  Sample = names(num_na),
  nas = as.numeric(num_na)
)

df_na$Diet <- ifelse(grepl("NF", df_na$Sample), "NF",
                     ifelse(grepl("SFA", df_na$Sample), "SFA",
                            ifelse(grepl("MUFA", df_na$Sample), "MUFA",
                                   ifelse(grepl("W3", df_na$Sample), "W3","other"))))

df_na$Diet <- factor(df_na$Diet, levels = c("NF", "SFA", "MUFA", "W3")) #Order by diet
df_na <- df_na[order(df_na$Diet), ]
df_na$Sample <- factor(df_na$Sample, levels = df_na$Sample)

my_colors <- c(
  "NF" = "grey70",       
  "SFA" = "salmon",      
  "MUFA" = "olivedrab",  
  "W3" = "lightblue")    


nas_plot <- ggplot(df_na, aes(x = Sample, y = nas, fill = Diet)) + #Plot details
  geom_col(color = "black", width = 0.75, linewidth = 0.3) + 
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(limits = c(0, 2400), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Missing values per sample",
    x = "Sample",
    y = "Number of Missing values",
    fill = "Diet Intervention" 
  ) +
  theme_classic(base_size = 12) + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8), 
    axis.title = element_text(face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    legend.position = "top" 
  )

print(nas_plot)

ggsave("./naperSample.svg", plot = nas_plot, width = 12, height = 6)


#Brief data distribution analysis. 

boxplot(clean_data, outline = F)

## Impute missing values (NaNs) replicating the Perseus algorithm (to perform PCA analysis)

#Nan -> # Value similar to the lowest  found in the sample
datoslog = log2(clean_data+1) #Perseus algorithm requieres normality

impute_perseus <- function(column, shift = 1.8, width = 0.3) {
  is_na <- is.na(column)
  if(!any(is_na)) return(column) # If there are not missing values, return the exact column. 
  
  # Compute the mean and standard deviation of detected proteins
  media <- mean(column, na.rm = TRUE)
  st_dev <- sd(column, na.rm = TRUE)
  
  # Shift the distribution downward (by 1.8 standard deviations)
  new_mean <- media - (shift * st_dev)
  new_dev <- width * st_dev
  
  # Fill missing values (NAs) with randomly generated low values
  column[is_na] <- rnorm(sum(is_na), mean = new_mean, sd = new_dev)
  return(column)
}

# lapply to apply the previous function  per column (sample by sample).
imputed_data <- as.data.frame(lapply(datoslog, impute_perseus))

rownames(imputed_data) <- rownames(datos_limpios)
colnames(imputed_data) <- colnames(datos_limpios)

# Reverse the log transformation
imputed_data = 2 ^ imputed_data

#write.csv(imputed_data, file = "./datos_generados/imputed_data.csv", row.names = T)




############################################################
#################      PCA        ##########################
############################################################

imputed_data <- read.table("../datos_generados/imputed_data.csv", header = TRUE, sep = ",")
colnames(imputed_data) = c("X",colnames(clean_data))

# Performing PCA analysis with imputed data to avoid predefined PCA imputation by the mean. 


# MUFA

data_mufa <- imputed_data %>%
  select(contains("MUFA"))

data_mufa = t(data_mufa) 


head(data_mufa)
dim(data_mufa)

data_mufa = as.data.frame(data_mufa)


data_mufa$qualysup <- as.factor(str_extract(rownames(data_mufa), "B$|P$|PP$"))

table(data_mufa$qualysup) #Check samples distribution

pca.mufa = PCA(data_mufa, scale.unit = T, graph = F, quali.sup = ncol(data_mufa)) #Perform PCA


MUFA_timepoint <- fviz_pca_ind( #Plotting
  pca.mufa, 
  axes = c(1, 2), 
  geom = "point",                   
  pointsize = 3,                    
  pointshape = 19,                  
  habillage = ncol(data_mufa),      
  palette = c("#0073C2", "#EFC000", "pink"), 
  addEllipses = TRUE, 
  ellipse.type = "confidence", 
  ellipse.level = 0.95,
  ellipse.alpha = 0.15,             
  mean.point = FALSE,               
  title = "MUFA Diet by Time Point", 
  legend.title = "Time Point"
) +
  
  theme_classic(base_size = 14) +   
  theme(
    legend.position = "top",        
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1) 
  )


print(MUFA_timepoint)

ggsave("./MUFA_timepoint.svg", plot = MUFA_timepoint, width = 8, height = 6)

## MUFA Group subject replicates


data_mufa$qualysup <- str_extract(rownames(data_mufa), "^\\d+")
data_mufa$qualysup <- as.factor(data_mufa$qualysup)

pca.mufa = PCA(data_mufa, scale.unit = T, graph = F, quali.sup = ncol(data_mufa)) 
table(data_mufa$qualysup)

colores <- c(
  "#1f77b4", 
  "#ff7f0e", 
  "#2ca02c", 
  "#d62728", 
  "#9467bd", 
  "#8c564b", 
  "#e377c2", 
  "#7f7f7f", 
  "#bcbd22", 
  "#17becf"  
)

MUFA_subject <- fviz_pca_ind(
  pca.mufa, 
  axes = c(1, 2), 
  geom = "point",                   
  pointsize = 3,                    
  pointshape = 19,                  
  habillage = data_mufa$qualysup,      
  palette = colores, 
  addEllipses = TRUE, 
  ellipse.type = "confidence", 
  ellipse.level = 0.95,
  ellipse.alpha = 0.15,             
  mean.point = FALSE,               
  title = "MUFA Diet by Subject", 
  legend.title = "Subject"
) +
  
  theme_classic(base_size = 14) +   
  theme(
    legend.position = "top",        
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1) 
  )


print(MUFA_subject)

ggsave("./MUFA_subject.svg", plot = MUFA_subject, width = 8, height = 6)

# NF

data_nf <- imputed_data %>%
  select(contains("NF"))

data_nf = t(data_nf) 

head(data_nf)
dim(data_nf)

data_nf = as.data.frame(data_nf)


data_nf$qualysup <- as.factor(str_extract(rownames(data_nf), "B$|P$|PP$"))
table(data_nf$qualysup)


pca.nf = PCA(data_nf, scale.unit = T, graph = F, quali.sup = ncol(data_nf))

NF_timepoint <- fviz_pca_ind(
  pca.nf, 
  axes = c(1, 2), 
  geom = "point",                   
  pointsize = 3,                    
  pointshape = 19,                  
  habillage = ncol(data_nf),      
  palette = c("#0073C2", "#EFC000", "pink"), 
  addEllipses = TRUE, 
  ellipse.type = "confidence", 
  ellipse.level = 0.95,
  ellipse.alpha = 0.15,             
  mean.point = FALSE,               
  title = "NF Diet by Time Point", 
  legend.title = "Time Point"
) +
  
  theme_classic(base_size = 14) +   
  theme(
    legend.position = "top",        
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1) 
  )


print(NF_timepoint)

ggsave("./NF_timepoint.svg", plot = NF_timepoint, width = 8, height = 6)


# NF Group subject replicates

data_nf$qualysup <- str_extract(rownames(data_nf), "^\\d+")
data_nf$qualysup <- as.factor(data_nf$qualysup) 
pca.nf = PCA(data_nf, scale.unit = T, graph = F, quali.sup = ncol(data_nf)) 
table(data_nf$qualysup)


nf_subject <- fviz_pca_ind(
  pca.nf, 
  axes = c(1, 2), 
  geom = "point",                   
  pointsize = 3,                    
  pointshape = 19,                  
  habillage = data_nf$qualysup,      
  palette = colores, 
  addEllipses = TRUE, 
  ellipse.type = "confidence", 
  ellipse.level = 0.95,
  ellipse.alpha = 0.15,             
  mean.point = FALSE,               
  title = "NF Diet by Subject", 
  legend.title = "Subject"
) +
  
  theme_classic(base_size = 14) +   
  theme(
    legend.position = "top",        
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1) 
  )


print(nf_subject)

ggsave("./NF_subject.svg", plot = nf_subject, width = 8, height = 6)


# SFA

data_sfa <- imputed_data %>%
  select(contains("SFA"))

data_sfa = t(data_sfa) 

head(data_sfa)
dim(data_sfa)

data_sfa = as.data.frame(data_sfa)


data_sfa$qualysup <- str_extract(rownames(data_sfa), "B$|P$|PP$")
table(data_sfa$qualysup)


pca.sfa = PCA(data_sfa, scale.unit = T, graph = F, quali.sup = ncol(data_sfa))

SFA_timepoint <- fviz_pca_ind(
  pca.sfa, 
  axes = c(1, 2), 
  geom = "point",                   
  pointsize = 3,                    
  pointshape = 19,                  
  habillage = ncol(data_sfa),      
  palette = c("#0073C2", "#EFC000", "pink"), 
  addEllipses = TRUE, 
  ellipse.type = "confidence", 
  ellipse.level = 0.95,
  ellipse.alpha = 0.15,             
  mean.point = FALSE,               
  title = "SFA Diet by Time Point", 
  legend.title = "Time Point"
) +
  
  theme_classic(base_size = 14) +   
  theme(
    legend.position = "top",        
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1) 
  )


print(SFA_timepoint)
ggsave("./SFA_timepoint.svg", plot = SFA_timepoint, width = 8, height = 6)


## SFA Group subject replicates

data_sfa$qualysup <- str_extract(rownames(data_sfa), "^\\d+")
data_sfa$qualysup <- as.factor(data_sfa$qualysup)
pca.sfa = PCA(data_sfa, scale.unit = T, graph = F, quali.sup = ncol(data_sfa)) 
table(data_sfa$qualysup)

sfa_subject <- fviz_pca_ind(
  pca.sfa, 
  axes = c(1, 2), 
  geom = "point",                   
  pointsize = 3,                    
  pointshape = 19,                  
  habillage = data_sfa$qualysup,      
  palette = colores, 
  addEllipses = TRUE, 
  ellipse.type = "confidence", 
  ellipse.level = 0.95,
  ellipse.alpha = 0.15,             
  mean.point = FALSE,               
  title = "SFA Diet by Subject", 
  legend.title = "Subject"
) +
  
  theme_classic(base_size = 14) +   
  theme(
    legend.position = "top",        
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1) 
  )


print(sfa_subject)

ggsave("./SFA_subject.svg", plot = sfa_subject, width = 8, height = 6)


# W3

data_w3 <- clean_data %>%
  select(contains("W3"))

data_w3 = t(data_w3) 

head(data_w3)
dim(data_w3)

data_w3 = as.data.frame(data_w3)


data_w3$qualysup <- str_extract(rownames(data_w3), "B$|P$|PP$")
table(data_w3$qualysup)


pca.w3 = PCA(data_w3, scale.unit = T, graph = F, quali.sup = ncol(data_w3))

w3_timepoint <- fviz_pca_ind(
  pca.w3, 
  axes = c(1, 2), 
  geom = "point",                   
  pointsize = 3,                    
  pointshape = 19,                  
  habillage = ncol(data_w3),      
  palette = c("#0073C2", "#EFC000", "pink"), 
  addEllipses = TRUE, 
  ellipse.type = "confidence", 
  ellipse.level = 0.95,
  ellipse.alpha = 0.15,             
  mean.point = FALSE,               
  title = "ω3 Diet by Time Point", 
  legend.title = "Time Point"
) +
  
  theme_classic(base_size = 14) +   
  theme(
    legend.position = "top",        
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1) 
  )


print(w3_timepoint)
ggsave("./w3_timepoint.svg", plot = w3_timepoint, width = 8, height = 6)

## w3 group subject replicates

data_w3$qualysup <- str_extract(rownames(data_w3), "^\\d+")
data_w3$qualysup <- as.factor(data_w3$qualysup)
pca.w3 = PCA(data_w3, scale.unit = T, graph = F, quali.sup = ncol(data_w3)) 
table(data_w3$qualysup)

w3_subject <- fviz_pca_ind(
  pca.w3, 
  axes = c(1, 2), 
  geom = "point",                   
  pointsize = 3,                    
  pointshape = 19,                  
  habillage = data_w3$qualysup,      
  palette = colores, 
  addEllipses = TRUE, 
  ellipse.type = "confidence", 
  ellipse.level = 0.95,
  ellipse.alpha = 0.15,             
  mean.point = FALSE,               
  title = "ω3 Diet by Subject", 
  legend.title = "Subject"
) +
  
  theme_classic(base_size = 14) +   
  theme(
    legend.position = "top",        
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1) 
  )


print(w3_subject)

ggsave("./w3_subject.svg", plot = w3_subject, width = 8, height = 6)


# P Time Point 

data_P <- imputed_data %>%
  select(ends_with("P"), -ends_with("PP"))

data_P = t(data_P) 

head(data_P)
dim(data_P)

data_P = as.data.frame(data_P)


data_P$qualysup <- str_extract(rownames(data_P), "MUFA|NF|SFA|W3")
table(data_P$qualysup)


pca.P = PCA(data_P, scale.unit = T, graph = F, quali.sup = ncol(data_P))


P_diet <- fviz_pca_ind(
  pca.P, 
  axes = c(1, 2), 
  geom = "point",                   
  pointsize = 3,                    
  pointshape = 19,                  
  habillage = ncol(data_P),      
  palette = c("#0073C2", "#EFC000", "pink","lightgreen"), 
  addEllipses = TRUE, 
  ellipse.type = "confidence", 
  ellipse.level = 0.95,
  ellipse.alpha = 0.15,             
  mean.point = FALSE,               
  title = "Dietary interventions in postprandial peak", 
  legend.title = "Time Point"
) +
  
  theme_classic(base_size = 14) +   
  theme(
    legend.position = "top",        
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1) 
  )


print(P_diet)
ggsave("./P_diet.svg", plot = P_diet, width = 8, height = 6)

## P time point group subject replicates

data_P$qualysup <- str_extract(rownames(data_P), "^\\d+")
data_P$qualysup <- as.factor(data_P$qualysup)
pca.P = PCA(data_P, scale.unit = T, graph = F, quali.sup = ncol(data_P)) 
table(data_P$qualysup)

P_subject <- fviz_pca_ind(
  pca.P, 
  axes = c(1, 2), 
  geom = "point",                   
  pointsize = 3,                    
  pointshape = 19,                  
  habillage = data_P$qualysup,      
  palette = colores, 
  addEllipses = TRUE, 
  ellipse.type = "confidence", 
  ellipse.level = 0.95,
  ellipse.alpha = 0.15,             
  mean.point = FALSE,               
  title = "Postprandial peak by subject", 
  legend.title = "Subject"
) +
  
  theme_classic(base_size = 14) +   
  theme(
    legend.position = "top",        
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1) 
  )


print(P_subject)

ggsave("./P_subject.svg", plot = P_subject, width = 8, height = 6)




# PP time point

data_PP <- imputed_data %>%
  select(contains("PP"))

data_PP = t(data_PP) 

head(data_PP)
dim(data_PP)

data_PP = as.data.frame(data_PP)


data_PP$qualysup <- str_extract(rownames(data_PP), "MUFA|NF|SFA|W3")
table(data_PP$qualysup)



pca.PP = PCA(data_PP, scale.unit = T, graph = F, quali.sup = ncol(data_PP))

PP_diet <- fviz_pca_ind(
  pca.PP, 
  axes = c(1, 2), 
  geom = "point",                   
  pointsize = 3,                    
  pointshape = 19,                  
  habillage = ncol(data_PP),      
  palette = c("#0073C2", "#EFC000", "pink","lightgreen"), 
  addEllipses = TRUE, 
  ellipse.type = "confidence", 
  ellipse.level = 0.95,
  ellipse.alpha = 0.15,             
  mean.point = FALSE,               
  title = "Dietary interventions in late postprandial phase", 
  legend.title = "Time Point"
) +
  
  theme_classic(base_size = 14) +   
  theme(
    legend.position = "top",        
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1) 
  )


print(PP_diet)
ggsave("./PP_diet.svg", plot = PP_diet, width = 8, height = 6)

## PP time point group subject replicates

data_PP$qualysup <- str_extract(rownames(data_PP), "^\\d+")
data_PP$qualysup <- as.factor(data_PP$qualysup)
pca.PP = PCA(data_PP, scale.unit = T, graph = F, quali.sup = ncol(data_PP)) 
table(data_PP$qualysup)

PP_subject <- fviz_pca_ind(
  pca.PP, 
  axes = c(1, 2), 
  geom = "point",                   
  pointsize = 3,                    
  pointshape = 19,                  
  habillage = data_PP$qualysup,      
  palette = colores, 
  addEllipses = TRUE, 
  ellipse.type = "confidence", 
  ellipse.level = 0.95,
  ellipse.alpha = 0.15,             
  mean.point = FALSE,               
  title = "Late postprandial phase by subject", 
  legend.title = "Subject"
) +
  
  theme_classic(base_size = 14) +   
  theme(
    legend.position = "top",        
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1) 
  )


print(PP_subject)

ggsave("./PP_subject.svg", plot = PP_subject, width = 8, height = 6)



#############################################################
#############################################################


# Correlation between baseline samples. 

baseline = seq(from = 1, to = 120, by = 3)
data_baseline = clean_data[,baseline]

# Selecting Baseline samples per condition
data_baseline_nf <- data_baseline %>%
  select(contains("NF"))
data_baseline_mufa <- data_baseline %>%
  select(contains("MUFA"))
data_baseline_sfa <- data_baseline %>%
  select(contains("SFA"))
data_baseline_w3 <- data_baseline %>%
  select(contains("W3"))
 
## Correrlation NF - MUFA

cor_nf_mufa <- cor(x = data_baseline_nf, 
                         y = data_baseline_mufa, 
                         use = "pairwise.complete.obs", #if use has the value "pairwise.complete.obs" 
                   #then the correlation or covariance between each pair of variables is computed using 
                   #all complete pairs of observations on those variables
                         method = "spearman") #Spearman correlation

color =  colorRampPalette(c("darkblue", "white", "darkred"))(400)

svglite("cor_nf_mufa.svg", width = 9, height = 9)

corr <- corrplot(cor_nf_mufa, 
                 method = "color", 
                 type = "full",       
                 addCoef.col = "black",  
                 number.cex = 1.5,      
                 tl.col = "black",    
                 tl.cex = 1.5,        
                 na.label = "NA",
                 col.lim = c(0, 1),      # No negative correlations were found
                 col = color,
                 is.corr = FALSE)

dev.off()

## Correrlation NF - SFA

cor_nf_sfa <- cor(x = data_baseline_nf, 
                   y = data_baseline_sfa, 
                   use = "pairwise.complete.obs",
                  method = "spearman")

svglite("cor_nf_sfa.svg", width = 9, height = 9)

corr <- corrplot(cor_nf_sfa, 
                 method = "color", 
                 type = "full",       
                 addCoef.col = "black",  
                 number.cex = 1.5,       
                 tl.col = "black",    
                 tl.cex = 1.5,        
                 na.label = "NA",
                 col.lim = c(0, 1),      
                 col = color,
                 is.corr = FALSE)

dev.off()

## Correrlation NF - w3

cor_nf_w3 <- cor(x = data_baseline_nf, 
                  y = data_baseline_w3, 
                  use = "pairwise.complete.obs",
                 method = "spearman")

svglite("cor_nf_w3.svg", width = 9, height = 9)

corr <- corrplot(cor_nf_w3, 
                 method = "color", 
                 type = "full",       
                 addCoef.col = "black",  
                 number.cex = 1.5,       
                 tl.col = "black",    
                 tl.cex = 1.5,        
                 na.label = "NA",
                 col.lim = c(0, 1),      
                 col = color,
                 is.corr = FALSE)

dev.off()


#########################################################

## Acquiring the names of the total protein groups -> Dataset functional enrichment

PG_names_clean <- unlist(strsplit(PG_names, ";"))
PG_names_clean <- unique(PG_names_clean)
PG_names_clean <- unique(PG_names_clean)

writeLines(PG_names_clean, "./data_generados/total_names_proteins.txt")

## EnrichR - Jensen Compartments

# Total protein list
ids <- readLines("total_names_proteins.txt")

# Formatting our protein group names
ids_clean <- ids[!grepl("^ENSEMBL:", ids)]       
ids_clean <- sub("-\\d+$", "", ids_clean)          
ids_clean <- trimws(ids_clean)
ids_clean <- unique(ids_clean[ids_clean != ""])

# Convert UniProt accession → Gene Symbol
up <- UniProt.ws(taxId = 9606)  # 9606 = Homo sapiens

mapping <- AnnotationDbi::select(up,
                  keys = ids_clean,
                  columns = c("gene_primary"),
                  keytype = "UniProtKB")

gene_symbols <- unique(na.omit(mapping$))

#  enrichR woth gene symbols
setEnrichrSite("Enrichr")
dbs <- c("Jensen_COMPARTMENTS")  

results <- enrichr(gene_symbols, dbs)

## Results
colores_low = c(rep("beige",5))
colores_high = c(rep("darkgrey",5))

df_cluster <- results[["Jensen_COMPARTMENTS"]] %>%
  
  # Significative filtering
  filter(Adjusted.P.value < 0.05) %>%
  
  # Calculating -log10(Adjusted.P.value)
  mutate(LogQ_pos = -log10(Adjusted.P.value)) %>%
  
  # Term name
  mutate(Term = str_wrap(Term, width = 40)) %>%
  
  # Top 10
  arrange(desc(LogQ_pos)) %>%
  head(10)

# Plotting

p <- ggplot(df_cluster, aes(x = reorder(Term, LogQ_pos), y = LogQ_pos, fill = LogQ_pos)) +
  geom_col(color = "black", width = 0.7, linewidth = 0.3) +
  coord_flip() +
  scale_fill_gradient(low = colores_low[1], high = colores_high[1],
                      name = "-log10(adj.P)", limits = c(0, NA)) +
  labs(
    x = "",
    y = "-log10(Adjusted P-value)",
    title = "Jensen Compartments - Functional Enrichment"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(color = "black", size = 9, lineheight = 0.7),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(p)
ggsave(filename = "jensen_compartments.svg", plot = p, device = "svg", height = 5, width = 7, path = "./")    


#########################################################

## Z-Score - Heatmap

##Reorganize data frame 

ordered_data <- imputed_data %>% #Ordering by dietary intervention
  select(
    contains("NF"),
    contains("SFA"),
    contains("MUFA"),
    contains("W3"),
    everything()  #Ensuring that if there are any extra columns 
    #(such as variable names), they are placed at the end and not lost.
  )

ordered_data_2 <- ordered_data %>% #Ordering by time point
  select(
    contains("B"),
    ends_with("P"), -ends_with("PP"),
    contains("PP"),
    everything()  
  ) 

ordered_data_2$X <- NULL #Eliminate column with protein names
rownames(ordered_data_2) = rownames(clean_data) #protein names as row names
dim(ordered_data_2)


## Z - Score calculation
matrix_scaled <- t(scale(t(matrix_heatmap)))
matrix_scaled[matrix_scaled >  3] <-  3
matrix_scaled[matrix_scaled < -3] <- -3

col_fun <- colorRamp2(
  c(-3, -1.5, 0, 1.5, 3),
  c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B")
)


## Dietary group
grupos <- dplyr::case_when(
  grepl("NF",   colnames(matrix_scaled)) ~ "NF",
  grepl("SFA",  colnames(matrix_scaled)) ~ "SFA",
  grepl("MUFA", colnames(matrix_scaled)) ~ "MUFA",
  grepl("W3",   colnames(matrix_scaled)) ~ "W3",
  TRUE ~ "Other"
)

group_colors <- c(
  "NF"   = "darkgrey",
  "SFA"  = "#E41A1C",
  "MUFA" = "olivedrab",
  "W3"   = "#377EB8"
)

## Time point
time <- dplyr::case_when(
  grepl("PP", colnames(matrix_scaled)) ~ "PP", # PP before P to avoid missmatches
  grepl("B",  colnames(matrix_scaled)) ~ "B",
  grepl("P",  colnames(matrix_scaled)) ~ "P",
  TRUE ~ "Other"
)

time_colors <- c(
  "B"  = "#8073AC", 
  "P"  = "#E08214",   
  "PP" = "pink"    
)

# Combined annotation
col_annotation <- HeatmapAnnotation(
  
  "Diet"   = grupos,
  "Time" = time,
  
  col = list(
    "Diet"   = group_colors,
    "Time" = time_colors
  ),
  
  annotation_height  = unit(c(4, 4), "mm"),   # height of the bars
  annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
  gap                = unit(1, "mm"),           # space between the  bars
  
  annotation_legend_param = list(
    "Diet" = list(
      title       = "Diet group",
      title_gp    = gpar(fontsize = 9, fontface = "bold"),
      labels_gp   = gpar(fontsize = 8),
      grid_height = unit(4, "mm"),
      grid_width  = unit(4, "mm")
    ),
    "Time" = list(
      title       = "Time point",
      title_gp    = gpar(fontsize = 9, fontface = "bold"),
      labels_gp   = gpar(fontsize = 8),
      grid_height = unit(4, "mm"),
      grid_width  = unit(4, "mm")
    )
  )
)

## Heatmap
ht <- Heatmap(
  matrix_scaled,
  col  = col_fun,
  name = NULL,
  
  cluster_columns          = FALSE,
  cluster_rows             = TRUE,
  clustering_distance_rows = "pearson",
  clustering_method_rows   = "ward.D2",
  
  show_row_dend    = TRUE,
  show_column_dend = FALSE,
  
  top_annotation = col_annotation,
  
  show_row_names    = FALSE,
  show_column_names = TRUE,
  column_names_gp   = gpar(fontsize = 2.5, fontface = "plain"),
  column_names_rot  = 45,
  column_names_side = "bottom",
  
  rect_gp = gpar(col = NA),
  
  heatmap_legend_param = list(
    title         = "Z-score",
    title_gp      = gpar(fontsize = 9, fontface = "bold"),
    labels_gp     = gpar(fontsize = 8),
    legend_height = unit(35, "mm"),
    at            = c(-3, -1.5, 0, 1.5, 3),
    labels        = c("−3", "−1.5", "0", "1.5", "3"),
    border        = "black"
  ),
  
  row_title    = NULL,
  column_title = NULL,
  width        = unit(14, "cm"),
  height       = unit(12, "cm")
)

## Export to SVG format
svglite("zscore_heatmap.svg", width = 18, height = 14)

draw(ht,
     heatmap_legend_side    = "right",
     annotation_legend_side = "right",
     padding = unit(c(5, 5, 5, 5), "mm")
)

dev.off()

