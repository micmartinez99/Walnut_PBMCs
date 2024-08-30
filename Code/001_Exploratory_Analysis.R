# The purpose of this code is to conduct an exploratory analysis for the Eve Technologies PBMC data collected by
# Yasheen Gao and Nuoxi Fan (Summer 2024)

# Load libraries
library(tidyverse)
library(dplyr)
library(ggpllot2)
library(ggrepel)
library(RColorBrewer)
library(lme4)
library(lmerTest)
library(factoextra)

# Initialize a function to generate new directories
newDir <- function(x) {
  if (!dir.exists(x)) {
    dir.create(x)
  }
}

# Generate an output directory for this script
opDir <- newDir("Outputs/001_Exploratory_Analysis")

################################################################################

# Read in the data
eve <- read.csv("Data/Input_Data/Formatted_Eve_Tech_Data.tsv", sep = "\t") %>% 
  column_to_rownames(var = "Sample_Name_R")


# ----- Calculation of the delta for each marker

# Subset just the V2 rows
V2 <- eve[eve$Visit == "V2",]
V2meta <- V2[,1:10]
V2[,1:10] <- NULL

# Subset just the V3 rows
V3 <- eve[eve$Visit == "V3",]
V3meta <- V3[,1:10]
V3[,1:10] <- NULL

# Calculate delta and re-append metadata
delta <- V3 - V2
delta <- cbind(V2meta, delta)

# Change "V2" to "Delta"
delta$Visit <- "Delta"

# Factor levels
delta$Obesity <- factor(delta$Obesity, levels = c("NonObese", "Obese"))
delta$Sex <- factor(delta$Sex, levels = c("M", "F"))

#----- Principal Component Analysis

# Run PCA
pca_result <- prcomp(delta[,11:30], center = TRUE, scale. = TRUE)

# View the PCA results
PCA_summary <- summary(pca_result) 

# Get the % variation explained by each PC
proportion_variance <- PCA_summary$importance[2, ]
percent_var <- proportion_variance * 100

# Convert the PCA results to a dataframe and append metadata
pca_df <- as.data.frame(pca_result$x)
pca_df <- cbind(V2meta, pca_df)
pca_df$Visit <- "Delta"

# Get PCA loadings
loadings <- pca_result$rotation

# Plot PCA plot
PCA_analysis <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Obesity, label = Sample)) +
  geom_point() +
  geom_text_repel() +
  stat_ellipse(aes(color = Obesity, fill = Obesity), level = 0.95) + 
  labs(x = "PC1: 36.1%",
       y = "PC2: 20.27%",
       color = "") +
  scale_color_manual(values = c("NonObese" = "cornflowerblue",
                                "Obese" = "firebrick")) +
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        legend.position = "bottom",
        legend.text = element_text(face = "bold"))
ggsave("Outputs/001_Exploratory_Analysis/PCA_Analysis.png", PCA_analysis, width = 6, height = 6)

pca_df$marker <- rownames(pca_df)

colors <- c("NonObese" = "cornflowerblue",
            "Obese" = "firebrick")
pca_df$Colors <- ifelse(pca_df$Obesity == "Obese", "firebrick", "cornflowerblue")

################################################################################

#----- Heatmap
# Extract the meta
eveMeta <- eve[,1:10]
eve[,1:10] <- NULL
eve$TNF.a <- NULL

# Add a small constant to each value to prepare to log transform
eveLog <- eve + 0.01

# Log1p transform
eveLog <- log1p(eveLog)

# Re-merge the metadata
eveLog <- cbind(eveMeta, eveLog)

# Separate by V2 and V3
v2 <- eveLog[eveLog$Visit == "V2",]
v2Meta <- v2[,1:10]
v3 <- eveLog[eveLog$Visit == "V3",]

# Take the delta and re-append metadata
deltaLog <- v3[,11:30] - v2[,11:30]
deltaLog <- cbind(v2Meta, deltaLog)

# Set annotation dataframe
anno <- as.data.frame(deltaLog$Obesity)
rownames(anno) <- deltaLog$Sample
rownames(deltaLog) <- deltaLog$Sample

# Split by obese
obese <- deltaLog[deltaLog$Obesity == "Obese",]
obese[,1:10] <- NULL

# Take colmeans for obese patients
obeseMeans <- as.data.frame(colMeans(obese))
colnames(obeseMeans) <- "Obese"

# Take colmeans for nonobese patients
no <- deltaLog[deltaLog$Obesity == "NonObese",]
no[,1:10] <- NULL
noMeans <- as.data.frame(colMeans(nonObese))
colnames(noMeans) <- "NonObese"

# Merge the means
means <- cbind(noMeans, obeseMeans)

# Plot heatmap
HM <- pheatmap(means, 
               scale = "column",
               cluster_cols = FALSE,
               fontsize_row = 14,
               fontsize_col = 14,
               angle_col = 0,
               color = turbo(100))
ggsave("Outputs/001_Exploratory_Analysis/Individual_Effects/Log1p_Delta_HM.png", HM)

################################################################################
#---- Delta T-test

# Get a vector of markers
markers <- colnames(delta[,11:ncol(delta)])

# Tests for normality
o <- delta[delta$Obesity == "Obese",]
no <- delta[delta$Obesity == "NonObese",]

for (i in markers) {
  res <- shapiro.test(delta[,i])
  print(i)
  print(res)
  print("##################################")
}

# Iterate through the markers and test
for (i in markers) {
  
  # Set a file name
  fileName <- paste(i, "Delta_Obesity_Boxplot.png", sep = "_")
  
  # For obesity
  obesity <- ggplot(delta, aes_string(x = "Obesity", y = i, fill = "Obesity")) +
    geom_boxplot(outliers = FALSE, alpha = 0.5, width = 0.4) +
    geom_point() +
    stat_compare_means(method = "wilcox.test", hjust = -0.5, size = 8) +
    labs(y = paste("Delta", i, "Levels (pg/ml)", sep = " "),
         x = "",
         fill = "") +
    scale_fill_manual(values = c("NonObese" = "cornflowerblue",
                                 "Obese" = "firebrick")) +
    theme_classic() +
    theme(axis.title.y = element_text(face = "bold", size = 20),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(face = "bold", size = 20),
          legend.position = "none")
  
  # ggsave
  ggsave(paste("Outputs/001_Exploratory_Analysis/Individual_Effects/Delta_TTest_Obesity", fileName, sep = "/"), obesity, width = 8, height = 8)
}


#----- Individual Paired plots

for (i in markers) {
  
  # Set filename
  fileName <- paste(i, "Individual_PairedPlot.png", sep = "_")
  
  # Plot
  pairedPlot <- ggplot(eve, aes_string(x = "Visit", y = i, color = "Visit", label = "Sample")) +
      geom_text_repel(size = 3, color = "black") +
      geom_point(size = 5, alpha = 0.7) +
      scale_color_manual(values = c("V2" = "grey75", 
                                  "V3" = "grey25")) +
      geom_line(aes(group = Sample), color = "darkgrey", linetype = "dashed") +  
      geom_vline(xintercept = 30, linetype = "dashed") +
      labs(x = "Visit",
           y = paste(i, "pg/ml", sep = " "),
           color = "") +
      facet_grid(~Obesity, scales = "free") +
      theme_classic() +
      theme(axis.title.x = element_text(face = "bold", size = 14),
            axis.title.y = element_text(face = "bold", size = 14),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            strip.text = element_text(face = "bold", size = 14),
            legend.position = "none")
  
  # Save
  ggsave(paste("Outputs/001_Exploratory_Analysis/Individual_Effects/Subject_Delta", fileName, sep = "/"), pairedPlot, width = 8, height = 8)
}


#----- Dunn Test

# Set a grouping variable
eve$Group <- paste(eve$Visit, eve$Obesity, sep = "-")
eve$Group <- factor(eve$Group, levels = c("V2-NonObese", "V3-NonObese", "V2-Obese", "V3-Obese"))

# Iterate through
for (i in markers) {
  
  # Set filename
  fileName <- paste(i, "Delta_LMPlot.png", sep = "_")
  
  # Plot
  delta_lm <- ggplot(eve, aes_string(x = "Group", y = i, fill = "Group", label = "Sample")) +
    geom_boxplot() +
    geom_point() +
    stat_pwc(method = "dunn_test") +
    geom_text_repel(size = 3, color = "black") +
    scale_fill_manual(values = c("V2-NonObese" = "cornflowerblue",
                                 "V3-NonObese" = "cornflowerblue",
                                 "V2-Obese" = "firebrick",
                                 "V3-Obese" = "firebrick"))

  # Save
  ggsave(paste("Outputs/001_Exploratory_Analysis/Individual_Effects/Dunn_Test", fileName, sep = "/"), delta_lm, width = 8, height = 8)
}

#-----Correlations

IL17cor <- ggplot(delta, aes(x = BMI, y = Il.17a)) +
  geom_point() +
  stat_cor(inherit.aes = FALSE, aes(x = BMI, y = Il.17a), method = "spearman", vjust = -1, size = 4) +
  geom_smooth(method = lm, se = TRUE) +
  geom_vline(xintercept = 30, linetype = "dashed", color = "darkgrey") +
  labs(y = "Delta IL-17a",
       x = "BMI") +
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold", size = 20),
        axis.title.y = element_text(face = "bold", size = 20),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none")
ggsave("Outputs/001_Exploratory_Analysis/IL17a_BMI_Correlation.png", IL17cor, width = 6, height = 6)

IL21cor <- ggplot(delta, aes(x = BMI, y = Il.21)) +
  geom_point() +
  stat_cor(inherit.aes = FALSE, aes(x = BMI, y = Il.21), method = "spearman", vjust = -1, size = 4) +
  geom_smooth(method = lm, se = TRUE) +
  geom_vline(xintercept = 30, linetype = "dashed", color = "darkgrey") +
  labs(y = "Delta IL-21",
       x = "BMI") +
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold", size = 20),
        axis.title.y = element_text(face = "bold", size = 20),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none")
ggsave("Outputs/001_Exploratory_Analysis/IL21_BMI_Correlation.png", IL21cor, width = 6, height = 6)


newDir("Outputs/001_Exploratory_Analysis/Individual_Effects/Correlations")
for (i in markers) {
  
  fileName <- paste(i, "Correlation.png", sep = "_")
  
  
  corPlot <- ggplot(delta, aes_string(x = "BMI", y = i)) +
    geom_point() +
    stat_cor(inherit.aes = FALSE, aes_string(x = "BMI", y = i), method = "spearman", size = 4) +
    geom_smooth(method = lm, se = TRUE) +
    geom_vline(xintercept = 30, linetype = "dashed", color = "darkgrey") +
    labs(y = paste("Delta", i, sep = " "),
         x = "BMI") +
    theme_classic() +
    theme(axis.title.x = element_text(face = "bold", size = 20),
          axis.title.y = element_text(face = "bold", size = 20),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.position = "none")
  ggsave(paste("Outputs/001_Exploratory_Analysis/Individual_Effects/Correlations", fileName, sep = "/"), corPlot, width = 6, height = 6)
}




         











  
  



    







