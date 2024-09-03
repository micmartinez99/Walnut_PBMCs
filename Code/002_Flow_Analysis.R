# The purpose of this code is to analyze the intracellular staining data for Yasheen Gao and Nuoxi Fan. 

# Load libraries
library(tidyverse)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)

# Initialize a new function to generate output directories
newDir <- function(x) {
  if (!dir.exists(x)) {
    dir.create(x)
  }
}

# Generate an output directory
opDir <- newDir("Outputs/002_Flow_Analysis")

# Generate an output directory for files to input to prism
groupFiles <- newDir("Outputs/002_Flow_Analysis/Processed_Data")

################################################################################

# Read in the flow cytometry data
flow <- read.csv("Data/Flow_Cytometry_Data/Walnut_Flow_Cytometry.csv")


# Take the rows that indicate how many live cells we have (all statistics are
# normalized to the number of live cells already)
liveCells <- flow[flow$Depth == "> > > ",]

# Create a column for treatment
flow <- flow %>%
  mutate(treatment = str_extract(Name, "^[^\\.]+"),
         Name = str_replace(Name, "^.*?\\.fcs", ""))

# Create replicate information
flow$Exp <- ifelse(flow$treatment == "C1 Stim 1" | flow$treatment == "C3 Stim 3", "Stimulated",
                   ifelse(flow$treatment == "C2 Unstim 1" | flow$treatment == "C4 unstim 2", "Unstimulated",
                          ifelse(flow$treatment == "C5 uro A 5uM 1" | flow$treatment == "C6 uro A 5uM 2", "Uro A 5uM",
                                 ifelse(flow$treatment == "C7 uro A 10uM -1" | flow$treatment == "C8 uro A 10uM - 2", "Uro A 10uM",
                                        ifelse(flow$treatment == "D1 DMSO 1" | flow$treatment == "D2 DMSO 2", "DMSO Cntl",
                                               ifelse(flow$treatment == "D3 DMSO 1 no IC" | flow$treatment == "D4 DMSO 2 no IC", "DMSO IC-",
                                                      ifelse(flow$treatment == "D5 uro A 10uM -2 no IC" | flow$treatment == "E4 uroA 10uM 1 - no IC", "Uro A 10uM IC-",
                                                             ifelse(flow$treatment == "D6 uro A 5 uM - 2 no IC" | flow$treatment == "D7 uro A 5 uM - 1 no IC", "Uro A 5uM IC-",
                                                                    ifelse(flow$treatment == "D8 Stim 1 no IC" | flow$treatment == "E2 Stim 3 no IC", "Stim IC-", "Unstim IC-")))))))))




################################################################################

#----- Loop to process the data

# Remove NA
flow <- na.omit(flow)

# Only take experimental groups
keep <- c("Stimulated", "Unstimulated", "Uro A 5uM", "Uro A 10uM")
flow <- flow[flow$Exp %in% keep,]

# We are not interested in the scatter or single cell columns
remove <- c("Scatter", "Single_Cells")
flow <- flow[!flow$Group %in% remove,]

# Get a vector of unique gated groups
gates <- unique(flow$Group)

for (i in gates) {
  
  # Subset the gate
  subset <- flow[flow$Group == i,]
  
  print(subset)

  # Order
  subset <- subset[order(subset$Exp, decreasing = TRUE),]
  
  # Add a column for replicate
  subset$replicate <- rep(1:2, length.out = nrow(subset))
  
  # Keep only these rows
  keep <- c("Group", "Exp", "Statistic", "replicate")
  subsetNew <- subset[,colnames(subset) %in% keep]
  
  # Calculate the mean and SEM
  summaryData <- subsetNew %>%
    group_by(Exp) %>%
    summarize(
      Mean = mean(Statistic),
      SEM = sd(Statistic) / sqrt(n())
    )
  
  # Pivot wider
  subsetNew <- subsetNew %>%
    pivot_wider(names_from = Exp, values_from = Statistic)
  
  subsetNew <- subsetNew[,c(1,2,5,6,3,4)]

  fileName <- paste(i, "Processed.csv", sep = "_")
  
  # Set colots
  colors <- c("Stimulated" = "firebrick",
              "Unstimulated" = "grey",
              "Uro A 5uM" = "goldenrod",
              "Uro A 10uM" = "cornflowerblue")
  
  # Factor
  summaryData$Exp <- factor(summaryData$Exp, levels = c("Unstimulated", "Stimulated", "Uro A 5uM", "Uro A 10uM"))
  
  
  # Plot
  plot <- ggplot(summaryData, aes(x = Exp, y = Mean, fill = Exp)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.2) +
    labs(x = "",
         y = "Mean Percent Cells", 
         title = i,
         fill = "") +
    scale_fill_manual(values = colors) +
    theme_classic() +
    theme(axis.title.y = element_text(face = "bold", size = 20),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(face = "bold", size = 14),
          legend.position = "none") +
    stat_pwc(method = "tukey.hsd")
  ggsave(paste("Outputs/002_Flow_Analysis/Processed_Data", paste(i, ".png"), sep = "/"), plot, width = 8, height = 8)
  
  

  
  # Write as a csv
  write.csv(subsetNew, file = paste("Outputs/002_Flow_Analysis/Processed_Data", paste(i, ".csv"), sep = "/"))
}






