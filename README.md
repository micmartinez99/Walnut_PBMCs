# Walnut_PBMCs

This analysis was conducted on samples from a clinical study (clinical trial ID: NCT05195970).
Individuals at high risk for colon cancer were subjected to a week 1 wash out for all ellagic acid-containing foods, followed by 3 weeks of eating walnuts. Whole blood from these patients at visit 2 (directly after wash out) and visit 3 (after walnut supplementation) was collected. Whole blood samples were processed to collect PBMCs.

PBMCs were sent to Eve Technologies for a 21-plex inflammatory marker panel. This respository holds the R code used to analyze the data.

# Analysis Methods
The delta (visit 3 - visit 2) was taken for each marker. TNF-a was removed from the analysis completely do to most samples being over-saturated. From these delta values, principal component analysis was conducted. 

Raw data was then log transformed, the delta of the log values was taken and the average marker expression level was taken across obese and non obese patients to generate a heatmap (scaled by column.)

Normality of data was assessed using Shapiro-Wilk test, demonstrating most markers followed a non-normal distribution. As a result, delta values for each patient, stratified by obese vs non-obese were compared using Wilcox Test. 

Correlations between delta marker levels and BMI were assessed using the non-parametric Spearman test.

# Intracellular Staining
PBMCs (n of 2) were obtained from ATCC and stained for intracellular targets (IL-1b, IL-6, TNF-a, FoxP3), as well as surface markers for B and T cells and Monocytes.
Gating was performed in FlowJo by Dr. Evan Jellison (UCHC Flow Cytometry Core Director). An excel file with cell popualtion percentages was provided. Code 002 is the code to process this data and export smaller csvs for use in GraphPad.
Statistical analysis was conducted using ANOVA and Tukey HSD. 
Monocytes were dropped from analysis due to low sample presence (as per the advise of Dr. Jellison.)



