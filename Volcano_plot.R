
############ Volcano plot of the methanol water shared features ################

# Loading necessary packages
library(notame)
library(dplyr)
library(readxl)
library(scales)
library(gridExtra)
library(ggplot2)
library(cowplot)

######################### Positive polarity ####################################
# Loading the feature list
data_vol_pos <- read_from_excel(file = "Data/Feat_list_Volcano_MeOH_H2O_pos_2notame.xlsx",
                                sheet = 1, corner_row = 4, corner_column = "F",
                                split_by = c("Column", "Ion_Mode"))
# Create a MetaboSet  
modes_vol_pos <- construct_metabosets(exprs = data_vol_pos$exprs,
                                      pheno_data = data_vol_pos$pheno_data,
                                      feature_data = data_vol_pos$feature_data,
                                      group_col = "Group")
# Extract the positive ionization dataset
mode_vol_pos <- modes_vol_pos$RP_Positive
# Change the features with value equal to 0 to NA
mode_vol_pos <- mark_nas(mode_vol_pos, value = 0)
# Loading a specific Water/Methanol metabolites subset using the MolNetInvert tool
specific_pos <- read_excel("Data/Metadata_MN_specific_H2O_MeOH_metabolites.xlsx")
# Remove the Water/Methanol-specific features
mode_vol_pos <- mode_vol_pos[!mode_vol_pos@featureData@data$Mzmine_ID %in%
                               specific_pos$`Mzmine ID`,]
# Flag the features with a low detection rate
mode_vol_pos <- flag_detection(mode_vol_pos, qc_limit = 0, group_limit = 0.5)
# Features clustering
clustered_vol_pos <- cluster_features(mode_vol_pos, rt_window = 1/60,
                                      all_features = FALSE, corr_thresh = 0.90,
                                      d_thresh = 0.85)
compressed_vol_pos <- compress_clusters(clustered_vol_pos)
# Data imputation
set.seed(5646)
imputed_vol_pos <- impute_rf(compressed_vol_pos)
# Extract clean data
noflag_vol_pos <- drop_flagged(imputed_vol_pos)
# Drop QC
noqc_vol_pos <- drop_qcs(noflag_vol_pos)
# Performing homoscedasticity test
Levene_vol_pos <- perform_homoscedasticity_tests(noqc_vol_pos,
                                                 formula_char = "Feature ~ Group")
# Adding homoscedasticity results to notame MetaboSet
noqc_vol_pos <- join_fData(noqc_vol_pos, Levene_vol_pos)
# Extracting features with equal variance
vol_pos_equal_var <- noqc_vol_pos[noqc_vol_pos@featureData@data$Levene_P_FDR > 0.05,]
# Extracting features with unequal variance
vol_pos_unequal_var <- noqc_vol_pos[noqc_vol_pos@featureData@data$Levene_P_FDR < 0.05,]
# Fold change between Methanol vs. Water
vol_pos_fc <- fold_change(noqc_vol_pos, group = "Group")
# The two-sample t-test performing
vol_pos_tt <- perform_t_test(vol_pos_equal_var,
                             formula_char = "Feature ~ Group",
                             var.equal = TRUE)
# Adding a tag for t-test results
vol_pos_tt <- vol_pos_tt %>%
  mutate(Statistic_test = case_when(
    Feature_ID != 0 ~ "Student's t-test"))
# The Welch's t-test performing
vol_pos_wtt <- perform_t_test(vol_pos_unequal_var,
                              formula_char = "Feature ~ Group",
                              var.equal = FALSE)
# Adding a tag for Welch's t-test results
vol_pos_wtt <- vol_pos_wtt %>%
  mutate(Statistic_test = case_when(
    Feature_ID != 0 ~ "Welch's t-test"))
# Merge Welch's t-test and Student's t-test
vol_pos_pvalue <- rbind(vol_pos_tt, vol_pos_wtt)
# Adding the fold change to statistical results
vol_pos_data <- left_join(vol_pos_pvalue, vol_pos_fc)
# Log-transform for visualization
vol_pos_data$logP <- -log10(vol_pos_data$H2O_vs_MEOH_t_test_P_FDR)
vol_pos_data$log2_fc <- log2(vol_pos_data$MEOH_vs_H2O_FC)
# Determine point colors based on significance and fold change
vol_pos_data <- vol_pos_data %>%
  mutate(point_color = case_when(
    H2O_vs_MEOH_t_test_P_FDR < 0.05 & log2_fc < -1 ~ "Down", # significantly down
    H2O_vs_MEOH_t_test_P_FDR < 0.05 & log2_fc > 1 ~ "Up"))   # significantly up
vol_pos_data$point_color[is.na(vol_pos_data$point_color)] <- "Not sig"
# Extracting feature identified
metab_data_pos <- noqc_vol_pos[!is.na(noqc_vol_pos@featureData@data$Metabolite),]
# Extracting metabolite table
meta_table_pos <- metab_data_pos@featureData@data
# Creating a new small table of the annotated compounds
# Keep metabolites with p-value < 0.05
volc_compouds_pos <- subset(vol_pos_data, H2O_vs_MEOH_t_test_P_FDR < 0.05 &
                              point_color != "Not sig")
volc_compouds_pos <- left_join(meta_table_pos, volc_compouds_pos)
# Volcano plot
vol_pos_plot <- ggplot(vol_pos_data, aes(log2_fc, logP)) +
  scale_colour_manual(values = c("Not sig" = "#E5E5E5",
                                 "Down" = "#3366FF",
                                 "Up" = "#009933")) +
  theme_classic() +
  geom_point(data = vol_pos_data,
             aes(shape = Statistic_test,
                 color = point_color),
             size = 2.5,
             alpha = 0.4)  +
  guides(shape = "none") +
  labs(shape = 'Statistical test',
       color = 'Significance') +
  ggrepel::geom_label_repel(data = volc_compouds_pos,
                            aes(label = meta_table_pos$Metabolite),
                            color = "black",
                            box.padding = 0.37,
                            label.padding = 0.22,
                            label.r = 0.30,
                            cex = 2.5,
                            max.overlaps = 20,
                            min.segment.length = 0,
                            seed = 42) +
  xlab(bquote(Log[2](Methanol/Water))) + ylab(bquote(-Log[10](p-value))) +
  scale_y_continuous(limits = c(0,10), labels = label_scientific(),
                     expand=c(0.003, 0.003)) +
  theme(legend.position = c(0.85, 0.90),
        legend.background = element_rect(fill = "white", color = "black")) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill= "transparent")) +
  geom_vline(xintercept = 0, linetype = "longdash", colour="gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "longdash", colour="gray")
vol_pos_plot

######################### Negative polarity ####################################
# Loading the feature list
data_vol_neg <- read_from_excel(file = "Data/Feat_list_Volcano_MeOH_H2O_neg_2notame.xlsx",
                                sheet = 1, corner_row = 4, corner_column = "F",
                                split_by = c("Column", "Ion_Mode"))
# Create a MetaboSet  
modes_vol_neg <- construct_metabosets(exprs = data_vol_neg$exprs,
                                      pheno_data = data_vol_neg$pheno_data,
                                      feature_data = data_vol_neg$feature_data,
                                      group_col = "Group")
# Extract the positive ionization dataset
mode_vol_neg <- modes_vol_neg$RP_Negative
# Change the features with value equal to 0 to NA
mode_vol_neg <- mark_nas(mode_vol_neg, value = 0)
# Loading a specific Water/Methanol metabolites subset using the MolNetInvert tool
specific_neg <- read_excel("Data/Metadata_MN_specific_H2O_MeOH_metabolites_neg.xlsx")
# Remove the Water/Methanol-specific features
mode_vol_neg <- mode_vol_neg[!mode_vol_neg@featureData@data$Mzmine_ID %in%
                               specific_neg$`Mzmine ID`,]
# Flag the features with a low detection rate
mode_vol_neg <- flag_detection(mode_vol_neg, qc_limit = 0, group_limit = 0.5)
# Features clustering
clustered_vol_neg <- cluster_features(mode_vol_neg, rt_window = 1/60,
                                      all_features = FALSE, corr_thresh = 0.90,
                                      d_thresh = 0.85)
compressed_vol_neg <- compress_clusters(clustered_vol_neg)
# Data imputation
set.seed(1562)
imputed_vol_neg <- impute_rf(compressed_vol_neg)
# Extract clean data
noflag_vol_neg <- drop_flagged(imputed_vol_neg)
# Drop QC
noqc_vol_neg <- drop_qcs(noflag_vol_neg)
# Performing homoscedasticity test
Levene_vol_neg <- perform_homoscedasticity_tests(noqc_vol_neg,
                                                 formula_char = "Feature ~ Group")
# Adding homoscedasticity results to notame MetaboSet
noqc_vol_neg <- join_fData(noqc_vol_neg, Levene_vol_neg)
# Extracting features with equal variance
vol_neg_equal_var <- noqc_vol_neg[noqc_vol_neg@featureData@data$Levene_P_FDR > 0.05,]
# Extracting features with unequal variance
vol_neg_unequal_var <- noqc_vol_neg[noqc_vol_neg@featureData@data$Levene_P_FDR < 0.05,]
# Fold change between Methanol vs. Water
vol_neg_fc <- fold_change(noqc_vol_neg, group = "Group")
# The two-sample t-test performing
vol_neg_tt <- perform_t_test(vol_neg_equal_var,
                             formula_char = "Feature ~ Group",
                             var.equal = TRUE)
# Adding a tag for t-test results
vol_neg_tt <- vol_neg_tt %>%
  mutate(Statistic_test = case_when(
    Feature_ID != 0 ~ "Student's t-test"))
# The Welch's t-test performing
vol_neg_wtt <- perform_t_test(vol_neg_unequal_var,
                              formula_char = "Feature ~ Group",
                              var.equal = FALSE)
# Adding a tag for Welch's t-test results
vol_neg_wtt <- vol_neg_wtt %>%
  mutate(Statistic_test = case_when(
    Feature_ID != 0 ~ "Welch's t-test"))
# Merge Welch's t-test and Student's t-test
vol_neg_pvalue <- rbind(vol_neg_tt, vol_neg_wtt)
# Adding the fold change to statistical results
vol_neg_data <- left_join(vol_neg_pvalue, vol_neg_fc)
# Log-transform for visualization
vol_neg_data$logP <- -log10(vol_neg_data$H2O_vs_MEOH_t_test_P_FDR)
vol_neg_data$log2_fc <- log2(vol_neg_data$MEOH_vs_H2O_FC)
# Determine point colors based on significance and fold change
vol_neg_data <- vol_neg_data %>%
  mutate(point_color = case_when(
    H2O_vs_MEOH_t_test_P_FDR < 0.05 & log2_fc < -1 ~ "Down", # significantly down
    H2O_vs_MEOH_t_test_P_FDR < 0.05 & log2_fc > 1 ~ "Up"))   # significantly up
vol_neg_data$point_color[is.na(vol_neg_data$point_color)] <- "Not sig"
# Extracting feature identified
metab_data_neg <- noqc_vol_neg[!is.na(noqc_vol_neg@featureData@data$Metabolite),]
# Extracting metabolite table
meta_table_neg <- metab_data_neg@featureData@data
# Creating a new small table of the annotated compounds
# Keep metabolites with p-value < 0.05
volc_compouds_neg <- subset(vol_neg_data, H2O_vs_MEOH_t_test_P_FDR < 0.05)
volc_compouds_neg <- left_join(meta_table_neg, volc_compouds_neg)
# Volcano plot
vol_neg_plot <- ggplot(vol_neg_data, aes(log2_fc, logP)) +
  scale_colour_manual(values = c("Not sig" = "#E5E5E5",
                                 "Down" = "#3366FF",
                                 "Up" = "#009933")) +
  theme_classic() +
  geom_point(data = vol_neg_data,
             aes(shape = Statistic_test,
                 color = point_color),
             size = 2.5,
             alpha = 0.4)  +
  guides(shape = "none") +
  labs(shape = 'Statistical test',
       color = 'Significance') +
  ggrepel::geom_label_repel(data = volc_compouds_neg,
                            aes(label = meta_table_neg$Metabolite),
                            color = "black",
                            box.padding = 0.37,
                            label.padding = 0.22,
                            label.r = 0.30,
                            cex = 2.5,
                            max.overlaps = 20,
                            min.segment.length = 0,
                            seed = 42) +
  xlab(bquote(Log[2](Methanol/Water))) + ylab(bquote(-Log[10](p-value))) +
  scale_y_continuous(limits = c(0,10), labels = label_scientific(),
                     expand=c(0.003, 0.003)) +
  theme(legend.position = c(0.85, 0.90),
        legend.background = element_rect(fill = "white", color = "black")) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill= "transparent")) +
  geom_vline(xintercept = 0, linetype = "longdash", colour="gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "longdash", colour="gray")
vol_neg_plot

########################## Join the volcano plots ##############################

# Figure matrix
volcano_plot <- arrangeGrob(vol_neg_plot,
                            vol_pos_plot,
                        layout_matrix = rbind(c(1, 2)))
# Adding label to the figures
volcano_one <- ggpubr::as_ggplot(volcano_plot) +
  draw_plot_label(label = LETTERS[1:2],
                  x = c(0, 0.5),
                  y = c(.99, .99))
# Exporting (*.pdf) file
ggsave(filename = "../Vibrio-diabolicus-Ili-/Plots/pdf/Volcano_plot.pdf", plot = volcano_one,
      width = 100, height = 70, units = "mm", dpi = 300, scale = 2.5)
# Exporting (*.png) file
ggsave(filename = "../Vibrio-diabolicus-Ili-/Plots/png/Volcano_plot.png", plot = volcano_one,
      width = 100, height = 70, units = "mm", dpi = 300, scale = 2.5)
# Exporting (*.jpg) file
#ggsave(filename = "Result/notame_Result/figure_1_glog.jpg", plot = volcano_one,
#      width = 175, height = 120, units = "mm", dpi = 300, scale = 2.5)








