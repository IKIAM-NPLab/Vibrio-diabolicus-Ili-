
############ Volcano plot of the methanol water shared features ################

# Loading necessary packages
library(notame)
library(dplyr)

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
                             var.equal = TRUE, paired = FALSE)
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
volc_compouds_pos <- subset(vol_pos_data, H2O_vs_MEOH_t_test_P_FDR < 0.05 |
                             point_color == "Fold change")
volc_compouds <- left_join(meta_table, volc_compouds)
# Volcano plot
vol_pos_plot <- ggplot(vol_pos_data, aes(log2_fc, logP)) +
  scale_colour_manual(values = c("Not sig" = "#808080",
                                 "Down" = "#4F9D4EFF",
                                 "Up" = "#ff80ff")) +
  theme_classic() +
  geom_point(data = vol_pos_data,
             aes(shape = Statistic_test,
                 color = point_color),
             size = 2.5,
             alpha = 0.4)  +
  labs(shape = 'Statistical test',
       color = 'Significance') +
  ggrepel::geom_label_repel(data = sp_volc_compouds,
                            aes(label = meta_table$Metabolite),
                            color = "black",
                            box.padding = 0.37,
                            label.padding = 0.22,
                            label.r = 0.30,
                            cex = 2.5,
                            max.overlaps = 20,
                            min.segment.length = 0,
                            seed = 42)
vol_pos_plot











