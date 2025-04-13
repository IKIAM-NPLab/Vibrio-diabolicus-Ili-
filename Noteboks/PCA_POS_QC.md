## Introduction

The present document aims to record the procedure given for the
statistical analysis of secondary metabolites present in the different
conditions of *bioles*. For each step a brief explanation, the code and
graphics obtained are included.

The workflow used is taken from the paper [“notame”: Workflow for
Non-Targeted LC–MS Metabolic
Profiling](https://doi.org/10.3390/metabo10040135). Which offers a wide
variety of functions to perform metabolomic profile analysis.

## Before to start

The “notame” package accepts as input a feature table that can be
obtained through software such as MZMine, MSDial, among others. In this
case, the table was obtained with the help of MZmine. The (\*.txt) file
was slightly modified to obtain the feature table.

Modifications made to the raw (\*.txt) file can be summarized in adding
and renaming columns. The added columns “Column” and “Ion Mode” allow to
analyze samples with different types of columns and with different
ionization modes respectively. Also, the cells corresponding to mass and
retention time must be renamed so that the package can detect and
process it.

## Notame workflow

As a first step for the analysis, all the necessary libraries were
installed and loaded in Rstudio.

``` r
# if (!requireNamespace("devtools", quietly = TRUE)) {
# install.packages("devtools")

devtools::install_github("antonvsdata/notame")
library(notame)
library(Biobase)
library(BiocGenerics)
library(futile.logger)
library(ggplot2)
library(magrittr)
library(foreach)
library(iterators)
library(parallel)
library(doParallel)
# if (!require("BiocManager", quietly = TRUE))
 #  install.packages("BiocManager")
# BiocManager::install("pcaMethods")
library(pcaMethods)
library(patchwork)
library(cowplot)
library(Rtsne)
library(ggdendro)
library(dplyr)
library(readxl)
library(ggsci)
# install.packages("igraph")
```

Then, a log system was added to have a record of each process executed.

``` r
init_log(log_file = "Results/Results_POS_QC.txt")
```

    ## INFO [2025-04-10 18:21:33] Starting logging

Next, the MZmine suitable feature list was imported.

``` r
data <- read_from_excel(file = "Data/POS_QC.xlsx", sheet = 1, 
                        corner_row = 6, corner_column = "F", 
                        split_by = c("Column", "Ion Mode"))
```

Once the data is read, the next step was to create a MetaboSet in order
to obtain a specific R object.

``` r
modes_pos <- construct_metabosets(exprs = data$exprs, 
                              pheno_data = data$pheno_data, 
                              feature_data = data$feature_data,
                              group_col = "Group")
```

We can visualize the raw data in order to inspect the processing
routines.

``` r
mode_pos <- modes_pos$RP_POS
Prueba_mode_pos <- modes_pos$RP_POS
POS_raw_sambx <- plot_sample_boxplots(Prueba_mode_pos, order_by = "Group", fill_by = "Group")
POS_raw_pca <- plot_pca(Prueba_mode_pos,
                        color = "Culture_media")
POS_raw_pca + POS_raw_sambx
```

![](PCA_POS_QC_files/figure-markdown_github/unnamed-chunk-5-1.png)

## Preprocessing

The first step of the preprocessing is to change the features with value
equal to 0 to NA.

``` r
mode_pos <- mark_nas(mode_pos, value = 0)
```

Then, features with low detection rate are first flagged and then will
be removed. The notame package employs two criteria to select this
features. First, is the feature presence in a percentage of QC
injections, and then the feature presence in a percentage within a
sample group or class.

``` r
mode_pos <- flag_detection(mode_pos, qc_limit = 0.70, group_limit = 0.75)
```

With these values, features which that were not detected in the 70% of
the QC injections and 75% of sample groups will be discarded.

The next step for preprocessing correspond to drift correction. The
drift correction can be applied by smoothed cubic spline regression.

``` r
corrected_pos <- correct_drift(mode_pos)
corrected_pos <- flag_quality(corrected_pos, condition = "RSD_r < 0.27 & D_ratio_r < 0.56")
```

Then we can visualize the data after drift correction.

``` r
EI_corr_sambx_pos <- plot_sample_boxplots(corrected_pos, order_by = "Group", fill_by = "Group")
EI_corr_pca_pos <- plot_pca(corrected_pos, center = T)
EI_corr_pca_pos + EI_corr_sambx_pos
```

![](PCA_POS_QC_files/figure-markdown_github/unnamed-chunk-9-1.png)

``` r
# Tercer método de corrección del efecto batch
 #Instalación del paquete RUVSeq
#if (!require("BiocManager", quietly = TRUE))
    #install.packages("BiocManager")
#BiocManager::install("RUVSeq")

#library(RUVSeq)

#Corrección del efecto batch mediante método RUVSeq (Removal of Unwanted variation)
#replicates <- list(which(corrected$QC == "QC"))
#batch_corrected <- ruvs_qc(corrected[["exprs"]], batch = "Batch")
#corrected_batch <- normalize_batches(corrected, batch = "Batch", group = "QC", ref_label = "QC")
```

The next step is feature clustering. This step helps us reduce the
number of features of the same molecule that were split due to
ionization behavior (In-source fragmentation for example).

``` r
clustered_pos <- cluster_features(corrected_pos,
                              rt_window = 1/60,
                              all_features = FALSE,
                              corr_thresh = 0.90,
                              d_thresh = 0.85
                              #plotting = TRUE,
                              #prefix = paste0(ppath, "Cluster/LCMS/LCMS_Cluster")
                              )
```

    ## INFO [2025-04-10 18:21:46] Identified m/z column mass and retention time column RT
    ## INFO [2025-04-10 18:21:46] 
    ## Starting feature clustering at 2025-04-10 18:21:46.445006
    ## INFO [2025-04-10 18:21:46] Finding connections between features in RP_POS
    ## [1] 100
    ## [1] 200
    ## [1] 300
    ## INFO [2025-04-10 18:21:53] Found 2245 connections in RP_POS
    ## INFO [2025-04-10 18:21:53] Found 2245 connections
    ## 46 components found
    ## 
    ## 19 components found
    ## 
    ## 2 components found
    ## 
    ## 2 components found
    ## 
    ## INFO [2025-04-10 18:21:53] Found 54 clusters of 2 or more features, clustering finished at 2025-04-10 18:21:53.679663

``` r
compressed_pos <- compress_clusters(clustered_pos)
```

    ## INFO [2025-04-10 18:21:53] Clusters compressed, left with 179 features

We can inspect data after clustering algorithm.

``` r
# PCA
compr_pca_pos <- plot_pca(compressed_pos,
                        center = TRUE,
                        shape = "Group",
                        color = "Group") 

compr_pca_pos
```

![](PCA_POS_QC_files/figure-markdown_github/unnamed-chunk-11-1.png)

The next step imputes the data.

``` r
# To clean data
set.seed(35)
imputed_pos <- impute_rf(compressed_pos)
```

    ## INFO [2025-04-10 18:21:54] 
    ## Starting random forest imputation at 2025-04-10 18:21:54.194337
    ## INFO [2025-04-10 18:22:00] Out-of-bag error in random forest imputation: 0.097
    ## INFO [2025-04-10 18:22:00] Random forest imputation finished at 2025-04-10 18:22:00.44602

We can inspect PCA plot after imputation.

``` r
# Boxplot
imp_bp_pos <- plot_sample_boxplots(imputed_pos,
                               order_by = "Group",
                               fill_by = "Group")
# PCA
imp_pca_pos <- plot_pca(imputed_pos,
                    center = TRUE,
                    shape = "Group",
                    color = "Group")
# Plot
imp_pca_pos + imp_bp_pos
```

![](PCA_POS_QC_files/figure-markdown_github/unnamed-chunk-13-1.png)

# Second PCA and loading plot

Droping flagged features

``` r
# Extract clean data
no_flag_pos <- drop_flagged(imputed_pos)
# Extracting feature height table
peak_height_pos <- exprs(no_flag_pos)
# Extracting samples information
pheno_data <- no_flag_pos@phenoData@data
# Extracting feature information
feat_data_pos <- no_flag_pos@featureData@data
```

Preparing data and transposing feature table.

``` r
# Transposing feature height table
transp_table_pos  <- t(peak_height_pos)

# Changing NA to 0 
transp_table_pos[is.na(transp_table_pos)]=0

# Centering and Scaling features
pos_pca <- prcomp(transp_table_pos, center = TRUE, scale. = TRUE)
```

Plotting PCA results.

``` r
# PCA scores
scores_pos <- pos_pca$x %>%                  # Get PC coordinates
  data.frame %>%                         # Convert to data frames
  mutate(Sample_ID = rownames(.)) %>%    # Create a new column with the sample names
  left_join(pheno_data)                 # Adding metadata

# PCA plot
ggplot(scores_pos,
       aes(PC1, PC2, shape = Group, color = Group)) +
  geom_point(size = 3) +
  guides(x=guide_axis(title = "PC1 (33.1 %)"),
         y=guide_axis(title = "PC2 (23.5 %)")) +
  #geom_text(label=pheno_data$Group,
  #          nudge_x = 1,
  #          nudge_y = 1,
  #          show.legend = FALSE) +
  theme_classic() +
#  theme(legend.text = element_text(face="italic")) +
  theme(legend.position = c(0.9, 0.7),
        legend.background = element_rect(fill = "white", color = "black")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  geom_vline(xintercept = 0, linetype = "longdash", colour="gray") +
  geom_hline(yintercept = 0, linetype = "longdash", colour="gray") 
```

![](PCA_POS_QC_files/figure-markdown_github/unnamed-chunk-16-1.png)

``` r
# Save plot
ggsave('Plots/PCA_POS_QC.pdf',
       width = 7, height = 6, device='pdf', dpi="print")
ggsave('Plots/PCA_POS_QC.png',
       width = 7, height = 6, device='png', dpi="print")
```

Plotting loading results.

``` r
loadings_pos <- pos_pca$rotation %>%    # Extract loadings
  data.frame(Cluster_ID = rownames(.))  # New column with feat name
```

Creating an artificial table with feature ID and compound name column.

``` r
# Exporting notame output to find and filter identified metabolites
write_to_excel(clustered_pos, "Results/Firts_features_name_POS_QC.xlsx")
# Load a metabolite name table
metab_name_pos <- readxl::read_excel("Data/Firts_POS_QC_Metabolites.xlsx", 1)
# Creating a new small table of the annotated compounds
pos_compouds <- left_join(metab_name_pos, loadings_pos)
```

``` r
ggplot(loadings_pos, aes(PC1, PC2)) +
  geom_point(alpha = 0.2) +
  theme_classic() +  # Aplicar tema base primero
  ggrepel::geom_label_repel(data = pos_compouds,
                            aes(label = Metabolite_name),
                            box.padding = 0.4,
                            label.padding = 0.1,
                            label.r = 0.1,
                            max.overlaps = 50,
                            # Cambiado 'cex' por 'size' y ajustado el valor
                            # Prueba con valores como 4, 5, 6, etc. hasta obtener el tamaño deseado
                            size = 5.23) +
  guides(x=guide_axis(title = "PC1 (33.1 %)"),
         y=guide_axis(title = "PC2 (23.5 %)")) +
  ggsci::scale_color_aaas() +
  # Añadir la capa theme() para modificar el tamaño de los títulos de los ejes
  theme(
    axis.title.x = element_text(size = 16), # Ajusta este valor (ej: 14, 16, 18)
    axis.title.y = element_text(size = 16)  # Ajusta este valor (ej: 14, 16, 18)
    # Opcionalmente, puedes ajustar también el tamaño del texto de las marcas de los ejes:
    # axis.text.x = element_text(size = 12),
    # axis.text.y = element_text(size = 12)
  ) + 
  scale_x_continuous(limits = c(-0.16, 0.16)) + # Reemplaza con tus valores deseados para X
  scale_y_continuous(limits = c(-0.16, 0.19)) +  # Reemplaza con tus valores deseados para Y
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  geom_vline(xintercept = 0, linetype = "longdash", colour="gray") +
  geom_hline(yintercept = 0, linetype = "longdash", colour="gray") 
```

![](PCA_POS_QC_files/figure-markdown_github/unnamed-chunk-19-1.png)

``` r
#Save plot
ggsave('Plots/LOADING_POS_QC.pdf',
       width = 15, height = 7, device='pdf', dpi="print")
ggsave('Plots/LOADING_POS_QC.png',
       width = 15, height = 7, device='png', dpi="print")
```

# Heat map plot

ComplexHeatmap package and dependency installation.

``` r
# ComplexHeatmap package installation
# if (!requireNamespace("BiocManager", quietly=TRUE))
 # install.packages("BiocManager")
# BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

# ColorRamp2 package installation
# if (!requireNamespace("devtools", quietly = TRUE)) 
  # install.packages("devtools")

# devtools::install_github("jokergoo/colorRamp2")
library(colorRamp2)

# Cowplot package installation
# install.packages("cowplot")
library(cowplot)
```

Extracting and loaded of identified metabolites abundance.

``` r
#"notame" output to find and filter height of identified metabolites
write_to_excel(no_flag_pos, "Results/parametaboanalyst_POS_QC.xlsx")
```

Extracting and loaded of identified metabolites abundance.

``` r
# Metabolite name table
metab_name_hm_pos <- readxl::read_excel("Results/IDLOADING_POS_QC.xlsx", 1)

hm_scl_pos <- metab_name_hm_pos[, 5:19] %>% as.matrix %>% log10()
rownames(hm_scl_pos) <- metab_name_hm_pos$Metabolite_name

# Metabolite classification
metab_class_pos <- metab_name_hm_pos %>% select(Class = Class, Metabolite = Metabolite_name)

# Metabolite class to HeatMap anotation
met_class_annotation_pos <-  metab_class_pos %>% select(Class) %>% 
  as.matrix()
rownames(met_class_annotation_pos) <- metab_class_pos$Metabolite

# Top information
top_info_pos <- data.frame(Species = c( rep("Pep"),
                                    rep("Yeast"),
                                    rep("Star_pep"),
                                    rep("Star_yeast"),
                                    rep("Try_yeast"),
                                    rep("Pep"),
                                    rep("Yeast"),
                                    rep("Star_pep"),
                                    rep("Star_yeast"),
                                    rep("Try_yeast"),
                                    rep("Pep"),
                                    rep("Yeast"),
                                    rep("Star_pep"),
                                    rep("Star_yeast"),
                                    rep("Try_yeast"))) 

rownames(top_info_pos) <- paste(top_info_pos$Species, rep(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)))
top_info_pos <- as.matrix(top_info_pos)
```

Scaling, row and top heatmap anotation.

``` r
set.seed(2024)
# Metabolite class color
cols_metclass_pos <- c("Steroids and steroid derivatives" = "#800000FF",
                   "Naphthopyrans" = "#EE6A50",
                   "Organonitrogen compounds" = "#FFF0F5",
                   "Lipids and lipid-like molecules" = "#97FFFF", 
                   "Benzopyrans" = "#FFDEAD", 
                   "Diazanaphthalenes" = "#C16622FF",
                   "Indoles and derivatives" = "#FFB5C5",
                   "Glycerolipids" = "#91cf60",
                   "Fatty Acyls" = "#87CEFF",
                   "Glycerophospholipids" = "#e6550d",
                   "Cinnamic acids and derivatives" = "#76EEC6",
                   "Pyrrolidines" = "#777")
                   
# Add row anotation to HeatMap
hm_row_ann_pos <- rowAnnotation(Metabolite = met_class_annotation_pos,
                            col = list(Metabolite = cols_metclass_pos),
                            show_annotation_name = T,
                            show_legend= F)

cols_species <- c("Pep" = "#FFDAB9",
                 "Yeast" = "#FFB5C5",
                 "Star_pep" = "#B9D3EE",
                 "Star_yeast" = "#EEE685",
                 "Try_yeast" = "#76EEC6",
                 "Pep" = "#FFDAB9",
                 "Yeast" = "#FFB5C5",
                 "Star_pep" = "#B9D3EE",
                 "Star_yeast" = "#EEE685",
                 "Try_yeast" = "#76EEC6",
                 "Pep" = "#FFDAB9",
                 "Yeast" = "#FFB5C5",
                 "Star_pep" = "#B9D3EE",
                 "Star_yeast" = "#EEE685",
                 "Try_yeast" = "#76EEC6", 
                 "Pep" = "#FFDAB9",
                 "Yeast" = "#FFB5C5",
                 "Star_pep" = "#B9D3EE",
                 "Star_yeast" = "#EEE685",
                 "Try_yeast" = "#76EEC6")


# Add top anotation to HeatMap

top_info_ann_pos <- HeatmapAnnotation(Species = top_info_pos,
                                  col = list(Species = cols_species),
                                  show_annotation_name = T,
                                  show_legend = F, 
                                  border = TRUE)
# Color scale

mycol_pos <- colorRamp2(c(-2, 5, 10),
                    c("blue", "white", "red"))

# Heatmap matrix plotting

hm_plot_pos <- Heatmap(hm_scl_pos,
        col = mycol_pos,
        border_gp = grid::gpar(col = "black", lty = 0.05),
        rect_gp = grid::gpar(col = "black", lwd = 0.75),
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "complete",
        top_annotation = top_info_ann_pos,
        right_annotation = hm_row_ann_pos,
        show_heatmap_legend = F,
        row_km = 3, column_km = 2)
hm_plot_pos
```

![](PCA_POS_QC_files/figure-markdown_github/unnamed-chunk-23-1.png)

Adding legends to heatmap.

``` r
# Color scale legend
lgd1_pos <- Legend(col_fun = mycol_pos,
               title = "log10 abundance",
               at = seq(6),
               direction = "horizontal" )
# Group legend
lgd2_pos <- Legend(labels = gt_render(c("Pep",
                                    "Yeast",
                                    "Star_pep",
                                    "Star_yeast",
                                    "Try_yeast")),
              legend_gp = gpar(fill = cols_species),
              title = 'Group', ncol = 1,
              direction = "horizontal" )
              
# Metabolite class Legend
lgd3_pos <- Legend(labels = c(unique(metab_class_pos$Class)) ,
               legend_gp = gpar(fill = cols_metclass_pos), 
               title = 'Metabolite class', ncol = 3,
               direction = "horizontal" )
```

ComplexHeatmap plot

``` r
set.seed(1540)
# Converting to ggplot
gg_heatmap_pos <- grid.grabExpr(draw(hm_plot_pos))
gg_heatmap_pos <- ggpubr::as_ggplot(gg_heatmap_pos)

# Legends
all_legends_pos <- packLegend(lgd1_pos, lgd2_pos, lgd3_pos, direction = "horizontal")
gg_legend_pos <- grid.grabExpr(draw(all_legends_pos))
gg_legend_fn_pos <- ggpubr::as_ggplot(gg_legend_pos)

# Heatmap plot
posqc_hm_pos <- plot_grid(gg_legend_fn_pos,
          gg_heatmap_pos, ncol = 1,
          rel_heights = c(0.195, 0.88))
posqc_hm_pos
```

![](PCA_POS_QC_files/figure-markdown_github/unnamed-chunk-25-1.png)

``` r
# Save heatmap plot
ggsave(filename = "Plots/Firts_QC_POS_Heatmap.pdf", plot = posqc_hm_pos,
      width = 9, height = 15, units = "in", dpi = 600, scale = 2)

ggsave(filename = "Plots/Firts_QC_POS_Heatmap.png", plot = posqc_hm_pos,
      width = 9, height = 15, units = "in", dpi = 600, scale = 2)
```

Finish a record.

``` r
finish_log()
```

    ## INFO [2025-04-10 18:22:29] Finished analysis. Thu Apr 10 18:22:29 2025
    ## Session info:
    ## 
    ## INFO [2025-04-10 18:22:29] R version 4.4.2 (2024-10-31 ucrt)
    ## INFO [2025-04-10 18:22:29] Platform: x86_64-w64-mingw32/x64
    ## INFO [2025-04-10 18:22:29] Running under: Windows 11 x64 (build 26100)
    ## INFO [2025-04-10 18:22:29] 
    ## INFO [2025-04-10 18:22:29] Matrix products: default
    ## INFO [2025-04-10 18:22:29] 
    ## INFO [2025-04-10 18:22:29] 
    ## INFO [2025-04-10 18:22:29] locale:
    ## INFO [2025-04-10 18:22:29] [1] LC_COLLATE=Spanish_Ecuador.utf8  LC_CTYPE=Spanish_Ecuador.utf8   
    ## INFO [2025-04-10 18:22:29] [3] LC_MONETARY=Spanish_Ecuador.utf8 LC_NUMERIC=C                    
    ## INFO [2025-04-10 18:22:29] [5] LC_TIME=Spanish_Ecuador.utf8    
    ## INFO [2025-04-10 18:22:29] 
    ## INFO [2025-04-10 18:22:29] time zone: America/Guayaquil
    ## INFO [2025-04-10 18:22:29] tzcode source: internal
    ## INFO [2025-04-10 18:22:29] 
    ## INFO [2025-04-10 18:22:29] attached base packages:
    ## INFO [2025-04-10 18:22:29] [1] grid      parallel  stats     graphics  grDevices utils     datasets 
    ## INFO [2025-04-10 18:22:29] [8] methods   base     
    ## INFO [2025-04-10 18:22:29] 
    ## INFO [2025-04-10 18:22:29] other attached packages:
    ## INFO [2025-04-10 18:22:29]  [1] colorRamp2_0.1.0      ComplexHeatmap_2.22.0 ggsci_3.2.0          
    ## INFO [2025-04-10 18:22:29]  [4] readxl_1.4.5          dplyr_1.1.4           ggdendro_0.2.0       
    ## INFO [2025-04-10 18:22:29]  [7] Rtsne_0.17            cowplot_1.1.3         patchwork_1.3.0      
    ## INFO [2025-04-10 18:22:29] [10] pcaMethods_1.98.0     doParallel_1.0.17     iterators_1.0.14     
    ## INFO [2025-04-10 18:22:29] [13] foreach_1.5.2         notame_0.3.2          magrittr_2.0.3       
    ## INFO [2025-04-10 18:22:29] [16] ggplot2_3.5.1         futile.logger_1.4.3   Biobase_2.66.0       
    ## INFO [2025-04-10 18:22:29] [19] BiocGenerics_0.52.0  
    ## INFO [2025-04-10 18:22:29] 
    ## INFO [2025-04-10 18:22:29] loaded via a namespace (and not attached):
    ## INFO [2025-04-10 18:22:29]   [1] RColorBrewer_1.1-3   rstudioapi_0.17.1    shape_1.4.6.1       
    ## INFO [2025-04-10 18:22:29]   [4] magick_2.8.5         farver_2.1.2         rmarkdown_2.29      
    ## INFO [2025-04-10 18:22:29]   [7] GlobalOptions_0.1.2  fs_1.6.5             ragg_1.3.3          
    ## INFO [2025-04-10 18:22:29]  [10] vctrs_0.6.5          memoise_2.0.1        askpass_1.2.1       
    ## INFO [2025-04-10 18:22:29]  [13] rstatix_0.7.2        htmltools_0.5.8.1    usethis_3.1.0       
    ## INFO [2025-04-10 18:22:29]  [16] itertools_0.1-3      missForest_1.5       lambda.r_1.2.4      
    ## INFO [2025-04-10 18:22:29]  [19] curl_6.2.1           broom_1.0.7          cellranger_1.1.0    
    ## INFO [2025-04-10 18:22:29]  [22] Formula_1.2-5        htmlwidgets_1.6.4    futile.options_1.0.1
    ## INFO [2025-04-10 18:22:29]  [25] cachem_1.1.0         commonmark_1.9.2     igraph_2.1.4        
    ## INFO [2025-04-10 18:22:29]  [28] mime_0.12            lifecycle_1.0.4      pkgconfig_2.0.3     
    ## INFO [2025-04-10 18:22:29]  [31] R6_2.6.1             fastmap_1.2.0        shiny_1.10.0        
    ## INFO [2025-04-10 18:22:29]  [34] clue_0.3-66          digest_0.6.37        colorspace_2.1-1    
    ## INFO [2025-04-10 18:22:29]  [37] S4Vectors_0.44.0     pkgload_1.4.0        textshaping_1.0.0   
    ## INFO [2025-04-10 18:22:29]  [40] ggpubr_0.6.0         labeling_0.4.3       randomForest_4.7-1.2
    ## INFO [2025-04-10 18:22:29]  [43] abind_1.4-8          compiler_4.4.2       rngtools_1.5.2      
    ## INFO [2025-04-10 18:22:29]  [46] remotes_2.5.0        withr_3.0.2          backports_1.5.0     
    ## INFO [2025-04-10 18:22:29]  [49] carData_3.0-5        pkgbuild_1.4.6       ggsignif_0.6.4      
    ## INFO [2025-04-10 18:22:29]  [52] MASS_7.3-65          openssl_2.3.2        sessioninfo_1.2.3   
    ## INFO [2025-04-10 18:22:29]  [55] rjson_0.2.23         tools_4.4.2          zip_2.3.2           
    ## INFO [2025-04-10 18:22:29]  [58] httpuv_1.6.15        glue_1.8.0           promises_1.3.2      
    ## INFO [2025-04-10 18:22:29]  [61] gridtext_0.1.5       cluster_2.1.8        generics_0.1.3      
    ## INFO [2025-04-10 18:22:29]  [64] gtable_0.3.6         tidyr_1.3.1          car_3.1-3           
    ## INFO [2025-04-10 18:22:29]  [67] xml2_1.3.7           ggrepel_0.9.6        pillar_1.10.1       
    ## INFO [2025-04-10 18:22:29]  [70] markdown_1.13        stringr_1.5.1        later_1.4.1         
    ## INFO [2025-04-10 18:22:29]  [73] circlize_0.4.16      tidyselect_1.2.1     miniUI_0.1.1.1      
    ## INFO [2025-04-10 18:22:29]  [76] knitr_1.49           IRanges_2.40.1       stats4_4.4.2        
    ## INFO [2025-04-10 18:22:29]  [79] xfun_0.51            devtools_2.4.5       credentials_2.0.2   
    ## INFO [2025-04-10 18:22:29]  [82] matrixStats_1.5.0    stringi_1.8.4        yaml_2.3.10         
    ## INFO [2025-04-10 18:22:29]  [85] evaluate_1.0.3       codetools_0.2-20     tibble_3.2.1        
    ## INFO [2025-04-10 18:22:29]  [88] cli_3.6.4            xtable_1.8-4         systemfonts_1.2.1   
    ## INFO [2025-04-10 18:22:29]  [91] munsell_0.5.1        Rcpp_1.0.14          gert_2.1.4          
    ## INFO [2025-04-10 18:22:29]  [94] png_0.1-8            ellipsis_0.3.2       doRNG_1.8.6.1       
    ## INFO [2025-04-10 18:22:29]  [97] profvis_0.4.0        urlchecker_1.0.1     viridisLite_0.4.2   
    ## INFO [2025-04-10 18:22:29] [100] scales_1.3.0         openxlsx_4.2.8       purrr_1.0.4         
    ## INFO [2025-04-10 18:22:29] [103] crayon_1.5.3         GetoptLong_1.0.5     rlang_1.1.5         
    ## INFO [2025-04-10 18:22:29] [106] formatR_1.14

``` r
#Finish a record.
```
