#-------------------------------------------------------------------------------
# 1. GLOBAL SETUP, LIBRARIES & DIRECTORY STRUCTURE
#-------------------------------------------------------------------------------
cat("\n--- 1. INITIALIZING ANALYSIS PIPELINE ---")

# --- 1a. Library Loading & Parallel Setup ---
cat("\n- Loading required packages...")
required_packages <- c(
  "data.table", "ape", "bayesplot", "cluster", "DEqMS", "dplyr", "factoextra", 
  "foreach", "future", "future.apply", "gdata", "ggplot2", "ggrepel", 
  "gridExtra", "itertools", "missForest", "pheatmap", "purrr", 
  "randomForest", "readr", "reshape2", "rstanarm", "stringr", 
  "tibble", "tidyr", "vegan", "VennDiagram", "WGCNA",
  "AnnotationDbi", "Biostrings", "clusterProfiler", "enrichplot", 
  "GO.db", "impute", "imputeLCMD", "limma", "MSstats", 
  "pcaMethods", "preprocessCore"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    cat(paste0("\n  > Installing: ", pkg, "..."))
    if (pkg %in% c("DEqMS", "limma", "clusterProfiler", "Biostrings", "AnnotationDbi", 
                   "enrichplot", "impute", "imputeLCMD", "MSstats", 
                   "pcaMethods", "preprocessCore")) {
      if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg, dependencies = TRUE)
    }
    library(pkg, character.only = TRUE)
  }
}

if (require(WGCNA)) allowWGCNAThreads()
plan(multisession, workers = parallel::detectCores() - 1)

# --- 1b. Main Project Path & Dynamic Output Folder ---
project_dir <- "."   # folder containing the input and configuration files
exp_id      <- "proteomics_analysis"   # identifier used for output folder and FASTA names
output_main <- file.path(project_dir, paste0(exp_id, "_Results"))

# --- 1c. Directory Generation ---
cat("\n- Generating directory structure...")
subfolders <- c(
  "1_Cleaned_File", "2_QC_correction", "3_EDA_Pre_Normalization",
  "4_Normalization", "5_Missing_Values", "6_DEA",
  "7_CoExpression_Analysis", "8_EDA_and_CoExpression",
  "9_CBR_Filter", "10_ORA", "11_GSEA", "12_Validation_RF"
)

if(!dir.exists(output_main)) dir.create(output_main, recursive = TRUE)
for (f in subfolders) {
  path <- file.path(output_main, f)
  if(!dir.exists(path)) dir.create(path, recursive = TRUE)
}

log_file <- file.path(output_main, paste0(date_suffix, "_Analysis_Log.txt"))
write(paste("--- ANALYSIS LOG -", Sys.time(), "---"), file = log_file)
write("===========================================\nEXECUTION REPORT\n===========================================", file = log_file, append = TRUE)

# --- 1d. Smart File Detection ---
find_file <- function(pattern) {
  f <- list.files(project_dir, pattern = pattern, full.names = TRUE)
  if(length(f) == 0) stop(paste("[ERROR] Missing file for pattern:", pattern))
  return(sort(f, decreasing = TRUE)[1])
}

cat("\n- Detecting configuration files...")
path_metadata  <- find_file("experimental_metadata.*\\.csv$")
path_contrasts <- find_file("experimental_contrasts.*\\.csv$")
path_params    <- find_file("analysis_params.*\\.csv$")
path_desc      <- find_file("protein_descriptions.*\\.tsv$")
path_fasta     <- find_file("protein_database.*\\.fasta$")
input_raw      <- find_file("proteomic_results.*\\.tsv$")

# --- 1e. Load Configuration & Parameters ---
metadata        <- as_tibble(fread(path_metadata))
params_df       <- as_tibble(fread(path_params))
uniprot_db      <- as_tibble(fread(path_desc)) %>% rename(Protein_ID = Entry)
contrast_matrix <- as_tibble(fread(path_contrasts)) %>% filter(toupper(Active) == "YES")
params          <- setNames(params_df$Value, params_df$Parameter)

# --- DYNAMIC CPU CALCULATION ---
total_physical_cores <- parallel::detectCores(logical = FALSE)
n_cores_final        <- max(1, floor(total_physical_cores * as.numeric(params["cpu_usage_ratio"])))

cat(paste0("\n[SYSTEM] Unifying pipeline to use: ", n_cores_final, " cores."))

# --- 1f. Auto-Detecting TSV Columns ---
headers <- colnames(fread(input_raw, nrow = 0))
pick_col <- function(candidates, available) {
  match <- intersect(candidates, available); if(length(match) > 0) return(match[1]) else return(NULL)
}

protein_id_col <- pick_col(c("PG.ProteinAccessions", "Majority protein IDs", "Protein.Ids"), headers)
intensity_col  <- pick_col(c("PG.Quantity", "Intensity", "Raw.Intensity", "Area"), headers)
# Crucial: Finding the link to Raw_FileName from metadata
raw_file_col   <- pick_col(c("R.FileName", "Raw.File", "File.Name"), headers)
peptide_col    <- pick_col(c("PEP.StrippedSequence", "Peptides", "Unique Peptides"), headers)

# --- 1g. Mapping Thresholds ---
min_rep_pct <- as.numeric(params["min_rep_pct"]); pep_min <- as.numeric(params["pep_min"])

# --- 1h. Professional EDA Function (The "Master Suite" with UpSet & Venn) ---
generate_full_eda_analysis <- function(wide_data, subfolder, suffix, title_prefix) {
  cat(paste0("\n    > Processing Visual Master Suite: ", title_prefix))
  target_dir <- file.path(output_main, subfolder)
  if(!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)
  
  # --- 1. PREP DATA ---
  all_int_cols <- grep("^Int_", colnames(wide_data), value = TRUE)
  data_long <- wide_data %>%
    select(all_of(protein_id_col), all_of(all_int_cols)) %>%
    pivot_longer(-!!sym(protein_id_col), names_to = "Full_ID", values_to = "Intensity") %>%
    mutate(SampleID = gsub("Int_", "", Full_ID), 
           Condition = sub("_R[0-9]+$", "", SampleID)) %>%
    filter(!is.na(Intensity) & Intensity > 0)
  
  sample_metadata <- data_long %>% select(SampleID, Condition) %>% distinct() %>% as.data.frame()
  rownames(sample_metadata) <- sample_metadata$SampleID
  
  # Matrix for multivariate (Log2 transformed, omitting NAs)
  wide_mat <- data_long %>% mutate(Log2Int = log2(Intensity)) %>%
    select(!!sym(protein_id_col), SampleID, Log2Int) %>%
    pivot_wider(names_from = SampleID, values_from = Log2Int) %>%
    column_to_rownames(protein_id_col)
  mat_no_na <- na.omit(wide_mat)
  
  # --- A. BOXPLOT ---
  p_box <- ggplot(data_long, aes(x = SampleID, y = log2(Intensity), fill = Condition)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5) + theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste(title_prefix, "- Intensity Distribution"), y = "Log2 Intensity")
  ggsave(file.path(target_dir, paste0("1_Boxplot_", suffix, ".tiff")), p_box, width=12, height=6, dpi=300, compression="lzw")
  
  # --- B. CV PLOT (Technical Consistency) ---
  cv_data <- data_long %>%
    group_by(!!sym(protein_id_col), Condition) %>%
    summarise(CV = (sd(Intensity, na.rm=T) / mean(Intensity, na.rm=T)) * 100, .groups = "drop")
  
  p_cv <- ggplot(cv_data, aes(x = CV, color = Condition, fill = Condition)) +
    geom_density(alpha = 0.2) + theme_minimal() + xlim(0, 100) +
    labs(title = paste(title_prefix, "- Coefficient of Variation"), x = "CV (%)", y = "Density")
  ggsave(file.path(target_dir, paste0("2_CV_Plot_", suffix, ".tiff")), p_cv, width=8, height=6, dpi=300, compression="lzw")
  
  # --- C. PCA + PERMANOVA + PERMDISP ---
  if(nrow(mat_no_na) > 5) {
    pca_res <- prcomp(t(mat_no_na), scale. = TRUE)
    pc_vars <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)
    dist_eucl <- dist(t(mat_no_na))
    perm_pca <- vegan::adonis2(dist_eucl ~ Condition, data = sample_metadata)
    disp_pca <- vegan::permutest(vegan::betadisper(dist_eucl, sample_metadata$Condition), pairwise=F)$tab$`Pr(>F)`[1]
    
    pca_df <- as.data.frame(pca_res$x) %>% rownames_to_column("SampleID") %>% left_join(sample_metadata, by = "SampleID")
    p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition, label = SampleID)) +
      geom_point(size = 4) + ggrepel::geom_text_repel(size = 2.5) + theme_bw() +
      labs(title = paste(title_prefix, "- PCA"), 
           subtitle = paste0("PERMANOVA p=", round(perm_pca$`Pr(>F)`[1],4), " | PERMDISP p=", round(disp_pca,4)),
           x = paste0("PC1 (", pc_vars[1], "%)"), y = paste0("PC2 (", pc_vars[2], "%)"))
    ggsave(file.path(target_dir, paste0("3_PCA_", suffix, ".tiff")), p_pca, width=9, height=7, dpi=300, compression="lzw")
  }
  
  # --- D. PCoA (Bray-Curtis) ---
  if(nrow(mat_no_na) > 5) {
    dist_b <- vegan::vegdist(t(mat_no_na), method = "bray")
    pcoa_res <- cmdscale(dist_b, k = 2, eig = TRUE)
    pcoa_vars <- round(pcoa_res$eig / sum(pcoa_res$eig) * 100, 1)
    perm_pcoa <- vegan::adonis2(dist_b ~ Condition, data = sample_metadata)
    
    pcoa_df <- as.data.frame(pcoa_res$points) %>% setNames(c("D1", "D2")) %>% 
      rownames_to_column("SampleID") %>% left_join(sample_metadata, by = "SampleID")
    p_pcoa <- ggplot(pcoa_df, aes(x = D1, y = D2, color = Condition, label = SampleID)) +
      geom_point(size = 4) + ggrepel::geom_text_repel(size = 2.5) + theme_bw() +
      labs(title = paste(title_prefix, "- PCoA (Bray-Curtis)"), 
           subtitle = paste0("PERMANOVA p=", round(perm_pcoa$`Pr(>F)`[1], 4)),
           x = paste0("PCoA1 (", pcoa_vars[1], "%)"), y = paste0("PCoA2 (", pcoa_vars[2], "%)"))
    ggsave(file.path(target_dir, paste0("4_PCoA_", suffix, ".tiff")), p_pcoa, width=9, height=7, dpi=300, compression="lzw")
  }
  
  # --- E. HEATMAP (Top 50 Variantes) ---
  top_var <- mat_no_na[order(apply(mat_no_na, 1, sd), decreasing = TRUE)[1:min(50, nrow(mat_no_na))], ]
  pheatmap::pheatmap(top_var, scale = "row", show_colnames = T, show_rownames = F,
                     main = paste(title_prefix, "- Top 50 Most Variable Proteins"),
                     filename = file.path(target_dir, paste0("5_Heatmap_Top50_", suffix, ".tiff")), 
                     width = 10, height = 8)
  
  # --- F. VENN & UPSET (Overlap Analysis) ---
  list_overlap <- lapply(unique(sample_metadata$Condition), function(cond){
    cols <- grep(paste0("Int_", cond, "_"), colnames(wide_data), value = TRUE)
    wide_data[[protein_id_col]][rowSums(!is.na(wide_data[,cols, drop=F])) > 0]
  })
  names(list_overlap) <- unique(sample_metadata$Condition)
  
  # 1. Venn Diagram (Only for 5 groups or less)
  if(length(list_overlap) <= 5) {
    tiff(file.path(target_dir, paste0("6_Venn_", suffix, ".tiff")), width=2400, height=2400, res=300)
    grid::grid.newpage()
    venn_p <- VennDiagram::venn.diagram(x = list_overlap, filename = NULL, fill = rainbow(length(list_overlap)), alpha = 0.5)
    grid::grid.draw(venn_p)
    dev.off()
  }
  
  # 2. UpSet Plot (Universal, works best for > 5 groups)
  if(require(UpSetR) && length(list_overlap) >= 2) {
    tiff(file.path(target_dir, paste0("6_UpSet_", suffix, ".tiff")), width=3000, height=2000, res=300)
    print(UpSetR::upset(fromList(list_overlap), order.by = "freq", nsets = length(list_overlap), 
                        main.bar.color = "steelblue", sets.bar.color = "darkred"))
    dev.off()
  }
  
  # --- G. MA-PLOT ---
  if(length(unique(sample_metadata$Condition)) >= 2) {
    g1 <- unique(sample_metadata$Condition)[1]
    g2 <- unique(sample_metadata$Condition)[2]
    ma_df <- wide_data %>%
      mutate(M1 = rowMeans(select(., starts_with(paste0("Int_", g1))), na.rm=T),
             M2 = rowMeans(select(., starts_with(paste0("Int_", g2))), na.rm=T),
             A = (log2(M1+1) + log2(M2+1))/2, M = log2(M1+1) - log2(M2+1)) %>% filter(!is.na(M))
    
    p_ma <- ggplot(ma_df, aes(x = A, y = M)) + geom_point(alpha=0.3) + geom_hline(yintercept=0, color="red") +
      theme_minimal() + labs(title = paste(title_prefix, "- MA Plot (", g1, " vs ", g2, ")"))
    ggsave(file.path(target_dir, paste0("7_MA_Plot_", suffix, ".tiff")), p_ma, width=8, height=6, dpi=300)
  }
  
  cat(paste0("\n    [SUCCESS] Figures saved in: ", subfolder))
  return(invisible(NULL))
}
#-------------------------------------------------------------------------------
# 2. DATA LOADING AND METADATA-BASED RENAMING
#-------------------------------------------------------------------------------
cat("\n--- 2. DATA LOADING AND METADATA-BASED RENAMING ---")

# 1. Load Raw Data
raw_data <- as_tibble(fread(input_raw)) %>% mutate(across(where(is.character), trimws))

# 2. TRACEABILITY: Initial Count
initial_proteins_raw <- n_distinct(raw_data[[protein_id_col]])
cat("\n- Initial proteins detected in raw file:", initial_proteins_raw)

# 3. INTEGRAL RENAMING (METADATA-DRIVEN)
# We use Metadata as the bridge.
cat("\n- Mapping Raw File Names to Sample IDs...")

# Create mapping dictionary: Raw_FileName -> Sample_ID
mapping_dict <- metadata %>%
  select(Raw_FileName, Sample_ID) %>%
  deframe()

# Filter and Rename based on metadata entries
raw_mapped <- raw_data %>%
  filter(!!sym(raw_file_col) %in% names(mapping_dict)) %>%
  mutate(Mapped_SampleID = mapping_dict[!!sym(raw_file_col)])

# Verify mapping success
cat("\n- Samples successfully mapped from metadata:", n_distinct(raw_mapped$Mapped_SampleID))

# 4. Intensity & Peptide Consolidation
cat("\n- Expanding protein groups and calculating means...")
data_expanded <- raw_mapped %>%
  separate_rows(!!sym(protein_id_col), sep = ";") %>%
  mutate(!!sym(protein_id_col) := trimws(!!sym(protein_id_col))) %>%
  group_by(!!sym(protein_id_col), Mapped_SampleID) %>%
  summarise(
    Val_Int = mean(!!sym(intensity_col), na.rm = TRUE),
    Val_Pep = n_distinct(!!sym(peptide_col)), 
    .groups = "drop"
  )

# 5. Pivoting to Wide Format (Using prefixes Int_ and Peps_)
intensity_wide <- data_expanded %>%
  select(!!sym(protein_id_col), Mapped_SampleID, Val_Int) %>%
  pivot_wider(names_from = Mapped_SampleID, values_from = Val_Int, names_prefix = "Int_")

peptide_counts_wide <- data_expanded %>%
  select(!!sym(protein_id_col), Mapped_SampleID, Val_Pep) %>%
  pivot_wider(names_from = Mapped_SampleID, values_from = Val_Pep, names_prefix = "Peps_")

# Create Master Unfiltered Table
master_table_unfiltered <- left_join(intensity_wide, peptide_counts_wide, by = protein_id_col)
n_initial_prots <- nrow(master_table_unfiltered)

cat("\n- Master table successfully created with", ncol(master_table_unfiltered), "columns.")
cat("\n- Total proteins after expansion:", n_initial_prots)

# Save the raw master file
write_csv(master_table_unfiltered, file.path(output_main, "1_Cleaned_File", "1_Master_Unfiltered_Renamed.csv"))

cat("\n[SUCCESS] Data renamed according to metadata. Ready for filtering.")

#-------------------------------------------------------------------------------
# 3. DATA CLEANING, BIOLOGICAL FILTERING & EDA
#-------------------------------------------------------------------------------
cat("\n--- 3. DATA CLEANING AND BIOLOGICAL FILTERING ---")

# 3a. Initial Clean: Contaminants & Reverse entries
cat("\n- Removing contaminants and reverse sequences...")
clean_step1 <- master_table_unfiltered %>% 
  filter(!grepl("CON__|REV__", !!sym(protein_id_col)))

n_after_contaminants <- nrow(clean_step1)
n_elim_contaminants <- n_initial_prots - n_after_contaminants
cat(paste0("\n  > Contaminants/Reverse removed: ", n_elim_contaminants))

# 3b. Peptide Filter (Metadata-driven column matching)
cat(paste0("\n- Applying Peptide Filter (min ", pep_min, " peptides per protein)..."))
all_int_cols <- grep("^Int_", colnames(clean_step1), value = TRUE)

for(int_col in all_int_cols) {
  # Match Int_SampleID with Peps_SampleID
  pep_col <- gsub("Int_", "Peps_", int_col)
  if(pep_col %in% colnames(clean_step1)) {
    clean_step1[[int_col]] <- ifelse(clean_step1[[pep_col]] < pep_min | is.na(clean_step1[[pep_col]]), 
                                     NA, clean_step1[[int_col]])
  }
}

#-------------------------------------------------------------------------------
# 3c. Presence Filter (Consistency within biological groups) - FIXED DIMENSIONS
#-------------------------------------------------------------------------------
noise_cutoff <- 0.20   
min_rep_pct  <- 0.50   

cat(paste0("\n- Cleaning noise (< ", noise_cutoff*100, "%) and applying consistency filter (", min_rep_pct*100, "%)..."))

bio_groups <- setdiff(unique(gsub("Int_|_R[0-9]+$", "", all_int_cols)), "QC")

# --- STEP 1: CELL-LEVEL NOISE REMOVAL ---
# We clean the cells of the CURRENT clean_step1
for (g in bio_groups) {
  group_cols <- grep(paste0("^Int_", g, "_R[0-9]+$"), colnames(clean_step1), value = TRUE)
  
  if(length(group_cols) > 0) {
    # We use rowSums on the current version of the matrix
    presence_ratio <- rowSums(!is.na(clean_step1[, group_cols])) / length(group_cols)
    noise_rows <- which(presence_ratio > 0 & presence_ratio < noise_cutoff)
    
    if(length(noise_rows) > 0) {
      clean_step1[noise_rows, group_cols] <- NA
    }
  }
}

# --- STEP 2: ROW-LEVEL CONSISTENCY FILTER ---
# Crucial: We calculate keep_idx based ONLY on the current clean_step1 dimensions
current_intensities <- clean_step1[, all_int_cols]

keep_idx <- apply(current_intensities, 1, function(row) {
  any(sapply(bio_groups, function(g) {
    # We find columns within the 'row' vector itself to avoid indexing errors
    group_cols_idx <- grep(paste0("^Int_", g, "_R[0-9]+$"), names(row))
    if(length(group_cols_idx) == 0) return(FALSE)
    
    valid_count <- sum(!is.na(row[group_cols_idx]))
    (valid_count / length(group_cols_idx)) >= min_rep_pct
  }))
})

# Now dimensions will match perfectly (1816 vs 1816)
clean_step1 <- clean_step1[keep_idx, ]

cat(paste0("\n  [DONE] Filtered to ", nrow(clean_step1), " consistent proteins."))

#-------------------------------------------------------------------------------
# 3d. Pre-Normalization Exploratory Data Analysis (EDA)
#-------------------------------------------------------------------------------
# First, create the "bridge" variable so the rest of your script doesn't break
intensity_matrix_final <- clean_step1
n_final_cleaned <- nrow(intensity_matrix_final)

cat("\n- Generating Pre-Normalization EDA Figures (PCA, PCoA, Stats)...")

# EDA figures
eda_raw_results <- generate_full_eda_analysis(
  wide_data = master_table_unfiltered, 
  subfolder = "1_Cleaned_File/EDA_Initial", 
  suffix = "RawData", 
  title_prefix = "Initial Raw Data"
)

#-------------------------------------------------------------------------------
# 3e. Final Exports & Traceability Log
#-------------------------------------------------------------------------------
# Using your original target directory for the matrix export
target_dir_clean <- file.path(output_main, "1_Cleaned_File")
if(!dir.exists(target_dir_clean)) dir.create(target_dir_clean, recursive = TRUE)

write_csv(intensity_matrix_final, file.path(target_dir_clean, "Filtered_Protein_Matrix.csv"))

# Update log for the automated report
# (Ensuring all variables like n_elim_contaminants exist or are calculated)
n_elim_total <- n_initial_prots - n_final_cleaned

write(paste0("\n--- STEP 4: DATA CLEANING REPORT ---",
             "\n- Initial Proteins: ", n_initial_prots,
             "\n- Total Removed Proteins: ", n_elim_total,
             "\n  (Including Noise Cleaned < ", noise_cutoff*100, "% and Consistency < ", min_rep_pct*100, "%)",
             "\n- Final Cleaned Proteins: ", n_final_cleaned), 
      file = log_file, append = TRUE)

cat(paste0("\n[SUCCESS] Cleaning complete. Total proteins removed: ", n_elim_total))
cat(paste0("\nFinal count: ", n_final_cleaned, " proteins. Files saved in '1_Cleaned_File'.\n"))

#-------------------------------------------------------------------------------
# 4. SIGNAL DRIFT CORRECTION, BIAS TEST & ADVANCED EDA
#-------------------------------------------------------------------------------
cat("\n--- 4. QC CORRECTION, BIAS TEST & ADVANCED EDA ---")

# Pull values from the parameters list
spearman_cutoff <- as.numeric(params["spearman_cutoff"])
drift_dir       <- file.path(output_main, "2_QC_correction")
eda_pre_dir     <- file.path(output_main, "3_EDA_Pre_Normalization")

# --- 4a. Data Preparation (Long Format) ---
data_long_pre_qc <- intensity_matrix_final %>%
  select(!!sym(protein_id_col), starts_with("Int_")) %>%
  pivot_longer(-!!sym(protein_id_col), names_to = "Full_ID", values_to = "Intensity") %>%
  mutate(
    CleanID = gsub("Int_", "", Full_ID),
    Condition = sub("_R[0-9]+$", "", CleanID)
  )

# --- 4b. Signal Drift Correction using QC samples ---
qc_samples <- data_long_pre_qc %>% filter(Condition == "QC")

if(nrow(qc_samples) > 0) {
  cat("\n- QC samples detected. Calculating correction factors...")
  qc_medians <- qc_samples %>% 
    group_by(CleanID) %>%
    summarise(m = median(Intensity, na.rm = TRUE), .groups = "drop")
  
  qc_target <- mean(qc_medians$m, na.rm = TRUE)
  
  corrected_data_long <- data_long_pre_qc %>%
    left_join(qc_medians %>% mutate(f = qc_target/m) %>% select(CleanID, f), by = "CleanID") %>%
    mutate(Intensity = Intensity * ifelse(is.na(f), 1, f)) %>%
    select(-f)
  
  cat(paste0("\n  [INFO] Drift correction applied. QC Target Mean: ", round(qc_target, 2)))
} else {
  cat("\n- [WARNING] No QC samples found. Skipping drift correction.")
  corrected_data_long <- data_long_pre_qc
}

# Re-pivot and attach peptide counts
corrected_data_wide <- corrected_data_long %>%
  mutate(SampleID_Col = paste0("Int_", CleanID)) %>%
  select(!!sym(protein_id_col), SampleID_Col, Intensity) %>%
  pivot_wider(names_from = SampleID_Col, values_from = Intensity)

full_corrected_matrix <- left_join(corrected_data_wide, 
                                   intensity_matrix_final %>% select(!!sym(protein_id_col), starts_with("Peps_")), 
                                   by = protein_id_col)

write_csv(full_corrected_matrix, file.path(drift_dir, "2_Corrected_Data_QC_Signal.csv"))

# --- 4c. Spearman Test (Bias Diagnosis) ---
cat("\n- Running Spearman bias test...")
bio_conds <- setdiff(unique(data_long_pre_qc$Condition), "QC")

if(length(bio_conds) >= 2) {
  c1 <- bio_conds[1]; c2 <- bio_conds[2]
  ma_prep <- corrected_data_long %>%
    filter(Condition %in% c(c1, c2)) %>%
    group_by(!!sym(protein_id_col), Condition) %>%
    summarise(m = mean(log2(Intensity + 1), na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Condition, values_from = m) %>%
    mutate(A = (get(c1) + get(c2))/2, M = get(c1) - get(c2)) %>%
    filter(is.finite(M) & is.finite(A))
  
  rho_val <- as.numeric(cor.test(ma_prep$A, ma_prep$M, method = "spearman", exact = FALSE)$estimate)
} else {
  rho_val <- 0 
}

norm_decision <- if(abs(rho_val) > spearman_cutoff) "LOESS Normalization" else "Median Normalization"
cat(paste0("\n  > Spearman Rho: ", round(rho_val, 4), " -> Recommended: ", norm_decision))

# --- 4d. EDA POST-CORRECTION (In Folder 3) ---
# Aquí es donde ocurre la magia: llamamos a la función de arriba
cat("\n- Generating Post-Correction EDA Figures (Folder 3)...")

eda_results_qc <- generate_full_eda_analysis(
  wide_data = full_corrected_matrix, 
  subfolder = "3_EDA_Pre_Normalization", 
  suffix = "PostQC", 
  title_prefix = "Post-QC Signal Correction"
)

# --- 4e. Update the Global Log ---
report_step4 <- paste0(
  "\n--- STEP 5: COMPREHENSIVE QC & EDA REPORT ---",
  "\n- Signal Drift Correction: Applied based on QC samples",
  "\n- Spearman Rho Bias:      ", round(rho_val, 4),
  "\n- Normalization Choice:   ", norm_decision,
  "\n- PERMANOVA p-value (PCA): ", ifelse(is.null(eda_results_qc), "NA", round(eda_results_qc$p_permanova_pca, 5)),
  "\n---------------------------------------------\n"
)
write(report_step4, file = log_file, append = TRUE)
cat(report_step4)

#-------------------------------------------------------------------------------
# 5. NORMALIZATION & POST-EDA REPORT
#-------------------------------------------------------------------------------
cat("\n--- 5. NORMALIZATION & POST-EDA REPORT ---")
norm_dir <- file.path(output_main, "4_Normalization")
if(!dir.exists(norm_dir)) dir.create(norm_dir, recursive = TRUE)

# --- 5a. Data Consolidation ---
data_to_norm <- full_corrected_matrix %>%
  select(all_of(protein_id_col), starts_with("Int_")) %>%
  pivot_longer(-!!sym(protein_id_col), names_to = "SampleID", values_to = "Intensity") %>%
  mutate(Condition = sub("^Int_", "", SampleID) %>% sub("_R[0-9]+$", "", .)) %>%
  filter(!is.na(Intensity) & Intensity > 0)

# --- 5b. Normalization Decision Logic (Based on Spearman Rho from Step 4) ---
if (abs(rho_val) > spearman_cutoff) {
  cat(paste0("\n- DECISION: High bias (Rho=", round(rho_val,3), "). Applying LOESS..."))
  
  wide_loess <- data_to_norm %>%
    mutate(l2 = log2(Intensity)) %>% 
    select(!!sym(protein_id_col), SampleID, l2) %>%
    pivot_wider(names_from = SampleID, values_from = l2) %>%
    column_to_rownames(protein_id_col)
  
  # Usamos limma para un loess rápido y robusto
  loess_mat <- limma::normalizeCyclicLoess(as.matrix(wide_loess), method = "fast")
  
  normalized_wide <- as.data.frame(loess_mat) %>%
    mutate(across(everything(), ~2^.x)) %>% 
    rownames_to_column(protein_id_col)
  
  norm_method <- "LOESS"
} else {
  cat(paste0("\n- DECISION: Low bias (Rho=", round(rho_val,3), "). Applying Median Normalization..."))
  
  normalized_wide <- data_to_norm %>%
    group_by(SampleID) %>%
    mutate(Median_Val = median(Intensity, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(Norm_Int = (Intensity / Median_Val) * mean(Median_Val)) %>%
    select(!!sym(protein_id_col), SampleID, Norm_Int) %>%
    pivot_wider(names_from = SampleID, values_from = Norm_Int)
  
  norm_method <- "Median"
}

# --- 5c. POST-NORMALIZATION EDA (Folder 4) ---
cat("\n- Running Post-Normalization EDA Figures...")
norm_eda_results <- generate_full_eda_analysis(
  wide_data = normalized_wide, 
  subfolder = "4_Normalization", 
  suffix = paste0("Post_", norm_method), 
  title_prefix = paste("Post-", norm_method)
)

# --- 5d. PEPTIDE RECOVERY ---
peptide_summary <- full_corrected_matrix %>%
  select(all_of(protein_id_col), starts_with("Peps_")) %>%
  pivot_longer(-!!sym(protein_id_col), values_to = "Peps") %>%
  group_by(!!sym(protein_id_col)) %>%
  summarise(Max_Unique_Peptides = max(Peps, na.rm = TRUE), .groups = "drop")

final_matrix_traceable <- normalized_wide %>%
  left_join(peptide_summary, by = protein_id_col) %>%
  select(all_of(protein_id_col), Max_Unique_Peptides, everything())

# --- 5e. FINAL TRACEABILITY REPORT ---
total_final_prots <- nrow(final_matrix_traceable)
removed_total <- n_initial_prots - total_final_prots

# We check if p-values exist for the report
p_pca_report <- if(!is.null(norm_eda_results$p_permanova_pca)) round(norm_eda_results$p_permanova_pca, 5) else "N/A"
p_pcoa_report <- if(!is.null(norm_eda_results$p_permanova_pcoa)) round(norm_eda_results$p_permanova_pcoa, 5) else "N/A"

traceability_report <- paste0(
  "\n===========================================================\n",
  "          AUTOMATED PROTEOMICS TRACEABILITY REPORT          \n",
  "===========================================================\n",
  "Project: ", exp_id, "\n",
  "Date:                           ", Sys.Date(), "\n",
  "- TOTAL PROTEINS IN RAW DATA:   ", n_initial_prots, "\n",
  "- PROTEINS REMOVED (FILTERS):   ", removed_total, "\n",
  "  (Contaminants, <", pep_min, " Peps, <", noise_cutoff*100, "% Noise)\n",
  "  (Consistency check <", min_rep_pct*100, "% Presence)\n",
  "-----------------------------------------------------------\n",
  "- FINAL PROTEINS FOR ANALYSIS:  ", total_final_prots, "\n",
  "- NORMALIZATION METHOD APPLIED: ", norm_method, "\n",
  "- POST-NORM PCA PERMANOVA P:    ", p_pca_report, "\n",
  "- POST-NORM PCoA PERMANOVA P:   ", p_pcoa_report, "\n",
  "===========================================================\n"
)

cat(traceability_report)
write(traceability_report, file = log_file, append = TRUE)
writeLines(traceability_report, file.path(norm_dir, "Traceability_Summary.txt"))

write_csv(final_matrix_traceable, file.path(norm_dir, "Final_Traceable_Matrix_Normalized.csv"))
cat("\n[SUCCESS] Step 5 complete. Data normalized and report generated.")

#-------------------------------------------------------------------------------
# 6. MISSING VALUES CLASSIFICATION & IMPUTATION (BY GROUP LOGIC)
#-------------------------------------------------------------------------------
cat("\n--- 6. MISSING VALUES CLASSIFICATION & IMPUTATION ---")
impute_dir <- file.path(output_main, "5_Missing_Values")
if(!dir.exists(impute_dir)) dir.create(impute_dir, recursive = TRUE)

# --- 6a. Preparation ---
impute_prep <- final_matrix_traceable %>%
  select(all_of(protein_id_col), starts_with("Int_")) %>%
  column_to_rownames(protein_id_col)

# Log2 transformation for stability (+1 to avoid zero issues)
impute_mat_log2 <- log2(impute_prep + 1)
impute_mat_log2 <- as.matrix(impute_mat_log2)
impute_mat_log2[is.infinite(impute_mat_log2)] <- NA

sample_names <- colnames(impute_mat_log2)
sample_groups <- sub("^Int_", "", sample_names) %>% sub("_R[0-9]+$", "", .)
group_info <- data.frame(SampleID = sample_names, Condition = sample_groups)

# --- 6b. Classification and Imputation by Group ---
cat("\n- Starting Hybrid Imputation (POV: Random Forest | MEC: QRILC) by GROUP...")

imputed_log2_final <- impute_mat_log2
total_pov_count <- 0
total_mec_count <- 0

# Loop through each condition to ensure biological context is preserved
for(cond in unique(group_info$Condition)) {
  cols_group <- group_info$SampleID[group_info$Condition == cond]
  mat_group <- imputed_log2_final[, cols_group, drop = FALSE]
  
  # Count NAs per row within this specific group
  nas_per_row <- rowSums(is.na(mat_group))
  n_reps <- length(cols_group)
  
  # 1. IDENTIFY POV (MAR): Missing in some replicates, but present in at least one
  is_pov <- nas_per_row > 0 & nas_per_row < n_reps
  total_pov_count <- total_pov_count + sum(is_pov)
  
  if(any(is_pov)) {
    cat(paste0("\n  > Group ", cond, ": Imputing ", sum(is_pov), " POV proteins (missForest)..."))
    set.seed(123)
    # RF imputation constrained to the group's own variance
    rf_res <- missForest::missForest(t(mat_group[is_pov, , drop=FALSE]), ntree = 20, maxiter = 5)
    imputed_log2_final[rownames(mat_group)[is_pov], cols_group] <- t(rf_res$ximp)
  }
  
  # 2. IDENTIFY MEC (MNAR): Missing in ALL replicates of the group
  is_mec <- nas_per_row == n_reps
  total_mec_count <- total_mec_count + sum(is_mec)
}

# 3. MEC IMPUTATION (QRILC): For proteins completely absent in a condition
cat(paste0("\n- Imputing MEC gaps with QRILC (Low-intensity tail distribution)..."))
nas_remaining_pre <- sum(is.na(imputed_log2_final))

if(nas_remaining_pre > 0) {
  set.seed(123)
  # QRILC simulates the limit of detection for proteins missing in entire groups
  qrilc_res <- imputeLCMD::impute.QRILC(imputed_log2_final, tune.sigma = 1)
  imputed_log2_final <- qrilc_res[[1]]
}

# --- 6c. Back-transformation & Recovery ---
final_imputed_intensity <- 2^imputed_log2_final - 1
final_imputed_wide <- as.data.frame(final_imputed_intensity) %>%
  rownames_to_column(protein_id_col)

final_imputed_complete <- final_imputed_wide %>%
  left_join(peptide_summary, by = protein_id_col) %>%
  select(all_of(protein_id_col), Max_Unique_Peptides, everything())

# --- 6d. POST-IMPUTATION EDA ---
# Using the master function defined in step 1h
post_imp_eda <- generate_full_eda_analysis(
  wide_data = final_imputed_complete, 
  subfolder = "5_Missing_Values", 
  suffix = "PostImputation", 
  title_prefix = "Post-Imputation (POV+MEC)"
)

# --- 6e. Final Logging (Automated Traceability Report) ---
report_impute <- paste0(
  "\n===========================================================\n",
  "          STEP 6: MISSING VALUES & IMPUTATION REPORT          \n",
  "===========================================================\n",
  "Date:                           ", Sys.Date(), "\n",
  "- POV (MAR) Proteins Imputed:     ", total_pov_count, " (Random Forest by Group)\n",
  "- MEC (MNAR) Proteins Imputed:    ", total_mec_count, " (QRILC - Low Tail)\n",
  "- Total Gaps Filled (Data Points):", nas_remaining_pre + total_pov_count, "\n",
  "-----------------------------------------------------------\n",
  "- FINAL PROTEINS FOR ANALYSIS:  ", nrow(final_imputed_complete), "\n",
  "===========================================================\n"
)

cat(report_impute)
write(report_impute, file = log_file, append = TRUE)
writeLines(report_impute, file.path(impute_dir, "Imputation_Summary.txt"))

write_csv(final_imputed_complete, file.path(impute_dir, "Final_Imputed_Matrix.csv"))
cat("\n[SUCCESS] Imputation complete. Ready for Statistical Contrast.")

#-------------------------------------------------------------------------------
# 7. DIFFERENTIAL EXPRESSION ANALYSIS (DEqMS)
#-------------------------------------------------------------------------------
cat("\n--- 7. DIFFERENTIAL EXPRESSION ANALYSIS (DEqMS) ---")

dea_dir <- file.path(output_main, "6_DEA")
if(!dir.exists(dea_dir)) dir.create(dea_dir, recursive = TRUE)

# --- 7a. Data Preparation: Universal Mapping ---
# Remove QC samples for the statistical contrast
expr_data_filtered <- final_imputed_complete %>%
  select(all_of(protein_id_col), starts_with("Int_")) %>%
  select(-contains("QC")) %>%
  column_to_rownames(protein_id_col)

# Prepare Log2 Matrix
expr_mat_log2 <- log2(as.matrix(expr_data_filtered) + 1)

# Mapping Sample Names to Universal IDs (cond1, cond2...)
current_samples <- colnames(expr_data_filtered)
clean_labels <- sub("^Int_", "", current_samples) %>% sub("_R[0-9]+$", "", .)

# Create Design Matrix based on Universal IDs from Metadata
mapping_vector <- metadata$Universal_ID
names(mapping_vector) <- metadata$Short_Label
translated_groups <- mapping_vector[clean_labels]

group_factor <- factor(translated_groups)
design <- model.matrix(~ 0 + group_factor)
colnames(design) <- levels(group_factor)

cat("\n[UNIVERSAL] Design Matrix established: ", paste(colnames(design), collapse = ", "))

# --- 7b. DEqMS Core Pipeline ---
fit <- limma::lmFit(expr_mat_log2, design)

# Get active contrasts from your matrix
experimental_contrasts_list <- contrast_matrix %>%
  filter(toupper(Active) == "YES") %>%
  mutate(Formula = paste(Treat_ID, "-", Ctrl_ID)) %>%
  select(Contrast_ID, Formula) %>%
  deframe()

all_dea_results <- list()
dea_universe <- final_imputed_complete 
sig_counts <- list()

# Record count before hypothesis contrast
prots_before_dea <- nrow(dea_universe)

for(comp_name in names(experimental_contrasts_list)) {
  current_formula <- experimental_contrasts_list[comp_name]
  cat(paste0("\n- Processing Contrast: ", comp_name, " (", current_formula, ")..."))
  
  # 1. Create Contrast Matrix
  cont_matrix <- limma::makeContrasts(contrasts = current_formula, levels = design)
  fit2 <- limma::contrasts.fit(fit, cont_matrix)
  
  # 2. eBayes with safety check for variance
  fit2$sigma[fit2$sigma == 0] <- 1e-9 
  fit2 <- limma::eBayes(fit2)
  
  # 3. DEqMS Weighting (Variance correction by Peptide Count)
  # IMPORTANT: match Max_Unique_Peptides to the fit object
  fit2$count <- final_imputed_complete$Max_Unique_Peptides[match(rownames(fit2), final_imputed_complete[[protein_id_col]])]
  fit3 <- DEqMS::spectraCounteBayes(fit2)
  
  # 4. Extract Results (sca.adj.pval is the DEqMS adjusted p-value)
  res_full <- DEqMS::outputResult(fit3, coef = 1) %>% 
    rownames_to_column(protein_id_col)
  
  all_dea_results[[comp_name]] <- res_full
  sig_counts[[comp_name]] <- sum(res_full$sca.adj.pval < 0.05, na.rm = TRUE)
  
  # --- INDIVIDUAL VOLCANO PLOT ---
  tiff(file.path(dea_dir, paste0("Volcano_", comp_name, ".tiff")), 
       width = 2000, height = 2000, res = 300, compression = "lzw")
  p_volc <- EnhancedVolcano::EnhancedVolcano(res_full, 
                                             lab = res_full[[protein_id_col]],
                                             x = 'logFC', y = 'sca.adj.pval',
                                             title = paste("Volcano:", comp_name),
                                             subtitle = paste("Formula:", current_formula),
                                             pCutoff = 0.05, FCcutoff = 0.5,
                                             legendPosition = 'bottom')
  print(p_volc)
  dev.off()
  
  # --- CONSOLIDATE MASTER DATA ---
  stats_to_merge <- res_full %>%
    select(all_of(protein_id_col), logFC, sca.P.Value, sca.adj.pval) %>%
    rename(!!paste0("logFC_", comp_name)     := logFC,
           !!paste0("P.Value_", comp_name)   := sca.P.Value,
           !!paste0("adj.P.Val_", comp_name) := sca.adj.pval)
  
  dea_universe <- left_join(dea_universe, stats_to_merge, by = protein_id_col)
  
  # Save individual result
  write_csv(res_full %>% left_join(final_imputed_complete, by = protein_id_col), 
            file.path(dea_dir, paste0("DEA_Full_Data_", comp_name, ".csv")))
  cat(" [OK]")
}

# --- 7c. GLOBAL MULTIVARIATE EDA (Post-DEA, No QC) ---
cat("\n- Generating Global Multivariate Figures (Excluding QC)...")

# Filter out QC columns for the final biological visualization
final_data_no_qc <- final_imputed_complete %>%
  select(all_of(protein_id_col), starts_with("Int_"), Max_Unique_Peptides) %>%
  select(-contains("QC"))

# Using our Master Function from 1h to see the real biological variance
generate_full_eda_analysis(
  wide_data = final_data_no_qc, 
  subfolder = "6_DEA/Global_EDA_NoQC", 
  suffix = "Final_Biological_Only", 
  title_prefix = "Final Biological Variance (No QC)"
)

# --- 7d. OVERLAP & CORRELATION ---
sig_list <- lapply(all_dea_results, function(df) as.character(df[[protein_id_col]][df$sca.adj.pval < 0.05]))
sig_list <- sig_list[sapply(sig_list, length) > 0]

if(length(sig_list) > 1) {
  cat("\n- Generating UpSet Plot for Significant Overlap...")
  tiff(file.path(dea_dir, "Upset_Significant_Overlap.tiff"), width = 3000, height = 2000, res = 300, compression = "lzw")
  print(UpSetR::upset(UpSetR::fromList(sig_list), order.by = "freq", main.bar.color = "darkred"))
  dev.off()
  
  cat("\n- Generating Log2FC Correlation Heatmap...")
  logfc_data <- dea_universe %>% select(starts_with("logFC_"))
  colnames(logfc_data) <- sub("logFC_", "", colnames(logfc_data))
  
  cor_matrix <- cor(logfc_data, use = "pairwise.complete.obs")
  tiff(file.path(dea_dir, "Log2FC_Biological_Correlation_Heatmap.tiff"), 
       width = 2400, height = 2200, res = 300, compression = "lzw")
  pheatmap::pheatmap(cor_matrix, main = "Correlation of Biological Effects (Log2FC)",
                     display_numbers = TRUE, color = colorRampPalette(c("blue", "white", "red"))(100))
  dev.off()
}

# --- 7e. Final Traceability Logging (automated report) ---
write_csv(dea_universe, file.path(dea_dir, "DEA_Master_Universe.csv"))

sig_report_text <- ""
for(n in names(sig_counts)) { 
  sig_report_text <- paste0(sig_report_text, "  - ", n, ": ", sig_counts[[n]], " significant proteins (adj.p < 0.05)\n") 
}

report_dea <- paste0(
  "\n===========================================================\n",
  "          STEP 7: DIFFERENTIAL EXPRESSION REPORT          \n",
  "===========================================================\n",
  "Date:                           ", Sys.Date(), "\n",
  "- Proteins before DEA:          ", prots_before_dea, "\n",
  "- Total contrasts analyzed:     ", length(all_dea_results), "\n",
  "- Significant results per contrast:\n", sig_report_text,
  "-----------------------------------------------------------\n"
)

cat(report_dea)
write(report_dea, file = log_file, append = TRUE)

cat("\n[SUCCESS] Statistical Analysis complete. Results in Folder 6.")

#-------------------------------------------------------------------------------
# 8. WGCNA CO-EXPRESSION ANALYSIS
#-------------------------------------------------------------------------------
cat("\n--- 8. WGCNA CO-EXPRESSION ANALYSIS ---")
wgcna_dir <- file.path(output_main, "7_CoExpression_Analysis")
if(!dir.exists(wgcna_dir)) dir.create(wgcna_dir, recursive = TRUE)

WGCNA::allowWGCNAThreads()

# --- 8a. Data Preparation & Cleaning ---
# WGCNA requires samples in rows and proteins in columns
datExpr <- dea_universe %>%
  dplyr::select(!!sym(protein_id_col), starts_with("Int_")) %>%
  dplyr::select(-contains("QC")) %>%
  tibble::column_to_rownames(protein_id_col) %>%
  t() 

# Log2 transformation
datExpr <- log2(datExpr + 1)

# Check for proteins/samples with too many missing values
gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  cat("\n[WARNING] Some proteins/samples were removed due to excessive missing values.")
}

# --- 8b. Soft Threshold Selection ---
cat("\n- Calculating Soft Threshold (Scale-Free Topology)...")
powers <- c(1:10, seq(12, 20, by = 2))
sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)

# Diagnosis Plot: Scale Free Topology
tiff(file.path(wgcna_dir, "WGCNA_ScaleFree_Topology_Diagnosis.tiff"), 
     width = 2400, height = 1200, res = 300, compression = "lzw")
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2", 
     type="n", main="Scale Independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, col="red")
abline(h=0.80, col="blue", lty=2)

plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)", 
     ylab="Mean Connectivity", type="n", main="Mean Connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, col="red")
dev.off()

softPower <- ifelse(is.na(sft$powerEstimate), 8, sft$powerEstimate)
cat(paste0("\n- Selected Soft Threshold Power: ", softPower))

# --- 8c. Network Construction ---
cat("\n- Constructing Co-expression Network (this may take a while)...")
# BACK TO YOUR ORIGINAL LOGIC: Temporary re-assignment of cor
temp_cor <- cor
cor <- WGCNA::cor 

net <- WGCNA::blockwiseModules(datExpr, 
                               power = softPower,
                               TOMType = "unsigned", 
                               minModuleSize = 20,
                               reassignThreshold = 0, 
                               mergeCutHeight = 0.25,
                               numericLabels = FALSE, 
                               pamRespectsDendro = FALSE,
                               corType = "pearson", 
                               maxBlockSize = 5000,
                               saveTOMs = FALSE, 
                               verbose = 3)

cor <- temp_cor # Restore original cor function immediately

# --- 8d. Intra-modular Connectivity (kME) ---
moduleColors <- net$colors
MEs_raw <- WGCNA::moduleEigengenes(datExpr, moduleColors)$eigengenes
kME_table <- WGCNA::signedKME(datExpr, MEs_raw)

# --- 8e. Synchronizing with Master Universe ---
wgcna_results_clean <- data.frame(
  Protein_ID_Temp = colnames(datExpr), 
  WGCNA_Module = moduleColors
)
colnames(wgcna_results_clean)[1] <- protein_id_col

dea_universe <- dea_universe %>% 
  dplyr::select(-any_of(c("WGCNA_Module", colnames(kME_table)))) %>%
  dplyr::left_join(wgcna_results_clean, by = protein_id_col) %>%
  dplyr::left_join(tibble::rownames_to_column(kME_table, protein_id_col), by = protein_id_col)

# --- 8f. Summary Table Generation ---
cat("\n- Generating Readable Summary Table...")
module_names_prefixed <- paste0("kME", dea_universe$WGCNA_Module)
kME_data <- dea_universe %>% dplyr::select(starts_with("kME"))

assigned_kME <- sapply(1:nrow(dea_universe), function(i) {
  col_name <- module_names_prefixed[i]
  if (col_name %in% colnames(kME_data)) {
    return(kME_data[i, col_name])
  } else {
    return(NA)
  }
})

wgcna_summary <- data.frame(
  Protein_ID = dea_universe[[protein_id_col]],
  Module_Color = dea_universe$WGCNA_Module,
  kME_Membership = assigned_kME
) %>% 
  dplyr::filter(!is.na(Module_Color)) %>%
  dplyr::arrange(Module_Color, desc(kME_Membership))

write_csv(wgcna_summary, file.path(wgcna_dir, "WGCNA_Module_Assignments_Summary.csv"))

# --- 8g. Cluster Dendrogram Plot ---
tiff(file.path(wgcna_dir, "WGCNA_Cluster_Dendrogram.tiff"), 
     width = 2400, height = 1800, res = 300, compression = "lzw")
WGCNA::plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]], 
                           "Module Color", dendroLabels = FALSE, hang = 0.03, 
                           addGuide = TRUE, guideHang = 0.05)
dev.off()

# --- 8h. Logging for Automatic Report ---
wgcna_assigned <- sum(moduleColors != "grey")
wgcna_total <- length(moduleColors)

report_wgcna <- paste0(
  "\n===========================================================\n",
  "          STEP 8: WGCNA CO-EXPRESSION REPORT          \n",
  "===========================================================\n",
  "Date:                           ", Sys.Date(), "\n",
  "- Total proteins analyzed:      ", wgcna_total, "\n",
  "- Modules identified:           ", length(unique(moduleColors)), "\n",
  "- Proteins assigned to modules: ", wgcna_assigned, " (", round(wgcna_assigned/wgcna_total*100, 2), "%)\n",
  "- Unassigned proteins (grey):   ", sum(moduleColors == "grey"), "\n",
  "-----------------------------------------------------------\n"
)

cat(report_wgcna)
write(report_wgcna, file = log_file, append = TRUE)

saveRDS(dea_universe, file.path(output_main, "6_DEA", "Master_Universe_DEA_WGCNA.rds"))
cat("\n[SUCCESS] WGCNA complete. Master Matrix and Visualizations generated.")

#-------------------------------------------------------------------------------
# 9. INTEGRATED ANALYSIS: DEqMS & WGCNA FUSION (Param-Driven)
#-------------------------------------------------------------------------------
cat("\n--- 9. STARTING INTEGRATED ANALYSIS AND PRECISION FILTERING ---")

# 1. Load thresholds from parameters
p_cut_param   <- as.numeric(params["adj_pvalue_cutoff"])
kME_cut_param <- as.numeric(params["wgcna_kME_threshold"])

output_8 <- file.path(output_main, "8_EDA_and_CoExpression")
if(!dir.exists(output_8)) dir.create(output_8, recursive = TRUE)

# --- 9.1. GLOBAL SCALE SYNCHRONIZATION (The Fix) ---
# We transform the intensities to Log2 BEFORE fusion to ensure 
# Stan (Bayesian) and DEqMS talk the same language.
cat("\n- Synchronizing Intensity Scales (Ensuring Log2 for Bayesian/Downstream)...")

int_cols_master <- grep("^Int_", colnames(final_imputed_complete), value = TRUE)

# Apply Log2(x+1) only if the data is currently in linear scale (max > 50)
if(max(final_imputed_complete[, int_cols_master], na.rm = TRUE) > 50) {
  final_imputed_complete <- final_imputed_complete %>%
    mutate(across(all_of(int_cols_master), ~ log2(.x + 1)))
  cat("\n  [OK] Linear intensities detected and converted to Log2.")
} else {
  cat("\n  [SKIP] Intensities already appear to be in Log2 scale.")
}

# --- 9.2. DATA FUSION (Master Universe) ---
cat("\n- Synchronizing Master Table with Param-Driven Thresholds...")

master_integrated <- dea_universe %>%
  # Update intensities in the universe from the newly logged master table
  select(-all_of(int_cols_master)) %>%
  left_join(final_imputed_complete %>% select(all_of(protein_id_col), all_of(int_cols_master)), 
            by = protein_id_col) %>%
  mutate(
    Is_Significant_DEA = if_any(matches("adj\\.P\\.Val_"), ~ . <= p_cut_param),
    Is_Functional_WGCNA = !is.na(WGCNA_Module) & WGCNA_Module != "grey",
    Category = case_when(
      Is_Significant_DEA & Is_Functional_WGCNA  ~ "Both",
      Is_Significant_DEA & !Is_Functional_WGCNA ~ "DEA Only",
      !Is_Significant_DEA & Is_Functional_WGCNA ~ "WGCNA Only",
      TRUE                                      ~ "None"
    )
  )

# --- 9.3. FILTERING CANDIDATES ---
res_master_filtered <- master_integrated %>%
  filter(Category != "None")

cat(paste0("\n- Filtering complete (P-adj <= ", p_cut_param, "):"))
cat(paste0("\n  > Total Proteins in Universe: ", nrow(master_integrated)))
cat(paste0("\n  > Candidates for Bayesian: ", nrow(res_master_filtered)))

# --- 9.4. GLOBAL MULTIVARIATE ANALYSIS (Master Function 1h) ---
cat("\n- Generating Global Multivariate Analysis (PCA/PCoA/PERMANOVA)...")

data_for_eda <- final_imputed_complete %>%
  select(all_of(protein_id_col), starts_with("Int_")) %>%
  select(-contains("QC"))

global_multivar_results <- generate_full_eda_analysis(
  wide_data = data_for_eda, 
  subfolder = "8_EDA_and_CoExpression", 
  suffix = "Global_Integration", 
  title_prefix = "Global Integrated Dataset (No QC)"
)

# --- 9.5. GENERATION OF CONTRAST-SPECIFIC TABLES ---
cat("\n- Generating precision-filtered tables per contrast...")
output_8_contrasts <- file.path(output_8, "Contrast_Specific_Tables")
if(!dir.exists(output_8_contrasts)) dir.create(output_8_contrasts, recursive = TRUE)

active_contrasts <- contrast_matrix %>% filter(toupper(Active) == "YES") %>% pull(Contrast_ID)

for(cont_name in active_contrasts) {
  adj_p_col <- paste0("adj.P.Val_", cont_name)
  logfc_col <- paste0("logFC_", cont_name)
  
  tab_cont <- res_master_filtered %>%
    filter((!!sym(adj_p_col) <= p_cut_param) | Is_Functional_WGCNA) %>%
    dplyr::select(
      all_of(protein_id_col), 
      WGCNA_Module, 
      Category,
      !!sym(logfc_col), 
      !!sym(adj_p_col),
      starts_with("Int_")
    ) %>%
    mutate(Comparison = cont_name)
  
  write_csv(tab_cont, file.path(output_8_contrasts, paste0("Table_Integrated_", cont_name, ".csv")))
}

# --- 9.6. DATA VISUALIZATION ---
cat("\n- Generating integration-specific plots...")

# A. Venn Diagram
dea_prots <- master_integrated %>% filter(Is_Significant_DEA) %>% pull(!!sym(protein_id_col))
wgcna_prots <- master_integrated %>% filter(Is_Functional_WGCNA) %>% pull(!!sym(protein_id_col))
venn_list <- list(DEqMS_Significant = dea_prots, WGCNA_Functional = wgcna_prots)

tiff(file.path(output_8, "Venn_Integration_Methods.tiff"), width = 2000, height = 2000, res = 300, compression = "lzw")
grid::grid.newpage()
venn_plot <- VennDiagram::venn.diagram(
  x = venn_list, filename = NULL, fill = c("steelblue", "orange"), alpha = 0.5, main = "DEA vs WGCNA Integration"
)
grid::grid.draw(venn_plot)
dev.off()

# B. Module Membership
mod_stats <- res_master_filtered %>% group_by(WGCNA_Module, Category) %>% tally()
p_mod_dist <- ggplot(mod_stats, aes(x = reorder(WGCNA_Module, -n), y = n, fill = Category)) +
  geom_bar(stat = "identity") + theme_minimal() + scale_fill_brewer(palette = "Set2") +
  labs(title = "Protein Distribution per Module", x = "WGCNA Module", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_8, "Module_Significance_Distribution.tiff"), p_mod_dist, width = 10, height = 6, dpi = 300)

# --- 9.7. PERSISTENCE ---
write_csv(res_master_filtered, file.path(output_8, "Master_Integrated_Filtered_for_Bayesian.csv"))
saveRDS(res_master_filtered, file.path(output_8, "Master_Integrated_Filtered_for_Bayesian.rds"))

# --- 9.8. AUTOMATED REPORTING (automated precision report) ---
p_pca_val  <- as.numeric(global_multivar_results$p_permanova_pca)
p_pcoa_val <- as.numeric(global_multivar_results$p_permanova_pcoa)

report_integration <- paste0(
  "\n===========================================================\n",
  "          STEP 9: INTEGRATED ANALYSIS REPORT          \n",
  "===========================================================\n",
  "Date:                           ", Sys.Date(), "\n",
  "- Scale:                        Log2 Normalized (Transformed at Step 9)\n",
  "- Total Proteins in Universe:   ", nrow(master_integrated), "\n",
  "- Candidate Proteins (Selected):", nrow(res_master_filtered), "\n",
  "  > DEA Only:        ", sum(res_master_filtered$Category == "DEA Only"), "\n",
  "  > WGCNA Only:      ", sum(res_master_filtered$Category == "WGCNA Only"), "\n",
  "  > Both Methods:    ", sum(res_master_filtered$Category == "Both"), "\n",
  "-----------------------------------------------------------\n",
  "- Global PCA PERMANOVA P:        ", ifelse(is.na(p_pca_val), "N/A", round(p_pca_val, 5)), "\n",
  "- Global PCoA PERMANOVA P:       ", ifelse(is.na(p_pcoa_val), "N/A", round(p_pcoa_val, 5)), "\n",
  "===========================================================\n"
)

cat(report_integration)
write(report_integration, file = log_file, append = TRUE)

cat("\n[SUCCESS] Integrated Analysis complete. Candidate list is now SCALE-SYNCHRONIZED.")

#-------------------------------------------------------------------------------
# 10 - STEP 1: SYSTEMATIC EXTRACTION (DEA + WGCNA HUB TRACEABILITY)
#-------------------------------------------------------------------------------
cat("\n--- 10.1: GENERATING TABLES WITH PARAMETER-BASED TRACEABILITY ---")

kME_threshold_val   <- as.numeric(params["wgcna_kME_threshold"])
adj_p_threshold_val <- as.numeric(params["adj_pvalue_cutoff"])
protein_id_col_val  <- protein_id_col 

output_9 <- file.path(output_main, "9_CBR_Filter")
if(!dir.exists(output_9)) dir.create(output_9, recursive = TRUE)

active_contrasts <- contrast_matrix %>% dplyr::filter(toupper(Active) == "YES")

# A vector is created to save the paths in paths
contrast_folders <- c()

for (i in 1:nrow(active_contrasts)) {
  current_contrast <- active_contrasts$Contrast_ID[i]
  treat_label <- translate_condition(active_contrasts$Treat_ID[i], metadata)
  ctrl_label  <- translate_condition(active_contrasts$Ctrl_ID[i], metadata)
  
  contrast_dir <- file.path(output_9, current_contrast)
  if(!dir.exists(contrast_dir)) dir.create(contrast_dir, recursive = TRUE)
  contrast_folders <- c(contrast_folders, contrast_dir)
  
  pval_col  <- paste0("adj.P.Val_", current_contrast)
  logfc_col <- paste0("logFC_", current_contrast)
  
  t_reps <- colnames(res_master_filtered)[grepl(paste0("^Int_", treat_label, "_R"), colnames(res_master_filtered))]
  c_reps <- colnames(res_master_filtered)[grepl(paste0("^Int_", ctrl_label, "_R"), colnames(res_master_filtered))]
  
  cat(paste0("\n>>> Processing Contrast [", i, "/", nrow(active_contrasts), "]: ", current_contrast))
  
  # kME_Value is calculate, if it does not exist, to filter
  temp_subset <- res_master_filtered 
  temp_subset$kME_Value <- sapply(1:nrow(temp_subset), function(idx) {
    mod <- temp_subset$WGCNA_Module[idx]
    if(is.na(mod) || mod == "grey") return(0)
    col_name <- paste0("kME", mod)
    if(col_name %in% colnames(temp_subset)) return(as.numeric(temp_subset[idx, col_name])) else return(0)
  })
  
  filtered_subset <- temp_subset %>%
    dplyr::filter((!!sym(pval_col) <= adj_p_threshold_val) | (kME_Value >= kME_threshold_val)) %>%
    dplyr::select(all_of(protein_id_col_val), !!sym(logfc_col), !!sym(pval_col), 
                  all_of(t_reps), all_of(c_reps), WGCNA_Module, kME_Value)
  
  if(nrow(filtered_subset) > 0) {
    output_path <- file.path(contrast_dir, paste0("Input_Bayesian_", current_contrast, ".csv"))
    readr::write_csv(filtered_subset, output_path)
    cat(paste0(" -> Candidates: ", nrow(filtered_subset)))
  }
}

#-------------------------------------------------------------------------------
# 10 - STEP 2: STABLE BAYESIAN MCMC
#-------------------------------------------------------------------------------
cat("\n\n--- 10.2: STARTING PARAMETER-DRIVEN BAYESIAN INFERENCE ---")

cbr_fc_target <- as.numeric(params["log2fc_cutoff"]) 

# Función MCMC pura (sin detecciones raras, ya viene en Log2 del paso 9)
run_bayesian_mcmc_stable <- function(y_treat, y_ctrl, fc_target) {
  dat_b <- data.frame(
    val = c(y_ctrl, y_treat), 
    grp = factor(c(rep("Ctrl", length(y_ctrl)), rep("Treat", length(y_treat))))
  )
  dat_b$grp <- relevel(dat_b$grp, ref = "Ctrl")
  
  if(var(y_ctrl, na.rm=T) == 0) dat_b$val[1:length(y_ctrl)] <- dat_b$val[1:length(y_ctrl)] + rnorm(length(y_ctrl), 0, 1e-5)
  if(var(y_treat, na.rm=T) == 0) dat_b$val[(length(y_ctrl)+1):nrow(dat_b)] <- dat_b$val[(length(y_ctrl)+1):nrow(dat_b)] + rnorm(length(y_treat), 0, 1e-5)
  
  fit <- tryCatch({
    rstanarm::stan_glm(val ~ grp, data = dat_b, family = gaussian(),
                       chains = as.numeric(params["mcmc_chains"]), 
                       iter = as.numeric(params["mcmc_iter"]), 
                       refresh = 0, seed = 123)
  }, error = function(e) { return(NULL) })
  
  if(is.null(fit)) return(list(prob = NA, fc = NA, status = "Fit_Error"))
  
  post <- as.matrix(fit)
  mcmc_fc <- mean(post[, "grpTreat"])
  p_rel <- if(mcmc_fc > 0) mean(post[, "grpTreat"] > fc_target) else mean(post[, "grpTreat"] < -fc_target)
  return(list(prob = p_rel, fc = mcmc_fc, status = "Success"))
}

for (f_path in contrast_folders) {
  contrast_id <- basename(f_path)
  input_file  <- file.path(f_path, paste0("Input_Bayesian_", contrast_id, ".csv"))
  output_file <- file.path(f_path, paste0("Full_Bayesian_Results_", contrast_id, ".csv"))
  if(!file.exists(input_file)) next
  
  cat(paste0("\n>>> MCMC Inference: ", contrast_id))
  input_data <- readr::read_csv(input_file, show_col_types = FALSE)
  parts <- unlist(strsplit(contrast_id, "_vs_"))
  t_cols <- grep(paste0("Int_", parts[1], "_"), colnames(input_data), value = TRUE)
  c_cols <- grep(paste0("Int_", parts[2], "_"), colnames(input_data), value = TRUE)
  
  if(file.exists(output_file)) file.remove(output_file)
  
  pb <- txtProgressBar(min = 0, max = nrow(input_data), style = 3)
  for(j in 1:nrow(input_data)) {
    res <- run_bayesian_mcmc_stable(as.numeric(input_data[j, t_cols]), as.numeric(input_data[j, c_cols]), cbr_fc_target)
    res_row <- input_data[j, ] %>% dplyr::mutate(prob_biol_relevance = res$prob, log2FC_bayes = res$fc, status_bayesian = res$status)
    readr::write_csv(res_row, output_file, append = file.exists(output_file))
    setTxtProgressBar(pb, j)
  }
  close(pb)
}

#-------------------------------------------------------------------------------
# 10.3: ELITE FILTERING, MULTI-LEVEL EDA (GLOBAL & LOCAL) & VISUALS
#-------------------------------------------------------------------------------
cat("\n\n--- 10.3: GENERATING ELITE RESULTS AND FULL EDA SUITE ---")

prob_threshold <- as.numeric(params["prob_rel_cutoff"])
desc_clean     <- readr::read_tsv(path_desc, show_col_types = FALSE) %>%
  dplyr::select(Entry = 1, Protein_Names = 3, EC_Number = 6, Subcellular_Location = 7) %>%
  dplyr::distinct(Entry, .keep_all = TRUE)

# --- A. GLOBAL EDA (Project Level: All Samples, All Candidates) ---
# This generates PCA, PCoA, PERMANOVA and Correlation for the entire project
cat("\n- Executing Global EDA for all Bayesian candidates (Project Overview)...")
generate_full_eda_analysis(
  wide_data = res_master_filtered, 
  subfolder = "9_CBR_Filter", 
  suffix = "Global_Candidates_Project", 
  title_prefix = "Project Overview: All Candidates"
)

# --- B. PER-CONTRAST LOOP (Local Level: Contrast-Specific Elite Proteins) ---
elite_counts <- list()

for (f_path in contrast_folders) {
  contrast_id <- basename(f_path)
  results_file <- file.path(f_path, paste0("Full_Bayesian_Results_", contrast_id, ".csv"))
  if(!file.exists(results_file)) next
  
  all_res <- readr::read_csv(results_file, show_col_types = FALSE)
  
  # 1. Biological Relevance Filtering (Elite Proteins)
  res_elite <- all_res %>%
    dplyr::left_join(desc_clean, by = c("PG.ProteinAccessions" = "Entry")) %>%
    dplyr::filter(prob_biol_relevance >= prob_threshold) %>%
    dplyr::arrange(desc(prob_biol_relevance))
  
  elite_counts[[contrast_id]] <- nrow(res_elite)
  readr::write_excel_csv(res_elite, file.path(f_path, paste0("Final_Elite_Proteins_", contrast_id, ".csv")))
  
  # 2. LOCAL EDA SUITE (Full 1h Suite for this specific comparison)
  # Generates PCA/PCoA/PERMANOVA/Correlation ONLY for this contrast's Elite proteins
  if(nrow(res_elite) >= 5) {
    cat(paste0("\n- Generating Local EDA Suite for contrast: ", contrast_id))
    generate_full_eda_analysis(
      wide_data = res_elite, 
      subfolder = file.path("9_CBR_Filter", contrast_id), 
      suffix = paste0("Elite_Local_", contrast_id), 
      title_prefix = paste0("Contrast Elite: ", contrast_id)
    )
    
    # 3. SPECIFIC ELITE HEATMAP (Clustered view of the comparison)
    parts  <- unlist(strsplit(contrast_id, "_vs_"))
    t_cols <- grep(paste0("Int_", parts[1], "_"), colnames(res_elite), value = TRUE)
    c_cols <- grep(paste0("Int_", parts[2], "_"), colnames(res_elite), value = TRUE)
    plot_mat <- as.matrix(res_elite[, c(t_cols, c_cols)])
    rownames(plot_mat) <- res_elite[[protein_id_col]]
    
    tiff(file.path(f_path, paste0("1_Heatmap_Elite_Only_", contrast_id, ".tiff")), width = 2400, height = 3000, res = 300)
    pheatmap::pheatmap(plot_mat, scale = "row", main = paste("Elite Proteins Log2 (n =", nrow(plot_mat), ")"),
                       color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
    dev.off()
  }
  
  # 4. BAYESIAN VOLCANO PLOT
  p_volc <- ggplot(all_res, aes(x = log2FC_bayes, y = prob_biol_relevance)) +
    geom_point(aes(color = prob_biol_relevance >= prob_threshold), alpha = 0.5) +
    scale_color_manual(values = c("grey70", "firebrick")) +
    geom_hline(yintercept = prob_threshold, linetype = "dashed", color = "blue") +
    theme_minimal() + 
    labs(title = paste("Bayesian Volcano Plot:", contrast_id), 
         subtitle = paste("Prob. Threshold:", prob_threshold),
         x = "Log2 Fold Change (Bayesian)", 
         y = "Probability of Biological Relevance")
  
  ggsave(file.path(f_path, paste0("3_Volcano_Bayesian_", contrast_id, ".tiff")), p_volc, width = 8, height = 6, dpi = 300)
}

# --- C. AUTOMATED TRACEABILITY REPORT ---
elite_report_text <- ""
for(n in names(elite_counts)) { 
  elite_report_text <- paste0(elite_report_text, "  - ", n, ": ", elite_counts[[n]], " Elite Proteins\n") 
}

report_final <- paste0(
  "\n===========================================================\n",
  "          FINAL PROTEOMICS ELITE TRACEABILITY          \n",
  "===========================================================\n",
  "Date:                           ", Sys.Date(), "\n",
  "- Bayesian Prob. Threshold:      ", prob_threshold, "\n",
  "- Log2FC Target (CBR):          ", as.numeric(params["log2fc_cutoff"]), "\n",
  "-----------------------------------------------------------\n",
  "ELITE PROTEINS AFTER BIOLOGICAL RELEVANCE FILTER:\n",
  elite_report_text,
  "===========================================================\n"
)

cat(report_final)
write(report_final, file = log_file, append = TRUE)

cat("\n[SUCCESS] Universal Step 10 complete. Final Elite results and plots ready.")

#-------------------------------------------------------------------------------
# 11. UNIVERSAL FASTA GENERATION (WITH EXPERIMENT IDENTITY)
#-------------------------------------------------------------------------------
cat("\n--- PHASE 11: PREPARING FASTA FOR EXTERNAL ANNOTATION ---")

# 1. Filename with Experiment ID
# We use exp_id (e.g., "proteomics_analysis") so the FASTA name stays consistent
fasta_filename <- paste0(exp_id, ".fasta")

# 2. SMART SEARCH FOR THE UNIVERSE FILE
cat("\n- Searching for DEA_Master_Universe.csv...")
universe_file_path <- list.files(path = output_main, pattern = "DEA_Master_Universe.csv", recursive = TRUE, full.names = TRUE)

if(length(universe_file_path) == 0) stop("Critical Error: DEA_Master_Universe.csv not found.")

# 3. LOAD & FILTER
dea_master <- readr::read_csv(universe_file_path[1], show_col_types = FALSE)
unique_ids <- unique(dea_master[[1]])
unique_ids <- unique_ids[!is.na(unique_ids)]

# 4. SUBSET ORIGINAL FASTA
# path_fasta should be your original large database file
full_fasta <- Biostrings::readAAStringSet(path_fasta)
fasta_ids <- sapply(strsplit(names(full_fasta), "\\|"), function(x) if(length(x) >= 2) x[2] else x[1])
universe_fasta <- full_fasta[fasta_ids %in% unique_ids]

# 5. SAVE WITH UNIQUE NAME
path_universe_fasta <- file.path(output_main, fasta_filename)
Biostrings::writeXStringSet(universe_fasta, path_universe_fasta)

cat(paste0("\n[SUCCESS] Universal FASTA saved: ", path_universe_fasta))

#-------------------------------------------------------------------------------
# 12. WSL BRIDGE (DYNAMIC PROJECT SYNC)
#-------------------------------------------------------------------------------
cat("\n--- PHASE 12: WSL BRIDGE ---")

# Identify file and folder names from global variables
current_fasta_name <- basename(path_universe_fasta) # e.g., proteomics_analysis.fasta
exp_folder_name    <- basename(output_main)        # e.g., proteomics_analysis_Results

# 1. DETECT WSL HOME
wsl_home_win <- system('wsl.exe wslpath -w ~', intern = TRUE)

# 2. CREATE PROJECT FOLDER IN WSL
# This creates ~/FASTA_Annotations/<exp_folder_name>/ inside WSL
system(paste0('wsl.exe mkdir -p ~/FASTA_Annotations/', exp_folder_name))

# 3. TRANSFER FASTA TO WSL
from_file_win <- normalizePath(path_universe_fasta)
copy_cmd <- paste0('wsl.exe cp "$(wslpath \'', from_file_win, '\')" ~/FASTA_Annotations/', exp_folder_name, '/')
system(copy_cmd)

# 4. VERIFICATION
check_cmd <- paste0('wsl.exe ls ~/FASTA_Annotations/', exp_folder_name, '/', current_fasta_name)
check_result <- system(check_cmd, intern = TRUE, ignore.stderr = TRUE)

if(length(check_result) > 0) {
  cat(paste0("\n[OK] Bridged: ", current_fasta_name))
} else {
  stop("\n[ERROR] Transfer failed for: ", current_fasta_name)
}

# 5. EXECUTE ANNOTATOR
# Passing both: Experiment Folder ($1) AND Dynamic Cores ($2)
system(paste0('wsl.exe bash ~/ElAnotaor.sh ', exp_folder_name, ' ', n_cores_final))

# 6. RETRIEVE ANNOTATION RESULTS
# Results are pulled back into a dedicated subfolder in Windows
win_results_dest <- file.path(output_main, "13_Annotation_Results")
if(!dir.exists(win_results_dest)) dir.create(win_results_dest, recursive = TRUE)

return_cmd <- paste0('wsl.exe cp -r ~/FASTA_Annotations/', exp_folder_name, '/5_Annotation_results/* ',
                     ' "$(wslpath \'', win_results_dest, '\')"')
system(return_cmd)

cat("\n[SUCCESS] WSL Processing complete for Experiment: ", exp_id)

  #-------------------------------------------------------------------------------
  # 13. FUNCTIONAL ENRICHMENT MODULE (ORA & GSEA)
  #-------------------------------------------------------------------------------
  cat("\n--- STARTING FUNCTIONAL ENRICHMENT MODULE ---")
  
  # 13a. LOAD ANNOTATION DATA
  # ------------------------------------------------------------------------------
  path_egg_res <- list.files(win_results_dest, pattern = "\\.emapper\\.annotations$", full.names = TRUE)
  path_ips_res <- list.files(win_results_dest, pattern = "InterProScan.*\\.tsv$", full.names = TRUE)
  
  egg_df_cleaned <- readr::read_tsv(path_egg_res, skip = 4, show_col_types = FALSE, comment = "##")
  ips_df_cleaned <- readr::read_tsv(path_ips_res, col_names = FALSE, show_col_types = FALSE)
  
  # 13b. ROBUST ID NORMALIZATION
  # ------------------------------------------------------------------------------
  normalize_ids <- function(id_vec) {
    ifelse(grepl("\\|", id_vec), 
           gsub("^.*\\|(.*)\\|.*$", "\\1", id_vec), 
           id_vec)
  }
  
  egg_df_cleaned[[1]] <- normalize_ids(egg_df_cleaned[[1]])
  ips_df_cleaned[[1]] <- normalize_ids(ips_df_cleaned[[1]])
  
  # 13c. ANNOTATION SET PREPARATION
  # ------------------------------------------------------------------------------
  process_go_fixed <- function(df, col_idx, root_go = NULL) {
    df %>%
      dplyr::select(Term = !!col_idx, gene = 1) %>%
      dplyr::filter(!is.na(Term) & Term != "-" & Term != "") %>%
      dplyr::mutate(Term = gsub("\\([^)]+\\)", "", Term)) %>%
      tidyr::separate_rows(Term, sep = "[\\|,]") %>%
      dplyr::mutate(Term = trimws(Term)) %>%
      dplyr::filter(grepl("^GO:[0-9]{7}$", Term)) %>%
      dplyr::distinct()
  }
  
  # 13d. PREPARING IPS SETS
  cat("\n- Preparing InterProScan sets...")
  ips_go_bp <- process_go_fixed(ips_df_cleaned, 14)
  ips_go_mf <- process_go_fixed(ips_df_cleaned, 14)
  ips_go_cc <- process_go_fixed(ips_df_cleaned, 14)
  
  ips_paths_raw <- ips_df_cleaned %>%
    dplyr::select(gene = 1, Path = 15) %>%
    dplyr::filter(!is.na(Path) & Path != "-") %>%
    tidyr::separate_rows(Path, sep = "\\|") %>%
    dplyr::mutate(Path = trimws(Path))
  
  ips_reactome <- ips_paths_raw %>% dplyr::filter(grepl("^Reactome:", Path)) %>% dplyr::select(Term = Path, gene)
  ips_metacyc  <- ips_paths_raw %>% dplyr::filter(grepl("^MetaCyc:", Path)) %>% dplyr::select(Term = Path, gene)
  
  ips_signatures <- ips_df_cleaned %>%
    dplyr::select(gene = 1, Source = 4, Signature = 5) %>%
    dplyr::filter(!is.na(Signature) & Signature != "-") %>%
    dplyr::mutate(Term = paste0(Source, ":", Signature)) %>%
    dplyr::select(Term, gene) %>%
    dplyr::distinct()
  
  # 13e. PREPARING EGGNOG SETS
  cat("\n- Preparing EggNOG sets...")
  egg_go_all <- egg_df_cleaned %>%
    dplyr::select(Term = 10, gene = 1) %>%
    dplyr::filter(!is.na(Term) & Term != "-" & Term != "") %>%
    tidyr::separate_rows(Term, sep = ",") %>%
    dplyr::mutate(Term = trimws(Term)) %>%
    dplyr::filter(grepl("^GO:[0-9]{7}$", Term))
  
  egg_kegg <- egg_df_cleaned %>%
    dplyr::select(Term = 12, gene = 1) %>%
    dplyr::filter(!is.na(Term) & Term != "-" & Term != "") %>%
    tidyr::separate_rows(Term, sep = ",") %>%
    dplyr::mutate(Term = gsub("ko:", "", Term), Term = trimws(Term))
  
  egg_cazy <- egg_df_cleaned %>%
    dplyr::select(Term = 19, gene = 1) %>%
    dplyr::filter(!is.na(Term) & Term != "-" & Term != "") %>%
    tidyr::separate_rows(Term, sep = ",") %>%
    dplyr::mutate(Term = trimws(Term))
  
  egg_pfam <- egg_df_cleaned %>%
    dplyr::select(Term = 21, gene = 1) %>%
    dplyr::filter(!is.na(Term) & Term != "-" & Term != "") %>%
    tidyr::separate_rows(Term, sep = ",") %>%
    dplyr::mutate(Term = trimws(Term))
  
  # 13f. CONSOLIDATED ANALYSIS LIST
  all_analyses <- list(
    "IPS_GO_BP"    = ips_go_bp, "IPS_GO_MF" = ips_go_mf, "IPS_GO_CC" = ips_go_cc,
    "IPS_Reactome" = ips_reactome, "IPS_MetaCyc" = ips_metacyc, "IPS_Signatures" = ips_signatures,
    "EGG_GO_All"   = egg_go_all, "EGG_KEGG" = egg_kegg, "EGG_CAZy" = egg_cazy, "EGG_PFAM" = egg_pfam
  )
  
  # 13g. OVER-REPRESENTATION ANALYSIS (ORA)
  # ------------------------------------------------------------------------------
  cat("\n--- EXECUTING ORA & SAVING DIRECTIONAL PROTEIN LISTS ---")
  
  elite_files <- list.files(file.path(output_main, "9_CBR_Filter"), 
                            pattern = "Final_Elite_Proteins_.*\\.csv$", 
                            recursive = TRUE, full.names = TRUE)
  
  master_univ_path <- file.path(output_main, "6_DEA", "DEA_Master_Universe.csv")
  master_univ <- read.csv(master_univ_path)
  
  # Universal Column Detection for Master Universe
  id_col_master <- grep("ProteinAccessions|ProteinID|Accession|Protein", colnames(master_univ), value = TRUE)[1]
  master_univ[[id_col_master]] <- normalize_ids(master_univ[[id_col_master]])
  all_quant_prot <- unique(master_univ[[id_col_master]])
  
  for (f_path in elite_files) {
    contrast_name <- basename(dirname(f_path))
    cat(paste0("\n\n> Processing Contrast: ", contrast_name))
    
    df_elite <- read.csv(f_path)
    
    # Universal Column Detection for Elite Files
    id_col_elite <- grep("ProteinAccessions|ProteinID|Accession|Protein", colnames(df_elite), value = TRUE)[1]
    lfc_col <- grep("log2FC_bayes|Log2FC_Bayes|logFC|Log2FC", colnames(df_elite), value = TRUE)[1]
    
    df_elite[[id_col_elite]] <- normalize_ids(df_elite[[id_col_elite]])
    
    dir_ora_contrast <- file.path(output_main, "10_ORA", contrast_name)
    if(!dir.exists(dir_ora_contrast)) dir.create(dir_ora_contrast, recursive = TRUE)
    
    for (direction in c("UP", "DOWN")) {
      subset_elite <- if(direction == "UP") df_elite[df_elite[[lfc_col]] > 0, ] else df_elite[df_elite[[lfc_col]] < 0, ]
      if(nrow(subset_elite) == 0) next
      
      write.csv(subset_elite, file.path(dir_ora_contrast, paste0("Target_Proteins_", direction, ".csv")), row.names = FALSE)
      
      target_ids <- subset_elite[[id_col_elite]]
      cat(paste0("\n  - ", direction, ": ", length(target_ids), " proteins exported."))
      
      for (db_name in names(all_analyses)) {
        t2g <- all_analyses[[db_name]]
        univ_source <- intersect(all_quant_prot, unique(t2g$gene))
        genes_to_test <- intersect(target_ids, t2g$gene)
        
        if (length(genes_to_test) < 2) next
        
        res_ora <- tryCatch({
          clusterProfiler::enricher(gene = genes_to_test, universe = univ_source, TERM2GENE = t2g, pvalueCutoff = 0.05)
        }, error = function(e) NULL)
        
        if (!is.null(res_ora) && nrow(as.data.frame(res_ora)) > 0) {
          fname <- paste0("ORA_", db_name, "_", direction)
          write.csv(as.data.frame(res_ora), file.path(dir_ora_contrast, paste0(fname, ".csv")), row.names = FALSE)
          p <- enrichplot::dotplot(res_ora, showCategory = 15) + ggplot2::ggtitle(paste(contrast_name, "-", fname)) + ggplot2::theme_bw()
          ggplot2::ggsave(file.path(dir_ora_contrast, paste0(fname, ".tiff")), plot = p, width = 10, height = 8, dpi = 300)
        }
      }
    }
  }
  
  # ==============================================================================
  # 13h. GSEA EXECUTION (WITH PATH TRACEABILITY & ANNOTATION)
  # =================-------------------------------------------------------------
  cat("\n--- STARTING GSEA MODULE ---")
  
  # SAFETY: Confirm the base output directory
  # normalizePath converts relative paths to absolute Windows paths
  base_output <- normalizePath(output_main, mustWork = FALSE)
  cat("\n[CHECK] Analysis will be saved to: ", base_output)
  
  # 1. Preparation of Functional Dictionary for Annotation
  dict_extra <- egg_df_cleaned %>% 
    dplyr::select(ID = 12, Desc = 13) %>% 
    dplyr::bind_rows(egg_df_cleaned %>% dplyr::select(ID = 21, Desc = 21)) %>% 
    dplyr::filter(!is.na(ID) & ID != "-") %>%
    tidyr::separate_rows(ID, sep = ",") %>% 
    dplyr::mutate(ID = gsub("ko:", "", ID)) %>% 
    dplyr::distinct(ID, .keep_all = TRUE)
  
  logfc_cols <- grep("^logFC_", names(master_univ), value = TRUE)
  
  for (lfc_col in logfc_cols) {
    contrast_id <- gsub("logFC_", "", lfc_col)
    
    # Define and normalize the contrast folder path
    dir_gsea_contrast <- normalizePath(file.path(output_main, "11_GSEA", contrast_id), mustWork = FALSE)
    
    if(!dir.exists(dir_gsea_contrast)) {
      dir.create(dir_gsea_contrast, recursive = TRUE)
      cat(paste0("\n\n[FOLDER CREATED] -> ", dir_gsea_contrast))
    }
    
    # Normalize IDs to match functional databases
    clean_names <- normalize_ids(as.character(master_univ$PG.ProteinAccessions))
    gene_list <- sort(setNames(master_univ[[lfc_col]], clean_names), decreasing = TRUE)
    
    # Remove NAs and duplicates
    gene_list <- gene_list[!is.na(gene_list)]
    gene_list <- gene_list[!duplicated(names(gene_list))]
    
    for (db_name in names(all_analyses)) {
      cat(paste0("\n  > Running GSEA: ", contrast_id, " | Database: ", db_name))
      
      res_gsea <- tryCatch({
        clusterProfiler::GSEA(geneList     = gene_list, 
                              TERM2GENE    = all_analyses[[db_name]], 
                              pvalueCutoff = 0.05, 
                              eps          = 0,
                              seed         = TRUE)
      }, error = function(e) {
        cat(" | Error: ", e$message)
        return(NULL)
      })
      
      if (!is.null(res_gsea) && nrow(as.data.frame(res_gsea)) > 0) {
        fname <- paste0("GSEA_", db_name)
        df_res <- as.data.frame(res_gsea)
        
        # --- ANNOTATION ---
        df_res$Description_Full <- "No Description Found"
        try({
          if (grepl("GO", db_name)) {
            term_map <- AnnotationDbi::select(GO.db, keys = as.character(df_res$ID), columns = "TERM", keytype = "GOID")
            if (!is.null(term_map)) df_res$Description_Full <- term_map$TERM[match(df_res$ID, term_map$GOID)]
          } else {
            clean_ids <- gsub("ko:", "", as.character(df_res$ID))
            matched_desc <- dict_extra$Desc[match(clean_ids, dict_extra$ID)]
            if (any(!is.na(matched_desc))) df_res$Description_Full <- ifelse(is.na(matched_desc), "No Description Found", matched_desc)
          }
        }, silent = TRUE)
        
        # Save CSV with Full Path Feedback
        csv_out <- file.path(dir_gsea_contrast, paste0(fname, "_Annotated.csv"))
        write.csv(df_res, csv_out, row.names = FALSE)
        cat(paste0(" | [SAVED] ", basename(csv_out)))
        
        # --- PLOTS ---
        p1 <- enrichplot::dotplot(res_gsea, showCategory = 20, split = ".sign") + 
          ggplot2::facet_grid(.~.sign) + ggplot2::theme_bw()
        ggplot2::ggsave(file.path(dir_gsea_contrast, paste0(fname, "_dotplot.tiff")), plot = p1, width = 12, height = 8, dpi = 300)
        
        try({
          p2 <- enrichplot::gseaplot2(res_gsea, geneSetID = 1:min(3, nrow(res_gsea)), 
                                      title = paste(contrast_id, "-", db_name))
          ggplot2::ggsave(file.path(dir_gsea_contrast, paste0(fname, "_mountains.tiff")), plot = p2, width = 12, height = 8, dpi = 300)
        }, silent = TRUE)
        
        cat(" -> SUCCESS!")
      } else {
        cat(" -> No significant results.")
      }
    }
  }
  
  cat("\n\n--- FUNCTIONAL ENRICHMENT MODULE COMPLETED ---")
  
  # ==============================================================================
  # 14. FINAL INTEGRATION: ANNOTATION & GLOBAL REPORTING
  # =================-------------------------------------------------------------
  cat("\n--- STARTING FINAL POST-PROCESSING (GLOBAL INTEGRATION) ---")
  
  # 14a. Dictionary building for Annotation
  # ------------------------------------------------------------------------------
  dict_extra <- egg_df_cleaned %>% 
    dplyr::select(ID = 12, Desc = 13) %>% 
    dplyr::bind_rows(egg_df_cleaned %>% dplyr::select(ID = 21, Desc = 21)) %>% 
    dplyr::filter(!is.na(ID) & ID != "-") %>%
    tidyr::separate_rows(ID, sep = ",") %>% 
    dplyr::mutate(ID = gsub("ko:", "", ID)) %>% 
    dplyr::distinct(ID, .keep_all = TRUE)
  
  # 14b. CREATE ANNOTATED FILES (ORA & GSEA Mapping)
  # ------------------------------------------------------------------------------
  result_paths <- c(file.path(output_main, "10_ORA"), file.path(output_main, "11_GSEA"))
  all_files <- list.files(path = result_paths, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  # Skip files that are already annotated
  all_files <- all_files[!grepl("_Annotated.csv$", all_files)]
  
  cat("\n- Mapping descriptions to enrichment files...")
  
  for (f in all_files) {
    df_res <- read.csv(f); if (nrow(df_res) == 0) next
    df_res$Description_Full <- "No Description Found"
    try({
      if (grepl("GO", f)) {
        term_map <- AnnotationDbi::select(GO.db, keys = as.character(df_res$ID), columns = "TERM", keytype = "GOID")
        if (!is.null(term_map)) df_res$Description_Full <- term_map$TERM[match(df_res$ID, term_map$GOID)]
      } else {
        clean_ids <- gsub("ko:", "", as.character(df_res$ID))
        matched_desc <- dict_extra$Desc[match(clean_ids, dict_extra$ID)]
        if (any(!is.na(matched_desc))) df_res$Description_Full <- ifelse(is.na(matched_desc), "No Description Found", matched_desc)
      }
    }, silent = TRUE)
    write.csv(df_res, gsub("\\.csv$", "_Annotated.csv", f), row.names = FALSE)
  }
  
  # 14c. Global Protein Flow Statistics
  # ------------------------------------------------------------------------------
  cat("\n- Calculating protein flow statistics...")
  
  get_count <- function(folder_name) {
    f_path <- file.path(output_main, folder_name)
    csv_files <- list.files(f_path, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
    if(length(csv_files) > 0) return(nrow(read.csv(csv_files[1]))) else return(0)
  }
  
  n_raw      <- if(exists("raw_prot_data")) nrow(raw_prot_data) else get_count("1_Cleaned_File")
  n_cleaned  <- if(exists("prot_data_cleaned")) nrow(prot_data_cleaned) else get_count("2_QC_correction")
  n_master   <- get_count("6_DEA")
  
  # Identification of unique Elite Protein IDs
  elite_path <- file.path(output_main, "9_CBR_Filter")
  elite_files <- list.files(elite_path, pattern = "Final_Elite", recursive = TRUE, full.names = TRUE)
  n_elite_ids <- if(length(elite_files) > 0) {
    unique(unlist(lapply(elite_files, function(x) as.character(read.csv(x)[[1]]))))
  } else { character(0) }
  n_elite <- length(n_elite_ids)
  
  report_final <- data.frame(
    Pipeline_Stage = c("1. Initial Proteins", "2. Proteins Removed", "3. Total Quantified", "4. Significant (DEA)", "5. Elite (Bio-Relevance)"),
    Count = c(n_raw, n_raw - n_cleaned, n_cleaned, n_master, n_elite),
    Percentage = paste0(round(c(n_raw, n_raw - n_cleaned, n_cleaned, n_master, n_elite) / n_raw * 100, 2), "%")
  )
  write.csv(report_final, file.path(output_main, "Final_Pipeline_Report.csv"), row.names = FALSE)
  
  # 14d. Visualizing Method Intersection (Stacked Bar Plot)
  # ------------------------------------------------------------------------------
  cat("\n- Generating method intersection visuals...")
  master_path <- file.path(output_main, "8_EDA_and_CoExpression", "Master_Integrated_Filtered_for_Bayesian.csv")
  
  if(file.exists(master_path)) {
    master_int <- read.csv(master_path)
    elite_distribution <- master_int %>%
      dplyr::filter(PG.ProteinAccessions %in% n_elite_ids) %>%
      dplyr::group_by(Category) %>%
      dplyr::tally() %>%
      dplyr::mutate(Percentage = n / sum(n) * 100)
    
    p_overlap <- ggplot(elite_distribution, aes(x = "Elite Proteins", y = n, fill = Category)) +
      geom_bar(stat = "identity", width = 0.5) +
      geom_text(aes(label = paste0(n, "\n(", round(Percentage, 1), "%)")), 
                position = position_stack(vjust = 0.5), color = "white", fontface = "bold") +
      scale_fill_manual(values = c("Both" = "#2ecc71", "DEA Only" = "#3498db", "WGCNA Only" = "#f1c40f")) +
      theme_minimal() +
      labs(title = "Origin of Elite Proteins",
           subtitle = "Contribution of DEA (DEqMS) and Co-expression (WGCNA)",
           x = "", y = "Protein Count", fill = "Method Source") +
      theme(axis.text.x = element_blank())
    
    ggsave(file.path(output_main, "Elite_Proteins_Method_Overlap.png"), p_overlap, width = 6, height = 8)
  }
  
  # 14e. Global Flow Summary Plot
  # ------------------------------------------------------------------------------
  viz_data <- report_final
  viz_data$Stage <- factor(viz_data$Pipeline_Stage, levels = rev(viz_data$Pipeline_Stage))
  p_flow <- ggplot(viz_data, aes(x = Stage, y = Count, fill = Stage)) +
    geom_bar(stat = "identity", width = 0.7, show.legend = FALSE) +
    geom_text(aes(label = paste0(Count, " (", Percentage, ")")), hjust = -0.1) +
    coord_flip() + scale_fill_brewer(palette = "Blues") + theme_minimal() +
    labs(title = "Protein Flow Summary", x = "", y = "Protein Count")
  ggsave(file.path(output_main, "Protein_Flow_Summary.png"), p_flow, width = 8, height = 5)
  
  cat("\n[SUCCESS] Block 14 Completed. Reports and Visuals saved.")
  
  # ==============================================================================
  # 15. RANDOM FOREST VALIDATION (THE COMPLETE TRILOGY)
  # ==============================================================================
  cat("\n[STEP 15] Starting Full Validation: Importance + Robustness + Permutation...")
  
  active_contrasts <- contrast_matrix$Contrast_ID[toupper(contrast_matrix$Active) == "YES"]
  
  for (current_contrast in active_contrasts) {
    cat(paste0("\n\n>>> FULL VALIDATION FOR: ", current_contrast))
    
    # 1. Load Data
    elite_path <- file.path(input_9, current_contrast, paste0("Final_Elite_Proteins_", current_contrast, ".csv"))
    if(!file.exists(elite_path)) next
    elite_ids <- readr::read_csv(elite_path, show_col_types = FALSE)[[protein_id_col]]
    
    rf_input <- master_data %>%
      mutate(
        is_elite = ifelse(.[[protein_id_col]] %in% elite_ids, "Elite", "Non_Elite"),
        Log2FC = .[[paste0("logFC_", current_contrast)]],
        AdjP   = .[[paste0("adj.P.Val_", current_contrast)]],
        Relevant_kME = mapply(function(m, row_idx) {
          if(is.na(m) || m == "grey") return(0)
          col_name <- paste0("kME", m)
          return(as.numeric(master_data[row_idx, col_name]))
        }, WGCNA_Module, row_number())
      ) %>%
      select(is_elite, Log2FC, AdjP, Relevant_kME) %>%
      na.omit() %>%
      mutate(across(c(Log2FC, AdjP, Relevant_kME), ~as.vector(scale(.))))
    
    # --- GRAPH 1: FEATURE IMPORTANCE ---
    set.seed(123)
    rf_real <- randomForest(as.factor(is_elite) ~ ., data = rf_input, ntree = 500, importance = TRUE)
    acc_real <- 1 - mean(rf_real$err.rate[,1])
    
    imp_df <- as.data.frame(importance(rf_real))
    imp_df$Variable <- rownames(imp_df)
    
    p1 <- ggplot(imp_df, aes(x = reorder(Variable, MeanDecreaseGini), y = MeanDecreaseGini, fill = MeanDecreaseGini)) +
      geom_bar(stat = "identity") + coord_flip() + theme_minimal() +
      scale_fill_gradient(low = "#D4E6F1", high = "#1A5276") +
      labs(title = paste("Predictor Importance:", current_contrast), subtitle = paste("Accuracy:", round(acc_real*100,2), "%"))
    
    ggsave(file.path(output_12, paste0("1_Importance_", current_contrast, ".tiff")), plot = p1, device = "tiff", dpi = 300)
    
    # --- GRAPH 2: ROBUSTNESS BAR CHART ---
    set.seed(42)
    perm_labels_single <- sample(rf_input$is_elite)
    rf_single_perm <- randomForest(as.factor(perm_labels_single) ~ ., data = rf_input, ntree = 200)
    acc_perm_single <- 1 - mean(rf_single_perm$err.rate[,1])
    
    df_rob <- data.frame(Type = c("Random Chance", "Real Model"), Accuracy = c(acc_perm_single, acc_real))
    p2 <- ggplot(df_rob, aes(x = Type, y = Accuracy*100, fill = Type)) +
      geom_bar(stat = "identity", width = 0.5) + geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
      theme_minimal() + scale_fill_manual(values = c("Real Model" = "#1A5276", "Random Chance" = "#D5DBDB")) +
      labs(title = paste("Robustness Check:", current_contrast), y = "Accuracy (%)")
    
    ggsave(file.path(output_12, paste0("2_Robustness_", current_contrast, ".tiff")), plot = p2, device = "tiff", dpi = 300)
    
    # --- GRAPH 3: NULL DISTRIBUTION (1000 PERMUTATIONS) ---
    cat("\n   Simulating 1000 permutations...")
    null_dist <- numeric(1000)
    for (i in 1:1000) {
      p_labels <- sample(rf_input$is_elite)
      rf_tmp <- randomForest(as.factor(p_labels) ~ ., data = rf_input, ntree = 50) # ntree bajo para velocidad
      null_dist[i] <- 1 - mean(rf_tmp$err.rate[,1])
    }
    
    p3 <- ggplot(data.frame(Acc = null_dist), aes(x = Acc*100)) +
      geom_histogram(bins = 30, fill = "#AEB6BF", color = "white") +
      geom_vline(xintercept = acc_real*100, linetype = "dashed", color = "#1A5276", size = 1) +
      theme_minimal() + labs(title = paste("Permutation Distribution:", current_contrast), x = "Accuracy (%)", y = "Frequency")
    
    ggsave(file.path(output_12, paste0("3_NullDist_", current_contrast, ".tiff")), plot = p3, device = "tiff", dpi = 300)
    
    # --- DATA TABLE ---
    write_csv(imp_df, file.path(output_12, paste0("Data_Detail_", current_contrast, ".csv")))
  }
  cat("\n\n--- [TRILOGY COMPLETED] ---")