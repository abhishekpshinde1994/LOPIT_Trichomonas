###############################################################################
## LOPIT Analysis: SVM for 3 TMT Replicates
## Outputs: marker profiles, F1 scores, UMAP 3D plots (PNG + HTML), CSV predictions
## Author: Abhishek Prakash Shinde
###############################################################################

# ---- Libraries ----
library(tidyverse)
library(readxl)
library(MSnbase)
library(pRoloc)
library(pRolocdata)
library(plotly)
library(htmlwidgets)
library(uwot)
library(BiocParallel)
library(scales)

set.seed(42)

# ---- Configuration ----
data_file   <- "TMT_replicates.xlsx"
marker_file <- "markers_v6.xlsx"
out_base    <- "LOPIT_results"

# Abundance columns 20:33 are consistent across all 4 sheets (14 fractions)
# ReplicateA/B/C: TVAG ID extracted via GN= regex from column 5 (Description)
# Pilot:          TVAG ID is directly in column 112 (TVAG column)
# SVM cutoff 0.50; 3D UMAP used for visualisation
sheet_info <- list(
  #ReplicateA = list(sheet = "ReplicateA", abund_cols = 20:33, id_col = 5,   id_mode = "regex"),
  ReplicateB = list(sheet = "ReplicateB", abund_cols = 20:33, id_col = 5,   id_mode = "regex"),
  ReplicateC = list(sheet = "ReplicateC", abund_cols = 20:33, id_col = 5,   id_mode = "regex"),
  Pilot      = list(sheet = "Pilot",      abund_cols = 20:33, id_col = 112, id_mode = "direct")
)

desc_col            <- 5     # Description column (all sheets)
svm_t_cutoff        <- 0.50
svm_per_compartment <- TRUE
umap_n_neighbors    <- 15    # UMAP: local neighbourhood size (increase for more global structure)
umap_min_dist       <- 0.1   # UMAP: compactness of clusters (lower = tighter)
umap_metric         <- "euclidean"

dir.create(out_base, showWarnings = FALSE, recursive = TRUE)

###############################################################################
## HELPER FUNCTIONS
###############################################################################

#' Extract TVAG IDs depending on sheet type
extract_ids <- function(raw, id_col, id_mode) {
  if (id_mode == "regex") {
    # Column 5 contains full description; pull TVAG_XXXXXX from GN= field
    trimws(str_extract(as.character(raw[[id_col]]), "(?<=GN=)TVAG_\\d+"))
  } else {
    # Column 112 in Pilot has clean TVAG_XXXXXX directly
    trimws(as.character(raw[[id_col]]))
  }
}


#' Load and prepare a single replicate as MSnSet
prepare_replicate <- function(data_file, sheet_name, abund_cols, id_col,
                              id_mode, desc_col, markers_df) {
  
  cat("  Loading sheet:", sheet_name, "\n")
  raw <- read_xlsx(data_file, sheet = sheet_name)
  
  ids   <- extract_ids(raw, id_col, id_mode)
  descs <- as.character(raw[[desc_col]])
  abund <- as.data.frame(raw[, abund_cols])
  
  # Drop rows with no extractable TVAG ID
  has_id <- !is.na(ids) & ids != "" & grepl("^TVAG_", ids)
  ids    <- ids[has_id];  descs <- descs[has_id];  abund <- abund[has_id, ]
  
  # Convert to numeric
  abund <- mutate_all(abund, as.double)
  
  # Keep only rows with complete abundance across all 14 fractions
  complete_rows <- complete.cases(abund)
  ids   <- ids[complete_rows]
  descs <- descs[complete_rows]
  abund <- abund[complete_rows, ]
  
  # Remove duplicate accessions (keep first)
  keep  <- !duplicated(ids)
  ids   <- ids[keep];  descs <- descs[keep];  abund <- abund[keep, ]
  
  # Row-normalise (standard LOPIT)
  abund_norm         <- abund / rowSums(abund, na.rm = TRUE)
  expr_mat           <- as.matrix(abund_norm)
  rownames(expr_mat) <- ids
  colnames(expr_mat) <- paste0("F", seq_len(ncol(expr_mat)))
  
  # Feature data with marker assignment
  fdata <- data.frame(
    Accession   = ids,
    Description = descs,
    markers     = "unknown",
    stringsAsFactors = FALSE,
    row.names   = ids
  )
  m <- match(ids, markers_df$ID)
  fdata$markers[!is.na(m)] <- markers_df$markers[m[!is.na(m)]]
  
  cat("    Proteins:", nrow(expr_mat), "| Markers matched:", sum(!is.na(m)), "\n")
  
  pdata <- data.frame(
    sampleNames = colnames(expr_mat),
    Fraction    = seq_len(ncol(expr_mat)),
    row.names   = colnames(expr_mat)
  )
  
  MSnSet(
    exprs = expr_mat,
    fData = new("AnnotatedDataFrame", data = fdata),
    pData = new("AnnotatedDataFrame", data = pdata)
  )
}


#' Run SVM optimisation + classification
run_svm_analysis <- function(msnset, out_dir, suffix) {
  
  cat("\n=== SVM:", suffix, "===\n")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  n_markers <- sum(fData(msnset)$markers != "unknown")
  if (n_markers == 0) stop("No marker proteins found for '", suffix, "' — check ID matching.")
  cat("  Marker proteins:", n_markers, "\n")
  
  # Inverse-frequency class weights to handle imbalance
  marker_tbl <- table(getMarkers(msnset))
  weights    <- 1 / marker_tbl[names(marker_tbl) != "unknown"]
  
  cat("  Optimising SVM hyperparameters...\n")
  params <- svmOptimisation(msnset, fcol = "markers", class.weights = weights)
  save(params, file = file.path(out_dir, paste0("svm_params_", suffix, ".Rdata")))
  
  png(file.path(out_dir, paste0("F1_optimisation_", suffix, ".png")),
      width = 2000, height = 1500, res = 300)
  plot(params); dev.off()
  
  png(file.path(out_dir, paste0("F1_levelplot_", suffix, ".png")),
      width = 2000, height = 1500, res = 300)
  print(levelPlot(params)); dev.off()
  
  cat("  Running SVM classification...\n")
  svmres <- svmClassification(msnset, fcol = "markers", assessRes = params)
  save(svmres, file = file.path(out_dir, paste0("SVMres_", suffix, ".Rdata")))
  
  # Per-compartment score thresholds
  if (svm_per_compartment) {
    thresholds  <- orgQuants(svmres, fcol = "svm", scol = "svm.scores",
                             mcol = "markers", t = svm_t_cutoff)
    svmres_pred <- getPredictions(svmres, fcol = "svm", scol = "svm.scores",
                                  mcol = "markers", t = thresholds)
  } else {
    svmres_pred <- getPredictions(svmres, fcol = "svm", t = svm_t_cutoff)
  }
  
  # Score distribution plot
  p_box <- ggplot(fData(svmres), aes(x = svm, y = svm.scores)) +
    geom_boxplot(fill = "lightblue") +
    geom_hline(yintercept = svm_t_cutoff, color = "red", linetype = 2) +
    coord_flip() + theme_minimal(base_size = 12) +
    labs(x = "Compartment", y = "SVM Score",
         title = paste0("SVM Scores - ", suffix))
  ggsave(file.path(out_dir, paste0("SVM_scores_", suffix, ".png")),
         p_box, width = 10, height = 6, dpi = 300)
  
  cat("  SVM complete.\n")
  list(svmres = svmres, svmres_pred = svmres_pred, params = params)
}


#' Plot organellar marker profiles
plot_marker_profiles <- function(msnset, out_dir, suffix) {
  
  fdat           <- fData(msnset)
  expr           <- exprs(msnset)
  marker_classes <- unique(fdat$markers[fdat$markers != "unknown"])
  
  if (length(marker_classes) == 0) {
    cat("  No markers — skipping profile plots.\n"); return(invisible(NULL))
  }
  
  marker_dir <- file.path(out_dir, "marker_profiles")
  dir.create(marker_dir, showWarnings = FALSE, recursive = TRUE)
  
  all_data <- lapply(marker_classes, function(mc) {
    idx <- which(fdat$markers == mc)
    as.data.frame(expr[idx, , drop = FALSE]) %>%
      mutate(Accession = rownames(expr)[idx]) %>%
      pivot_longer(-Accession, names_to = "Fraction", values_to = "Abundance") %>%
      mutate(Fraction = factor(Fraction, levels = colnames(expr)), Compartment = mc)
  })
  
  # Individual compartment plots
  for (i in seq_along(marker_classes)) {
    mc        <- marker_classes[i]
    safe_name <- str_replace_all(mc, "[^A-Za-z0-9_]", "_")
    p <- ggplot(all_data[[i]], aes(x = Fraction, y = Abundance, group = Accession)) +
      geom_line(alpha = 0.4, color = "steelblue") +
      geom_point(size = 1, alpha = 0.4, color = "steelblue") +
      stat_summary(aes(group = 1), fun = mean, geom = "line",  color = "red", linewidth = 1.2) +
      stat_summary(aes(group = 1), fun = mean, geom = "point", color = "red", size = 2) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste0(mc, " - ", suffix), x = "Fraction", y = "Normalised Abundance")
    ggsave(file.path(marker_dir, paste0(safe_name, "_", suffix, ".png")),
           p, width = 8, height = 5, dpi = 300)
  }
  
  # Combined summary
  summary_df <- bind_rows(all_data) %>%
    group_by(Compartment, Fraction) %>%
    summarise(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop")
  
  p_all <- ggplot(summary_df, aes(x = Fraction, y = mean_abund,
                                  color = Compartment, group = Compartment)) +
    geom_line(linewidth = 1) + geom_point(size = 2) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0("All Marker Profiles - ", suffix),
         x = "Fraction", y = "Mean Normalised Abundance", color = "Compartment")
  ggsave(file.path(out_dir, paste0("all_marker_profiles_", suffix, ".png")),
         p_all, width = 12, height = 7, dpi = 300)
  
  cat("  Marker profiles saved.\n")
}


#' UMAP 3D + SVM prediction overlay -> PNG, interactive HTML, CSV
generate_umap_outputs <- function(msnset, out_dir, suffix) {
  
  cat("  Generating 3D UMAP...\n")
  expr_mat <- exprs(msnset)
  fdat     <- fData(msnset)
  
  # 3D UMAP embedding
  set.seed(42)
  coords <- umap(expr_mat,
                 n_components = 3,
                 n_neighbors  = umap_n_neighbors,
                 min_dist     = umap_min_dist,
                 metric       = umap_metric,
                 verbose      = FALSE)
  
  # Use thresholded predictions (svm.pred); fall back to raw svm if absent
  pred_col  <- if ("svm.pred"   %in% colnames(fdat)) "svm.pred"   else "svm"
  score_col <- if ("svm.scores" %in% colnames(fdat)) "svm.scores" else NULL
  
  plot_df <- tibble(
    Accession   = rownames(expr_mat),
    UMAP1       = coords[, 1],
    UMAP2       = coords[, 2],
    UMAP3       = coords[, 3],
    Marker      = fdat$markers,
    Prediction  = fdat[[pred_col]],
    Score       = if (!is.null(score_col)) round(as.numeric(fdat[[score_col]]), 4) else NA_real_,
    Description = fdat$Description
  ) %>%
    mutate(label = paste0(
      "ID: ",          Accession,
      "\nPrediction: ", Prediction,
      "\nScore: ",      Score,
      "\nMarker: ",     Marker,
      "\n",             str_trunc(Description, 60)
    ))
  
  # Colour scheme: one colour per compartment, grey for unknown
  compartments <- sort(unique(fdat$markers[fdat$markers != "unknown"]))
  all_levels   <- c(compartments, "unknown")
  colors       <- setNames(c(hue_pal()(length(compartments)), "grey70"), all_levels)
  
  plot_df$Prediction <- factor(plot_df$Prediction, levels = all_levels)
  marker_df          <- filter(plot_df, Marker != "unknown")
  
  # ---- Static 2D PNG (UMAP1 vs UMAP2 projection, for publication) ----
  p2d <- ggplot() +
    geom_point(data = plot_df,
               aes(x = UMAP1, y = UMAP2, color = Prediction, text = label),
               size = 1.5, alpha = 0.5) +
    geom_point(data = marker_df,
               aes(x = UMAP1, y = UMAP2),
               color = "black", size = 3.5) +
    geom_point(data = marker_df,
               aes(x = UMAP1, y = UMAP2, color = factor(Marker, levels = all_levels)),
               size = 2.2) +
    scale_color_manual(values = colors, drop = FALSE, name = "Compartment") +
    theme_minimal(base_size = 14) +
    labs(title = paste0("UMAP SVM Predictions (2D projection) - ", suffix))
  
  ggsave(file.path(out_dir, paste0("UMAP_SVM_", suffix, ".png")),
         p2d, width = 12, height = 8, dpi = 300)
  
  # ---- Interactive 3D HTML ----
  unknown_df <- filter(plot_df, Marker == "unknown")
  known_df   <- filter(plot_df, Marker != "unknown")
  
  fig3d <- plot_ly() %>%
    
    # Layer 1: predicted-unknown proteins (semi-transparent background)
    add_trace(
      data       = unknown_df,
      type       = "scatter3d",
      mode       = "markers",
      x          = ~UMAP1,
      y          = ~UMAP2,
      z          = ~UMAP3,
      color      = ~Prediction,
      colors     = colors,
      text       = ~label,
      hoverinfo  = "text",
      marker     = list(size = 3, opacity = 0.45),
      showlegend = FALSE,
      name       = ~as.character(Prediction)
    ) %>%
    
    # Layer 2: all proteins with legend
    add_trace(
      data       = plot_df,
      type       = "scatter3d",
      mode       = "markers",
      x          = ~UMAP1,
      y          = ~UMAP2,
      z          = ~UMAP3,
      color      = ~Prediction,
      colors     = colors,
      text       = ~label,
      hoverinfo  = "text",
      marker     = list(size = 3, opacity = 0.50),
      showlegend = TRUE,
      name       = ~as.character(Prediction)
    ) %>%
    
    # Layer 3: confirmed marker proteins (larger, black border)
    add_trace(
      data       = known_df,
      type       = "scatter3d",
      mode       = "markers",
      x          = ~UMAP1,
      y          = ~UMAP2,
      z          = ~UMAP3,
      color      = ~factor(Marker, levels = all_levels),
      colors     = colors,
      text       = ~label,
      hoverinfo  = "text",
      marker     = list(
        size    = 7,
        opacity = 1.0,
        line    = list(color = "black", width = 1.5)
      ),
      showlegend = FALSE,
      name       = ~paste0("Marker: ", Marker)
    ) %>%
    
    layout(
      title  = list(
        text = paste0("<b>3D UMAP — SVM Predictions: ", suffix, "</b>"),
        font = list(size = 16)
      ),
      scene  = list(
        xaxis = list(title = "UMAP 1", backgroundcolor = "rgb(240,240,240)",
                     gridcolor = "white", showbackground = TRUE),
        yaxis = list(title = "UMAP 2", backgroundcolor = "rgb(240,240,240)",
                     gridcolor = "white", showbackground = TRUE),
        zaxis = list(title = "UMAP 3", backgroundcolor = "rgb(240,240,240)",
                     gridcolor = "white", showbackground = TRUE),
        camera = list(eye = list(x = 1.5, y = 1.5, z = 1.2))
      ),
      legend = list(
        title         = list(text = "<b>Compartment</b>"),
        itemsizing    = "constant",
        tracegroupgap = 2
      ),
      margin = list(l = 0, r = 0, b = 0, t = 60)
    )
  
  saveWidget(fig3d,
             file          = file.path(normalizePath(out_dir),
                                       paste0("UMAP3D_SVM_", suffix, ".html")),
             selfcontained = TRUE)
  
  # ---- Predictions CSV (includes all 3 UMAP coords) ----
  plot_df %>%
    select(Accession, Description, Marker, Prediction, Score, UMAP1, UMAP2, UMAP3) %>%
    write_csv(file.path(out_dir, paste0("predictions_SVM_", suffix, ".csv")))
  
  cat("  UMAP 2D PNG + 3D HTML saved. Predictions CSV exported.\n")
}


###############################################################################
## MAIN PIPELINE
###############################################################################

# ---- Load markers ----
cat("Loading markers...\n")
markers_raw <- read_xlsx(marker_file, col_names = TRUE)
markers_df  <- tibble(
  ID      = trimws(as.character(markers_raw[[1]])),
  markers = trimws(as.character(markers_raw[[2]]))
) %>% filter(!is.na(ID), !is.na(markers), ID != "NA", ID != "ID")

cat("Markers loaded:", nrow(markers_df), "proteins across",
    n_distinct(markers_df$markers), "compartments\n")
cat("Compartments:", paste(sort(unique(markers_df$markers)), collapse = ", "), "\n\n")


# ---- Process each replicate ----
replicate_results <- list()

for (rep_name in names(sheet_info)) {
  
  cat("\n", strrep("=", 70), "\n")
  cat("PROCESSING:", rep_name, "\n")
  cat(strrep("=", 70), "\n")
  
  info    <- sheet_info[[rep_name]]
  rep_dir <- file.path(out_base, rep_name)
  dir.create(rep_dir, showWarnings = FALSE, recursive = TRUE)
  
  msnset <- prepare_replicate(
    data_file  = data_file,
    sheet_name = info$sheet,
    abund_cols = info$abund_cols,
    id_col     = info$id_col,
    id_mode    = info$id_mode,
    desc_col   = desc_col,
    markers_df = markers_df
  )
  
  plot_marker_profiles(msnset, rep_dir, rep_name)
  svm_result <- run_svm_analysis(msnset, rep_dir, rep_name)
  generate_umap_outputs(svm_result$svmres_pred, rep_dir, rep_name)
  
  replicate_results[[rep_name]] <- list(msnset = msnset, svm_result = svm_result)
  cat(rep_name, "COMPLETE\n")
}


###############################################################################
## OVERLAP ANALYSIS
###############################################################################

cat("\n", strrep("=", 70), "\n")
cat("OVERLAP ANALYSIS\n")
cat(strrep("=", 70), "\n\n")

overlap_dir <- file.path(out_base, "Overlap")
dir.create(overlap_dir, showWarnings = FALSE, recursive = TRUE)

# Proteins present in all 4 replicates
all_ids       <- lapply(replicate_results, function(r) rownames(exprs(r$msnset)))
core_proteins <- Reduce(intersect, all_ids)

cat("Proteins per replicate:\n")
for (nm in names(all_ids)) cat(" ", nm, ":", length(all_ids[[nm]]), "\n")
cat("Core overlap (all 4):", length(core_proteins), "\n\n")

write_csv(tibble(Accession = sort(core_proteins)),
          file.path(overlap_dir, "core_overlapping_proteins.csv"))

# Average expression across replicates then re-normalise
expr_list <- lapply(replicate_results,
                    function(r) exprs(r$msnset)[core_proteins, , drop = FALSE])
avg_expr  <- Reduce("+", expr_list) / length(expr_list)
avg_expr  <- avg_expr / rowSums(avg_expr)

first_fdata   <- fData(replicate_results[[1]]$msnset)
overlap_fdata <- data.frame(
  Accession   = core_proteins,
  Description = first_fdata[core_proteins, "Description"],
  markers     = "unknown",
  stringsAsFactors = FALSE,
  row.names   = core_proteins
)
m <- match(core_proteins, markers_df$ID)
overlap_fdata$markers[!is.na(m)] <- markers_df$markers[m[!is.na(m)]]
cat("Markers in overlap set:", sum(!is.na(m)), "\n")

overlap_msnset <- MSnSet(
  exprs = avg_expr,
  fData = new("AnnotatedDataFrame", data = overlap_fdata),
  pData = new("AnnotatedDataFrame", data = data.frame(
    sampleNames = colnames(avg_expr),
    Fraction    = seq_len(ncol(avg_expr)),
    row.names   = colnames(avg_expr)
  ))
)

plot_marker_profiles(overlap_msnset, overlap_dir, "Overlap")
svm_overlap <- run_svm_analysis(overlap_msnset, overlap_dir, "Overlap")
generate_umap_outputs(svm_overlap$svmres_pred, overlap_dir, "Overlap")


# ---- Consensus across replicates ----
cat("\nBuilding SVM consensus...\n")

consensus_df <- tibble(Accession = core_proteins)
for (rep_name in names(replicate_results)) {
  sfd <- fData(replicate_results[[rep_name]]$svm_result$svmres)
  consensus_df[[paste0("SVM_", rep_name)]] <- sfd[core_proteins, "svm"]
}

consensus_df <- consensus_df %>%
  rowwise() %>%
  mutate(
    consensus_pred  = { p <- c_across(starts_with("SVM_")); names(which.max(table(p))) },
    consensus_count = { p <- c_across(starts_with("SVM_")); max(table(p)) }
  ) %>%
  ungroup() %>%
  mutate(Description = first_fdata[core_proteins, "Description"]) %>%
  left_join(markers_df %>% rename(Accession = ID, MarkerLabel = markers),
            by = "Accession") %>%
  mutate(MarkerLabel = replace_na(MarkerLabel, "unknown"))

write_csv(consensus_df, file.path(overlap_dir, "SVM_consensus_predictions.csv"))

cat("Consensus 4/4:", sum(consensus_df$consensus_count == 4), "\n")
cat("Consensus 3/4:", sum(consensus_df$consensus_count == 3), "\n")
cat("Consensus 2/4:", sum(consensus_df$consensus_count == 2), "\n")
cat("Consensus 1/4:", sum(consensus_df$consensus_count == 1), "\n")


# ---- QSep spatial resolution ----
cat("\nRunning QSep...\n")
qsep_dir <- file.path(overlap_dir, "spatial_resolution")
dir.create(qsep_dir, showWarnings = FALSE, recursive = TRUE)

tryCatch({
  qsep_res <- QSep(object = overlap_msnset, fcol = "markers")
  
  for (nm in c("boxplot", "boxplot_norm", "heatmap", "heatmap_norm")) {
    norm <- grepl("norm", nm)
    png(file.path(qsep_dir, paste0("QSep_", nm, ".png")),
        width = 2000, height = if (grepl("heatmap", nm)) 2000 else 1500, res = 300)
    if (grepl("heatmap", nm)) levelPlot(qsep_res, norm = norm)
    else plot(qsep_res, norm = norm)
    dev.off()
  }
  cat("QSep plots saved.\n")
}, error = function(e) cat("QSep failed:", e$message, "\n"))


###############################################################################
## DONE
###############################################################################
cat("\n", strrep("=", 70), "\n")
cat("ANALYSIS COMPLETE — outputs in:", out_base, "\n")
cat(strrep("=", 70), "\n")