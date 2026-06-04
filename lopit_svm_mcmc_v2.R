library(tidyverse)    # data parsing
library(plotly)       # data visualization
library(scales)       # data scaling
library(coda)         # for gelman.diag
library(pRoloc)       # LOPIT analysis
library(xlsx)
library(pheatmap)     # for heatmap
library(dbscan)
library(pRolocGUI)
library(bandle)
library(ggforce)      # library for plotting circles in ggplot

#######################################################
## If chain-failing problems repeat:                 ##
#######################################################
# Temporarily replace MulticoreParam with SerialParam() to check whether
# the problem disappears (issue with chains failing).
# If yes, this confirms it is parallelization/environment related.

###########################
## Functions for calling ##
###########################

# Function that creates a consensus prediction from multiple LOPITs
lopit_consensus <- function(data,
                            smv_pred_pos,
                            consensusCount,
                            compartments) {
  # data            = data with the svm.predictions
  # smv_pred_pos    = position of the svm.pred column from the lopit function
  # consensusCount  = how many times a protein must be seen in a given compartment to achieve consensus
  # compartments    = a vector of all compartment names

  countedCompartments_T <- map(compartments,
                               ~{
                                 # select predictions
                                 allPredTien2_noNA <- data %>%
                                   dplyr::select(all_of(smv_pred_pos))

                                 # count predictions per LOPIT per row
                                 countedComp <- apply(allPredTien2_noNA, 1, str_count) %>%
                                   t() %>%
                                   apply(., 1, sum, na.rm = T) %>%
                                   unlist()

                                 # make a tibble from it
                                 countedComp_T <- tibble(countedComp)
                                 colnames(countedComp_T) <- .x

                                 return(countedComp_T)

                               }) %>%
    purrr::reduce(cbind)

  # create the consensus compartment
  consensusCompartment_Sch <- map(array_branch(countedCompartments_T, 1),
                                  ~{
                                    # get which LOPITs have seen the protein in the same compartment more than twice
                                    seenMoreThan2 <- which(.x >= consensusCount)

                                    if (length(seenMoreThan2) == 0) {
                                      "unknown"
                                    } else {
                                      paste0(names(seenMoreThan2), collapse = "|")
                                    }
                                  }) %>%
    unlist()

  return(consensusCompartment_Sch)

}

# Function that takes a vector of IDs and always returns the same ID as key and
# other IDs for the group (based on my matching of IDs for T. vaginalis)
lopit_get_uniqueID <- function(tibble,
                               IDcol,
                               key_tibble,
                               sep = ";") {
  # tibble      = tibble with ID column where the spectrometer thinks the IDs form a single group;
  #               expects a vector of separated IDs
  # IDcol       = position of the ID column
  # key_tibble  = tibble with ID matches among TVAG, TVAGG3 and UniProt (general datasets file)
  # from        = mapping direction: either from TVAG, TVAGG3 or Entry (UniProt), e.g. "TVAG"
  # sep         = if a single string with IDs separated by a delimiter is provided,
  #               provide the delimiter so it can be split

  # split the ID column
  tibble <- tibble %>%
    separate_rows(all_of(IDcol), sep = sep)

  # get the IDs from the dataset (possibly mixed ID types)
  IDs_tibble_T <- tibble(ID = str_remove(unname(unlist(tibble[, IDcol])), "-1-p1$")) %>%
    # remove trailing "-1-p1" often found in IDs from proteomics analysis
    mutate(rowN = row_number())

  # find where these IDs appear in all possible IDs from the key file
  key_tibble_sub <- key_tibble %>%
    filter(Entry %in% IDs_tibble_T$ID
           | TVAG %in% IDs_tibble_T$ID
           | TVAGG3 %in% IDs_tibble_T$ID)

  # join by every column, then paste, unique, and repeat
  keyCols <- colnames(key_tibble_sub)

  for (i in 1:length(keyCols)) {

    colnames(IDs_tibble_T)[1] <- keyCols[i] # rename column for merging

    # join
    IDs_tibble_T <- left_join(IDs_tibble_T, key_tibble_sub, by = keyCols[i])

    # remove duplicates and repeat
    ExtraID_col_pos <- c(ncol(IDs_tibble_T) - 1, ncol(IDs_tibble_T))
    colnames(IDs_tibble_T)[ExtraID_col_pos] <- c("ID1", "ID2")

    IDs_tibble_T <- IDs_tibble_T %>%
      group_by(rowN) %>%
      mutate(OtherIDs = paste0(unique(c(ID1, ID2)), collapse = ";")) %>%
      dplyr::select(1, 2, matches("OtherIDs"))

    colnames(IDs_tibble_T)[(3 + (i - 1))] <- paste0(colnames(IDs_tibble_T)[(3 + (i - 1))], i)
  }

  # rename the ID column used for joining
  colnames(IDs_tibble_T)[1] <- "oriID"

  # consolidate all other IDs
  IDs_tibble_T <- IDs_tibble_T %>%
    group_by(rowN) %>%
    mutate(OtherIDs = paste0(unique(c(unlist(str_split(OtherIDs1, ";")),
                                      unlist(str_split(OtherIDs2, ";")),
                                      unlist(str_split(OtherIDs3, ";")))
    ),
    collapse = ";"),
    OtherIDs = str_remove_all(OtherIDs, "(NA;|;NA)") # remove NAs
    ) %>%
    dplyr::select(1, 2, OtherIDs) %>%
    unique()

  # create TVAGG3 key
  IDs_tibble_T <- IDs_tibble_T %>%
    group_by(rowN) %>%
    mutate(key = paste0(unlist(str_extract_all(OtherIDs, "TVAGG3_\\d{7}")), collapse = ";")
    ) %>%
    separate_rows(key, sep = ";") %>%
    arrange(key) %>% # arrange so the first ID (lowest number) is always first
    mutate(key = if_else(key == "",
                         oriID,
                         key)
    ) %>%
    mutate(keyPos = row_number()) %>%
    filter(keyPos == 1)


  # join with the original data
  colnames(IDs_tibble_T)[1] <- colnames(tibble)[IDcol]
  tibble_res <- left_join(IDs_tibble_T, tibble)

  # create joined IDs
  return(tibble_res)

}

# Function that projects a given dataset onto a t-SNE or PCA plot
lopit_project_datasets <- function(lopit_analyse_results,
                                   otherData,
                                   IDcols,
                                   showCol,
                                   file_out) {
  # lopit_analyse_results = TSV returned by lopit_analyse that contains reduced-dimensionality coordinates (t-SNE or PCA)
  # otherData             = additional data to overlay
  # IDcols                = positions of the ID columns in the first and second dataset for merging, e.g. c(1,2)
  # showCol               = column to use for colouring; should be character or factor; must be in otherData
  # file_out              = path to output PNG

  # rename columns for easy handling
  colnames(lopit_analyse_results)[IDcols[1]] <- "ID"
  colnames(otherData)[IDcols[2]] <- "ID"
  colnames(otherData)[showCol] <- "colour"

  # join the data
  plot_data <- left_join(lopit_analyse_results, otherData, by = "ID")
  onlySelected_T <- plot_data %>%
    filter(!is.na(colour))

  ## plot the non-interactive plot ##
  tSNE_plot <- ggplot() +
    geom_jitter(data = plot_data, aes(Dim1, Dim2, color = colour), size = 2, alpha = 0.5) +
    geom_point(data = onlySelected_T, aes(Dim1, Dim2), color = "black", size = 3) +
    geom_point(data = onlySelected_T, aes(Dim1, Dim2, color = colour), size = 2) +
    theme_minimal(base_size = 20) +
    labs(color = "Localization")

  png(file_out,
      width = 4000,
      height = 2500,
      res = 300
  )

  plot(tSNE_plot)

  dev.off()
}

lopit_plot_markers <- function(pred_res,
                               DirOut,
                               quantiCols,
                               replicate_count) {
  # pred_res        = full path to the TSV output from the lopit_analyse function
  # DirOut          = directory to save marker plots
  # quantiCols      = positions of the quantification columns in the output dataset
  # replicate_count = number of replicates in the file

  # load the data
  res_T <- read_tsv(pred_res)

  # get the markers
  res_markers_T <- res_T %>%
    dplyr::filter(markers != "unknown")

  # factor the replicates
  levels_rep <- colnames(res_markers_T)[quantiCols]

  # create replicate groups
  fract_per_replicate <- length(levels_rep) / replicate_count
  replicate_groups_Sch <- paste0(levels_rep, "_x_", rep(1:replicate_count, each = fract_per_replicate))
  colnames(res_markers_T)[quantiCols] <- replicate_groups_Sch

  # plot the pattern of each marker group
  map(unique(res_markers_T$markers),
      ~{
        # remove any special characters from compartment name for use in filename
        outName <- str_remove_all(.x, "[^a-zA-Z0-9_]")

        # filter the data
        res_markers_map_T <- res_markers_T %>%
          filter(markers == .x) %>%
          gather(all_of(quantiCols), key = "fraction", value = "intensity") %>%
          mutate(intensity = as.double(intensity))

        # make fraction a factor
        res_markers_map_T$fraction <- factor(res_markers_map_T$fraction, levels = replicate_groups_Sch)

        # make group per ID per replicate
        res_markers_map_T <- res_markers_map_T %>%
          group_by(ID, fraction) %>%
          mutate(rep = str_split(fraction, "_x_", simplify = T)[, 2],
                 group = paste0(ID, "_", rep)
          )

        # plot the data
        plot <- ggplot(res_markers_map_T, aes(fraction, intensity, group = group, colour = rep)) +
          geom_point() +
          geom_line() +
          theme_minimal(base_size = 20) +
          labs(title = .x) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))

        png(paste0(DirOut, "/", outName, ".png"), width = 6000, height = 3000, res = 300)

        plot(plot)

        dev.off()

      })

}

lopit_GUI_visualisation <- function(pred_output,
                                    method = "t-SNE",
                                    colour_col,
                                    filter_colour_col = NULL,
                                    filter_cutoff = NULL,
                                    seed = 42) {
  # pred_output       = path to the prediction results from SVM or MCMC,
  #                     named SVMres_your_suffix.Rdata or MCMCres_your_suffix.Rdata
  # method            = dimensionality reduction method; options: "t-SNE", "PCA", "UMAP", "lda"
  # colour_col        = column to use for colouring;
  #                     if NULL, the type of result (MCMC or SVM) is detected automatically
  #                     and colouring is applied using the provided cutoff
  # filter_colour_col = column to use for filtering predictions (global cutoff);
  #                     note: the main function now applies per-compartment cutoffs for SVM;
  #                     TAGM-MAP and MCMC use global cutoffs
  # filter_cutoff     = probability cutoff (minimum 0.99) for MCMC, or t cutoff (minimum 0.5) for SVM

  # load the data
  load(pred_output)

  # find the variable name
  var_name <- ls()[str_detect(ls(), "res")][1]  # [1] takes first match if multiple

  # rename to a consistent variable
  assign("pred_res", get(var_name))
  rm(list = var_name)

  if (is.null(filter_colour_col) & is.null(filter_cutoff)) {
    # if NULL, colour by the provided column directly
    set.seed(seed)
    pRolocVis(object = pred_res,
              method = method,
              fcol = colour_col)
    cat("here")
  } else if (!is.null(filter_colour_col) & !is.null(filter_cutoff)) {
    # get the data
    data_T <- fData(pred_res)

    # rename the filter column
    filter_col <- which(colnames(data_T) == filter_colour_col)
    colnames(data_T)[filter_col] <- "filter_col"

    # rename the prediction column
    pred_col <- which(colnames(data_T) == colour_col)
    colnames(data_T)[pred_col] <- "pred_col"

    # create the filtered prediction column
    data_T <- data_T %>%
      mutate(pred_filter = if_else(filter_col >= filter_cutoff,
                                   pred_col,
                                   NA_character_)
      )

    # add the new filtered prediction
    fData(pred_res)$pred_filter <- data_T$pred_filter

    # plot
    set.seed(seed)
    pRolocVis(object = pred_res,
              method = method,
              fcol = "pred_filter")
  } else {
    errorCondition("Either provide both filter_colour_col and filter_cutoff, or provide only the colour column.")
  }

}

lopit_resolution <- function(pred_output,
                             marker_col = "markers",
                             DirOut) {
  # pred_output = path to prediction results from SVM or MCMC,
  #               named SVMres_your_suffix.Rdata or MCMCres_your_suffix.Rdata
  # marker_col  = name of the markers column
  # DirOut      = directory to save plots

  # load the data
  load(pred_output)

  # find the variable name
  var_name <- ls()[str_detect(ls(), "(svmres|tagm_map_res|mcmcRes)")][1]

  # rename to a consistent variable
  assign("pred_res", get(var_name))
  rm(list = var_name)

  # QSep function
  qsep_res <- QSep(object = pred_res, fcol = marker_col)

  # create the output directory
  dir.create(DirOut)

  ## Boxplot of raw QSep scores
  png(paste0(DirOut, "/", "box_plot.png"))
  plot(x = qsep_res, norm = FALSE)
  dev.off()

  ## Boxplot of normalised QSep scores
  png(paste0(DirOut, "/", "box_plot_norm.png"))
  plot(x = qsep_res, norm = TRUE)
  dev.off()

  ## Level plot of raw QSep values between each compartment
  png(paste0(DirOut, "/", "heatmap.png"))
  plot <- levelPlot(object = qsep_res, norm = FALSE)
  plot(plot)
  dev.off()

  ## Level plot of normalised QSep values between each compartment
  png(paste0(DirOut, "/", "heatmap_norm.png"))
  plot <- levelPlot(object = qsep_res, norm = TRUE)
  plot(plot)
  dev.off()

}

# Function that merges datasets similarly to lopit_merge_replicates, but importantly
# returns an MSnSet object (required for BANDLE MCMC prediction).
# All datasets must use the same ID type (e.g. UniProt)!
lopit_combined_data_for_bundle <- function(list_data_path,
                                           data_names,
                                           outDir,    # full path to output directory "my/output/dir"
                                           outSuffix  # suffix for output files
) {
  # list_data_path = list of paths to datasets to join; datasets are in MSnSet format,
  #                  as output by lopit_analyse when no prediction is run
  # data_names     = names of the loaded data objects (as named in the R environment
  #                  when loaded via load())

  data_list <- vector(mode = "list", length = length(list_data_path))

  # load the datasets
  for (i in 1:length(list_data_path)) {
    load(list_data_path[i])

    var_name <- ls()[str_detect(ls(), paste0("^", data_names[i], "$"))]

    assign(paste0("data_", i), get(var_name))
    rm(list = var_name)

    data_list[[i]] <- get(paste0("data_", i))
  }

  # create the BANDLE input by extracting features common across all datasets
  data <- commonFeatureNames(data_list)

  # save the combined input
  save(data,
       file = paste0(outDir, "/Bundle_joined_data_", outSuffix, ".Rdata")
  )

}

# Function to test different perplexity values for t-SNE
lopit_try_perplexities <- function(data,
                                   perplexities = seq(10, 100, by = 10),
                                   outDir,    # full path to output directory "my/output/dir"
                                   outSuffix, # suffix for output files
                                   seed = 42
) {
  # data = path to an MSnSet dataset prepared for prediction

  # load the data
  load(data)

  # create plots for different perplexities
  map(perplexities,
      ~{
        set.seed(seed)

        # run t-SNE
        pred_filtered_plot <- plot2D(data,
                                     method = "t-SNE",
                                     methargs = list(perplexity = .x)
        )

        # prepare data for the plot
        pred_filtered_plot_T <- tibble(ID = rownames(pred_filtered_plot),
                                       Dim1 = pred_filtered_plot[, 1],
                                       Dim2 = pred_filtered_plot[, 2])

        ## non-interactive plot ##
        tSNE_plot <- ggplot(pred_filtered_plot_T, aes(Dim1, Dim2)) +
          geom_jitter() +
          theme_minimal(base_size = 20) +
          labs(title = paste0("perplexity:", .x))

        png(paste0(outDir, "/t-SNE", "_perplexity", .x, outSuffix, ".png"),
            width = 4000,
            height = 2500,
            res = 300
        )

        plot(tSNE_plot)

        dev.off()

      })

}

# NOTE: Tested and working
lopit_mcmc_show_localization_distribution <- function(mcmcRes_converged_pooled_path,
                                                      ID,
                                                      outDir) {
  # mcmcRes_converged_pooled_path = full path to MCMC results file beginning with
  #                                  mcmcRes_converged_pooled_for_violin_plots
  # ID     = the protein ID to display
  # outDir = full path to the directory for output

  load(mcmcRes_converged_pooled_path)

  png(paste0(outDir, "/", ID, ".png"))

  plot(mcmcRes_converged_pooled, ID)

  dev.off()
}

# Function that selects points within a circle drawn on the t-SNE plot
lopit_points_in_circle <- function(data_path,
                                   center_x,
                                   center_y,
                                   radius,
                                   x_col = "Dim1",
                                   y_col = "Dim2",
                                   output_path) {
  # output_path = output path for the XLSX file and the graph (no file extension)

  # load the TSV output from lopit_analyse
  data_T <- read_tsv(data_path)

  # rename the dimension variables
  x_pos <- which(colnames(data_T) == x_col)
  y_pos <- which(colnames(data_T) == y_col)

  colnames(data_T[c(x_pos, y_pos)]) <- c("Dim1", "Dim2")

  # filter points inside the circle using vectorised distance calculation
  dist <- sqrt((data_T[[x_col]] - center_x)^2 + (data_T[[y_col]] - center_y)^2)
  inside <- data_T[dist <= radius, , drop = FALSE]

  # create plot with all points, highlighted interior points, and circle overlay
  p <- ggplot(data_T, aes(Dim1, Dim2)) +
    geom_point(alpha = 0.6, color = "gray50") +
    geom_point(data = inside, aes(Dim1, Dim2),
               color = "red", size = 2) +
    geom_circle(aes(x0 = center_x, y0 = center_y, r = radius),
                fill = NA, color = "blue", linewidth = 1, inherit.aes = FALSE) +
    coord_fixed() +  # equal aspect ratio for accurate circle representation
    labs(title = paste("Points inside circle (center:", center_x, ",", center_y,
                       "radius:", radius, ")"),
         subtitle = paste("Red points:", nrow(inside), "inside")) +
    theme_minimal()

  print(p)

  png(paste0(output_path, ".png"),
      width = 4000,
      height = 2500,
      res = 300
  )

  plot(p)

  dev.off()

  writexl::write_xlsx(inside, paste0(output_path, ".xlsx"))

  return(invisible(inside))  # returns filtered data invisibly to avoid cluttering the console
}

# MAIN LOPIT ANALYSIS FUNCTION
lopit_analyse <- function(data = NULL,     # full path to data input "this/is/my/data.xlsx";
                                            # must be an xlsx file, or an Rdata file output
                                            # by this function when data preparation was run previously
                          descLopit_data = NULL, # if providing pre-prepared Rdata, also provide
                                                  # the path to the descLopit_data TSV created during preparation
                          markers = NULL,  # full path to the markers file (Excel);
                                            # format: 3 columns named "ID", "markers", "colors"
                                            # (first: ID, second: marker localisation, third: coloured cells);
                                            # keep compartment names consistent!
                                            # small complexes should not be used as markers but can be
                                            # overlaid on the reduced-dimensionality plot later to check
                                            # whether they co-localise as expected
                          MS = "TNT",      # "TNT" or "LFQ"; use LFQ if combining LFQ and TMT data
                          filterRows = 18, # if MS is LFQ: rows with >= this many missing values are removed
                          minTNTval = 50,  # for TNT: intensities below this value are treated as NA
                          transform = NULL, # "log2" or "scale"
                          scale = c(0, 1000), # if transform = "scale": lower and upper bounds for rescaling
                          impute = F,      # whether to impute missing values;
                                            # if FALSE, zeros replace NAs;
                                            # if TRUE, imputation uses provided mean and SD
                          imputeMean = NULL, # mean for imputation normal distribution;
                                              # if not provided, defaults to mean shifted by -1.8 SD
                          imputeSd = NULL,   # SD for imputation normal distribution;
                                              # if not provided, defaults to 0.3 * data SD
                          outDir,          # full path to output directory "my/output/dir"
                          outSuffix,       # suffix for output files
                          IDsCol,          # position of the ID column in the data
                          quantiCols,      # positions of the intensity columns: c(5:15) or c(1,2,3,...)
                          quantiColsGroups = NULL, # if multiple datasets are used, group them for
                                                    # per-dataset scaling rather than global scaling;
                                                    # e.g. three datasets with 3 columns each: c(1,1,1,2,2,2,3,3,3)
                          descCol,         # position of the protein description column
                                            # (added to the interactive graph)
                          metaCols,        # additional metadata columns to include in the output
                          IDs_sep = ";",   # character separating IDs in the ID column
                          # normalisation is recommended when combining datasets; improves results
                          norm = NULL,     # options: "sum", "max", "quantiles", "quantiles.robust",
                                            # "center.mean", "center.median", "div.mean", "div.median",
                                            # "diff.median", "vsn"
                          predType = "SVM", # prediction method: "SVM", "MCMC", "TAGM-MAP" (supervised),
                                             # "PHEATMAP", "HDBSCAN" (unsupervised), or "NO" (prepare data only)

                          ##########################
                          ## Unsupervised methods ##
                          ##########################
                          # Variables for pheatmap clustering
                          # Results should not exceed ~15 clusters; more compartments become indistinguishable
                          max_cluster_size = 500, # maximum number of proteins per cluster
                          iteration = 100,        # iterations through the dendrogram;
                                                   # increase if the function reports this is too low

                          # Variables for HDBSCAN clustering
                          minPts = 10,    # typically 10-15

                          #######################
                          ## Supervised methods ##
                          #######################
                          # Variables for SVM
                          sigma = NULL,   # SVM kernel parameter; estimated if not provided
                          cost = NULL,    # SVM cost parameter; estimated if not provided
                          t_cutoff = 0.75, # SVM score cutoff to consider a protein as predicted
                                           # to a given compartment (must exceed this to be assigned)
                          t_cutoff_per_compartment = T, # if TRUE, compartment-specific t thresholds
                                                         # are used; t_cutoff then means quantile cutoff,
                                                         # not an absolute t value

                          # Variables for TAGM-MAP
                          numIter_tagmmap = 100,  # number of posterior sampling steps
                          p_cutoff_tagmmap = 0.99, # mean probability threshold for colouring predictions;
                                                    # probabilities below this are shown in grey

                          # Variables for MCMC
                          # Class imbalance can cause over-classification;
                          # partly addressed by the overall_probability composite score
                          # (product of localisation probability, outlier probability, and Shannon entropy)
                          numIter = 20000,  # number of posterior sampling steps
                          burnin = 10000,   # number of initial steps to discard
                                             # (the chain needs time to reach the posterior from a random start)
                          thin = 20,        # thinning interval to reduce autocorrelation
                                             # (e.g. 20 = keep every 20th step)
                          numChains = 6,    # number of chains to run
                          MCMCres_path = NULL, # path to existing MCMC results to skip re-running;
                                               # must also provide path to the prepared data Rdata file
                          p_cutoff = 0.99,  # overall probability threshold for predictions
                          hyppar = c(0.5, 3, 100), # BANDLE: priors for the hyper parameters
                          dirPrior = 0.0005,        # BANDLE: Dirichlet prior matrix;
                                                     # reflects expected probability that a protein
                                                     # has different localisation between datasets
                          conditions, # BANDLE: which datasets correspond to which condition;
                                       # same format as quantiColsGroups (e.g. c(1,1,2,2))

                          ###########################################################
                          # Chain selection strategy                                #
                          ###########################################################
                          # Options for MCMCres_keepChains:
                          # "geweke"           - keep chains passing Geweke test
                          # "gelman"           - keep chains passing Gelman test
                          # "geweke_and_gelman" - keep chains passing both
                          # c(1,2,3)           - manually specify chain numbers to keep

                          MCMCres_keepChains = "gelman",

                          extra_thin = NULL,    # MCMC only (not BANDLE): additional thinning of posterior samples
                          extra_burnin = NULL,  # MCMC only (not BANDLE): additional steps removed from chain start
                          geweke_pval = 0.05,   # Geweke test significance threshold; tests whether the first 10%
                                                 # of the chain is statistically equivalent to the last 50%;
                                                 # significant result means the chain did not converge

                          ##############################
                          # Visualisation parameters   #
                          ##############################
                          seed = 42,        # seed for reproducible t-SNE or PCA plots
                          dimRed = "t-SNE", # dimensionality reduction method: "PCA" or "t-SNE"
                          perplexity = 30,  # t-SNE perplexity parameter; roughly the number of
                                             # effective neighbours per point;
                                             # good practice: try several values (10, 20, 30, 40, 50)
                          dims = c(1, 2)    # which PCs to use for PCA (default: first two)
) {

  if (str_detect(data, pattern = "\\.xlsx")) {
    cat("Excel file provided, running the analysis from scratch.\n")
    ####################################
    #### change inputs and set seed ####
    ####################################
    # remove trailing backslash from path if present
    outDir <- str_remove(outDir, "/$")
    transform <- ifelse(is.null(transform),
                        "no",
                        transform
    )

    #### read in the data ####
    lopit_T <- readxl::read_xlsx(data) %>%
      mutate_at(quantiCols, as.double)

    #########################################
    #### separate the ID column by IDs_sep ####
    #########################################
    colnames(lopit_T)[IDsCol] <- "ID"

    lopit_T <- lopit_T %>%
      mutate(rowN = row_number()) %>%
      separate_rows(ID, sep = IDs_sep)

    # create description table
    descLopit_data <- lopit_T[, c(IDsCol, descCol, ncol(lopit_T))]

    # create intensities table
    intensities <- lopit_T[, quantiCols]

    #### for TNT data, replace low-intensity values with NA
    if (MS == "TNT") {
      intensities <- intensities %>%
        mutate_all(~if_else(. <= minTNTval,
                            NA_real_,
                            .)
        )
    }

    # count NAs per row and add to description
    NAs_n <- apply(intensities, 1, sumNAs)
    descLopit_data$NAs_n <- NAs_n

    colnames(descLopit_data) <- c("ID", "desc", "rowN", "NAs_n")

    ########################################
    #### keep only the first five words of the description ####
    ########################################
    keep5_words <- function(vec) {
      fristdesc <- str_split(vec, ";", simplify = T)[, 1]
      words <- unlist(str_split(fristdesc, "\\s"))

      if (length(words) > 5) {
        words <- words[1:5]
      }

      words <- paste0(words, collapse = " ")

      return(words)
    }

    descLopit_data$desc2 <- lapply(descLopit_data$desc, keep5_words) %>%
      unlist()

    # restore group IDs in the description dataset
    descLopit_data <- descLopit_data %>%
      group_by(rowN) %>%
      mutate(ID = paste0(ID, collapse = IDs_sep)) %>%
      sample_n(size = 1)

    write_tsv(descLopit_data, paste0(outDir, "/", "data_for_prediction_descLopit_data", outSuffix, ".tsv"))

    #####################################################################
    ### filter out rows with too many NAs (controlled by filterRows) ####
    #####################################################################
    lopit_T <- lopit_T[NAs_n < filterRows, ]

    if (transform == "log2") {
      ##################################
      ## log2-transform the data     ##
      ##################################
      lopit_T <- lopit_T %>%
        mutate_at(quantiCols, ~if_else(is.na(.) | . == 0,
                                       NA_real_,
                                       log2(.))
        )

      cat("Data were log2 transformed.\n")
    } else if (transform == "scale" & length(scale) == 2 & all(is.numeric(scale))) {

      ########################################
      ## row-wise rescaling                 ##
      ########################################
      intensities <- lopit_T[, quantiCols]

      if (is_null(quantiColsGroups)) {
        # if no groups provided, scale all intensities together
        scaled_intensities <- apply(intensities, 1, rescale_to_xx, scale) %>%
          t() %>%
          as_tibble()

      } else if (typeof(quantiColsGroups) == "double") {
        # if groups provided, scale each group separately
        quantiColsGroups_L <- map(unique(quantiColsGroups),
                                  ~{
                                    which(quantiColsGroups == .x)
                                  })

        scaled_intensities <- map(quantiColsGroups_L,
                                  ~{
                                    scaled_map_T <- apply(intensities[, .x], 1, rescale_to_xx, scale) %>%
                                      t() %>%
                                      as_tibble()

                                    return(scaled_map_T)
                                  }) %>%
          purrr::reduce(cbind)

      } else {
        errorCondition("Wrong input for quantiColsGroups.\nMust be NULL or a vector of numbers splitting the intensities into groups:\ne.g. c(1,1,1,2,2,2,3,3,3)")
      }

      lopit_T[, quantiCols] <- scaled_intensities

      cat("Data were scaled/transformed.\n")
    } else if (transform == "no") {
      cat("Data were NOT transformed.\n")
    } else {
      stop('Wrong input for transformation. Possible inputs are NULL, \"log2\" or \"scale\".\nIf using \"scale\", the scale argument must be a vector of 2 numbers.')
    }

    ################
    ## Imputation ##
    ################
    if (impute == F) {
      # replace NAs with 0
      lopit_T <- lopit_T %>%
        mutate_at(quantiCols, ~if_else(is.na(.),
                                       0,
                                       .),
        )

      cat("Zeros were placed instead of NAs.\n")
    } else if (impute == T & !is.null(imputeMean) & !is.null(imputeSd)) {
      ## impute with a normal distribution using user-specified parameters
      lopit_T <- lopit_T %>%
        group_by(ID) %>%
        mutate_at(quantiCols, ~if_else(is.na(.),
                                       rnorm(1, imputeMean, imputeSd),
                                       .)
        )

      cat("Missing values were imputed with a normal distribution using user-specified parameters.\n")
    } else if (impute == T & xor(is.null(imputeMean), is.null(imputeSd))) {
      stop("You provided only mean or SD for imputation. Provide both and rerun the function!\n")
    } else if (impute == T & is.null(imputeMean) & is.null(imputeSd)) {

      ########################################################################
      ## impute with default parameters:                                    ##
      ## mean shifted by -1.8 SD of the data distribution;                 ##
      ## SD = 0.3 * data SD                                                 ##
      ########################################################################
      meanData <- mean(unlist(lopit_T[, quantiCols]), na.rm = T)
      sdData <- sd(unlist(lopit_T[, quantiCols]), na.rm = T)

      lopit_T <- lopit_T %>%
        group_by(ID) %>%
        mutate_at(quantiCols, ~if_else(is.na(.),
                                       rnorm(1, (meanData - 1.8 * sdData), sdData * 0.3),
                                       .)
        )
      cat("Missing values were imputed with a normal distribution using default parameters:\n--mean shifted by -1.8 SD of the data\n--SD = 0.3 * data SD.\n")
    }

    #####################
    ## Add the markers ##
    #####################
    if (str_detect(predType, "(HDBSCAN|PHEATMAP|NO)")) {
      # if no prediction is requested, create mock markers
      markers_T <- tibble(ID = sample(lopit_T$ID, size = 100, replace = F),
                          markers = "random") %>%
        unique()

    } else if (str_detect(predType, "(MCMC|SVM|TAGM-MAP)")) {
      markers_T <- readxl::read_xlsx(markers)

      ########################################################
      ## check for colour column in the marker file (col 3) ##
      ########################################################
      if (any(is.na(colnames(markers_T)[3] != "colors") | colnames(markers_T)[3] != "colors", na.rm = T)) {
        # create levels and colours if no colors column is provided
        levels_T <- markers_T %>%
          dplyr::select(markers) %>%
          unique()

        levels <- c(levels_T$markers, "unknown")
        colors <- rainbow(length(levels) - 1)
        colors <- c(colors, "grey70")
      } else {
        cat("Using user-provided colours from the third column of the marker list!\n")

        ## read colours from the Excel cell fill colour
        markers_job <- loadWorkbook(markers)
        markers_sheet1 <- getSheets(markers_job)[[1]]
        rows <- getRows(markers_sheet1)
        cells <- getCells(rows)
        styles <- sapply(cells, getCellStyle)
        colours_cols_L <- styles[str_detect(names(styles), "\\.3$")]
        colours_cols_L <- colours_cols_L[-1]
        colors_sch <- sapply(colours_cols_L, cell_color)

        # add to markers
        markers_T$colors <- colors_sch

        # create levels and associated colours
        levels_T <- markers_T %>%
          dplyr::select(markers, colors) %>%
          unique()

        levels <- c(levels_T$markers, "unknown")
        colors <- unname(c(levels_T$colors, "grey70"))

      }

    } else {
      errorCondition("Wrong predType provided. Options are: \"MCMC\", \"SVM\", \"TAGM-MAP\", \"HDBSCAN\", or \"PHEATMAP\"!\n")

    }

    # add the markers
    lopit_T <- left_join(lopit_T, markers_T, by = "ID")
    ################################################################################
    # ensure that group members inherit the marker assignment of any labelled row  #
    ################################################################################
    lopit_T <- lopit_T %>%
      group_by(rowN) %>%
      mutate(markers = paste0(unique(markers[!is.na(markers)]), collapse = ","))

    # return an error if any markers span more than two compartments
    multimarker_protein <- unique(str_count(lopit_T$markers, ","))

    if (length(multimarker_protein) > 1) {
      dupID <- markers_T$ID[duplicated(markers_T$ID)]
      cat("These markers are duplicated:\n",
          paste0(dupID,
                 collapse = "\n"),
          "\n",
          sep = "")
      stop("Error: Some markers have more than two compartments. Check your marker list!\n")
    }

    # report markers not found in the dataset
    present_marker <- !is.na(lopit_T$markers)
    present_marker <- lopit_T$ID[present_marker]

    cat("These markers were not found in your dataset:\n",
        paste0(markers_T$ID[!(markers_T$ID %in% present_marker)],
               collapse = "\n"),
        "\n",
        sep = "")

    ## replace empty strings with "unknown"
    lopit_T$markers <- if_else(lopit_T$markers == "",
                               "unknown",
                               lopit_T$markers)

    ## restore grouped IDs to their original form
    lopit_T <- lopit_T %>%
      group_by(rowN) %>%
      mutate(ID = paste0(ID, collapse = IDs_sep)) %>%
      sample_n(size = 1)

    ## extract quantification columns
    expression_T <- lopit_T[, c(IDsCol, quantiCols)]
    colnames(expression_T) <- c("ID", paste0("quan", 1:length(quantiCols)))

    ## extract metadata columns
    meta_T <- lopit_T[, c(IDsCol, metaCols)]
    meta_T$markers <- lopit_T$markers

    ## create fraction metadata table
    fraction_T <- tibble(sampleNames = colnames(expression_T)[2:ncol(expression_T)],
                         Fractions = 1:length(colnames(expression_T)[2:ncol(expression_T)]))

    ## write TSVs for readMSnSet
    outQuan <- paste0(outDir, "/", "Quan_data", outSuffix, ".tsv")
    outMeta <- paste0(outDir, "/", "Meta_data", outSuffix, ".tsv")
    outFrac <- paste0(outDir, "/", "Frac_data", outSuffix, ".tsv")

    write_tsv(expression_T, outQuan)
    write_tsv(meta_T, outMeta)
    write_tsv(fraction_T, outFrac)

    ## read data into the MSnSet structure
    data <- readMSnSet(exprsFile = outQuan,
                       featureDataFile = outMeta,
                       phenoDataFile = outFrac,
                       sep = "\t")

    #######################
    #### Normalisation  ####
    #######################
    if (!is.null(norm)) {
      data <- normalise(data, method = norm)

      cat("Normalisation by", norm, "was performed.\n")
    } else {
      cat("Normalisation was NOT performed.\n")
    }

    ############################
    ## save the prepared data ##
    ############################
    save(data, file = paste0(outDir, "/", "data_for_prediction", outSuffix, ".Rdata"))

  } else if (str_detect(data, "\\.Rdata") & !is.null(descLopit_data)) {
    cat("Pre-prepared data and descLopit were provided; loading the data.\n")
    load(data)
    descLopit_data <- read_tsv(descLopit_data)
  } else {
    errorCondition("Neither an Excel file nor both data Rdata and descLopit_data were provided.\n")
  }

  ###########################
  ## Run selected prediction ##
  ###########################
  if (predType == "SVM") {
    cat("Running SVM predictions.\n")

    run_SVM(data = data,
            sigma = sigma,
            cost = cost,
            seed = seed,
            t_cutoff = t_cutoff,
            t_cutoff_per_compartment = t_cutoff_per_compartment,
            dimRed = dimRed,
            perplexity = perplexity,
            dims = dims,
            levels = levels,
            colors = colors,
            outDir = outDir,
            outSuffix = outSuffix,
            descLopit_data = descLopit_data)

  } else if (predType == "TAGM-MAP") {

    run_TAGMMAP(data = data,
                numIter_tagmmap = numIter_tagmmap,
                seed = seed,
                p_cutoff_tagmmap = p_cutoff_tagmmap,
                dimRed = dimRed,
                perplexity = perplexity,
                dims = dims,
                levels = levels,
                colors = colors,
                outDir = outDir,
                outSuffix = outSuffix,
                descLopit_data = descLopit_data)

  } else if (predType == "MCMC") {
    cat("Running MCMC predictions.\n")

    run_MCMC(data = data,
             numIter = numIter,
             burnin = burnin,
             thin = thin,
             numChains = numChains,
             MCMCres_path = MCMCres_path,
             MCMCres_keepChains = MCMCres_keepChains,
             extra_thin = extra_thin,
             extra_burnin = extra_burnin,
             geweke_pval = geweke_pval,
             seed = seed,
             dimRed = dimRed,
             perplexity = perplexity,
             dims = dims,
             levels = levels,
             colors = colors,
             p_cutoff = p_cutoff,
             outDir = outDir,
             outSuffix = outSuffix,
             descLopit_data = descLopit_data)

  } else if (predType == "PHEATMAP") {
    cat("Running unsupervised clustering with hclust via pheatmap.\n")
    run_pheatmap_clustering(data = data,
                            descLopit_data = descLopit_data,
                            max_cluster_size = max_cluster_size,
                            iteration = iteration,
                            outDir = outDir,
                            outSuffix = outSuffix,
                            seed = seed,
                            dimRed = dimRed,
                            perplexity = perplexity,
                            dims = dims
    )

  } else if (predType == "HDBSCAN") {
    cat("Running unsupervised clustering with HDBSCAN.\n")
    run_hdbscan_clustering(data = data,
                           descLopit_data = descLopit_data,
                           minPts = minPts,
                           outDir = outDir,
                           outSuffix = outSuffix,
                           seed = seed,
                           dimRed = dimRed,
                           perplexity = perplexity,
                           dims = dims
    )

  } else if (predType == "NO") {
    cat("No predType was selected. Not performing any further analysis!\n")

  } else {
    errorCondition("Wrong option for predType. Options are:\n\"SVM\", \"MCMC\", \"TAGM-MAP\" (supervised), \"PHEATMAP\", \"HDBSCAN\" (unsupervised), or \"NO\" (prepare data only)")
  }
}

###############################################
## Helper functions called within lopit_analyse
###############################################

# Function to extract fill colour from an Excel cell style
cell_color <- function(style) {
  fg <- style$getFillForegroundXSSFColor()
  hex <- tryCatch(fg$getRgb(), error = function(e) NULL)
  hex <- paste0("#", paste(hex, collapse = ""))
  tint <- tryCatch(fg$getTint(), error = function(e) NULL)

  if (!is.null(tint) & !is.null(hex)) {
    rgb_col <- col2rgb(col = hex)
    if (tint < 0) rgb_col <- (1 - abs(tint)) * rgb_col
    if (tint > 0) rgb_col <- rgb_col + (255 - rgb_col) * tint
    hex <- rgb(red = rgb_col[1, 1], green = rgb_col[2, 1], blue = rgb_col[3, 1], maxColorValue = 255)
  }
  return(hex)
}

# Function for row-wise rescaling
rescale_to_xx <- function(vec, scale) {
  rescale(vec,
          from = c(min(vec, na.rm = T),
                   max(vec, na.rm = T)
          ),
          to = scale
  )
}

# Function to count NAs in a vector
sumNAs <- function(vec) {
  sum(is.na(vec))
}

# Function to run Gelman diagnostics on MCMC chains
run_gelman <- function(inputMatrices,
                       names,
                       subset = NULL
) {
  # inputMatrices = list of results from the mcmc_get_ functions
  # names         = names of the matrices in a vector
  # subset        = a vector of chain indices to subset from the mcmc_get_ results

  gelman_diag_res_T <- map2(inputMatrices,
                            names,
                            ~{
                              if (is.null(subset)) {
                                subset <- 1:length(inputMatrices[[1]])
                              }

                              gelman_diag_map_res <- gelman.diag(.x[subset], transform = FALSE)
                              gelman_diag_map_res_T <- as_tibble(gelman_diag_map_res[1]$psrf)

                              gelman_diag_map_res_T$data <- .y

                              return(gelman_diag_map_res_T)
                            }) %>%
    purrr::reduce(rbind)

  colnames(gelman_diag_res_T) <- str_remove_all(colnames(gelman_diag_res_T),
                                                 "\\.") %>%
    str_replace_all("\\s", "_")

  return(gelman_diag_res_T)
}

# Creates diagnostic summaries and trace plots for a given MCMC matrix
make_MCMC_stat_plots <- function(matrix,
                                 output,
                                 nChains
) {
  # output  = output file path without extension
  # matrix  = matrix returned by mcmc_get_ functions
  # nChains = number of chains that were run for the MCMC prediction

  # create a summary for each chain and save it
  for (i in seq_len(nChains)) {
    sink(paste0(output, ".txt"))
    print(summary(matrix[[i]]))
    sink()
  }

  # create trace plots and histograms
  pdf(file = paste0(output, ".pdf"))

  for (i in seq_len(nChains)) {
    plot(matrix[[i]], main = paste("Chain", i), auto.layout = FALSE, col = i)
  }

  dev.off()
}

make_NA_plot <- function(pred_filtered_plot_T,
                         outDir,
                         dimRed,
                         outSuffix) {
  # pred_filtered_plot_T = data created in the SVM or MCMC function

  ## non-interactive NA count plot ##
  tSNE_plot_NAs <- ggplot(pred_filtered_plot_T, aes(colour = NAs_n)) +
    geom_jitter(data = pred_filtered_plot_T, aes(Dim1, Dim2), size = 2, alpha = 0.5) +
    theme_minimal(base_size = 20) +
    labs(color = "Number of NAs")

  png(paste0(outDir, "/", dimRed, "_NAplot_", outSuffix, ".png"),
      width = 4000,
      height = 2500,
      res = 300
  )

  plot(tSNE_plot_NAs)

  dev.off()

}

# Function that clusters using HDBSCAN
run_hdbscan_clustering <- function(data,
                                   descLopit_data,
                                   minPts,
                                   outDir,
                                   outSuffix,
                                   seed,
                                   dimRed,
                                   perplexity,
                                   dims
) {
  ## Run HDBSCAN with the specified minimum cluster size
  hdb_results <- data %>%
    exprs() %>%
    hdbscan(minPts = minPts)  # BUG FIX: was hardcoded to 10; now uses the minPts argument

  ## Add HDBSCAN results to the MSnSet
  fData(data)$hdb_cluster_id <- hdb_results$cluster
  fData(data)$hdb_cluster_prob <- hdb_results$membership_prob

  ## Check how many proteins are in each cluster
  data %>%
    fData() %>%
    pull(hdb_cluster_id) %>%
    table()

  ## run dimensionality reduction
  set.seed(seed)

  if (dimRed == "t-SNE") {
    pred_filtered_plot <- plot2D(data,
                                 fcol = "hdb_cluster_id",
                                 method = "t-SNE",
                                 methargs = list(perplexity = perplexity)
    )
  } else if (dimRed == "PCA") {
    pred_filtered_plot <- plot2D(data,
                                 fcol = "hdb_cluster_id",
                                 dims = dims,
                                 method = "PCA")
  } else {
    stop('Error: no dimensionality reduction method provided. Choose "PCA" or "t-SNE" for the dimRed argument.\n')
  }

  pred_filtered_plot_T <- tibble(ID = rownames(pred_filtered_plot),
                                 Dim1 = pred_filtered_plot[, 1],
                                 Dim2 = pred_filtered_plot[, 2])

  # add HDBSCAN cluster assignments
  hdbscan_clusters_T <- fData(data) %>%
    dplyr::select(hdb_cluster_id, hdb_cluster_prob) %>%
    mutate(ID = rownames(.))

  ## add cluster columns to the coordinates ##
  pred_filtered_plot_T <- left_join(pred_filtered_plot_T, hdbscan_clusters_T, by = "ID")

  ## add protein descriptions ##
  pred_filtered_plot_T <- left_join(pred_filtered_plot_T, descLopit_data, by = "ID")

  ## create interactive label ##
  pred_filtered_plot_T <- pred_filtered_plot_T %>%
    group_by(ID) %>%
    mutate(label = paste0(ID, ", dbscan_cluster: ", hdb_cluster_id, " - ", desc2, collapse = ""))

  ## convert cluster IDs to factor ##
  pred_filtered_plot_T$hdb_cluster_id <- factor(pred_filtered_plot_T$hdb_cluster_id, levels = 1:length(unique(pred_filtered_plot_T$hdb_cluster_id)))

  ## Generate distinct colours ##
  colors <- rainbow(length(unique(pred_filtered_plot_T$hdb_cluster_id)))

  ## non-interactive plot ##
  tSNE_plot <- ggplot(pred_filtered_plot_T, aes(label = label)) +
    geom_jitter(data = pred_filtered_plot_T, aes(Dim1, Dim2, color = hdb_cluster_id), size = 2, alpha = 0.5) +
    theme_minimal(base_size = 20) +
    labs(color = "Localization") +
    scale_color_manual(values = colors)

  png(paste0(outDir, "/", dimRed, outSuffix, "_hdbscan_cluster.png"),
      width = 4000,
      height = 2500,
      res = 300
  )

  plot(tSNE_plot)

  dev.off()

  write_tsv(pred_filtered_plot_T,
            paste0(outDir,
                   "/",
                   dimRed,
                   "_plot_data_hdbscan_cluster",
                   outSuffix,
                   ".tsv")
  )

}

# Function that takes a pheatmap object and returns clusters plotted on a dimensionality reduction plot
run_pheatmap_clustering <- function(data,
                                    descLopit_data,
                                    max_cluster_size = 250,
                                    iteration,
                                    outDir,
                                    outSuffix,
                                    seed,
                                    dimRed,
                                    perplexity,
                                    dims
) {
  # max_cluster_size = maximum number of proteins allowed in a cluster
  # iteration        = number of steps through the dendrogram to extract clusters;
  #                    increase if a higher value is needed

  #################################
  ## Cluster using pheatmap      ##
  #################################
  pheatmap_res <- pheatmap(mat = exprs(data),
                           cluster_rows = TRUE,
                           cluster_cols = FALSE,
                           show_rownames = FALSE,
                           filename = paste0(outDir, "/", dimRed, outSuffix, "_heatmap.pdf")
  )

  # extract clusters
  row_hcluts <- pheatmap_res

  # extract ordered IDs
  row_ordered_IDs <- row_hcluts$tree_row$labels

  # create the initial single-cluster assignment
  cluster_T <- tibble(ID = row_ordered_IDs,
                      cluster1 = 1
  )

  # iterate over increasing numbers of clusters and retain those exceeding max_cluster_size
  otherClusters_T <- map(2:iteration,
                         ~{
                           row_clusters <- cutree(row_hcluts$tree_row, k = .x)

                           clusters_T <- tibble(ID = row_ordered_IDs,
                                                cluster = row_clusters) %>%
                             group_by(cluster) %>%
                             mutate(keep = n() > max_cluster_size) %>%
                             filter(keep) %>%
                             dplyr::select(-keep)

                           colnames(clusters_T)[2] <- paste0(colnames(clusters_T)[2], "_", .x)

                           return(clusters_T)
                         }) %>%
    purrr::reduce(full_join, by = "ID")

  # check whether the last iteration has resolved all large clusters
  all_is_NA <- nrow(otherClusters_T) - sum(is.na(otherClusters_T[, ncol(otherClusters_T[, ])]))

  if (all_is_NA != 0) {
    cat("Higher iteration is needed to resolve all groups!\n")
    return(NULL)
  }

  all_clusters_T <- left_join(cluster_T, otherClusters_T, by = "ID")

  # assign final cluster based on number of NAs (used as a proxy for cut depth)
  all_clusters_T$total_clusters <- apply(all_clusters_T[2:ncol(all_clusters_T)], 1, sumNAs)

  # renumber clusters into sequential integers
  cluster_heatmap_T <- all_clusters_T %>%
    dplyr::select(total_clusters) %>%
    as.data.frame()

  factor_clusters <- cluster_heatmap_T %>%
    unique() %>%
    arrange(total_clusters)

  factor_clusters$new_cluster <- 1:nrow(factor_clusters)
  cluster_heatmap_T <- left_join(cluster_heatmap_T, factor_clusters, by = "total_clusters")

  cluster_heatmap_T <- cluster_heatmap_T %>%
    dplyr::select(new_cluster) %>%
    mutate(new_cluster = as.character(new_cluster))

  rownames(cluster_heatmap_T) <- all_clusters_T$ID

  pheatmap_res <- pheatmap(mat = exprs(data),
                           cluster_rows = TRUE,
                           cluster_cols = FALSE,
                           show_rownames = FALSE,
                           annotation_row = cluster_heatmap_T,
                           filename = paste0(outDir, "/", dimRed, outSuffix, "_heatmap_clusters.pdf")
  )

  # add IDs as a column for joining to coordinates
  cluster_heatmap_T$ID <- rownames(cluster_heatmap_T)

  ###################################################################
  #### Dimensionality reduction and plot                          ####
  ###################################################################
  set.seed(seed)

  if (dimRed == "t-SNE") {
    pred_filtered_plot <- plot2D(data,
                                 method = "t-SNE",
                                 methargs = list(perplexity = perplexity)
    )
  } else if (dimRed == "PCA") {
    pred_filtered_plot <- plot2D(data,
                                 dims = dims,
                                 method = "PCA")
  } else {
    stop('Error: no dimensionality reduction method provided. Choose "PCA" or "t-SNE" for the dimRed argument.\n')
  }

  pred_filtered_plot_T <- tibble(ID = rownames(pred_filtered_plot),
                                 Dim1 = pred_filtered_plot[, 1],
                                 Dim2 = pred_filtered_plot[, 2])

  ## add cluster columns to coordinates ##
  pred_filtered_plot_T <- left_join(pred_filtered_plot_T, cluster_heatmap_T, by = "ID")

  ## add protein descriptions ##
  pred_filtered_plot_T <- left_join(pred_filtered_plot_T, descLopit_data, by = "ID")

  ## create interactive label ##
  pred_filtered_plot_T <- pred_filtered_plot_T %>%
    group_by(ID) %>%
    mutate(label = paste0(ID, ", pHeatMap_cluster: ", new_cluster, " - ", desc2, collapse = ""))

  ## convert cluster ID to factor ##
  pred_filtered_plot_T$new_cluster <- factor(pred_filtered_plot_T$new_cluster,
                                             levels = as.character(1:length(unique(pred_filtered_plot_T$new_cluster))))

  ## Generate distinct colours ##
  colors <- rainbow(length(unique(pred_filtered_plot_T$new_cluster)))

  ## non-interactive plot ##
  tSNE_plot <- ggplot(pred_filtered_plot_T, aes(label = label)) +
    geom_jitter(data = pred_filtered_plot_T, aes(Dim1, Dim2, color = new_cluster), size = 2, alpha = 0.5) +
    theme_minimal(base_size = 20) +
    labs(color = "Localization") +
    scale_color_manual(values = colors)

  png(paste0(outDir, "/", dimRed, outSuffix, "_pheatmap_cluster.png"),
      width = 4000,
      height = 2500,
      res = 300
  )

  plot(tSNE_plot)

  dev.off()

  write_tsv(pred_filtered_plot_T,
            paste0(outDir,
                   "/",
                   dimRed,
                   "_plot_data_pheatmap_cluster",
                   outSuffix,
                   ".tsv")
  )

}

# SVM prediction sub-function
run_SVM <- function(data,
                    sigma,
                    cost,
                    seed,
                    t_cutoff,
                    t_cutoff_per_compartment,
                    dimRed,
                    perplexity,
                    dims,
                    levels,
                    colors,
                    outDir,
                    outSuffix,
                    descLopit_data
) {

  #######################################################
  #### Predict protein localisations using SVM        ####
  #######################################################
  if (is.null(sigma) & is.null(cost)) {
    cat("Sigma and cost not provided. Performing optimisation then predicting localisations.\n")

    ## Get markers and calculate inverse-frequency class weights
    ## (prevents bias towards compartments with many markers)
    marker_tbl <- data %>%
      getMarkers() %>%
      table()

    weights <- 1 / marker_tbl[names(marker_tbl) != "unknown"]

    ## optimisation ##
    params <- svmOptimisation(data,
                              fcol = "markers",
                              class.weights = weights)

    ## optimisation plot: distribution of F1 scores ##
    png(paste0(outDir, "/", "Dist_F1score_best_sigma_cost", outSuffix, ".png"),
        width = 1500,
        height = 1500,
        res = 300)

    plot1 <- plot(params)
    print(plot1)

    dev.off()

    ## optimisation plot: averaged F1 scores across sigma/cost grid ##
    png(paste0(outDir, "/", "Avereged_F1score_all_sigma_cost.", outSuffix, ".png"),
        width = 1500,
        height = 1500,
        res = 300)

    plot2 <- levelPlot(params)
    print(plot2)

    dev.off()

    ## save optimisation results ##
    save(params, file = paste0(outDir, "/", "paraOptimization", outSuffix, ".Rdata"))

    ## classify using best sigma and cost ##
    svmres <- svmClassification(data, fcol = "markers",
                                assessRes = params)

  } else if (!is.null(sigma) & !is.null(cost)) {
    cat("Using user-specified sigma and cost to predict localisations.\n")

    marker_tbl <- data %>%
      getMarkers() %>%
      table()

    weights <- 1 / marker_tbl[names(marker_tbl) != "unknown"]

    svmres <- svmClassification(data,
                                fcol = "markers",
                                sigma = sigma,
                                cost = cost,
                                class.weights = weights)

  } else {
    stop("Error: either sigma or cost was not provided.")
  }

  # save SVM results
  save(svmres, file = paste0(outDir, "/SVMres", outSuffix, ".Rdata"))

  ## print processing summary ##
  proccessing <- processingData(svmres)
  cat("Data processing summary:", proccessing@processing, sep = "\n")

  #######################
  #### Visualisation  ####
  #######################
  ## distribution of SVM scores relative to cutoff ##
  pred_T <- fData(svmres)

  plotF1_dist <- ggplot(pred_T, aes(svm, svm.scores)) +
    geom_boxplot() +
    geom_hline(yintercept = t_cutoff, color = "red", linetype = 2) +
    theme_minimal(base_size = 10) +
    coord_flip() +
    ylab("prediction") +
    xlab("prediction F1 score")

  png(paste0(outDir, "/", "Dist_F1score_prediction_with_cut_off", outSuffix, ".png"),
      width = 1500,
      height = 1500,
      res = 300)

  plot(plotF1_dist)

  dev.off()

  cat("Creating graphical outputs.\n")

  ## apply cutoff to get final predictions ##
  if (t_cutoff_per_compartment == F) {
    # single global cutoff for all compartments
    pred_filtered <- getPredictions(svmres,
                                    fcol = "svm",
                                    t = t_cutoff)
    pred_filtered_T <- fData(pred_filtered)
    pred_filtered_T$ID <- rownames(pred_filtered_T)

  } else if (t_cutoff_per_compartment == T) {

    ## per-compartment quantile-based SVM score thresholds
    score_thresholds <- orgQuants(object = svmres,
                                  fcol = "svm",
                                  scol = "svm.scores",
                                  mcol = "markers",
                                  t = t_cutoff)

    pred_filtered <- getPredictions(object = svmres,
                                    fcol = "svm",
                                    scol = "svm.scores",
                                    mcol = "markers",
                                    t = score_thresholds)
    pred_filtered_T <- fData(pred_filtered)
    pred_filtered_T$ID <- rownames(pred_filtered_T)

  } else {
    errorCondition("t_cutoff_per_compartment must be TRUE or FALSE!\n")
  }


  ###################################################################
  #### Dimensionality reduction and plot                          ####
  ###################################################################
  set.seed(seed)

  if (dimRed == "t-SNE") {
    pred_filtered_plot <- plot2D(pred_filtered,
                                 fcol = "svm.pred",
                                 method = "t-SNE",
                                 methargs = list(perplexity = perplexity)
    )
  } else if (dimRed == "PCA") {
    pred_filtered_plot <- plot2D(pred_filtered,
                                 fcol = "svm.pred",
                                 dims = dims,
                                 method = "PCA")
  } else {
    stop('Error: no dimensionality reduction method provided. Choose "PCA" or "t-SNE" for the dimRed argument.\n')
  }

  pred_filtered_plot_T <- tibble(ID = rownames(pred_filtered_plot),
                                 Dim1 = pred_filtered_plot[, 1],
                                 Dim2 = pred_filtered_plot[, 2])

  ## add prediction columns to coordinates ##
  pred_filtered_plot_T <- left_join(pred_filtered_plot_T, pred_filtered_T, by = "ID")

  ## add protein descriptions ##
  pred_filtered_plot_T <- left_join(pred_filtered_plot_T, descLopit_data, by = "ID")

  ## create interactive label ##
  pred_filtered_plot_T <- pred_filtered_plot_T %>%
    group_by(ID) %>%
    mutate(label = paste0(ID, ", marker: ", markers, ", svm: ", svm, ", svm.score: ", round(svm.scores, 2), " - ", desc2, collapse = ""))

  ## convert svm.pred to factor ##
  pred_filtered_plot_T$svm.pred <- factor(pred_filtered_plot_T$svm.pred, levels = levels)

  ## filter to marker proteins only for overlay ##
  markers_plot_T <- pred_filtered_plot_T %>%
    filter(markers != "unknown")

  ## non-interactive plot ##
  tSNE_plot <- ggplot(pred_filtered_plot_T, aes(label = label)) +
    geom_jitter(data = pred_filtered_plot_T, aes(Dim1, Dim2, color = svm.pred), size = 2, alpha = 0.5) +
    geom_point(data = markers_plot_T, aes(Dim1, Dim2), color = "black", size = 3) +
    geom_point(data = markers_plot_T, aes(Dim1, Dim2, color = svm.pred), size = 2) +
    theme_minimal(base_size = 20) +
    labs(color = "Localization") +
    scale_color_manual(values = colors)

  png(paste0(outDir, "/", dimRed, outSuffix, ".png"),
      width = 4000,
      height = 2500,
      res = 300
  )

  plot(tSNE_plot)

  dev.off()

  ## NA count plot
  make_NA_plot(pred_filtered_plot_T = pred_filtered_plot_T,
               dimRed = dimRed,
               outSuffix = outSuffix,
               outDir = outDir
  )

  ## interactive plotly HTML
  plotly <- ggplotly(tSNE_plot)
  htmlwidgets::saveWidget(as_widget(plotly),
                          paste0(outDir, "/", dimRed, "_plot_interactive", outSuffix, ".html")
  )

  write_tsv(pred_filtered_plot_T,
            paste0(outDir,
                   "/",
                   dimRed,
                   "_plot_data",
                   outSuffix,
                   ".tsv")
  )

  ########################
  ## return the dataset ##
  ########################
  returnList <- list(data = pred_filtered_plot_T)

  return(returnList)

  cat("Everything is done.\n")
}

# TAGM-MAP prediction sub-function
run_TAGMMAP <- function(data,
                        numIter_tagmmap,
                        seed,
                        p_cutoff_tagmmap,
                        dimRed,
                        perplexity,
                        dims,
                        levels,
                        colors,
                        outDir,
                        outSuffix,
                        descLopit_data
) {

  ## train TAGM-MAP model ##
  params <- tagmMapTrain(data,
                         numIter = numIter_tagmmap,
                         fcol = "markers")

  ## convergence plot ##
  png(paste0(outDir, "/", "Convergence_graph", outSuffix, ".png"),
      width = 1500,
      height = 1500,
      res = 300)

  params %>%
    logPosteriors %>%
    plot(type = "b", col = "blue", cex = 0.3,
         ylab = "log-posterior",
         xlab = "Iteration")

  dev.off()

  ## predict localisations using TAGM-MAP ##
  tagm_map_res <- tagmMapPredict(data,
                                 fcol = "markers",
                                 params = params)

  # save results
  save(tagm_map_res, file = paste0(outDir, "/TAGM-MAPres", outSuffix, ".Rdata"))

  #######################
  #### Visualisation  ####
  #######################
  pred_T <- fData(tagm_map_res)

  plot_prob_dist <- ggplot(pred_T, aes(tagm.map.allocation, tagm.map.probability)) +
    geom_boxplot() +
    geom_hline(yintercept = p_cutoff_tagmmap, color = "red", linetype = 2) +
    theme_minimal(base_size = 10) +
    coord_flip() +
    xlab("prediction") +
    ylab("prediction probability")

  png(paste0(outDir, "/", "Dist_prob_prediction_with_cut_off", outSuffix, ".png"),
      width = 1500,
      height = 1500,
      res = 300)

  plot(plot_prob_dist)

  dev.off()

  cat("Creating graphical outputs.\n")

  ## composite probability: localisation probability * (1 - outlier probability)
  tagm_prob <- fData(tagm_map_res)[, "tagm.map.probability"]
  tagm_out <- 1 - fData(tagm_map_res)[, "tagm.map.outlier"]

  fData(tagm_map_res)[, "overall_prob"] <- tagm_prob * tagm_out

  pred_filtered <- getPredictions(tagm_map_res,
                                  fcol = "tagm.map.allocation",
                                  scol = "overall_prob",
                                  t = p_cutoff_tagmmap)

  pred_filtered_T <- fData(pred_filtered)
  pred_filtered_T$ID <- rownames(pred_filtered_T)

  ###################################################################
  #### Dimensionality reduction and plot                          ####
  ###################################################################
  set.seed(seed)

  if (dimRed == "t-SNE") {
    pred_filtered_plot <- plot2D(pred_filtered,
                                 fcol = "tagm.map.allocation.pred",
                                 method = "t-SNE",
                                 methargs = list(perplexity = perplexity)
    )
  } else if (dimRed == "PCA") {
    pred_filtered_plot <- plot2D(pred_filtered,
                                 fcol = "tagm.map.allocation.pred",
                                 dims = dims,
                                 method = "PCA")
  } else {
    stop('Error: no dimensionality reduction method provided. Choose "PCA" or "t-SNE" for the dimRed argument.\n')
  }

  pred_filtered_plot_T <- tibble(ID = rownames(pred_filtered_plot),
                                 Dim1 = pred_filtered_plot[, 1],
                                 Dim2 = pred_filtered_plot[, 2])

  pred_filtered_plot_T <- left_join(pred_filtered_plot_T, pred_filtered_T, by = "ID")
  pred_filtered_plot_T <- left_join(pred_filtered_plot_T, descLopit_data, by = "ID")

  pred_filtered_plot_T <- pred_filtered_plot_T %>%
    group_by(ID) %>%
    mutate(label = paste0(ID, ", marker: ", markers, ", allocation: ", tagm.map.allocation, ", prob: ", round(overall_prob, 4), " - ", desc2, collapse = ""))

  pred_filtered_plot_T$tagm.map.allocation.pred <- factor(pred_filtered_plot_T$tagm.map.allocation.pred, levels = levels)

  markers_plot_T <- pred_filtered_plot_T %>%
    filter(markers != "unknown")

  ## non-interactive plot ##
  tSNE_plot <- ggplot(pred_filtered_plot_T, aes(label = label)) +
    geom_jitter(data = pred_filtered_plot_T, aes(Dim1, Dim2, color = tagm.map.allocation.pred), size = 2, alpha = 0.5) +
    geom_point(data = markers_plot_T, aes(Dim1, Dim2), color = "black", size = 3) +
    geom_point(data = markers_plot_T, aes(Dim1, Dim2, color = tagm.map.allocation.pred), size = 2) +
    theme_minimal(base_size = 20) +
    labs(color = "Localization") +
    scale_color_manual(values = colors)

  png(paste0(outDir, "/", dimRed, outSuffix, ".png"),
      width = 4000,
      height = 2500,
      res = 300
  )

  plot(tSNE_plot)

  dev.off()

  ## NA count plot
  make_NA_plot(pred_filtered_plot_T = pred_filtered_plot_T,
               dimRed = dimRed,
               outSuffix = outSuffix,
               outDir = outDir
  )

  ## interactive plotly HTML
  plotly <- ggplotly(tSNE_plot)
  htmlwidgets::saveWidget(as_widget(plotly),
                          paste0(outDir, "/", dimRed, "_plot_interactive", outSuffix, ".html")
  )

  write_tsv(pred_filtered_plot_T,
            paste0(outDir,
                   "/",
                   dimRed,
                   "_plot_data",
                   outSuffix,
                   ".tsv")
  )

  ########################
  ## return the dataset ##
  ########################
  returnList <- list(data = pred_filtered_plot_T)

  return(returnList)

  cat("Everything is done.\n")
}

# MCMC prediction sub-function
# Currently no manual prior specification (see https://f1000research.com/articles/8-446/v1
# if prior customisation is needed in future via the S0 argument)
run_MCMC <- function(data,
                     numIter,
                     burnin,
                     thin,
                     numChains,
                     MCMCres_path,
                     MCMCres_keepChains,
                     extra_thin,
                     extra_burnin,
                     geweke_pval,
                     seed,
                     dimRed,
                     perplexity,
                     dims,
                     levels,
                     colors,
                     p_cutoff,
                     outDir,
                     outSuffix,
                     descLopit_data
) {

  #######################
  ## Run the MCMC      ##
  #######################
  if (is.null(MCMCres_path)) {
    cat("Running MCMC from scratch (no existing results provided).\n")
    n_workers <- as.integer(Sys.getenv("PBS_NCPUS", unset = "4"))

    multicoreParam <- MulticoreParam(workers = n_workers)
    register(multicoreParam)

    max_tries <- 3
    attempt <- 1
    success <- FALSE
    mcmcRes <- NULL

    while (attempt <= max_tries && !success) {
      tryCatch({
        message(sprintf("Attempt %d of %d", attempt, max_tries))
        mcmcRes <- tagmMcmcTrain(data,
                                 numIter = numIter,
                                 burnin = burnin,
                                 thin = thin,
                                 numChains = numChains,
                                 BPPARAM = multicoreParam)

        success <- TRUE
      }, error = function(e) {
        message("Error in MCMC run: ", e$message)
        attempt <<- attempt + 1
        if (attempt > max_tries) {
          stop("Maximum retry attempts reached. MCMC failed.")
        } else {
          message("Retrying...")
        }
      })
    }

    if (success) {
      save(mcmcRes, file = paste0(outDir, "/MCMCres", outSuffix, ".Rdata"))
    }

  } else {
    cat("Loading existing MCMC results.\n")
    load(MCMCres_path)
  }

  if (!is.null(extra_thin) && !is.null(extra_burnin)) {
    mcmcRes <- mcmc_thin_chains(mcmcRes, freq = extra_thin)
    mcmcRes <- mcmc_burn_chains(mcmcRes, n = extra_burnin)
    cat("Applied additional thinning and additional burn-in removal.\n")
  } else if (!is.null(extra_thin)) {
    mcmcRes <- mcmc_thin_chains(mcmcRes, freq = extra_thin)
    cat("Applied additional thinning.\n")
  } else if (!is.null(extra_burnin)) {
    mcmcRes <- mcmc_burn_chains(mcmcRes, n = extra_burnin)
    cat("Applied additional burn-in removal from the start of chains.\n")
  }

  ###########################
  ## Convergence checking  ##
  ###########################
  nChains <- length(mcmcRes)
  selectedChains <- 1:nChains

  ######################################
  ## Diagnostic plots and statistics  ##
  ######################################
  out <- mcmc_get_outliers(mcmcRes)
  meanoutProb <- mcmc_get_meanoutliersProb(mcmcRes)
  meanAlloc <- mcmc_get_meanComponent(mcmcRes)

  statistics_L <- list(out, meanoutProb, meanAlloc)
  names(statistics_L) <- c("Outliers_number", "Outliers_probability", "Mean_component_allocation")

  map2(statistics_L,
       names(statistics_L),
       ~{
         make_MCMC_stat_plots(matrix = .x,
                              output = paste0(outDir,
                                              "/",
                                              .y,
                                              outSuffix
                              ),
                              nChains = nChains
         )
       })

  ############################
  ## Gelman diagnostics     ##
  ############################
  if (nChains == 1) {
    cat("Skipping Gelman diagnostics. More than one chain is required!\n")

  } else if (nChains == 2) {
    cat("Running Gelman diagnostics for two chains.\n")

    gelman_diag_res_T <- run_gelman(list(out, meanAlloc, meanoutProb),
                                    c("out", "meanAlloc", "meanoutProb")
    )

    write_tsv(gelman_diag_res_T,
              paste0(outDir,
                     "/gelman_diagnostics",
                     outSuffix,
                     ".tsv")
    )

    if (max(gelman_diag_res_T$Upper_CI) <= 1.2) {
      cat("Both chains converged according to Gelman diagnostics.\n")
      selectedChains_gelman <- selectedChains
    } else {
      cat("The combination of the two chains failed Gelman diagnostics.\n")
      selectedChains_gelman <- vector()
    }

  } else {
    cat("Running Gelman diagnostics for all chains.\n")

    gelman_diag_res_T <- run_gelman(list(out, meanAlloc, meanoutProb),
                                    c("out", "meanAlloc", "meanoutProb")
    )

    write_tsv(gelman_diag_res_T,
              paste0(outDir,
                     "/gelman_diagnostics",
                     outSuffix,
                     ".tsv")
    )

    upper_CI <- max(gelman_diag_res_T[, 2])

    ########################
    ## Select good chains ##
    ########################
    if (upper_CI <= 1.2) {
      cat("All chains converged; keeping all.\n")
      selectedChains_gelman <- selectedChains
    } else {
      cat("Gelman diagnostics indicate not all chains converged.\nRunning all chain combinations and selecting the largest converged subset.\n")

      combinations_L <- map(2:(nChains - 1),
                            ~{
                              combinations_map <- combn(1:nChains, .x)

                              combinations_map_L <- map(array_branch(combinations_map, 2),
                                                        ~{
                                                          .x
                                                        })

                              return(combinations_map_L)
                            }) %>%
        flatten()

      combinations_gelmanDia_res_T <- map(combinations_L,
                                          ~{
                                            gelman_diag_map_res_T <- run_gelman(list(out, meanAlloc, meanoutProb),
                                                                                c("out", "meanAlloc", "meanoutProb"),
                                                                                subset = .x)

                                            gelman_diag_map_res_T$nChain <- length(.x)
                                            gelman_diag_map_res_T$ChainCombList <- list(.x)
                                            gelman_diag_map_res_T$ChainCombChar <- paste0(.x, collapse = ",")
                                            return(gelman_diag_map_res_T)
                                          }) %>%
        purrr::reduce(rbind)

      # select the largest chain combination that passes Gelman diagnostics
      combinations_gelmanDia_res_T <- combinations_gelmanDia_res_T %>%
        group_by(ChainCombChar) %>%
        mutate(max_upper_CI = max(Upper_CI)) %>%
        filter(max_upper_CI <= 1.2) %>%
        mutate(pass = n() == 3) %>%
        filter(pass) %>%
        arrange(desc(nChain), desc(Upper_CI))

      selectedChains_gelman <- combinations_gelmanDia_res_T$ChainCombList[1][[1]]

      combinations_gelmanDia_res_T <- combinations_gelmanDia_res_T %>%
        dplyr::select(-ChainCombList)

      write_tsv(combinations_gelmanDia_res_T,
                paste0(outDir,
                       "/gelman_diagnostics_comb",
                       outSuffix,
                       ".tsv")
      )

      cat(paste0("Chains selected based on Gelman diagnostics: "),
          paste0(selectedChains_gelman, collapse = ","),
          "\n")
    }
  }

  #######################################################################################################
  ## Geweke test: compare first 10% of chain with last 50% to assess within-chain convergence         ##
  ## (can be affected by burn-in if the chain converges later than the discarded period)              ##
  #######################################################################################################
  FailedChains <- map(statistics_L,
                      ~{
                        geweke_test_res <- geweke_test(.x)

                        geweke_test_res_T <- as_tibble(geweke_test_res)
                        write_tsv(geweke_test_res_T,
                                  paste0(outDir,
                                         "/geweke_test_res",
                                         outSuffix,
                                         ".tsv")
                        )

                        failTest <- geweke_test_res[2, ] <= geweke_pval
                        failChains <- names(failTest[failTest])

                        if (!is.null(failChains)) {
                          failChainsNumber <- parse_number(failChains)
                          return(failChainsNumber)
                        } else {
                          failChainsNumber <- vector()
                          return(failChainsNumber)
                        }

                      }) %>%
    unlist() %>%
    unique()

  ################################################################
  ## Remove chains that failed the Geweke test                  ##
  ################################################################
  if (is_empty(FailedChains)) {
    cat("No chain was removed during the Geweke test.\n")
    selectedChains_geweke <- selectedChains
  } else {
    removedChain <- selectedChains[selectedChains == FailedChains]
    cat("Based on the Geweke test, Chain ",
        paste0(removedChain, collapse = ", "),
        " was removed.\n",
        sep = "")
    selectedChains_geweke <- selectedChains[selectedChains != FailedChains]
  }

  ##############################################################################################################
  ## Apply chain selection strategy                                                                          ##
  ##############################################################################################################
  if (typeof(MCMCres_keepChains) == "double") {
    cat("Keeping the manually specified chains.\n")
    selectedChains <- MCMCres_keepChains
  } else if (MCMCres_keepChains == "geweke_and_gelman") {
    # intersection: chains that pass both tests
    selectedChains <- selectedChains_gelman[selectedChains_gelman %in% selectedChains_geweke]
    cat("Keeping chains that pass both Gelman and Geweke criteria.\n")
  } else if (MCMCres_keepChains == "geweke") {
    selectedChains <- selectedChains_geweke
    cat("Keeping chains based on the Geweke test only.\n")
  } else if (MCMCres_keepChains == "gelman") {
    selectedChains <- selectedChains_gelman
    cat("Keeping chains based on the Gelman test only.\n")
  } else {
    errorCondition("Wrong value for MCMCres_keepChains. Options are: \"geweke\", \"gelman\", \"geweke_and_gelman\", or a numeric vector of chain indices.")
  }

  ################################
  ## Keep the converged chains  ##
  ################################
  if (is_empty(selectedChains)) {
    errorCondition("No chains remain after selection. Check the supplied chain indices or statistical tests.\n")
  }

  mcmcRes_converged <- mcmcRes[selectedChains]

  ########################
  ## Pool the chains    ##
  ########################
  # Pooling unconverged chains with converged ones leads to poor results; only converged chains are pooled
  mcmcRes_converged_pooled <- mcmc_pool_chains(mcmcRes_converged)

  # BUG FIX: original code passed the object as second positional argument instead of using file=
  save(mcmcRes_converged_pooled,
       file = paste0(outDir,
                     "/mcmcRes_converged_pooled_for_violin_plots",
                     outSuffix,
                     ".Rdata")
  )

  #####################################
  ## Add predictions to the data     ##
  #####################################
  mcmcRes_converged_pooled <- tagmMcmcProcess(mcmcRes_converged_pooled)

  mcmcRes_converged_pooled_MSdata <- tagmPredict(object = data,
                                                 params = mcmcRes_converged_pooled,
                                                 probJoint = TRUE)

  ## composite probability: localisation probability * (1 - outlier probability) * (1 - Shannon entropy)
  tagm_prob <- fData(mcmcRes_converged_pooled_MSdata)[, "tagm.mcmc.probability"]
  tagm_out <- 1 - fData(mcmcRes_converged_pooled_MSdata)[, "tagm.mcmc.outlier"]
  shan_out <- 1 - fData(mcmcRes_converged_pooled_MSdata)[, "tagm.mcmc.mean.shannon"]

  fData(mcmcRes_converged_pooled_MSdata)[, "overall_prob"] <- tagm_prob * tagm_out * shan_out

  ## apply probability threshold
  mcmcRes_converged_pooled_MSdata <- getPredictions(mcmcRes_converged_pooled_MSdata,
                                                    fcol = "tagm.mcmc.allocation",
                                                    scol = "overall_prob",
                                                    t = p_cutoff)

  ###################
  ## Visualisation ##
  ###################
  tibble_data_T <- fData(mcmcRes_converged_pooled_MSdata)
  colnames(tibble_data_T)[1] <- "ID"

  set.seed(seed)

  if (dimRed == "t-SNE") {
    pred_filtered_plot <- plot2D(mcmcRes_converged_pooled_MSdata,
                                 fcol = "tagm.mcmc.probability",
                                 cex = fData(mcmcRes_converged_pooled_MSdata)$tagm.mcmc.probability,
                                 method = "t-SNE",
                                 methargs = list(perplexity = perplexity),
                                 main = "TAGM MCMC allocations"
    )
  } else if (dimRed == "PCA") {
    pred_filtered_plot <- plot2D(mcmcRes_converged_pooled_MSdata,
                                 fcol = "tagm.mcmc.probability",
                                 cex = fData(mcmcRes_converged_pooled_MSdata)$tagm.mcmc.probability,
                                 method = "PCA",
                                 dims = dims,
                                 main = "TAGM MCMC allocations"
    )
  } else {
    stop('Error: no dimensionality reduction method provided. Choose "PCA" or "t-SNE" for the dimRed argument.\n')
  }

  pred_filtered_plot_T <- tibble(ID = rownames(pred_filtered_plot),
                                 Dim1 = pred_filtered_plot[, 1],
                                 Dim2 = pred_filtered_plot[, 2])

  pred_filtered_plot_T <- left_join(pred_filtered_plot_T, tibble_data_T, by = "ID")
  pred_filtered_plot_T <- left_join(pred_filtered_plot_T, descLopit_data, by = "ID")

  pred_filtered_plot_T <- pred_filtered_plot_T %>%
    group_by(ID) %>%
    mutate(label = paste0(ID,
                          ", marker: ", markers,
                          ", localization: ", tagm.mcmc.allocation,
                          ", probability: ", round(overall_prob, 2),
                          " - ",
                          desc2,
                          collapse = "")
    )

  # BUG FIX: original code factored overall_prob (a continuous variable) using compartment levels,
  # which would produce NA factors and break the colour scale.
  # The correct column to factor for colouring is tagm.mcmc.allocation.pred.
  pred_filtered_plot_T$tagm.mcmc.allocation.pred <- factor(
    fData(mcmcRes_converged_pooled_MSdata)[rownames(fData(mcmcRes_converged_pooled_MSdata)) %in% pred_filtered_plot_T$ID, "tagm.mcmc.allocation.pred"],
    levels = levels
  )

  markers_plot_T <- pred_filtered_plot_T %>%
    filter(markers != "unknown")

  ## non-interactive plot ##
  tSNE_plot <- ggplot(pred_filtered_plot_T, aes(label = label)) +
    geom_jitter(data = pred_filtered_plot_T, aes(Dim1, Dim2, color = tagm.mcmc.allocation.pred), size = 2, alpha = 0.5) +
    geom_point(data = markers_plot_T, aes(Dim1, Dim2), color = "black", size = 3) +
    geom_point(data = markers_plot_T, aes(Dim1, Dim2, color = tagm.mcmc.allocation.pred), size = 2) +
    theme_minimal(base_size = 20) +
    labs(color = "Localization") +
    scale_color_manual(values = colors)

  png(paste0(outDir, "/", dimRed, outSuffix, ".png"),
      width = 4000,
      height = 2500,
      res = 300
  )

  plot(tSNE_plot)

  dev.off()

  ## NA count plot
  make_NA_plot(pred_filtered_plot_T = pred_filtered_plot_T,
               dimRed = dimRed,
               outSuffix = outSuffix,
               outDir = outDir
  )

  ## interactive plotly HTML
  plotly <- ggplotly(tSNE_plot)
  htmlwidgets::saveWidget(as_widget(plotly),
                          paste0(outDir, "/", dimRed, "_plot_interactive", outSuffix, ".html")
  )

  #######################
  ## Shannon entropy plot ##
  #######################
  tSNE_shannon_plot <- ggplot(pred_filtered_plot_T, aes(label = label)) +
    geom_jitter(data = pred_filtered_plot_T, aes(Dim1, Dim2,
                                                 color = tagm.mcmc.allocation,
                                                 size = tagm.mcmc.mean.shannon),
                alpha = 0.5) +
    theme_minimal(base_size = 20) +
    labs(color = "Localization", size = "Shannon entropy") +
    scale_color_manual(values = colors)

  png(paste0(outDir, "/", dimRed, "_shannon", outSuffix, ".png"),
      width = 4000,
      height = 2500,
      res = 300
  )

  plot(tSNE_shannon_plot)

  dev.off()

  ## interactive Shannon entropy plot
  plotly_shannon <- ggplotly(tSNE_shannon_plot)
  htmlwidgets::saveWidget(as_widget(plotly_shannon),
                          paste0(outDir, "/", dimRed, "_shannon_plot_interactive", outSuffix, ".html")
  )

  ###########################
  ## Write the final results ##
  ###########################
  # expand per-compartment probability columns
  all_loc_prob <- as_tibble(pred_filtered_plot_T$tagm.mcmc.joint)
  colnames(all_loc_prob) <- paste0("prob_", colnames(all_loc_prob))

  pred_filtered_plot_T <- pred_filtered_plot_T %>%
    dplyr::select(-tagm.mcmc.joint) %>%
    cbind(all_loc_prob)

  write_tsv(pred_filtered_plot_T,
            paste0(outDir,
                   "/",
                   dimRed,
                   "_plot_data",
                   outSuffix,
                   ".tsv")
  )

  ########################
  ## return the dataset ##
  ########################
  returnList <- list(data = pred_filtered_plot_T)

  return(returnList)

  cat("Everything is done.\n")

}

run_MCMC_bandle <- function(data,
                            numIter,
                            burnin,
                            thin,
                            numChains,
                            MCMCres_path,
                            MCMCres_keepChains,
                            extra_thin,
                            extra_burnin,
                            geweke_pval,
                            seed,
                            dimRed,
                            perplexity,
                            dims,
                            p_cutoff,
                            outDir,
                            outSuffix,
                            descLopit_data,
                            hyppar,    # BANDLE: priors for hyper parameters
                            dirPrior,  # BANDLE: Dirichlet prior matrix;
                                       # reflects expected probability of differential localisation between datasets
                            conditions
) {

  ############################################
  ## Prepare hyper-parameter priors         ##
  ############################################
  mrkCl <- getMarkerClasses(data[[1]], fcol = "markers")

  set.seed(seed)

  K <- length(mrkCl)
  pc_prior <- matrix(NA, ncol = 3, K)
  pc_prior[seq.int(1:K), ] <- matrix(rep(hyppar, each = K), ncol = 3)

  ## Fit GP priors to each dataset separately
  set.seed(seed)
  gpParams <- lapply(data,
                     function(x) fitGPmaternPC(x, hyppar = pc_prior))

  ## Plot GP prior fits overlaid on marker data for each dataset
  map(1:length(data),
      ~{
        png(paste0(outDir,
                   "/",
                   "BANDLE_fitted_hyperpara_markers_fit_dataset_",
                   .x,
                   outSuffix,
                   ".png"
        )
        )

        par(mfrow = c(4, 3))
        plotGPmatern(data[[.x]], gpParams[[.x]])
        dev.off()
      })

  #################################
  ## Prepare the Dirichlet prior ##
  #################################
  dirPrior <- diag(rep(1, K)) + matrix(dirPrior, nrow = K, ncol = K)

  map(1:length(data),
      ~{
        png(paste0(outDir,
                   "/",
                   "BANDLE_fitted_dirichlet_prior_determined_prob_other_localization_",
                   .x,
                   outSuffix,
                   ".png"
        )
        )

        ## Determine the prior probability of more than 15 differential localisation events
        predDirPrior <- prior_pred_dir(object = data[[.x]],
                                       dirPrior = dirPrior,
                                       q = 15)

        ## Histogram of differential localisation probability based on priors
        hist(predDirPrior$priornotAlloc, col = getStockcol()[1])

        dev.off()
      })

  ## Create lists of replicate MSnSets per condition
  conditions_L <- map(unique(conditions),
                      ~{
                        datasetGroups <- which(conditions == .x)

                        map(datasetGroups,
                            ~{
                              data[[.x]]
                            })
                      })

  ###################################
  ## Run the BANDLE MCMC algorithm ##
  ###################################
  if (is.null(MCMCres_path)) {
    cat("Running BANDLE MCMC from scratch (no existing results provided).\n")
    n_workers <- as.integer(Sys.getenv("PBS_NCPUS", unset = "4"))

    multicoreParam <- MulticoreParam(workers = n_workers)
    register(multicoreParam)

    max_tries <- 3
    attempt <- 1
    success <- FALSE
    mcmcRes <- NULL

    while (attempt <= max_tries && !success) {
      tryCatch({
        message(sprintf("Attempt %d of %d", attempt, max_tries))

        badleres <- bandle(objectCond1 = conditions_L[[1]],
                           objectCond2 = conditions_L[[2]],
                           numIter = numIter,
                           burnin = burnin,
                           thin = thin,
                           gpParams = gpParams,
                           numChains = numChains,
                           dirPrior = dirPrior,
                           pcPrior = pc_prior,
                           seed = seed,
                           BPPARAM = multicoreParam
        )

        success <- TRUE
      }, error = function(e) {
        message("Error in BANDLE MCMC run: ", e$message)
        attempt <<- attempt + 1
        if (attempt > max_tries) {
          stop("Maximum retry attempts reached. BANDLE MCMC failed.")
        } else {
          message("Retrying...")
        }
      })
    }

    if (success) {
      save(badleres, file = paste0(outDir, "/BANDLEres", outSuffix, ".Rdata"))
    }

  } else {
    cat("Loading existing BANDLE results.\n")
    load(MCMCres_path)
  }

  ######################
  ## Check the chains ##
  ######################
  Gelman_res <- bandleres %>%
    calculateGelman()

  Gelman_res_T <- map2(Gelman_res,
                       names(Gelman_res),
                       ~{
                         df_map_T <- as_tibble(.x)
                         df_map_T$gelmanDiag <- rownames(.x)
                         df_map_T$condition <- .y

                         return(df_map_T)

                       }) %>%
    purrr::reduce(rbind)

  write_tsv(Gelman_res_T, paste0(outDir, "/BANDLEres_gelman_diagnostics", outSuffix, ".tsv"))


  #############################
  ## Outlier diagnostic plot ##
  #############################
  png(paste0(outDir,
             "/",
             "BANDLEres_chain_outliers_",
             outSuffix,
             ".png"
  ), width = 7500,
  height = 7500,
  res = 300
  )

  bandleres %>%
    plotOutliers()

  dev.off()

  ################################
  ## Optionally subset chains   ##
  ################################
  if (!is_empty(MCMCres_keepChains)) {
    cat("Keeping the specified chains.\n")
    bandleres_converged <- bandleres[MCMCres_keepChains]
  } else {
    cat("Keeping all chains. Provide a numeric vector to MCMCres_keepChains to select specific chains.\n")
    bandleres_converged <- bandleres
  }

  ## Process BANDLE results
  params <- bandleProcess(bandleres_converged)

  ## Predict localisations for both conditions
  res <- bandlePredict(objectCond1 = conditions_L[[1]],
                       objectCond2 = conditions_L[[2]],
                       params = params,
                       fcol = "markers")

  res_cond1 <- res[[1]][[1]]
  res_cond2 <- res[[2]][[1]]

  res_per_cond <- c(res_cond1, res_cond2)

  ## apply probability threshold
  res_per_cond <- map(res_per_cond,
                      ~{
                        map_data <- .x
                        class_prob <- fData(map_data)$bandle.probability
                        out_prob <- 1 - fData(map_data)$bandle.outlier

                        fData(map_data)$bandle.probability.overall <- class_prob * out_prob

                        map_data <- getPredictions(map_data,
                                                   fcol = "bandle.allocation",
                                                   scol = "bandle.probability.overall",
                                                   mcol = "markers",
                                                   t = p_cutoff)

                        return(map_data)
                      })

  ## Visualise BANDLE results per condition
  map2(res_per_cond,
       c("cond1", "cond2"),
       ~{
         map_data <- .x
         # BUG NOTE: original code referenced res_unstim_rep1 which is undefined in this scope;
         # replaced with map_data
         res_no_mrk <- unknownMSnSet(map_data, fcol = "markers")

         alloc <- map_data %>%
           fData() %>%
           pull(bandle.allocation.pred)

         pe <- res_no_mrk %>%
           fData() %>%
           pull(bandle.probability)

         ## bar plot of predicted localisations
         png(paste0(outDir,
                    "/",
                    "BANDLEres_predictions_with_prob_cutoff_",
                    .y,
                    outSuffix,
                    ".png"
         )
         )

         barplot(alloc %>% table,
                 las = 2, main = "Predicted location",
                 ylab = "Number of proteins")

         dev.off()

         ## box plot of posterior distributions by localisation
         png(paste0(outDir,
                    "/",
                    "BANDLEres_predictions_posterior_prob_",
                    .y,
                    outSuffix,
                    ".png"
         )
         )

         boxplot(pe ~ alloc, las = 2, main = "Posterior",
                 ylab = "Probability")

         dev.off()

       })

}
