library(tidyverse) #data prsing
library(plotly) #data visualization
library(scales) #data scaling
library(coda) #for gelman.diag
library(pRoloc) #lopit anna
library(xlsx)
library(pheatmap) #for heat map 
library(dbscan)
library(pRolocGUI)
library(bandle)
library(ggforce) #library for ploting  circle in ggplot

#######################################################
## pokud se bude opakovat problem with chain failing ##
#######################################################
#Temporarily replace MulticoreParam with SerialParam() to see if the problem disappears (problem with chains failing).
#If yes, this confirms it's parallelization/environment related.

###########################
## functions for calling ##
###########################
#function that makes consensus from more lopits
lopit_consensus <- function(data,
                            smv_pred_pos,
                            consensusCount,
                            compartments) {
  #data = data with the svm.predictions
  #smv_pred_pos = position of the svm.pred from the lopit function
  #consensusCount = how many times the protein must be seen in the given compartment to have consensus
  #compartments = a vector of all compartments names
  countedCompartments_T <- map(compartments,
                               ~{
                                 #select predictions
                                 allPredTien2_noNA <- data %>% 
                                   dplyr::select(all_of(smv_pred_pos))
                                 
                                 #count the predictions per lopit per line
                                 countedComp <- apply(allPredTien2_noNA, 1 , str_count) %>% 
                                   t() %>%
                                   apply(., 1, sum, na.rm = T) %>% 
                                   unlist()
                                 
                                 #make atibble out of it
                                 countedComp_T <- tibble(countedComp)
                                 colnames(countedComp_T) <- .x
                                 
                                 return(countedComp_T)
                                 
                               }) %>% 
    purrr::reduce(cbind)
  
  #create the consensus comparment
  consensusCompartment_Sch <- map(array_branch(countedCompartments_T, 1),
                                  ~{
                                    #get which lopits have seen the portein in the same compartment more than twice
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

#function that takes a vector of IDs and always returns the same ID as key and other IDs, for the group (base on my matching of IDs for TV)
lopit_get_uniqueID <- function(tibble, 
                               IDcol,
                               key_tibble, 
                               sep = ";") {
  #tibble = with ID column were spectrometer thinks the IDs to be a single group, excpects vector of separated IDs
  #IDcol = position of the ID column 
  #key_tibble = tibble with IDs matches among TVAG, TVAGG3 and Uniprot (my file in general datasets)
  #from = mapping either from TVAG, TVAGG3 or Entry (which is uniprot), e.g. "TVAG"
  #sep = if a single string with IDs with separator is provided provide the separate so it can be split
  
  #split the IDs column
  tibble <- tibble %>% 
    separate_rows(all_of(IDcol), sep = sep)
  
  #get the IDs from your dataset, possibly mixed
  IDs_tibble_T <- tibble(ID = str_remove(unname(unlist(tibble[,IDcol])), "-1-p1$")) %>% #just remove often found in ID from proteomic analysis
    mutate(rowN = row_number())
  
  #where are these IDs found in all possible IDs from my key file
  key_tibble_sub <- key_tibble %>% 
    filter(Entry %in% IDs_tibble_T$ID
           | TVAG %in% IDs_tibble_T$ID
           | TVAGG3 %in% IDs_tibble_T$ID)
  
  #join by every column and then paste and unique an repeat
  keyCols <- colnames(key_tibble_sub)
  
  for (i in 1:length(keyCols)) {
    
    colnames(IDs_tibble_T)[1] <- keyCols[i] #rename the colname for merging
    
    #now join
    IDs_tibble_T <- left_join(IDs_tibble_T, key_tibble_sub, by = keyCols[i])
    
    #remove replicates and repeat
    ExtraID_col_pos <- c(ncol(IDs_tibble_T) -1, ncol(IDs_tibble_T))
    colnames(IDs_tibble_T)[ExtraID_col_pos] <- c("ID1", "ID2")
    
    IDs_tibble_T <- IDs_tibble_T %>% 
      group_by(rowN) %>% 
      mutate(OtherIDs = paste0(unique(c(ID1, ID2)), collapse = ";")) %>% 
      dplyr::select(1,2, matches("OtherIDs"))
    
    colnames(IDs_tibble_T)[(3+(i-1))] <- paste0(colnames(IDs_tibble_T)[(3+(i-1))], i)
  }
  
  #chnage the name of the IDs column used for joining
  colnames(IDs_tibble_T)[1] <- "oriID"
  
  #correct the other IDs
  IDs_tibble_T <- IDs_tibble_T %>% 
    group_by(rowN) %>% 
    mutate(OtherIDs = paste0(unique(c(unlist(str_split(OtherIDs1, ";")),
                                      unlist(str_split(OtherIDs2, ";")),
                                      unlist(str_split(OtherIDs3, ";")))
    ),
    collapse = ";"), #collpase all Other IDs
    OtherIDs = str_remove_all(OtherIDs, "(NA;|;NA)")#remove the NAs
    ) %>% 
    dplyr::select(1, 2 , OtherIDs) %>% 
    unique()
  
  #create TVAGG3 key
  IDs_tibble_T <- IDs_tibble_T %>% 
    group_by(rowN) %>% 
    mutate(key = paste0(unlist(str_extract_all(OtherIDs, "TVAGG3_\\d{7}")), collapse = ";")
    ) %>% 
    separate_rows(key, sep = ";") %>% 
    arrange(key) %>% #arrange to have the first ID (lewest number always first)
    mutate(key = if_else(key == "",
                         oriID,
                         key)
    ) %>% 
    mutate(keyPos = row_number()) %>% 
    filter(keyPos == 1)
  
  
  #join with the original data
  colnames(IDs_tibble_T)[1] <- colnames(tibble)[IDcol]
  tibble_res <- left_join(IDs_tibble_T, tibble)
  
  #create joined IDs
  return(tibble_res)
  
}


################
#NOT USED USING THE BUNDLE JOINING
###############
# #funkce co merguje data pro to abychom jeli replikaty dohromady
# #this merges data into a tsv
# lopit_merge_replicates <- function(replicatesPath,
#                                    IDcolumn, 
#                                    IDsep = ";",
#                                    outputPath) {
#   #replicatesPath = path were the replicates are (full paths)
#   #IDcolumn = position of the ID columbn for each dataset
#   #IDsep = pattern that separates indistinguishable IDs
#   #outputPath = ful output path for the joined dataset
#   
#   #load and join the data
#   allJoined <- pmap(list(replicatesPath,
#                          paste0("rep", 1:length(replicatesPath)),
#                          IDcolumn),
#                     ~{
#                       #load the data
#                       if (str_detect(..1, pattern = "\\.tsv")) {
#                         df_map <- read_tsv(..1) 
#                         
#                         #rename the column to be able to join it
#                         colnames(df_map)[..3] <- "ID"
#                         
#                         #reparate for joining
#                         df_map <- df_map %>% 
#                           mutate(rowN = row_number()) %>% 
#                           separate_rows(ID, sep = IDsep)
#                         
#                       } else if (str_detect(..1, pattern = "\\.xlsx")) {
#                         df_map <- readxl::read_xlsx(..1) 
#                         
#                         #rename the column to be able to join it
#                         colnames(df_map)[..3] <- "ID"
#                         
#                         #reparate for joining
#                         df_map <- df_map %>% 
#                           mutate(rowN = row_number()) %>% 
#                           separate_rows(ID, sep = IDsep)
#                       }
#                       
#                       colnames(df_map)[-..3] <- paste0(colnames(df_map)[-..3], "_", ..2)
#                       
#                       return(df_map)
#                       
#                       
#                     }) %>%
#     purrr::reduce(full_join, by = "ID")
#   
#   #filter out proteins that are not found in all repliacates
#   data_sumNAs <- allJoined %>% 
#     dplyr::select(matches("rowN"))
#   summedNAs <- apply(data_sumNAs, 1, sumNAs)
#   found_in_all_replicates <- summedNAs == 0
#   
#   allJoined <- allJoined[found_in_all_replicates, ]
# 
#   write_tsv(allJoined, outputPath)
#   
# }
# 
# #funkce ktera po lopit merge da ID k sobe co maji stejny intensity
# lopit_merge_replicates_clean <- function(inputMerged, 
#                                          intesity_cols,
#                                          outputPath) {
#   #inputMerged = data created by lopit_merge_replicates
#   #intesity_cols = the column numbers of the merged dataset
#   #outputPath = ful output path for the joined dataset
#   
#   inputMerged_T <- read_tsv(inputMerged)
#   
#   inputMerged_T <- inputMerged_T %>% 
#     group_by(across(all_of(names(inputMerged_T)[intesity_cols]))) %>% 
#     mutate(ID = paste0(ID, collapse = ";")) %>% 
#     unique()
#   
#   write_tsv(inputMerged_T, outputPath)
# }

#funkce ktera promitne dany dataset na plot tsne nebo pca
lopit_project_datasets <- function(lopit_analyse_results, 
                                   otherData,
                                   IDcols, 
                                   showCol, 
                                   file_out) {
  #lopit_analyse_results = the tsv returned by the lopit_analyse function that contains data for the reduced dataset (t-sne or pca)
  #otherData = the datat to show there
  #IDcols = positions of columns in the first and second dataset to know how to merdge e.g. c(1,2), first number is for the lopit_analyse_results and the second for otherData
  #showCol = by which column to colour, should be characters or factor, this column is in the otherData
  #file_out = path to output the png
  
  #rename the columns for easy handling
  colnames(lopit_analyse_results)[IDcols[1]] <- "ID"
  colnames(otherData)[IDcols[2]] <- "ID"
  colnames(otherData)[showCol] <- "colour"
  
  #join the data
  plot_data <- left_join(lopit_analyse_results, otherData, by = "ID")
  onlySelected_T <- plot_data %>% 
    filter(!is.na(colour))
  
  ## plot the non-interactive plot ##
  tSNE_plot <- ggplot() +
    geom_jitter(data = plot_data, aes(Dim1, Dim2, color = colour), size = 2, alpha = 0.5) +
    geom_point(data = onlySelected_T, aes(Dim1, Dim2), color = "black", size = 3)+
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
  #pred_res = full path to the tsv that is outputed by the lopit_analyse function
  #DirOut = Directory where to save the plots of all the markers
  #quantiCols = the position of the quantity columns in the output dataset
  #the number of replicates that is in the file
  
  #load the data
  res_T <- read_tsv(pred_res)
  
  #get the markers
  res_markers_T <- res_T %>% 
    dplyr::filter(markers != "unknown")
  
  #factor the replicates
  levels_rep <- colnames(res_markers_T)[quantiCols]
  
  #create repliecate groups
  fract_per_replicate <- length(levels_rep)/replicate_count
  replicate_groups_Sch <- paste0(levels_rep, "_x_", rep(1:replicate_count, each = fract_per_replicate))
  colnames(res_markers_T)[quantiCols] <- replicate_groups_Sch
  
  #plot the pattern of each marjer group
  map(unique(res_markers_T$markers),
      ~{
        #remove any special characters
        outName <- str_remove_all(.x, "[^a-zA-Z0-9_]")
        
        #filter the data
        res_markers_map_T <- res_markers_T %>% 
          filter(markers == .x) %>% 
          gather(all_of(quantiCols), key = "fraction", value = "intensity") %>% 
          mutate(intensity = as.double(intensity))
        
        #make fraction into a factor
        res_markers_map_T$fraction <- factor(res_markers_map_T$fraction, levels = replicate_groups_Sch )
        
        #make group per ID per replicate 
        res_markers_map_T <- res_markers_map_T %>% 
          group_by(ID,fraction ) %>% 
          mutate(rep = str_split(fraction, "_x_", simplify = T)[,2],
                 group = paste0(ID, "_", rep)
                 )
        
        #plot the data
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
  #pred_output = path to the perditopn results from the SVM or MCMC prediction, SVMres_your_suffix.Rdata or MCMCres_your_suffix.Rdata
  #method = method to be used for dimensional reduction, possible: "t-SNE", "PCA", "UMAP" and "lda
  #colour_col = by which coloumn to colour the dataset, if NULL it is recognized wheter the results are MCMC or SVM and the proper colouring is calculated with the provided cutoff
  #filter_colour_col = which column to use for filtering of the prediction, ale pozor je to global cutoff hlavni funkce ted dela per compartment v pripade SVM, ostatni jsou global (TAGM-MAP a MCMC)
  #filter_cutoff = either the probability cutoff (minimum 0.99) if MCMC or t cutoff (minimum 0.5) if SVM
  
  #load the data
  load(pred_output)

  # Find the variable name
  var_name <- ls()[str_detect(ls(), "res")][1]  # [1] takes first match if multiple
  
  # Rename it (e.g., to "my_results")
  assign("pred_res", get(var_name))
  rm(list = var_name)  # Remove original
  
  if (is.null(filter_colour_col) & is.null(filter_cutoff)) {
    #if null just colour by the provided column
    set.seed(seed)
    pRolocVis(object = pred_res,
              method = method,
              fcol = colour_col)
    cat("here")
  } else if (!is.null(filter_colour_col) & !is.null(filter_cutoff)) {
    # get the data
   data_T <- fData(pred_res)
   
   # renname the filter column
   filter_col <- which(colnames(data_T) == filter_colour_col)
   colnames(data_T)[filter_col] <- "filter_col"
   
   #prediction col rename
   pred_col<- which(colnames(data_T) == colour_col)
   colnames(data_T)[pred_col] <- "pred_col"
    
   #create the prediction filter
   data_T <- data_T %>% 
     mutate(pred_filter = if_else(filter_col >= filter_cutoff,
                                 pred_col,
                                 NA_character_)
            )
   
   # add the new filtered prediction
   fData(pred_res)$pred_filter <- data_T$pred_filter
   
   # make the graph
   set.seed(seed)
   pRolocVis(object = pred_res,
             method = method,
             fcol = "pred_filter")
  } else {
    errorCondition("Ether you provied both filter_colour_col and filter_cutoff or provied only the column by which to colour.")
  }
  

}

lopit_resolution <- function(pred_output,
                             marker_col = "markers",
                             DirOut) {
  #pred_output = path to the perditopn results from the SVM or MCMC prediction, SVMres_your_suffix.Rdata or MCMCres_your_suffix.Rdata
  #marker_col = name of the markers column
  #DirOut = Directory where to save the plots
  
  #load the data
  load(pred_output)
  
  # Find the variable name
  var_name <- ls()[str_detect(ls(), "(svmres|tagm_map_res|mcmcRes)")][1]  # [1] takes first match if multiple
  
  # Rename it (e.g., to "my_results")
  assign("pred_res", get(var_name))
  rm(list = var_name)  # Remove original
  
  #qsep function
  qsep_res <- QSep(object = pred_res, fcol = marker_col)
  
  #create the directory
  dir.create(DirOut)
  
  ## Plot a boxplot of raw QSep scores
  png(paste0(DirOut, "/", "box_plot.png"))
  plot(x = qsep_res, norm = FALSE)
  dev.off()
  
  ## Plot a boxplot of normalised QSep scores
  png(paste0(DirOut, "/", "box_plot_norm.png"))
  plot(x = qsep_res, norm = TRUE)
  dev.off()
  
  ## Plot a level plot of raw QSep values between each compartment
  png(paste0(DirOut, "/", "heatmap.png"))
  plot <- levelPlot(object = qsep_res, norm = FALSE)
  plot(plot)
  dev.off()
  
  ## Plot a level plot of normalised QSep values between each compartment
  png(paste0(DirOut, "/", "heatmap_norm.png"))
  plot <- levelPlot(object = qsep_res, norm = TRUE)
  plot(plot)
  dev.off()
  
}

#this funtion merges datasets a lopit_merge_replicates, but importantly returns MSnSet dataset (important for bundle MCMC prediction)
#the data mast have the same ID types e.g. uniprot!!!!
lopit_combined_data_for_bundle <- function(list_data_path,
                                           data_names,
                                           outDir, #full path to the output dir "my/output/dir"
                                           outSuffix #suffix of the output files
                                           ) {
  #list_data_path = list of paths of datasets that will be joined, the datasets are in MSnSET format, results for lopit_analyse when you run no prediction
  #data_names = the names (a vecore) of thel oaded data, (how is it named in the R enviroment when load with the load() function)
  data_list <- vector(mode = "list", length = length(list_data_path))
  
  #load the datasets
  for (i in 1:length(list_data_path)) {
    # Load the data
    load(list_data_path[i])
    
    # Find the variable name
    var_name <- ls()[str_detect(ls(), paste0("^", data_names[i], "$"))]  # Takes first match if multiple
    
    # Rename it (e.g., to "data_1", "data_2", etc.)
    assign(paste0("data_", i), get(var_name))
    rm(list = var_name)  # Remove original
    
    # assing it to the list
    data_list[[i]] <- get(paste0("data_", i))
  }

  #make the bundle input
  data <- commonFeatureNames(data_list)
  
  #return the input
  save(data, 
       file = paste0(outDir, "/Bundle_joined_data_", outSuffix, ".Rdata")
       )
  
}

#function to test different perplexities
lopit_try_perplexities <- function(data,
                                   perplexities = seq(10, 100, by = 10),
                                   outDir, #full path to the output dir "my/output/dir"
                                   outSuffix, #suffix of the output files,
                                   seed = 42
                                   )  {
  #data path to a MSnSet dataset, the data prepared for prediction

  #load the data
  load(data)
  
  #create diffferent perplexiti plots
  map(perplexities,
      ~{
        #set seed
        set.seed(seed)
        
        #run t-SNE
        pred_filtered_plot <- plot2D(data,
                                     method = "t-SNE",
                                     methargs = list(perplexity = .x)
                                     )
        
        #prepare data for the plot
        pred_filtered_plot_T <- tibble(ID = rownames(pred_filtered_plot),
                                       Dim1 = pred_filtered_plot[,1],
                                       Dim2 = pred_filtered_plot[,2])
        
        
       
        ## plot the non-interactive plot ##
        tSNE_plot <- ggplot(pred_filtered_plot_T, aes(Dim1, Dim2)) +
          geom_jitter() +
          theme_minimal(base_size = 20) +
          labs(title = paste0("perplexity:", .x))
        
        #dopsat at se to uklada kam chces
        png(paste0(outDir, "/t-SNE", "_perplexity", .x, outSuffix, ".png"),
            width = 4000,
            height = 2500,
            res = 300
            )
        
        plot(tSNE_plot)
        
        dev.off()
        
      })
  
}


# odzkouset !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# odzkouset !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# odzkouset !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# odzkouset !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
lopit_mcmc_show_localization_distribution <- function(mcmcRes_converged_pooled_path,
                                                      ID,
                                                      outDir) {
  #mcmcRes_converged_pooled_path = full path to the results from the MCMC results that have mcmcRes_converged_pooled_for_violin_plots at the begininng
  #ID = the ID of the protein to show
  #outDir = full path of the directory where to put the results
  #load the data
  load(mcmcRes_converged_pooled_path)
  
  #plot the data
  png(paste0(outDir, "/", ID, ".png"))
  
  plot(mcmcRes_converged_pooled, ID)
  
  dev.off()
}
# odzkouset !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# odzkouset !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# odzkouset !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# odzkouset !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# a function the selects point in a circle you chose from the t-sne plot
lopit_points_in_circle <- function(data_path, 
                                   center_x, 
                                   center_y,
                                   radius,
                                   x_col = "Dim1", 
                                   y_col = "Dim2",
                                   output_path) {
  
  #output_path = output path of the xlsx file and the graph, no suffix (like no .xlsx)
 
  # load the tsv data you get from lopit_analyse
  data_T <- read_tsv(data_path)
  
  # rename the dimensions variables
  x_pos <- which(colnames(data_T) == x_col)
  y_pos <- which(colnames(data_T) == y_col)
  
  colnames(data_T[c(x_pos, y_pos)]) <- c("Dim1", "Dim2")
  
  # Filter points inside the circle using vectorized distance calculation
  dist <- sqrt((data_T[[x_col]] - center_x)^2 + (data_T[[y_col]] - center_y)^2)
  inside <- data_T[dist <= radius, , drop = FALSE]
  
  # Create plot with all points, highlighted inside points, and circle
  p <- ggplot(data_T, aes(Dim1, Dim2)) +
    geom_point(alpha = 0.6, color = "gray50") +
    geom_point(data = inside, aes(Dim1, Dim2), 
               color = "red", size = 2) +
    geom_circle(aes(x0 = center_x, y0 = center_y, r = radius), 
                fill = NA, color = "blue", linewidth = 1, inherit.aes = FALSE) +
    coord_fixed() +  # Equal aspect ratio for accurate circle
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
  
  return(invisible(inside))  # Returns filtered data (invisible to avoid cluttering console)
}

#FUNCTION FOR LOPIT ANAYLSIS
lopit_analyse <- function(data = NULL, #full path to the data input "this/is/my/data.xlsx, must be a xlsx file or Rdata file of the output of the data prepared for the analysis
                          descLopit_data = NULL, #if you provide the data preparaded for analysis you have to also provided the path for descLopit_data (description created during the data preparation)
                          markers = NULL, #full path to the markers file, in excel, format: 3 columns, with names "ID", "markers" and "colors" (first is ID, second the marker localization, 3rd is colours, simply coloured cells), keep the names of the localisation constant!
                                          #small complexes should not be used as markers but they should be later added to the dimensionaly reduced data to see if they are next to each other as expected
                          MS = "TNT", # the possible values are "TNT" or "LFQ", use LFQ if combining LFQ and TMT data
                          filterRows = 18, #if MS is LFQ then if the provided number is higher or equal to the number of missing values in a row (per protein) the row is removed
                          minTNTval = 50, #if the value of TNT intensity is lower then this it is treated as NA and filtered out based on the filterRows variable
                          transform = NULL, #possible is: "log2" or "scale"
                          scale = c(0, 1000), #if the previous parameter is scale you need to include the upper and lower bond of the scale to transform to
                          impute = F, #T or F whether to impute or not, if F, 0 is put in the place of all NAs, if T, imputation with provided mean and sd is done (missing data must be NA or NaN)
                          imputeMean = NULL, #impute the NAs by a xxx distribution with given mean and sd, if not provided the data are imputed by normal distribution with sd 0.3 of the data distribution and with the mean shifted by 1.8 sd of the data distribution
                          imputeSd = NULL, #impute the NAs by a xxx distribution with given mean and sd, if not provided the data are imputed by normal distribution with sd 0.3 of the data distribution and with the mean shifted by 1.8 sd of the data distribution
                          outDir, #full path to the output dir "my/output/dir"
                          outSuffix, #suffix of the output files
                          IDsCol, #position of the column with IDs in the data
                          quantiCols, #positions of the intensitiy columns, either as: c(5:15), or c(1,2,3,4,5,6, ...)
                          quantiColsGroups = NULL, #if more datasets are used, you can create groups of them so when you apply scaling it is done per dataset not per all of the intensities, e.g. you have 3 datasets with 3 intesities column each, than you have to input: c(1,1,1,2,2,2,3,3,3), the numbers say what LFQs are in what groups
                          descCol, #posistion of descreption of the proteins (will be added to the interactive graph)
                          metaCols, #other metacolumns you want to add and will be present in the output file
                          IDs_sep = ";", #character separating IDs in the ID column
                          #normalizace je vhodna kdyz kombnujes data zlepsuje vysledky
                          norm = NULL, #possible is: "sum", "max", "quantiles", "quantiles.robust", "center.mean", "center.median", "div.mean", "div.median", "diff.median", "quantiles", "vsn",
                          predType = "SVM", #use SVM, MCMC, TAGM-MAP (supervised with markers), PHEATMAP, HDBSCAN (unsupervised without markers) prediction or "NO" (only the data are prepared)
                          
                          ##########################
                          ## unsupervised methods ##
                          ###################################
                          ## varible for pheatmap clusters ##
                          ###################################
                          #the resulrs should not give more than 15 clusters usualy you cannot s´distiquidśh more comparmtents
                          max_cluster_size = 500, #maximim number of proteins in the group
                          iteration = 100, #iteration that run through the dendogram, no need to change unless it complains it is too low than it needs to ne icreased
                          
                          ##################################
                          ## varible for hdbscan clusters ##
                          ##################################
                          minPts = 10, #normaly between 10-15
                          
                          #######################
                          ## suprvised methods ##
                          #######################
                          # variable for SVM #
                          ####################
                          sigma = NULL, #parameter for prediction, estimated if not provided
                          cost = NULL, #parameter for prediction,  estimated if not provided,
                          t_cutoff = 0.75, #cut off of the svm.score to consider a protein as predicted to be in certain compartment (if higher it is predicted to be in the given compartment)
                          t_cutoff_per_compartment = T, #if True compartment specific t treeshold is used, if true t_cutoff means qunatile cutoff not specific t value
                         
                          #########################
                          # variable for TAGM-MAP #
                          #########################
                          numIter_tagmmap = 100, #the number of times the programme samples from posterior
                          p_cutoff_tagmmap = 0.99, #them mean probability used to colour the prediction if the probability is higher than the treshold it is coloured otherwise it is grey
                          
                          ####################
                          # variable forMCMC #
                          ####################
                          #Class imbalance can lead to over-classification = trosku osetrono ze overall_probability je slozeni pravdepodobnosti lokalizace, outlier a shannon entropy
                          numIter = 20000, #the number of times the programme samples from posterior
                          burnin = 10000, #the number of steps from the begining toi discard (it takes time until it gets to the posterior from the renadom defualt, so you have to discard the begining)
                          thin = 20, #remove autocorelation (I think it means take only evry 20th step, since teps right after one another correlate)
                          numChains = 6, #how many predictions (chains) to run
                          MCMCres_path = NULL, #put the path to the result of the MCMC analysis and it will continue with visualization and statistics, MUST then also provide path to the data (the data variable) ready for MCMC analysis output by the programme (names as data_for_prediction_{suffix}.Rdata)
                          p_cutoff = 0.99, #them mean probability used to colour the prediction if the probability is higher than the treshold it is coloured otherwise it is grey
                          hyppar = c(0.5, 3, 100), #for BUNDLE priors for the hyper paramteres
                          dirPrior = 0.0005, #for BUNDLE Dirichlet prior matrix, how much you expect that 1 protein will have different localization between datsets, not much
                          conditions,#for BUNDLE, which datasets are which condition, based on how the data were merge with lopit_combined_data_for_bundle, the same input as quantiColsGroups (e.g. c(1,1,2,2))
                          
                          ###########################################################
                          # what chains to keep = important for downstream analysis #
                          ###########################################################
                          #you can select to keep the chains according to the "geweke", "gelman", "geweke_and_gelman" (these are statistics)
                          #alternatively, look at the result of these statistics and select manualy by providing a vector of chains (e.g. c(1,2,3)) to keep
                          
                          MCMCres_keepChains = "gelman", 
                          
                          extra_thin = NULL, #MCMC only (not Bandle), a number that does extra thining of the sample from the posterio, reduces auto correlation
                          extra_burnin = NULL, #MCMC only (not Bandle), a number of cycles to be removed from the beging of the chains, extra to the burnin already provided
                          geweke_pval = 0.05, #testing whether in a single hcain the first 10 percent is the same as the last 50, assesses converdgence, if significant the chains did not converdge and the chain is removed
                          
                          ##############################
                          # variable for visualization #
                          ##############################
                          seed = 42, #to produce always the same graphs from tsne or pca,
                          dimRed = "t-SNE",#what dimensionality reduction algorithm to use "PCA" or "t-SNE"
                          perplexity = 30, #, t-SNE parameter, which can be though of as the number of effective neighbors of every point
                                          # it is goo practise to try several perplexity values, 10, 20, 30, 40 and 50 at least
                          dims = c(1, 2) #which PC for PCA to use default is the first 2
                          ) {
  
  if (str_detect(data, pattern = "\\.xlsx")) {
    cat("Excell provided running the analysis from scratch.\n")
    ####################################
    #### change inputs and set seed ####
    ####################################
    #### remove backslash at the and of the path if present ####
    outDir <- str_remove(outDir, "/$")
    transform <- ifelse(is.null(transform),
                        "no",
                        transform
                        )
    
    #### read in the data ####
    lopit_T <- readxl::read_xlsx(data) %>% 
      mutate_at(quantiCols, as.double)
    
    #########################################
    #### separate the ID cols by IDs_sep ####
    #########################################
    colnames(lopit_T)[IDsCol] <- "ID"
 
    lopit_T <- lopit_T %>% 
      mutate(rowN = row_number()) %>% 
      separate_rows(ID, sep = IDs_sep)
    
    #create description
    descLopit_data <- lopit_T[,c(IDsCol,descCol, ncol(lopit_T))]

    #create intensities
    intensities <- lopit_T[,quantiCols]

    #### if TNT put the NA for the proteins with low intensity
    if (MS == "TNT") {
      intensities <- intensities %>% 
        mutate_all(~if_else(. <= minTNTval,
                            NA_real_,
                            .)
                   )
    }
    
    #sum NAs and put it to the description
    NAs_n <- apply(intensities, 1, sumNAs)
    descLopit_data$NAs_n <- NAs_n
    
    colnames(descLopit_data) <- c("ID", "desc", "rowN", "NAs_n")
    
    ########################################
    #### keep prvnich pet slov z popisu ####
    ########################################
    keep5_words <- function(vec) {
      fristdesc <- str_split(vec, ";", simplify = T)[,1]
      words <- unlist(str_split(fristdesc, "\\s"))
      
      if (length(words) > 5) {
        words <- words[1:5]
      }
      
      words <- paste0(words, collapse = " ")
      
      return(words)
    }
    
    descLopit_data$desc2 <- lapply(descLopit_data$desc, keep5_words) %>% 
      unlist()
    
    #put back the IDs of the desc dataset
    descLopit_data <- descLopit_data %>% 
      group_by(rowN) %>% 
      mutate(ID = paste0(ID, collapse = IDs_sep)) %>% 
      sample_n(size = 1)
    
    write_tsv(descLopit_data, paste0(outDir, "/", "data_for_prediction_descLopit_data", outSuffix, ".tsv") )
    
    #####################################################################
    ### filter out the rows with to many NAs controled by filterRows ####
    #####################################################################
    lopit_T <- lopit_T[NAs_n < filterRows,]
    
    if (transform ==  "log2") {
      ##################################
      ## transforming the data by log ##
      ##################################
      lopit_T <- lopit_T %>%
        mutate_at(quantiCols, ~if_else(is.na(.) | . == 0,
                                       NA_real_,
                                       log2(.)
                                       )
                  )
      
      cat("Data were log2 transformed.\n")
    } else if (transform == "scale" & length(scale) == 2 & all(is.numeric(scale))) {
      
      ########################################
      ## make a rescale function for apply ###
      #######################################
      intensities <- lopit_T[,quantiCols]
      
      ##########################
      #### rescale the data ####
      ##########################
      if (is_null(quantiColsGroups)) {
        #if groups not provided scale all intensities as one
        scaled_intensities <- apply(intensities, 1, rescale_to_xx, scale) %>% 
          t() %>% 
          as_tibble()
        
        } else if (typeof(quantiColsGroups) == "double") {
          #if group provided do the scaling in the groups separatelly
        quantiColsGroups_L <- map(unique(quantiColsGroups),
                                  ~{
                                    which(quantiColsGroups == .x)
                                  })
      
        scaled_intensities <- map(quantiColsGroups_L,
                                  ~{
                                    scaled_map_T <- apply(intensities[,.x], 1, rescale_to_xx, scale) %>% 
                                      t() %>% 
                                      as_tibble()
                                    
                                    return(scaled_map_T)
                                  }) %>% 
          purrr::reduce(cbind)
        
        } else {
        errorCondition("Wrong input for quantiColsGroups.\nMust be NULL or a vector of numbers spliting the intesities into groups:\ne.g. c(1,1,1,2,2,2,3,3,3)")
        }
      
      lopit_T[,quantiCols] <- scaled_intensities
      
      cat("Data were scaled/transformed.\n")
    } else if (transform == "no") {
      cat("Data were NOT transformed.\n")
    } else {
      stop('You used wrong input for transformation. Possible inputs are NULL, \"log2\" or \"scale\".\nPossibly, the input for scale is wrong, it has to be vector of 2 numbers')
    }
    
    ################
    ## imputation ##
    ################
    if (impute == F) {
      
      ########################
      ## replace NAs with 0 ##
      ########################
      lopit_T <- lopit_T %>% 
        mutate_at(quantiCols, ~if_else(is.na(.),
                                       0,
                                       .),
                  )
      
      cat("Zeros were placed instead of NAs.\n")
    } else if (impute == T & !is.null(imputeMean) & !is.null(imputeSd)) {
      ## impute with normal distribution
      lopit_T <- lopit_T %>%
        group_by(ID) %>%
        mutate_at(quantiCols, ~if_else(is.na(.),
                                       rnorm(1, imputeMean, imputeSd),
                                       .)
                  )
      
      cat("Missing values were imputed with normal distribtuion with user specified parameters.\n")
    } else if (impute == T & xor(is.null(imputeMean), is.null(imputeSd))){
      stop("You provided only mean or sd for the imputation. Provided both and rerun the function!\n")
    } else if (impute == T & is.null(imputeMean) & is.null(imputeSd)) {
      
      #####################################################################################################################################
      ## impute with normal distribution with default para: mean shifted by 1.8 sd of the data distribution abd sd is 0.3 of the data sd ##
      #####################################################################################################################################
      ## calculate mean and sd of the data ##
      meanData <- mean(unlist(lopit_T[,quantiCols]), na.rm = T)
      sdData <- sd(unlist(lopit_T[,quantiCols]), na.rm = T)
      
      ## impute ##
      lopit_T <- lopit_T %>%
        group_by(ID) %>%
        mutate_at(quantiCols, ~if_else(is.na(.),
                                       rnorm(1, (meanData -1.8*sdData), sdData*0.3),
                                       .)
                  )
      cat("Missing values were imputed with normal distribtuin with default parameters:\n--mean shifted by 1.8 sd of the data\n--sd is 0.3 of the data sd.\n")
    }
    
    #####################
    ## add the markers ##
    #####################
    if (str_detect(predType, "(HDBSCAN|PHEATMAP|NO)")) {
      #if no prediction wanted create mock markers
      markers_T <- tibble(ID = sample(lopit_T$ID, size =  100, replace = F),
                          markers = "random") %>% 
        unique()
      
    } else if (str_detect(predType, "(MCMC|SVM|TAGM-MAP)")){
      markers_T <- readxl::read_xlsx(markers)
      
      ########################################################
      ## check for colour in the marker dataset, 3rd column ##
      ########################################################
      if (any(is.na(colnames(markers_T)[3] != "colors") | colnames(markers_T)[3] != "colors", na.rm = T)){
        #create the levels of markers and their associated colors
        levels_T <- markers_T %>% 
          dplyr::select(markers) %>% 
          unique()
        
        ## Generate distinct colors if not colors column provided ##
        levels <- c(levels_T$markers, "unknown")
        colors <- rainbow(length(levels)-1 )
        colors <- c(colors, "grey70")
      } else {
        cat("Using the user provided colors from the 3rd column of marker list!\n")
        
        ## get the colours
        markers_job <- loadWorkbook(markers)
        markers_sheet1 <- getSheets(markers_job)[[1]]
        rows <- getRows(markers_sheet1)
        cells <- getCells(rows)
        styles <- sapply(cells, getCellStyle)
        colours_cols_L <- styles[str_detect(names(styles), "\\.3$")]
        colours_cols_L <- colours_cols_L[-1]
        colors_sch <- sapply(colours_cols_L, cell_color)
        
        #add it to the markers
        markers_T$colors <- colors_sch
        
        #create the levels of markers and their associated colors
        levels_T <- markers_T %>% 
          dplyr::select(markers, colors) %>% 
          unique()
        
        levels <- c(levels_T$markers, "unknown")
        colors <- unname(c(levels_T$colors, "grey70"))
        
      }
      
      
    } else {
      errorCondition("Wrong predType, provided. Options are: \"MCMC\", \"SVM\", \"TAGM-MAP\", \"HDBSCAN\", or \"PHEATMAP\"!\n")
      
    }
    
    #add the markers
    lopit_T <- left_join(lopit_T, markers_T, by = "ID")
    #################################################################################
    # make sure that the groups that have marker that it is propagated in the group #
    #################################################################################
    lopit_T <- lopit_T %>% 
      group_by(rowN) %>% 
      mutate(markers = paste0(unique(markers[!is.na(markers)]), collapse = ","))
    
    #retrun Error if some markers have more than 2 compartments
    multimarker_protein <- unique(str_count(lopit_T$markers, ","))

    if(length(multimarker_protein) > 1) {
      #return the duplicated ID
      
    dupID <- markers_T$ID[duplicated(markers_T$ID)]
      cat("These markers are duplicated:\n",
          paste0(dupID,
                 collapse = "\n"),
          "\n",
          sep = "")
      stop("Error: Some markers have more than 2 compartments. Check your marker list!\n")
    }

    #report how many markers is missing
    present_marker <- !is.na(lopit_T$markers)
    present_marker <- lopit_T$ID[present_marker]
    
    cat("These markers were not found in your dataset:\n",
        paste0(markers_T$ID[!(markers_T$ID %in% present_marker)],
               collapse = "\n"),
        "\n",
        sep = "")
    
    ## change the empty strings into unknown ##
    lopit_T$markers <- if_else(lopit_T$markers == "",
                               "unknown",
                               lopit_T$markers)
    
    ## return the IDs into the previous form ##
    lopit_T <- lopit_T %>% 
      group_by(rowN) %>% 
      mutate(ID = paste0(ID, collapse = IDs_sep)) %>% 
      sample_n(size = 1)
    
    ## get the quantification columns ##
    expression_T <- lopit_T[,c(IDsCol, quantiCols)]
    colnames(expression_T) <- c("ID", paste0("quan", 1:length(quantiCols)))
    
    ## get the metadata ##
    meta_T <- lopit_T[,c(IDsCol, metaCols)] 
    meta_T$markers <- lopit_T$markers
    
    ## create the table for what fraction is in what quantification ##
    fraction_T <- tibble(sampleNames = colnames(expression_T)[2:ncol(expression_T)],
                         Fractions = 1:length(colnames(expression_T)[2:ncol(expression_T)]))#jen automaticky tam neco davam, ale neni to potreba vubec
    
    ## write csvs so they can be read in by the readMSnSet function ##
    outQuan <- paste0(outDir, "/", "Quan_data", outSuffix, ".tsv")
    outMeta <- paste0(outDir, "/", "Meta_data", outSuffix, ".tsv")
    outFrac <- paste0(outDir, "/", "Frac_data", outSuffix, ".tsv")
    
    write_tsv(expression_T, outQuan)
    write_tsv(meta_T, outMeta)
    write_tsv(fraction_T, outFrac)
    
    ## read in the data in correct data structure ##
    data <- readMSnSet(exprsFile = outQuan,
                       featureDataFile = outMeta,
                       phenoDataFile = outFrac,
                       sep = "\t")
    
    #######################
    #### normalization ####
    #######################
    if (!is.null(norm)) {
      data <- normalise(data, method = norm)
      
      cat("Normalization by", norm, "was performed.\n")
    } else {
      cat("Normalization was NOT performed.\n")
    }
    
    ############################
    ## save the prepared data ##
    ############################
    save(data, file = paste0(outDir, "/", "data_for_prediction", outSuffix, ".Rdata"))

    } else if (str_detect(data, "\\.Rdata") & !is.null(descLopit_data)) {
      cat("Data and DescLopit were provided were provided, loading the data.\n")
      load(data)
      descLopit_data <- read_tsv(descLopit_data)
    } else {
      errorCondition("Data and DescLopit or an excel were not provided.\n")
    }

  ###########################
  ## what prediction to do ##
  ###########################
  if (predType == "SVM") {
    cat("Running SVM predictions.\n")
    
    #############
    ## run SVM ##
    #############
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
    ##############
    ## run MCMC ##
    ##############
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
      cat("Running unsuprived clustering with hclust via pheatmap.\n")
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
      cat("Running unsuprived clustering with hdbscan.\n")
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
      errorCondition("Wrong option for predType selected. Possiblilities are:\n\"SVM\", \"MCMC\", \"TAGM-MAP\" (supervised with markers), \"PHEATMAP\", \"HDBSCAN\" (unsupervised without markers) or \"NO\" (just prepare data)")
    }
  }

###############################################
## functions called within the lopit_analyse ##
###############################################
#function to get colors from xlsx
cell_color <- function(style) {
  fg <- style$getFillForegroundXSSFColor()
  hex <- tryCatch(fg$getRgb(), error = function(e) NULL)
  hex <- paste0("#", paste(hex, collapse = ""))
  tint <- tryCatch(fg$getTint(), error = function(e) NULL)
  
  if (!is.null(tint) & !is.null(hex)) {
    rgb_col <- col2rgb(col = hex)
    if (tint < 0) rgb_col <- (1 - abs(tint)) * rgb_col
    if (tint > 0) rgb_col <- rgb_col + (255 - rgb_col) * tint
    hex <- rgb(red = rgb_col[1,1], green = rgb_col[2,1], blue = rgb_col[3,1], maxColorValue = 255)
  }
  return(hex)
}


#function for scaling rowwise
rescale_to_xx <- function(vec, scale) {
  rescale(vec,
          from = c(min(vec, na.rm = T),
                   max(vec, na.rm = T)
          ),
          to = scale
  )
}

#function to sum NAs
sumNAs <- function(vec) {
  sum(is.na(vec))
}

#function to run gelma diagnostics on MCMC chains
run_gelman <- function(inputMatrices,
                       names, 
                       subset = NULL
                       ) {
  #inputMatrices = list of the results from the mcmc_get_ functions
  #names = names of the matrices in a vector
  #subset = a vector of numbers to use for subset of the mcmc_get_ function results
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

#makes summury of statistics and plots of the the given matrix
make_MCMC_stat_plots <- function(matrix,
                                 output,
                                 nChains
                                 ) {
  #output = output without extension
  #matrix = matrix returned by mcmc_get_ functions
  #nChains = number of chains that were run for the MCMC prediction
  
  #create sumary for each chain and save it
  for (i in seq_len(nChains)) {
    sink(paste0(output, ".txt"))
    print(summary(matrix[[i]]))
    sink()
  }
  
  #create plot for outliers
  pdf(file = paste0(output, ".pdf"))
  
  ## Using coda S3 objects to produce trace plots and histograms
  for (i in seq_len(nChains)) {
    
    plot(matrix[[i]], main = paste("Chain", i), auto.layout = FALSE, col = i)
  }
  
  dev.off()
}


make_NA_plot <- function(pred_filtered_plot_T,
                         outDir,
                         dimRed,
                         outSuffix) {
  # pred_filtered_plot_T = data created in the svm or MCMC function
  
  
  ## plot the non-interactive plot ##
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


#functions that clusters based on hdbscan
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
  ## Run HDBSCAN using optimal minimum cluster size - here 10
  hdb_results <- data %>%
    exprs() %>%
    hdbscan(minPts = 10)
  
  ## Add the results of HDBSCAN to our MSnSet
  fData(data)$hdb_cluster_id <- hdb_results$cluster
  fData(data)$hdb_cluster_prob <- hdb_results$membership_prob
  
  ## Check how many proteins are in each cluster
  data %>%
    fData() %>%
    pull(hdb_cluster_id) %>%
    table()
  
  ## run the dim reduction
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
    stop('Error: no method for dimensional reduction provided. Provided one of these: "PCA", "t-SNE" as the dimRed variable.\n')
  }
  
  pred_filtered_plot_T <- tibble(ID = rownames(pred_filtered_plot),
                                 Dim1 = pred_filtered_plot[,1],
                                 Dim2 = pred_filtered_plot[,2])
  
  #add the hdbscan clusters
  hdbscan_clusters_T <- fData(data) %>% 
    dplyr::select(hdb_cluster_id, hdb_cluster_prob) %>% 
    mutate(ID = rownames(.))

  ## add columns to the locations ##
  pred_filtered_plot_T <- left_join(pred_filtered_plot_T, hdbscan_clusters_T, by = "ID")
  
  ## add desc to the data ##
  pred_filtered_plot_T <- left_join(pred_filtered_plot_T, descLopit_data, by = "ID")
  
  ## create label ##
  pred_filtered_plot_T <- pred_filtered_plot_T %>% 
    group_by(ID) %>% 
    mutate(label = paste0(ID, ", dbscan_cluster: ", hdb_cluster_id, " - ", desc2, collapse = ""))
  
  ## change the clusters to factor ##
  pred_filtered_plot_T$hdb_cluster_id  <- factor(pred_filtered_plot_T$hdb_cluster_id , levels = 1:length(unique(pred_filtered_plot_T$hdb_cluster_id)))
  
  ## create only marker dataset ##
  markers_plot_T <- pred_filtered_plot_T
  
  ## Generate distinct colors ##
  colors <- rainbow(length(unique(pred_filtered_plot_T$hdb_cluster_id)), )
  
  ## plot the non-interactive plot ##
  tSNE_plot <- ggplot(pred_filtered_plot_T, aes(label = label)) +
    geom_jitter(data = pred_filtered_plot_T, aes(Dim1, Dim2, color = hdb_cluster_id), size = 2, alpha = 0.5) +
    theme_minimal(base_size = 20) +
    labs(color = "Localization")+
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


# funtion that takes in pheatmap object and returns the clusters ploted against a dim reduction plot
run_pheatmap_clustering <- name <- function(data, 
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
  #max_cluster_size = what is the maximum size of a cluster
  #iteration = how many steps done the dendogramt the programme should go to extract clusters
  #output_pdf = path where to output the pdf of the heatmap, do not add pdf it is added outomaticaly

  #################################
  ## make the cluster by heatmap ##
  #################################
  
  #run the pheatmap
  pheatmap_res <- pheatmap(mat = exprs(data),
                           cluster_rows = TRUE,
                           cluster_cols = FALSE,
                           show_rownames = FALSE,
                           filename = paste0(outDir, "/", dimRed, outSuffix, "_heatmap.pdf")
                           )
  
  #extract clusters
  row_hcluts <- pheatmap_res
  
  #extract the ordered IDs
  row_ordered_IDs <- row_hcluts$tree_row$labels
  
  #create the first cluster
  cluster_T <- tibble(ID = row_ordered_IDs,
                      cluster1 = 1
                      )
  
  #map over the clusters size and return all the possible cluster with wheter to keep the clusters or not
  otherClusters_T <- map(2:iteration,
      ~{
        #get the new cluster
        row_clusters <- cutree(row_hcluts$tree_row, k = .x)
        
        #create the current cluster with the given number of clusters
        clusters_T <- tibble(ID = row_ordered_IDs,
                             cluster = row_clusters) %>% 
          group_by(cluster) %>% 
          mutate(keep = n() > max_cluster_size) %>%  #decide wheter to keep certain clusters or not
          filter(keep) %>% 
          dplyr::select(-keep)
        
        colnames(clusters_T)[2] <- paste0(colnames(clusters_T)[2], "_", .x)
        
        return(clusters_T)
      }) %>% 
    purrr::reduce(full_join, by = "ID")
  
  #check the last iteration if all is na if not return that higher iteration is needed
  all_is_NA <- nrow(otherClusters_T) -sum(is.na(otherClusters_T[,ncol(otherClusters_T[,])]))
  
  if (all_is_NA != 0) {
    cat("Higher iteration is needed to resolve the groups!\n")
    return(NULL)
  }
  
  all_clusters_T <- left_join(cluster_T, otherClusters_T, by = "ID")
  
  #cluster based on the number of NAs
  all_clusters_T$total_clusters <- apply(all_clusters_T[2:ncol(all_clusters_T)], 1, sumNAs)
  
  #make the clusters into sensible numbers
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
                           annotation_row = cluster_heatmap_T,#colour the heatmap by clusters
                           filename = paste0(outDir, "/", dimRed, outSuffix, "_heatmap_clusters.pdf")
                           )

  #add ids a columns to join it for reduction
  cluster_heatmap_T$ID <- rownames(cluster_heatmap_T)

  ###################################################################
  #### create the grahp and get the location of the reduced data ####
  ###################################################################
  ## set seed ##
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
    stop('Error: no method for dimensional reduction provided. Provided one of these: "PCA", "t-SNE" as the dimRed variable.\n')
  }
  
  pred_filtered_plot_T <- tibble(ID = rownames(pred_filtered_plot),
                                 Dim1 = pred_filtered_plot[,1],
                                 Dim2 = pred_filtered_plot[,2])
  
  ## add columns to the locations ##
  pred_filtered_plot_T <- left_join(pred_filtered_plot_T, cluster_heatmap_T, by = "ID")
  
  ## add desc to the data ##
  pred_filtered_plot_T <- left_join(pred_filtered_plot_T, descLopit_data, by = "ID")
  
  ## create label ##
  pred_filtered_plot_T <- pred_filtered_plot_T %>% 
    group_by(ID) %>% 
    mutate(label = paste0(ID, ", pHeatMap_cluster: ", new_cluster, " - ", desc2, collapse = ""))

  ## change the clusters to factor ##
  pred_filtered_plot_T$new_cluster <- factor(pred_filtered_plot_T$new_cluster, levels = as.character(1:length(unique(pred_filtered_plot_T$new_cluster))))
  
  ## create only marker dataset ##
  markers_plot_T <- pred_filtered_plot_T
  
  ## Generate distinct colors ##
  colors <- rainbow(length(unique(pred_filtered_plot_T$new_cluster)), )
  
  ## plot the non-interactive plot ##
  tSNE_plot <- ggplot(pred_filtered_plot_T, aes(label = label)) +
    geom_jitter(data = pred_filtered_plot_T, aes(Dim1, Dim2, color = new_cluster), size = 2, alpha = 0.5) +
    theme_minimal(base_size = 20) +
    labs(color = "Localization")+
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

#SVM prediction
run_SVM <- function(data,
                    sigma, #parameter for prediction, estimated if not provided
                    cost, #parameter for prediction,  estimated if not provided,
                    seed, #to produce always the same graphs from tsne or pca,
                    t_cutoff, #cut off of the svm.score to consider a protein as predicted to be in certain compartment (if higher it is predicted to be in the given compartment)
                    t_cutoff_per_compartment,
                    dimRed,#what dimensionality reduction algorithm to use "PCA" or "t-SNE"
                    perplexity, #, t-SNE parameter, which can be though of as the number of effective neighbors of every point
                    dims,
                    levels, #levels of the markers so the colours provided by user match (created before this funtion is run in the main lopit_analyses function)
                    colors, #colours as found in the 3rd column of markers
                    outDir, #full path to the output dir "my/output/dir"
                    outSuffix, #suffix of the output files
                    descLopit_data #passed as a dataframe from the analyseLopit function
                    ) {
  
  #######################################################
  #### predict other proteins from markers using SVM ####
  #######################################################
  ## prediction using the default best sigma and cost (should be checked and added manually if needed) ##
  if (is.null(sigma) & is.null(cost)) {
    cat("Sigma and cost not provided. Performing their optimisation and subsequent localization prediction.\n")
    
    ## Set class weights as inverse of class frequencies
     # makes sure that when we have too many markers for 1 compartment it will not be predicting false those that have more markers
    
    ## Get markers ##
    marker_tbl <- data %>%
      getMarkers() %>%
      table()
    
    ## calculate the weights ##
    weights <- 1 / marker_tbl[names(marker_tbl) != "unknown"]
    
    ## optimization ##
    params <- svmOptimisation(data, 
                              fcol = "markers",
                              class.weights = weights)
    
    ## optimization graphs ##
    png(paste0(outDir, "/", "Dist_F1score_best_sigma_cost", outSuffix, ".png"),
        width = 1500, 
        height = 1500,
        res = 300)
    
    plot1 <- plot(params)
    print(plot1)
    
    dev.off()
    
    png(paste0(outDir, "/", "Avereged_F1score_all_sigma_cost.", outSuffix, ".png"),
        width = 1500, 
        height = 1500,
        res = 300)
    
    plot2 <- levelPlot(params)
    print(plot2)
    
    dev.off()
    
    ## save the data of optimization ##
    save(params, file = paste0(outDir, "/", "paraOptimization", outSuffix, ".Rdata"))
    
    ## clasification by default ##
    svmres <- svmClassification(data, fcol = "markers",
                                assessRes = params)
    
  } else if (!is.null(sigma) & !is.null(cost)) {
    cat("Using user specified sigma and cost to predict localization.\n")
    
    ## Get markers ##
    marker_tbl <- data %>%
      getMarkers() %>%
      table()
    
    ## calculate the weights ##
    weights <- 1 / marker_tbl[names(marker_tbl) != "unknown"]
    
    ## use custom provided sigma and cost ##
    svmres <- svmClassification(data, 
                                fcol = "markers",
                                sigma = sigma,
                                cost = cost, 
                                class.weights = weights)
    
  } else {
    stop("Error: either sigma or cost was not provided.")
  }
  
  #save the svm res
  save(svmres, file = paste0(outDir, "/SVMres", outSuffix, ".Rdata"))
  
  ## return what was done with the dataset (mainly sigma and cost) ##
  proccessing <- processingData(svmres)
  cat("Dataprocessing summary:", proccessing@processing, sep = "\n")
  
  #######################
  #### visualization ####
  #######################
  ## show how many prediction are returned with given cut off ##
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
  
  #########################
  ### supervised method ###
  #########################
  cat("Creating graphical outputs.\n")
  
  ## get the predictions with cutoff ##
  if (t_cutoff_per_compartment == F) {
    ## for all compartment a single t cutoff
    pred_filtered <- getPredictions(svmres,
                                    fcol = "svm", 
                                    t = t_cutoff)
    pred_filtered_T <- fData(pred_filtered)
    pred_filtered_T$ID <- rownames(pred_filtered_T)
    
  } else if (t_cutoff_per_compartment == T) {
    
    ## Get organelle-specific quantile SVM scores
    score_thresholds <- orgQuants(object = svmres,
                                  fcol = "svm",
                                  scol = "svm.scores",
                                  mcol = "markers",
                                  t = t_cutoff)
    
    ## Use organelle-specific quantiles to get thresholded localisation predictions
    pred_filtered <- getPredictions(object = svmres,
                                    fcol = "svm",
                                    scol = "svm.scores",
                                    mcol = "markers",
                                    t = score_thresholds)
    pred_filtered_T <- fData(pred_filtered)
    pred_filtered_T$ID <- rownames(pred_filtered_T)
    
  } else {
    errorCondition("t_cutoff_per_compartment must be T or F!\n")
  }
 
  
  ###################################################################
  #### create the grahp and get the location of the reduced data ####
  ###################################################################
  ## set seed ##
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
    stop('Error: no method for dimensional reduction provided. Provided one of these: "PCA", "t-SNE" as the dimRed variable.\n')
  }
  
  pred_filtered_plot_T <- tibble(ID = rownames(pred_filtered_plot),
                                 Dim1 = pred_filtered_plot[,1],
                                 Dim2 = pred_filtered_plot[,2])
  
  ## add columns to the locations ##
  pred_filtered_plot_T <- left_join(pred_filtered_plot_T, pred_filtered_T, by = "ID")
  
  ## add desc to the data ##
  pred_filtered_plot_T <- left_join(pred_filtered_plot_T, descLopit_data, by = "ID")
  
  ## create label ##
  pred_filtered_plot_T <- pred_filtered_plot_T %>% 
    group_by(ID) %>% 
    mutate(label = paste0(ID, ", marker: ", markers, ", svm: ", svm, ", svm.score: ", round(svm.scores, 2), " - ", desc2, collapse = ""))
  
  ## change the svm.pred to factor ##
  pred_filtered_plot_T$svm.pred <- factor(pred_filtered_plot_T$svm.pred, levels = levels)
  
  ## create only marker dataset ##
  markers_plot_T <- pred_filtered_plot_T %>% 
    filter(markers != "unknown")
  
  ## plot the non-interactive plot ##
  tSNE_plot <- ggplot(pred_filtered_plot_T, aes(label = label)) +
    geom_jitter(data = pred_filtered_plot_T, aes(Dim1, Dim2, color = svm.pred), size = 2, alpha = 0.5) +
    geom_point(data = markers_plot_T, aes(Dim1, Dim2), color = "black", size = 3)+
    geom_point(data = markers_plot_T, aes(Dim1, Dim2, color = svm.pred), size = 2) +
    theme_minimal(base_size = 20) +
    labs(color = "Localization")+
    scale_color_manual(values = colors)
  
  png(paste0(outDir, "/", dimRed, outSuffix, ".png"),
      width = 4000,
      height = 2500,
      res = 300
  )
  
  plot(tSNE_plot)
  
  dev.off()
  
  ## plot the NA plot
  make_NA_plot(pred_filtered_plot_T = pred_filtered_plot_T,
               dimRed = dimRed,
               outSuffix = outSuffix,
               outDir = outDir
               )
  
  ## plot the data interactively in plotly ##
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
  proteinsN <- length(pred_filtered_plot_T$svm)
  proteins_pred_with_t_cutoff <- pred_filtered_plot_T$svm.pred %>% 
    table() %>% 
    .[names(.) != "unknown"] %>% 
    sum()
  
  returnList <- list(data  = pred_filtered_plot_T)
  
  return(returnList)
  
  cat("Everything is done.\n")
}

#TAGM-MAP preidiction
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
  
  
  #######################################################
  #### predict other proteins from markers using SVM ####
  #######################################################
  ## prediction using the default best sigma and cost (should be checked and added manually if needed) ##
  ## optimization ##
  params <- tagmMapTrain(data,
                         numIter = numIter_tagmmap,
                         fcol = "markers")
  
  ## convergence graphs ##
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
  
  ## clasification by tagm-map ##
  tagm_map_res <- tagmMapPredict(data, 
                                 fcol = "markers",
                                 params = params)
    

  #save the svm res
  save(tagm_map_res, file = paste0(outDir, "/TAGM-MAPres", outSuffix, ".Rdata"))
  
  #######################
  #### visualization ####
  #######################
  ## show how many prediction are returned with given cut off ##
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
  
  #########################
  ### supervised method ###
  #########################
  cat("Creating graphical outputs.\n")
  
  ## Store classification and outlier probabilities
  tagm_prob <- fData(tagm_map_res)[, "tagm.map.probability"]
  tagm_out <- 1 - fData(tagm_map_res)[, "tagm.map.outlier"]
  
  ## Create new column containing overall probability
  fData(tagm_map_res)[, "overall_prob"] <- tagm_prob * tagm_out
  ## Set prediction thresholds on overall probability
  pred_filtered <- getPredictions(tagm_map_res,
                                  fcol = "tagm.map.allocation",
                                  scol = "overall_prob",
                                  t = p_cutoff_tagmmap)

  pred_filtered_T <- fData(pred_filtered)
  pred_filtered_T$ID <- rownames(pred_filtered_T)
  
  ###################################################################
  #### create the grahp and get the location of the reduced data ####
  ###################################################################
  ## set seed ##
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
    stop('Error: no method for dimensional reduction provided. Provided one of these: "PCA", "t-SNE" as the dimRed variable.\n')
  }
  
  pred_filtered_plot_T <- tibble(ID = rownames(pred_filtered_plot),
                                 Dim1 = pred_filtered_plot[,1],
                                 Dim2 = pred_filtered_plot[,2])
  
  ## add columns to the locations ##
  pred_filtered_plot_T <- left_join(pred_filtered_plot_T, pred_filtered_T, by = "ID")
  
  ## add desc to the data ##
  pred_filtered_plot_T <- left_join(pred_filtered_plot_T, descLopit_data, by = "ID")
  
  ## create label ##
  pred_filtered_plot_T <- pred_filtered_plot_T %>% 
    group_by(ID) %>% 
    mutate(label = paste0(ID, ", marker: ", markers, ", allocation: ", tagm.map.allocation, ", prob: ", round(overall_prob, 4), " - ", desc2, collapse = ""))

  ## change the tagm.map.allocation.pred to factor ##
  pred_filtered_plot_T$tagm.map.allocation.pred <- factor(pred_filtered_plot_T$tagm.map.allocation.pred, levels = levels)
  
  ## create only marker dataset ##
  markers_plot_T <- pred_filtered_plot_T %>% 
    filter(markers != "unknown")

  ## plot the non-interactive plot ##
  tSNE_plot <- ggplot(pred_filtered_plot_T, aes(label = label)) +
    geom_jitter(data = pred_filtered_plot_T, aes(Dim1, Dim2, color = tagm.map.allocation.pred), size = 2, alpha = 0.5) +
    geom_point(data = markers_plot_T, aes(Dim1, Dim2), color = "black", size = 3)+
    geom_point(data = markers_plot_T, aes(Dim1, Dim2, color = tagm.map.allocation.pred), size = 2) +
    theme_minimal(base_size = 20) +
    labs(color = "Localization")+
    scale_color_manual(values = colors)
  
  png(paste0(outDir, "/", dimRed, outSuffix, ".png"),
      width = 4000,
      height = 2500,
      res = 300
  )
  
  plot(tSNE_plot)
  
  dev.off()
  
  ## plot the NA plot
  make_NA_plot(pred_filtered_plot_T = pred_filtered_plot_T,
               dimRed = dimRed,
               outSuffix = outSuffix,
               outDir = outDir
               )
  
  ## plot the data interactively in plotly ##
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
  proteinsN <- length(pred_filtered_plot_T$tagm.map.allocation)
  proteins_pred_with_p_cutoff <- pred_filtered_plot_T$tagm.map.allocation.pred %>% 
    table() %>% 
    .[names(.) != "unknown"] %>% 
    sum()
  
  returnList <- list(data  = pred_filtered_plot_T)
  
  return(returnList)
  
  cat("Everything is done.\n")
}

#MCMC function
#zatim bez moznosti ovlivnovat prior (pokud to chces dodelat precti si https://f1000research.com/articles/8-446/v1)
#pokud budu chtit nekdy dodavat prior je to pomoci S0 argument
#DODELAT VIZULALIZACI
run_MCMC <- function(data,
                     numIter, #the number of times the programme samples from posterior
                     burnin, #the number of steps from the begining toi discard (it takes time until it gets to the posterior from the renadom defualt, so you have to discard the begining)
                     thin, #remove autocorelation (if e.g. 20 then every 20th step is kept), it is to remove autocorrelation since steps right after one another correlate)
                     numChains, #number of chains to run
                     MCMCres_path, #put the path to the result of the MCMC analysis and it will continue with visualization and statistics,
                     MCMCres_keepChains, #a vector of chains (e.g. c(1,2,3)) to keep and skip the statistical calculation of the chain convergance             
                     extra_thin, #a number that does extra thining of the sample from the posterio, reduces auto correlation
                     extra_burnin, #a number of cycles to be removed from the beging of the chains, extra to the burnin already provided
                     geweke_pval, #testing whether in a single hcain the first 10 percent is the same as the last 50, assesses converdgence, if significant the chains did not converdge and the chain is removed
                     seed,
                     dimRed,#what dimensionality reduction algorithm to use "PCA" or "t-SNE"
                     perplexity, #, t-SNE parameter, which can be though of as the number of effective neighbors of every point
                     dims,
                     levels,
                     colors, #colours as found in the 3rd column of markers
                     p_cutoff,
                     outDir, #full path to the output dir "my/output/dir"
                     outSuffix, #suffix of the output files
                     descLopit_data #passed as a dataframe from the analyseLopit function
                     ) {
  
  #######################
  ## run the MCMC pred ##
  #######################
  if (is.null(MCMCres_path)) {
    cat("Running the MCMC from scratch since no results of the MCMC run were provided.\n")
    # prepare mulcicore for parallel computing
    n_workers <- as.integer(Sys.getenv("PBS_NCPUS", unset = "4"))  # default 4 if PBS_NCPUS not set (it should always be, it is the number of CPUS allocated)
    
    multicoreParam <- MulticoreParam(workers = n_workers) #set number of cpus??
    register(multicoreParam)
    
    #run MCMC with tryCatch to rerun if chain fails (possible it is just random the fails)
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
        
        success <- TRUE  # if this runs without error, mark success
      }, error = function(e) {
        message("Error occurred in MCMC run: ", e$message)
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
    cat("Results of the MCMC run were provided, loading them now.\n")
    load(MCMCres_path)
  }
  
  if (!is.null(extra_thin) && !is.null(extra_burnin)) {
    mcmcRes <- mcmc_thin_chains(mcmcRes, freq = extra_thin)
    mcmcRes <- mcmc_burn_chains(mcmcRes, n = extra_burnin)
    cat("Extra thinning and extra removal from the beginning of the chains was done\n")
  } else if (!is.null(extra_thin)) {
    mcmcRes <- mcmc_thin_chains(mcmcRes, freq = extra_thin)
    cat("Extra thinning was done\n")
  } else if (!is.null(extra_burnin)) {
    mcmcRes <- mcmc_burn_chains(mcmcRes, n = extra_burnin)
    cat("Extra removal from the beginning of the chains was done\n")
  }
  
  ###########################
  ## check the convergenge ##
  ###########################
  # Get number of chains
  nChains <- length(mcmcRes)
  selectedChains <- 1:nChains #keep all chains as default
  
  ######################################
  ## plot and check the chains course ##
  ######################################
  #####################################################
  ## outliers analysis and mean component allocation ##
  #####################################################
  ## Convergence diagnostic to see if we need to discard any
  out <- mcmc_get_outliers(mcmcRes) #get data for gelman diagnostics, outliers
  meanoutProb <- mcmc_get_meanoutliersProb(mcmcRes) #outliers probability
  meanAlloc <- mcmc_get_meanComponent(mcmcRes)#mean component allocation
  
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
  ## run gelman diganostics ##
  ############################
  #only makes sense to look at this when the number of chain is 3 or higher, since than you can look at combinations of chains and kick out the bad ones
  if (nChains == 1) {
    cat("Skiping gelman diagnostics. You need more than a single chain!\n")
    
  } else if (nChains == 2) {
    cat("Running gelman diagnostics for the two chains only.\n")
    
    gelman_diag_res_T <- run_gelman(list(out, meanAlloc, meanoutProb),
                                    c("out", "meanAlloc", "meanoutProb")
                                    )
    
    #write the gelman diagnostics ouput
    write_tsv(gelman_diag_res_T,
              paste0(outDir,
                     "/gelman_diagnostics", 
                     outSuffix, 
                     ".tsv")
              )
    
    #keep the 2 chains if the statistics is good otherwise not
    if (max(gelman_diag_res_T$Upper_CI) <= 1.2) {
      cat("Both Chains are good according to the gelman statistics.\n")
      selectedChains_gelman <- selectedChains
    } else {
      cat("The combination of the 2 chains is bad according to the gelman statistics.\n")
      
      selectedChains_gelman <- vector()
    }
   
    
  } else {
    cat("Running gelman diagnostics for all chains.\n")
    
    gelman_diag_res_T <- run_gelman(list(out, meanAlloc, meanoutProb),
                                    c("out", "meanAlloc", "meanoutProb")
                                    )
    
    #write the gelman diagnostics ouput
    write_tsv(gelman_diag_res_T,
              paste0(outDir,
                     "/gelman_diagnostics", 
                     outSuffix, 
                     ".tsv")
              )
    
    #calculate the maximum of the interval
    upper_CI <- max(gelman_diag_res_T[,2]) #maximum must not be more than 1.2 (idicates problem in at least one of the dignostics)
    
    ########################
    ## select good chains ##
    ########################
    #if upper CI of gelman diagnostics is higher than 1.2 and you have more than 2 chains
    #run all combination of chains and take the one with the most chains and with the diagnostics below 1.2
    if (upper_CI <= 1.2) {
      cat("All chains converged, keeping all.\n")
      
      selectedChains_gelman <- selectedChains
    } else {
      cat("Gelman diagnostics suggest that not all chains converged.\nRunning combination of all chains and keeping the one with highest number of chain and with gelman diagnostics bellow 1.2.\n")
      
      #create all combinations
      combinations_L <- map(2:(nChains -1), #take combinations with no rep (from 2 up to number of chain minus 1)
                            ~{
                              combinations_map <- combn(1:nChains, .x) #create the combinations
                              
                              #keep them in the form of a vector in a list
                              combinations_map_L <- map(array_branch(combinations_map, 2),
                                                        ~{
                                                          .x
                                                        })
                              
                              return(combinations_map_L)
                            }) %>% 
        flatten()
      
      #do the diagnostics for all
      combinations_gelmanDia_res_T <- map(combinations_L,
                                          ~{
                                            gelman_diag_map_res_T <- run_gelman(list(out, meanAlloc, meanoutProb),
                                                                                c("out", "meanAlloc", "meanoutProb"),
                                                                                subset = .x)
                                            
                                            gelman_diag_map_res_T$nChain <- length(.x)#add the number of chains
                                            gelman_diag_map_res_T$ChainCombList <- list(.x)#add the chain combination as alist
                                            gelman_diag_map_res_T$ChainCombChar <- paste0(.x, collapse = ",")#add the chain combination as character pasted
                                            return(gelman_diag_map_res_T)
                                          }) %>% 
        purrr::reduce(rbind)
      
      #select the highest number of chains with good diagnostic
      combinations_gelmanDia_res_T <- combinations_gelmanDia_res_T %>% 
        group_by(ChainCombChar) %>% 
        mutate(max_upper_CI = max(Upper_CI)) %>%  #make the maximum for the combination of chains from all the matrices
        filter(max_upper_CI <= 1.2) %>% 
        mutate(pass = n() == 3) %>% #calculate the number of matrices left if not for now 3 then soome had the statistics worse than 1.2
        filter(pass) %>% #if so remove those that do not have the number 3 (most of them are anyways remove already with the max_upper_CI, since if 1 matric is bad the other as well it seems)
        arrange(desc(nChain), desc(Upper_CI))
      
      #select the chains
      selectedChains_gelman <- combinations_gelmanDia_res_T$ChainCombList[1][[1]]
      
      #write the results
      combinations_gelmanDia_res_T <- combinations_gelmanDia_res_T %>% 
        dplyr::select(-ChainCombList)
     
       write_tsv(combinations_gelmanDia_res_T,
                paste0(outDir,
                       "/gelman_diagnostics_comb", 
                       outSuffix, 
                       ".tsv")
                )
      
      cat(paste0("Based on Gelman diagnostics selected Chains are:"),
          paste0(selectedChains_gelman, collapse = ","),
          "\n")
      }
    }
    
  #######################################################################################################
  ## run diagnostics on each chain and if some end up below still remove them from the selected chains ##
  #######################################################################################################
  #geweke_test = porovna prvnich 10% chainu s poslednimy 50, pokud tam je rozdil asi neconverdged
  #ale muze byt ovlivneno pomoci burn in parametru pac proste convergoval pozdeji nez jme dali burn in
  FailedChains <- map(statistics_L,
                      ~{
                        #run the test
                        geweke_test_res <- geweke_test(.x)
                        
                        #write the results
                        geweke_test_res_T <- as_tibble(geweke_test_res)
                        write_tsv(geweke_test_res_T,
                                  paste0(outDir,
                                         "/geweke_test_res", 
                                         outSuffix, 
                                         ".tsv")
                                  )
                        
                        #flag bad chains
                        failTest <- geweke_test_res[2,] <= geweke_pval
                        failChains <- names(failTest[failTest])
                        
                        if (!is.null(failChains)) {
                          failChainsNumber <- parse_number(failChains) #take these that fail test
                          return(failChainsNumber)
                          
                        } else {
                          failChainsNumber <- vector()
                          return(failChainsNumber)
                        }
                        
                      }) %>% 
    unlist() %>% 
    unique()
  
  ################################################################
  ## remove the chains that failed according ti the geweke test ##
  ################################################################
  if (is_empty(FailedChains)) {
    cat("No chain was removed during the Geweke test.\n")
    
    selectedChains_geweke <- selectedChains
  } else {
    
    removedChain <- selectedChains[selectedChains == FailedChains]
    cat("Based on the Geweke test Chain ", 
        paste0(removedChain, collapse = ", "),
        " was removed.\n",
        sep = "")
    
    selectedChains_geweke <- selectedChains[selectedChains != FailedChains]
  }
  
  ##############################################################################################################
  ## if chains to be selected where provided change the selected chains otherwise use the selected statistics ##
  ##############################################################################################################
  if (typeof(MCMCres_keepChains) == "double") {
    cat("Keeping the provided Chains!\n")
    selectedChains <- MCMCres_keepChains 
  } else if(MCMCres_keepChains == "geweke_and_gelman") {
    #since gelman ti rika o kombinaci a to je dulezitejsi, tak jejich kombinace bude ze vezmes vse z gelmana co proslo gewekem
    selectedChains <-  selectedChains_gelman[selectedChains_gelman %in% selectedChains_geweke]
    cat("Keeping the chains based on the intersection of the gelman and geweke statistics!!\n")
  } else if(MCMCres_keepChains == "geweke") {
    selectedChains <- selectedChains_geweke
    cat("Keeping the chains based on the geweke statistics only!\n")
  } else if(MCMCres_keepChains == "gelman") {
    selectedChains <- selectedChains_gelman
    cat("Keeping the chains based on the gelman statistics only!\n")
  } else {
  errorCondition("Wrong condition for MCMCres_keepChains was provided, shoule be: \"geweke\", \"gelman\", \"geweke_and_gelman\" or a vector of numbers")
  }

  ################################
  ## keep the converged chains ##
  ################################
  #before running check at least 1 chain was left7
  if (is_empty(selectedChains)) {
    errorCondition("No Chain was left after chain selection. Probably a problem with your suplied chain numbers\nor that no chained passed your selected statistical test.\n")
  }
  
  mcmcRes_converged <- mcmcRes[selectedChains]
  
  ########################
  ## pooling the chains ##
  ########################
  #if covered and UNconverged samples are pooled = can llead to bad results
  mcmcRes_converged_pooled <- mcmc_pool_chains(mcmcRes_converged)
  
  save(mcmcRes_converged_pooled,
       paste0(outDir,
              "/mcmcRes_converged_pooled_for_violin_plots", 
              outSuffix, 
              ".Rdata")
       )
  
  #####################################
  ## add the predictions to the data ##
  ####################################
  mcmcRes_converged_pooled <- tagmMcmcProcess(mcmcRes_converged_pooled) 
  
  mcmcRes_converged_pooled_MSdata <- tagmPredict(object = data,
                                                 params = mcmcRes_converged_pooled,
                                                 probJoint = TRUE)
  
  ## Store allocation and outlier probabilities
  ##!!!!!!!!!!!!!!!! somehow also include shanon.entropy !!!!!!!!!!!!!!!!!!!!##
  ##!!!!!!!!!!!!!!!! somehow also include shanon.entropy !!!!!!!!!!!!!!!!!!!!##
  ##!!!!!!!!!!!!!!!! somehow also include shanon.entropy !!!!!!!!!!!!!!!!!!!!##
  ##!!!!!!!!!!!!!!!! somehow also include shanon.entropy !!!!!!!!!!!!!!!!!!!!##
  ##!!!!!!!!!!!!!!!! somehow also include shanon.entropy !!!!!!!!!!!!!!!!!!!!##
  ##!!!!!!!!!!!!!!!! somehow also include shanon.entropy !!!!!!!!!!!!!!!!!!!!##
  
  tagm_prob <- fData(mcmcRes_converged_pooled_MSdata)[, "tagm.mcmc.probability"]
  tagm_out <- 1 - fData(mcmcRes_converged_pooled_MSdata)[, "tagm.mcmc.outlier"]
  shan_out <- 1 - fData(mcmcRes_converged_pooled_MSdata)[,"tagm.mcmc.mean.shannon"]
  
  ## Create a new column containing overall probability
  fData(mcmcRes_converged_pooled_MSdata)[, "overall_prob"] <- tagm_prob * tagm_out * shan_out
  
  ## Set prediction thresholds based on overall probability
  mcmcRes_converged_pooled_MSdata <- getPredictions(mcmcRes_converged_pooled_MSdata,
                                                    fcol = "tagm.mcmc.allocation",
                                                    scol = "overall_prob",
                                                    t = p_cutoff)
  
  ###################
  ## visualization ##
  ###################
  tibble_data_T <- fData(mcmcRes_converged_pooled_MSdata)
  colnames(tibble_data_T)[1] <- "ID"
  
  ## set seed ##
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
    stop('Error: no method for dimensional reduction provided. Provided one of these: "PCA", "t-SNE" as the dimRed variable.\n')
  }
  
  ###################
  ## visualization ##
  ###################
  #get the position from the dim reduction
  pred_filtered_plot_T <- tibble(ID = rownames(pred_filtered_plot),
                                 Dim1 = pred_filtered_plot[,1],
                                 Dim2 = pred_filtered_plot[,2])

  ## add columns to the locations ##
  pred_filtered_plot_T <- left_join(pred_filtered_plot_T, tibble_data_T, by = "ID")

  ## add desc to the data ##
  pred_filtered_plot_T <- left_join(pred_filtered_plot_T, descLopit_data, by = "ID")
  
  ## create label ##
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

  ## change the svm.pred to factor ##
  pred_filtered_plot_T$overall_prob <- factor(pred_filtered_plot_T$overall_prob, levels = levels)

  ## create only marker dataset ##
  markers_plot_T <- pred_filtered_plot_T %>%
    filter(markers != "unknown")
  
  ## plot the non-interactive plot ##
  tSNE_plot <- ggplot(pred_filtered_plot_T, aes(label = label)) +
    geom_jitter(data = pred_filtered_plot_T, aes(Dim1, Dim2, color = overall_prob), size = 2, alpha = 0.5) +
    geom_point(data = markers_plot_T, aes(Dim1, Dim2), color = "black", size = 3)+
    geom_point(data = markers_plot_T, aes(Dim1, Dim2, color = overall_prob), size = 2) +
    theme_minimal(base_size = 20) +
    labs(color = "Localization")+
    scale_color_manual(values = colors)
  
  png(paste0(outDir, "/", dimRed, outSuffix, ".png"),
      width = 4000,
      height = 2500,
      res = 300
      )
  
  plot(tSNE_plot)
  
  dev.off()
  
  ## make the NA plot
  make_NA_plot(pred_filtered_plot_T = pred_filtered_plot_T,
               dimRed = dimRed,
               outSuffix = outSuffix,
               outDir = outDir
               )
  
  ## plot the data interactively in plotly ##
  plotly <- ggplotly(tSNE_plot)
  htmlwidgets::saveWidget(as_widget(plotly),
                          paste0(outDir, "/", dimRed, "_plot_interactive", outSuffix, ".html")
                          )
  
  #######################
  ## visualize shannon ##
  #######################
  ## plot the non-interactive plot ##
  tSNE_shannon_plot <- ggplot(pred_filtered_plot_T, aes(label = label)) +
    geom_jitter(data = pred_filtered_plot_T, aes(Dim1, Dim2, 
                                                 color = tagm.mcmc.allocation, 
                                                 size = tagm.mcmc.mean.shannon),
                alpha = 0.5) +
    theme_minimal(base_size = 20) +
    labs(color = "Localization", size = "Shannon entropy" ) +
    scale_color_manual(values = colors)
  
  png(paste0(outDir, "/", dimRed, "_shannon", outSuffix, ".png"),
      width = 4000,
      height = 2500,
      res = 300
      )
  
  plot(tSNE_shannon_plot)
  
  dev.off()
  
  ## plot the data interactively in plotly ##
  plotly <- ggplotly(tSNE_shannon_plot)
  htmlwidgets::saveWidget(as_widget(plotly),
                          paste0(outDir, "/", dimRed, "_shannon_plot_interactive", outSuffix, ".html")
                          )
  
  ###########################
  # write the final results #
  ###########################
  #get the probability for each compartment (asi mean nebo median)
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
  returnList <- list(data  = pred_filtered_plot_T)
  
  return(returnList)
  
  cat("Everything is done.\n")

}

run_MCMC_bandle <- function(data,
                            numIter, #the number of times the programme samples from posterior
                            burnin, #the number of steps from the begining toi discard (it takes time until it gets to the posterior from the renadom defualt, so you have to discard the begining)
                            thin, #remove autocorelation (if e.g. 20 then every 20th step is kept), it is to remove autocorrelation since steps right after one another correlate)
                            numChains, #number of chains to run
                            MCMCres_path, #put the path to the result of the MCMC analysis and it will continue with visualization and statistics,
                            MCMCres_keepChains, #a vector of chains (e.g. c(1,2,3)) to keep and skip the statistical calculation of the chain convergance             
                            extra_thin, #a number that does extra thining of the sample from the posterio, reduces auto correlation
                            extra_burnin, #a number of cycles to be removed from the beging of the chains, extra to the burnin already provided
                            geweke_pval, #testing whether in a single hcain the first 10 percent is the same as the last 50, assesses converdgence, if significant the chains did not converdge and the chain is removed
                            seed,
                            dimRed,#what dimensionality reduction algorithm to use "PCA" or "t-SNE"
                            perplexity, #, t-SNE parameter, which can be though of as the number of effective neighbors of every point
                            dims,
                            p_cutoff,
                            outDir, #full path to the output dir "my/output/dir"
                            outSuffix, #suffix of the output files
                            descLopit_data, #passed as a dataframe from the analyseLopit function,
                            hyppar, #for BUNDLE priors for the hyper paramteres
                            dirPrior, #for BUNDLE Dirichlet prior matrix, how much you expect that 1 protein will have different localization between datsets, not much
                            conditions
                            ) {
 
  ############################################
  ## prepare the prior for hyper parameters ##
  ############################################
  ## Extract subcellular classes
  mrkCl <- getMarkerClasses(data[[1]], fcol = "markers")
 
   ## Construct a pc_prior
  set.seed(seed)
  
  K <- length(mrkCl) # K = number of subcellular classes
  pc_prior <- matrix(NA, ncol = 3, K) # Initiate empty matrix
  pc_prior[seq.int(1:K), ] <- matrix(rep(hyppar, each = K), ncol = 3)
  
  ## Fit GP priors to the data - each sample separately
  set.seed(seed)
  gpParams <- lapply(data,
                     function(x) fitGPmaternPC(x, hyppar = pc_prior))
  
  ## Using plotGPmatern to overlay the predictives for all datasets
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
  ## prepare the Dirichlet prior ##
  #################################
 
  ## Set up Dirichlet prior matrix
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
        
        ## Determine prior probability of >15 differential localisation events
        predDirPrior <- prior_pred_dir(object = data[[.x]],
                                       dirPrior = dirPrior,
                                       q = 15)
        
        
        ## Plot histogram of differential localisation probability based on priors
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
  ## run the bandle MCMC algorithm ##
  ###################################
  if (is.null(MCMCres_path)) {
    cat("Running the MCMC with BANDLE from scratch since no results of the MCMC run were provided.\n")
    # prepare mulcicore for parallel computing
    n_workers <- as.integer(Sys.getenv("PBS_NCPUS", unset = "4"))  # default 4 if PBS_NCPUS not set (it should always be, it is the number of CPUS allocated)
    
    multicoreParam <- MulticoreParam(workers = n_workers) #set number of cpus??
    register(multicoreParam)
    
    #run MCMC with tryCatch to rerun if chain fails (possible it is just random the fails)
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
        
        
        success <- TRUE  # if this runs without error, mark success
      }, error = function(e) {
        message("Error occurred in MCMC run: ", e$message)
        attempt <<- attempt + 1
        if (attempt > max_tries) {
          stop("Maximum retry attempts reached. MCMC failed.")
        } else {
          message("Retrying...")
        }
      })
    }
    
    if (success) {
      save(badleres, file = paste0(outDir, "/BANDLEres", outSuffix, ".Rdata"))
    }
    
  } else {
    cat("Results of the MCMC run were provided, loading them now.\n")
    load(MCMCres_path)
  
  }
  
  ######################
  ## chech the chains ##
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
    
  write_tsv(Gelman_res_T,  paste0(outDir, "/BANDLEres_gelman_diagnostics", outSuffix, ".tsv"))
  
  
  #############################
  ## plot the outliers plots ##
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
  ## Remove chains if provided  ##
  ################################
  if (!is_empty(MCMCres_keepChains)) {
    cat("Keeping the custome chains!\n")
    
    bandleres_converged <- bandleres[MCMCres_keepChains]
  } else {
    cat("Allchains were kept, if you want to kee specific chains, plese provided them in a vector to MCMCres_keepChains!\n")
    
    bandleres_converged <- bandleres
   
  }

  ## Process BANDLE results
  params <- bandleProcess(bandleres_converged)
 
  ## Append BANDLE results to first entry in list of each condition
  res <- bandlePredict(objectCond1 = conditions_L[[1]],
                       objectCond2 = conditions_L[[2]],
                       params = params,
                       fcol = "markers")
  
  ## get the results for both conditions
  res_cond1 <- res[[1]][[1]]
  res_cond2 <- res[[2]][[1]]
  
  res_per_cond <- c(res_cond1, res_cond2)
  
  ## threshold the prediction
  ## add filtering by shanon entropy here if wnated
  res_per_cond <- map(res_per_cond,
                      ~{
                        map_data <- .x
                        class_prob <- fData(map_data)$bandle.probability
                        out_prob <- 1 - fData(map_data)$bandle.outlier
                        
                        fData(map_data)$bandle.probability.overall <- class_prob * out_prob
                        
                        ## Threshold BANDLE localisation predictions
                        map_data <- getPredictions(map_data,
                                                   fcol = "bandle.allocation",
                                                   scol = "bandle.probability.overall",
                                                   mcol = "markers",
                                                   t = p_cutoff)
                        
                        return(map_data)
                      })
  
  ## Extract BANDLE results minus marker proteins
  map2(res_per_cond,
      c("cond1", "cond2"),
      ~{
        #get the predictions wuthout markers
        map_data <- .x
        res_no_mrk <- unknownMSnSet(res_unstim_rep1, fcol = "markers")
        
        ## Pull column containing BANDLE localisation predictions in each condition
        alloc <- map_data %>%
          fData() %>%
          pull(bandle.allocation.pred)
        
        ## Pull column containing BANDLE localisation posterior estimates
        pe <- res_no_mrk %>%
          fData() %>%
          pull(bandle.probability)
        
        ## make the plots bar plot of allocation
        png(paste0(outDir,
                   "/",
                   "BANDLEres_predictions_with_prob_cutoff",
                   outSuffix,
                   ".png"
                   )
            )
        
        barplot(alloc %>% table,
                las = 2, main = "Predicted location",
                ylab = "Number of proteins")
        
        dev.off()
        
        ## make the plots of posterior destribution
        png(paste0(outDir,
                   "/",
                   "BANDLEres_predictions_posterio_prob",
                   outSuffix,
                   ".png"
                   )
            )
        
        boxplot(pe ~ alloc, las = 2, main = "Posterior",
                ylab = "Probability")
        
        dev.off()
        
      })
  
}

#########################
  ## not working for now, melo by patrit do analyseLoopt ##
  #########################
  # #### semi-supervised method ####
  # if (all(table(pred_filtered_T$svm.pred) > 6) & semiSup == T) {
  #   
  #   cat("Performing semi supervised method\n")
  #   
  #   pheno_pred <- phenoDisco(pred_filtered,
  #                            GS = 10,
  #                            times = 100,
  #                            fcol = "svm.pred")
  #   
  #   ## save the data ##
  #   save(pheno_pred, file = paste0(outDir, "/", "pheno_pred", outSuffix, ".Rdata"))
  #   
  #   plot2D(pdres, fcol = "pd")
  #   addLegend(pdres, fcol = "pd", ncol = 2,
  #             where = "bottomright",
  #             cex = .5)
  #   
  # } else {
  #   cat("Everything is done.\n")
  # }
  # return(returnList)

##############
## dodat prikald analyzy jak jsem to delal u michala 
##############
