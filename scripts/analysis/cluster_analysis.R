###################################################################
# cluster_analysis.R
#
# Utlities for wrangling and analysing cluster data
###################################################################

###################################################################
# Common path utilities
###################################################################

# Cached directory to script
# https://stackoverflow.com/questions/3452086/getting-path-of-an-r-script
SCRIPT_DIR = getSrcDirectory(function(x) {
  x
})

# @return full path to directory containing this script
get_script_dir <- function() {
  return(SCRIPT_DIR)
}

# @return root directory of project
get_root_dir <- function() {
  file.path(get_script_dir(), '..', '..')
}

# @return base output directory of project
get_base_output_dir <- function() {
  file.path(get_root_dir(), "ompg", "output")
}

# @return base cluster directory of project
get_clust_dir <- function() {
  file.path(get_base_output_dir(), "multi_loop", "clust")
}

# @param sim_id - simulation identifier
# @return directory containing pca plotting coordinates
get_pca_plots_dir <- function(sim_id = "wt_2iwv") {
  file.path(get_clust_dir(), "g_pca_for_plots", sim_id)
}

# @param sim_id - simulation identifier
# @return path to csv containing pca plot data
get_pca_plots_csv_path <- function(sim_id = "wt_2iwv") {
  file.path(get_pca_plots_dir(sim_id),
            paste0(sim_id, ".pca.plot.csv"))
}

# @param sim_id - simulation identifier
# @return directory containing cluster results
get_clust_results_dir <- function(sim_id = "wt_2iwv") {
  file.path(get_clust_dir(), "d_clust_pdbs", sim_id)
}

#@param sim_id - simulation identifier
#@param dat_id - data identifier, must be one of {c2s|ex|s2c}
#   - c2s: cluster to sample mapping
#   - ex -
get_clust_results_csv_path <- function(sim_id = "wt_2iwv",
                                       dat_id = c("c2s", "ex", "s2c")[1]) {
  stopifnot(dat_id %in% c("c2s", "ex", "s2c"))
  file.path(get_clust_results_dir(sim_id),
            paste0(sim_id, ".pca.bin.", dat_id, ".csv"))
}

###################################################################
# Data set utilities
###################################################################

source(file.path(get_script_dir(), "load_data_utils.R"))

# Utility loads data sets necessary for clustering
# @param sim_id - simulation identifier
# @param overwrite - If TRUE, overwrites all cached data
# @return list with fields:
#   $mat - pca matrix where each sample is a column
#   $s2c - sample to 1-based cluster id mapping
#   $ex - 1-based indices of exemplars (gives column in mat)
#   $plot_base_path - Path to output plot data sans extension
load_clust_info <- function(sim_id = "wt_2iwv", overwrite = FALSE) {
  # Determine base plot path
  pca_plots_csv_path = get_pca_plots_csv_path(sim_id)
  plot_base_path = gsub("\\.csv$", "", pca_plots_csv_path)
  
  # Load data sets
  
  mat = load_csv_cached(
    csv_path = pca_plots_csv_path,
    header = FALSE,
    to_matrix = TRUE,
    overwrite = overwrite
  )
  
  s2c = as.integer(
    load_csv_cached(
      csv_path = get_clust_results_csv_path(sim_id = sim_id,
                                            dat_id = "s2c"),
      header = FALSE,
      to_matrix = FALSE,
      overwrite = overwrite
    )[, 1]
  )
  
  ex = as.integer(
    load_csv_cached(
      csv_path = get_clust_results_csv_path(sim_id = sim_id,
                                            dat_id = "ex"),
      header = FALSE,
      to_matrix = FALSE,
      overwrite = overwrite
    )[, 1]
  )
  
  # Convert from 0-based to 1-based indexing
  s2c = s2c + 1
  ex = ex + 1
  
  return(list(
    mat = mat,
    s2c = s2c,
    ex = ex,
    plot_base_path = plot_base_path
  ))
}

###################################################################
# Globals
###################################################################

# 2iwv cluster info
WT_2IWV_CLUST_INFO = load_clust_info(sim_id = "wt_2iwv")

# 2iww cluster info
WT_2IWW_CLUST_INFO = load_clust_info(sim_id = "wt_2iww")

###################################################################
# Plots
###################################################################

# Default scatter plot title which reports number of samples and PDB
# (note: does not report mutant)
# @param sim_id - simulation identifier
# @param sim_info - simulation info list with members $mat, $s2c
# @return Scatter plot title derived from simulation identifier and data
get_cpp_clusters_plot_title <-
  function(sim_id = "wt_2iwv",
           sim_info = load_clust_info(sim_id)) {
    pdb_id = strsplit(sim_id, "_")[[1]][2]
    title_ = paste0(
      "Clustering of ",
      ncol(sim_info$mat),
      " OmpG multi-loop samples\n(",
      length(unique(sim_info$s2c)),
      " clusters, PDB template: ",
      pdb_id,
      ")"
    )
    return(title_)
  }

# Plots first two principle components for each sample and color codes by
# cluster assignment.
# @param mat_pca - each column is a sample in principal component space
# @param cpp_clusters - vector mapping sample to 1-based cluster id
# @param title_ - title of scatter plot
# @param ex_indices - [optional] vector containing sample ids (1-based indices) of exemplars
#   may be set to NULL
plot_cpp_clusters <- function(mat_pca = WT_2IWV_CLUST_INFO$mat,
                              cpp_clusters = WT_2IWV_CLUST_INFO$s2c,
                              title_ = get_cpp_clusters_plot_title(),
                              ex_indices = WT_2IWV_CLUST_INFO$ex) {
  # Plot first two principal components of dataset
  x = mat_pca[1,] # PC1
  y = mat_pca[2,] # PC2
  nclust = length(unique(cpp_clusters))
  clust_colors = rainbow(n = nclust)
  # Plot first two principle components as 2-d scatter plot
  plot(
    x,
    y,
    col = clust_colors[cpp_clusters],
    xlab = "PC1",
    xlim = range(x),
    ylab = "PC2",
    ylim = range(y),
    main = title_,
    pch = 19,
    cex = 0.8
  )
  
  if (!is.null(ex_indices) && (length(ex_indices) > 0)) {
    # Connect samples to their exemplars
    segments(
      x0 = x,
      y0 = y,
      x1 = x[ex_indices[cpp_clusters]],
      y1 = y[ex_indices[cpp_clusters]],
      col = clust_colors[cpp_clusters]
    )
    
    # Plot exemplars as hollow black squares
    points(
      x = x[ex_indices],
      y = y[ex_indices],
      col = "black",
      type = "p",
      pch = 22,
      cex = 1.5
    )
  }
}

###################################################################
# Shell
###################################################################

# Generates plots for each simulation identifier
do_plot_cpp_clusters <-
  function(sim_ids = c("wt_2iwv", "wt_2iww")) {
    # Utility to export plot of specific file type
    do_plot <- function(ext, nfo, title_, ...) {
      plot_path = paste0(nfo$plot_base_path, ext)
      
      print(paste0("Generating plot: ", plot_path))
      
      plot_fn = function() {
        plot_cpp_clusters(
          mat_pca = nfo$mat,
          cpp_clusters = nfo$s2c,
          title_ = title_,
          ex_indices = nfo$ex
        )
      }
      
      # Special case: powerpoint
      if (ext == ".pptx") {
        # Increase java heap size
        # http://www.bramschoenmakers.nl/en/node/726 
        options(java.parameters = "-Xmx32g")
        library(ReporteRs)
        doc = pptx()
        doc = addSlide(doc, slide.layout = "Title and Content")
        doc = addPlot(doc,
                      fun = plot_fn,
                      vector.graphic = TRUE)
        doc = addTitle(doc, "Cluster Scatter Plot")
        writeDoc(doc, file = plot_path)
        return(TRUE)
      }
      
      if (ext == ".pdf")
      {
        pdf(file = plot_path, ...)
      } else if (ext == ".png") {
        png(filename = plot_path,
            width = 1080,
            height = 810,
            ...)
      } else {
        print(paste0("Unrecognized graphics device extension: ", ext))
        return(FALSE)
      }
      plot_fn()
      dev.off()
      return(TRUE)
    }
    
    for (sim_id in sim_ids) {
      nfo = load_clust_info(sim_id = sim_id)
      title_ = get_cpp_clusters_plot_title(sim_id = sim_id, sim_info = nfo)
      
      do_plot(ext = ".pptx",
              nfo = nfo,
              title_ = title_)
      do_plot(ext = ".pdf",
              nfo = nfo,
              title_ = title_)
      do_plot(ext = ".png",
              nfo = nfo,
              title_ = title_)
    }
  }
