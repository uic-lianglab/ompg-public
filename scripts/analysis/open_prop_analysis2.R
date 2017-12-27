###################################################################
# open_prop_analysis2.R
#
# Analysis relating to open vs closed conformations of OMPG
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

###################################################################
# Data set utilities
###################################################################

source(file.path(get_script_dir(), "electrostat_merge_utilities.R"))

###################################################################
# Globals
###################################################################

# Probe radius in Angstroms for open/close analysis
# 2.75 Angstroms = K
# 1.75 Angstroms = Cl
PROBE = 2.75

# Switch to toggle overwrite behavior
SHOULD_OVERWRITE = FALSE

###################################################################
# Datasets
###################################################################

SIM_WT = load_sim(base_output_dir = get_relax_dir(),
                  sim_prefix = "wt",
                  probe = PROBE,
                  max_energy_rank = 5000,
                  should_load_elec = FALSE,
                  overwrite = SHOULD_OVERWRITE)

###################################################################
# Enery rank plots
###################################################################

# Obtain cumulative open proportion
get_cum_open_prop <- function(oc) {
  # Cumulative open count
  cum_n_open = as.numeric(cumsum(oc))
  
  # Cumulative length - same as (energy) rank
  ones = rep(1.0, length(oc))
  cum_n_all = as.numeric(cumsum(ones))
  
  # Cumulative proportion of open samples
  cum_open_prop = 100.0 * (cum_n_open / cum_n_all)
  return(list(cum_open_prop = cum_open_prop,
              cum_n_all = cum_n_all))
}

# Plot % open versus energy rank:
# Example usage:
#   - pH 5, top 5k (wildtype):
#       plot_open_prop(oc=SIM_WT$pH5$merge$oc[1:5000], sub="wt, pH=5")
#   - pH 7, top 5k (wildtype):
#       plot_open_prop(oc=SIM_WT$pH7$merge$oc[1:5000], sub="wt, pH=7")
plot_open_prop <-
  function(oc = SIM_WT$pH5$merge$oc[1:5000],
           sub = "wt, pH=5") {
    ls.res = get_cum_open_prop(oc)
    cum_n_all = ls.res$cum_n_all
    cum_open_prop = ls.res$cum_open_prop
    plot(
      x = cum_n_all,
      y = cum_open_prop,
      main = paste0("% open vs enery rank\n", sub),
      xlab = "energy rank",
      ylab = "% open",
      xlim = range(cum_n_all),
      ylim = range(cum_open_prop)
    )
    
    return(cum_open_prop)
  }

plot_open_prop_w_barrel <-
  function(oc = SIM_WT$pH5$merge$oc,
           sub = "wt, pH=5") {
    # Compute total open proportion among all barrel types
    ls.res.oc = get_cum_open_prop(oc)
    # get logical for entries belonging to 2iwv (pH7, open barrel)
    logi.2iwv = grepl("2iwv", names(oc))
    # Obtain cumulative proportions of 2iwv and 2iww barrels
    ls.res.2iwv = get_cum_open_prop(logi.2iwv)
    ls.res.2iww = list(
      cum_n_all = ls.res.oc$cum_n_all,
      cum_open_prop = 100.0 - ls.res.2iwv$cum_open_prop
    )
    # Determine open proportion attributed to barrel type
    ls.res.2iwv$cum_open_prop = ls.res.2iwv$cum_open_prop / 100.0 * ls.res.oc$cum_open_prop
    ls.res.2iww$cum_open_prop = ls.res.2iww$cum_open_prop / 100.0 * ls.res.oc$cum_open_prop
    
    # Plot 2iwv contribution
    plot(
      x = ls.res.2iwv$cum_n_all,
      y = ls.res.2iwv$cum_open_prop,
      col = "red",
      main = paste0("% open vs enery rank\n", sub),
      xlab = "energy rank",
      ylab = "% open",
      xlim = range(ls.res.oc$cum_n_all),
      ylim = range(
        ls.res.oc$cum_open_prop,
        ls.res.2iwv$cum_open_prop,
        ls.res.2iww$cum_open_prop
      )
    )
    
    # Plot 2iww contribution
    points(x = ls.res.2iww$cum_n_all,
           y = ls.res.2iww$cum_open_prop,
           col = "blue")
    
    # Plot total contribution
    lines(x = ls.res.oc$cum_n_all,
          y = ls.res.oc$cum_open_prop,
          col = "black")
    
    # Turn off clipping
    # http://stackoverflow.com/questions/3932038/plot-a-legend-outside-of-the-plotting-area-in-base-graphics
    par(xpd = TRUE)
    legend(
      x = 'topright',
      legend = c("2iwv", "2iww", "All"),
      col = c("red", "blue", "black"),
      pch = 1,
      inset = c(0, 0.3)
    )
    
    return(ls.res.oc$cum_open_prop)
  }

# Example usage:
#   - pH 5, 2iww, top 5k (wildtype):
#       plot_barrel_open_rate(oc=SIM_WT$pH5$b2iww$oc[1:5000], barrel = "2iww", sub="wt, pH=5")
#   - pH 5, 2iwv, top 5k (wildtype):
#       plot_barrel_open_rate(oc=SIM_WT$pH5$b2iwv$oc[1:5000], barrel = "2iwv", sub="wt, pH=5")
#   - pH 7, 2iww, top 5k (wildtype):
#       plot_barrel_open_rate(oc=SIM_WT$pH7$b2iww$oc[1:5000], barrel = "2iww", sub="wt, pH=7")
#   - pH 7, 2iwv, top 5k (wildtype):
#       plot_barrel_open_rate(oc=SIM_WT$pH7$b2iwv$oc[1:5000], barrel = "2iwv", sub="wt, pH=7")
plot_barrel_open_rate <- function(oc = SIM_WT$pH5$merge$oc,
                                  barrel = "2iww",
                                  sub = paste0(barrel, " wt, pH=5")) {
  # get logical for entries belonging to barrel type
  logi.barrel = grepl(barrel, names(oc))
  # extract open/close status for barrel type
  oc.barrel = oc[logi.barrel]
  plot_open_prop(oc = oc.barrel,
                 sub = sub)
}

# Plot weighted % open vs energy rank
# Example usage:
#   - pH 5, top 5k:
#       plot_weighted_open_prop(oc=SIM_WT$pH5$merge$oc[1:5000],
#                               ener=SIM_WT$pH5$merge$ener[1:5000],
#                               sub="wt, pH=5")
#   - pH 7, top 5k:
#       plot_weighted_open_prop(oc=SIM_WT$pH7$merge$oc[1:5000],
#                               ener=SIM_WT$pH7$merge$ener[1:5000],
#                               sub="wt, pH=5")
plot_weighted_open_prop <-
  function(oc = SIM_WT$pH5$merge$oc[1:5000],
           # Energy vector in kcal/mol (NAMD units)
           ener = SIM_WT$pH5$merge$ener[1:5000],
           sub = "wt, pH=5",
           # Boltzmann constant in kcal/(mol*K) (NAMD units)
           # found in NAMD source file: colvarproxy_namd.h
           k = 0.001987191,
           # Temperature in Kelvin (NAMD units)
           temp = 300) {
    # Compute log of Boltzmann factor
    log_boltzf = -ener / (k * temp)
    # Determine normalization value
    norm_log_boltzf = max(log_boltzf)
    # Set all values relative to normalization factor
    log_boltzf = log_boltzf - norm_log_boltzf
    # Compute Boltzmann factor
    boltzf = exp(log_boltzf)
    
    # Compute cumulative Boltzmann factors for open conformation
    open_boltzf = boltzf
    open_boltzf[!oc] = 0.0
    cum_open = cumsum(open_boltzf)
    
    # Compute cumulative partition constant
    cum_z = cumsum(boltzf)
    
    # Cumulative proportion of open samples
    cum_open_prop = 100.0 * (cum_open / cum_z)
    
    # Cumulative length - same as (energy) rank
    ones = rep(1.0, length(oc))
    cum_length = as.numeric(cumsum(ones))
    
    plot(
      x = cum_length,
      y = cum_open_prop,
      main = paste0("weighted % open vs enery rank\n", sub),
      xlab = "energy rank",
      ylab = "% open",
      xlim = range(cum_length),
      ylim = range(cum_open_prop)
    )
    
    return(cum_open_prop)
  }

# Adjust weight by sampling energy bias
# Example usage:
#   - pH 5, top 5k:
#       plot_density_weighted_open_prop(oc=SIM_WT$pH5$merge$oc[1:5000],
#                                       ener=SIM_WT$pH5$merge$ener[1:5000],
#                                       sub="wt, pH=5")
#   - pH 7, top 5k:
#       plot_density_weighted_open_prop(oc=SIM_WT$pH7$merge$oc[1:5000],
#                                       ener=SIM_WT$pH7$merge$ener[1:5000],
#                                       sub="wt, pH=5")
plot_density_weighted_open_prop <-
  function(oc = SIM_WT$pH5$merge$oc[1:5000],
           # Energy vector in kcal/mol (NAMD units)
           ener = SIM_WT$pH5$merge$ener[1:5000],
           sub = "wt, pH=5",
           # Boltzmann constant in kcal/(mol*K) (NAMD units)
           # found in NAMD source file: colvarproxy_namd.h
           k = 0.001987191,
           # Temperature in Kelvin (NAMD units)
           temp = 300) {
    ###############################################
    # Determine sampling distribution
    ###############################################
    
    # Compute energy density estimate (bias term)
    # Consider KernSmooth package's bkde method instead
    samp_ener_dens = density(x = ener)
    
    # Create function to interpolate density estimates
    samp_ener_pdf = approxfun(
      x = samp_ener_dens$x,
      y = samp_ener_dens$y,
      method = "linear",
      yleft = 0.0,
      yright = 0.0
    )
    
    # Evaluate pdf at each sample input
    samp_prob = samp_ener_pdf(ener)
    
    ###############################################
    # Determine target distribution
    ###############################################
    
    # Compute log of Boltzmann factor
    targ_log_boltzf = -ener / (k * temp)
    # Determine normalization value
    targ_norm_log_boltzf = max(targ_log_boltzf)
    # Set all values relative to normalization factor
    targ_log_boltzf = targ_log_boltzf - targ_norm_log_boltzf
    # Compute Boltzmann factor
    targ_boltzf = exp(targ_log_boltzf)
    
    ###############################################
    # Compute self-normalized weighted expectation
    ###############################################
    
    # Compute weight of each sample as
    #   target_probability / sample_probability
    #     = targ_boltzf / samp_prob
    w = targ_boltzf / samp_prob
    
    f_open = w
    f_open[!oc] = 0.0
    
    cum_open = cumsum(f_open)
    cum_z = cumsum(w)
    
    cum_open_prop = 100.0 * (cum_open / cum_z)
    
    ones = rep(1.0, length(oc))
    cum_length = as.numeric(cumsum(ones))
    
    plot(
      x = cum_length,
      y = cum_open_prop,
      main = paste0("weighted % open vs enery rank (corrected)\n", sub),
      xlab = "energy rank",
      ylab = "% open",
      xlim = range(cum_length),
      ylim = range(cum_open_prop)
    )
    
    return(cum_open_prop)
  }
