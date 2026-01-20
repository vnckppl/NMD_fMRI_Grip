#!/usr/bin/env Rscript
## * Environment
base <- "/datadisk/tmp/20251119_NMD_ClenchTask"
idir <- file.path(base, "data")
odir <- file.path(base, "20251123_1stlevel")
dir.create(odir, showWarnings = FALSE, recursive = TRUE)


## * Libraries
library(INLA)
library(ciftiTools)
ciftiTools.setOption(
  "wb_path", "/home/vincent/Software/workbench/2.1.0/bin_linux64"
)
library(BayesfMRI)
library(hrf)
library(ggplot2)


## * Function for design plot SPM-style with rotated x-axis labels
plot_spm_design <- function(spikes, nuisances, filename) {

  ## ** Create single data matrix
  x_mat <- as.matrix(cbind(nuisances, spikes))
  nr <- nrow(x_mat)
  nc <- ncol(x_mat)

  ## ** Create long-format data frame
  df <- data.frame(
    x = rep(colnames(x_mat), each = nr),  # column names repeated for all rows
    y = rep(nr:1, times = nc),            # reverse rows so top = first row
    value = as.vector(x_mat)
  )

  ## ** Make x a factor to preserve column order
  df$x <- factor(df$x, levels = colnames(x_mat))

  ## ** Plot
  p <- ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "black") +   # grayscale
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      axis.ticks.length = unit(0.2, "cm"),
      axis.ticks.x = element_line(color = "black")
    )

  ## ** Save to file
  ggsave(filename, p, width = 10, height = 10)

  # return(p)
}


## * Function for exporting results to cifti, gifti, and nifti
replace_na_with_label <-
  function(xii, NA_value = -1, NA_label = "NA", NA_color = "#FFFFFFFF") {
    stopifnot(is.xifti(xii))

    # Make new row for NA label
    col_table <- col2rgb(NA_color, alpha = TRUE) / 255
    col_table <- rbind(NA_value, col_table)
    rownames(col_table) <- c("Key", "Red", "Green", "Blue", "Alpha")
    col_table <- as.data.frame(t(col_table))
    rownames(col_table) <- NA_label

    # Add this row to each label table
    nc <- ncol(xii)
    for (cc in seq(nc)) {
      stopifnot(!(NA_value %in% c(xii$meta$cifti$labels[[cc]]$Key)))
      xii$meta$cifti$labels[[cc]] <- rbind(
        col_table, xii$meta$cifti$labels[[cc]]
      )
    }

    # Replace NA values with the new value
    q <- as.matrix(xii)
    q[is.na(q)] <- NA_value
    xii <- newdata_xifti(xii, q)

    stopifnot(is.xifti(xii))
    xii
  }


## * Setup design matrix
## ** Timings and duration for hand clenching
clench <- data.frame(
  start = c(10, 55, 100, 145, 190, 235),
  duration = rep(30, 6),
  idk = rep(1, 6)
)

## ** Timings and duration for rest (just looking at crosshair)
crosshair <- data.frame(
  start = c(0, 40, 85, 130, 175, 220, 265),
  duration = rep(15, 7),
  idk = rep(1, 7)
)

## ** Combine into list
events <- list(clench = clench, crosshair = crosshair)


## * List all input files
files <- Sys.glob(file.path(idir, "*.nii"))
subs <- gsub(paste0(idir, "/"), "", files) |> substr(1, 10)


## * Loop over subjects
for (sub in subs) {
#for (sub in "sub-NMD001") {

  ## ** Announce
  vkr::h2(paste("Working on:", sub))

  ## ** Environment
  ## *** Folder for design related files
  odir_d <- file.path(odir, sub, "01_Design")
  dir.create(odir_d, showWarnings = FALSE, recursive = TRUE)
  ## *** Folder for BayesfMRI first level output
  odir_b <- file.path(odir, sub, "02_BayesfMRI")
  dir.create(odir_b, showWarnings = FALSE, recursive = TRUE)
  ## *** Folder for Significant Activation output
  odir_a <- file.path(odir, sub, "03_Activation")
  dir.create(odir_a, showWarnings = FALSE, recursive = TRUE)


  ## ** Prepare Regressors
  regressors <- read.csv(
    file.path(
      idir, paste0(sub, "_ses-01_task-handdom_desc-confounds_timeseries.tsv")),
    header = TRUE, sep = "	"
  )

  ## *** Nuisance Regressors
  ## **** Select nuisance regressors
  reg_cols <- c(
    "trans_x", "trans_y", "trans_z",
    "rot_x", "rot_y", "rot_z",
    "t_comp_cor_00", "t_comp_cor_01", "t_comp_cor_02",
    "c_comp_cor_00", "c_comp_cor_01", "c_comp_cor_02",
    "w_comp_cor_00", "w_comp_cor_01", "w_comp_cor_02"
  )

  ## **** Subset regressor data frame
  nuisances <- regressors[, reg_cols]

  ## **** Save output file
  write.table(
    nuisances, file = file.path(odir_d, "nuisance_regressors.txt"), sep = " ",
    col.names = FALSE, row.names = FALSE, quote = FALSE
  )

  ## *** Spike Regressors
  ## **** Select spike regressors
  spike_cols <- c(
    names(regressors)[grepl("^non_steady_state_outlier", names(regressors))],
    names(regressors)[grepl("^motion", names(regressors))]
  )

  ## **** Subset regressor data frame
  spikes <- regressors[, spike_cols]

  ## **** Remove duplicates
  # In some cases, the steady-state outliers and the motion outliers may be the
  # same. Make sure that no duplicate regressors are in the matrix, because that
  # will make it impossible to invert.
  spikes <- t(unique(t(spikes)))

  ## **** Save output file
  write.table(
    spikes, file = file.path(odir_d, "spike_regressors.txt"), sep = " ",
    col.names = FALSE, row.names = FALSE, quote = FALSE
  )


  ## ** Load data
  ## *** Define file paths
  fname_bold <- file.path(
    idir, paste0(
            sub, "_ses-01_task-handdom_space-fsLR_den-91k_bold.dtseries.nii"
          ))

  ## *** Read in data
  # Also resmple the data to 10K vertices
  (bold <- read_cifti(
     cifti_fname = fname_bold,
     brainstructures = "all",
     resamp_res = 10000))
  (reptime <- bold$meta$cifti$time_step)
  (vols <- ncol(bold))


  ## ** Build design
  # We can optionally set dHRF = 1 to include the temporal derivative of the
  # HRF. This allows for small shifts in the timing of the HRF. Mandy suggested
  # to set this to 0 for our design.

  # This function checks for collinearity between the predictors by returning
  # and printing the maximum pairwise correlation and the variance inflation
  # factor (VIF). If the VIF exceeds five, it is important to examine the design
  # matrix for sources of multicollinearity.

  # In addition to the design matrix itself, this function returns additonal
  # information, including the hemodynamic response function (HRF) convolved
  # with the stimulus response function.
  design <- make_design(EVs = events, nTime = vols, TR = reptime, dHRF = 0)

  ## *** Plot design and save to output file
  filename <- file.path(odir_d, "design_hrf.pdf")
  pdf(filename, width = 10, height = 10)
  plot(design, colors = c("red", "black", "pink", "grey"))
  (dev.off())

  ## *** Check the design list
  names(design)
  head(design$design)

  ## *** Save the design file (spm design matrix style image)
  filename <- file.path(odir_d, "design_spm.pdf")
  plot_spm_design(
    spikes = spikes,
    nuisances = nuisances,
    filename = filename
  )

  ## *** Save the design file (line plot)
  filename <- file.path(odir_d, "design_line.pdf")
  pdf(filename, width = 10, height = 10)
  plotdata <- as.data.frame(cbind(nuisances, spikes))
  names(plotdata) <- paste0("V", seq_len(ncol(plotdata)))
  plot_design(plotdata, method = "lineplot")
  (dev.off())


  ## ** Build contrast
  # We are not interested in average activation during a condition, but in
  # in the contarst 'squeeze' vs. 'crosshair'. Create this design here.
  ## *** Backup the old design
  design_main_effects <- design

  ## *** Create contrasts
  # We are going to contrast the clenc and crosshair columns and store these in
  # the design object. The code calculates a model for both positive and
  # negative values, so no need to calculate two contrasts.
  contrasts <- c("CL_gt_CH")
  design$design[, 1] <-
    (design_main_effects$design[, 1] - design_main_effects$design[, 2]) * 0.5
  # Remove the second column (or actually, just keep only the first)
  design$design <- design$design[, 1, drop = FALSE]
  colnames(design$design) <- contrasts
  design$field_names <- contrasts
  head(design$design)

  ## *** Plot design and save to output file
  filename <- file.path(odir_d, "design_hrf_contrast.pdf")
  pdf(filename, width = 10, height = 10)
  plot(design, colors = c("red", "black", "pink", "grey"))
  (dev.off())


  ## ** List of subcortical ROIs
  subcortical_rois <- c(
    "Accumbens-L", "Accumbens-R", "Amygdala-L", "Amygdala-R", "Brain Stem",
    "Caudate-L", "Caudate-R", "Cerebellum-L", "Cerebellum-R",
    "Diencephalon-L", "Diencephalon-R", "Hippocampus-L", "Hippocampus-R",
    "Pallidum-L", "Pallidum-R", "Putamen-L", "Putamen-R",
    "Thalamus-L", "Thalamus-R"
  )


  ## ** Fit the model with BayesfMRI
  # This took 10 minutes with 10 cores.
  system.time(
    bglm <- BayesGLM(
      BOLD = bold,
      design = design$design,
      brainstructures = c("sub", "left", "right"),
      subROI = subcortical_rois,
      TR = reptime,
      nuisance = nuisances,
      scrub = grepl(1, rowSums(spikes)),
      scale_BOLD = "mean",
      surfL = "fs_LR",
      surfR = "fs_LR",
      hpf = .01,
      nbhd_order = 1,
      ar_order = 3,
      ar_smooth = 0,
      Bayes = TRUE,
      verbose = 1,
      meanTol = 1,
      n_threads = 10
    ))


  ## ** Save result object
  save(bglm, file = file.path(odir_b, "bglm.Rdata"))


  ## ** Show results
  ## *** Bayesian
  for (c in contrasts) {
    ofile <- file.path(odir_b, paste0("estimates_", c, ".html"))
    plot(bglm, idx = c, fname = ofile, zlim = c(-2, 2))
  }

  ## *** Classical
  for (c in contrasts) {
    ofile <- file.path(odir_b, paste0("estimates_glassical_glm_", c, ".html"))
    plot(bglm, idx = c, fname = ofile, Bayes = FALSE, zlim = c(-2, 2))
  }


  ## ** Export results to cifti, gifti, and nifti files
  # Potentially, I could aslo use write_xifti2() for this.
  ## *** Export cifti file
  ofile_cii <- file.path(odir_b, paste0(sub, "_bglm_1st_level.dscalar.nii"))
  write_cifti(
    bglm$estimate_xii$Bayes$single_sess,
    cifti_fname = ofile_cii,
    verbose = TRUE
  )

  ## *** Split cifti in gifti and nifti file
  separate_cifti(
    cifti_fname = ofile_cii,
    brainstructures = "all",
    write_dir = odir_b
  )

  ## *** Resample gifti surfaces used by CiftiTools
  # Because the fMRI data was resampled to 10k, I need to do the same
  # for the surface models. This is quick and the files are small, so I
  # will just do this for each subject.
  ciftitools_dir <- file.path(.libPaths()[1], "ciftiTools/extdata")

  ## **** Left Hemisphere
  resample_gifti(
    original_fname = file.path(
      ciftitools_dir, "S1200.L.inflated_MSMAll.32k_fs_LR.surf.gii"),
    target_fname = file.path(
      odir_b, "S1200.L.inflated_MSMAll.10k_fs_LR.surf.gii"),
    hemisphere = "left",
    resamp_res = 10000)

  ## **** Right Hemisphere
  resample_gifti(
    original_fname = file.path(
      ciftitools_dir, "S1200.R.inflated_MSMAll.32k_fs_LR.surf.gii"),
    target_fname = file.path(
      odir_b, "S1200.R.inflated_MSMAll.10k_fs_LR.surf.gii"),
    hemisphere = "right",
    resamp_res = 10000)


  ## ** Caculate significance maps
  act <- activations(
    bglm,
    Bayes = TRUE,
    gamma = 0,
    alpha = 0.05,
    verbose = 1
  )
  print(act)


  ## ** Save result object with significant activations
  save(act, file = file.path(odir_a, "activation.Rdata"))


  ## ** Plot significant activations
  for (c in contrasts) {
    ofile <- file.path(odir_a, paste0("activation_", c, ".html"))
    plot(act, fname = ofile, idx = c, title = "spatial GLM")
  }

  ## *** Export cifti file
  ofile_cii <- file.path(odir_a, paste0(sub, "_bglm_1st_level.dlabel.nii"))
  write_cifti(
    replace_na_with_label(act$activations_xii$single_sess),
    cifti_fname = ofile_cii,
    verbose = TRUE
  )

  ## *** Split cifti in gifti and nifti file
  separate_cifti(
    cifti_fname = ofile_cii,
    brainstructures = "all",
    write_dir = odir_a
  )

}
