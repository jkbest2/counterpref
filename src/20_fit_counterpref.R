## if script is run in the REPL interactively, use the local `spatq` package and
## manually set a range of replicates to fit, because they were probably not
## passed as command line arguments when R was started. Otherwise (typically on
## Hyak via SLURM), use the installed version of `spatq` and read in the
## replicate numbers from the command line arguments.
if (interactive()) {
  devtools::load_all("~/dev/spatq", helpers = FALSE)
  repl_arg <- c(1, 5)
} else {
  library(spatq)
  repl_arg <- as.numeric(commandArgs(trailingOnly = TRUE))
}
library(tidyverse)

## Which simulation study are we fitting?
study <- "counterpref"
## What range of replicates are going to be fit?
repls <- repl_arg[1]:repl_arg[2]
## How many years to fit?
max_T <- 15
## Tune the optimization routine
optcontrol <- list(eval.max = 1000L, iter.max = 750L)
## Names of the operating models
opmods <- 1:5
## Names of the estimation models
estmods <- c("comm500x", "comm1000x", "comm2000x",
             "comm500", "comm1000", "comm2000")

## List all possible combinations of OM/EM in given replicate range
specify_fits <- function(study, repls, opmods, estmods, root_dir = "sims") {
  create_res_dir(study, repls, root_dir = root_dir)

  res_paths <- all_res_file_paths(study, repls, opmods, estmods, root_dir)

  df <- cross_df(list(estmod = estmods, opmod = opmods, repl = repls)) %>%
    mutate(Rdata = res_paths$rdata,
           index = res_paths$indexcsv,
           sub_df = map(estmod, specify_subset),
           estd = map(estmod, estmod_pars))
  df
}

## Which still need to be fit?
fits_todo <- function(fit_spec, result_root = "sims") {
  fit_spec %>%
    filter(!file.exists(Rdata))
}

## How many observations to use from each ?
specify_subset <- function(estmod) {
  sub_df <- switch(estmod,
                   comm500x = data.frame(vessel_idx = c(1, 2), n = c(0, 500)),
                   comm1000x = data.frame(vessel_idx = c(1, 2), n = c(0, 1000)),
                   comm2000x = data.frame(vessel_idx = c(1, 2), n = c(0, 2000)),
                   comm500 = data.frame(vessel_idx = 2, n = 500),
                   comm1000 = data.frame(vessel_idx = 2, n = 1000),
                   comm2000 = data.frame(vessel_idx = 2, n = 2000))
  sub_df
}

## Specify which parameters to estimate for each estimation model; don't
## estimate catchability parameters if using a single survey vessel.
estmod_pars <- function(estmod) {
  specify_estimated(beta = TRUE,
                    gamma = FALSE,
                    omega = list(omega_n = TRUE,
                                 omega_w = FALSE),
                    epsilon = FALSE,
                    phi = FALSE,
                    psi = FALSE,
                    kappa_map =
                      c(1, NA, NA, NA, NA, NA, NA, NA),
                    obs_lik = 1L)
}

fit_list <- fits_todo(specify_fits(study = study,
                                   repls = repls,
                                   opmods = opmods,
                                   estmods = estmods,
                                   root_dir = "sims"))

## Iterate over rows to fit each model
for (idx in seq_len(nrow(fit_list))) {
  spec <- fit_list[idx, ]

  setup <- spatq_simsetup(repl = spec$repl,
                          study,
                          spec$opmod,
                          spec$sub_df[[1]],
                          max_T = max_T,
                          root_dir = "sims",
                          index_step = 1,
                          spec_estd = spec$estd[[1]])
  obj <- spatq_obj(setup,
                   runSymbolicAnalysis = TRUE,
                   normalize = TRUE,
                   silent = TRUE)

  fit <- tryCatch({
    ## Fit with large number of iterations and do it twice so more likely to
    ## reach optimum. Previous fits have ended early and/or with large
    ## gradient components, and many of these did not have PD Hessians
    fit <- spatq_fit(obj = obj, control = optcontrol)
    fit <- spatq_fit(obj = obj, fit = fit, control = optcontrol)
    fit},
    error = function(e) list(fail = TRUE))
  lpb <- tryCatch(
    gather_nvec(obj$env$last.par.best),
    error = function(e) list(fail = TRUE))
  rep <- tryCatch(
    report_spatq(obj),
    error = function(e) list(fail = TRUE))
  sdr <- tryCatch(
    sdreport_spatq(obj),
    error = function(e) list(fail = TRUE))

  saveRDS(list(spec = spec, fit = fit, lpb = lpb, rep = rep, sdr = sdr),
          spec$Rdata)

  ## Read true population state and calculate index
  true_index <- read_popstate(study = study,
                              repl = spec$repl,
                              opmod = spec$opmod,
                              root_dir = "sims") %>%
    rename(year = time,
           raw_true = pop) %>%
    filter(year <= max_T) %>%
    mutate(index_true = rescale_index(raw_true)$index)

  if (!("fail" %in% names(sdr))) {
    ## Organize details for estimated index
    which_index <- which(names(sdr$value) == "Index")
    est_index <- tibble(repl = spec$repl,
                        opmod = spec$opmod,
                        estmod = spec$estmod,
                        year = 1:max_T,
                        raw_est = sdr$value[which_index],
                        index_est = rescale_index(raw_est)$index,
                        raw_unb = sdr$unbiased$value[which_index],
                        index_unb = rescale_index(raw_unb)$index,
                        raw_sd = sdr$sd[which_index],
                        index_sd = rescale_index(raw_est, raw_sd)$sd,
                        raw_unb_sd = sdr$unbiased$sd,
                        unb_sd = rescale_index(raw_unb, raw_unb_sd)$sd)
  } else {
    est_index <- tibble(repl = spec$repl,
                        opmod = spec$opmod,
                        estmod = spec$estmod,
                        year = 1:max_T,
                        raw_est = rep(NA, max_T),
                        index_est = rep(NA, max_T),
                        raw_unb = rep(NA, max_T),
                        index_unb = rep(NA, max_T),
                        raw_sd = rep(NA, max_T),
                        index_sd = rep(NA, max_T),
                        raw_unb_sd = rep(NA, max_T),
                        unb_sd = rep(NA, max_T))

  }
  ## Join and write to CSV file
  index_df <- left_join(est_index, true_index, by = "year")
  write_csv(index_df, spec$index)
}
