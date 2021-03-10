## if script is run in the REPL interactively, use the local `spatq` package and
## manually set a range of replicates to fit, because they were probably not
## passed as command line arguments when R was started. Otherwise (typically on
## Hyak via SLURM), use the installed version of `spatq` and read in the
## replicate numbers from the command line arguments.
if (interactive()) {
  devtools::load_all("~/dev/spatq", helpers = FALSE)
  repl_arg <- c(1, 25)
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
## Names of the operating models
opmods <- 1:5
## Names of the estimation models
estmods <- c("comm500",
             "comm1000",
             "comm2000")

results <- lapply(study,
                  all_res_file_paths,
                  repls = repls,
                  opmods = opmods,
                  estmods = estmods,
                  root_dir = "sims")
names(results) <- study

get_convcode <- function(fitlist) {
  pluck(fitlist, "fit", "convergence")
}

get_mgc <- function(fitlist) {
  pluck(fitlist, "fit", attr_getter("mgc"))
}

read_index_csv <- function(filename) {
  read_csv(filename,
           col_types = cols(
             ## repl = col_factor(levels = repls),
             repl = col_integer(),
             ## opmod = col_factor(levels = opmods),
             opmod = col_integer(),
             estmod = col_factor(levels = estmods),
             year = col_integer(),
             raw_est = col_double(),
             index_est = col_double(),
             raw_unb = col_double(),
             index_unb = col_double(),
             raw_sd = col_double(),
             index_sd = col_double(),
             raw_unb_sd = col_double(),
             unb_sd = col_double(),
             raw_true = col_double(),
             index_true = col_double())) %>%
    mutate(repl = factor(repl, levels = repls),
           opmod = factor(opmod, levels = opmods))
}

read_all_indices <- function(csv_list) {
  map_df(csv_list, read_index_csv)
}

evaluate_bias <- function(index_df) {
  index_df %>%
    group_by(opmod, estmod) %>%
    nest() %>%
    mutate(mod = map(data, ~ lm(log(raw_unb) ~ repl + log(raw_true), data = .x)),
           coef = map(mod, coef),
           delta = map_dbl(coef, pluck, "log(raw_true)")) %>%
    select(opmod, estmod, delta)
}

evaluate_rmse <- function(index_df) {
  index_df %>%
    mutate(sq_err = (index_unb - index_true)^2) %>%
    group_by(opmod, estmod) %>%
    summarize(rmse = sqrt(mean(sq_err))) %>%
    select(opmod, estmod, rmse)
}

evaluate_calibration <- function(index_df) {
  index_df %>%
    mutate(pnorm = pnorm(index_true, index_unb, unb_sd),
           n_survey = c(50, 100, 200, 300, 400)[opmod],
           n_comm = as.numeric(substr(estmod, 5, 10))) %>%
    ggplot(aes(x = pnorm, y = stat(density), fill = estmod)) +
    geom_histogram(breaks = seq(0.0, 1.0, 0.1)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    facet_grid(n_survey ~ n_comm) +
    labs(x = "Quantile") +
    guides(fill = FALSE) +
    theme_minimal() +
    theme(axis.line.x = element_line(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank())
}

plot_index_devs <- function(index_df) {
  index_df %>%
    mutate(dev = index_est - index_true) %>%
    ggplot(aes(x = year, y = dev, color = estmod, group = repl)) +
    geom_line(alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = 2) +
    facet_grid(estmod ~ opmod) +
    guides(color = FALSE)
}

plot_bias <- function(index_df) {
  index_df %>%
    ggplot(aes(x = index_true, y = index_est, color = estmod, group = repl)) +
    geom_point(alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    facet_grid(estmod ~ opmod) +
    coord_fixed() +
    guides(color = FALSE)
}

eval_dir <- "evaluation"
pdonly <- FALSE

## Allow postprocessing fits when fitting is still in progress
res_df <- tibble(rdata = results[[study]]$rdata,
                 index_csv = results[[study]]$indexcsv,
                 rdata_exists = map_lgl(rdata, file.exists),
                 index_exists = map_lgl(index_csv, file.exists)) %>%
  filter(rdata_exists, index_exists)

index_df <- map_df(res_df$index_csv, read_index_csv)

## index_df <- res_df %>%
##   map_read_all_indices(results[[study]]$indexcsv) %>%
##   filter(complete.cases(.))

convg_df <- map_dfr(res_df$rdata,
                    ~ readRDS(.)$spec) %>%
  select(estmod, opmod, repl, Rdata) %>%
  mutate(repl = factor(repl, levels = 1:50),
          convcode = map_dbl(res_df$rdata,
                            ~ readRDS(.)$fit$convergence),
          outer_mgc = map_dbl(res_df$rdata,
                              ~ max(readRDS(.)$fit$grad)),
          pdhess = map_lgl(res_df$rdata,
                          function(fn) {
            pdhess <- readRDS(fn)$sdr$pdHess
            if(is.null(pdhess)) pdhess <- FALSE
            return(pdhess)
          })) %>%
  select(estmod, opmod, repl, pdhess) %>%
  mutate(opmod = factor(opmod, levels = opmods),
          estmod = factor(estmod, levels = estmods))

index_df <- left_join(index_df, convg_df,
                      by = c("estmod", "opmod", "repl"))

pdhess_df <- index_df %>%
  group_by(estmod, opmod, repl) %>%
  summarize(pdhess = all(pdhess)) %>%
  summarize(pdhess = sum(pdhess)) %>%
  pivot_wider(names_from = opmod, values_from = pdhess)

if (pdonly) {
    index_df <- filter(index_df, pdhess)
}

if (!file.exists(eval_dir))
  dir.create(eval_dir)

if (!file.exists(file.path(eval_dir, study)))
  dir.create(file.path(eval_dir, study))

bias_df <- evaluate_bias(index_df)
bias_wide <- bias_df %>%
  pivot_wider(names_from = opmod,
              values_from = delta)
rmse_df <- evaluate_rmse(index_df)
rmse_wide <- rmse_df %>%
  pivot_wider(names_from = opmod,
              values_from = rmse)

bias_plot <- plot_bias(index_df)
calibration_plot <- evaluate_calibration(index_df)
index_devs <- plot_index_devs(index_df)

write_csv(bias_df, "figs/bias.csv")
write_csv(bias_wide, "figs/bias_wide.csv")
write_csv(rmse_df, "figs/rmse.csv")
write_csv(rmse_wide, "figs/rmse_wide.csv")
ggsave("figs/bias_plot.svg", bias_plot, scale = 1/2)
ggsave("figs/calibration.svg", calibration_plot, scale = 1/2)
ggsave("figs/index_devs.svg", index_devs, scale = 1/2)
ggsave("figs/bias_plot.png", width = 7, height = 7, bias_plot, scale = 1/2)
ggsave("figs/calibration.png", width = 7, height = 7, calibration_plot, scale = 1/2)
ggsave("figs/index_devs.png", width = 7, height = 7, index_devs, scale = 1/2)

par_df <- res_df %>%
  mutate(rd = map(rdata, readRDS))


par_df2 <- map_df(par_df$rd, ~ tibble(repl = pluck(., "spec", "repl"),
                                      opmod = pluck(., "spec", "opmod"),
                                      estmod = pluck(., "spec", "estmod"),
                                      par = names(pluck(., "sdr", "par.fixed")),
                                      val = pluck(., "sdr", "par.fixed"),
                                      var = diag(pluck(., "sdr", "cov.fixed")),
                                      sd = sqrt(var))) %>%
  mutate(n_comm = as.numeric(substr(estmod, 5, 99)),
         n_survey = c(50, 100, 200, 300, 400)[opmod])

par_df2 %>%
  filter(par == "lambda_n") %>%
  ggplot(aes(x = val, color = factor(n_comm))) +
  geom_density() +
  facet_wrap(~ factor(n_survey)) +
  labs(x = "λ",
       color = "Commercial\nobservations") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  theme(axis.line.x = element_line(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave("figs/lambda_val.svg", scale = 1/2)

par_df2 %>%
  filter(par == "lambda_n") %>%
  ggplot(aes(x = sd, color = factor(n_comm))) +
  geom_density() +
  facet_wrap(~ factor(n_survey)) +
  labs(x = "Standard error of λ",
       color = "Commercial\nobservations") +
  theme_minimal() +
  theme(axis.line.x = element_line(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave("figs/lambda_sd.svg", scale = 1/2)

## Bias plots
bias_df %>%
  mutate(n_survey = c(50, 100, 200, 300, 400)[opmod],
         n_comm = as.numeric(substr(estmod, 5, 10))) %>%
  ggplot(aes(x = n_survey, y = delta, color = factor(n_comm))) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x = "Number of survey stations",
       y = "Bias measure",
       color = "Commercial\neffort") +
  theme_minimal()
ggsave("figs/bias_plot2.svg", scale = 1/2)

rmse_df %>%
  mutate(n_survey = c(50, 100, 200, 300, 400)[opmod],
         n_comm = as.numeric(substr(estmod, 5, 10))) %>%
  ggplot(aes(x = n_survey, y = rmse, color = factor(n_comm))) +
  geom_point() +
  geom_line() +
  labs(x = "Number of survey stations",
       y = "RMSE",
       color = "Commercial\neffort") +
  theme_minimal()
ggsave("figs/rmse_plot.svg", scale = 1/2)

rep_df <- map_df(par_df$rd, ~ tibble(repl = pluck(., "spec", "repl"),
                                     opmod = pluck(., "spec", "opmod"),
                                     estmod = pluck(., "spec", "estmod"),
                                     par = "rho_sp",
                                     val = pluck(., "rep", "rho_sp", 1)))
