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
library(hdf5r)
library(patchwork)

prep_file <- "~/gscratch/spatq_sims/prep.h5"

study <- "counterpref"
repl <- 1
opmods <- 1:5
max_T <- 15

## preph5 <- h5file(prep_file, "r")
poph5_files <- map_chr(opmods, ~ sim_file_paths(study, repl, ., "sims")$poph5)

poph5_path <- sim_file_paths(study, repl, 1, "sims")$poph5
poph5 <- h5file(poph5_path)
init_pop <- poph5[["popstate"]][, , seq_len(max_T)]

pop_df <- cross_df(list(s1 = 0.5:99.5,
                        s2 = 0.5:99.5,
                        time = seq_len(max_T))) %>%
  mutate(pop = as.vector(init_pop))

pop_df %>%
  filter(time == 1) %>%
  ggplot(aes(x = s1, y = s2, fill = pop)) +
  geom_raster() +
  scale_fill_viridis_c() +
  coord_equal() +
  lims(x = c(0, 100), y = c(0, 100)) +
  labs(fill = "Abundance") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank())
ggsave("figs/initial_abundance.svg", scale = 1/2)

pop_df %>%
  ggplot(aes(x = s1, y = s2, fill = pop)) +
  geom_raster() +
  scale_fill_viridis_c() +
  coord_equal() +
  facet_wrap(~ time, nrow = 3) +
  lims(x = c(0, 100), y = c(0, 100)) +
  labs(fill = "Abundance") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank())
ggsave("figs/abundance_series.svg", scale = 1/2)

catch1_df <- read_catch(study, repl, 5, "sims") %>%
  filter(time <= 15)
pop1_df <- read_popstate(study, repl, 5, "sims") %>%
  filter(time <= 15)

tot_catch_df <- catch1_df %>%
  group_by(time, vessel_idx) %>%
  summarize(tot_catch = sum(catch_biomass)) %>%
  mutate(Vessel = factor(vessel_idx,
                         labels = c("Survey", "Commercial")))

cplot <- ggplot(tot_catch_df,
                mapping = aes(x = time,
                              y = tot_catch,
                              color = Vessel)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0) +
  guides(color = FALSE) +
  facet_grid(Vessel ~ ., scales = "free_y") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  ## lims(y = c(0, NA)) +
  labs(x = "Year", y = "Total catch") +
  theme_minimal() +
  theme(axis.line = element_line())
cplot

pplot <- ggplot(pop1_df, aes(x = time, y = pop)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  ## lims(y = c(0, NA)) +
  labs(x = "Year", y = "Population biomass") +
  theme_minimal() +
  theme(axis.line = element_line())
pplot

pplot / cplot
ggsave("figs/abundcatch_tots.svg", scale = 1/2)

## -------
## Habitat
preph5 <- h5file("../spatq_sims/prep.h5")
hab1 <- preph5[["habitat"]][, , 1]
hab_df <- cross_df(list(s1 = 0.5:99.5,
                     s2 = 0.5:99.5)) %>%
  mutate(habitat = as.vector(hab1))

hab_df %>%
  ggplot(aes(x = s1, y = s2, fill = habitat)) +
  geom_raster() +
  scale_fill_viridis_c(option = "E") +
  coord_equal() +
  lims(x = c(0, 100), y = c(0, 100)) +
  labs(fill = "Habitat\ncovariate") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank())
ggsave("figs/habitat.svg", scale = 1/2)
