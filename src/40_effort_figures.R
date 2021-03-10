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

study <- "counterpref"
repl <- 1
estmods <- c("comm500",
             "comm1000",
             "comm2000")
opmods <- 1:5

catch_df <- read_catch(study, repl, opmods[1], root_dir = "sims")

catches <- tibble(opmod = opmods) %>%
  mutate(catch_df = map(.$opmod, ~ read_catch(study, repl, .x, "sims"))) %>%
  unnest(cols = catch_df) %>%
  filter(time == 1) %>%
  group_by(opmod, time, vessel_idx, s1, s2) %>%
  summarize(eff = sum(effort)) %>%
  mutate(n_survey = c(50, 100, 200, 300, 400)[opmod])

## Plot commercial and survey together; difficult to read
ggplot(mapping = aes(x = s1, y = s2)) +
  geom_tile(aes(fill = eff),
            filter(catches, vessel_idx == 2)) +
            geom_point(data = filter(catches, vessel_idx == 1),
                       size = 1) +
  scale_fill_viridis_b(option = "B") +
  coord_equal() +
  facet_wrap(~ n_survey, nrow = 2) +
  lims(x = c(0, 100), y = c(0, 100)) +
  labs(fill = "Commercial\neffort") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank())

## Plot only commercial effort
catches %>%
  filter(vessel_idx == 2,
         opmod == 1) %>%
ggplot(mapping = aes(x = s1, y = s2, fill = eff)) +
  geom_tile() +
  scale_fill_viridis_b(option = "B") +
  coord_equal() +
  lims(x = c(0, 100), y = c(0, 100)) +
  labs(fill = "Commercial\neffort") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank())
ggsave("figs/commercial_effort.svg", scale = 1/2)

## Plot survey stations, fade commercial effort
ggplot(mapping = aes(x = s1, y = s2)) +
  geom_tile(aes(fill = eff),
            filter(catches, vessel_idx == 2),
            alpha = 0.5) +
  geom_point(data = filter(catches, vessel_idx == 1),
                       size = 1) +
  scale_fill_viridis_b(option = "B") +
  coord_equal() +
  facet_wrap(~ n_survey, nrow = 2) +
  lims(x = c(0, 100), y = c(0, 100)) +
  labs(fill = "Commercial\neffort") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank())

## Plot each survey effort distribution
walk(opmods,
     function (om) {
       df <- filter(catches, opmod == om)
       ## Plot survey stations, fade commercial effort
       ggplot(mapping = aes(x = s1, y = s2)) +
         geom_tile(aes(fill = eff),
                   filter(df, vessel_idx == 2),
                   alpha = 0.5) +
         geom_point(data = filter(df, vessel_idx == 1),
                    size = 1.5) +
         scale_fill_viridis_b(option = "B", guide = FALSE) +
         coord_equal() +
         ## facet_wrap(~ n_survey, nrow = 2) +
         lims(x = c(0, 100), y = c(0, 100)) +
         theme_minimal() +
         theme(axis.title = element_blank(),
               axis.ticks = element_blank(),
               axis.text = element_blank(),
               panel.grid = element_blank())
       fn <- paste0("figs/survey_effort_", df$n_survey[1], ".svg")
       ggsave(fn, scale = 1/2)
     })

preph5 <- h5file("../spatq_sims/prep.h5")
pop <- preph5[["spatial_eq"]][, , 1]
binpop_df <- cross_df(list(s1 = 0.5:99.5,
                           s2 = 0.5:99.5)) %>%
  mutate(pop = as.vector(pop)) %>%
  mutate(g1 = cut(s1, seq(0, 100, 5)),
         g2 = cut(s2, seq(0, 100, 5))) %>%
  group_by(g1, g2) %>%
  summarize(tot_pop = sum(pop),
            s1 = mean(s1),
            s2 = mean(s2),
            .groups = "drop") %>%
  mutate(iwt1 = 2 * mean(tot_pop) - tot_pop,
         iwt2 = iwt1 - min(iwt1),
         iwt3 = iwt2 + 0.1 * iwt2,
         st_prob = iwt3 / sum(iwt3)) %>%
  select(s1, s2, st_prob)

ggplot(binpop_df, aes(x = s1, y = s2, fill = st_prob)) +
  geom_raster() +
  geom_point() +
  coord_equal() +
  lims(x = c(0, 100), y = c(0, 100)) +
  labs(fill = "Station\nselection\nprobability") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank())
ggsave("figs/station_prob.svg", scale = 1/2)
