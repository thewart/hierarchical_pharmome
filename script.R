library(tidyverse)
library(tidybayes)
library(rstan)
obs <- read_csv("downloadable/observed_points_table.csv")
eve_res <- read_csv("downloadable/detailed_result_table.csv") 
obs_screen_antag <- filter(obs, Phase == "Screening (3pt)" & Mode == "Antagonist")


to_drop <- obs_screen_antag |>
  filter(Activity > 150 | Activity < -50) |> distinct(Compound_ID, Target_ID) |> 
  reframe(str_c(Compound_ID, Target_ID)) |> pull()

obs_final <- obs_screen_antag |> filter(!(str_c(Compound_ID, Target_ID) %in% to_drop))

names_to_groups <- obs_final |>
  reframe(Compound_ID, Target_ID, group = str_c(Compound_ID, Target_ID)) |>
  mutate(group = factor(group) |> as.numeric())

stan_data <- list(
  N = nrow(obs_final),
  M = obs_final |> distinct(Compound_ID, Target_ID) |> nrow(),
  y = obs_final |> pull(Activity)/100,
  x = obs_final |> pull(Concentration) |> dense_rank() - 2,
  group = obs_final |> pull(group)
)

model_doublet <- stan_model("model.stan")
fit <- sampling(model_doublet, stan_data, iter = 300, chains = 4, warmup = 200,
                sample_file = "doublet_samples_12tree", control = list(max_treedepth = 12))

samps <- rstan::extract(fit, pars = c(
  "alpha", "beta", "sigma_alpha", "sigma_beta", "lognu_alpha"
  ))

nrep <- 10
x <- c(-1, 0, 1)
trueff <- list(vector(), vector(), vector())
for (i in 1:10) {
  for (j in 1:3) {
    beta <- with(samps, beta + sigma_beta * rnorm(length(beta)))
    alpha <- with(samps, alpha + sigma_alpha * rt(length(alpha), exp(lognu_alpha)))
    tmp <- pnorm(alpha + x[j] * beta) * 100
    trueff[[j]] <- c(trueff[[j]], tmp)
  }
}

mybreaks <- c(0, 1, 5, 15, 30, 50, 75, 100)
trueff_prop <- predprop <- sapply(trueff, \(x) {
    x |> cut(mybreaks) |> table() |> prop.table()
  }) |> as_tibble(rownames = "Activity") |> 
  pivot_longer(starts_with("V"), values_to = "Proportion", names_to = "Concentration") |> 
  mutate(Concentration = case_match(Concentration, "V1" ~ 6.25e-7, "V2" ~ 2.5e-6, "V3" ~ 1e-5) |> ordered()) |>
  mutate(Activity = ordered(Activity, levels = unique(Activity)))

mybreaks <- seq(0, 100, 5)
indests <- spread_draws(fit, alpha, sigma_alpha, beta, sigma_beta, alpha_z[i], beta_z[i])
indests <- indests |> reframe(.draw, i, alpha = alpha + sigma_alpha * alpha_z, beta = beta + sigma_beta * beta_z) |> 
  reframe(.draw, i, 
          act1 = (pnorm(alpha - beta) * 100) |> cut(mybreaks),
          act2 = (pnorm(alpha) * 100) |> cut(mybreaks),
          act3 = (pnorm(alpha + beta) * 100) |> cut(mybreaks)) |>
  pivot_longer(starts_with("act")) |> group_by(i, name) |> 
  count(value) |> mutate(p = n/sum(n)) |> ungroup()

eve_res <- eve_res |> filter(Mode == "Antagonist" & !(str_c(Compound_ID, Target_ID) %in% to_drop))

indests <- indests |>
  mutate(name = case_match(name, "act1" ~ 6.25e-7, "act2" ~ 2.5e-6, "act3" ~ 1e-5) |> ordered()) |> 
  rename(group = "i", Concentration = "name", Activity = "value") |> 
  left_join(names_to_groups |> distinct()) |>
  left_join(reframe(everes, Target_ID, Compound_ID, Scr_Category_Simplified, Compound))

everes |> count(Scr_Category_Simplified) |> mutate(p = n/sum(n))
