# test GAMs

rm(list = ls())

if(!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}
library(here)

source(here::here("utils", "functions.R"))

install_load_packages(c(
  "tidyverse",
  "tictoc",
  "grf",
  "kernelshap",
  "hstats",
  "glue",
  "shapviz",
  "patchwork",
  "cobalt",
  "purrr",
  "corrplot",
  "gridExtra",
  "rlang",
  "arrow",
  "tictoc",
  "mgcv",
  "gratia"))

dir_figs <- here::here('figs', 'exploratory')
dir_derived <- here::here('data', 'derived')
dir_raw <- here::here('data', 'raw')
dir_manual <- here::here('data', 'manual')
dir_ensure(c(dir_figs,
             dir_derived,
             dir_manual,
             dir_raw))

seed = 1234
set.seed(seed)

dats_long <- arrow::read_parquet(here(dir_derived, "test_dats_long.parquet"))

x <- dats_long |>
  filter(fire == 1 & uniqueID == 89) |>
  filter(year_from_fire_index >= -20) |>
  select(rap_tree, year)



fit <- gam(rap_tree ~ s(year, k = 8),
           family = gaussian(),
           data = x,
           method = "REML")
summary(fit)


b0 <- coef(fit)[["(Intercept)"]]

draw(fit) +
  ggplot2::scale_y_continuous(
    name = "rap_tree (fitted)",
    labels = function(z) z + b0
  ) + geom_vline(xintercept = 2000)

gam.check(fit)
k.check(fit)



year_grid <- data.frame(
  year = seq(min(x$year), max(x$year), length.out = 200)
)

pred <- predict(fit, newdata = year_grid, se.fit = TRUE, type = "response")

year_grid$fit <- pred$fit
year_grid$lwr <- pred$fit - 1.96 * pred$se.fit
year_grid$upr <- pred$fit + 1.96 * pred$se.fit

ggplot(x, aes(year, rap_tree)) +
  geom_point(alpha = 0.35) +
  geom_ribbon(
    data = year_grid,
    aes(x = year, ymin = lwr, ymax = upr),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  geom_line(
    data = year_grid,
    aes(x = year, y = fit),
    inherit.aes = FALSE
  ) +
  geom_vline(xintercept = 2000)

