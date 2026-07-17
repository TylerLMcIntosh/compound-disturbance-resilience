here()

agg_obj <- readRDS(here("results-exo/v8/tables/by_run/nfg_aspen_birch__rap_tree__sixyr_b90global_pdsisum90global__sunab_twfe_unweighted_6yr/agg_event_study__conley_75km_5km.rds"))

names(agg_obj)

agg_obj$coef
g <- agg_obj$groups
vcov <- agg_obj$vcov


# check wald computations
years <- c(4:6)
group_a <- "cd_f"
group_b <- "cd_bf"

idx_a <- which(
  agg_obj$groups$dummy_group == group_a &
    agg_obj$groups$event_time %in% years
)

idx_b <- which(
  agg_obj$groups$dummy_group == group_b &
    agg_obj$groups$event_time %in% years
)

L <- numeric(length(agg_obj$coef))

L[idx_a] <- 1 / length(idx_a)
L[idx_b] <- -1 / length(idx_b)

est <- sum(L * agg_obj$coef)


se <- sqrt(
  as.numeric(
    t(L) %*% agg_obj$vcov %*% L
  )
)


test <- t(L) %*% agg_obj$vcov
tt <- test %*% L


# simulate data
set.seed(123)

n <- 100

x1 <- rnorm(n, mean = 10, sd = 2)
x2 <- rnorm(n, mean = 5, sd = 1.5)

y <- 3 + 2 * x1 - 1.5 * x2 + rnorm(n, mean = 0, sd = 2)


fit <- lm(y ~ x1 + x2)

summary(fit)
vcov(fit)


# manually derive vcov
X <- model.matrix(fit)

sigma2 <- sum(residuals(fit)^2) / df.residual(fit)

vcov_manual <- sigma2 * solve(t(X) %*% X)


# Compare naive and robust vcov

# Coefficients and variance-covariance matrix
b <- coef(fit)
V <- vcov(fit)

# Contrast: beta_x1 - beta_x2
contrast_estimate <- b["x1"] - b["x2"]

# Correct variance and SE, including covariance
contrast_variance_correct <-
  V["x1", "x1"] +
  V["x2", "x2"] -
  2 * V["x1", "x2"]

contrast_se_correct <- sqrt(contrast_variance_correct)

# Naive variance and SE, incorrectly assuming independence
contrast_variance_naive <-
  V["x1", "x1"] +
  V["x2", "x2"]

contrast_se_naive <- sqrt(contrast_variance_naive)

# Compare results
data.frame(
  contrast = "x1 - x2",
  estimate = unname(contrast_estimate),
  covariance_x1_x2 = V["x1", "x2"],
  se_correct = unname(contrast_se_correct),
  se_assumed_independent = unname(contrast_se_naive)
)










