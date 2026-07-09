


# A function to use sufrep::make_encoder such that the encoded variables are added to the dataframe
sufrep_encode_append <- function(dats,
                                 method = 'means',
                                 X_cols,
                                 G_col,
                                 remove_G_col = TRUE,
                                 ...) {
  
  #Select variables for encoding
  x_df <- dats |> 
    select(all_of(X_cols)) |>
    as.data.frame()
  g_df <- dats |>
    pull(G_col)
  
  #Encode
  encoder <- sufrep::make_encoder(method = method, X = x_df, G = g_df, ...)
  dats_encoded <- encoder(x_df, G = g_df)
  
  #Make interpretable encoding columns
  encoded_cols <- dats_encoded |> select(starts_with("ENC"))
  names(encoded_cols) <- paste0(names(encoded_cols), "_", G_col)
  
  dats_out <- dats |>
    cbind(encoded_cols)
  
  if(remove_G_col) {
    dats_out <- dats_out |>
      select(-all_of(G_col))
  }
  
  return(dats_out)
}



create_ml_data <- function(dats, scale = FALSE,
                           cont_covars, cat_covars, outcome, trts, cluster) {
  
  all_cols <- c(cont_covars, cat_covars, outcome, trts, cluster)
  
  # Scale all continuous covariates
  # filter out anything with NAs
  # ensure categorical covariates are factors
  # encode categorical as one-hot dummies or sufficient representations
  dats_ml <- dats |>
    select(all_of(all_cols)) |>
    filter(dplyr::if_all(everything(), ~ !is.na(.x)))
  
  if(scale) {
    dats_ml <- dats_ml |>
      mutate(
        across(
          all_of(cont_covars),
          \(x) as.numeric(scale(x))
        )
      )
  }
  dats_ml <- dats_ml |>
    mutate(
      across(
        all_of(cat_covars),
        \(x) forcats::fct_drop(as.factor(x))
      )
    )
  
  # for(var in cat_covars) {
  #   dats_ml <- dats_ml |>
  #     sufrep_encode_append(method = 'low_rank', # Use low-rank sufficient representation encoding to reduce to 2 dimensions
  #                          X_cols = cont_covars,
  #                          G_col = var,
  #                          num_components = 2)
  # }
  # 
  
  for(var in cat_covars) {
    
    n_levels <- dplyr::n_distinct(dats_ml[[var]], na.rm = TRUE)
    
    if(n_levels <= 1) {
      message("Skipping ", var, ": only ", n_levels, " level after filtering.")
      
      dats_ml <- dats_ml |>
        dplyr::select(-dplyr::all_of(var))
      
      next
    }
    
    encoded_successfully <- FALSE
    
    for(n_comp_use in rev(seq_len(min(2, n_levels - 1)))) {
      
      message("Trying to encode ", var, " with ", n_comp_use, " component(s).")
      
      dats_try <- tryCatch(
        {
          sufrep_encode_append(
            dats = dats_ml,
            method = 'low_rank',
            X_cols = cont_covars,
            G_col = var,
            num_components = n_comp_use
          )
        },
        error = function(e) {
          message("Failed: ", conditionMessage(e))
          NULL
        }
      )
      
      if(!is.null(dats_try)) {
        dats_ml <- dats_try
        encoded_successfully <- TRUE
        break
      }
    }
    
    if(!encoded_successfully) {
      message("Skipping ", var, ": encoding failed for all attempted component counts.")
      
      dats_ml <- dats_ml |>
        dplyr::select(-dplyr::all_of(var))
    }
  }
  
  enc_xvars <- dats_ml |>
    select(starts_with("ENC")) |>
    names()
  
  return(list('dats_ml' = dats_ml,
              'enc_xvars' = enc_xvars))
  
}





cf_analysis <- function(dats,
                        yvar,
                        xvars,
                        trt,
                        data_nm,
                        seed = 1234,
                        dir_out,
                        ...) {
  # Filter out rows with NA in the outcome variable
  dats_filt <- dats %>%
    filter(!is.na(.data[[yvar]]))
  
  n_total <- nrow(dats_filt)
  n_trt <- nrow(dats_filt |>
                  filter(.data[[trt]] == 1))
  
  # Select covariates and outcome
  covariates <- dats_filt %>%
    select(all_of(xvars))
  outcome <- dats_filt[[yvar]]
  treatment <- dats_filt[[trt]]
  
  # Fit causal forest
  fit <- grf::causal_forest(
    X = covariates,
    Y = outcome,
    W = treatment,
    seed = seed,
    ...
  )
  
  # Variable importance - ensure normalized
  imp_raw <- sort(setNames(variable_importance(fit), xvars))
  imp_normalized <- imp_raw / sum(imp_raw)
  
  imp_df <- tibble(
    variable = names(imp_normalized),
    importance = as.numeric(imp_normalized)
  ) %>%
    arrange(importance)  # So the smallest is at the bottom of the horizontal bar plot
  
  imp_p <- ggplot(imp_df, aes(x = importance, y = reorder(variable, importance))) +
    geom_col(fill = "orange") +
    labs(
      x = "Variable Importance (Normalized)",
      y = NULL,
      title = glue("Causal Forest Variable Importance, yvar = {yvar}"),
      subtitle = data_nm,
      caption = glue("Total sample size: {n_total}\nTreated sample size: {n_trt}")
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.y = element_text(hjust = 1),
      plot.margin = margin(10, 20, 5, 5)
    )
  ggsave(plot = imp_p,
         filename = here(dir_out, glue("var_importance_{data_nm}_{yvar}.png")),
         bg = "white")
  
  pred_fun <- function(object, newdata, ...) {
    predict(object, newdata, ...)$predictions
  }
  
  xvars_sorted <- rev(names(imp_normalized))
  
  pdps <- lapply(xvars_sorted, function(v) {
    pd <- partial_dep(fit, v, X = covariates, pred_fun = pred_fun)
    plot(pd) +
      theme(axis.title.y = element_blank())
  })
  
  #For adding y axis label to overall plot
  ylab_plot <- ggplot() +
    theme_void() +
    annotate("text", x = 0, y = 0.5, label = "Treatment effect", angle = 90, hjust = 0.5, size = 4)
  
  pdps_p <- wrap_plots(pdps, guides = "collect", ncol = 3) #&
  #ylab("Treatment effect")
  
  # Add format & add titles
  pdps_p <- ylab_plot + pdps_p +
    plot_layout(widths = c(0.05, 1)) +
    plot_annotation(
      title = glue("Partial Dependence Plots: trt = {yvar}"),
      subtitle = data_nm,
      caption = glue("NOTES:\nNegative values mean a stronger (positive) treatment effect (i.e. fire is causing a greater difference in forest cover change)\nPDPs are in order of variable importance\nTotal sample size: {n_total}\nBurned sample size: {n_trt}"))
  
  ggsave(filename = here(dir_out, glue("partial_dependence_{data_nm}_{yvar}.png")),
         plot = pdps_p,
         units = 'px',
         width = 2000,
         height = 2000)
  
  return(list(fit = fit,
              imp_p = imp_p,
              pdps = pdps))
}



