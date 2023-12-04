library(minpack.lm)
library(dplyr)

fit_models <- function(df) {
  df_total <- data.frame(YTA = numeric(), AR = numeric(), a = numeric(), b = numeric())
  unique_YTA <- unique(df$YTA)
  unique_AR <- unique(df$AR)
  
  # Initialize an empty data frame to store the results
  df_ar <- data.frame(YTA = numeric(), AR = numeric(), a = numeric(), b = numeric())
  
  # Loop to fit models for AR
  for (ar in unique_AR) {
    da_ar <- subset(df, AR == ar)
    mod_ar <- nlsLM(H ~ a * D^b, start = c(a = 2, b = 0.5), data = da_ar)
    
    if (!inherits(mod_ar, "try-error")) {
      co_ar <- data.frame(AR = ar, a = coef(mod_ar)[1], b = coef(mod_ar)[2])
      df_ar <- rbind(df_ar, co_ar, row.names = FALSE) %>% filter(AR > 0)
      rownames(df_ar) <- NULL
    } else {
      cat("Error fitting model for AR =", ar, ": ", conditionMessage(mod_ar), "\n")
      co_ar <- NULL
    }
  }
  
  df_ar_yta <- data.frame(YTA = numeric(), AR = numeric(), a = numeric(), b = numeric())
  
  for (yta in unique_YTA) {
    for (ar in unique_AR) {
      da <- subset(df, YTA == yta & AR == ar)
      
      tryCatch({
        mod <- nlsLM(H ~ a * D^b, start = c(a = 2, b = 0.5), data = da)
        
        if (!inherits(mod, "try-error")) {
          co <- data.frame(YTA = yta, AR = ar, a = coef(mod)[1], b = coef(mod)[2])
          df_ar_yta <- rbind(df_ar_yta, co, row.names = FALSE) %>% filter(AR > 0)
          rownames(df_ar_yta) <- NULL
        } else {
          cat("Model did not converge for YTA =", yta, "and AR =", ar, "\n")
          cat("Retrying with all observations for AR =", ar, "\n")
        }
      }, error = function(e) {
        cat("Error fitting model for YTA =", yta, "and AR =", ar, ": ", conditionMessage(e), "\n")
      })
    }
  }
  
  # Identify missing combinations of YTA and AR
  missing_combinations <- expand.grid(unique(df_ar_yta$YTA), unique(df_ar$AR))
  missing_combinations <- missing_combinations[!do.call(paste, missing_combinations) %in% paste(df_ar_yta$YTA, df_ar_yta$AR), ]
  
  # Add missing rows from df_ar to df_ar_yta
  new_rows <- lapply(seq_len(nrow(missing_combinations)), function(i) {
    yta <- missing_combinations$Var1[i]
    ar <- missing_combinations$Var2[i]
    row_df_ar <- df_ar[df_ar$AR == ar, ]
    row_df_ar$YTA <- yta
    return(row_df_ar)
  })
  
  # Combine the original df_ar_yta and the new rows
  ab_total <- bind_rows(df_ar_yta, new_rows) %>% arrange(YTA, AR)
  
  return(ab_total)
}

df_total_result <- fit_models(df)
