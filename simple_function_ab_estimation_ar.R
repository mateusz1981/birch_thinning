fit_simple_model <- function(df)  {df_total <- data.frame(AGE = numeric(), a = numeric(), b = numeric())

unique_age <- unique(df$AGE)


# Initialize an empty data frame to store the results
df_age <- data.frame(AGE = numeric(), a = numeric(), b = numeric())

# Loop to fit models for AR
for (age in unique_age) {
  da_age <- subset(df, AGE == age)
  mod_age <- nlsLM(H ~ a * D^b, start = c(a = 2, b = 0.5), data = da_age)

    co_age <- data.frame(AGE = age, a = coef(mod_age)[1], b = coef(mod_age)[2])
    df_age <- rbind(df_age, co_age, row.names = FALSE) %>% filter(AGE > 0)
    rownames(df_age) <- NULL
 
  }
  aa <- df_age
  return(df_age)
}



