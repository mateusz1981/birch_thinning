setwd("C:/Users/mali/OneDrive - Skogforsk/BJÖRK/Björk_gallrig")

exp = "F1248"


library(readxl)
library(tidyverse)
stems <- read_excel("DB_S1325_SödraVi.xlsx", sheet = "Data", na = ".") %>%
  filter(Exp == exp)
df <- stems



gal <- read_excel("DB_S1325_SödraVi.xlsx", sheet = "Försök", na = ".") %>%
  filter(EXP == exp) %>%
  select(YTA, AR, UTF)  

beh <- read_excel("DB_S1325_SödraVi.xlsx", sheet = "Försök", na = ".") %>% 
  filter(EXP == exp) %>%
  select(YTA, BEH) %>%
  unique()  
  

area <- read_excel("DB_S1325_SödraVi.xlsx", sheet = "Försök", na = ".") %>% 
  filter(EXP %in% exp) %>%
  select(YTA, AREA) %>% 
  unique() %>%
  mutate(AREA = 10000/AREA)  
  



df_total <- data.frame(YTA = numeric(), AR = numeric(), a = numeric(), b = numeric())
unique_YTA <- unique(df$YTA)
unique_AR <- unique(df$AR)


library(minpack.lm)

# Initialize an empty data frame to store the results
df_ar <- data.frame(YTA = numeric(), AR = numeric(), a = numeric(), b = numeric())

# Get unique values of  AR
unique_AR <- unique(df$AR)

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
  }}
  


df_ar_yta <- data.frame(YTA = numeric(), AR = numeric(), a = numeric(), b = numeric())
unique_YTA <- unique(df$YTA)
unique_AR <- unique(df$AR)

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

# Assuming your dataframes are named df_ar and df_ar_yta

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
ab_total <- bind_rows(df_ar_yta, new_rows) %>% arrange(YTA,AR)






df1 <- left_join(df, df_total, by = c("YTA", "AR"))
df2 <- left_join(df1, beh, by = c("YTA")) %>%
  mutate(hest = a*D^b) %>% select(-a, -b)




birch_volume <- function(D, H, na.rm = FALSE){
  
  if (na.rm) {
    D <- na.omit(D)
    H <- na.omit(H)
  }
  
  V <- (10^(-0.89363) * D^2.23818 * (D + 20.0)^(-1.06930) * H^6.02015 * (H - 1.3)^(-4.51472))/1000
  return(V)
}


basalA <- function(D, na.rm = FALSE){
  if(na.rm) {
    D <- na.omit(D)
  }
  
  BA <- pi*(D/2000)^2
}

df2 <- df2 %>% mutate(vol = birch_volume(D/10, hest/10), ba = basalA(D))

df2 <- left_join(df2, area, by = "YTA")

head(df2)





kvar <- df2 %>% group_by(YTA, BEH,  AR) %>% filter(!is.na(D) | MORT3 == 1) %>%
  summarise(sumvol = sum(vol, na.rm = T)*mean(AREA), sumba = sum(ba, na.rm = T)*mean(AREA)) 

kvar

utgall <- df2 %>% group_by(YTA, BEH, AR) %>% filter(MORT3 == 0) %>%
  summarise(sumvolGall = sum(vol, na.rm = T)*mean(AREA), sumbaGall = sum(ba, na.rm = T)*mean(AREA))


dod <- df2 %>% filter(MORT3 == 5) %>% group_by(YTA, BEH, AR)  %>% 
  summarise(sumvolDod = sum(vol, na.rm = T)*mean(AREA), sumbaDod = sum(ba, na.rm = T)*mean(AREA))




volume <- left_join(kvar, utgall, by = c("YTA", "BEH", "AR")) %>% ungroup() %>%
  mutate(sumvolGall  = ifelse(is.na(sumvolGall ), 0, sumvolGall )) %>%
  mutate(kvarvol = sumvol - sumvolGall) %>% group_by(YTA, BEH) %>%
  mutate(utgall = cumsum(sumvolGall)) %>%
  mutate(totvol = kvarvol + utgall) %>%
  filter(!is.na(totvol) & sumvol != 0)



volume %>% filter(YTA == 5)

### test anova and figures
res <- df2 %>% filter(AR == 2020) %>% filter(GALL == 1) %>% group_by(YTA, BEH) %>%
  summarise(mh = mean(H, na.rm = T), D = mean(D, na.rm = T))



ggplot(aes(x = AR, y = totvol, color = BEH), data = volume) + geom_line() +
  facet_wrap(YTA~BEH)

ggplot(aes(x = AR, y = totvol), data = volume) + 
  geom_line(aes(group = YTA, color = BEH), linewidth = 1) +
  theme_bw()

ggplot(aes(x = AR, y = totvol), data = volume) + 
  geom_bar(stat = "identity", position = "dodge", aes(group = YTA, fill = BEH)) +
  theme_bw() +
  facet_wrap(YTA~BEH)


ggplot(aes(x = AR, y = totvol, group = BEH), data = volume %>% filter(AR != 2013)) + 
  stat_summary(fun.y = "mean", aes(group = BEH, fill = BEH ), geom = "bar", position = "dodge")+
  theme_bw() +
  facet_wrap(YTA~BEH)



#############basal area##############################
basalarea <- left_join(kvar %>% select(-sumvol), utgall %>% select(-sumvolGall), by = c("YTA", "BEH", "AR")) %>% ungroup() %>%
  mutate(sumbaGall  = ifelse(is.na(sumbaGall ), 0, sumbaGall )) %>%
  mutate(kvarba = sumba - sumbaGall) %>% group_by(YTA, BEH) %>%
  mutate(utgallba = cumsum(sumbaGall)) %>%
  mutate(totba = kvarba + utgallba) %>%
  filter(!is.na(totba) & sumba != 0)



basalarea <- basalarea[rep(seq_len(nrow(basalarea)), each = 2), ]
basalarea[is.na(basalarea)] <- 0

for (i in seq(1, length(basalarea$AR), 2)){
  basalarea$kvarba[i] = basalarea$kvarba[i] + basalarea$sumbaGall [i]
}
basalarea
ggplot(aes(x = AR, y = kvarba), data = basalarea) + geom_line() + facet_wrap(YTA~BEH )

ggplot(aes(x = AR, y = kvarba, group = YTA, color = factor(BEH)), data = basalarea) + geom_line() + geom_point() + facet_wrap(~BEH )

##############volym figure###########################
volume
volume <- volume[rep(seq_len(nrow(volume)), each = 2), ]
volume[is.na(volume)] <- 0
volume
for (i in seq(1, length(volume$AR), 2)){
  volume$kvarvol[i] = volume$kvarvol[i] + volume$sumvolGall[i]
}
volume
ggplot(aes(x = AR, y = kvarvol), data = volume) + geom_line() + facet_wrap(YTA~BEH )



##############MAI################
ggplot(aes(x = AR-1999+2, y = totvol/(AR-2000+2)), data = volume) + 
  geom_line(aes(group = YTA, color = BEH), linewidth = 1) +
  xlab("Ålder") + ylab("MAI (m3 ha år)") +
  theme_bw()

unique(df2$MORT)

head(df2)
antal <- df2 %>% group_by(YTA,BEH, AR, MORT3)  %>%
  summarise(n= round(sum(!is.na(D), na.rm = T)*mean(AREA), 0)) %>% 
  filter(!is.na(n) & n != 0)



########## Antal rätt################
antal1 <- df2 %>% filter(!is.na(MORT3)) %>% group_by(YTA, BEH, AR, MORT3) %>%
  summarise(li = round(sum(!is.na(MORT3)) * mean(AREA, na.rm = T), 0))
antal1

antal1 %>% filter(YTA == 11)

antal1 <- reshape2::dcast(antal1, YTA + BEH + AR ~ MORT3, value.var = "li")
antal1
names(antal1)[names(antal1)=="1"] <- "n"
names(antal1)[names(antal1)== "0"] <- "utgall_n"





antal1 <- antal1[rep(seq_len(nrow(antal1)), each = 2), ]
antal1[is.na(antal1)] <- 0
antal1
for (i in seq(1, length(antal1$AR), 2)){
  antal1$n[i] = antal1$n[i] + antal1$utgall_n[i]
}
antal1 %>% filter(YTA == 10)
ggplot(aes(x = AR, y = n), data = antal1) + geom_line() + facet_wrap(YTA~BEH )
#################################
