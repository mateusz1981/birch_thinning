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
gal
beh <- read_excel("DB_S1325_SödraVi.xlsx", sheet = "Försök", na = ".") %>% 
  filter(EXP == exp) %>%
  select(YTA, BEH) %>%
  unique()  
  

area <- read_excel("DB_S1325_SödraVi.xlsx", sheet = "Försök", na = ".") %>% 
  filter(EXP %in% exp) %>%
  select(YTA, AREA) %>% 
  unique() %>%
  mutate(AREA = 10000/AREA)  
  


source("function_ab_estimation.R")
df_total <- fit_models(df)





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

tail(df2)


####volume calculations


kvar <- df2 %>% group_by(YTA, BEH,  AR) %>% filter(!is.na(D) | MORT3 == 1) %>%
  summarise(sumvol = sum(vol, na.rm = T)*mean(AREA), sumba = sum(ba, na.rm = T)*mean(AREA)) 



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

volume[is.na(volume)] <- 0


#to result table##########################
res_volume <- volume %>% select(YTA, BEH, AR, sumvol, sumvolGall, kvarvol, utgall, totvol) %>%
  rename(volfg = sumvol, volut = sumvolGall, voleg = kvarvol, volut_cum = utgall, Totvol = totvol)

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


#####################################################################################

#############basal area##############################
basalarea <- left_join(kvar %>% select(-sumvol), utgall %>% select(-sumvolGall), by = c("YTA", "BEH", "AR")) %>% ungroup() %>%
  mutate(sumbaGall  = ifelse(is.na(sumbaGall ), 0, sumbaGall )) %>%
  mutate(kvarba = sumba - sumbaGall) %>% group_by(YTA, BEH) %>%
  mutate(utgallba = cumsum(sumbaGall)) %>%
  mutate(totba = kvarba + utgallba) %>%
  filter(!is.na(totba) & sumba != 0)

basalarea[is.na(basalarea)] <- 0


res_ba <- basalarea %>% select(YTA, BEH, AR, sumba, sumbaGall, kvarba, utgallba, totba) %>%
  rename(bafg = sumba, baut = sumbaGall, baeg = kvarba, baut_cum = utgallba, Totba = totba)

basalarea <- basalarea[rep(seq_len(nrow(basalarea)), each = 2), ]



##BA figure##################
for (i in seq(1, length(basalarea$AR), 2)){
  basalarea$kvarba[i] = basalarea$kvarba[i] + basalarea$sumbaGall [i]
}
basalarea
ggplot(aes(x = AR, y = kvarba), data = basalarea) + geom_line() + facet_wrap(YTA~BEH )

ggplot(aes(x = AR, y = kvarba, group = YTA, color = factor(BEH)), data = basalarea) + geom_line() + geom_point() + facet_wrap(~BEH )




##############MAI################
ggplot(aes(x = AR-1995+2, y = totvol/(AR-1995+2)), data = volume) + 
  geom_line(aes(group = YTA, color = BEH), linewidth = 1) +
  xlab("Ålder") + ylab("MAI (m3 ha år)") +
  theme_bw()




########## Antal ################
antal1 <- df2 %>% filter(!is.na(MORT3)) %>% group_by(YTA, BEH, AR, MORT3) %>%
  summarise(li = round(sum(!is.na(MORT3)) * mean(AREA, na.rm = T), 0))



antal1[is.na(antal1)] <- 0
antal1
library(reshape2)
library(dplyr)

process_antal_data <- function(antal1) {
  antal1 <- reshape2::dcast(antal1, YTA + BEH + AR ~ MORT3, value.var = "li")
  
  names(antal1)[names(antal1)== "0"] <- "utgall_n" #removed in thinning
  names(antal1)[names(antal1)=="1"] <- "n" # left in stand
  names(antal1)[names(antal1)== "5"] <- "dod_n" # natural mortality
  if (!all(c("utgall_n", "n", "dod_n") %in% names(antal1))) {
    antal1 <- antal1 %>% mutate(dod_n = 0)
  }
  antal1[is.na(antal1)] <- 0
  res_antal <- antal1 %>% group_by(YTA, BEH) %>%
     mutate(nfg = n + utgall_n + dod_n, nut = utgall_n, neg = n) %>%
    select(-utgall_n, -n) %>%
    mutate(nutcum = cumsum(nut))
  
  return(res_antal)
  
}

# Example usage:
antal1 <- process_antal_data(antal1)

res_antal <- antal1


antal1 <- antal1[rep(seq_len(nrow(antal1)), each = 2), ]
antal1[is.na(antal1)] <- 0




for (i in seq(1, length(antal1$AR), 2)) {
  antal1$neg[i] = antal1$neg[i] + antal1$nut[i] + antal1$dod_n[i] }


ggplot(aes(x = AR, y = neg), data = antal1) + geom_line() + facet_wrap(YTA~BEH )


#################################


############################################################
### results diameter_height##############################
##########################################################
#before thinning
res_height_diameter_fg <- df2 %>% group_by(YTA, BEH, AR) %>%
  summarise(mh_fg = mean(hest, na.rm = T), md_fg = mean(D, na.rm = T)) 


#thinned = 0 and unthinnded = 1
res_height_diameter<- df2 %>% group_by(YTA, BEH, AR, MORT3) %>%
  summarise(mh = round(mean(hest, na.rm = T), 0), md = round(mean(D, na.rm = T), 0)) %>%
  filter(!is.na(MORT3))

#height
hh <- reshape2::dcast(res_height_diameter, YTA + BEH + AR ~ MORT3, value.var = "mh")

names(hh)[names(hh)== "0"] <- "mh_ut" #removed in thinning
names(hh)[names(hh)=="1"] <- "mh_eg" # left in stand'
hh[is.na(hh)] <- "0"
hh
#diameter
dd <- reshape2::dcast(res_height_diameter, YTA + BEH + AR ~ MORT3, value.var = "md")
names(dd)[names(dd)== "0"] <- "md_ut" #removed in thinning
names(dd)[names(dd)=="1"] <- "md_eg" # left in stand'

dd[is.na(dd)] <- "0"

hd <- left_join(res_height_diameter_fg, hh, by = c("YTA", "BEH", "AR"))
hd <- left_join(hd, dd, by = c("YTA", "BEH", "AR")) %>% 
  select("YTA", "BEH", "AR", "mh_fg", "mh_ut", "mh_eg", "md_fg", "md_ut", "md_eg")
hd

rm(pri)
pri <- left_join(res_antal, res_volume, by = c("YTA", "BEH", "AR"))
pri <- left_join(pri, hd, by = c("YTA", "BEH", "AR"))
pri <- left_join(pri, res_ba, by = c("YTA", "BEH", "AR")) %>% mutate(GALL = ifelse(baut > 1, 1, 0)) %>%
  select(YTA, BEH, AR, GALL, everything())


write.csv(pri %>% mutate_if(is.numeric, list(~format(., nsmall = 1))), paste(exp, "_results.csv", sep = ""), row.names = F, fileEncoding = "UTF-8" )
