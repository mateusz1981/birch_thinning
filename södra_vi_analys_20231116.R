setwd("C:/Users/mali/OneDrive - Skogforsk/BJÖRK/Björk_gallrig")

exp = "Test"


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

#to result table
res_volume <- volume %>% select(YTA, BEH, AR, sumvol, sumvolGall, kvarvol, utgall, totvol) %>%
  rename(volfg = sumvol, volut = sumvolGall, voleg = kvarvol, volut_cum = utgall, Totvol = totvol)


### results diameter_height
res_height_diameter <- df2 %>% group_by(YTA, BEH, AR, MORT3) %>%
  summarise(mh = mean(hest, na.rm = T), md = mean(D, na.rm = T)) %>%
  filter(!is.na(MORT3))

# height/diameter before, after treamtent, and diameter/height of thinned stems


library(reshape2)
height <- dcast(res_height_diameter, YTA + BEH + AR ~ MORT3, value.var = "mh")
diameter <- dcast(res_height_diameter, YTA + BEH + AR ~ MORT3, value.var = "md")



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
ggplot(aes(x = AR-1995+2, y = totvol/(AR-1995+2)), data = volume) + 
  geom_line(aes(group = YTA, color = BEH), linewidth = 1) +
  xlab("Ålder") + ylab("MAI (m3 ha år)") +
  theme_bw()




########## Antal rätt################
antal1 <- df2 %>% filter(!is.na(MORT3)) %>% group_by(YTA, BEH, AR, MORT3) %>%
  summarise(li = round(sum(!is.na(MORT3)) * mean(AREA, na.rm = T), 0))

antal1 <- reshape2::dcast(antal1, YTA + BEH + AR ~ MORT3, value.var = "li")

names(antal1)[names(antal1)== "0"] <- "utgall_n" #removed in thinning
names(antal1)[names(antal1)=="1"] <- "n" # left in stand
names(antal1)[names(antal1)== "5"] <- "dod_n" # natural mortality




antal1 <- antal1[rep(seq_len(nrow(antal1)), each = 2), ]
antal1[is.na(antal1)] <- 0




for (i in seq(1, length(antal1$AR), 2)) {
  if (exists("antal1$n") && exists("antal1$utgall_n") && exists("antal1$dod_n")) {
    antal1$n[i] = antal1$n[i] + antal1$utgall_n[i] + antal1$dod_n[i]
  } else {
    antal1$n[i] = antal1$n[i] + antal1$utgall_n[i]
}}
ggplot(aes(x = AR, y = n), data = antal1) + geom_line() + facet_wrap(YTA~BEH )


#################################


