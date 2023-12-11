#height
process_height_data <- function(df) {
  
  
  
  
  #thinned = 0 and unthinnded = 1
  res_height_diameter<- df2 %>% group_by(YTA, BEH, AGE, MORT3) %>%
    summarise(mh = round(mean(hest, na.rm = T), 0), md = round(mean(D, na.rm = T), 0)) %>%
    filter(!is.na(MORT3))
  
  #height
  hh <- reshape2::dcast(res_height_diameter, YTA + BEH + AGE ~ MORT3, value.var = "mh")
  
  names(hh)[names(hh)== "0"] <- "mh_ut" #removed in thinning
  names(hh)[names(hh)=="1"] <- "mh_eg" # left in stand'
  names(hh)[names(hh)=="5"] <- "mh_dod" # left in stand'
  
  
  if (!all(c("utgall_n", "n", "dod_n") %in% names(hh))) {
    if (!"mh_ut" %in% names(hh)) {
      hh <- hh %>% mutate(mh_ut = 0)
    }
    if (!"mh_dod" %in% names(hh)) {
      hh <- hh %>% mutate(mh_dod = 0)
    }
  }
  
  hh[is.na(hh)] <- 0
  
  return(hh)
}

#process_height_data(df2)
#diameter
process_diameter_data <- function(df) {
  
  
  
  
  #thinned = 0 and unthinnded = 1
  res_height_diameter<- df %>% group_by(YTA, BEH, AGE, MORT3) %>%
    summarise(mh = round(mean(hest, na.rm = T), 0), md = round(mean(D, na.rm = T), 0)) 
  #height
  hh <- reshape2::dcast(res_height_diameter, YTA + BEH + AGE ~ MORT3, value.var = "md")
  
  names(hh)[names(hh)== "0"] <- "md_ut" #removed in thinning
  names(hh)[names(hh)=="1"] <- "md_eg" # left in stand'
  names(hh)[names(hh)=="5"] <- "md_dod" # left in stand'
  
  
  if (!all(c("md_ut", "md_ut", "md_dod") %in% names(hh))) {
    if (!"md_ut" %in% names(hh)) {
      hh <- hh %>% mutate(md_ut = 0)
    }
    if (!"md_dod" %in% names(hh)) {
      hh <- hh %>% mutate(md_dod = 0)
    }
  }
  
  hh[is.na(hh)] <- 0
  
  return(hh)
}

#process_diameter_data(df2)
###volume
process_volume_data <- function(df) {
  
  
  
  
  #thinned = 0 and unthinnded = 1
  res_height_diameter<- df %>% group_by(YTA, BEH,  AGE, MORT3) %>% filter(!is.na(D) ) %>%
    summarise(sumvol = round(sum(vol, na.rm = T)*mean(AREA), 1), sumba = round(sum(ba, na.rm = T)*mean(AREA), 1)) %>%
    filter(!is.na(MORT3))
  
  #height
  hh <- reshape2::dcast(res_height_diameter, YTA + BEH + AGE ~ MORT3, value.var = "sumvol")
  
  names(hh)[names(hh)== "0"] <- "utgall" #removed in thinning
  names(hh)[names(hh)=="1"] <- "kvarvol" # left in stand'
  names(hh)[names(hh)=="5"] <- "dodvol" # left in stand'
  
  #hh <- hh %>% mutate(dodvol = numeric(dodvol))
  
  if (!all(c("kvarvol", "utgall", "dodvol") %in% names(hh))) {
    if (!"utgall" %in% names(hh)) {
      hh <- hh %>% mutate(utgall = 0)
    }
    if (!"dodvol" %in% names(hh)) {
      hh <- hh %>% mutate(dodvol = 0)
    }
  }
  
  hh[is.na(hh)] <- 0
  
  return(hh)
}
#process_volume_data(df2)

#basal area
process_ba_data <- function(df) {
  
  #thinned = 0 and unthinnded = 1
  res_height_diameter<- df %>% group_by(YTA, BEH,  AGE, MORT3) %>% filter(!is.na(D) ) %>%
    summarise(sumba = round(sum(ba, na.rm = T)*mean(AREA), 1)) %>%
    filter(!is.na(MORT3))
  
  #height
  hh <- reshape2::dcast(res_height_diameter, YTA + BEH + AGE ~ MORT3, value.var = "sumba")
  
  names(hh)[names(hh)== "0"] <- "utgall_ba" #removed in thinning
  names(hh)[names(hh)=="1"] <- "kvarba" # left in stand'
  names(hh)[names(hh)=="5"] <- "dodba" # left in stand'
  
  hh[is.na(hh)] <- 0
  if (!all(c("kvarba", "utgall_ba", "dodba") %in% names(hh))) {
    if (!"utgall_ba" %in% names(hh)) {
      hh <- hh %>% mutate(utgall_ba = 0)
    }
    if (!"dodba" %in% names(hh)) {
      hh <- hh %>% mutate(dodba = 0)
    }
  }

  
  return(hh)
}

#process_ba_data(df2)

#antal
process_antal_data <- function(df) {
  
  antal1 <- df %>% filter(!is.na(MORT3)) %>% group_by(YTA, BEH, AGE, MORT3) %>%
  summarise(li = round(sum(!is.na(MORT3)) * mean(AREA, na.rm = T), 0))
  
  
  antal1 <- reshape2::dcast(antal1, YTA + BEH + AGE ~ MORT3, value.var = "li")
  
  names(antal1)[names(antal1) == "0"] <- "utgall_n" # removed in thinning
  names(antal1)[names(antal1) == "1"] <- "n"       # left in stand
  names(antal1)[names(antal1) == "5"] <- "dod_n"   # natural mortality
  
  if (!all(c("utgall_n", "n", "dod_n") %in% names(antal1))) {
    if (!"utgall_n" %in% names(antal1)) {
      antal1 <- antal1 %>% mutate(utgall_n = 0)
    }
    if (!"dod_n" %in% names(antal1)) {
      antal1 <- antal1 %>% mutate(dod_n = 0)
    }
  }
  
  antal1[is.na(antal1)] <- 0
  
  res_antal <- antal1 %>% group_by(YTA, BEH) %>%
    mutate(nfg = n + utgall_n + dod_n, nut = utgall_n, neg = n) %>%
    select(-utgall_n, -n) %>%
    mutate(nutcum = cumsum(nut) + cumsum(dod_n))
  
  return(res_antal)
}
#process_antal_data(df2)
