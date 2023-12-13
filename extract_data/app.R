# Load required libraries
library(shiny)
library(readxl)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(minpack.lm)

source("simple_function_ab_estimation_ar.R")
source("site_index_functions.R")
source("process_diameter_height_functions.R")
source("site_index_functions.R")


# ... (other libraries used in your script)

# Define UI
ui <- fluidPage(
  titlePanel("Extract summary data "),
 
  
  sidebarLayout(
    sidebarPanel(
      p("Tillgängliga försök i databas (trädslag/nummer/namn):"),
      
      uiOutput("lista_med_fsk"),
      
      # Input for experiment ID
      textInput("exp_id", "Write the Experiment ID ex:", value = "S1325"),
      
      # Button to trigger analysis
      actionButton("run_analysis", "Run Analysis"),
      
      #tags$img(src = "skogforsk logo.jfif"),
      
      
      # Text about version
      
      
      downloadButton("downloadCSV", "Download CSV"),
      p(),
      p("Version: 1.0"),
      
      # Text about developer
      p("Developer: Mateusz Liziniewicz"),
      
    ),
    
    mainPanel(
      # Output for displaying results
      textOutput("results"), 
      
      tableOutput("pri_table"),
      
      plotOutput("my_plot")
    )
  )
)

# Define server
server <- function(input, output) {
  
 
  
  output$lista_med_fsk <- renderUI({
    lista_fsk <- read_excel("DB_S1325_SödraVi.xlsx", sheet = "Försökmeta", na = ".") %>%
      mutate(fsk = paste(Trsl, Fsk_nummer, Name, sep = " / ")) %>%
      select(fsk) %>%
      distinct()
    # Replace this with your actual list or dynamic list generation logic
    my_list <- c(lista_fsk$fsk)
    
    
    # Render the list as an HTML list
    tags$ul(
      lapply(my_list, function(item) {
        tags$li(item)
      })
    )
  })
  
  
  
  observeEvent(input$run_analysis, {
    # Read experiment ID from input
    exp <- input$exp_id
    
    lista_fsk <- read_excel("DB_S1325_SödraVi.xlsx", sheet = "Försökmeta", na = ".") %>%
      mutate(fsk = paste(Fsk_nummer, Name, sep = " / ")) %>%
      select(fsk) %>%
      distinct()
      
    print(lista_fsk)
   
    year1 <- read_excel("DB_S1325_SödraVi.xlsx", sheet = "Försökmeta", na = ".") %>%
      filter(Exp == exp)
    
    year <- year1$Planteringsar
    experiment_name <-paste(year1$Fsk_nummer, " / ", year1$Name, sep = "")
    print(cat("Analyseras: ", experiment_name))
    
    library(readxl)
    library(tidyverse)
    df <- read_excel("DB_S1325_SödraVi.xlsx", sheet = "Data", na = ".") %>%
      filter(Exp == exp) %>% mutate(AGE = AR - year + 1)
    
    
    
    #####judge species######################################
    get_most_common_value_numeric <- function(data_frame, variable_name, result_variable_name) {
      # Extract the specified variable
      selected_variable <- data_frame[[variable_name]]
      
      # Use table() to get frequency of each unique value
      frequency_table <- table(selected_variable)
      
      # Find the most common value and convert to numeric
      most_common_value_numeric <- as.numeric(names(frequency_table)[which.max(frequency_table)])
      
      # Assign the result to a variable in the global environment
      assign(result_variable_name, most_common_value_numeric, envir = .GlobalEnv)
      return(cat("Dominant species is: ", most_common_value_numeric))
      
    }
    
    # Example usage with a data frame named 'your_data_frame' and a variable named 'your_variable'
    
    get_most_common_value_numeric(df, "TRSL", "species")
    #########################################################################################
    
    
    
    gal <- read_excel("DB_S1325_SödraVi.xlsx", sheet = "Försök", na = ".") %>%
      filter(EXP == exp) %>%
      select(YTA, AR, UTF) %>% mutate(AGE = AR - year + 1) %>% select(-AR)
    
    beh <- read_excel("DB_S1325_SödraVi.xlsx", sheet = "Försök", na = ".") %>% 
      filter(EXP == exp) %>%
      select(YTA, BEH) %>%
      unique()  
    
    
    area <- read_excel("DB_S1325_SödraVi.xlsx", sheet = "Försök", na = ".") %>% 
      filter(EXP %in% exp) %>%
      select(YTA, AREA) %>% 
      unique() %>%
      mutate(AREA = 10000/AREA)  
    
    
  
    source("simple_function_ab_estimation_ar.R")
    df_total <- fit_simple_model(df)
    
    
    
    
    
    df1 <- left_join(df, df_total, by = c("AGE"))
    df2 <- left_join(df1, beh, by = c("YTA")) %>%
      mutate(hest = a*D^b) %>% select(-a, -b)
    
    
    
    
    vol_calc <- function(D, H, TRSL = 10, na.rm = FALSE) {
      
      if (na.rm) {
        D <- na.omit(D)
        H <- na.omit(H)
      }
      
      vol <- ifelse(TRSL == 30,
                    (10^(-0.89363) * D^2.23818 * (D + 20.0)^(-1.06930) * H^6.02015 * (H - 1.3)^(-4.51472))/1000,
                    ifelse(TRSL %in% c(20, 21),
                           10^(-1.02039) * D^2.00128 * (D + 20)^(-0.47473) * H^2.87128 * (H - 1.3)^(-1.61083)/1000,
                           NA))
      
      return(vol)
    }
    
    basalA <- function(D, na.rm = FALSE){ #basal area function
      if(na.rm) {
        D <- na.omit(D)
      }
      
      BA <- pi*(D/2000)^2
    }
    
    df2 <- df2 %>% mutate(vol = vol_calc(D/10, hest/10, TRSL = TRSL), ba = basalA(D))
    
    
    df2 <- left_join(df2, area, by = "YTA")
    
    
    ####volume calculations
    
    
    # kvar <- df2 %>% group_by(YTA, BEH,  AGE) %>% filter(!is.na(D) | MORT3 == 1) %>%
    #   summarise(sumvol = sum(vol, na.rm = T)*mean(AREA), sumba = sum(ba, na.rm = T)*mean(AREA)) 
    
    
    # 
    # utgall <- df2 %>% group_by(YTA, BEH, AGE) %>% filter(MORT3 == 0) %>%
    #   summarise(sumvolGall = sum(vol, na.rm = T)*mean(AREA), sumbaGall = sum(ba, na.rm = T)*mean(AREA))
    # 
    # 
    # dod <- df2 %>% filter(MORT3 == 5) %>% group_by(YTA, BEH, AGE)  %>% 
    #   summarise(sumvolDod = sum(vol, na.rm = T)*mean(AREA), sumbaDod = sum(ba, na.rm = T)*mean(AREA))
    # 
    # 
    # df2 %>% group_by(YTA, BEH,  AGE, MORT3) %>% filter(!is.na(D) ) %>%
    #   summarise(sumvol = sum(vol, na.rm = T)*mean(AREA), sumba = sum(ba, na.rm = T)*mean(AREA)) 
    # 
    # 
    # process_volume_data(df2)
    # 
    # 
    # volume <- left_join(kvar, utgall, by = c("YTA", "BEH", "AGE")) %>% ungroup() %>%
    #   mutate(sumvolGall  = ifelse(is.na(sumvolGall ), 0, sumvolGall )) %>%
    #   mutate(kvarvol = sumvol - sumvolGall) %>% group_by(YTA, BEH) %>%
    #   mutate(utgall = cumsum(sumvolGall)) %>%
    #   mutate(totvol = kvarvol + utgall) %>%
    #   filter(!is.na(totvol) & sumvol != 0)
    
   
    
    volume <- process_volume_data(df2) 
    
    
    
    #to result table##########################
    res_volume <- volume %>% mutate(sumvol = kvarvol + as.numeric(dodvol) + utgall) %>% 
      group_by(YTA, BEH, AGE) %>% mutate(cumutgall = cumsum(utgall) + cumsum(dodvol), Totprod = sumvol + cumutgall) %>%
      select(YTA, BEH, AGE, kvarvol , utgall, dodvol, sumvol, cumutgall, Totprod) %>%
      rename(volfg = sumvol, volut = utgall, voleg = kvarvol, volut_cum = cumutgall) %>% mutate(MAI = Totprod/AGE) %>%
      select(YTA, BEH, AGE, volfg, volut, dodvol, voleg, volut_cum, Totprod, MAI)
    
    ##############volym figure###########################
    volume
    volume <- volume[rep(seq_len(nrow(volume)), each = 2), ]
    
    for (i in seq(1, length(volume$AGE), 2)){
      volume$kvarvol[i] = volume$kvarvol[i] + volume$utgall[i] + volume$dodvol[i]
    }
    
    ggplot(aes(x = AGE, y = kvarvol), data = volume) + geom_line() + geom_point() + facet_wrap(YTA~BEH )
    
    
    
    #####################################################################################
    
    #### top height #############
    
    
    htoh <- df2 %>% group_by(YTA, AGE) %>% top_n(5, wt = D) %>%
      summarise(HtOH = mean(hest, na.rm = T))
    
    
    
    #############basal area##############################
  
    
    basalarea <- process_ba_data(df2)
    
    
    res_ba <- basalarea %>% mutate(sumba = kvarba + as.numeric(dodba) + utgall_ba) %>% 
      group_by(YTA, BEH, AGE) %>% mutate(cumutgall_ba = cumsum(utgall_ba) + cumsum(dodba), Totprod_ba = sumba + cumutgall_ba) %>%
      select(YTA, BEH, AGE, kvarba , utgall_ba, dodba, sumba, cumutgall_ba, Totprod_ba) %>%
      rename(bafg = sumba, baut = utgall_ba, baeg = kvarba, baut_cum = cumutgall_ba) %>%
      select(YTA, BEH, AGE, bafg, baut, dodba, baeg, baut_cum, Totprod_ba)
    
    
    basalarea <- basalarea[rep(seq_len(nrow(basalarea)), each = 2), ]
    
  
    
    ##BA figure##################
    for (i in seq(1, length(basalarea$AGE), 2)){
      basalarea$kvarba[i] = basalarea$kvar[i] + basalarea$utgall_ba [i] + basalarea$dodba[i]
    }
    basalarea
    ggplot(aes(x = AGE, y = kvarba), data = basalarea) + geom_line() + geom_point() + facet_wrap(YTA~BEH )
    
    ggplot(aes(x = AGE, y = kvarba, group = YTA, color = factor(BEH)), data = basalarea) + geom_line() + geom_point() + facet_wrap(~BEH )
    
    
    
    
    
    
    
    
    ########## Antal ################
    
    
    
    
    
    
    
    # Example usage:
    antal1 <- process_antal_data(df2)
    
    res_antal <- antal1
    antal1 <- antal1[rep(seq_len(nrow(antal1)), each = 2), ]
    antal1[is.na(antal1)] <- 0
    
    
    
    
    for (i in seq(1, length(antal1$AGE), 2)) {
      antal1$neg[i] = antal1$neg[i] + antal1$nut[i] + antal1$dod_n[i] }
    
    
    ggplot(aes(x = AGE, y = neg), data = antal1) + geom_line() + geom_point() + facet_wrap(YTA~BEH )
    
    
    #################################
    
    
    ############################################################
    ### results diameter_height##############################
    ##########################################################
    #before thinning
    res_height_diameter_fg <- df2 %>% group_by(YTA, BEH, AGE) %>%
      summarise(mh_fg = mean(hest, na.rm = T), md_fg = mean(D, na.rm = T)) 
    
    
    #thinned = 0 and unthinnded = 1
    res_height_diameter<- df2 %>% group_by(YTA, BEH, AGE, MORT3) %>%
      summarise(mh = round(mean(hest, na.rm = T), 0), md = round(mean(D, na.rm = T), 0)) %>%
      filter(!is.na(MORT3))
    
    
    
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
      
      hh[is.na(hh)] <- "0"
      
      return(hh)
    }
    
    hh <- process_height_data(df2)
    
    #diameter
    process_diameter_data <- function(df) {
      
      
      
      
      #thinned = 0 and unthinnded = 1
      res_height_diameter<- df %>% group_by(YTA, BEH, AGE, MORT3) %>%
        summarise(mh = round(mean(hest, na.rm = T), 0), md = round(mean(D, na.rm = T), 0)) %>%
        filter(!is.na(MORT3))
      
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
      
      hh[is.na(hh)] <- "0"
      
      return(hh)
    }
    
    
    dd <- process_diameter_data(df2)
    
    
    rm(pri)
    #merging results
    hd <- left_join(res_height_diameter_fg, hh, by = c("YTA", "BEH", "AGE"))
    hd <- left_join(hd, dd, by = c("YTA", "BEH", "AGE")) %>% 
      select("YTA", "BEH", "AGE", "mh_fg", "mh_ut", "mh_eg", "mh_dod", "md_fg", "md_ut", "md_eg", "md_dod")
    
    
    pri <- left_join(res_antal, res_volume, by = c("YTA", "BEH", "AGE"))
    pri
    pri <- left_join(pri, hd, by = c("YTA", "BEH", "AGE"))
    
    
    pri <- left_join(pri, res_ba, by = c("YTA", "BEH", "AGE")) %>% mutate(GALL = ifelse(volut > 1, 1, 0)) %>%
      select(YTA, BEH, AGE, GALL, everything()) %>% 
      left_join(., htoh, by = c("YTA", "AGE")) %>%
      mutate(SI = calculate_si(HtOH/10, AGE, species))
    
    
    pri
    write.csv(pri %>% mutate_if(is.numeric, list(~format(., nsmall = 1))), paste(exp, "_results.csv", sep = ""), row.names = F, fileEncoding = "UTF-8" )
    
    
    output$pri_table <- renderTable({
      pri %>% mutate(AGE = factor(AGE), YTA = factor(YTA), GALL = factor(GALL)) %>% mutate(across(where(is.numeric), ~round(., digits = 1)))
    })
    
    
    
    
    #table to save as a result
    output$results <- renderText({
      write.csv(
        pri %>% mutate_if(is.numeric, list(~format(., nsmall = 1))),
        paste(exp, "_results.csv", sep = ""),
        row.names = FALSE,
        fileEncoding = "UTF-8"
      )
      print(paste(experiment_name, " / Analysis completed. Results saved as CSV.", sep = ""))
    })
    
    
    output$my_plot <- renderPlot({
      ggplot(aes(x = AGE, y = neg), data = antal1) + geom_line() + facet_wrap(YTA~BEH ) +
        theme_bw() +
        ggtitle("Number of stems over time")})

      
    output$downloadCSV <- downloadHandler(
        filename = function() {
          paste(input$exp_id, "_results.csv", sep = "")
        },
        content = function(file) {
          write.csv(pri, file, row.names = FALSE)
        }
      )      
    
    
  })
}

# Run the app
shinyApp(ui, server)
