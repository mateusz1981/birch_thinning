
calculate_si <- function(input_height, input_age, species_code) {
  co_si <- read.csv("coefficients_main.csv", header = TRUE, sep = ",")

  
  cal_si <- function(species_code) {
    da <- subset(co_si, co_si$TRSL == species_code)
    
    if (species_code != 20) {
      h = (input_height + da$beta / da$asi^da$b2 + sqrt((input_height - da$beta / da$asi^da$b2)^2 + (4 * da$beta * input_height) / input_age^da$b2)) / 
        (2 + 4 * da$beta / ((da$ref.age^da$b2) * (input_height - da$beta / da$asi^da$b2 + sqrt((input_height - da$beta / da$asi^da$b2)^2 + (4 * da$beta * input_height) / input_age^da$b2))))
      return(round(h, 2))
    } else if (species_code == 20) {
      hg = (input_height + da$beta / da$asi^da$b2 + sqrt((input_height - da$beta / da$asi^da$b2)^2 + (4 * da$beta * input_height) / (input_age - 3)^da$b2)) / 
        (2 + 4 * da$beta / (((da$ref.age - 3)^da$b2) * (input_height - da$beta / da$asi^da$b2 + sqrt((input_height - da$beta / da$asi^da$b2)^2 + (4 * da$beta * input_height) / (input_age - 3)^da$b2))))
      return(round(hg, 2))
    }
  }
  
  result <- cal_si(species_code)
  return(result)
}

calculate_si(20, 34, 10)



