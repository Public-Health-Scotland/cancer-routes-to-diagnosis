#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1_housekeeping
# Calum Purdie
# 03/11/2021
# Housekeeping script for use in later scripts
# Written/run on Posit Workbench
# R version 4.1.2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### 1 Housekeeping ----

# Load packages

library(odbc)
library(dplyr)
library(haven)
library(janitor)
library(tidyr)
library(lubridate)
library(ggplot2)
library(here)
library(purrr)
library(glue)
library(phsverse)
library(stringr)
library(openxlsx)
library(sjlabelled)
library(PHEindicatormethods)
library(forcats)
library(ggrepel)
library(tidylog)

# Connect to SMRA tables using odbc connection
# The suppressWarnings function prevents your password from appearing in the
# console if the connection is unsuccessful

channel <- suppressWarnings(dbConnect(
  odbc(),
  dsn = "SMRA",
  uid = .rs.askForPassword("What is your user ID?"),
  pwd = .rs.askForPassword("What is your LDAP password?")))

# Define dates
# These are the start and end dates for admissions and cancer diagnoses
# Admission start date should be at least 30 days before cancer start date to
# allow for emergencies 30 days before admission to be identified

adm_start <- "2014-12-01"
adm_end <- "2022-12-31"
cancer_start <- "2015-01-01"
cancer_end <- "2022-12-31"

# Define years from start and end dates

start <- str_sub(cancer_start, 1, 4)
end <- str_sub(cancer_end, 1, 4)

model_start <- 2015
model_end <- 2019

# Stop scientific notation for small numbers

options(scipen = 999)



### 2 Functions ----

# Define function for calculating age groups for standard populations

standard_pop_age_groups <- function(current_col){
  
  # Define an age group for each current_col value
  
  case_when(current_col >= 0 & current_col <= 4 ~ 0, 
            current_col >= 5 & current_col <= 9 ~ 1, 
            current_col >= 10 & current_col <= 14 ~ 2, 
            current_col >= 15 & current_col <= 19 ~ 3, 
            current_col >= 20 & current_col <= 24 ~ 4, 
            current_col >= 25 & current_col <= 29 ~ 5, 
            current_col >= 30 & current_col <= 34 ~ 6, 
            current_col >= 35 & current_col <= 39 ~ 7, 
            current_col >= 40 & current_col <= 44 ~ 8, 
            current_col >= 45 & current_col <= 49 ~ 9, 
            current_col >= 50 & current_col <= 54 ~ 10, 
            current_col >= 55 & current_col <= 59 ~ 11, 
            current_col >= 60 & current_col <= 64 ~ 12, 
            current_col >= 65 & current_col <= 69 ~ 13, 
            current_col >= 70 & current_col <= 74 ~ 14, 
            current_col >= 75 & current_col <= 79 ~ 15, 
            current_col >= 80 & current_col <= 84 ~ 16, 
            current_col >= 85 & current_col <= 89 ~ 17, 
            current_col >= 90 ~ 18)
  
}

# Define function for creating fitted poisson regressions
# Separate this into two functions, one including sex in model and one without

fitted_model_without_sex <- function(df, ...){
  
  # Filter for years within start and end
  # Group by specified columns and nest data
  # Add a column for regression model for each row by using map()
  # Include year, age_group and offset population within model and use a poisson
  # regression model
  
  # I tried using nest_by() rather than group_by() and nest() but this didn't
  # work when mapping the model
  
  new_df <- df |>  
    filter(year >= model_start & year <= model_end) |>  
    group_by(...) |> 
    nest() |> 
    mutate(mod = map(data, ~ glm(n ~ year + as.factor(age_group) + 
                                   offset(log(pop)), 
                                 data = ., family = poisson(link = "log"))
                     )
           )
  
}

fitted_model_with_sex <- function(df, ...){
  
  # Filter for years within start and end
  # Group by specified columns and nest data
  # Add a column for regression model for each row by using map()
  # Include year, age_group, sex and offset population within model and use a 
  # poisson regression model
  
  new_df <- df |>  
    filter(year >= model_start & year <= model_end) |>  
    group_by(...) |> 
    nest() |> 
    mutate(mod = map(data, ~ glm(n ~ year + as.factor(age_group) + 
                                   as.factor(sex) + offset(log(pop)), 
                                 data = ., family = poisson(link = "log"))
                     )
           )
  
} 

# Define function to calculate incidence rates and standardised incidence ratios

calculate_incidence <- function(df, grouping_variables, standard_pop,
                                events, population, confidence_level){
  
  # Group across all grouping variables
  # Calculate the sum of the standard population squared multiplied by number
  # of events, divided by actual population squared
  # This gives us the numerator for our square root
  # Calculate the standard incidence as the sum of the standard population
  # multiplied by number of events, divided by actual population squared and
  # then divide this number by the sum of the standard population
  # Calculate the sum of events and the squared sum of the standard population
  # Define lower and upper chi-square quantile distributions
  # Define the lower and upper multiplicands based on their inverse chi-square
  # divided by 2 minus the sum of events
  
  # We can then summarise this data by taking n as the sum of events and inc as 
  # the standard incidence multiplied by 100,000
  # We can then derive the lower and upper confidence intervals using the 
  # standardised incidenc added to the square root value multiplied by the 
  # multiplicand, which is then multiplied by 100,000
  
  df_new <- df |>  
    group_by(across(all_of(grouping_variables))) |> 
    mutate(ci_sqrt_num = sum(({{standard_pop}} ^ 2 * {{events}})/({{population}} ^ 2)),
           inc_std = sum({{standard_pop}} * {{events}} / {{population}}) / sum({{standard_pop}}),
           sum_events = sum({{events}}),
           sum_standard_pop_all_squared = sum({{standard_pop}}) ^ 2,
           ci_sqrt = sqrt(ci_sqrt_num / (sum_events *
                                           sum_standard_pop_all_squared)),
           l_ci_inv_chi_sq = qchisq((1 - confidence_level) / 2, 2 * sum_events),
           u_ci_inv_chi_sq = qchisq(confidence_level + (1 - confidence_level)/2, 
                                    (2 * sum_events) + 2),
           l_ci_multiplicand = l_ci_inv_chi_sq / 2 - sum_events,
           u_ci_multiplicand = u_ci_inv_chi_sq / 2 - sum_events) |> 
    summarise(n = sum({{events}}), 
              inc = inc_std * 100000, 
              l_ci = (inc_std + (ci_sqrt * l_ci_multiplicand)) * 100000,
              u_ci = (inc_std + (ci_sqrt * u_ci_multiplicand)) * 100000) |> 
    distinct()
  
}

calculate_ratio <- function(df, grouping_variables, observed, expected, 
                            confidence_level, type){
  
  # Group across all grouping variables
  
  df_grouped <- df |>  
    group_by(across(all_of(grouping_variables)))
  
  # Define the ratio name based on the input type
  
  ratio_name <- deparse(substitute(type))
  
  # Use phe_rate to define the rate ratios
  # Input the data, observed values, expected values and confidence level
  # Multiplier can just be 1 here as our data is already per 100,000
  # Drop unnecessary columns and rename ratios
  
  ratio_data <- phe_rate(data = df_grouped, x = {{observed}}, n = {{expected}}, 
                         confidence = confidence_level, multiplier = 1) |>  
    select(-c(confidence, statistic, method)) |>  
    rename(ratio = value, 
           ratio_l_ci = lowercl, 
           ratio_u_ci = uppercl) |>  
    rename_with(~ tolower(gsub("ratio", ratio_name, .x, fixed = TRUE)))
  
  # Define the output by joining the initial grouped data with the ratio data
  
  output <- df_grouped |>  
    full_join(ratio_data) |>  
    ungroup()
  
}

# Define function for creating a chart 

create_chart <- function(df, presentation_type, line_data){
  
  # Filter for data where emergency_flag is presentation_type
  # Define year as character and set line_data variable as an ordered factor
  # Define ggplot aesthetics and add line using line_data
  # Add an error bar for 95% confidence intervals
  # Add theme details to adjust text - this will make the chart look odd here
  # but looks fine once saved out
  
  df |>  
    filter(emergency_flag == presentation_type) |>  
    mutate(year = as.character(year), 
           !!as.name(line_data) := fct_reorder(!!as.name(line_data), obs_inc, 
                                               tail, n = 1, .desc = TRUE)) |>        
    ggplot(aes(x = year, y = obs_inc, group = !!as.name(line_data))) +
    geom_line(aes(colour = !!as.name(line_data)), size = 2) + 
    geom_errorbar(aes(ymin = obs_l_ci, ymax = obs_u_ci), 
                  width = .1) +
    theme(plot.title = element_text(hjust = 0.5, size = 40),
          axis.text.x = element_text(vjust = 0.5, size = 30),
          axis.text.y = element_text(size = 30),
          axis.title = element_text(size = 35),
          legend.title = element_text(size = 30),
          legend.key.size = unit(1.5, "cm"), 
          legend.text = element_text(size = 30))
    
}

# Define function for creating a chart with regression line

create_regression_chart <- function(df, presentation_type, line_data){
  
  # Filter for data where emergency_flag is presentation_type
  # Define year as character and set line_data variable as an ordered factor
  # Define ggplot aesthetics and add line using line_data
  # Add a line for expected value from regression model
  # Add an error bar for 95% confidence intervals
  # Add theme details to adjust text - this will make the chart look odd here
  # but looks fine once saved out
  
  df |>  
    filter(emergency_flag == presentation_type) |>  
    mutate(year = as.character(year), 
           !!as.name(line_data) := fct_reorder(!!as.name(line_data), obs_inc, 
                                               tail, n = 1, .desc = TRUE)) |>     
    ggplot(aes(x = year, group = !!as.name(line_data), 
               colour = !!as.name(line_data))) +
    geom_line(aes(y = obs_inc, colour = !!as.name(line_data)), size = 2.5) + 
    geom_line(aes(y = exp_inc), size = 1, linetype = "dashed") + 
    geom_errorbar(aes(ymin = obs_l_ci, ymax = obs_u_ci), 
                  width = .1) +
    theme(plot.title = element_text(hjust = 0.5, size = 40),
          axis.text.x = element_text(vjust = 0.5, size = 30),
          axis.text.y = element_text(size = 30),
          axis.title = element_text(size = 35),
          legend.title = element_text(size = 30),
          legend.key.size = unit(1.5, "cm"), 
          legend.text = element_text(size = 30))
  
}

# Define functions for Wilson confidence intervals
# These are used for the simple proportion calculations

Wilson_lowerCI <- function(kpi_p, alpha, n){
  ((kpi_p + qchisq((1 - alpha), df = 1) / (2 * n) 
    - qnorm(1 - (alpha / 2) , mean = 0, sd = 1)
    * sqrt((kpi_p * (1 - kpi_p) + qchisq((1 - alpha), df = 1) 
            / (4 * n)) / n))) / (1 + qchisq((1 - alpha), df = 1) / n)   
}

Wilson_upperCI <- function(kpi_p, alpha, n){
  ((kpi_p + qchisq((1 - alpha), df = 1) / (2 * n) 
    + qnorm(1 - (alpha / 2), mean = 0, sd = 1)
    * sqrt((kpi_p * (1 - kpi_p) + qchisq((1 - alpha), df = 1) 
            / (4 * n)) / n))) / (1 + qchisq((1 - alpha), df = 1) / n)   
}
