#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1_housekeeping
# Calum Purdie
# 03/11/2021
# Sets up initial joined SMR01/SMR06 extract for analysis
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
adm_end <- "2021-12-31"
cancer_start <- "2015-01-01"
cancer_end <- "2021-12-31"

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
    geom_line(aes(y = obs_inc, colour = !!as.name(line_data)), size = 2) + 
    geom_line(aes(y = exp_inc), size = 0.5, linetype = "dashed") + 
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



### 3 Data Extraction ----

# Create SMR06 query
# This selects relevant variables from SMR06 between the specified cancer dates

smr06_query <-
  paste(
    "SELECT INCIDENCE_DATE, SEX, UPI_NUMBER, ICD10S_CANCER_SITE, ", 
    "AGE_IN_YEARS, POSTCODE, METHOD_1ST_DETECTION, ETHNIC_GROUP, TUMOUR_NO, ",
    "ENCR_INCIDENCE_DATE, OUT_OF_SCOTLAND ", 
    "FROM ANALYSIS.SMR06_PI",
    "WHERE INCIDENCE_DATE >= TO_DATE('", cancer_start, "', 'yyyy-mm-dd') AND",
    "INCIDENCE_DATE <= TO_DATE('", cancer_end, "', 'yyyy-mm-dd')"
  )
 
# Extract data from SMR06 using above query
# Clean names and remove records with blank upi_number
# Remove records where out_of_scotland is not blank
# Create 3 character cancer_site code based on icd10s_cancer_site
# Define incidence type, format date and calculate incidence year
# Remove men with breast or cervical cancer and women with prostate cancer

smr06_data <- as_tibble(dbGetQuery(channel, statement = smr06_query)) |> 
  clean_names() |> 
  filter(!is.na(upi_number)) |>  
  filter(is.na(out_of_scotland)) |> 
  mutate(cancer_site = substr(icd10s_cancer_site, 1, 3), 
         incidence_type = case_when(cancer_site %in% c("C18", "C19", "C20") ~ 
                                      "Colorectal", 
                                    cancer_site %in% c("C33", "C34") ~ "Lung",
                                    cancer_site %in% c("C50") ~ "Breast", 
                                    cancer_site %in% c("C61") ~ "Prostate", 
                                    cancer_site %in% c("C00", "C01", "C02", 
                                                       "C03", "C04", "C05", 
                                                       "C06", "C07", "C08", 
                                                       "C09", "C10", "C11", 
                                                       "C12", "C13", "C14", 
                                                       "C30", "C31", "C32") ~ 
                                      "Head and Neck", 
                                    cancer_site %in% c("C53") ~ "Cervical", 
                                    cancer_site %in% c("C15", "C16", "C17", 
                                                       "C22", "C23", "C24", 
                                                       "C25", "C26") ~ "Upper GI"),
         old_incidence_date = as.Date(incidence_date),
         encr_incidence_date = as.Date(encr_incidence_date),
         incidence_date = case_when(old_incidence_date < "2019-01-01" ~ old_incidence_date,
                                    encr_incidence_date >= "2019-01-01" ~ encr_incidence_date
         )) |>  
  mutate(incidence_year = year(incidence_date)) |>  
  filter(between(incidence_year, as.numeric(start), as.numeric(end))) |>  
  filter(!(incidence_type == "Breast" & sex == "1") & 
         !(incidence_type == "Prostate" & sex == "2") & 
         !(incidence_type == "Cervical" & sex == "1"))

# Create SMR01 query for emergency admissions
# This selects variables from SMR01 between the specified admission dates

query_admissions <-
  paste(
    "SELECT UPI_NUMBER, ADMISSION_DATE, HBTREAT_CURRENTDATE, DOB,
    MAIN_CONDITION, OTHER_CONDITION_1, OTHER_CONDITION_2,
    OTHER_CONDITION_3, OTHER_CONDITION_4, OTHER_CONDITION_5,
    CIS_MARKER, DISCHARGE_DATE, MANAGEMENT_OF_PATIENT, LENGTH_OF_STAY, 
    ADMISSION_TYPE",
    "FROM ANALYSIS.SMR01_PI",
    "WHERE ADMISSION_DATE >= TO_DATE('", adm_start, "', 'yyyy-mm-dd') AND ",
    "ADMISSION_DATE <= TO_DATE('", adm_end, "', 'yyyy-mm-dd')"
  )

# Extract data from SMR01 using above query
# Clean names and filter for UPIs in SMR06 extract
# Define emergency and elective admissions

smr01_adm <- as_tibble(dbGetQuery(channel, statement = query_admissions)) |> 
  clean_names() |>  
  filter(upi_number %in% smr06_data$upi_number) |>  
  mutate(admission = case_when(admission_type %in% c("20", "21", "22", "30", 
                                                     "31", "32", "33", "34", 
                                                     "35", "36", "38", 
                                                     "39") ~ "emergency", 
                               admission_type %in% c("10", "11", "12", "18", 
                                                     "19") ~ "elective"))

# Free unused memory

gc()

# Read in postcode simd lookup file
# Select columns and recode lower and upper quintile values

pc_simd <- readRDS(glue("/conf/linkage/output/lookups/Unicode/Deprivation/", 
                        "postcode_2023_2_simd2020v2.rds")) |>  
  select(postcode = pc7, simd2020v2_sc_quintile) |>  
  mutate(simd2020v2_sc_quintile = as.character(simd2020v2_sc_quintile)) |>  
  mutate(simd2020v2_sc_quintile = recode(simd2020v2_sc_quintile, 
                                         "1" = "1 (most deprived)", 
                                         "5" = "5 (least deprived)"))

# Read in data zone population for incidence and select columns
# Group data by year, sex and simd2020v2_sc_quintile and sum population by age
# Filter for data between specified years
# Pivot data to longer format, with age group columns all falling into one
# larger age_group column and their values going to pop
# Group data by year, sex and simd2020v2_sc_quintile and recalculate age_group
# as the row number - 1, this gives us the 0-18 groupings used elsewhere
# Recode sex to numeric for matching and add most/least deprived to SIMD names

simd_pop <- readRDS(glue("/conf/linkage/output/lookups/Unicode/Populations/", 
                         "Estimates/", 
                         "DataZone2011_pop_est_5year_agegroups_2011_2021.rds")) |> 
  select(year, datazone2011, sex, starts_with("ageg"), 
         simd2020v2_sc_quintile) |>  
  group_by(year, sex, simd2020v2_sc_quintile) |>  
  summarise(across(starts_with("ageg"), ~sum(.x))) |> 
  filter(between(year, as.numeric(start), as.numeric(end))) |>  
  pivot_longer(cols = starts_with("ageg"), names_to = "age_group", 
               values_to = "pop") |>  
  group_by(year, sex, simd2020v2_sc_quintile) |>  
  mutate(age_group = row_number() - 1) |>  
  ungroup() |>  
  mutate(sex = case_when(sex %in% c("m", "M") ~ "1", 
                         sex %in% c("f", "F") ~ "2"), 
         simd2020v2_sc_quintile = recode(as.character(simd2020v2_sc_quintile), 
                                         "1" = "1 (most deprived)", 
                                         "5" = "5 (least deprived)"))

# Calculate the Scotland level population
# Group by year, sex and age_group and summarise pop within each group

scot_pop <- simd_pop |>  
  group_by(year, sex, age_group) |> 
  summarise(pop = sum(pop))

# Read in European standard populations by sex and drop SPSS formatting
# Align age group numbering with other data and select columns required

esp2013 <- read_sav(glue("/conf/linkage/output/lookups/Unicode/Populations/", 
                         "Standard/ESP2013_by_sex.sav")) |> 
  zap_formats() |> 
  zap_widths() |> 
  remove_all_labels() |>  
  clean_names() |> 
  mutate(age_group = agegroup - 1, 
         sex = as.character(sex)) |> 
  select(sex, age_group, esp2013)



### 4 Data Manipulation ----

# Join data and calculate time between admission and incidence
# Join on SIMD data
# Define age groups and sex name

joined_data <- smr06_data |>  
  left_join(smr01_adm) |>  
  mutate(time_diff = time_length(admission_date %--% incidence_date, "days")) |>  
  left_join(pc_simd) |>  
  mutate(age_group = create_age_groups(age_in_years, 0, 90, 5)) |>  
  mutate(age_group = if_else(age_in_years >= 0 & age_in_years <= 39, 
                             "0-39", age_group)) |>  
  mutate(sex_name = recode(sex, "1" = "Male", "2" = "Female"))

# Arrange data by upi, admission date and admission type
# This orders the data by upi, date and then takes emergency over elective
# We take emergency over elective for admissions where someone is admitted as
# an emergency and then transferred as an elective on the same day
# Group by upi and cis marker
# Calculate the admission type for the first stay of that cis
# Flag people who had an emergency admission within 30 days prior to diagnosis
# This is any case where the first admission type is emergency and the admission
# is within 30 days, or any emergency admission within 30 days

joined_data <- joined_data |>  
  arrange(upi_number, admission_date, desc(admission)) |>  
  group_by(upi_number, cis_marker) |>  
  mutate(first_admission_type = first(admission)) |>  
  ungroup() |>  
  mutate(emergency_flag = case_when(time_diff >= 0 & time_diff <= 30
                                    & first_admission_type == "emergency" ~ "Emergency",
                                    time_diff >= 0 & time_diff <= 30 & 
                                      admission == "emergency" ~ "Emergency", 
                                    TRUE ~ "Non-Emergency"))

# Arrange by tumour_no, UPI, emergency_flag, admission and incidence date
# This orders the data by tumour_no and then upi, prioritises records with yes 
# for emergency_flag, prioritises the most recent admission and then sorts by 
# incidence_date
# Group data by tumour_no and take first row for each
# This means we only keep one row per person and cancer

joined_data <- joined_data |>  
  arrange(tumour_no, upi_number, emergency_flag, desc(admission), incidence_date) |>  
  group_by(tumour_no) |>  
  slice(1) |>  
  ungroup()

# Tidy environment

rm(smr01_adm, smr06_data, pc_simd)
gc()

# saveRDS(joined_data, "Temp/20230923_data.rds")
