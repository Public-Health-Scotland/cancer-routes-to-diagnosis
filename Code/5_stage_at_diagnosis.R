#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5_stage_at_diagnosis
# Calum Purdie
# 27/05/2022
# Calculates standardised incidence rates and ratios by stage at diagnosis
# Written/run on Posit Workbench
# R version 4.1.2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



### 1 Housekeeping ----

source(here::here("Code/1_housekeeping.R"))



### 2 Extract CAD data ----

# Read in CAD - only keep relevant columns
# Removing SPSS formatting and names
# Drop rows with blank tnm_stage_2 and define tnm_stage_2 as character

cad <- readRDS(glue("/PHI_conf/CancerGroup1/Topics/CancerStatistics/Data/",
                    "Registry/CAD/Cancer_Analysis_Database.rds")) |> 
  select(tumourid, tnm_stage_2, figo_stage_2, site10) |> 
  zap_formats() |>
  zap_widths() |>
  remove_all_labels() |>
  mutate(tnm_stage_2 = case_when(site10 == "C53" ~ figo_stage_2,
                                 TRUE ~ tnm_stage_2)) |>
  filter(tnm_stage_2 != "") |> 
  mutate(tnm_stage_2 = as.character(tnm_stage_2))

# Join cad data onto analysis data

joined_data <- joined_data |>
  left_join(cad, by = c("tumour_no" = "tumourid")) |>
  mutate(tnm_stage_2 = case_when(is.na(tnm_stage_2) ~ "X",
                                 TRUE ~ tnm_stage_2))



### 3 Overall Analysis ----

# Count data by incidence_year, emergency_flag and tnm_stage_2
# Group by incidence_year and emergency_flag and calculate percentages by group

stage_at_diagnosis <- joined_data |> 
  count(incidence_year, emergency_flag, tnm_stage_2) |> 
  group_by(incidence_year, emergency_flag) |> 
  mutate(pc = round_half_up(n * 100 / sum(n), 1)) |> 
  ungroup()



### 4 Site Level Analysis ----

# Count data by incidence_year, emergency_flag, incidence_type and tnm_stage_2
# Group by incidence_year and emergency_flag and calculate percentages by group

stage_at_diagnosis_site <- joined_data |> 
  count(incidence_year, emergency_flag, incidence_type, tnm_stage_2) |> 
  group_by(incidence_year, emergency_flag, incidence_type) |> 
  mutate(pc = round_half_up(n * 100 / sum(n), 1)) |> 
  ungroup()



### 5 Regression Models and Incidence Rates ----

# Get a list of all cancer types in analysis

stages <- joined_data |> 
  distinct(tnm_stage_2) |> 
  pull()

# Create a data frame for matching to joined_data
# This is every combination of year, sex, age group, flag and incidence_type
# This ensures all populations are counted for the incidence rates, even where
# there are no admissions for a specific category, e.g. if there are no elective
# admissions for men in age group 0 in 2019, joining on this data frame will 
# create a row for this group

match_df <- expand.grid(year = c(start:end), sex = c("1", "2"), 
                        age_group = c(0:18), 
                        emergency_flag = c("Non-Emergency", "Emergency"), 
                        tnm_stage_2 = c(stages))

# Define age groups for admission data
# Count data to get total number of admissions by year, sex, age_group,  
# emergency_flag and incidence_type
# Join on match_df_site to add any groups with no admissions and set number of
# admissions for these groups to 0, as they appear as NAs when joined on
# Join on Scotland populations and also ESP and arrange data

inc_rate_data <- joined_data |> 
  mutate(age_group = standard_pop_age_groups(age_in_years)) |> 
  count(year = incidence_year, sex, age_group, emergency_flag, tnm_stage_2) |> 
  full_join(match_df) |> 
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .))) |> 
  left_join(scot_pop) |> 
  left_join(esp2013) |> 
  arrange(year, sex, age_group, emergency_flag, tnm_stage_2)

# Define fitted model for all cancers
# Group data by emergency_flag, incidence_type and tnm_stage_2
# Drop data column

reg_mod <- fitted_model_with_sex(inc_rate_data, emergency_flag, tnm_stage_2) |> 
  select(-data)

# Group by emergency_flag and tnm_stage_2 and nest data
# Do an inner join with site_sex_reg_mod
# Create column for expected values by predicting the model in the mod column
# using the data from the data column
# Unnest data to get data back into standard format

reg_mod_output <- inc_rate_data |> 
  group_by(emergency_flag, tnm_stage_2) |>
  nest() |> 
  inner_join(reg_mod) |> 
  mutate(exp = purrr::map2(mod, data, predict, type = "response")) |> 
  unnest(c(data, exp))

# Calculate incidence rates by year, emergency_flag and tnm_stage_2 for
# observed and expected values
# Rename columns to indicate observed or expected

obs_inc_output <- calculate_incidence(reg_mod_output, 
                                      c("year", "emergency_flag", 
                                        "tnm_stage_2"),
                                      esp2013, n, pop, 0.95) |> 
  rename_with(~ paste0("obs_", .), c("n", "inc", "l_ci", "u_ci"))

exp_inc_output <- calculate_incidence(reg_mod_output, 
                                      c("year", "emergency_flag", 
                                        "tnm_stage_2"),
                                      esp2013, exp, pop, 0.95) |> 
  rename_with(~ paste0("exp_", .), c("n", "inc", "l_ci", "u_ci"))

# Join observed and expected together

inc_output <- full_join(obs_inc_output, exp_inc_output)

# Calculate standardised incidence ratio and standardised rate ratio

inc_output_sir <- calculate_ratio(inc_output, 
                                  c("year", "emergency_flag", 
                                    "tnm_stage_2"), 
                                  obs_n, exp_n, 0.95, sir)

# Calculate p-values for SIRs
# This only works for intergers, so use obs_n and exp_n

inc_output_sir <- inc_output_sir |> 
  mutate(sir_p_value = map2_dbl(obs_n, exp_n, ~ 
                                  poisson.test(x = .x, `T` = .y, 
                                               alternative = "t", 
                                               conf.level = 0.95)$p.value))

inc_rate_output <- calculate_ratio(inc_output_sir, 
                                   c("year", "emergency_flag", 
                                     "tnm_stage_2"), 
                                   obs_inc, exp_inc, 0.95, srr)

# Tidy environment

rm(match_df, reg_mod, reg_mod_output, inc_rate_data, obs_inc_output, 
   exp_inc_output, inc_output, inc_output_sir)



### 6 Regression Models and Incidence Rates by Site ----

# Get a list of all cancer types in analysis

cancers <- joined_data |> 
  distinct(incidence_type) |> 
  pull()

# Get a list of all stages in analysis

stages <- joined_data |> 
  distinct(tnm_stage_2) |> 
  pull()

# Create a data frame for matching to joined_data
# This is every combination of year, sex, age group, flag and incidence_type
# This ensures all populations are counted for the incidence rates, even where
# there are no admissions for a specific category, e.g. if there are no elective
# admissions for men in age group 0 in 2019, joining on this data frame will 
# create a row for this group
# Remove rows for men with breast/cervical cancer and women with prostate cancer

match_df_site <- expand.grid(year = c(start:end), sex = c("1", "2"), 
                             age_group = c(0:18), 
                             emergency_flag = c("Non-Emergency", "Emergency"), 
                             tnm_stage_2 = c(stages), 
                             incidence_type = c(cancers)) |> 
  filter(!(incidence_type == "Breast" & sex == "1") & 
           !(incidence_type == "Prostate" & sex == "2") & 
           !(incidence_type == "Cervical" & sex == "1"))

# Define age groups for admission data
# Count data to get total number of admissions by year, sex, age_group,  
# emergency_flag and incidence_type
# Join on match_df_site to add any groups with no admissions and set number of
# admissions for these groups to 0, as they appear as NAs when joined on
# Join on Scotland populations and also ESP and arrange data

inc_rate_site_data <- joined_data |> 
  mutate(age_group = standard_pop_age_groups(age_in_years)) |> 
  count(year = incidence_year, sex, age_group, emergency_flag, tnm_stage_2, 
        incidence_type) |> 
  full_join(match_df_site) |> 
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .))) |> 
  left_join(scot_pop) |> 
  left_join(esp2013) |> 
  arrange(year, sex, age_group, incidence_type, emergency_flag, tnm_stage_2)

# Define fitted model for breast, cervical and prostate cancers
# Group data by emergency_flag, incidence_type and tnm_stage_2

site_reg_mod_one_sex <- fitted_model_without_sex(
  inc_rate_site_data |> 
    filter(incidence_type %in% c("Breast", "Cervical", "Prostate")), 
  emergency_flag, incidence_type, tnm_stage_2)

# Define fitted model for colorectal, head and neck, lung and upper GI cancers
# Group data by emergency_flag, incidence_type and tnm_stage_2

site_reg_mod_two_sex <- fitted_model_with_sex(
  inc_rate_site_data |> 
    filter(incidence_type %in% c("Colorectal", "Head and Neck", "Lung", 
                                 "Upper GI")), 
  emergency_flag, incidence_type, tnm_stage_2)

# Bind site_reg_mod_one_sex and site_reg_mod_two_sex and drop data column

site_reg_mod <- bind_rows(site_reg_mod_one_sex, site_reg_mod_two_sex) |> 
  select(-data)

# Group by emergency_flag, incidence_type and tnm_stage_2 and nest data
# Do an inner join with site_reg_mod
# Create column for expected values by predicting the model in the mod column
# using the data from the data column
# Unnest data to get data back into standard format

site_reg_mod_output <- inc_rate_site_data |> 
  group_by(emergency_flag, incidence_type, tnm_stage_2) |>
  nest() |> 
  inner_join(site_reg_mod) |> 
  mutate(exp = purrr::map2(mod, data, predict, type = "response")) |> 
  unnest(c(data, exp))

# Calculate incidence rates by year, emergency_flag, incidence_type and sex for
# observed and expected values
# Rename columns to indicate observed or expected

obs_inc_site_output <- calculate_incidence(site_reg_mod_output, 
                                           c("year", "emergency_flag", 
                                             "incidence_type", "tnm_stage_2"),
                                           esp2013, n, pop, 0.95) |> 
  rename_with(~ paste0("obs_", .), c("n", "inc", "l_ci", "u_ci"))

exp_inc_site_output <- calculate_incidence(site_reg_mod_output, 
                                           c("year", "emergency_flag", 
                                             "incidence_type", "tnm_stage_2"),
                                           esp2013, exp, pop, 0.95) |> 
  rename_with(~ paste0("exp_", .), c("n", "inc", "l_ci", "u_ci"))

# Join observed and expected together

inc_site_output <- full_join(obs_inc_site_output, 
                             exp_inc_site_output)

# Calculate standardised incidence ratio and standardised rate ratio

inc_site_output_sir <- calculate_ratio(inc_site_output, 
                                       c("year", "emergency_flag", 
                                         "incidence_type", "tnm_stage_2"), 
                                       obs_n, exp_n, 0.95, sir)

# Calculate p-values for SIRs
# This only works for intergers, so use obs_n and exp_n

inc_site_output_sir <- inc_site_output_sir |> 
  mutate(sir_p_value = map2_dbl(obs_n, exp_n, ~ 
                                  poisson.test(x = .x, `T` = .y, 
                                               alternative = "t", 
                                               conf.level = 0.95)$p.value))

inc_rate_site_output <- calculate_ratio(inc_site_output_sir, 
                                        c("year", "emergency_flag", 
                                          "incidence_type", "tnm_stage_2"), 
                                        obs_inc, exp_inc, 0.95, srr)

# Tidy environment

rm(site_reg_mod_one_sex, site_reg_mod_two_sex, site_reg_mod, 
   site_reg_mod_output, obs_inc_site_output, exp_inc_site_output, 
   inc_site_output, inc_site_output_sir)



### 7 Charts ----

### 7.1 Total ----

# Create chart for emergency presentations by tnm_stage_2
# Define colours based on phsstyles package
# Add axis labels and set y axis to start from 0 and go to specified length

inc_rate_em_plot <- create_chart(inc_rate_output, "Emergency", "tnm_stage_2") +
  scale_colour_manual("TNM Stage",
                      values = phs_colours(c("phs-purple", "phs-magenta",
                                             "phs-blue", "phs-green", 
                                             "phs-teal"))) +
  ylab("Age-Sex Standardised Rate per 100,000") + 
  xlab("Year") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 53))

# # Save chart
# 
# ggsave(here::here("Charts/incidence_rate_emergency_line_plot.png"),
#        plot = incidence_rate_em_plot,
#        height = 14,
#        width = 23,
#        dpi = 500)


# Same as above but with linear regression

inc_rate_em_plot_reg <- create_regression_chart(inc_rate_output, 
                                                "Emergency", "tnm_stage_2") + 
  scale_colour_manual("TNM Stage",
                      values = phs_colours(c("phs-purple", "phs-magenta",
                                             "phs-blue", "phs-green", 
                                             "phs-teal"))) +
  ylab("Age-Sex Standardised Rate per 100,000") + 
  xlab("Year") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 63))

# Save chart

ggsave(here(glue("Charts/{end}/incidence_rate_emergency_stage_line_regression_plot.png")),
       plot = inc_rate_em_plot_reg,
       height = 14,
       width = 23,
       dpi = 500)


# Create chart for non-emergency presentations by tnm_stage_2
# Define colours based on phsstyles package
# Add axis labels and set y axis to start from 0 and go to specified length

inc_rate_el_plot <- create_chart(inc_rate_output, "Non-Emergency", 
                                 "tnm_stage_2") +
  scale_colour_manual("TNM Stage",
                      values = phs_colours(c("phs-purple", "phs-magenta",
                                             "phs-blue", "phs-green", 
                                             "phs-teal"))) +
  ylab("Age-Sex Standardised Rate per 100,000") + 
  xlab("Year") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 95))

# Save chart

ggsave(here(glue("Charts/{end}/incidence_rate_non_emergency_stage_line_plot.png")),
       plot = inc_rate_el_plot,
       height = 14,
       width = 23,
       dpi = 500)


# Same as above but with linear regression

inc_rate_el_plot_reg <- create_regression_chart(inc_rate_output, 
                                                "Non-Emergency", 
                                                "tnm_stage_2") + 
  scale_colour_manual("TNM Stage",
                      values = phs_colours(c("phs-purple", "phs-magenta",
                                             "phs-blue", "phs-green", 
                                             "phs-teal"))) +
  ylab("Age-Sex Standardised Rate per 100,000") + 
  xlab("Year") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 95))

# Save chart

ggsave(here(glue("Charts/{end}/incidence_rate_non_emergency_stage_line_regression_plot.png")),
       plot = inc_rate_el_plot_reg,
       height = 14,
       width = 23,
       dpi = 500)


### 7.2 Site Level ----

# Create regression chart for non-emergency breast presentations by tnm_stage_2
# Define colours based on phsstyles package
# Add axis labels and set y axis to start from 0 and go to specified length

# inc_rate_site_el_breast_plot <- create_regression_chart(
#   inc_rate_site_output |> filter(incidence_type == "Breast"), 
#   "Non-Emergency", "tnm_stage_2") + 
#   scale_colour_manual("TNM Stage",
#                       values = phs_colours(c("phs-purple", "phs-magenta",
#                                              "phs-blue", "phs-green", 
#                                              "phs-teal"))) +
#   ylab("Age-Sex Standardised Rate per 100,000") + 
#   xlab("Year") + 
#   scale_y_continuous(expand = c(0, 0), limits = c(0, 63))
# 
# # Save chart
# 
# ggsave(here::here("Charts/incidence_rate_non_emergency_breast_stage_line_plot.png"),
#        plot = inc_rate_site_el_breast_plot,
#        height = 14,
#        width = 23,
#        dpi = 500)



### 8 Output ----

### 8.1 Incidence Rates ----

# Create a workbook

wb <- createWorkbook()

# Define a header style for workbook

hs <- createStyle(fontColour = "#ffffff", fgFill = "#0078D4",
                  halign = "center", valign = "center", 
                  textDecoration = "bold", border = "TopBottomLeftRight")

# Add sheet for Stage At Diagnosis and write stage_at_diagnosis data
# Rename relevant columns and set column widths to auto for this sheet

addWorksheet(wb, "Stage At Diagnosis", gridLines = FALSE)

writeData(wb, sheet = "Stage At Diagnosis", 
          stage_at_diagnosis |> 
            rename(Year = incidence_year, 
                   Type = emergency_flag, 
                   "TNM Stage" = tnm_stage_2, 
                   "Number of Cancers" = n, 
                   "Percentage" = pc), 
          borders = "all", headerStyle = hs)

setColWidths(wb, sheet = "Stage At Diagnosis", cols = 1:5, widths = "auto")

# Add sheet for Stage At Diagnosis by Site and write stage_at_diagnosis_site
# Rename relevant columns and set column widths to auto for this sheet

addWorksheet(wb, "Stage At Diagnosis by Site", gridLines = FALSE)

writeData(wb, sheet = "Stage At Diagnosis by Site", 
          stage_at_diagnosis_site |> 
            rename(Year = incidence_year, 
                   Type = emergency_flag, 
                   Cancer = incidence_type, 
                   "TNM Stage" = tnm_stage_2, 
                   "Number of Cancers" = n, 
                   "Percentage" = pc), 
          borders = "all", headerStyle = hs)

setColWidths(wb, sheet = "Stage At Diagnosis by Site", cols = 1:6, 
             widths = "auto")

# Add sheet for Incidence Rates and write inc_rate_output data
# Round relevant columns to two decimal places
# Select and rename relevant columns
# Set column widths to auto for this sheet

addWorksheet(wb, "Incidence Rates", gridLines = FALSE)

writeData(wb, sheet = "Incidence Rates", 
          inc_rate_output |> 
            mutate(across(c(obs_inc, obs_l_ci, obs_u_ci, 
                            exp_inc, exp_l_ci, exp_u_ci, 
                            sir, sir_l_ci, sir_u_ci, sir_p_value, 
                            srr, srr_l_ci, srr_u_ci, 
                            exp_n),
                          ~ round_half_up(., 2))) |>
            select(Year = year, 
                   Type = emergency_flag, 
                   "TNM Stage" = tnm_stage_2, 
                   "Observed Number" = obs_n, 
                   "Observed Rate" = obs_inc, 
                   "Observed L 95 CI" = obs_l_ci, 
                   "Observed U 95 CI" = obs_u_ci, 
                   "Expected Number" = exp_n, 
                   "Expected Rate" = exp_inc, 
                   "Expected L 95 CI" = exp_l_ci, 
                   "Expected U 95 CI" = exp_u_ci, 
                   "SIR" = sir, 
                   "SIR L 95 CI" = sir_l_ci, 
                   "SIR U 95 CI" = sir_u_ci, 
                   "SIR P Value" = sir_p_value, 
                   "SRR" = srr, 
                   "SRR L 95 CI" = srr_l_ci, 
                   "SRR U 95 CI" = srr_u_ci), 
          borders = "all", headerStyle = hs, startCol = 1, startRow = 1)

setColWidths(wb, sheet = "Incidence Rates", cols = 1:18, widths = "auto")

# Add sheet for Site Incidence Rates and write inc_rate_site_output data
# Round relevant columns to two decimal places
# Select and rename relevant columns
# Set column widths to auto for this sheet

addWorksheet(wb, "Site Incidence Rates", gridLines = FALSE)

writeData(wb, sheet = "Site Incidence Rates", 
          inc_rate_site_output |> 
            mutate(across(c(obs_inc, obs_l_ci, obs_u_ci, 
                            exp_inc, exp_l_ci, exp_u_ci, 
                            sir, sir_l_ci, sir_u_ci, sir_p_value, 
                            srr, srr_l_ci, srr_u_ci, 
                            exp_n),
                          ~ round_half_up(., 2))) |>
            select(Year = year, 
                   Type = emergency_flag, 
                   "TNM Stage" = tnm_stage_2, 
                   "Cancer" = incidence_type, 
                   "Observed Number" = obs_n, 
                   "Observed Rate" = obs_inc, 
                   "Observed L 95 CI" = obs_l_ci, 
                   "Observed U 95 CI" = obs_u_ci, 
                   "Expected Number" = exp_n, 
                   "Expected Rate" = exp_inc, 
                   "Expected L 95 CI" = exp_l_ci, 
                   "Expected U 95 CI" = exp_u_ci, 
                   "SIR" = sir, 
                   "SIR L 95 CI" = sir_l_ci, 
                   "SIR U 95 CI" = sir_u_ci, 
                   "SIR P Value" = sir_p_value, 
                   "SRR" = srr, 
                   "SRR L 95 CI" = srr_l_ci, 
                   "SRR U 95 CI" = srr_u_ci), 
          borders = "all", headerStyle = hs, startCol = 1, startRow = 1)

setColWidths(wb, sheet = "Site Incidence Rates", cols = 1:19, widths = "auto")

# Save workbook

saveWorkbook(wb, here(glue("Data/{end}/stage_at_diagnosis_{start}_{end}.xlsx")), 
             overwrite = TRUE)



### 8.2 Tables for Report ----

# Create a workbook

report_wb <- createWorkbook()

# Define a header style for workbook

hs <- createStyle(fontColour = "#ffffff", fgFill = "#0078D4",
                  halign = "center", valign = "center", 
                  textDecoration = "bold", border = "TopBottomLeftRight")

# Add sheet for Site Incidence Rates and write inc_rate_site_output data
# Round relevant columns to two decimal places
# Select and rename relevant columns
# Set column widths to auto for this sheet

addWorksheet(report_wb, "SIR Data", gridLines = FALSE)

writeData(report_wb, sheet = "SIR Data", 
          inc_rate_output |> 
            filter(year %in% c(2020, 2021)) |> 
            mutate(across(c(obs_inc, obs_l_ci, obs_u_ci, 
                            exp_inc, exp_l_ci, exp_u_ci, 
                            sir, sir_l_ci, sir_u_ci, sir_p_value, 
                            srr, srr_l_ci, srr_u_ci, 
                            exp_n),
                          ~ round_half_up(., 2))) |>
            mutate(sir_ci = paste0("(", sir_l_ci, ", ", sir_u_ci, ")")) |> 
            select(Year = year, 
                   Type = emergency_flag, 
                   "Observed Number" = obs_n, 
                   "Expected Number" = exp_n, 
                   "SIR" = sir, 
                   "SIR CI" = sir_ci, 
                   "SIR P Value" = sir_p_value), 
          borders = "all", headerStyle = hs)

setColWidths(report_wb, sheet = "SIR Data", cols = 1:7, widths = "auto")

# Add sheet for Sex Incidence Rates and write inc_rate_site_sex_output data
# Round relevant columns to two decimal places
# Select and rename relevant columns
# Set column widths to auto for this sheet

addWorksheet(report_wb, "SRR Data", gridLines = FALSE)

writeData(report_wb, sheet = "SRR Data", 
          inc_rate_site_output |> 
            filter(year %in% c(2020, 2021)) |> 
            mutate(across(c(obs_inc, obs_l_ci, obs_u_ci, 
                            exp_inc, exp_l_ci, exp_u_ci, 
                            sir, sir_l_ci, sir_u_ci, 
                            srr, srr_l_ci, srr_u_ci, 
                            exp_n),
                          ~ round_half_up(., 2))) |>
            mutate(srr_ci = paste0("(", srr_l_ci, ", ", srr_u_ci, ")")) |> 
            select(Year = year, 
                   Type = emergency_flag, 
                   "Observed Rate" = obs_inc, 
                   "Expected Rate" = exp_inc, 
                   "SRR" = srr, 
                   "SRR CI" = srr_ci), 
          borders = "all", headerStyle = hs, startCol = 1, startRow = 1)

setColWidths(report_wb, sheet = "SRR Data", cols = 1:6, widths = "auto")


# Save workbook

saveWorkbook(report_wb, here(glue("Data/{end}/report_tables_stage.xlsx")), 
             overwrite = TRUE)


### 8.3 Tables for ENCR ----

# Create a workbook

report_wb <- createWorkbook()

# Define a header style for workbook

hs <- createStyle(fontColour = "#ffffff", fgFill = "#0078D4",
                  halign = "center", valign = "center", 
                  textDecoration = "bold", border = "TopBottomLeftRight")

# Add sheet for Site Incidence Rates and write inc_rate_site_output data
# Round relevant columns to two decimal places
# Select and rename relevant columns
# Set column widths to auto for this sheet

addWorksheet(report_wb, "SIR Data", gridLines = FALSE)

writeData(report_wb, sheet = "SIR Data", 
          inc_rate_output |> 
            filter(year %in% c(2020, 2021)) |> 
            mutate(across(c(obs_inc, obs_l_ci, obs_u_ci, 
                            exp_inc, exp_l_ci, exp_u_ci, 
                            sir, sir_l_ci, sir_u_ci, sir_p_value, 
                            srr, srr_l_ci, srr_u_ci, 
                            exp_n),
                          ~ round_half_up(., 2))) |>
            select(Year = year, 
                   Type = emergency_flag, 
                   Stage = tnm_stage_2, 
                   "SIR" = sir) |> 
            pivot_wider(names_from = "Year", values_from = "SIR", 
                        names_glue = "SIR_{Year}"), 
          borders = "all", headerStyle = hs)

setColWidths(report_wb, sheet = "SIR Data", cols = 1:7, widths = "auto")

# Save workbook

saveWorkbook(report_wb, here(glue("Data/{end}/report_tables_stage_encr.xlsx")), 
             overwrite = TRUE)

