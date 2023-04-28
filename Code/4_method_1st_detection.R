#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4_method_1st_detection
# Calum Purdie
# 11/05/2022
# Data extraction/preparation
# Written/run on R Studio Server
# R version 3.6.1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



### 1 Housekeeping ----

source(here::here("Code/1_housekeeping.R"))



### 2 Method 1st Detection ----

# Recode method_1st_detection to add descriptions

joined_data <- joined_data %>% 
  mutate(method_1st_detection = case_when(
    method_1st_detection == 1 ~ "1 (Screening Examination)", 
    method_1st_detection == 2 ~ "2 (Incidental Finding)", 
    method_1st_detection == 3 ~ "3 (Clinical Presentation)", 
    method_1st_detection == 4 ~ "4 (Incidental Finding at Autopsy)", 
    method_1st_detection == 5 ~ "5 (Interval Cancer)", 
    method_1st_detection == 8 ~ "8 (Other)", 
    method_1st_detection == 9 ~ "9 (Not Known)", 
    TRUE ~ "Error"
  ))

# Count data by incidence_year, emergency_flag and method_1st_detection
# Complete data to get all method_1st_detection values for incidence_year and 
# emergency_flag, filling in any blanks with 0
# Pivot data wider, taking new names from incidence_year and values from n
# Select columns

method_1st_detection_n <- joined_data %>% 
  count(incidence_year, emergency_flag, method_1st_detection) %>% 
  complete(method_1st_detection, nesting(incidence_year, emergency_flag), 
           fill = list(n = 0)) %>% 
  pivot_wider(names_from = incidence_year, values_from = n) %>% 
  select("Presentation" = emergency_flag, 
         "Method 1st Detection" = method_1st_detection, 
         everything())

# Count data by incidence_year, emergency_flag and method_1st_detection
# Complete data to get all method_1st_detection values for incidence_year and 
# emergency_flag, filling in any blanks with 0
# Group by incidence_year and emergency_flag
# Calculate percentages for each group and drop the n column
# Pivot data wider, taking new names from incidence_year and values from pc
# Select columns

method_1st_detection_pc <- joined_data %>% 
  count(incidence_year, emergency_flag, method_1st_detection) %>% 
  complete(method_1st_detection, nesting(incidence_year, emergency_flag), 
           fill = list(n = 0)) %>% 
  group_by(incidence_year, emergency_flag) %>% 
  mutate(pc = round_half_up(n * 100 / sum(n), 1)) %>% 
  ungroup() %>% 
  select(-n) %>% 
  pivot_wider(names_from = incidence_year, values_from = pc) %>% 
  select("Presentation" = emergency_flag, 
         "Method 1st Detection" = method_1st_detection, 
         everything())



### 3 Incidence Rates ----

# Get a list of all cancer types in analysis

method <- joined_data %>% 
  distinct(method_1st_detection) %>% 
  pull()

# Create a data frame for matching to joined_data
# This is every combination of year, sex, age group, flag and incidence_type
# This ensures all populations are counted for the incidence rates, even where
# there are no admissions for a specific category, e.g. if there are no elective
# admissions for men in age group 0 in 2019, joining on this data frame will 
# create a row for this group

match_df_site <- expand.grid(year = c(start:end), sex = c("1", "2"), 
                             age_group = c(0:18), 
                             emergency_flag = c("Non-Emergency", "Emergency"), 
                             method_1st_detection = c(method))

# Define age groups for admission data
# Count data to get total number of admissions by year, sex, age_group,  
# emergency_flag and incidence_type
# Join on match_df_site to add any groups with no admissions and set number of
# admissions for these groups to 0, as they appear as NAs when joined on
# Join on Scotland populations and also ESP and arrange data

inc_rate_m1d_data <- joined_data %>% 
  mutate(age_group = standard_pop_age_groups(age_in_years)) %>% 
  count(year = incidence_year, sex, age_group, emergency_flag, 
        method_1st_detection) %>% 
  full_join(match_df_site) %>% 
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .))) %>% 
  left_join(scot_pop) %>% 
  left_join(esp2013) %>% 
  arrange(year, sex, age_group, emergency_flag, method_1st_detection)

# Define fitted model for all cancers
# Group data by emergency_flag and method_1st_detection
# Drop data column

m1d_reg_mod <- fitted_model_with_sex(inc_rate_m1d_data, emergency_flag, 
                                     method_1st_detection) %>% 
  select(-data)

# Group by emergency_flag and method_1st_detection and nest data
# Do an inner join with m1d_reg_mod
# Create column for expected values by predicting the model in the mod column
# using the data from the data column
# Unnest data to get data back into standard format

m1d_reg_mod_output <- inc_rate_m1d_data %>% 
  group_by(emergency_flag, method_1st_detection) %>%
  nest() %>% 
  inner_join(m1d_reg_mod) %>% 
  mutate(exp = purrr::map2(mod, data, predict, type = "response")) %>% 
  unnest(c(data, exp))

# Calculate incidence rates by year, emergency_flag and method_1st_detection for
# observed and expected values
# Rename columns to indicate observed or expected

obs_inc_m1d_output <- calculate_incidence(m1d_reg_mod_output, 
                                           c("year", "emergency_flag", 
                                             "method_1st_detection"),
                                           esp2013, n, pop, 0.95) %>% 
  rename_with(~ paste0("obs_", .), c("n", "inc", "l_ci", "u_ci"))

exp_inc_m1d_output <- calculate_incidence(m1d_reg_mod_output, 
                                           c("year", "emergency_flag", 
                                             "method_1st_detection"),
                                           esp2013, exp, pop, 0.95) %>% 
  rename_with(~ paste0("exp_", .), c("n", "inc", "l_ci", "u_ci"))

# Join observed and expected together

inc_m1d_output <- full_join(obs_inc_m1d_output, exp_inc_m1d_output)

# Calculate standardised incidence ratio and standardised rate ratio

inc_m1d_output_sir <- calculate_ratio(inc_m1d_output, 
                                       c("year", "emergency_flag", 
                                         "method_1st_detection"), 
                                       obs_n, exp_n, 0.95, sir)

inc_rate_m1d_output <- calculate_ratio(inc_m1d_output_sir, 
                                        c("year", "emergency_flag", 
                                          "method_1st_detection"), 
                                        obs_inc, exp_inc, 0.95, srr)



### 4 Charts ----

# Create chart for emergency presentations by method_1st_detection
# Define colours based on phsstyles package
# Add axis labels and set y axis to start from 0 and go to specified length

inc_rate_em_plot <- create_chart(inc_rate_m1d_output, 
                                 "Emergency", "method_1st_detection") +
  scale_colour_manual("Method 1st Detection",
                      values = phs_colours(c("phs-graphite", "phs-magenta",
                                             "phs-blue", "phs-green", 
                                             "phs-rust", "phs-teal", 
                                             "phs-purple", "phs-liberty"))) +
  ylab("Age-Sex Standardised Rate per 100,000") + 
  xlab("Year")

# # Save chart
# 
# ggsave(here::here("Charts/incidence_rate_emergency_line_plot.png"),
#        plot = incidence_rate_em_plot,
#        height = 14,
#        width = 23,
#        dpi = 500)


# Same as above but with linear regression

inc_rate_em_plot_reg <- create_regression_chart(inc_rate_m1d_output, 
                                                   "Emergency", 
                                                   "method_1st_detection") +
  scale_colour_manual("Method 1st Detection",
                      values = phs_colours(c("phs-graphite", "phs-magenta",
                                             "phs-blue", "phs-green", 
                                             "phs-rust", "phs-teal", 
                                             "phs-purple", "phs-liberty"))) +
  ylab("Age-Sex Standardised Rate per 100,000") + 
  xlab("Year") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 113))

# Save chart

ggsave(here(glue("Charts/{end}/incidence_rate_emergency_method_line_regression_plot.png")),
       plot = inc_rate_em_plot_reg,
       height = 14,
       width = 23,
       dpi = 500)

# Create chart for non-emergency presentations by method_1st_detection
# Define colours based on phsstyles package
# Add axis labels and set y axis to start from 0 and go to specified length

inc_rate_el_plot <- create_chart(inc_rate_m1d_output, 
                                 "Non-Emergency", "method_1st_detection") +
  scale_colour_manual("Method 1st Detection",
                      values = phs_colours(c("phs-graphite", "phs-magenta",
                                             "phs-blue", "phs-green", 
                                             "phs-rust", "phs-teal", 
                                             "phs-purple", "phs-liberty"))) +
  ylab("Age-Sex Standardised Rate per 100,000") + 
  xlab("Year") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 293))


# Same as above but with linear regression

inc_rate_el_plot_reg <- create_regression_chart(inc_rate_m1d_output, 
                                                       "Non-Emergency", 
                                                       "method_1st_detection") +
  scale_colour_manual("Method 1st Detection",
                      values = phs_colours(c("phs-graphite", "phs-magenta",
                                             "phs-blue", "phs-green", 
                                             "phs-rust", "phs-teal", 
                                             "phs-purple", "phs-liberty"))) +
  ylab("Age-Sex Standardised Rate per 100,000") + 
  xlab("Year") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 293))

# Save chart

ggsave(here(glue("Charts/{end}/incidence_rate_non_emergency_method_line_regression_plot.png")),
       plot = inc_rate_el_plot_reg,
       height = 14,
       width = 23,
       dpi = 500)



### 5 Output ----

# Create a workbook

wb_m1d <- createWorkbook()

# Define a header style for workbook

hs <- createStyle(fontColour = "#ffffff", fgFill = "#0078D4",
                  halign = "center", valign = "center", 
                  textDecoration = "bold", border = "TopBottomLeftRight")

# Add sheet for Method 1st detectiom and write data
# Select and rename relevant columns
# Set column widths to auto for this sheet

addWorksheet(wb_m1d, "Method 1st Detection", gridLines = FALSE)

writeData(wb_m1d, "Method 1st Detection", 
          "Number of Cancers  by Presentation Type for 2015-2020 by Method 1st Detection", 
          startRow = 1, startCol = 1)

writeData(wb_m1d, sheet = "Method 1st Detection", 
          method_1st_detection_n, 
          borders = "all", headerStyle = hs, 
          startRow = 3, startCol = 1)

writeData(wb_m1d, "Method 1st Detection", 
          "Percentage of Cancers by Presentation Type for 2015-2020 by Method 1st Detection", 
          startRow = 19, startCol = 1)

writeData(wb_m1d, sheet = "Method 1st Detection", 
          method_1st_detection_pc, 
          borders = "all", headerStyle = hs, 
          startRow = 21, startCol = 1)

setColWidths(wb_m1d, sheet = "Method 1st Detection", cols = 2:8, widths = "auto")

# Add sheet for incidence rates and write inc_rate_m1d_output data
# Round relevant columns to two decimal places
# Select and rename relevant columns
# Set column widths to auto for this sheet

addWorksheet(wb_m1d, "Incidence Rates", gridLines = FALSE)

writeData(wb_m1d, sheet = "Incidence Rates", 
          inc_rate_m1d_output %>% 
            mutate(across(c(obs_inc, obs_l_ci, obs_u_ci, 
                            exp_inc, exp_l_ci, exp_u_ci, 
                            sir, sir_l_ci, sir_u_ci, 
                            srr, srr_l_ci, srr_u_ci, 
                            exp_n),
                          ~ round_half_up(., 2))) %>%
          select(Year = year, 
                 Type = emergency_flag, 
                 "Method 1st Detection" = method_1st_detection, 
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
                 "SRR" = srr, 
                 "SRR L 95 CI" = srr_l_ci, 
                 "SRR U 95 CI" = srr_u_ci), 
          borders = "all", headerStyle = hs)

setColWidths(wb_m1d, sheet = "Incidence Rates", cols = 1:17, widths = "auto")

# Save workbook

saveWorkbook(wb_m1d, here(glue("Data/{end}/method_1st_detection_{start}_{end}.xlsx")), 
             overwrite = TRUE)
