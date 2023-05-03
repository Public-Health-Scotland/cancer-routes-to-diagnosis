#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3_incidence_rates
# Calum Purdie
# 03/11/2021
# Data extraction/preparation
# Written/run on R Studio Server
# R version 3.6.1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



### 1 Housekeeping ----

source(here::here("Code/1_housekeeping.R"))



### 2 Regression Models and Incidence Rates ----

### 2.1 Cancer Site Split ----

# Get a list of all cancer types in analysis

cancer_types <- joined_data %>% 
  distinct(incidence_type) %>% 
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
                             incidence_type = c(cancer_types)) %>% 
  filter(!(incidence_type == "Breast" & sex == "1") & 
           !(incidence_type == "Prostate" & sex == "2") & 
           !(incidence_type == "Cervical" & sex == "1"))

# Define age groups for admission data
# Count data to get total number of admissions by year, sex, age_group,  
# emergency_flag and incidence_type
# Join on match_df_site to add any groups with no admissions and set number of
# admissions for these groups to 0, as they appear as NAs when joined on
# Join on Scotland populations and also ESP and then arrange data

inc_rate_site_data <- joined_data %>% 
  mutate(age_group = standard_pop_age_groups(age_in_years)) %>% 
  count(year = incidence_year, sex, age_group, emergency_flag, incidence_type) %>% 
  full_join(match_df_site) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .))) %>% 
  left_join(scot_pop) %>% 
  left_join(esp2013) %>% 
  arrange(year, sex, age_group, incidence_type, emergency_flag)

# inc_rate_site_data <- inc_rate_site_data %>% 
#   mutate(age_group = case_when(age_group >= 0 & age_group <= 6 ~ 6, 
#                                TRUE ~ age_group)) %>% 
#   group_by(year, sex, age_group, emergency_flag, incidence_type) %>% 
#   summarise(n = sum(n), 
#             pop = sum(pop), 
#             esp2013 = sum(esp2013))

# Define fitted model for breast, cervical and prostate cancers
# Group data by emergency_flag and incidence_type

site_reg_mod_without_sex <- fitted_model_without_sex(
  inc_rate_site_data %>% 
    filter(incidence_type %in% c("Breast", "Cervical", "Prostate")), 
  emergency_flag, incidence_type)

# Define fitted model for colorectal, head and neck, lung and upper GI cancers
# Group data by emergency_flag and incidence_type

site_reg_mod_with_sex <- fitted_model_with_sex(
  inc_rate_site_data %>% 
    filter(incidence_type %in% c("Colorectal", "Head and Neck", "Lung", 
                                 "Upper GI")), 
  emergency_flag, incidence_type)

# Bind site_reg_mod_one_sex and site_reg_mod_two_sex and drop data column

site_reg_mod <- bind_rows(site_reg_mod_without_sex, site_reg_mod_with_sex) %>% 
  select(-data)

# Group by emergency_flag and incidence_type and nest data
# Do an inner join with site_reg_mod
# Create column for expected values by predicting the model in the mod column
# using the data from the data column
# Unnest data to get data back into standard format

site_reg_mod_output <- inc_rate_site_data %>% 
  group_by(emergency_flag, incidence_type) %>%
  nest() %>% 
  inner_join(site_reg_mod) %>% 
  mutate(exp = purrr::map2(mod, data, predict, type = "response")) %>% 
  unnest(c(data, exp))

# Calculate goodness of fit for Poisson regression models
# Calculate the squared sum of residuals using the Pearson method
# Use sum_res in pchisq with the residuals from the model as degrees of freedom

goodness_of_fit <- site_reg_mod_output %>% 
  mutate(sum_res = sum(residuals(mod[[1]], type="pearson")^2), 
         p_chisq = pchisq(sum_res, mod[[1]]$df.residual, lower.tail = FALSE))

goodness_of_fit %>% count(emergency_flag, incidence_type, sum_res, p_chisq)

# Calculate incidence rates by year, emergency_flag and incidence_type for
# observed and expected values
# Rename columns to indicate observed or expected

obs_inc_site_output <- calculate_incidence(site_reg_mod_output, 
                                           c("year", "emergency_flag", 
                                             "incidence_type"),
                                           esp2013, n, pop, 0.95) %>% 
  rename_with(~ paste0("obs_", .), c("n", "inc", "l_ci", "u_ci"))

exp_inc_site_output <- calculate_incidence(site_reg_mod_output, 
                                           c("year", "emergency_flag", 
                                             "incidence_type"),
                                           esp2013, exp, pop, 0.95) %>% 
  rename_with(~ paste0("exp_", .), c("n", "inc", "l_ci", "u_ci"))

# Join observed and expected together

inc_site_output <- full_join(obs_inc_site_output, exp_inc_site_output)

# Calculate standardised incidence ratio

inc_site_output_sir <- calculate_ratio(inc_site_output, 
                                       c("year", "emergency_flag", 
                                         "incidence_type"), 
                                       obs_n, exp_n, 0.95, sir)

# Calculate p-values for SIRs
# This only works for intergers, so use obs_n and exp_n

inc_site_output_sir <- inc_site_output_sir %>% 
  mutate(sir_p_value = map2_dbl(obs_n, exp_n, ~ 
                                  poisson.test(x = .x, `T` = .y, 
                                               alternative = "t", 
                                               conf.level = 0.95)$p.value))

# Calculate standardised rate ratio

inc_rate_site_output <- calculate_ratio(inc_site_output_sir, 
                                        c("year", "emergency_flag", 
                                          "incidence_type"), 
                                        obs_inc, exp_inc, 0.95, srr)



# Tidy environment

rm(site_reg_mod_without_sex, site_reg_mod_with_sex, site_reg_mod, 
   site_reg_mod_output, obs_inc_site_output, exp_inc_site_output, 
   inc_site_output, inc_site_output_sir)


### 2.2 Cancer Site and Sex Split ----

# Define fitted model for all cancers
# Group data by emergency_flag, incidence_type and sex
# Drop data column

site_sex_reg_mod <- fitted_model_without_sex(inc_rate_site_data,
                                             emergency_flag, incidence_type, 
                                             sex) %>% 
  select(-data)

# Group by emergency_flag, incidence_type and sex and nest data
# Do an inner join with site_sex_reg_mod
# Create column for expected values by predicting the model in the mod column
# using the data from the data column
# Unnest data to get data back into standard format

site_sex_reg_mod_output <- inc_rate_site_data %>% 
  group_by(emergency_flag, incidence_type, sex) %>%
  nest() %>% 
  inner_join(site_sex_reg_mod) %>% 
  mutate(exp = purrr::map2(mod, data, predict, type = "response")) %>% 
  unnest(c(data, exp))

# Calculate incidence rates by year, emergency_flag, incidence_type and sex for
# observed and expected values
# Rename columns to indicate observed or expected

obs_inc_site_sex_output <- calculate_incidence(site_sex_reg_mod_output, 
                                               c("year", "emergency_flag", 
                                                 "incidence_type", "sex"),
                                               esp2013, n, pop, 0.95) %>% 
  rename_with(~ paste0("obs_", .), c("n", "inc", "l_ci", "u_ci"))

exp_inc_site_sex_output <- calculate_incidence(site_sex_reg_mod_output, 
                                               c("year", "emergency_flag", 
                                                 "incidence_type", "sex"),
                                               esp2013, exp, pop, 0.95) %>% 
  rename_with(~ paste0("exp_", .), c("n", "inc", "l_ci", "u_ci"))

# Join observed and expected together

inc_site_sex_output <- full_join(obs_inc_site_sex_output, 
                                 exp_inc_site_sex_output)

# Calculate standardised incidence ratio

inc_site_sex_output_sir <- calculate_ratio(inc_site_sex_output, 
                                           c("year", "emergency_flag", 
                                             "incidence_type", "sex"), 
                                           obs_n, exp_n, 0.95, sir)

# Calculate p-values for SIRs
# This only works for intergers, so use obs_n and exp_n

inc_site_sex_output_sir <- inc_site_sex_output_sir %>% 
  mutate(sir_p_value = map2_dbl(obs_n, exp_n, ~ 
                                  poisson.test(x = .x, `T` = .y, 
                                               alternative = "t", 
                                               conf.level = 0.95)$p.value))

# Calculate standardised rate ratio

inc_rate_site_sex_output <- calculate_ratio(inc_site_sex_output_sir, 
                                            c("year", "emergency_flag", 
                                              "incidence_type", "sex"), 
                                            obs_inc, exp_inc, 0.95, srr)

# Tidy environment

rm(match_df_site, inc_rate_site_data, site_sex_reg_mod, site_sex_reg_mod_output, 
   obs_inc_site_sex_output, exp_inc_site_sex_output, inc_site_sex_output, 
   inc_site_sex_output_sir)


### 2.3 Cancer Site and SIMD Split ----

# Get a list of all SIMD quintiles in analysis
# Exclude blank quintiles

simd_types <- joined_data %>% 
  distinct(simd2020v2_sc_quintile) %>% 
  filter(!is.na(simd2020v2_sc_quintile)) %>% 
  pull()

# Create a data frame for matching to joined_data
# This is every possible combination of year, sex, age group, emergency_flag, 
# incidence_type and simd2020v2_sc_quintile
# This ensures all populations are counted for the incidence rates, even where
# there are no admissions for a specific category, e.g. if there are no elective
# admissions for men in age group 0 in 2019, joining on this data frame will 
# create a row for this group
# Remove rows for men with breast/cervical cancer and women with prostate cancer

match_df_simd <- expand.grid(year = c(start:end), sex = c("1", "2"), 
                             age_group = c(0:18), 
                             emergency_flag = c("Non-Emergency", "Emergency"), 
                             incidence_type = c(cancer_types), 
                             simd2020v2_sc_quintile = c(simd_types)) %>% 
  filter(!(incidence_type == "Breast" & sex == "1") & 
           !(incidence_type == "Prostate" & sex == "2") & 
           !(incidence_type == "Cervical" & sex == "1"))

# Exclude any rows with blank SIMD quintile
# Define age groups for admission data
# Count data to get total number of admissions by year, sex, age_group,  
# emergency_flag, incidence_type and simd2020v2_sc_quintile
# Join on match_df_simd to add any groups with no admissions and set number of
# admissions for these groups to 0, as they appear as NAs when joined on
# Join on Scotland populations and also ESP and arrange data

inc_rate_simd_data <- joined_data %>% 
  filter(!is.na(simd2020v2_sc_quintile)) %>% 
  mutate(age_group = standard_pop_age_groups(age_in_years)) %>% 
  count(year = incidence_year, sex, age_group, emergency_flag, incidence_type, 
        simd2020v2_sc_quintile) %>% 
  full_join(match_df_simd) %>% 
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .))) %>% 
  left_join(simd_pop) %>% 
  left_join(esp2013) %>% 
  arrange(year, sex, age_group, incidence_type, emergency_flag, 
          simd2020v2_sc_quintile)

# Define fitted model for breast, cervical and prostate cancers
# Group data by emergency_flag, incidence_type and simd2020v2_sc_quintile

simd_reg_mod_without_sex <- fitted_model_without_sex(
  inc_rate_simd_data %>% 
    filter(incidence_type %in% c("Breast", "Cervical", "Prostate")), 
  emergency_flag, incidence_type, simd2020v2_sc_quintile)

# Define fitted model for colorectal, head and neck, lung and upper GI cancers
# Group data by emergency_flag, incidence_type and simd2020v2_sc_quintile

simd_reg_mod_with_sex <- fitted_model_with_sex(
  inc_rate_simd_data %>% 
    filter(incidence_type %in% c("Colorectal", "Head and Neck", "Lung", 
                                 "Upper GI")), 
  emergency_flag, incidence_type, simd2020v2_sc_quintile)

# Bind site_reg_mod_one_sex and site_reg_mod_two_sex and drop data column

simd_reg_mod <- bind_rows(simd_reg_mod_without_sex, simd_reg_mod_with_sex) %>% 
  select(-data)

# Group by emergency_flag and incidence_type and nest data
# Do an inner join with site_reg_mod
# Create column for expected values by predicting the model in the mod column
# using the data from the data column
# Unnest data to get data back into standard format

simd_reg_mod_output <- inc_rate_simd_data %>% 
  group_by(emergency_flag, incidence_type, simd2020v2_sc_quintile) %>%
  nest() %>% 
  inner_join(simd_reg_mod) %>% 
  mutate(exp = purrr::map2(mod, data, predict, type = "response")) %>% 
  unnest(c(data, exp))

# Calculate incidence rates by year, emergency_flag, incidence_type and 
# simd2020v2_sc_quintile for observed and expected values
# Rename columns to indicate observed or expected

obs_inc_site_simd_output <- calculate_incidence(simd_reg_mod_output, 
                                                c("year", "emergency_flag", 
                                                  "incidence_type", 
                                                  "simd2020v2_sc_quintile"),
                                                esp2013, n, pop, 0.95) %>% 
  rename_with(~ paste0("obs_", .), c("n", "inc", "l_ci", "u_ci"))

exp_inc_site_simd_output <- calculate_incidence(simd_reg_mod_output, 
                                                c("year", "emergency_flag", 
                                                  "incidence_type", 
                                                  "simd2020v2_sc_quintile"),
                                                esp2013, exp, pop, 0.95) %>% 
  rename_with(~ paste0("exp_", .), c("n", "inc", "l_ci", "u_ci"))

# Join observed and expected together

inc_site_simd_output <- full_join(obs_inc_site_simd_output, 
                                  exp_inc_site_simd_output)

# Calculate standardised incidence ratio

inc_site_simd_output_sir <- calculate_ratio(inc_site_simd_output, 
                                            c("year", "emergency_flag", 
                                              "incidence_type", 
                                              "simd2020v2_sc_quintile"), 
                                            obs_n, exp_n, 0.95, sir)

# Calculate p-values for SIRs
# This only works for intergers, so use obs_n and exp_n

inc_site_simd_output_sir <- inc_site_simd_output_sir %>% 
  mutate(sir_p_value = map2_dbl(obs_n, exp_n, ~ 
                                  poisson.test(x = .x, `T` = .y, 
                                               alternative = "t", 
                                               conf.level = 0.95)$p.value))

# Calculate standardised rate ratio

inc_rate_site_simd_output <- calculate_ratio(inc_site_simd_output_sir, 
                                             c("year", "emergency_flag", 
                                               "incidence_type", 
                                               "simd2020v2_sc_quintile"), 
                                             obs_inc, exp_inc, 0.95, srr)

# Tidy environment

rm(match_df_simd, simd_reg_mod, simd_reg_mod_without_sex, simd_reg_mod_with_sex, 
   simd_reg_mod_output, inc_rate_simd_data, obs_inc_site_simd_output, 
   exp_inc_site_simd_output, inc_site_simd_output, inc_site_simd_output_sir)



### 3 Charts ----

### 3.1 Cancer Site Split ----

# Create chart for emergency presentations by cancer site
# Define colours based on phsstyles package
# Add axis labels and set y axis to start from 0 and go to specified length

inc_rate_em_plot <- create_chart(inc_rate_site_output, 
                                 "Emergency", "incidence_type") +
  scale_colour_manual("Cancer",
                      values = phs_colours(c("phs-purple", "phs-magenta",
                                             "phs-blue", "phs-green", 
                                             "phs-teal", "phs-liberty", 
                                             "phs-graphite"))) +
  ylab("Age-Sex Standardised Rate per 100,000") + 
  xlab("Year") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 48))

# Save chart

ggsave(here(glue("Charts/{end}/incidence_rate_emergency_line_plot.png")),
       plot = inc_rate_em_plot,
       height = 14,
       width = 23,
       dpi = 500)

# Same as above but with linear regression

inc_rate_em_plot_reg <- create_regression_chart(inc_rate_site_output, 
                                                "Emergency", 
                                                "incidence_type") + 
  scale_colour_manual("Cancer",
                      values = phs_colours(c("phs-purple", "phs-magenta",
                                             "phs-blue", "phs-green", 
                                             "phs-teal", "phs-liberty", 
                                             "phs-graphite"))) +
  ylab("Age-Sex Standardised Rate per 100,000") + 
  xlab("Year") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 48))

# Save chart

ggsave(here(glue("Charts/{end}/incidence_rate_emergency_line_regression_plot.png")),
       plot = inc_rate_em_plot_reg,
       height = 14,
       width = 23,
       dpi = 500)

# Create chart for non-emergency presentations by cancer site
# Define colours based on phsstyles package
# Add axis labels and set y axis to start from 0 and go to specified length

inc_rate_el_plot <- create_chart(inc_rate_site_output, 
                                 "Non-Emergency", "incidence_type") +
  scale_colour_manual("Cancer",
                      values = phs_colours(c("phs-purple", "phs-magenta",
                                             "phs-blue", "phs-green", 
                                             "phs-teal", "phs-liberty", 
                                             "phs-graphite"))) +
  ylab("Age-Sex Standardised Rate per 100,000") + 
  xlab("Year") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 183))

# Save chart

ggsave(here(glue("Charts/{end}/incidence_rate_non_emergency_line_plot.png")),
       plot = inc_rate_el_plot,
       height = 14,
       width = 23,
       dpi = 500)

# Same as above but with linear regression

# inc_rate_el_plot_reg <- create_regression_chart(inc_rate_site_output, 
#                                                        "Non-Emergency", 
#                                                        "incidence_type") + 
#   scale_colour_manual("Cancer",
#                       values = phs_colours(c("phs-purple", "phs-magenta",
#                                              "phs-blue", "phs-green", 
#                                              "phs-teal", "phs-liberty", 
#                                              "phs-graphite"))) +
#   ylab("Age-Sex Standardised Rate per 100,000") + 
#   xlab("Year") + 
#   scale_y_continuous(expand = c(0, 0), limits = c(0, 155))

# inc_rate_el_plot_reg <- create_regression_chart(
#   inc_rate_site_output %>% 
#     mutate(plot_label = if_else(year == max(year), 
#                                 incidence_type, NA_character_)), 
#                         "Non-Emergency", 
#                         "incidence_type") + 
#   scale_colour_manual("Cancer",
#                       values = phs_colours(c("phs-purple", "phs-magenta",
#                                              "phs-blue", "phs-green", 
#                                              "phs-teal", "phs-liberty", 
#                                              "phs-graphite"))) +
#   ylab("Age-Sex Standardised Rate per 100,000") + 
#   xlab("Year") + 
#   scale_y_continuous(expand = c(0, 0), limits = c(0, 185)) + 
#   geom_text_repel(aes(y = obs_inc,
#                       label = gsub("^.*$", " ", plot_label)), # This will force the correct position of the link's right end.
#                   segment.colour = NA,
#                   size = 6,
#                   box.padding = 0.1,
#                   point.padding = 0.6,
#                   nudge_x = 0.7,
#                   nudge_y = 1,
#                   force = 0.5,
#                   hjust = 0,
#                   # direction="y",
#                   na.rm = TRUE,
#                   # xlim = as.Date(c("2021-02-16", "2021-03-01")),
#                   # ylim = c(0,150),
#   ) +
#   geom_text_repel(data = . %>% filter(!is.na(plot_label)),
#                    aes(y = obs_inc, label = str_wrap(plot_label, 10), 
#                        color = factor(plot_label)),
#                    segment.alpha = 0, ## This will 'hide' the link
#                    # segment.curvature = -0.1,
#                    # segment.square = TRUE,
#                    size = 9, 
#                    fontface = "bold", 
#                    # segment.color = 'grey',
#                    box.padding = 0,
#                    point.padding = 0, 
#                    nudge_x = 0.4,
#                    nudge_y = 1,
#                    force = 0.5,
#                    hjust = 0,
#                    # direction = "y",
#                    na.rm = TRUE
#   ) + 
#   scale_fill_manual(values = phs_colours(c("phs-purple", "phs-graphite",
#                                            "phs-green", "phs-liberty", 
#                                            "phs-blue", "phs-magenta", 
#                                            "phs-teal")))

inc_rate_el_plot_reg <- create_regression_chart(inc_rate_site_output, 
                                                "Non-Emergency", "incidence_type") +
  scale_colour_manual("Cancer",
                      values = phs_colours(c("phs-purple", "phs-magenta",
                                             "phs-blue", "phs-green", 
                                             "phs-teal", "phs-liberty", 
                                             "phs-graphite"))) +
  ylab("Age-Sex Standardised Rate per 100,000") + 
  xlab("Year") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 183))

inc_rate_el_plot_reg

# Save chart

ggsave(here(glue("Charts/{end}/incidence_rate_non_emergency_line_regression_plot.png")),
       plot = inc_rate_el_plot_reg,
       height = 14,
       width = 23,
       dpi = 500)


### 3.2 Cancer Site and Sex Split ----

# Create chart for non-emergency presentations by sex
# Only include cancers with data for more than one sex in analysis
# Define year as character field and recode sex
# Define ggplot aesthetics and create a line chart with one line for each cancer
# Add 95% confidence interval error bars
# Define colours based on phsstyles package
# Add facet_wrap to get one chart for each cancer type
# Add axis labels and set y axis to start from 0 and go to specified length
# Add theme details to adjust text - this will make the chart look odd in R but
# looks fine once saved out

inc_rate_sex_el_plot <- inc_rate_site_sex_output %>% 
  filter(emergency_flag == "Non-Emergency" & 
           incidence_type %in% c("Colorectal", "Head and Neck", "Lung", 
                                 "Upper GI")) %>% 
  mutate(year = as.character(year), 
         sex = recode(sex, "1" = "Male", "2" = "Female")) %>% 
  ggplot(aes(x = year, y = obs_inc, group = sex)) +
  geom_line(aes(colour = sex), size = 2) + 
  geom_errorbar(aes(ymin = obs_l_ci, ymax = obs_u_ci), width = .1) +
  scale_color_manual("Sex",
                     values = phs_colours(c("phs-green", "phs-blue"))) + 
  facet_wrap(~ incidence_type) + 
  ylab("Age-Standardised Rate per 100,000") + 
  xlab("Year") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 75)) +
  theme(plot.title = element_text(hjust = 0.5, size = 40),
        axis.text.x = element_text(vjust = 0.5, size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 35),
        strip.text.x = element_text(size = 25), 
        legend.title = element_text(size = 30),
        legend.key.size = unit(1.5, "cm"), 
        legend.text = element_text(size = 30))

# Save chart

ggsave(here(glue("Charts/{end}/incidence_rate_sex_non_emergency_plot.png")),
       plot = inc_rate_sex_el_plot,
       height = 14,
       width = 23,
       dpi = 500)


# Same as above but with linear regression

inc_rate_sex_el_plot_reg <- inc_rate_site_sex_output %>% 
  filter(emergency_flag == "Non-Emergency" & 
           incidence_type %in% c("Colorectal", "Head and Neck", "Lung", 
                                 "Upper GI")) %>% 
  mutate(year = as.character(year), 
         sex = recode(sex, "1" = "Male", "2" = "Female")) %>% 
  ggplot(aes(x = year, group = sex, colour = sex)) +
  geom_line(aes(y = obs_inc, colour = sex), size = 2) + 
  geom_line(aes(y = exp_inc), size = 0.5, linetype = "dashed") + 
  geom_errorbar(aes(ymin = obs_l_ci, ymax = obs_u_ci), width = .1) +
  scale_color_manual("Sex",
                     values = phs_colours(c("phs-green", "phs-blue"))) + 
  facet_wrap(~ incidence_type) + 
  ylab("Age-Standardised Rate per 100,000") + 
  xlab("Year") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 75)) +
  theme(plot.title = element_text(hjust = 0.5, size = 40),
        axis.text.x = element_text(vjust = 0.5, size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 35),
        strip.text.x = element_text(size = 25), 
        legend.title = element_text(size = 30),
        legend.key.size = unit(1.5, "cm"), 
        legend.text = element_text(size = 30))

# Save chart

ggsave(here(glue("Charts/{end}/incidence_rate_sex_non_emergency_regression_plot.png")),
       plot = inc_rate_sex_el_plot_reg,
       height = 14,
       width = 23,
       dpi = 500)


# Same as above but for emergency presentations

inc_rate_sex_em_plot_reg <- inc_rate_site_sex_output %>% 
  filter(emergency_flag == "Emergency" & 
           incidence_type %in% c("Colorectal", "Head and Neck", "Lung", 
                                 "Upper GI")) %>% 
  mutate(year = as.character(year), 
         sex = recode(sex, "1" = "Male", "2" = "Female")) %>% 
  ggplot(aes(x = year, group = sex, colour = sex)) +
  geom_line(aes(y = obs_inc, colour = sex), size = 2) + 
  geom_line(aes(y = exp_inc), size = 0.5, linetype = "dashed") + 
  geom_errorbar(aes(ymin = obs_l_ci, ymax = obs_u_ci), width = .1) +
  scale_color_manual("Sex",
                     values = phs_colours(c("phs-green", "phs-blue"))) + 
  facet_wrap(~ incidence_type) + 
  ylab("Age-Standardised Rate per 100,000") + 
  xlab("Year") + 
  # ggtitle("Age-Standardised Incidence Rate per 100,000 for Cancers with an
  #         Emergency Route to Diagnosis by Sex") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 62)) +
  theme(plot.title = element_text(hjust = 0.5, size = 40),
        axis.text.x = element_text(vjust = 0.5, size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 35),
        strip.text.x = element_text(size = 25), 
        legend.title = element_text(size = 30),
        legend.key.size = unit(1.5, "cm"), 
        legend.text = element_text(size = 30))

# Save chart

ggsave(here(glue("Charts/{end}/incidence_rate_sex_emergency_regression_plot.png")),
       plot = inc_rate_sex_em_plot_reg,
       height = 14,
       width = 23,
       dpi = 500)



### 3.3 Cancer Site and SIMD Split ----

# Create chart for non-emergency lung presentations by SIMD quintile
# Define colours based on phsstyles package
# Add axis labels and set y axis to start from 0 and go to specified length

inc_rate_simd_el_plot <- create_chart(inc_rate_site_simd_output %>% 
                                        filter(incidence_type == "Lung"), 
                                      "Non-Emergency", "simd2020v2_sc_quintile") + 
  scale_colour_manual("SIMD Quintile",
                      values = phs_colours(c("phs-purple", "phs-magenta",
                                             "phs-blue", "phs-green", 
                                             "phs-teal"))) +
  ylab("Age-Sex Standardised Rate per 100,000") + 
  xlab("Year") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 126))

# Save chart

ggsave(here(glue("Charts/{end}/incidence_rate_simd_non_emergency_plot.png")),
       plot = inc_rate_simd_el_plot,
       height = 14,
       width = 23,
       dpi = 500)

# Same as above but with linear regression

inc_rate_simd_el_plot_reg <- create_regression_chart(
  inc_rate_site_simd_output %>% 
    filter(incidence_type == "Lung"), 
  "Non-Emergency", "simd2020v2_sc_quintile") + 
  scale_colour_manual("SIMD Quintile",
                      values = phs_colours(c("phs-purple", "phs-magenta",
                                             "phs-blue", "phs-green", 
                                             "phs-teal"))) +
  ylab("Age-Sex Standardised Rate per 100,000") + 
  xlab("Year") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 125))

# Save chart

ggsave(here(glue("Charts/{end}/incidence_rate_simd_non_emergency_regression_plot.png")),
       plot = inc_rate_simd_el_plot_reg,
       height = 14,
       width = 23,
       dpi = 500)



### 4 Output ----

### 4.1 Main Output ----

# Create a workbook

wb <- createWorkbook()

# Define a header style for workbook

hs <- createStyle(fontColour = "#ffffff", fgFill = "#0078D4",
                  halign = "center", valign = "center", 
                  textDecoration = "bold", border = "TopBottomLeftRight")

# Add sheet for Site Incidence Rates and write inc_rate_site_output data
# Round relevant columns to two decimal places
# Select and rename relevant columns
# Set column widths to auto for this sheet

addWorksheet(wb, "Site Incidence Rates", gridLines = FALSE)

writeData(wb, sheet = "Site Incidence Rates", 
          inc_rate_site_output %>% 
            mutate(across(c(obs_inc, obs_l_ci, obs_u_ci, 
                            exp_inc, exp_l_ci, exp_u_ci, 
                            sir, sir_l_ci, sir_u_ci, sir_p_value, 
                            srr, srr_l_ci, srr_u_ci, 
                            exp_n),
                          ~ round_half_up(., 2))) %>%
            select(Year = year, 
                   Type = emergency_flag, 
                   "Cancer Type" = incidence_type, 
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
          borders = "all", headerStyle = hs)

setColWidths(wb, sheet = "Site Incidence Rates", cols = 1:18, widths = "auto")

# Add sheet for Sex Incidence Rates and write inc_rate_site_sex_output data
# Round relevant columns to two decimal places
# Select and rename relevant columns
# Set column widths to auto for this sheet

addWorksheet(wb, "Sex Incidence Rates", gridLines = FALSE)

writeData(wb, sheet = "Sex Incidence Rates", 
          inc_rate_site_sex_output %>% 
            mutate(across(c(obs_inc, obs_l_ci, obs_u_ci, 
                            exp_inc, exp_l_ci, exp_u_ci, 
                            sir, sir_l_ci, sir_u_ci, sir_p_value, 
                            srr, srr_l_ci, srr_u_ci, 
                            exp_n),
                          ~ round_half_up(., 2))) %>%
            mutate(sex = recode(sex, "1" = "Male", "2" = "Female")) %>% 
            select(Year = year, 
                   Type = emergency_flag, 
                   "Cancer Type" = incidence_type, 
                   Sex = sex, 
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

setColWidths(wb, sheet = "Sex Incidence Rates", cols = 1:19, widths = "auto")

# Add sheet for SIMD Incidence Rates and write inc_rate_simd_output data
# Round relevant columns to two decimal places
# Select and rename relevant columns
# Set column widths to auto for this sheet

addWorksheet(wb, "SIMD Incidence Rates", gridLines = FALSE)

writeData(wb, sheet = "SIMD Incidence Rates", 
          inc_rate_site_simd_output %>% 
            mutate(across(c(obs_inc, obs_l_ci, obs_u_ci, 
                            exp_inc, exp_l_ci, exp_u_ci, 
                            sir, sir_l_ci, sir_u_ci, sir_p_value, 
                            srr, srr_l_ci, srr_u_ci, 
                            exp_n),
                          ~ round_half_up(., 2))) %>%
            select(Year = year, 
                   Type = emergency_flag, 
                   "Cancer Type" = incidence_type, 
                   SIMD = simd2020v2_sc_quintile, 
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

setColWidths(wb, sheet = "SIMD Incidence Rates", cols = 1:19, widths = "auto")

# Save workbook

saveWorkbook(wb, here(glue("Data/{end}/incidence_rates_{start}_{end}.xlsx")), 
             overwrite = TRUE)



### 4.2 Tables for Report ----

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
          inc_rate_site_output %>% 
            filter(year %in% c(2020, 2021)) %>% 
            mutate(across(c(obs_inc, obs_l_ci, obs_u_ci, 
                            exp_inc, exp_l_ci, exp_u_ci, 
                            sir, sir_l_ci, sir_u_ci, sir_p_value, 
                            srr, srr_l_ci, srr_u_ci, 
                            exp_n),
                          ~ round_half_up(., 2))) %>%
            mutate(sir_ci = paste0("(", sir_l_ci, ", ", sir_u_ci, ")")) %>% 
            select(Year = year, 
                   Type = emergency_flag, 
                   "Cancer Type" = incidence_type, 
                   "Observed Number" = obs_n, 
                   "Expected Number" = exp_n, 
                   "SIR" = sir, 
                   "SIR CI" = sir_ci, 
                   "SIR P Value" = sir_p_value), 
          borders = "all", headerStyle = hs)

setColWidths(report_wb, sheet = "SIR Data", cols = 1:8, widths = "auto")

# Add sheet for Sex Incidence Rates and write inc_rate_site_sex_output data
# Round relevant columns to two decimal places
# Select and rename relevant columns
# Set column widths to auto for this sheet

addWorksheet(report_wb, "SRR Data", gridLines = FALSE)

writeData(report_wb, sheet = "SRR Data", 
          inc_rate_site_output %>% 
            filter(year %in% c(2020, 2021)) %>% 
            mutate(across(c(obs_inc, obs_l_ci, obs_u_ci, 
                            exp_inc, exp_l_ci, exp_u_ci, 
                            sir, sir_l_ci, sir_u_ci, 
                            srr, srr_l_ci, srr_u_ci, 
                            exp_n),
                          ~ round_half_up(., 2))) %>%
            mutate(srr_ci = paste0("(", srr_l_ci, ", ", srr_u_ci, ")")) %>% 
            select(Year = year, 
                   Type = emergency_flag, 
                   "Cancer Type" = incidence_type, 
                   "Observed Rate" = obs_inc, 
                   "Expected Rate" = exp_inc, 
                   "SRR" = srr, 
                   "SRR CI" = srr_ci), 
          borders = "all", headerStyle = hs, startCol = 1, startRow = 1)

setColWidths(report_wb, sheet = "SRR Data", cols = 1:7, widths = "auto")


# Save workbook

saveWorkbook(report_wb, here(glue("Data/{end}/report_tables.xlsx")), 
             overwrite = TRUE)
