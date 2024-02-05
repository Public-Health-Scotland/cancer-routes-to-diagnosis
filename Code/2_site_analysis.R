#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2_site_analysis
# Calum Purdie
# 03/11/2021
# Calculates totals and proportions by cancer site
# Written/run on Posit Workbench
# R version 4.1.2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



### 1 Housekeeping ----

source(here::here("Code/1_housekeeping.R"))



### 2 Cancer Sites ----

# Count by incidence_year, incidence_type and cancer_site to get total cancers

site_totals <- joined_data |>   
  count(incidence_year, incidence_type, name = "Cancers")

# Count by incidence_year, incidence_type and emergency_flag to get total 
# presentations

site_admissions <- joined_data |> 
  count(incidence_year, incidence_type, emergency_flag, 
        name = "Presentations")

# Join totals and admissions and set any NAs to 0 for numeric columns
# Calculate percentage cancers which had an emergency presentation
# Calculate 95% confidence intervals for PC_Presentations and Cancers
# Round PC_Presentations, Lower_95_CI and Upper_95_CI to 1 decimal place and 
# multiply by 100 to get percentages
# Arrange data and select columns

site_output <- full_join(site_totals, site_admissions) |> 
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .))) |> 
  mutate(PC_Presentations = Presentations / Cancers,
         Lower_95_CI = Wilson_lowerCI(PC_Presentations, 0.05, Cancers), 
         Upper_95_CI = Wilson_upperCI(PC_Presentations, 0.05, Cancers)) |> 
  mutate(across(c(PC_Presentations, Lower_95_CI, Upper_95_CI), 
                ~ round_half_up(.x * 100, 1))) |> 
  arrange(incidence_type, incidence_year) |> 
  select(incidence_year, incidence_type, emergency_flag, Cancers, Presentations, 
         PC_Presentations, Lower_95_CI, Upper_95_CI)

# Get Scotland level totals
# Group by incidence_year and emergency_flag and sum Cancers and Presentations
# Calculate percentage presentations and 95% confidence intervals
# Round PC_Presentations, Lower_95_CI and Upper_95_CI to 1 decimal place and 
# multiply by 100 to get percentages

scotland_output <- site_output |> 
  group_by(incidence_year, emergency_flag) |> 
  summarise(Cancers = sum(Cancers), 
            Presentations = sum(Presentations)) |> 
  mutate(PC_Presentations = Presentations / Cancers,
         Lower_95_CI = Wilson_lowerCI(PC_Presentations, 0.05, Cancers), 
         Upper_95_CI = Wilson_upperCI(PC_Presentations, 0.05, Cancers)) |> 
  mutate(across(c(PC_Presentations, Lower_95_CI, Upper_95_CI), 
                ~ round_half_up(.x * 100, 1)))



### 3 Charts ----

# Create chart for non-emergency presentations by cancer site
# Filter for emergency presentations
# Select columns and define year as character field
# Define ggplot aesthetics and create a line chart with one line for each cancer
# Define colours based on phsstyles package
# Add axis labels and set y axis to start from 0 and go to specified length
# Add theme details to adjust text - this will make the chart look odd in R but
# looks fine once saved out

n_cancers_em_plot <- site_output |> 
  filter(emergency_flag == "Emergency") |> 
  select(incidence_year, incidence_type, Presentations) |> 
  mutate(incidence_year = as.character(incidence_year)) |> 
  ggplot(aes(x = incidence_year, y = Presentations, group = incidence_type)) +
  geom_line(aes(colour = incidence_type), size = 2) + 
  scale_colour_manual("Cancer",
                      values = phs_colours(c("phs-graphite", "phs-magenta",
                                             "phs-blue", "phs-green", 
                                             "phs-rust", "phs-teal", 
                                             "phs-purple"))) +
  ylab("Number of Emergency Cancer Presentations") + 
  xlab("Year") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2400)) +
  theme(plot.title = element_text(hjust = 0.5, size = 40),
        axis.text.x = element_text(vjust = 0.5, size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 35),
        legend.title = element_text(size = 30),
        legend.key.size = unit(1.5, "cm"), 
        legend.text = element_text(size = 30))

# Save chart

ggsave(here(glue("Charts/{end}/number_cancers_emergency_line_plot.png")),
       plot = n_cancers_em_plot,
       height = 14,
       width = 23,
       dpi = 500)

# Create chart for non-emergency presentations by cancer site
# Filter for non-emergency presentations
# Select columns and define year as character field
# Define ggplot aesthetics and create a line chart with one line for each cancer
# Define colours based on phsstyles package
# Add axis labels and set y axis to start from 0 and go to specified length
# Add theme details to adjust text - this will make the chart look odd in R but
# looks fine once saved out

n_cancers_el_plot <- site_output |> 
  filter(emergency_flag == "Non-Emergency") |> 
  select(incidence_year, incidence_type, Presentations) |> 
  mutate(incidence_year = as.character(incidence_year)) |> 
  ggplot(aes(x = incidence_year, y = Presentations, group = incidence_type)) +
  geom_line(aes(colour = incidence_type), size = 2) + 
  scale_colour_manual("Cancer",
                      values = phs_colours(c("phs-graphite", "phs-magenta",
                                             "phs-blue", "phs-green", 
                                             "phs-rust", "phs-teal", 
                                             "phs-purple"))) +
  ylab("Number of Non-Emergency Cancer Presentations") + 
  xlab("Year") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 5600)) +
  theme(plot.title = element_text(hjust = 0.5, size = 40),
        axis.text.x = element_text(vjust = 0.5, size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 35),
        legend.title = element_text(size = 30),
        legend.key.size = unit(1.5, "cm"), 
        legend.text = element_text(size = 30))

# Save chart

ggsave(here(glue("Charts/{end}/number_cancers_non_emergency_line_plot.png")),
       plot = n_cancers_el_plot,
       height = 14,
       width = 23,
       dpi = 500)

# Create chart for percentage presentations by year
# Define year as character field
# Define ggplot aesthetics and create a line chart with one line for each cancer
# Define colours based on phsstyles package
# Add axis labels and set y axis to start from 0 and go to specified length
# Add theme details to adjust text - this will make the chart look odd in R but
# looks fine once saved out

pc_cancers_scotland_plot <- scotland_output |> 
  mutate(incidence_year = as.character(incidence_year)) |> 
  ggplot(aes(x = incidence_year, y = PC_Presentations, fill = emergency_flag)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_errorbar(aes(ymin = Lower_95_CI, ymax = Upper_95_CI),
                width = .2,
                position = position_dodge(.9)) +
  scale_fill_manual("Presentation Type",
                    values = phs_colours(c("phs-green", "phs-blue"))) +
  
  ylab("Proportion of Cancer Presentations") + 
  xlab("Year") + 
  theme(plot.title = element_text(hjust = 0.5, size = 40),
        axis.text.x = element_text(vjust = 0.5, size = 30),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 35),
        legend.title = element_text(size = 30),
        legend.key.size = unit(1.5, "cm"), 
        legend.text = element_text(size = 30))

# Save chart

ggsave(here(glue("Charts/{end}/pc_cancers_scotland_bar_plot.png")),
       plot = pc_cancers_scotland_plot,
       height = 14,
       width = 23,
       dpi = 500)



### 4 Output ----

# Create a workbook

wb <- createWorkbook()

# Define a header style for workbook

hs <- createStyle(fontColour = "#ffffff", fgFill = "#0078D4",
                  halign = "center", valign = "center", 
                  textDecoration = "bold", border = "TopBottomLeftRight")

# Add sheet for Scotland Data and write scotland_output data
# Select and rename relevant columns
# Set column widths to auto for this sheet

addWorksheet(wb, "Scotland Data", gridLines = FALSE)

writeData(wb, sheet = "Scotland Data", 
          scotland_output |> 
            select(Year = incidence_year, 
                   "Presentation Type" = emergency_flag, 
                   Cancers, 
                   Presentations, 
                   "Percentage Presentations" = PC_Presentations, 
                   "Lower 95 CI" = Lower_95_CI, 
                   "Upper 95 CI" = Upper_95_CI), 
          borders = "all", headerStyle = hs)

setColWidths(wb, sheet = "Scotland Data", cols = 1:7, widths = "auto")

# Add sheet for Scotland Data and write scotland_output data
# Select and rename relevant columns
# Set column widths to auto for this sheet

addWorksheet(wb, "Site Totals", gridLines = FALSE)

writeData(wb, sheet = "Site Totals", 
          site_output |> 
            select(Year = incidence_year, 
                   Type = emergency_flag, 
                   "Cancer Type" = incidence_type, 
                   Cancers, 
                   Presentations, 
                   "Percentage Presentations" = PC_Presentations, 
                   "Lower 95 CI" = Lower_95_CI, 
                   "Upper 95 CI" = Upper_95_CI), 
          borders = "all", headerStyle = hs)

setColWidths(wb, sheet = "Site Totals", cols = 1:8, widths = "auto")

# Save workbook

saveWorkbook(wb, 
             here(glue("Data/{end}/emergency_presentations_{start}_{end}.xlsx")), 
             overwrite = TRUE)
