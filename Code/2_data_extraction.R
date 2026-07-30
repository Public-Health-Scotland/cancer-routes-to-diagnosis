#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2_data_extraction
# 03/11/2021
# Sets up initial joined SMR01/SMR06 extract for analysis
# Written/run on Posit Workbench
# R version 4.1.2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


### 1 Housekeeping ----

source(here::here("Code/1_housekeeping.R"))



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
         incidence_type = case_when(cancer_site %in% c("C00", "C01", "C02", 
                                                       "C03", "C04", "C05", 
                                                       "C06", "C07", "C08", 
                                                       "C09", "C10", "C11", 
                                                       "C12", "C13", "C14", 
                                                       "C30", "C31", "C32") ~ 
                                      "Head and Neck", 
                                    cancer_site %in% c("C15") ~ "Oesophagus", 
                                    cancer_site %in% c("C16") ~ "Stomach", 
                                    cancer_site %in% c("C18", "C19", "C20") ~ 
                                      "Colorectal", 
                                    cancer_site %in% c("C33", "C34") ~ "Lung",
                                    cancer_site %in% c("C43") ~ "Melanoma", 
                                    cancer_site %in% c("C50") ~ "Breast", 
                                    cancer_site %in% c("C53") ~ "Cervical", 
                                    cancer_site %in% c("C54") ~ "Corpus Uteri", 
                                    cancer_site %in% c("C56") ~ "Ovary", 
                                    cancer_site %in% c("C61") ~ "Prostate", 
                                    cancer_site %in% c("C64", "C65") ~ "Kidney",
                                    cancer_site %in% c("C67") ~ "Bladder",
                                    cancer_site %in% c("C70", "C71", "C72") |
                                      icd10s_cancer_site %in% c("C751", "C752", "C753") ~ "Brain", 
                                    cancer_site %in% c("C82", "C83", "C84", "C85") ~ "Non-Hodgkin's Lymphoma", 
                                    cancer_site %in% c("C91", "C92", "C93", "C94", "C95") ~ "Leukaemia"),
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
                        "postcode_2025_1_simd2020v2.rds")) |>  
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
                         "DataZone2011_pop_est_5year_agegroups_2011_2022.rds")) |> 
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
