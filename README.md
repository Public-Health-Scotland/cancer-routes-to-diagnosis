# cancer-routes-to-diagnosis
Emergency and Non-Emergency Routes to Cancer Diagnosis

This repository contains code for analysing emergency and non-emergency routes to cancer diagnosis in Scotland. This work involves data linkage between SMR01 (national inpatient dataset) and SMR06 (Scottish Cancer Registry) to identify emegrency and non-emergency routes to cancer diagnosis. There are five R scripts used within this analysis:

- `1_housekeeping` sets up the data extracts from SMR01 and SMR06, joins them together and defines dates and years to use for extracts and analysis. This script also contains various functions used throughout the other scripts for things like fitting Poisson regression models and calculating incidence rates. The joined extract defined within this script is used throughout the other scripts.
- `2_site_analysis` produced counts for number of diagnosis by route, year and cancer site.
- `3_incidence_rates` is the largest script and produces age-sex standardised incidence rates and ratios by route, year, cancer site, sex and deprivation (SIMD).
- `4_method_1st_detection` produces age-sex standardised incidence rates and ratios by route, year, cancer site and method of 1st detection.
- `5_stage_at_diagnosis` produces age-sex standardised incidence rates and ratios by route, year, cancer site and stage at diagnosis.

If you have any questions about this analysis or wish further information please contact [calum.purdie@phs.scot](mailto:calum.purdie@phs.scot)
