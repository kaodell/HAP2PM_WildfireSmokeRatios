##README - HAP2PM_WildfireSmokeRatios
# created by: Katelyn O'Dell
# date: 12/19/19
# contact: katelyn.odell 'at' colostate.edu

# This git repository contains python v2.7 scripts used to create a TOGA merge, estimate chemical 
# smoke age, anthropogenic influence, and HAP to VOC ratios in the WE-CAN campaign.

# The paper resulting from this analysis is currently in prep for ES&T. 

# Details and data from the WE-CAN campaign are publically available at the following link:
# https://www.eol.ucar.edu/field_projects/we-can

# These python scripts require data from the following instruments onboard the C130 during
# WE-CAN: PTR-MS (special data files of specific VOCs from Wade Permer), TOGA R2 merge, and 
# updated TOGA formaldehyde values.

# Details on each code:
# find_anthropogenic_tracers_NASAmrg.py - input TOGA species to use as anthropogenic tracers, identify and flag anthropogeic-influenced observations
# VOC-PM_TOGA_sensitivity_tests_NASAmrg.py - input CO, HCN, and CH3CN limits to use as smoke cut-offs, calculate chemical smoke age,       identify smoke and smoke-free observtions, calcualte non-smoke backgrounds of HAPs, place PTR data on TOGA timescale
# WE-CAN_paper_figures_final_NASAmrg.py - calculate HAP/PM ratios, HAPs-PM weighted risk, and create figures 1-3 from the paper
# VOC2PM_output4paper.py - calcualte HAP/PM ratios for each age category and create supplementary csv file from paper
# PM_and_VOC_maps4paper.py - calculate chronic HAPs exposure and risk and create figure 4 from paper
# check_age_calc.py - create figures S1.1 and S1.2
# check_age_distriubtions_wdilution.py - create figure S1.3

# To perform analysis presented in paper:
# 1) download TOGA merge from WE-CAN site and special PTR data from:
# 2) download table from paper appendix of HAPs names, molecular weights, reference concentrations, and unit risk estimates
# 2) run find_anthropogenic_tracers_NASAmrg.py to find and flag data with urban influence
# 3) run VOC-PM_TOGA_sensitivity_tests_NASAmrg.py to identify smoke-influenced points, calculate non-smoke mixing ratios, calculate       chemical smoke age and place PTR variables on TOGA timescale
# 4) run WE-CAN_paper_figures_final_NASAmrg.py to weight HAPs by PM and risk and create figures 1-3 and S2.1 in the paper
# 5) run VOC2PM_output4paper.py to calculte HAP/PM ratios from paper supplentary files
# 6) download kriged PM2.5 from Colorado State University Repositories, DOI(s): xxxxxx
# 7) run PM_and_VOC_maps4paper.py to calculate chronic HAPs exposure and risk, and create figure 4 in the paper
# ... for supplemental figures
# 8) run check_age_calc.py to create figures S1.1 and S1.2
# 9) run check_age_distributions_wdiltuion.py to create figure S1.3

# Version Details

# 12.19.19 - initial comit with code used to create data files of PM2VOC relationships
#	in young smoke for Kat Navarro. 

# 06.15.20 - final codes for paper analysis uploaded
