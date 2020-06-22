#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
VOC2PM_output4paper.py
    a python script to calculate VOC2PM ratios by smoke age for paper,
    adapted from OG version of this code "young_VOC2PM_data4Kat.py"
    written to calcualte VOC2PM ratios in young smoke for Kat Navarro.
Created on Thu Dec 19 10:46:54 2019
V2 created Thu Jan 7 2020
@author: kodell
"""
#%% user inputs
# this list is the full list of HAPs measured during WE-CAN
# file name for table of VOC names
HAPs_RFs_fn = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/health risk tables/HAPs_EFs_RFs_forcode_final_final_WadeUpdate.csv'

# specify which TOGA merge version to use
# for paper
mrg_name = 'TOGA_mrg_w_age__NASA_R2_anth_tracers_PTR16fix_R1TOGAupdate85.275_w200_CH3CN.csv'
voc_bk_name = 'VOC_bks__NASA_R2_anth_tracers_PTR16fix_R1TOGAupdate85.275_w200_CH3CN.csv'
out_fn = 'VOC2PM_ratios_4paper_NASAmrg_R1TOGAupdate.csv'

# detection limit file name
TOGA_LLOD_fn =  '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/WECAN_datafiles/TOGA_data_05.21.20/TOGA_LLOD_R1.csv'

# sensitivty tests
#mrg_name = 'TOGA_mrg_w_age__NASA_R2_anth_tracers_bk_nagebk_minus1sig85.275_w200_CH3CN.csv'
#voc_bk_name = 'VOC_bks__NASA_R2_anth_tracers_bk_nagebk_minus1sig85.275_w200_CH3CN.csv'
#out_fn = 'VOC2PM_ratios_4paper_NASAmrg_allbk_minus1sig.csv'

#out_fn = 'VOC2PM_ratios_4paper_NASAmrg_addPTRcheck_dl4bdl.csv'
#mrg_name = 'TOGA_mrg_w_age__NASA_R2_anth_tracers_bk_addPTRcheck_dl4bdl85.275_w200_CH3CN.csv'
#voc_bk_name = 'VOC_bks__NASA_R2_anth_tracers_bk_addPTRcheck_dl4bdl85.275_w200_CH3CN.csv'

# list of age names to loop through later
age_names = ['not smoke', 'young, < 0.5 days','med, 0.5-1.5 days', 'old, < 4 days','> 4 days']

#%% load modules
import pandas as pd
import numpy as np
import plotly.io as pio
pio.renderers.default = "chrome"

#%% load data
# TOGA merge
mrg_fn = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/data analysis/sensitivity_tests/sensitivity_data/' + mrg_name
TOGA_mrg = pd.read_csv(mrg_fn)

# detection limit file
TOGA_LLOD = pd.read_csv(TOGA_LLOD_fn)

# background estimates
bk_fn = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/data analysis/sensitivity_tests/sensitivity_data/' + voc_bk_name
VOC_bks = pd.read_csv(bk_fn)

# VOC names and chemical formulas
HAPs_RFs_raw = pd.read_csv(HAPs_RFs_fn,header=0,skiprows=[1,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,
                                                          49,50,51])
HAPs_RFs = HAPs_RFs_raw.set_index('TOGA name')
healthVOCs = HAPs_RFs.index

#%% create dataframe for ratios and R2
VOC_PM_ratio_1 = pd.DataFrame(data = {'variable name':HAPs_RFs.index})
VOC_PM_ratio = VOC_PM_ratio_1.set_index('variable name')
VOC_PM_ratio['VOC name'] = HAPs_RFs['Name']
VOC_PM_ratio['Chemical Formula'] = HAPs_RFs['Chemical Formula']
    

#%% create VOC_elv dataframe by subtracting off background
VOC_smoke_elv = pd.DataFrame(data={'  STARTTIME' : TOGA_mrg['  STARTTIME'],
                                   ' STOPTIME': TOGA_mrg[' STOPTIME'],
                                   ' FLIGHT':TOGA_mrg[' FLIGHT']})
# add PM and CO
VOC_smoke_elv['PM_elv'] = TOGA_mrg['PM1_ams_sp2'].values - VOC_bks['PM1_ams_sp2'].values
VOC_smoke_elv['CO_comb_elv'] = TOGA_mrg['CO_comb'].values - VOC_bks['CO_comb'].values

for VOC_name in healthVOCs:
    VOC_smoke_elv[' '+VOC_name + '_elv'] = TOGA_mrg[' '+VOC_name].values - VOC_bks[VOC_name].values
# if elv < 0, make zero.
# test not doing this to see if main conclusions will change. 
negs = np.where(VOC_smoke_elv <0.)
for i in range(len(negs[0])):
    VOC_smoke_elv.iloc[negs[0][i],negs[1][i]]=0
# add chem age
VOC_smoke_elv['chem_smoke_age'] = TOGA_mrg['chem_smoke_age']
# remove urban points
anth_inds = np.where(TOGA_mrg['anth_flag']>0)[0]
VOC_smoke_elv_nourban = VOC_smoke_elv.drop(index=anth_inds)
# drop does not reset index, so we need to reset this
VOC_smoke_elv_nourban.reset_index(inplace=True,drop=True)
TOGA_mrg_nourban = TOGA_mrg.drop(index=anth_inds)
TOGA_mrg_nourban.reset_index(inplace=True,drop=True)
#%% group smoke elevated concentrations by age
smoke_elv_age_groups = VOC_smoke_elv_nourban.groupby(by='chem_smoke_age')
TOGA_mrg_agegroups = TOGA_mrg_nourban.groupby(by='chem_smoke_age')
#%% stats for smoke
ai = 1
for age in ['young', 'medium', 'old', 'extra old']:
    VOC_PM_ratio[age+' VOC v PM R2'] = -999.
    VOC_PM_ratio[age+' VOC PM1 ratio, median [ppt/ugm-3]'] = -999.
    VOC_PM_ratio[age+' VOC PM1 ratio, mean [ppt/ugm-3]'] = -999.
    VOC_PM_ratio[age+' VOC PM1 ratio, standard deviation'] = -999.
    VOC_PM_ratio[age+' VOC PM1 ratio, n obs'] = -999.
    VOC_PM_ratio[age+' VOC PM1 ratio, % obs above DL'] = -999.
    
    TOGA_mrg_agegroup = TOGA_mrg_agegroups.get_group(age_names[ai])
    VOC_smoke_elv_agegroup = smoke_elv_age_groups.get_group(age_names[ai])
    VOC_PM_weighted_agegroup = VOC_smoke_elv_agegroup[VOC_smoke_elv_agegroup.columns[:4]].copy(deep=True)
    zeroPM = np.where(VOC_smoke_elv_agegroup['PM_elv']==0)
    nanPM = np.where(np.isnan(VOC_smoke_elv_agegroup['PM_elv']))
    rmv_inds = np.hstack([zeroPM,nanPM])
    i = 0
    for voc in VOC_PM_ratio.index:
        # R2
        data = VOC_smoke_elv_agegroup[[' '+voc + '_elv','PM_elv']]
        VOC_PM_ratio[age+' VOC v PM R2'].iloc[i] = round(data.corr('pearson').values[0,1]**2.,2)
        # median, mean, standard deviation
        VOC_PM_weighted_agegroup[' '+voc+'_elv'] = VOC_smoke_elv_agegroup[' '+voc + '_elv'].values/VOC_smoke_elv_agegroup['PM_elv'].values
        VOC_PM_weighted_agegroup[' '+voc+'_elv'].iloc[rmv_inds] = np.nan
        # for calculation above DL
        TOGA_mrg_agegroup[' '+voc].iloc[rmv_inds] = np.nan
        VOC_PM_ratio[age+' VOC PM1 ratio, median [ppt/ugm-3]'].iloc[i] = np.round(np.nanmedian(VOC_PM_weighted_agegroup[' '+voc+'_elv']),2)
        VOC_PM_ratio[age+' VOC PM1 ratio, mean [ppt/ugm-3]'].iloc[i] = np.round(np.nanmean(VOC_PM_weighted_agegroup[' '+voc+'_elv']),2)
        VOC_PM_ratio[age+' VOC PM1 ratio, standard deviation'].iloc[i] = np.round(np.nanstd(VOC_PM_weighted_agegroup[' '+voc+'_elv'],ddof=1),2)
        # n
        nobs = len(np.where(np.isfinite(VOC_PM_weighted_agegroup[' '+voc+'_elv']))[0])
        VOC_PM_ratio[age+' VOC PM1 ratio, n obs'].iloc[i] = nobs
        if voc[-4:] == 'TOGA':
            nobs_above_DL =  len(np.where(TOGA_mrg_agegroup[' '+voc]>TOGA_LLOD[voc].values[0])[0])
        else:
            nobs_above_DL =  len(np.where(TOGA_mrg_agegroup[' '+voc]>200.0)[0])
    
        VOC_PM_ratio[age+' VOC PM1 ratio, % obs above DL'].iloc[i] = 100.0*(np.float(nobs_above_DL)/nobs)
        print nobs_above_DL
        i += 1
        
    ai += 1
#%% save
VOC_PM_ratio.to_csv('/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/health risk tables/'+out_fn)

#%% make figure
'''
# test case for young smoke, only named species in the box plots
VOCs_plot = ['Benzene','Formaldehyde','Acetaldehyde','Acrylonitrile (2Propenenitrile)',
                 'Acrolein','Hydrogen Cyanide']
VOC_PM_ratio_plot = VOC_PM_ratio.set_index('VOC name')
VOC_PM_ratio_subset = VOC_PM_ratio_plot.loc[VOCs_plot]
fig = plotly.subplots.make_subplots(rows=2, cols=2,
                                    subplot_titles = ['young','medium','old','older'])
fig.add_trace(go.Scatter(x=VOC_PM_ratio_subset['young VOC v PM R2'],
                 y=VOC_PM_ratio_subset['young VOC PM1 ratio [ppt/ugm-3]'],
                 mode='markers+text',
                 marker_color='purple',
                 text = VOC_PM_ratio_subset.index,
                 textposition='top center'),1,1)
fig.add_trace(go.Scatter(x=VOC_PM_ratio_subset['medium VOC v PM R2'],
                 y=VOC_PM_ratio_subset['medium VOC PM1 ratio [ppt/ugm-3]'],
                 mode='markers+text',
                 marker_color='red',
                 text = VOC_PM_ratio_subset.index,
                 textposition='top center'),1,2)
fig.add_trace(go.Scatter(x=VOC_PM_ratio_subset['old VOC v PM R2'],
                 y=VOC_PM_ratio_subset['old VOC PM1 ratio [ppt/ugm-3]'],
                 mode='markers+text',
                 marker_color='orange',
                 text = VOC_PM_ratio_subset.index,
                 textposition='top center'),2,1)
fig.add_trace(go.Scatter(x=VOC_PM_ratio_subset['>1 week VOC v R2'],
                 y=VOC_PM_ratio_subset['>1week VOC PM1 ratio [ppt/ugm-3]'],
                 mode='markers+text',
                 marker_color='green',
                 text = VOC_PM_ratio_subset.index,
                 textposition='top center'),2,2)
fig.show()
'''





