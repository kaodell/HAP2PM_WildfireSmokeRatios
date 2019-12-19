#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
young_VOC2PM_data4Kat.py
    a python script to calcualte VOC2PM ratios in young smoke for Kat Navarro.
    'young' in this case is ~< 1 day of aging
Created on Thu Dec 19 10:46:54 2019
@author: kodell
"""
#%% user inputs
# list VOC variable names you would like to calcualte relationships for
# this list is the full list of HAPs measured during WE-CAN
VOC_names = ['CH3CCl3_TOGA','CH2ClCH2Cl_TOGA','x224TrimePentane_TOGA','Acetaldehyde_TOGA',
            'CH3CN_TOGA','Acrolein_TOGA','Acrylonitrile_TOGA','Benzene_TOGA','CHBr3_TOGA',
            'CH3Br_TOGA','CS2_TOGA','ClBenzene_TOGA','CH3Cl_TOGA','CH2Cl2_TOGA','EtBenzene_TOGA',
            'CH2O_TOGA','HCN_TOGA','CH3I_TOGA','CH3OH_TOGA','nHexane_TOGA','oXylene_TOGA',
            'Propanenitrile_TOGA','Propanal_TOGA','Styrene_TOGA','C2Cl4_TOGA','Toluene_TOGA',
            'CHCl3_TOGA','mpXylene_TOGA','MeAcrylonitrile_TOGA',
            'ACETAMIDE_C2H5NO_PTR','METHYL_METHACRYLATE_C5H8O2_PTR','PHENOL_C6H6O_PTR',
            'NITROBENZENE_C6H5NO2_PTR','NAPHTHALENE_C10H8_PTR',
            '1-3-BUTADIENE_1-2-BUTADIENE_C4H6_PTR',
            'METHYL-ISOCYANATE_HYDROXY-ACETONITRILE_C2H3NO_PTR',
            'QUINONE_C6H4O2_PTR','ISOCYANIC_ACID_HNCO_PTR']

# specify CO and HCN cutoffs for age ... this indicates which files to load
CO_cutoff = '85' #ppb
HCN_cutoff = '275' #ppt

# specify which TOGA merge version to use
mrg_version = '121619'
# this is the most recent version with the full, updated variable list

#%% load modules
import pandas as pd
import numpy as np

#%% load data
# TOGA merge
mrg_fn = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/data analysis/sensitivity_tests/sensitivity_data/'\
+ 'TOGA_mrg_w_age_ptrfix_both_ns_mean_ndc_anth_tracers_'\
+mrg_version+CO_cutoff+'.'+HCN_cutoff+'.csv'
TOGA_mrg = pd.read_csv(mrg_fn)

# background estimates
bk_fn = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/data analysis/sensitivity_tests/sensitivity_data/'\
+'VOC_bks_ptrfix_both_ns_mean_ndc_anth_tracers_'\
+mrg_version+CO_cutoff+'.'+HCN_cutoff+'.csv'
VOC_bks = pd.read_csv(bk_fn)

# VOC names and chemical formulas
HAPs_RFs_fn = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/health risk tables/HAPs_EFs_RFs_forcode_final.csv'
HAPs_RFs_raw = pd.read_csv(HAPs_RFs_fn,header=0,skiprows=[1,39,40,41,42,43,44,45,46,47,48,
                                                          49,50,51])
HAPs_RFs = HAPs_RFs_raw.set_index('TOGA name')


#%% create dataframe for ratios and R2
VOC_PM_ratio_1 = pd.DataFrame(data = {'variable name':HAPs_RFs.index})
VOC_PM_ratio = VOC_PM_ratio_1.set_index('variable name')
VOC_PM_ratio['VOC name'] = HAPs_RFs['Name']
VOC_PM_ratio['Chemical Formula'] = HAPs_RFs['Chemical Formula']
    

#%% create VOC_elv dataframe by subtracting off background
VOC_smoke_elv = pd.DataFrame(data={'Start_datetimeUTC' : TOGA_mrg['datetime_start'],
                                   'Stop_datetimeUTC': TOGA_mrg['datetime_stop'],
                                   'RF_number':TOGA_mrg['RF_number']})
# add PM
VOC_smoke_elv['PM_elv_ams'] = TOGA_mrg['PM_TOGA_time_ams'] - VOC_bks['PM_ams'].values
VOC_smoke_elv['PM_elv_csd'] = TOGA_mrg['PM_TOGA_time_csd'] - VOC_bks['PM_csd'].values


for VOC_name in VOC_names:
    VOC_smoke_elv[VOC_name + '_elv'] = TOGA_mrg[VOC_name] - VOC_bks[VOC_name].values

# if elv < 0, make zero.
negs = np.where(VOC_smoke_elv <0.)
for i in range(len(negs[0])):
    VOC_smoke_elv.iloc[negs[0][i],negs[1][i]]=0
# add age
VOC_smoke_elv['chem_smoke_age'] = TOGA_mrg['chem_smoke_age']

#%% remove urban points
# anth_flag is a numerical flag indicating the number of urban tracers that exceed an
# estimated background. Four tracers are included.
# If any of them are elevated, remove points.
# This removes 61 of the total 3311 TOGA observations
anth_inds = np.where(TOGA_mrg['anth_flag']>0)[0]
TOGA_mrg_nourban = TOGA_mrg.drop(index=anth_inds)
VOC_smoke_elv_nourban = VOC_smoke_elv.drop(index=anth_inds)

#%% pull young smoke
smoke_elv_age_groups = VOC_smoke_elv_nourban.groupby(by='chem_smoke_age')
VOC_smoke_elv_young = smoke_elv_age_groups.get_group('young, < 7 hours')

#%% calculate VOC 2 PM linear relationship and R2
VOC_PM_ratio['VOC v PM pearson correlation squared'] = -999.
i = 0
for voc in VOC_PM_ratio.index:
    data = VOC_smoke_elv_young[[voc + '_elv','PM_elv_ams']]
    VOC_PM_ratio['VOC v PM pearson correlation squared'].iloc[i] = np.float64(data.corr('pearson').values[0,1]**2.)
    i += 1

#%% calculate average VOC 2 PM ratios
VOC_PM_weighted = VOC_smoke_elv_young[VOC_smoke_elv_young.columns[:5]].copy(deep=True)
zeroPM = np.where(VOC_PM_weighted['PM_elv_ams']==0)
VOC_PM_ratio['VOC PM1 ratio [ppt/ugm-3]'] = -999.
i = 0
for voc in VOC_PM_ratio.index:
    VOC_PM_weighted[voc+'_elv'] = VOC_smoke_elv_young[voc + '_elv'].values/VOC_smoke_elv_young['PM_elv_ams'].values
    VOC_PM_weighted[voc+'_elv'].iloc[zeroPM] = 0
    VOC_PM_ratio['VOC PM1 ratio [ppt/ugm-3]'].iloc[i] = np.float64(np.nanmean(VOC_PM_weighted[voc+'_elv']))
    i += 1
    print voc
#%% check this matches previous code

#%% save
VOC_PM_ratio.to_csv('/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/collabroation_w_Kat/VOC2PM_ratios.csv')
