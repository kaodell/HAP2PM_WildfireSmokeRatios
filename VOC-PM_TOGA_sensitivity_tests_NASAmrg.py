#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
VOC-PM_TOGA_sensitivity_tests_NASAmrg.py
    python script to test VOC_PM sensitivity to:
        1) CO and HCN smoke cutoffs
        2) method of background estimation
        3) which PM data souce is used
Created on Tue Sep 17 09:04:38 2019
updated on 01.30.2020 to use NASA toga merge and renmaed with _NASAmrg
@author: kodell
"""
#%% Load modules
import numpy as np
import pandas as pd

#%% user inputs
# CO and HCN maxes to test
CO_maxes = [85]
HCN_maxes = [275]
CH3CN_max = 200 #ppt this is the value used in Singh et al., 2012
use_CH3CN = True
# whether or not to dilution correct
# don't do this; it doesn't work
dilution_correct = False

PTR_VOC_names = ['ACETAMIDE_C2H5NO','METHYL_METHACRYLATE_C5H8O2','PHENOL_C6H6O',
'QUINONE_C6H4O2','ISOCYANIC_ACID_HNCO', 'Formaldehyde_CH2O', 'METHANOL_C4HO']

# 5.1.2020 - add formaldyhyde and methanol to make sure PTR averaging matches 
# the TOGA merge values - update: it does!

# thsese should be in order [young, medium, old]
traceVOC_names = [' x2MeFuran_TOGA',' Acrolein_TOGA',' Acrylonitrile_TOGA']
age_names = ['not smoke', 'young, < 0.5 days','med, 0.5-1.5 days', 'old, < 4 days']
age_colors = ['blue','purple','red','orange','green']

# paths for figures and datafiles
fig_path = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/data analysis/sensitivity_tests/sensitivity_figures/'
outfn_path = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/data analysis/sensitivity_tests/sensitivity_data/'

#TOGA mrg to load and PTR filepath
mrg_fn = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/data analysis/processed datafiles/TOGA_mrg_NASA_R2_anth_tracers_R1TOGAupdate.csv'
PTR_fp = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/WECAN_datafiles/PTR_fromWade/Kate_Odell_HAPS/'

# this indictes method of background calculation and PM estimation portion indicating HCN and CO elevation done in loop
# nameing convetion in sensitivity_readme file
desc_name_pt1 = '_NASA_R2_anth_tracers_PTR16fix_R1TOGAupdate_95pct_50pctall_' 
#%% user defined functions
# here can change what is used for age-tracer and all VOCs background

# function for age tracers, use mean no smoke + 1 standard deviation
# update: post-review we use 95th percentile
def find_age_inds(mrg_data,tracer_name,nsmoke_inds,smoke_inds,dc):
    #age_bk = np.nanmean(mrg_data[tracer_name].iloc[nsmoke_inds]) + np.nanstd(mrg_data[tracer_name].iloc[nsmoke_inds],ddof=1)
    age_bk = np.round(np.nanpercentile(mrg_data[tracer_name].iloc[nsmoke_inds],95),1)
    if dc:
        age_ER = (mrg_data[tracer_name] - age_bk)/(mrg_data['CO_comb']-bk_CO)
        smoke_age_inds = np.where(age_ER[smoke_inds] > 0 )
        age_inds = smoke_inds[smoke_age_inds[0]]
    else:
        smoke_age_inds = np.where(mrg_data[tracer_name].iloc[smoke_inds] > age_bk)
        age_inds = smoke_inds[smoke_age_inds]

    return age_inds, age_bk

# function for other VOCs, use median no smoke
def calculate_bk(merge,var_name,bk_inds):
    var = np.copy(merge[var_name].values)
    #var_std = np.nanstd(var[bk_inds])
    var_bk = np.nanmedian(var[bk_inds])
    return var_bk 
#%% import data
# load data from 'find_anthropogenic_tracers_NASAmrg.py' output
TOGA_mrg = pd.read_csv(mrg_fn)
# load TOGA limit of detection
# TOGA limits of detection from TOGA R1 data
TOGA_LLOD_fn =  '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/WECAN_datafiles/TOGA_data_05.21.20/TOGA_LLOD_R1.csv'
TOGA_LLOD = pd.read_csv(TOGA_LLOD_fn)

# load HAPs spreadsheet for healthVOC names
HAPs_RFs_fn = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/health risk tables/HAPs_EFs_RFs_forcode_final_final_WadeUpdate_onlyREL4acute.csv'
HAPs_RFs_raw = pd.read_csv(HAPs_RFs_fn,header=0,skiprows=[1,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,
                                                          49,50,51])
HAPs_RFs = HAPs_RFs_raw.set_index('TOGA name')
healthVOCs = HAPs_RFs.index

#%% for species in TOGA merge, replace -9999999 with nans and -8888888 with 0.5*LOD
# new TOGA vars use -999,-888,-777; other WE-CAN vars use -9999999
TOGA_mrg.replace(to_replace = -999.0, value = np.nan, inplace = True)
TOGA_mrg.replace(to_replace = -9999999.0, value = np.nan, inplace = True)

TOGA_LLOD_flag = -888.0

# replace LLOD with limit of detection for all tracers
# test subbing 0 or the detection limit here to make sure it isn't sensitive
for VOC_name in TOGA_LLOD.keys()[3:]:
    TOGA_mrg[' ' + VOC_name].replace(to_replace = TOGA_LLOD_flag,value = 0.5*TOGA_LLOD[VOC_name][0],inplace=True)
    # sensitivity check for paper
    # TOGA_mrg[' ' + VOC_name].replace(to_replace = TOGA_LLOD_flag,value = TOGA_LLOD[VOC_name][0],inplace=True)

#PTR LOD replacement is handled later
    
#%% when co from the qcl is missing, use picarro, otherwise use qcl 
# call new combined CO CO_comb
TOGA_mrg['CO_comb'] = np.where(np.isnan(TOGA_mrg[' CO_QCL']),
                          TOGA_mrg[' CO_PICARRO'].values,
                          TOGA_mrg[' CO_QCL'].values)

#%% create a column that adds black carbon mass and total AMS mass (these are both in ug/m3 STP)
TOGA_mrg['PM1_ams_sp2'] = TOGA_mrg[' M_rBC_SDI_SP2'] + TOGA_mrg[' Total_AMS']

#%% load all PTR data and prep to average for TOGA flight
# define dates and RFs for each flight
dates = ['0724','0726','0730','0731','0802','0803','0806','0808','0809','0813','0815','0816','0820','0823','0826','0828','0906','0910','0913']
RFs = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19'] # 17, 18, and 19 are educational flights; no TOGA data
# load first flight to make array
PTR_filename = PTR_fp + '/WECAN-PTRTOFMS_C130_2018'+ dates[0] + '_R1.ict'
PTR_mrg_all_raw = pd.read_csv(PTR_filename,header=77)
PTR_mrg_all_raw['datetime_start'] = [int(RFs[0])*10**5.]*PTR_mrg_all_raw.shape[0] + PTR_mrg_all_raw['Start_UTC']
# in header Wade says LLODs are not flagged
# load all files
load_date_strs = dates[1:]
i = 1
for date_str in load_date_strs:
    if date_str != '0826': # no PTR data on this day
        PTR_filename = PTR_fp + '/WECAN-PTRTOFMS_C130_2018'+ date_str + '_R1.ict'
        PTR_mrg_raw = pd.read_csv(PTR_filename,header=77)
        PTR_mrg_raw['datetime_start'] = [int(RFs[i])*10**5.]*PTR_mrg_raw.shape[0] + PTR_mrg_raw['Start_UTC']
        PTR_mrg_all_raw = pd.concat([PTR_mrg_all_raw,PTR_mrg_raw])
        del PTR_mrg_raw
    i+=1
PTR_mrg_all_raw.reset_index(inplace=True)

# make -99999 nans
PTR_mrg_all = PTR_mrg_all_raw.replace(to_replace = -99999., value = np.nan)

# replace below detection limit values with  1/2 the limit of detection
# use LOD of 200 ppt (0.2 ppb) per Wade email to me on 11.05.19
for PTR_VOC_name in PTR_VOC_names:
    PTR_mrg_all[PTR_VOC_name] = np.where(PTR_mrg_all[PTR_VOC_name].values < 0.2, 
                   0.1, # ppb, 1/2 the detection limit of 0.2 ppb
                   PTR_mrg_all[PTR_VOC_name].values)
    # also tested sensitivity of this detection limit on results

#%% average PTR data to TOGA timescale
# add datetime start and stop to TOGA so these can be matched with PTR
TOGA_mrg['datetime_start'] = TOGA_mrg[' FLIGHT'].values*10**5. + TOGA_mrg['  STARTTIME']
TOGA_mrg['datetime_stop'] = TOGA_mrg[' FLIGHT'].values*10**5. + TOGA_mrg[' STOPTIME']

# now merge to TOGA timescale
ptr_start_times = np.array(PTR_mrg_all['datetime_start'])
ptr_stop_times = np.array(PTR_mrg_all['datetime_start'] + 1.0)
for var_name in PTR_VOC_names:
    TOGA_mrg[' '+var_name+'_PTR'] = [-555]*TOGA_mrg.shape[0]
for i in range(TOGA_mrg.shape[0]):
    toga_start_time = TOGA_mrg['datetime_start'].iloc[i]
    toga_stop_time = TOGA_mrg['datetime_stop'].iloc[i]

    ptr_ind_start = np.where(ptr_start_times == toga_start_time)[0]
    ptr_ind_stop = np.where(ptr_stop_times == toga_stop_time)[0] 

    if len(ptr_ind_start) < 1: 
        # no PTR data at start of TOGA, sometimes happens at start/end of flights
        for PTR_VOC_name in PTR_VOC_names:
             TOGA_mrg[' '+PTR_VOC_name+'_PTR'].iloc[i] = np.nan
        print 'ptr times dont allign @ start'
    elif len(ptr_ind_stop) < 1:
        # no PTR data at end of TOGA, sometimes happens at start/end of flights
        for PTR_VOC_name in PTR_VOC_names:
             TOGA_mrg[' '+PTR_VOC_name+'_PTR'].iloc[i] = np.nan
        print 'ptr times dont allign @ stop'
    else:   
        # assign value and convert to ppt
        for PTR_VOC_name in PTR_VOC_names:
            TOGA_mrg[' '+PTR_VOC_name+'_PTR'].iloc[i] = 1000.0*np.nanmean(PTR_mrg_all[PTR_VOC_name].iloc[ptr_ind_start[0]:ptr_ind_stop[0]+1],axis=0)          
    if i%100 == 0:
        print 100.*(float(i)/TOGA_mrg.shape[0]), '% done' 


#%% identify smoke/ nosmoke inds and calculate backgrounds for different sensitivitys
#   save data
# outter loops are CO/HCN maxes, for final paper version, these are just one value, defined above
for CO_max in CO_maxes:
    for HCN_max in HCN_maxes:
        # create description name for output file
        desc_name = desc_name_pt1 + str(CO_max) +'.' + str(HCN_max)
        
        # identify no smoke and smoke inds
        if use_CH3CN:
            smoke_inds1 = np.where(np.logical_and(TOGA_mrg['CO_comb']>CO_max,TOGA_mrg[' HCN_TOGA']>HCN_max))[0]
            not_smoke_inds1 = np.where(np.logical_and(TOGA_mrg['CO_comb']<=CO_max,TOGA_mrg[' HCN_TOGA']<=HCN_max))[0]
            rmv_smk_inds = np.where(TOGA_mrg[' CH3CN_TOGA'].iloc[smoke_inds1] <= CH3CN_max)[0]
            rmv_nsmk_inds = np.where(TOGA_mrg[' CH3CN_TOGA'].iloc[not_smoke_inds1]>CH3CN_max)
            smoke_inds = np.delete(smoke_inds1,rmv_smk_inds)
            not_smoke_inds = np.delete(not_smoke_inds1,rmv_nsmk_inds)
        print len(not_smoke_inds)    
        '''
        else:
            smoke_inds = np.where(np.logical_and(TOGA_mrg['CO_comb']>CO_max,TOGA_mrg[' HCN_TOGA']>HCN_max))[0]
            not_smoke_inds2 = np.where(np.logical_and(TOGA_mrg['CO_comb']>CO_max,TOGA_mrg[' HCN_TOGA']>HCN_max)==False)[0]
            # if any of the tracers are missing, remove these points
            CO_nans = np.where(np.isnan(TOGA_mrg['CO_comb'].iloc[not_smoke_inds2])==True)
            HCN_nans = np.where(np.isnan(TOGA_mrg[' HCN_TOGA'].iloc[not_smoke_inds2])==True)
            CH3CN_nans = np.where(np.isnan(TOGA_mrg[' CH3CN_TOGA'].iloc[not_smoke_inds2])==True)
            all_nans = np.hstack([CO_nans,HCN_nans,CH3CN_nans])
            unans = np.unique(all_nans)
            # 240 of these, remove them from not_smoke_inds
            not_smoke_inds = np.delete(not_smoke_inds2,all_nans)
        '''

        # inner section calculates different backgrounds and saves data
        bk_CO = calculate_bk(TOGA_mrg,'CO_comb',not_smoke_inds)
        PMbk = calculate_bk(TOGA_mrg,'PM1_ams_sp2',not_smoke_inds)
        
        #create array to save VOC bks
        VOC_bks = pd.DataFrame(data={' STARTTIME' : TOGA_mrg['  STARTTIME'],
                                           ' STOPTIME': TOGA_mrg[' STOPTIME'],
                                           ' FLIGHT':TOGA_mrg[' FLIGHT']})
        # add PM and CO
        VOC_bks['PM1_ams_sp2'] = PMbk
        VOC_bks['CO_comb'] = bk_CO
        
        for VOC_name in healthVOCs:
                VOC_bks[VOC_name] = calculate_bk(TOGA_mrg,' '+VOC_name,not_smoke_inds)
        
        ### now for age estimates ###
        
        # young
        young_inds, young_bk = find_age_inds(TOGA_mrg,traceVOC_names[0],not_smoke_inds,smoke_inds,dilution_correct)
        # med
        med_inds, med_bk = find_age_inds(TOGA_mrg,traceVOC_names[1],not_smoke_inds,smoke_inds,dilution_correct)
        # old
        old_inds, old_bk = find_age_inds(TOGA_mrg,traceVOC_names[2],not_smoke_inds,smoke_inds,dilution_correct)

        # now assign age to an array
        # note, order here matters! young smoke will also have an elevated 'old' tracer
        # need to assign young age last, or else the old inds would overwrite it
        chem_age = np.empty(TOGA_mrg.shape[0],dtype='|S32')
        chem_age[:] = 'NA'
        chem_age[not_smoke_inds] = 'not smoke'
        chem_age[smoke_inds] = '> 4 days'
        chem_age[old_inds] = 'old, < 4 days'
        chem_age[med_inds] = 'med, 0.5-1.5 days'
        chem_age[young_inds] = 'young, < 0.5 days'

        # remove points where any of the age tracers or smoke identifiers are nans
        nans1 = np.where(np.isnan(TOGA_mrg[traceVOC_names[0]]))
        nans2 = np.where(np.isnan(TOGA_mrg[traceVOC_names[1]]))
        nans3 = np.where(np.isnan(TOGA_mrg[traceVOC_names[2]]))
        allCO_nans = np.where(np.isnan(TOGA_mrg['CO_comb']))
        allHCN_nans = np.where(np.isnan(TOGA_mrg[' HCN_TOGA']))
        allCH3CN_nans = np.where(np.isnan(TOGA_mrg[' CH3CN_TOGA']))

        all_nans = np.hstack([nans1,nans1,nans3,allCO_nans,allHCN_nans,allCH3CN_nans])
        unans = np.unique(all_nans)
        chem_age[unans] = 'NA'

        # add age to TOGA dataframe and define colors for the ages
        TOGA_mrg['chem_smoke_age'] = chem_age

        # save out to csv
        TOGA_mrg.to_csv(path_or_buf= outfn_path + 'TOGA_mrg_w_age_'+desc_name+'_w200_CH3CN.csv')
        VOC_bks.to_csv(path_or_buf= outfn_path + 'VOC_bks_'+desc_name+'_w200_CH3CN.csv')
        fn = outfn_path + 'bk_inds' + desc_name + '_w200_CH3CN' +'.npz'
        np.savez(fn,not_smoke_inds=not_smoke_inds)


