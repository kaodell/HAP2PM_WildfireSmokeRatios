#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 11:32:31 2020
@author: kodell

check_age_distributions_wdilution.py
adapted from check_age_distributions

Python script to plot VOCs we use for age for each of the age categories
they allign to. This can test to see if this age binning is actually doing
what we want. 

written by: Katelyn O'Dell 
"""

#%% Load modules
import numpy as np
import pandas as pd
import plotly
import plotly.subplots
import plotly.graph_objs as go
import plotly.io as pio
pio.renderers.default = "chrome"
pio.templates.default = "seaborn"
# %% load data and set output location
file_path = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/data analysis/sensitivity_tests/sensitivity_data/'
TOGA_mrg = pd.read_csv(file_path + 'TOGA_mrg_w_age__NASA_R2_anth_tracers_PTR16fix_R1TOGAupdate_95pct_50pctall_85.275_w200_CH3CN.csv')

fig_path = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/paper_figures/final/SI_final_figures/'
#%% user-defined functions
# code to calc back
def calculate_bk(merge,var_name,bk_inds):
    var = np.copy(merge[var_name].values)
    var_bk = np.nanmedian(var[bk_inds])
    return var_bk 
def calculate_age_bk(merge,var_name,bk_inds):
    var = np.copy(merge[var_name].values)
    var_bk = np.round(np.nanpercentile(var[bk_inds],95),1)
    return var_bk 

#%% calc smoke-elevated values
# identify smoke and nosmoke inds
CO_max = 85.0
HCN_max = 275.0
CH3CN_max = 200.0
smoke_inds1 = np.where(np.logical_and(TOGA_mrg['CO_comb']>CO_max,TOGA_mrg[' HCN_TOGA']>HCN_max))[0]
not_smoke_inds1 = np.where(np.logical_and(TOGA_mrg['CO_comb']<=CO_max,TOGA_mrg[' HCN_TOGA']<=HCN_max))[0]
rmv_smk_inds = np.where(TOGA_mrg[' CH3CN_TOGA'].iloc[smoke_inds1] <= CH3CN_max)
rmv_nsmk_inds = np.where(TOGA_mrg[' CH3CN_TOGA'].iloc[not_smoke_inds1]>CH3CN_max)
smoke_inds = np.delete(smoke_inds1,rmv_smk_inds)
not_smoke_inds = np.delete(not_smoke_inds1,rmv_nsmk_inds)

# calc backgrounds
bk_xMeFuran = calculate_age_bk(TOGA_mrg,' x2MeFuran_TOGA',not_smoke_inds)
bk_Acrolein = calculate_age_bk(TOGA_mrg,' Acrolein_TOGA',not_smoke_inds)
bk_Acrylonitrile = calculate_age_bk(TOGA_mrg,' Acrylonitrile_TOGA',not_smoke_inds)
bk_CH3CN = calculate_bk(TOGA_mrg,' CH3CN_TOGA',not_smoke_inds)

# create dataframe of elevated values
VOC_smoke_elv = pd.DataFrame(data={
        'x2MeFuran_TOGA_elv':TOGA_mrg[' x2MeFuran_TOGA'].values - bk_xMeFuran,
        'Acrolein_TOGA_elv':TOGA_mrg[' Acrolein_TOGA'].values - bk_Acrolein,
        'Acrylonitrile_TOGA_elv':TOGA_mrg[' Acrylonitrile_TOGA'].values - bk_Acrylonitrile,
        'CH3CN_TOGA_elv':TOGA_mrg[' CH3CN_TOGA'].values - bk_CH3CN})

# replace negs with zeros
for VOC in ['CH3CN_TOGA_elv','x2MeFuran_TOGA_elv','Acrolein_TOGA_elv','Acrylonitrile_TOGA_elv']:
    VOC_smoke_elv[VOC]=np.where(VOC_smoke_elv[VOC]<0,0,VOC_smoke_elv[VOC])
    
VOC_smoke_elv['chem_smoke_age'] = TOGA_mrg['chem_smoke_age'].copy(deep=True)

#%% remove urban points
anth_inds = np.where(TOGA_mrg['anth_flag']>0)[0]
VOC_smoke_elv_nourban = VOC_smoke_elv.drop(index=anth_inds)
# drop does not reset index, so we need to reset this
VOC_smoke_elv_nourban.reset_index(inplace=True,drop=True)

#%% group by age
age_groups_VOC_elv = VOC_smoke_elv_nourban.groupby(by='chem_smoke_age')

young = age_groups_VOC_elv.get_group('young, < 0.5 days')
med = age_groups_VOC_elv.get_group('med, 0.5-1.5 days')
old = age_groups_VOC_elv.get_group('old, < 4 days')

#%% plot regular and dilutions-corrected distributions

fig = plotly.subplots.make_subplots(rows=3, cols=2,
                    subplot_titles=['&#916;2MethylFuran','&#916;2MethylFuran/&#916;Acetonitrile',
                                    '&#916;Acrolein','&#916;Acrolein/&#916;Acetonitrile',
                                    '&#916;Acrylonitrile','&#916;Acrylonitrile/&#916;Acetonitrile'])
# first plot non-dilution corrected
i=0
for VOC_plot in ['x2MeFuran_TOGA_elv','Acrolein_TOGA_elv','Acrylonitrile_TOGA_elv']:
    i+=1
    fig.add_trace(go.Box(x=old[VOC_plot],name='old',marker_color='orange',showlegend=False),row=i,col=1)
    fig.add_trace(go.Box(x=med[VOC_plot],name='medium',marker_color='red',showlegend=False),row=i,col=1)
    fig.add_trace(go.Box(x=young[VOC_plot],name='young',marker_color='purple',showlegend=False),row=i,col=1)
fig.update_xaxes(type='log',tickfont_size=20,title='mixing ratio [ppt]',row=3,col=1)
fig.update_yaxes(tickfont_size=20)
fig.update_layout(font_size=20)

# now plot dilution corrected values
i=0
showlegend_bool = [False,False,True]
for VOC_plot in ['x2MeFuran_TOGA_elv','Acrolein_TOGA_elv','Acrylonitrile_TOGA_elv']:
    i+=1
    fig.add_trace(go.Box(x=old[VOC_plot].values/old['CH3CN_TOGA_elv'].values,
                         name='old',marker_color='orange',showlegend=showlegend_bool[i-1]),row=i,col=2)
    fig.add_trace(go.Box(x=med[VOC_plot].values/med['CH3CN_TOGA_elv'].values,
                         name='medium',marker_color='red',showlegend=showlegend_bool[i-1]),row=i,col=2)
    fig.add_trace(go.Box(x=young[VOC_plot].values/young['CH3CN_TOGA_elv'].values,
                         name='young',marker_color='purple',showlegend=showlegend_bool[i-1]),row=i,col=2)
fig.update_xaxes(type='log',tickfont_size=20)
fig.update_yaxes(tickfont_size=20)
fig.update_layout(legend=dict(traceorder='reversed',font_size=20),plot_bgcolor='rgba(0,0,0,0)')
plotly.offline.plot(fig, filename= fig_path+'dilution_corrected_age_tacer_distributions.html')
fig.write_image(fig_path+'dilution_corrected_age_tacer_distributions.png',width=1200,height=900,scale=4)  