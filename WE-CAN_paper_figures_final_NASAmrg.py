#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
WE-CAN_paper_figures_final_NASAmrg.py
    script to create the final versions of WE-CAN paper figures, 
    adapted from WE-CAN_paper_figures.py
Created on Mon Nov 11 12:07:52 2019
updated 01.31.2020 to use NASA TOGA merge
@author: kodell
"""
#%% import modules
import numpy as np
import pandas as pd
import plotly
import plotly.graph_objs as go
import plotly.io as pio
pio.renderers.default = "chrome"
pio.templates.default = "seaborn"
mapbox_access_token='pk.eyJ1Ijoia2FvZGVsbCIsImEiOiJjanZza3k1bGkzMHZoNDhwYjdybTYyOTliIn0.KJyzHWVzu2U087Ps215_LA'
#%% user inputs

# for paper
fig_path = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/paper_figures/final/'
fig_desc = '_NASA_R2mrg_85.275_w200_CH3CN_PTRfix_R1TOGAupdate'
mrg_name = 'TOGA_mrg_w_age__NASA_R2_anth_tracers_PTR16fix_R1TOGAupdate85.275_w200_CH3CN.csv'
voc_bk_name = 'VOC_bks__NASA_R2_anth_tracers_PTR16fix_R1TOGAupdate85.275_w200_CH3CN.csv'

# sensitivity tests
#fig_path = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/data analysis/sensitivity_tests/sensitivity_figures/leave in negative enhancements/'
#fig_desc = '_NASA_R2mrg_85.275_w200_CH3CN_leave_enhance_negs'

#fig_path = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/data analysis/sensitivity_tests/sensitivity_figures/background_sensitivity/'
#fig_desc = '_NASA_R2mrg_85.275_w200_CH3CN_allbk_plusonesig'
#mrg_name = 'TOGA_mrg_w_age__NASA_R2_anth_tracers_bk_nagebk_plus1sig85.275_w200_CH3CN.csv'
#voc_bk_name = 'VOC_bks__NASA_R2_anth_tracers_bk_nagebk_plus1sig85.275_w200_CH3CN.csv'

#fig_path = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/data analysis/sensitivity_tests/sensitivity_figures/detection_limit_sensitivity/'
#fig_desc = '_NASA_R2mrg_85.275_w200_CH3CN_addPTRcheck_04bdl'
#mrg_name = 'TOGA_mrg_w_age__NASA_R2_anth_tracers_bk_addPTRcheck_04bdl85.275_w200_CH3CN.csv'
#voc_bk_name = 'VOC_bks__NASA_R2_anth_tracers_bk_addPTRcheck_04bdl85.275_w200_CH3CN.csv'


age_colors = ['grey','purple','red','orange']
age_names = ['not smoke', 'young, < 0.5 days','med, 0.5-1.5 days', 'old, < 4 days']
age_names_fig = np.array(['not smoke','young, < 1 day', 'medium, 1-3 days','old, > 3 days' ])
states_plot = ['Colorado','California','Idaho','Montana','Wyoming','Washington',
               'Oregon','Nevada','Utah']

#%% parameters and constants
P_std = 101325 #Pa
T_std = 273.15 #K
R = 8.314 #j/molK

#%% load data
mrg_fn = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/data analysis/sensitivity_tests/sensitivity_data/' + mrg_name
TOGA_mrg = pd.read_csv(mrg_fn)
TOGA_mrg.drop(columns=['Unnamed: 0', 'Unnamed: 0.1'], inplace=True)
# load calculated VOC background concentrations
bk_fn = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/data analysis/sensitivity_tests/sensitivity_data/' + voc_bk_name
VOC_bks = pd.read_csv(bk_fn)

# spreadsheet with HAPs, reference concentrations, and molecular weight
HAPs_RFs_fn = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/health risk tables/HAPs_EFs_RFs_forcode_final_final_WadeUpdate.csv'
HAPs_RFs_raw = pd.read_csv(HAPs_RFs_fn,header=0,skiprows=[1,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,
                                                          49,50,51])

#%% reset HAPs index and pull healthVOC names
HAPs_RFs = HAPs_RFs_raw.set_index('TOGA name')
healthVOCs = HAPs_RFs.index
healthVOCs_fignames = HAPs_RFs['Name'].values

#%% Convert units of chronic cancer risk ratios and noncancer hazard values to ppt
# Convert units of chronic cancer risk ratios from 1/(ug/m3) -> 1/ppt
# conversion factor converts convert 1/ug/m3 to 1/ppt
conv_factor = ((P_std)/(T_std*R))*(1.*10.**-12.)*(1.*10.**6)*(1.*10.**6) #second 10**6 so number is ppl/mil risk

crf_ppt = HAPs_RFs['chronic cancer risk factor']*HAPs_RFs['molecular weight']*conv_factor   
HAPs_RFs['crf_ppt'] = crf_ppt

# convert referce concentrations from mg/m3 to ppt
arl_ppt = HAPs_RFs['acute risk factor']*(1./(1000.*HAPs_RFs['molecular weight']))*(R*T_std/P_std)*(1.*10.**12)
HAPs_RFs['arl_ppt'] = arl_ppt

crl_ppt = HAPs_RFs['chronic noncancer risk factor']*(1./(1000.*HAPs_RFs['molecular weight']))*(R*T_std/P_std)*(1.*10.**12)
HAPs_RFs['crl_ppt'] = crl_ppt

#%% combine old and v old age categories
TOGA_mrg.loc[TOGA_mrg.chem_smoke_age =='> 4 days', 'chem_smoke_age'] = 'old, < 4 days'

#%% remove urban points
anth_inds = np.array(TOGA_mrg['anth_flag']==0)
TOGA_mrg_nourban = TOGA_mrg.iloc[anth_inds]
TOGA_mrg_nourban.reset_index(inplace=True,drop=True)
#%% remove TOGA inds where AMS is nans or all TOGA or all PTR are nans
all_nans = []
count = 0
for i in range(TOGA_mrg_nourban.shape[0]):
    if np.isnan(TOGA_mrg_nourban['PM1_ams_sp2'].iloc[i]):
        all_nans.append(i)
    elif np.all(np.isnan(TOGA_mrg_nourban[' '+ healthVOCs].iloc[i])):
        all_nans.append(i) # this happens ~ 30 times
        print i
        count += 1
unans = np.unique(all_nans)
TOGA_mrg_datause = TOGA_mrg_nourban.drop(unans,axis=0)
TOGA_mrg_datause.reset_index(inplace=True,drop=True)
#%% make df for fig names and colors
fig_names = pd.DataFrame(data = {'VOC_names' : ' '+healthVOCs,'fig_names' : healthVOCs_fignames})
fig_names = fig_names.set_index('VOC_names')
# add colors for plotting in box plots
chemical_colors =  fig_names.shape[0]*['grey']
fig_names['color'] = chemical_colors

i=0
VOCs_w_color = ['CH3CHO_TOGA','Acrolein_TOGA','Acrylonitrile_TOGA','Benzene_TOGA','CH2O_TOGA','HCN_TOGA']
colors = ['#ec7063','#5dade2','#f39c12','#45b39d','#1F618D','#884EA0 ']
for VOC_name in VOCs_w_color:
    fig_names['color'][' '+VOC_name] = colors[i]
    i += 1

#%% create VOC_elv dataframe by subtracting off background
VOC_smoke_elv = pd.DataFrame(data={'  STARTTIME' : TOGA_mrg['  STARTTIME'],
                                   ' STOPTIME': TOGA_mrg[' STOPTIME'],
                                   ' FLIGHT':TOGA_mrg[' FLIGHT']})
# add PM and CO
VOC_smoke_elv['PM_elv'] = TOGA_mrg['PM1_ams_sp2'] - VOC_bks['PM1_ams_sp2'].values
VOC_smoke_elv['CO_comb_elv'] = TOGA_mrg['CO_comb'] - VOC_bks['CO_comb'].values

for VOC_name in healthVOCs:
    VOC_smoke_elv[' '+VOC_name + '_elv'] = TOGA_mrg[' '+VOC_name] - VOC_bks[VOC_name].values
# if elv < 0, make zero.
# test not doing this to see if main conclusions will change... they don't. 
negs = np.where(VOC_smoke_elv <0.)
for i in range(len(negs[0])):
    VOC_smoke_elv.iloc[negs[0][i],negs[1][i]]=0
# age chem age
VOC_smoke_elv['chem_smoke_age'] = TOGA_mrg['chem_smoke_age'].copy(deep=True)
# remove urban points
anth_inds = np.where(TOGA_mrg['anth_flag']>0)[0]
VOC_smoke_elv_nourban = VOC_smoke_elv.drop(index=anth_inds)
VOC_smoke_elv_nourban.reset_index(inplace=True,drop=True)
####################################################
# Figures
####################################################  
#%% Figure 1: Flight tracks with data used in this study colored by chem age 
fig = go.Figure()
# make trace of all flight tracks
fig.add_trace(go.Scattergeo(
    lat=TOGA_mrg[' LATITUDE'],
    lon=TOGA_mrg[' LONGITUDE'],
    mode='lines',
    marker_color='grey',marker_size=1, name = 'WE-CAN flight track'))
i=0
age_groups_use = TOGA_mrg_datause.groupby(by='chem_smoke_age')
n = []
for age_name in age_names[1:]:
    group = age_groups_use.get_group(age_name)
    nobs = len(group)
    n.append(nobs)
    fig.add_trace(go.Scattergeo(lat=group[' LATITUDE'],lon=group[' LONGITUDE'],mode='markers',
        marker = dict(color=age_colors[i+1],size=4,opacity=0.7),
        hovertext = group[' FLIGHT'],name = age_names_fig[i+1] + ' (' + str(nobs)+')'))
    i += 1
fig.update_layout(
    geo = dict(scope = 'north america',showland = True,landcolor = "white",
        subunitcolor = "black",countrycolor = "black",showlakes = True,lakecolor = "rgb(255, 255, 255)",
        showsubunits = True,showcountries = True,resolution = 50,
        projection = dict(type = 'mercator',rotation_lon = -100),
        lonaxis = dict(showgrid = False,gridwidth = 0.5,range= [ -126.0, -108.0 ],dtick = 5),
        lataxis = dict(showgrid = False,gridwidth = 0.5,range= [ 37.0, 49.0 ],dtick = 5)),
    legend=dict(x=0.43,y=-0.08,traceorder="normal",
        font=dict(size=8.5,color="black"),
        bgcolor="White",bordercolor="Grey",borderwidth=1))
fig.show()
fig.write_image(fig_path+fig_desc+  'ChemSmokeAge_wout_vold.png',scale=8)  

#%% Figure 2: distribution of HAPs for each age with reference concentration
# first, sort VOC list by WE-CAN ER
sorted_healthVOCs = HAPs_RFs['WE-CAN ER'][healthVOCs].sort_values(ascending=True,
                            na_position='last').keys()
# also create a list with the space infront of the name to load from the mrg
HAPs_mrg_names = np.copy(sorted_healthVOCs)
for i in range(len(HAPs_mrg_names)):
    HAPs_mrg_names[i]=' '+HAPs_mrg_names[i]

fig = go.Figure()
fig.add_trace(go.Scatter(x=HAPs_RFs['arl_ppt'][sorted_healthVOCs],
                         y=fig_names['fig_names'][HAPs_mrg_names],
         mode = 'markers',marker_color='black',marker_size=5,
         marker=dict(symbol=2),name='acute reference concentration'))
i=0
for age_name in np.flip(age_names[1:],axis=0):
    names = []
    concs = []
    inds = np.where(TOGA_mrg_nourban['chem_smoke_age'] == age_name)
    for VOC_name in sorted_healthVOCs:
        conc = VOC_smoke_elv_nourban[' '+VOC_name + '_elv'].iloc[inds]
        name = [fig_names['fig_names'][' '+VOC_name]]*len(inds[0])
        concs = np.hstack([concs,conc])
        names = np.hstack([names,name])
    # remove 0's
    lowinds = np.where(concs == 0)
    concs = np.delete(concs,lowinds)
    names = np.delete(names,lowinds)
    names = np.array(names)
    concs = np.array(concs)
    fig.add_trace(go.Box(x = concs,y=names,name=np.flip(age_names_fig[1:],axis=0)[i],
                         boxpoints=False,marker_color=np.flip(age_colors[1:],axis=0)[i]))
    i += 1  
fig.update_layout(boxmode='group',
                  plot_bgcolor='white',
                  legend=dict(traceorder='reversed',x=0.21,y=0.03,font_size=18,
                              bordercolor="White",borderwidth=1))
fig.update_xaxes(title_text='background-corrected mixing ratio [pptv]',
                 type='log',range=[-1,8],tickfont=dict(size=16),titlefont=dict(size=20))             
fig.update_traces(orientation='h')
fig.update_yaxes(tickfont=dict(size=14))
plotly.offline.plot(fig, filename= fig_path+ fig_desc+ '_VOC_elv_boxplot_hazardv2')
fig.write_image(fig_path+fig_desc+ '_VOC_elv_boxplot_hazardv2.png',width=900,height=1200,scale=4)  
 

#%% Figure 3: bar charts of total health hazard for each age group, weighted by PM 
#           also supplemental figure of each weighted by CO for supplement

# first calc PM- and CO-weighted VOCs for each TOGA observation
VOC_PM_weighted = VOC_smoke_elv_nourban[VOC_smoke_elv_nourban.columns[:5]].copy(deep=True)
VOC_CO_weighted = VOC_smoke_elv_nourban[VOC_smoke_elv_nourban.columns[:5]].copy(deep=True)
zeroPM = np.where(VOC_smoke_elv_nourban['PM_elv']==0)
PM_nans = np.where(np.isnan(VOC_smoke_elv_nourban['PM_elv'])) # we will use this in the loop

for voc in VOC_smoke_elv_nourban.columns[5:-1]:
    VOC_PM_weighted[voc] = VOC_smoke_elv_nourban[voc].values/VOC_smoke_elv_nourban['PM_elv'].values
    # divide by zero gives inf, replace with nans
    VOC_PM_weighted[voc].iloc[zeroPM] = np.nan

    VOC_CO_weighted[voc] = VOC_smoke_elv_nourban[voc].values/VOC_smoke_elv_nourban['CO_comb_elv'].values
    # make CO nans here to so we are comparing the same data
    # CO_elv is never 0 in smoke, because the background CO calculation (~70 ppb) is less 
    # than the minimum required for smoke events (85 ppb), so it doesn't make a difference to replace CO infs
    VOC_CO_weighted[voc].iloc[zeroPM] = np.nan

VOC_PM_weighted['chem_smoke_age'] = VOC_smoke_elv_nourban['chem_smoke_age'].copy(deep=True)
VOC_CO_weighted['chem_smoke_age'] = VOC_smoke_elv_nourban['chem_smoke_age'].copy(deep=True)

# make figure and calc medians of PM- and CO- weighted VOCs for each age category
legend_bool = [False,False,True] # to only plot legend on final age group so it isn't repeated
f = 0 #indicate which figure is being plotted
for VOC_weighted in [VOC_PM_weighted,VOC_CO_weighted]:
    row = 1
    
    # create figure 
    fig = plotly.subplots.make_subplots(rows=3, cols=1,
    subplot_titles = ['(a) Acute','(b) Chronic Noncancer','(c) Cancer'])

    for risk in ['arl_ppt','crl_ppt','crf_ppt']:
        VOCs_w_risk = []
        VOC_risk_weighted = VOC_weighted[VOC_weighted.columns[:]].copy(deep=True)
        # now calculate risk per 10 ug PM
        if risk in ['arl_ppt','crl_ppt']:
            for voc in VOC_smoke_elv_nourban.columns[5:-1]:
                VOC_risk_weighted[voc] = 10.0*VOC_weighted[voc]/HAPs_RFs[risk][voc[1:-4]]
                if np.isfinite(HAPs_RFs[risk][voc[1:-4]]):
                    VOCs_w_risk.append(voc)
        if risk == 'crf_ppt':
            for voc in VOC_smoke_elv_nourban.columns[5:-1]:
                VOC_risk_weighted[voc] = 10.0*VOC_weighted[voc]*HAPs_RFs[risk][voc[1:-4]]
                if np.isfinite(HAPs_RFs[risk][voc[1:-4]]):
                    VOCs_w_risk.append(voc)

        age_groups = VOC_risk_weighted.groupby(by='chem_smoke_age')
        
        total_risk_methodA = [0,0,0] # need total risk via this method to add bars
        # add 'other' first so it's at the bottom of the colors
        other_med_risk = [0,0,0]
        for VOC_name in healthVOCs:
            if VOC_name in VOCs_w_color:
                continue
            other_med_risk = np.vstack([other_med_risk,
                                              age_groups.median()[' '+VOC_name + '_elv'][age_names[1:]].values])
        other_med_risk_tot = np.nansum(other_med_risk,axis=0)
        total_risk_methodA+= other_med_risk_tot

        fig.add_trace(go.Bar(x=age_names_fig[1:],y=other_med_risk_tot,name='other',
                             legendgroup='other',showlegend=legend_bool[row-1],
                             marker_color='#7F8C8D'),row=row,col=1)
        # add VOCs with color
        for VOC_name in VOCs_w_color:
                # calculate median for VOCs with color
                med_risk = age_groups.median()[' '+VOC_name + '_elv'][age_names[1:]].values
                # add to figure
                fig.add_trace(go.Bar(x=age_names_fig[1:], y=med_risk,
                         name=fig_names['fig_names'][' '+VOC_name],
                         legendgroup=fig_names['fig_names'][' '+VOC_name],
                         showlegend=legend_bool[row-1],
                         marker_color=fig_names['color'][' '+VOC_name]),row=row,col=1)
                total_risk_methodA += np.where(np.isnan(med_risk),0,med_risk)
            
        # add 'error' bars
        med_tot_risk = []
        lower_q_tot_risk = []
        upper_q_tot_risk = []

        # sum for total risk of individual points
        total_risk = np.nansum(VOC_risk_weighted.iloc[:,5:-1].values,axis=1)
        # make places where PM_elv is nans or zeros nans, in the sum they turn into zeros
        # PM_elv is only zero 3 times during the smoke points (twice in young, once in old)
        total_risk[PM_nans] = np.nan
        total_risk[zeroPM] = np.nan
        # also remove places where any of the VOCs with risk are nans, these will
        # be counted as zeros in the sum, but shouldn't be and will bias risk estiamte low
        rmv_inds = []
        for i in range(VOC_risk_weighted.shape[0]):
            if np.any(np.isnan(VOC_risk_weighted[VOCs_w_risk].iloc[i])):
                rmv_inds = np.append(rmv_inds,i)
        rmv_inds = np.array(rmv_inds,dtype=int)
        total_risk[rmv_inds] = np.nan

        # calculate 25th and 75th percentile for 'error' bars
        for age_name in age_names[1:]:
            age_inds = np.where(VOC_risk_weighted['chem_smoke_age']==age_name)
            total_risk_age = total_risk[age_inds]
            # print 'total risk length', age_name, len(np.where(np.isfinite(total_risk_age))[0])
            med_tot_risk = np.append(med_tot_risk,np.nanmedian(total_risk_age))
            lower_q_tot_risk = np.append(lower_q_tot_risk,np.nanpercentile(total_risk_age,25,interpolation='linear'))
            upper_q_tot_risk = np.append(upper_q_tot_risk,np.nanpercentile(total_risk_age,75,interpolation='linear'))
        print 'methodA, sum indv. med\n',total_risk_methodA, 'methodB, total med\n', med_tot_risk
        # plot the risk based on total_med risk in method A, this way
        # what is plotted is the 25th-75th percentile of the total risk, when it
        # can be estimated for individual points (all species w/ risk have data)
        fig.add_trace(go.Bar(x=age_names_fig[1:], y=[0,0,0],
             name='total_med',showlegend=False,
             error_y=dict(type='data', symmetric=False,
                           array=upper_q_tot_risk-total_risk_methodA,
                           arrayminus=total_risk_methodA-lower_q_tot_risk)),row=row,col=1)
        row += 1                     
    # clean up figure and add labels
    fig.update_xaxes(tickfont=dict(size=16),row=1,col=1)
    fig.update_xaxes(tickfont=dict(size=16),row=2,col=1)
    fig.update_xaxes(tickfont=dict(size=16),row=3,col=1)
    fig.update_layout(barmode='stack',plot_bgcolor='white',
                      xaxis=dict(zeroline=True),
                      legend=dict(traceorder='reversed',
                                  x=0.87,y=1.06),
                      margin=dict(l=80, r=80, t=20, b=20))
    if f==0:
        fig.update_yaxes(title="hazard index <br> 10 &mu;g m<sup>-3</sup> PM<sub>1</sub>",
                         color='black',
                         row=1, col=1)
        fig.update_yaxes(title="hazard index <br> 10 &mu;g m<sup>-3</sup> PM<sub>1</sub>",
                         color='black',
                         row=2, col=1)
        fig.update_yaxes(title="cancer risk per 10<sup>6</sup> <br> 10 &mu;g m<sup>-3</sup> PM<sub>1</sub>",
                         color='black',
                         row=3, col=1)
        plotly.offline.plot(fig, filename= fig_path+ fig_desc + 'barchart_PMweight.html')
        fig.write_image(fig_path+fig_desc+ '_barchartrisks_PMweight.png',width=900,height=700,scale=4) 
    if f==1:
        fig.update_yaxes(title="hazard index <br> 10 ppb CO",
                 color='black',
                 row=1, col=1)
        fig.update_yaxes(title="hazard index <br> 10 ppb CO",
                         color='black',
                         row=2, col=1)
        fig.update_yaxes(title="cancer risk per 10<sup>6</sup> <br> 10 ppb CO",
                         color='black',
                         row=3, col=1)
        plotly.offline.plot(fig, filename= fig_path+ fig_desc + 'barchart_COweight.html')
        fig.write_image(fig_path+fig_desc+ '_barchartrisks_COweight.png',width=900,height=700,scale=4) 

    f+= 1