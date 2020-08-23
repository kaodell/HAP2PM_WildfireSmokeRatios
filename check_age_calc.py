#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
check_age_calc.py
    a python script to check the background cutoffs I'm using for my chem age estiamtes
Created on Mon Oct 21 11:31:43 2019
@author: kodell
"""
#%% load modules
import numpy as np
import pandas as pd
import plotly
from plotly import subplots
import plotly.graph_objs as go
import plotly.io as pio
pio.renderers.default = "chrome"
pio.templates.default = "seaborn"

#%% user inputs
# reaction time with OH
# OH number concentration
OH_conc = 2.0*(10**6.)
acrolein_koh = 1.96*(10**-11.) #NIST, source: Atkinson et al., 1986 review paper
acrylonitrile_koh = 4.04*(10.**-12.) #NIST, source: Harris et al., 1981
x2mefuran_koh = 7.31*(10**-11.) # Aschmann, Nishino, Arey, and Atkinson et al., 2011
# updated figures to be good for explanations of things in talks, new fig path
# fig_path = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/conferences and presentations/'
# fig_path for paper
fig_path = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/paper_figures/final/SI_final_figures/'
#%% # make date string and RF string for files
RFs = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
dates = ['0724','0726','0730','0731','0802','0803','0806','0808','0809','0813','0815','0816','0820','0823','0826','0828','0906','0910','0913']

date_strs = []
RF_strs = []
for RF in RFs:
    date_strs.append('2018' + dates[RF-1])
    if len(str(RF))<2:
        RF_str = 'RF0' + str(RF)
        RF_strs.append(RF_str)
    else:
        RF_str = 'RF' + str(RF)
        RF_strs.append(RF_str)

#%% load files 
# load TOGA mrg file
file_path = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/data analysis/sensitivity_tests/sensitivity_data/'
mrg_fn = 'TOGA_mrg_w_age__NASA_R2_anth_tracers_PTR16fix_R1TOGAupdate_95pct_50pctall_85.275_w200_CH3CN.csv'
TOGA_mrg = pd.read_csv(file_path+mrg_fn)

#WE-CAN ERs from Wade
wadeERs = pd.read_csv('/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/WECAN_datafiles/ERs_fromWade/updated4allfires_03.24.20/Speciated_Emissions_wTime+GC_species.csv')

#%% calc back
# load bk inds
CO_max = 85. # ppb
HCN_max = 275. # ppt
CH3CN_max = 200. # ppt

# first indentify smoke inds for my merge
smoke_inds1 = np.where(np.logical_and(TOGA_mrg['CO_comb']>CO_max,TOGA_mrg[' HCN_TOGA']>HCN_max))[0]
not_smoke_inds1 = np.where(np.logical_and(TOGA_mrg['CO_comb']<=CO_max,TOGA_mrg[' HCN_TOGA']<=HCN_max))[0]
rmv_smk_inds = np.where(TOGA_mrg[' CH3CN_TOGA'].iloc[smoke_inds1] <= CH3CN_max)
rmv_nsmk_inds = np.where(TOGA_mrg[' CH3CN_TOGA'].iloc[not_smoke_inds1]>CH3CN_max)
smoke_inds = np.delete(smoke_inds1,rmv_smk_inds)
nsmoke_inds = np.delete(not_smoke_inds1,rmv_nsmk_inds)

# calculate backgrounds
#acrylonitrile_bk = np.nanmean(TOGA_mrg[' Acrylonitrile_TOGA'].iloc[nsmoke_inds]) + np.nanstd(TOGA_mrg[' Acrylonitrile_TOGA'].iloc[nsmoke_inds],ddof=1)
#acrolein_bk = np.nanmean(TOGA_mrg[' Acrolein_TOGA'].iloc[nsmoke_inds]) + np.nanstd(TOGA_mrg[' Acrolein_TOGA'].iloc[nsmoke_inds],ddof=1)
#x2mefuran_bk = np.nanmean(TOGA_mrg[' x2MeFuran_TOGA'].iloc[nsmoke_inds]) + np.nanstd(TOGA_mrg[' x2MeFuran_TOGA'].iloc[nsmoke_inds],ddof=1)
acrylonitrile_bk = np.round(np.nanpercentile(TOGA_mrg[' Acrylonitrile_TOGA'].iloc[nsmoke_inds],95),1) 
acrolein_bk = np.round(np.nanpercentile(TOGA_mrg[' Acrolein_TOGA'].iloc[nsmoke_inds],95),1) 
x2mefuran_bk = np.round(np.nanpercentile(TOGA_mrg[' x2MeFuran_TOGA'].iloc[nsmoke_inds],95),1)

#%% caclulate emitted concentrations of species for each emissions pass
# emission ratios are ppb VOC/ ppb CO per email from Wade 1/22/2020
acrylonitrile_inds = np.where(wadeERs['NMOG_contributor']=='Acrylonitrile')
acrylonitrile_emit = 1000.0*wadeERs['Emission_ratio'].iloc[acrylonitrile_inds].values*wadeERs['Integrated_CO_ppb'].iloc[acrylonitrile_inds].values

acrolein_inds = np.where(wadeERs['NMOG_contributor']=='Acrolein')
acrolein_emit = 1000.0*wadeERs['Emission_ratio'].iloc[acrolein_inds].values*wadeERs['Integrated_CO_ppb'].iloc[acrolein_inds].values

x2mefuran_inds = np.where(wadeERs['NMOG_contributor']=='2-Methyl furan, 3-methyl furan')
x2mefuran_emit = 0.85*1000.0*wadeERs['Emission_ratio'].iloc[x2mefuran_inds].values*wadeERs['Integrated_CO_ppb'].iloc[x2mefuran_inds].values
# multiply by 0.85 per Wade email on 3/11/2020

#%% calc t to get to bk I use
t_acrolein = (-1./(acrolein_koh*OH_conc))*np.log(acrolein_bk/acrolein_emit)/3600.
t_acrylonitrile = (-1./(acrylonitrile_koh*OH_conc))*np.log(acrylonitrile_bk/acrylonitrile_emit)/3600.
t_x2mefuran = (-1./(x2mefuran_koh*OH_conc))*np.log(x2mefuran_bk/x2mefuran_emit)/3600.

#%% calc loss rates and plot
acrylonitrile_life = acrylonitrile_emit*np.exp(-1.)
acrolein_life = acrolein_emit*np.exp(-1.)
x2mefuran_life = x2mefuran_emit*np.exp(-1.)

acrolein_tau = (1./(OH_conc*acrolein_koh))/3600.
acrylonitrile_tau = (1./(OH_conc*acrylonitrile_koh))/3600.
x2mefuran_tau = (1./(OH_conc*x2mefuran_koh))/3600.

#%% plot emission, concentration at t=tau, and background for each age tracer
fig = subplots.make_subplots(rows=3, cols=1,shared_xaxes=True,
                                    subplot_titles = ['(a) 2-methylfuran (young)',
                                                      '(b) acrolein (medium)',
                                                      '(c) acrylonitrile (old)'])
# this is for legend
fig.add_trace(go.Scatter(x=wadeERs['Fire'].iloc[acrylonitrile_inds][:-4],
                         y=acrylonitrile_emit[:-4]/1000.0,
                         mode = 'markers',
                         name = 'emit',
                         marker_color = 'black',showlegend=True),row=3,col=1)
fig.add_trace(go.Scatter(x=wadeERs['Fire'].iloc[acrylonitrile_inds][:-4],
                         y=acrylonitrile_life[:-4]/1000.0,
                         mode = 'markers',
                         name = r't=&#964;',
                         marker_color = 'black', marker = dict(symbol=101),
                         showlegend=True),row=3,col=1)
fig.add_trace(go.Scatter(x=wadeERs['Fire'].iloc[acrylonitrile_inds][:-4],
                         y=[acrylonitrile_bk/1000.0]*(len(acrylonitrile_inds[0])-4),
                         mode = 'lines',
                         name = 'background',
                         marker_color='black',showlegend=True),row=3,col=1)

# now actual plotting
VOCnames = ['acrylonitrile','acrolein','2methylfuran']
colors = ['orange','red','purple']
inds = [acrylonitrile_inds,acrolein_inds,x2mefuran_inds]
taus = [acrylonitrile_tau,acrolein_tau,x2mefuran_tau]
bks = [acrylonitrile_bk,acrolein_bk,x2mefuran_bk]
emits = [acrylonitrile_emit,acrolein_emit,x2mefuran_emit]
lifes = [acrylonitrile_life,acrolein_life,x2mefuran_life]
i = 0
for i in range(3):
    fig.add_trace(go.Scatter(x=wadeERs['Fire'].iloc[inds[i]][:-4],
                             y=emits[i][:-4]/1000.0,
                             mode = 'markers',
                             name = VOCnames[i]+' emit',
                             marker_color = colors[i],showlegend=False),row=3-i,col=1)
    fig.add_trace(go.Scatter(x=wadeERs['Fire'].iloc[inds[i]][:-4],
                             y=lifes[i][:-4]/1000.0,
                             mode = 'markers',
                             name = VOCnames[i]+' t=tau='+str(taus[i])[:3],
                             marker_color = colors[i], marker = dict(symbol=101),
                             showlegend=False),row=3-i,col=1)
    fig.add_trace(go.Scatter(x=wadeERs['Fire'].iloc[inds[i]][:-4],
                             y=[bks[i]/1000.0]*(len(inds[i][0])-4),
                             mode = 'lines',
                             name = VOCnames[i]+' background',
                             marker_color=colors[i],showlegend=False),row=3-i,col=1)

fig.update_layout(plot_bgcolor='white')
fig.update_xaxes(showline=True, linewidth=2, linecolor='black',
                 tickfont=dict(size=15),tickangle=45)
fig.update_yaxes(title_text= 'mixing ratio [ppb]',
                 title_font=dict(size=18),row=2,col=1)

plotly.offline.plot(fig, filename=fig_path + '_check_age_tracers_OHrxn' )
fig.write_image(fig_path+'all_age_desc_R2NASAmrg_reviewerupdate.png',width=900,height=900,scale=5)  

#%% plot time to background
fig = go.Figure()
fig.add_trace(go.Box(x=t_acrylonitrile[:-4],name='acrylonitrile (old)',marker_color='orange'))
fig.add_trace(go.Box(x=t_acrolein[:-4],name='acrolein (medium)',marker_color='red'))
fig.add_trace(go.Box(x=t_x2mefuran[:-4],name='2-methylfuran (young)',marker_color='purple'))
fig.update_layout( plot_bgcolor='white',legend = dict(traceorder='reversed',x=0.55,y=0.9))
fig.update_xaxes(showline=True, linewidth=2, linecolor='black', #range = (0,300),
                 title = 'tracer decay time [hours]')
fig.update_yaxes(showline=True, linewidth=2, linecolor='black')

plotly.offline.plot(fig, filename=fig_path + '_check_age_tracers_OHrxn_time2bk' )
fig.write_image(fig_path+'t2bk_4agetracers_w200CH3CN_4paper_R2NASAmrg_reviewerupdates.png',width=600,height=500,scale=5)  

#%% pick flight and plot time series of TOGA-averaged PM colored by age
'''
RF = 'RF15'
flights = TOGA_mrg.groupby(by='RF_number')

flight_data = flights.get_group(RF)
flight_ages = flight_data.groupby(by='chem_smoke_age')
fig = plotly.subplots.make_subplots(rows=2,cols=2,
                                    subplot_titles = ['2methylfuran','acrolein',
                                                      'acrylonitrile','CO'])
age_colors = ['grey','purple','red','orange','blue']
i = 0
for age_name in ['not smoke', 'young, < 7 hours', 'med, < 36 hours','old, < 1 week', '> 1 week']:
    if age_name in flight_ages.groups.keys():
        data_plot = flight_ages.get_group(age_name)
        fig.add_trace(go.Scatter(x=data_plot['Start_UTC'],
                             y = data_plot['x2MeFuran_TOGA'],
                             mode='markers',
                             name=age_name,
                             legendgroup = age_name,
                             marker_color = age_colors[i]),
                             row=1,col=1)
        fig.add_trace(go.Scatter(x=data_plot['Start_UTC'],
                             y = data_plot['Acrolein_TOGA'],
                             mode='markers',
                             name=age_name,
                             legendgroup = age_name,
                             marker_color = age_colors[i]),
                             row=1,col=2)
        fig.add_trace(go.Scatter(x=data_plot['Start_UTC'],
                             y = data_plot['Acrylonitrile_TOGA'],
                             mode='markers',
                             name=age_name,
                             legendgroup = age_name,
                             marker_color = age_colors[i]),
                             row=2,col=1)
        fig.add_trace(go.Scatter(x=data_plot['Start_UTC'],
                             y = data_plot['CO_comb'],
                             mode='markers',
                             name=age_name,
                             legendgroup = age_name,
                             marker_color = age_colors[i]),
                             row=2,col=2)
    i += 1

# and plot flight on map

plotly.offline.plot(fig,
                    filename='/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/data analysis/figures/check_age/3mf_'+RF)

mapbox_access_token='pk.eyJ1Ijoia2FvZGVsbCIsImEiOiJjanZza3k1bGkzMHZoNDhwYjdybTYyOTliIn0.KJyzHWVzu2U087Ps215_LA'
fig = go.Figure()
i=0
for age_name in ['not smoke', 'young, < 7 hours', 'med, < 36 hours','old, < 1 week', '> 1 week']:
    if age_name in flight_ages.groups.keys():
        group = flight_ages.get_group(age_name)
        fig.add_trace(go.Scattergeo(
            lat=group['LATITUDE_start'].values,
            lon=group['LONGITUDE_stop'].values,
            mode='markers',
            marker = dict(color=age_colors[i],
            size=4,
            opacity=0.7),
            hovertext = group['CO_comb'],
            name = age_name))
    i += 1

fig.update_layout(
    geo = dict(
        scope = 'north america',
        showland = True,
        landcolor = "rgb(212, 212, 212)",
        subunitcolor = "rgb(255, 255, 255)",
        countrycolor = "rgb(255, 255, 255)",
        showlakes = True,
        lakecolor = "rgb(255, 255, 255)",
        showsubunits = True,
        showcountries = True,
        resolution = 50,
        projection = dict(
            type = 'mercator',
            rotation_lon = -100
        ),
        lonaxis = dict(
            showgrid = False,
            gridwidth = 0.5,
            range= [ -126.0, -108.0 ],
            dtick = 5
        ),
        lataxis = dict (
            showgrid = False,
            gridwidth = 0.5,
            range= [ 37.0, 49.0 ],
            dtick = 5
        )
    ))
plotly.offline.plot(fig,filename='/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/data analysis/figures/check_age/R2smoke_age_map_'+RF)  
'''
