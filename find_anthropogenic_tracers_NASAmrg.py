#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
find_anthropogenic_tracers.py
    python script to find the best tracer to identify smoke mixed with urban air 
    and learning how to plot with mapbox maps. 
Created on Mon Nov  4 13:07:57 2019
updated on 01.30.2020 to use NASA toga merge and renmaed with _NASAmrg
06.03.20 - updated to load new TOGA data and combine with NASAmrg
@author: kodell
"""
#%% Load modules
import numpy as np
import pandas as pd
import plotly
import plotly.graph_objs as go
mapbox_access_token='pk.eyJ1Ijoia2FvZGVsbCIsImEiOiJjanZza3k1bGkzMHZoNDhwYjdybTYyOTliIn0.KJyzHWVzu2U087Ps215_LA'

# %% User inputs
# path to load files
# these files can be downloaded on the WECAN data archive:
# https://data.eol.ucar.edu/master_lists/generated/we-can/
# note: as of 06.03.2020 TOGA files with udpated CH2O are not available at this link
# and can be accessed by emailing Rebecca Hornbrook, rsh@ucar.edu
mrg_loadfn_path = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/WECAN_datafiles/WECAN_R2_merge/'
updated_TOGA_path = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/WECAN_datafiles/TOGA_data_05.21.20/'

# mrg filename to load
mrg_fn = mrg_loadfn_path + 'wecan-mrgTOGA-c130_merge_20180724_R2_thru20180828.ict'
# anth tracers according to Becky
anth_tracers = ['x224TrimePentane_TOGA','C2Cl4_TOGA','HFC134a_TOGA','HCFC22_TOGA']

# location and name to save out TOGA merge with anthropogenic flag
outfn_path = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/data analysis/processed datafiles/'
desc_name = 'NASA_R2_anth_tracers_R1TOGAupdate'

# paths for figures
fig_path = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/data analysis/figures/'

#%% load data
# WE-CAN TOGA merge R2 with all WE-CAN variables, but old TOGA names and HCHO numbers
TOGA_mrg_old = pd.read_csv(mrg_fn,header=340)

# TOGA limits of detection from TOGA R1 data
TOGA_LLOD_fn =  '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/WECAN_datafiles/TOGA_data_05.21.20/TOGA_LLOD_R1.csv'
TOGA_LLOD = pd.read_csv(TOGA_LLOD_fn)

# load updated TOGA R1 data and replace 
dates = ['0724','0726','0730','0731','0802','0803','0806','0808','0809','0813','0815','0816','0820','0823','0826','0828']
i=0
for RF in ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16']:
    if RF =='01':       
        updated_TOGA_fn = updated_TOGA_path + 'WECAN-TOGA_C130_2018'+dates[i]+'_R1.ict'
        TOGA_update_allflights=pd.read_csv(updated_TOGA_fn,header=103)
    else:
        updated_TOGA_fn = updated_TOGA_path + 'WECAN-TOGA_C130_2018'+dates[i]+'_R1.ict'
        updated_TOGA=pd.read_csv(updated_TOGA_fn,header=103)
        TOGA_update_allflights = pd.concat([TOGA_update_allflights,updated_TOGA],ignore_index=True)
    i+=1

# variable names change between the new and old TOGA data, so find both    
new_TOGA_vars = TOGA_update_allflights.keys()[2:]
old_TOGA_vars = TOGA_mrg_old.keys()[-88:-15]

# copy everything from old TOGA merge except old TOGA variables
TOGA_mrg = TOGA_mrg_old.drop(columns=old_TOGA_vars)

# add new TOGA data to the R2 WECAN TOGA merge
for TOGA_var in new_TOGA_vars:
    TOGA_mrg[' '+ TOGA_var] = TOGA_update_allflights[TOGA_var]
#%% replace LOD for anthropogenic tracers
LLOD_flag = -888.0

# replace LLOD with limit of detection for anthropogenic tracers
for VOC_name in anth_tracers:
    TOGA_mrg[' '+VOC_name].replace(to_replace = LLOD_flag,value = 0.5*TOGA_LLOD[VOC_name][0],inplace=True)
# also replace LLOD flag for 2methylpentane, which is used later
TOGA_mrg[' x2MePentane_TOGA'].replace(to_replace = LLOD_flag,value = 0.5*TOGA_LLOD['x2MePentane_TOGA'][0],inplace=True)

TOGA_mrg.replace(-9999999.0,np.nan, inplace=True)
TOGA_mrg.replace(-999.0,np.nan,inplace=True)
#%% plot anthropogenic tracers
for tracer in anth_tracers:
    fig = go.Figure(go.Scattermapbox(lat=TOGA_mrg[' LATITUDE'],lon=TOGA_mrg[' LONGITUDE'],
            mode='markers',
            marker=go.scattermapbox.Marker(color = TOGA_mrg[' ' + tracer],
                colorscale='Magma', size=10, colorbar=dict(thickness=20)),
            name = tracer,
            text = TOGA_mrg[' ' + tracer]))
    
    fig.update_layout(autosize=True,title_text = tracer,hovermode='closest',
        mapbox=go.layout.Mapbox(accesstoken=mapbox_access_token,
            bearing=0,center=go.layout.mapbox.Center(lat=45.92,lon=-115.07),
            pitch=0,zoom=5))    
    plotly.offline.plot(fig, filename=fig_path+  'anth_tracer_'+tracer  )

#%% now we need to find when these tracers are elevated
TOGA_mrg['anth_flag'] = [0]*TOGA_mrg[' x224TrimePentane_TOGA'].shape[0]

anth_ind1 = np.where(TOGA_mrg[' x224TrimePentane_TOGA'] > 20.0)
TOGA_mrg['anth_flag'].iloc[anth_ind1] += 1

anth_ind2 = np.where(TOGA_mrg[' C2Cl4_TOGA'] > 2.0)
TOGA_mrg['anth_flag'].iloc[anth_ind2] += 1

anth_ind3 = np.where(TOGA_mrg[' HFC134a_TOGA'] > 125.0)
TOGA_mrg['anth_flag'].iloc[anth_ind3] += 1

anth_ind4 = np.where(TOGA_mrg[' HCFC22_TOGA'] > 275.0)
TOGA_mrg['anth_flag'].iloc[anth_ind4] += 1

#%% plot anth flag
fig = go.Figure(go.Scattermapbox(
        lat=TOGA_mrg[' LATITUDE'],
        lon=TOGA_mrg[' LONGITUDE'],
        mode='markers',
        marker=go.scattermapbox.Marker(
            color = TOGA_mrg['anth_flag'],
            colorscale='Magma', size=10, colorbar=dict(thickness=20)),
        name = 'anth_flag'))

fig.update_layout(
    autosize=True,title_text = tracer,hovermode='closest',
    mapbox=go.layout.Mapbox(
        accesstoken=mapbox_access_token,bearing=0,
        center=go.layout.mapbox.Center(lat=45.92,lon=-115.07),pitch=0,zoom=5))    
plotly.offline.plot(fig, filename=fig_path+  'anth_tracer'  )

#%% plot 2methylpentane vs PM to see if this removed the clearly urban influced points 
#   in the VOC:PM relationship plot, just using AMS total PM is fine here, we will calc PM for the rest of the 
#   analysis in the next code
fig = go.Figure(go.Scatter(x=TOGA_mrg[' Total_AMS'],
                           y=TOGA_mrg[' x2MePentane_TOGA'],
                           mode='markers'))
fig.update_layout(title_text = 'with anthropogenic points')
fig.show()

anth_inds = np.where(TOGA_mrg['anth_flag']>0)[0]
TOGA_mrg_nourban = TOGA_mrg.drop(index=anth_inds)
fig = go.Figure(go.Scatter(x=TOGA_mrg_nourban[' Total_AMS'],
                           y=TOGA_mrg_nourban[' x2MePentane_TOGA'],
                           mode='markers'))
fig.update_layout(title_text = 'anthropogenic points removed')
fig.show()

#%% save file with anthropogenic flag
TOGA_mrg.to_csv(path_or_buf= outfn_path + 'TOGA_mrg_'+desc_name+'.csv')

