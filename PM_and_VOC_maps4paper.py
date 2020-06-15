#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
PM_and_VOC_maps4paper.py
    a python script to make map figures for paper using kriging data and ratios
Created on Mon Feb 10 09:24:06 2020
@author: kodell
"""
#%% import modules
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import matplotlib.pyplot as pl
import matplotlib as mplt
from mpl_toolkits.basemap import Basemap
import shapefile
import pylab as plb
from matplotlib.colors import BoundaryNorm

#%% user imputs
years = np.arange(2006,2019)
states_plot = ['Colorado','California','Idaho','Montana','Wyoming','Washington',
               'Oregon','Nevada','Utah']


# kriging data path, HMS files are also here
kdata_path = '/Users/kodell/Local Google Drive /CSU/Research/NASA_fires/kriging data/PM_4HAPs_paper/'

# ratio file name
ratio_data_fn = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/health risk tables/VOC2PM_ratios_4paper_NASAmrg_R1TOGAupdate.csv'

# file path for health risk factors and health VOC names
HAPs_RFs_fn = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/health risk tables/HAPs_EFs_RFs_forcode_final_final_WadeUpdate.csv'
HAPs_RFs_raw = pd.read_csv(HAPs_RFs_fn,header=0,skiprows=[1,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,
                                                          49,50,51,52])
HAPs_RFs = HAPs_RFs_raw.set_index('TOGA name')

healthVOCs = HAPs_RFs.index

output_figure_path = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/paper_figures/final/'
fig_desc = 'PTR16fix_R1TOGAupdate'

# sensitivity tests
# leave in negative enhancements
#output_figure_path = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/data analysis/sensitivity_tests/sensitivity_figures/leave in negative enhancements/'
#fig_desc = '_leave_negs'
#ratio_data_fn = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/health risk tables/VOC2PM_ratios_sensitivity_leave_negs.csv'

# background tests
#output_figure_path = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/data analysis/sensitivity_tests/sensitivity_figures/background_sensitivity/'
#fig_desc = '_allbk_minusonesig'
#ratio_data_fn = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/health risk tables/VOC2PM_ratios_4paper_NASAmrg_allbk_minus1sig.csv'

#output_figure_path = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/data analysis/sensitivity_tests/sensitivity_figures/detection_limit_sensitivity/'
#fig_desc = '_NASA_R2mrg_85.275_w200_CH3CN_addPTRcheck_04bdl'
#ratio_data_fn = '/Users/kodell/Local Google Drive /CSU/Research/NSF/WECAN PM_VOC health/health risk tables/VOC2PM_ratios_4paper_NASAmrg_addPTRcheck_04bdl.csv'

#%% user-defined functions
def make_map(states2plot):
    # draw base map
    m = Basemap(llcrnrlon=-125,llcrnrlat=32,urcrnrlon=-102,urcrnrlat=50,
        projection='merc',lat_1=33,lat_2=45,lon_0=-95,resolution='l')
    # draw states
    m.readshapefile('cb_2018_us_state_20m/cb_2018_us_state_20m','states',
                    drawbounds=False,linewidth=0.45,color='gray')
    for info, shape in zip(m.states_info,m.states):
        #print info['NAME']
        if info['NAME'] in states2plot:
            x,y = zip(*shape)
            m.plot(x,y,marker=None,color='gray',linewidth=0.45)
    pl.box(on=None)

    return m

#%% parameters and constants
P_std = 101325 #Pa
T_std = 273.15 #K
R = 8.314 #j/molK

#%% Load Files
# ratios of VOCs to PM
ratios = pd.read_csv(ratio_data_fn)
ratios.set_index('variable name',inplace=True)
# Load kPM and HMS
# let's load files and average in this loop
anavg_smokePM_allyears = np.empty([1,189,309])
anavg_smokePM_allyears[:] = -999

smokePM_daily_allyears = np.empty([1,366,189,309])
smokePM_daily_allyears[:] = -999

# when calculating smoke PM, check if it makes a significant difference
# to remove negative nosmoke background values. They are small negative 
# numbers < 1 micron, but could make a difference on the average and
# make the smoke artifically high

for year in years:
    # load kriging data
    krig_nc_fid = Dataset(kdata_path + 'krigedPM25_' + str(year) + '.nc')
    PM25 = krig_nc_fid['PM25'][:]
    print PM25.max()
    nosmokePM25 = krig_nc_fid['Background PM25'][:]
    glat = krig_nc_fid['lat'][:]
    glon = krig_nc_fid['lon'][:]
    krig_nc_fid.close()
    
    # load HMS to multiply by smoke PM so smoke PM is zero on non-smoke days
    HMS_fid = Dataset(kdata_path + '/HMS/hms_smoke_'+str(year)+'.nc')
    HMS_smoke = HMS_fid['HMS_Smoke'][:]
    HMS_fid.close()
    
    # remove mask (mask=False for all arrays so can just take values)
    PM25_vals = PM25.data    
    nosmokePM25_vals = nosmokePM25.data
    HMS_smoke_vals = HMS_smoke.data
    
    smokePM = HMS_smoke_vals*np.array([PM25_vals - nosmokePM25_vals])    
    
    if year in [2006,2007,2009,2010,2011,2013,2014,2015,2017,2018]:
        nan_array = np.zeros([1,189,309])
        nan_array[:] = np.nan
        smokePM_ly_adj = np.insert(smokePM,59,nan_array,axis=1)
        smokePM_stack = smokePM_ly_adj
    else:
        smokePM_stack = smokePM
        
    smokePM_daily_allyears = np.vstack([smokePM_daily_allyears,smokePM_stack])    
    
    anavg_smokePM = np.nanmean(smokePM, axis=1)
    # where annual average smokePM < 0, make 0 ... likely in places w few monitors,
    # little smoke impact, or smoke is often not at the surface
    neginds = np.where(anavg_smokePM < 0)
    print len(neginds[0])
    anavg_smokePM[neginds] = 0 
    anavg_smokePM_allyears = np.vstack([anavg_smokePM_allyears, anavg_smokePM])    

    # just to ensure we don't stack the same year twice 
    del smokePM_stack

# remove first array of -999s
smokePM_daily_allyears = smokePM_daily_allyears[1:,:,:,:]
anavg_smokePM_allyears = anavg_smokePM_allyears[1:,:]

# now, calculate average smokePM 2006-2018
decadal_avg_smokePM = np.nanmean(smokePM_daily_allyears,axis=(0,1))
# there are negatives here, but once we pull for the western US/over land, 
# they go away.

# also want 2018 as a case for 'extreme year' 
PM2018 = anavg_smokePM_allyears[-1,:,:]

#%% Convert units of chronic cancer risk ratios and noncancer hazard values to ppt
# Convert units of chronic cancer risk ratios from ug/m3 -> ppt
# conversion factor converts mw to g/molec and will convert 1/ug/m3 to ppt
conv_factor = ((P_std)/(T_std*R))*(1.*10.**-12.)*(1.*10.**6)*(1.*10.**6) #second 10**6 so number is ppl/mil risk

crf_ppt = HAPs_RFs['chronic cancer risk factor']*HAPs_RFs['molecular weight']*conv_factor   
HAPs_RFs['crf_ppt'] = crf_ppt

arl_ppt = HAPs_RFs['acute risk factor']*(1./(1000.*HAPs_RFs['molecular weight']))*(R*T_std/P_std)*(1.*10.**12)
HAPs_RFs['arl_ppt'] = arl_ppt

crl_ppt = HAPs_RFs['chronic noncancer risk factor']*(1./(1000.*HAPs_RFs['molecular weight']))*(R*T_std/P_std)*(1.*10.**12)
HAPs_RFs['crl_ppt'] = crl_ppt

  
#%% Prep data for plotting on maps
# create the figure and axes instances.
# adjust center lat/lons to llcorner lat lons for pcolormesh
# adjust xlon,xlat values so they represent corners of grid cells for mapping using pcolor
# calculate average between two points and reassign lat/lon pairs
# skip first row and column of data since there are no points outside of domain to average with

# switch these to lower left corner lat/lon (glat/glon are centers)
# kind of complicated with kriging grid
dx = 0.5*(glon[1:,1:]-glon[:-1,:-1])
dy = 0.5*(glat[1:,1:]-glat[:-1,:-1])

hfdx = np.vstack([dx[[0],:],dx])
fdx = np.hstack([hfdx[:,[0]],hfdx])

hfdy = np.vstack([dy[[0],:],dy])
fdy = np.hstack([hfdy[:,[0]],hfdy])

pglon = glon-fdx
pglat = glat-fdy

# mask data outside of specified states (this is extrapolation with kriging)
full_mask = np.zeros(glon.shape)
US_file = shapefile.Reader('cb_2018_us_state_20m/cb_2018_us_state_20m')
US_shp = US_file.shape(0)
state_records = US_file.records()
state_shapes = US_file.shapes()
for j in range(len(state_records)):
    #print state_records[j][-4]
    if state_records[j][-4] in states_plot:
        state_shp = state_shapes[j]
        for i in range(len(state_shp.parts)):
            i0 = state_shp.parts[i]
            if i < len(state_shp.parts)-1:
                i1 = state_shp.parts[i+1] - 1
            else:
                i1 = len(state_shp.points)
            seg = state_shp.points[i0:i1+1]
            mpath = mplt.path.Path(seg)
            points = np.array((glon.flatten(), glat.flatten())).T
            mask = mpath.contains_points(points).reshape(glon.shape)
            mask_int = np.array(mask,dtype=int)
            full_mask += mask_int
US_decadal_avg_smokePM = np.where(full_mask>0,decadal_avg_smokePM, np.nan)
US_PM2018 = np.where(full_mask>0,PM2018, np.nan)

#%% calculate HAPS hazard quotient and cancer risk for the decade and 2018
chronicHAPs_risk_decade = np.zeros(decadal_avg_smokePM.shape)
cancerHAPs_risk_decade = np.zeros(decadal_avg_smokePM.shape)
chronicHAPs_risk_2018 = np.zeros(US_PM2018.shape)
cancerHAPs_risk_2018 = np.zeros(US_PM2018.shape)

for VOC_name in healthVOCs:
    VOC2PMratio = ratios['medium VOC PM1 ratio, median [ppt/ugm-3]'].loc[VOC_name]
    VOC_decade = VOC2PMratio*US_decadal_avg_smokePM
    VOC_2018 = VOC2PMratio*US_PM2018

    if np.isfinite(HAPs_RFs['crl_ppt'][VOC_name]):
        VOC_risk_decade = VOC_decade/HAPs_RFs['crl_ppt'][VOC_name]
        chronicHAPs_risk_decade += VOC_risk_decade

        VOC_risk_2018 = VOC_2018/HAPs_RFs['crl_ppt'][VOC_name]
        chronicHAPs_risk_2018 += VOC_risk_2018


    if np.isfinite(HAPs_RFs['crf_ppt'][VOC_name]):
        VOC_cancer_risk_decade = VOC_decade*HAPs_RFs['crf_ppt'][VOC_name]
        cancerHAPs_risk_decade += VOC_cancer_risk_decade    

        VOC_risk_2018 = VOC_2018*HAPs_RFs['crf_ppt'][VOC_name]
        cancerHAPs_risk_2018 += VOC_risk_2018

#%% figure 2 or 3: full time period PM
#Hack to fix missing PROJ4 env var
import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

#subplot 1: average 2006-2018 annual PM from fires
pl.figure()
m = make_map(states_plot)
# plot data
cs = m.pcolormesh(pglon,pglat, US_decadal_avg_smokePM,
                  latlon=True,vmin=0,vmax=2,cmap='YlOrRd')
# add colorbar
m.colorbar(cs,"bottom", size="5%", pad='2%',extend='max',label = r'PM$_{2.5}$ [$\mu$g m$^{-3}$]')
pl.title('2006-2018 Average Smoke PM$_{2.5}$ ')
pl.savefig(output_figure_path +'PM_map_allyears'+fig_desc+'.png',
           dpi = 400)

#%% Figure 4
# create colormaps
cmapPM = plb.cm.get_cmap('Greys',16)
cmapPM.set_over('black')
cmapPM.set_under('white')
normPM = BoundaryNorm([0,0.25, 0.5, 0.75, 1,1.5,2,2.5,3,4,5], ncolors=cmapPM.N)

cmapA = plb.cm.get_cmap('PuRd',16)
cmapA.set_over('black')
cmapA.set_under('white')
normA = BoundaryNorm([0,0.025, 0.05, 0.075, 0.1,0.2,0.3,0.4,0.5,1], ncolors=cmapA.N)

cmapB = plb.cm.get_cmap('BuPu',8)
cmapB.set_under('white')
normB = BoundaryNorm([0,1,2,3,4,5,10,20], ncolors=cmapB.N)

# make figure
fig = pl.figure()
# a) decade PM
ax = fig.add_subplot(231)
ax.text(-0.15, 0.93, '(a)',transform=ax.transAxes,fontsize=10,va='top')
m = make_map(states_plot)
cs = m.pcolormesh(pglon,pglat, US_decadal_avg_smokePM,latlon=True,
                  cmap=cmapPM, norm=normPM)

# b) decade hazard index
ax = fig.add_subplot(232)
ax.text(-0.15, 0.93, '(b)',transform=ax.transAxes,fontsize=10,va='top')
ax.set_title('\n2006-2018 average exposure')
m = make_map(states_plot)
cs = m.pcolormesh(pglon,pglat, chronicHAPs_risk_decade,latlon=True,
                  vmax=1,cmap=cmapA, norm=normA)

# c) decade cancer risk
ax = fig.add_subplot(233)
ax.text(-0.15, 0.93, '(c)',transform=ax.transAxes,fontsize=10,va='top')
m = make_map(states_plot)
cs = m.pcolormesh(pglon,pglat, cancerHAPs_risk_decade,latlon=True,
                  vmax=25,cmap=cmapB,norm=normB)

# d) 2018 PM
ax = fig.add_subplot(234)
ax.text(-0.15, 0.93, '(d)',transform=ax.transAxes,fontsize=10,va='top')
m = make_map(states_plot)
cs = m.pcolormesh(pglon,pglat, US_PM2018,latlon=True,
                  cmap=cmapPM,norm=normPM)
m.colorbar(cs,"bottom", size="5%", pad='2%',extend='max',
           label = r'PM$_{2.5}$ [$\mu$g m$^{-3}$]')

# e) 2018 hazard index
ax = fig.add_subplot(235)
ax.text(-0.15, 0.93, '(e)',transform=ax.transAxes,fontsize=10,va='top')
ax.set_title('2018 as representative exposure')
m = make_map(states_plot)
cs = m.pcolormesh(pglon,pglat, chronicHAPs_risk_2018,latlon=True,
                  cmap=cmapA,norm=normA)
colorbar = m.colorbar(cs,"bottom", size="5%", pad='2%',extend='max',
                      label = 'hazard index')
colorbar.set_ticks([0,0.05,0.1,0.3,0.5,1])
colorbar.ax.set_xticklabels(['0','0.05','0.1','0.3','0.5','1']) 

# f) 2018 cancer risk
ax = fig.add_subplot(236)
ax.text(-0.15, 0.93, '(f)',transform=ax.transAxes,fontsize=10,va='top')
m = make_map(states_plot)
cs = m.pcolormesh(pglon,pglat, cancerHAPs_risk_2018,latlon=True,
                  vmax=25,cmap=cmapB,norm=normB)
m.colorbar(cs,"bottom", size="5%", pad='2%', 
           label='cancer risk per 10$^6$ exposed')
# adjust figure white space and save
pl.subplots_adjust(top=0.969,bottom=0.1,left=0.05,right=0.956,
                   hspace=0.0,wspace=0.4)
pl.savefig(output_figure_path + 'HAPs_hazard_and_risk'+fig_desc+'.png',dpi=300)
