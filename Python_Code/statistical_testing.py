#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import xarray as xr
from matplotlib import pyplot as plt
import matplotlib as mpl
import pandas as pd
import urllib
import cartopy.crs as ccrs
import os
import datetime
from matplotlib import ticker
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
from scipy.stats import t
GeoAxes._pcolormesh_patched = Axes.pcolormesh


# In[2]:


def remove_time_mean(x):
    return x - x.mean(dim='time')


# In[3]:


import os
os.chdir('/home/ivanov/matlab/Code/MJO_TW')

# Open OMI_index document and convert to a matrix
from numpy import genfromtxt
OMI_ind = genfromtxt('OMI_Index.txt', dtype = float)
OMI_amp = OMI_ind[:,6]

# Convert to polar coordinates
def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

[rho, phi] = cart2pol(OMI_ind[:,5], -OMI_ind[:,4])

# Determine MJO phase
angles = np.arange(-np.pi, 5*np.pi/4,np.pi/4)
amp_std = np.std(OMI_amp)
OMI_phase = 0*phi

for i in range(len(phi)):
    if OMI_amp[i] >= amp_std:
        if (phi[i] >= angles[0] and phi[i] < angles[1]):
            OMI_phase[i] = 1
        elif (phi[i] >= angles[1] and phi[i] < angles[2]):
            OMI_phase[i] = 2
        elif (phi[i] >= angles[2] and phi[i] < angles[3]):
            OMI_phase[i] = 3
        elif (phi[i] >= angles[3] and phi[i] < angles[4]):
            OMI_phase[i] = 4
        elif (phi[i] >= angles[4] and phi[i] < angles[5]):
            OMI_phase[i] = 5
        elif (phi[i] >= angles[5] and phi[i] < angles[6]):
            OMI_phase[i] = 6
        elif (phi[i] >= angles[6] and phi[i] < angles[7]):
            OMI_phase[i] = 7
        else:
            OMI_phase[i] = 8
    else:
        OMI_phase[i] = 0


# In[4]:


# Create datetime array
start = datetime.datetime.strptime("01-01-1979", "%d-%m-%Y")
end = datetime.datetime.strptime("08-02-2020", "%d-%m-%Y")
dates = [start + datetime.timedelta(days=x) for x in range(0, (end-start).days)]

# Create data arrays for amplitude and phase (OMI index)
OMIamp_DA = xr.DataArray(OMI_amp, coords=[dates], dims=["time"])
OMIamp_DS = xr.Dataset({"OMIamp": OMIamp_DA})

OMIphase_DA = xr.DataArray(OMI_phase, coords=[dates], dims=["time"])
OMIphase_DS = xr.Dataset({"OMIphase": OMIphase_DA})


# ## OLR

# In[ ]:


# Load OLR data
mtnlwrf_data = xr.open_mfdataset('/dx01/data/ERA5/top_net_longwave_radiation_dailymean/*.nc', parallel = True, combine = 'by_coords')

OLR_data = -1*mtnlwrf_data

# Calculate OLR anomalies
OLR_anom = OLR_data.groupby('time.month').apply(remove_time_mean)

# Add OMI amplitude and phase variables to OLR dataset
OLR_anom_OMI = xr.merge([OLR_anom,OMIamp_DS,OMIphase_DS], join='inner')

# Remove outliers
OLR_real_OMI = OLR_anom_OMI.sel(time = OLR_anom_OMI['OMIamp'] < 5)

# Group by season
MJ_OLRdata = OLR_real_OMI.where(OLR_real_OMI.time.dt.month.isin([5, 6]), drop=True)
JA_OLRdata = OLR_real_OMI.where(OLR_real_OMI.time.dt.month.isin([7, 8]), drop=True)

# Regional data
MJ_OLR_reg = MJ_OLRdata.sel(longitude = slice(40,100), latitude = slice(35,0))
JA_OLR_reg = JA_OLRdata.sel(longitude = slice(40,100), latitude = slice(35,0))

# Calculate average OLR anomaly in each phase, seasonally
MJ_OLR_mean = MJ_OLR_reg.groupby(MJ_OLR_reg.OMIphase).mean()
JA_OLR_mean = JA_OLR_reg.groupby(JA_OLR_reg.OMIphase).mean()

# Calculate STD in each phase, seasonally
MJ_OLR_std = MJ_OLR_reg.groupby(MJ_OLR_reg.OMIphase).std()
JA_OLR_std = JA_OLR_reg.groupby(JA_OLR_reg.OMIphase).std()


# #### MJ

# In[ ]:


# Calculate student t-test for each grid point within each phase

#pVal = 0.05

# Phase 0
#mean_0 = MJ_OLR_mean.sel(OMIphase = 0).mtnlwrf
#std_0 = MJ_OLR_std.sel(OMIphase = 0).mtnlwrf
#n0 = len(MJ_OLR_reg.sel(time = MJ_OLR_reg['OMIphase'] == 0).time)


# Iterate through each phase
#for phase in np.arange(1,9):
    
    ##### STUDENT T-TEST
    
   # mean_1 = MJ_OLR_mean.sel(OMIphase = phase).mtnlwrf
    #std_1 = MJ_OLR_std.sel(OMIphase = phase).mtnlwrf
    
    # Calculate number of days in each phase, seasonally
    #n1 = len(MJ_OLR_reg.sel(time = MJ_OLR_reg['OMIphase'] == phase).time)
    
    # Calculate t
    #sdelta = np.sqrt((std_1**2)/n1 + (std_0**2)/n0)
    #num = mean_1 - mean_0
    #t_stat = num/sdelta

    #df_num = ((std_1**2)/n1 + (std_0**2)/n0)**2
    #df_denom = (((std_1**2)/n1)**2)/(n1-1) + (((std_0**2)/n0)**2)/(n0-1)
    #df = df_num/df_denom
    
    #pvals = t.sf(np.abs(t_stat), df)*2 ## Two tailed t-test (since want pos and neg)
    
    ##### WILKS 2016 METHOD
    
    #N = np.sum(np.isfinite(pvals))
    
    #p_FDR = np.ravel(pvals)[np.ravel(np.isfinite(pvals))] #vectorize the field if it was a matrix
    #p_FDR.sort();p_FDR[p_FDR==1.0]=np.nan #sort the vector of p-values and reset nans (in my field p=1 was nan)

    #for i in range(N):
        #if (p_FDR[i]>((i+1)/np.float(N)*pVal*2)): #find the cut-point
            #p_FDR[i]=0  #set above the cut point to 0
            
    #p_FDR = np.nanmax(p_FDR) #find the max value ahead of the cut point, where all values were set to 0
    
    #msk = pvals < p_FDR
    
    #np.savetxt("MJ_OLR_phase" + np.str(phase) + "_sigmask.csv", msk, delimiter=",")


# #### JA

# In[ ]:


# Calculate student t-test for each grid point within each phase

pVal = 0.05

# Phase 0
mean_0 = JA_OLR_mean.sel(OMIphase = 0).mtnlwrf
std_0 = JA_OLR_std.sel(OMIphase = 0).mtnlwrf
n0 = len(JA_OLR_reg.sel(time = JA_OLR_reg['OMIphase'] == 0).time)


# Iterate through each phase
for phase in np.arange(1,9):
    
    ##### STUDENT T-TEST
    
    mean_1 = JA_OLR_mean.sel(OMIphase = phase).mtnlwrf
    std_1 = JA_OLR_std.sel(OMIphase = phase).mtnlwrf
    
    # Calculate number of days in each phase, seasonally
    n1 = len(JA_OLR_reg.sel(time = JA_OLR_reg['OMIphase'] == phase).time)
    
    # Calculate t
    sdelta = np.sqrt((std_1**2)/n1 + (std_0**2)/n0)
    num = mean_1 - mean_0
    t_stat = num/sdelta

    df_num = ((std_1**2)/n1 + (std_0**2)/n0)**2
    df_denom = (((std_1**2)/n1)**2)/(n1-1) + (((std_0**2)/n0)**2)/(n0-1)
    df = df_num/df_denom
    
    pvals = t.sf(np.abs(t_stat), df)*2 ## Two tailed t-test (since want pos and neg)
    
    ##### WILKS 2016 METHOD
    
    N = np.sum(np.isfinite(pvals))
    
    p_FDR = np.ravel(pvals)[np.ravel(np.isfinite(pvals))] #vectorize the field if it was a matrix
    p_FDR.sort();p_FDR[p_FDR==1.0]=np.nan #sort the vector of p-values and reset nans (in my field p=1 was nan)

    for i in range(N):
        if (p_FDR[i]>((i+1)/np.float(N)*pVal*2)): #find the cut-point
            p_FDR[i]=0  #set above the cut point to 0
            
    p_FDR = np.nanmax(p_FDR) #find the max value ahead of the cut point, where all values were set to 0
    
    msk = pvals < p_FDR
    
    np.savetxt("JA_OLR_phase" + np.str(phase) + "_sigmask.csv", msk, delimiter=",")
    
    print('JA OLR phase '+ np.str(phase) + ' saved.')


# ## Temperature

# In[5]:


# 2m air temperature
temp_data = xr.open_mfdataset('/dx01/data/ERA5/2mtemp_dailymean/*.nc', parallel = True, combine = 'by_coords', chunks = {'time':400})
temp_C = temp_data - 273.15 #convert from kelvin

temp_anom = temp_C.groupby('time.month').apply(remove_time_mean)

# Add amplitude and phase columns to temperature dataset
temp_anom = xr.merge([temp_anom,OMIamp_DS,OMIphase_DS], join= 'inner')

# Remove outliers
temp_real = temp_anom.sel(time = temp_anom['OMIamp'] < 5)

# Group by season and select region
MJ_T_data = temp_real.sel(longitude = slice(40,100), latitude = slice(35,0)).where(temp_real.time.dt.month.isin([5, 6]), drop=True)
JA_T_data = temp_real.sel(longitude = slice(40,100), latitude = slice(35,0)).where(temp_real.time.dt.month.isin([7, 8]), drop=True)

# Calculate average temp anomaly in each phase, seasonally
MJ_T_mean = MJ_T_data.groupby(MJ_T_data.OMIphase).mean()
JA_T_mean = JA_T_data.groupby(JA_T_data.OMIphase).mean()

# Calculate STD in each phase, seasonally
MJ_T_std = MJ_T_data.groupby(MJ_T_data.OMIphase).std()
JA_T_std = JA_T_data.groupby(JA_T_data.OMIphase).std()


# #### MJ

# In[ ]:


# Calculate student t-test for each grid point within each phase

pVal = 0.05

# Phase 0
mean_0 = MJ_T_mean.sel(OMIphase = 0).t2m
std_0 = MJ_T_std.sel(OMIphase = 0).t2m
n0 = len(MJ_T_data.sel(time = MJ_T_data['OMIphase'] == 0).time)


# Iterate through each phase
for phase in np.arange(1,9):
    
    ##### STUDENT T-TEST
    
    mean_1 = MJ_T_mean.sel(OMIphase = phase).t2m
    std_1 = MJ_T_std.sel(OMIphase = phase).t2m
    
    # Calculate number of days in each phase, seasonally
    n1 = len(MJ_T_data.sel(time = MJ_T_data['OMIphase'] == phase).time)
    
    # Calculate t
    sdelta = np.sqrt((std_1**2)/n1 + (std_0**2)/n0)
    num = mean_1 - mean_0
    t_stat = num/sdelta

    df_num = ((std_1**2)/n1 + (std_0**2)/n0)**2
    df_denom = (((std_1**2)/n1)**2)/(n1-1) + (((std_0**2)/n0)**2)/(n0-1)
    df = df_num/df_denom
    
    pvals = t.sf(np.abs(t_stat), df)*2 ## Two tailed t-test (since want pos and neg)
    
    ##### WILKS 2016 METHOD
    
    N = np.sum(np.isfinite(pvals))
    
    p_FDR = np.ravel(pvals)[np.ravel(np.isfinite(pvals))] #vectorize the field if it was a matrix
    p_FDR.sort();p_FDR[p_FDR==1.0]=np.nan #sort the vector of p-values and reset nans (in my field p=1 was nan)

    for i in range(N):
        if (p_FDR[i]>((i+1)/np.float(N)*pVal*2)): #find the cut-point
            p_FDR[i]=0  #set above the cut point to 0
            
    p_FDR = np.nanmax(p_FDR) #find the max value ahead of the cut point, where all values were set to 0
    
    msk = pvals < p_FDR
    
    np.savetxt("MJ_T_phase" + np.str(phase) + "_sigmask.csv", msk, delimiter=",")
    
    print('MJ temp phase '+ np.str(phase) + ' saved.')


# #### JA

# In[ ]:


# Calculate student t-test for each grid point within each phase

pVal = 0.05

# Phase 0
mean_0 = JA_T_mean.sel(OMIphase = 0).t2m
std_0 = JA_T_std.sel(OMIphase = 0).t2m
n0 = len(JA_T_data.sel(time = JA_T_data['OMIphase'] == 0).time)


# Iterate through each phase
for phase in np.arange(1,9):
    
    ##### STUDENT T-TEST
    
    mean_1 = JA_T_mean.sel(OMIphase = phase).t2m
    std_1 = JA_T_std.sel(OMIphase = phase).t2m
    
    # Calculate number of days in each phase, seasonally
    n1 = len(JA_T_data.sel(time = JA_T_data['OMIphase'] == phase).time)
    
    # Calculate t
    sdelta = np.sqrt((std_1**2)/n1 + (std_0**2)/n0)
    num = mean_1 - mean_0
    t_stat = num/sdelta

    df_num = ((std_1**2)/n1 + (std_0**2)/n0)**2
    df_denom = (((std_1**2)/n1)**2)/(n1-1) + (((std_0**2)/n0)**2)/(n0-1)
    df = df_num/df_denom
    
    pvals = t.sf(np.abs(t_stat), df)*2 ## Two tailed t-test (since want pos and neg)
    
    ##### WILKS 2016 METHOD
    
    N = np.sum(np.isfinite(pvals))
    
    p_FDR = np.ravel(pvals)[np.ravel(np.isfinite(pvals))] #vectorize the field if it was a matrix
    p_FDR.sort();p_FDR[p_FDR==1.0]=np.nan #sort the vector of p-values and reset nans (in my field p=1 was nan)

    for i in range(N):
        if (p_FDR[i]>((i+1)/np.float(N)*pVal*2)): #find the cut-point
            p_FDR[i]=0  #set above the cut point to 0
            
    p_FDR = np.nanmax(p_FDR) #find the max value ahead of the cut point, where all values were set to 0
    
    msk = pvals < p_FDR
    
    np.savetxt("JA_T_phase" + np.str(phase) + "_sigmask.csv", msk, delimiter=",")
    
    print('JA temp phase '+ np.str(phase) + ' saved.')


# ## Specific Humidity

# In[7]:


# Surface dewpoint
d2m_data = xr.open_mfdataset('/dx01/data/ERA5/2mdewpoint_dailymean/*.nc', parallel = True, combine = 'by_coords')
d2m_C = d2m_data.d2m - 273.15 #convert from kelvin

# Surface pressure
sp_data = xr.open_mfdataset('/dx01/data/ERA5/surface_pressure_dailymean/*.nc', parallel = True, combine = 'by_coords')
sp_mb = sp_data/100

# Specific humidity
vap_pres = 6.112*np.exp((17.67*d2m_C)/(d2m_C + 243.5))
q = (0.622 * vap_pres)/(sp_mb.sp - (0.378 * vap_pres))
q_derived = q.to_dataset(name = 'q')
q_derived

# Calculate specific humidity anomalies
q_anom_derived = q_derived.groupby('time.month').apply(remove_time_mean)

# Add to precip anom
q_anom_derived = xr.merge([q_anom_derived, OMIamp_DS, OMIphase_DS], join = 'inner')

# Remove outliers
q_real = q_anom_derived.sel(time = q_anom_derived['OMIamp'] < 5)

# Group by season and select region
MJ_qdata = q_real.sel(longitude = slice(40,100), latitude = slice(35,0)).where(q_real.time.dt.month.isin([5, 6]), drop=True)
JA_qdata = q_real.sel(longitude = slice(40,100), latitude = slice(35,0)).where(q_real.time.dt.month.isin([7, 8]), drop=True)

# Calculate average temperature anomaly in each phase, seasonally
MJ_q_mean = MJ_qdata.groupby(MJ_qdata.OMIphase).mean()
JA_q_mean = JA_qdata.groupby(JA_qdata.OMIphase).mean()

# Calculate STD in each phase, seasonally
MJ_q_std = MJ_qdata.groupby(MJ_qdata.OMIphase).std()
JA_q_std = JA_qdata.groupby(JA_qdata.OMIphase).std()


# #### MJ

# In[ ]:


# Calculate student t-test for each grid point within each phase

pVal = 0.05

# Phase 0
mean_0 = MJ_q_mean.sel(OMIphase = 0).q
std_0 = MJ_q_std.sel(OMIphase = 0).q
n0 = len(MJ_qdata.sel(time = MJ_qdata['OMIphase'] == 0).time)


# Iterate through each phase
for phase in np.arange(1,9):
    
    ##### STUDENT T-TEST
    
    mean_1 = MJ_q_mean.sel(OMIphase = phase).q
    std_1 = MJ_q_std.sel(OMIphase = phase).q
    
    # Calculate number of days in each phase, seasonally
    n1 = len(MJ_qdata.sel(time = MJ_qdata['OMIphase'] == phase).time)
    
    # Calculate t
    sdelta = np.sqrt((std_1**2)/n1 + (std_0**2)/n0)
    num = mean_1 - mean_0
    t_stat = num/sdelta

    df_num = ((std_1**2)/n1 + (std_0**2)/n0)**2
    df_denom = (((std_1**2)/n1)**2)/(n1-1) + (((std_0**2)/n0)**2)/(n0-1)
    df = df_num/df_denom
    
    pvals = t.sf(np.abs(t_stat), df)*2 ## Two tailed t-test (since want pos and neg)
    
    ##### WILKS 2016 METHOD
    
    N = np.sum(np.isfinite(pvals))
    
    p_FDR = np.ravel(pvals)[np.ravel(np.isfinite(pvals))] #vectorize the field if it was a matrix
    p_FDR.sort();p_FDR[p_FDR==1.0]=np.nan #sort the vector of p-values and reset nans (in my field p=1 was nan)

    for i in range(N):
        if (p_FDR[i]>((i+1)/np.float(N)*pVal*2)): #find the cut-point
            p_FDR[i]=0  #set above the cut point to 0
            
    p_FDR = np.nanmax(p_FDR) #find the max value ahead of the cut point, where all values were set to 0
    
    msk = pvals < p_FDR
    
    np.savetxt("MJ_q_phase" + np.str(phase) + "_sigmask.csv", msk, delimiter=",")
    
    print('MJ q phase '+ np.str(phase) + ' saved.')


# #### JA

# In[ ]:


# Calculate student t-test for each grid point within each phase

pVal = 0.05

# Phase 0
mean_0 = JA_q_mean.sel(OMIphase = 0).q
std_0 = JA_q_std.sel(OMIphase = 0).q
n0 = len(JA_qdata.sel(time = JA_qdata['OMIphase'] == 0).time)


# Iterate through each phase
for phase in np.arange(1,9):
    
    ##### STUDENT T-TEST
    
    mean_1 = JA_q_mean.sel(OMIphase = phase).q
    std_1 = JA_q_std.sel(OMIphase = phase).q
    
    # Calculate number of days in each phase, seasonally
    n1 = len(JA_qdata.sel(time = JA_qdata['OMIphase'] == phase).time)
    
    # Calculate t
    sdelta = np.sqrt((std_1**2)/n1 + (std_0**2)/n0)
    num = mean_1 - mean_0
    t_stat = num/sdelta

    df_num = ((std_1**2)/n1 + (std_0**2)/n0)**2
    df_denom = (((std_1**2)/n1)**2)/(n1-1) + (((std_0**2)/n0)**2)/(n0-1)
    df = df_num/df_denom
    
    pvals = t.sf(np.abs(t_stat), df)*2 ## Two tailed t-test (since want pos and neg)
    
    ##### WILKS 2016 METHOD
    
    N = np.sum(np.isfinite(pvals))
    
    p_FDR = np.ravel(pvals)[np.ravel(np.isfinite(pvals))] #vectorize the field if it was a matrix
    p_FDR.sort();p_FDR[p_FDR==1.0]=np.nan #sort the vector of p-values and reset nans (in my field p=1 was nan)

    for i in range(N):
        if (p_FDR[i]>((i+1)/np.float(N)*pVal*2)): #find the cut-point
            p_FDR[i]=0  #set above the cut point to 0
            
    p_FDR = np.nanmax(p_FDR) #find the max value ahead of the cut point, where all values were set to 0
    
    msk = pvals < p_FDR
    
    np.savetxt("JA_q_phase" + np.str(phase) + "_sigmask.csv", msk, delimiter=",")
    
    print('JA q phase '+ np.str(phase) + ' saved.')


# ## Wet Bulb

# In[9]:


# Surface wet bulb
TW_data = xr.open_dataset('/dx01/ivanov/data/ERA5/MJO_TW/submission1/regional_Tw_dailymean.nc')

# Calculate specific humidity anomalies
TW_anom = TW_data.groupby('time.month').apply(remove_time_mean)

# Add to precip anom
TW_anom = xr.merge([TW_anom, OMIamp_DS, OMIphase_DS], join = 'inner')

# Remove outliers
TW_real = TW_anom.sel(time = TW_anom['OMIamp'] < 5)

# Group by season
MJ_TWdata = TW_real.sel(longitude = slice(40,100), latitude = slice(35,0)).where(TW_real.time.dt.month.isin([5, 6]), drop=True)
JA_TWdata = TW_real.sel(longitude = slice(40,100), latitude = slice(35,0)).where(TW_real.time.dt.month.isin([7, 8]), drop=True)

# Calculate average temperature anomaly in each phase, seasonally
MJ_TW_mean = MJ_TWdata.groupby(MJ_TWdata.OMIphase).mean()
JA_TW_mean = JA_TWdata.groupby(JA_TWdata.OMIphase).mean()

# Calculate average temperature anomaly in each phase, seasonally
MJ_TW_std = MJ_TWdata.groupby(MJ_TWdata.OMIphase).std()
JA_TW_std = JA_TWdata.groupby(JA_TWdata.OMIphase).std()


# #### MJ

# In[ ]:


# Calculate student t-test for each grid point within each phase

pVal = 0.05

# Phase 0
mean_0 = MJ_TW_mean.sel(OMIphase = 0).Tw
std_0 = MJ_TW_std.sel(OMIphase = 0).Tw
n0 = len(MJ_TWdata.sel(time = MJ_TWdata['OMIphase'] == 0).time)


# Iterate through each phase
for phase in np.arange(1,9):
    
    ##### STUDENT T-TEST
    
    mean_1 = MJ_TW_mean.sel(OMIphase = phase).Tw
    std_1 = MJ_TW_std.sel(OMIphase = phase).Tw
    
    # Calculate number of days in each phase, seasonally
    n1 = len(MJ_TWdata.sel(time = MJ_TWdata['OMIphase'] == phase).time)
    
    # Calculate t
    sdelta = np.sqrt((std_1**2)/n1 + (std_0**2)/n0)
    num = mean_1 - mean_0
    t_stat = num/sdelta

    df_num = ((std_1**2)/n1 + (std_0**2)/n0)**2
    df_denom = (((std_1**2)/n1)**2)/(n1-1) + (((std_0**2)/n0)**2)/(n0-1)
    df = df_num/df_denom
    
    pvals = t.sf(np.abs(t_stat), df)*2 ## Two tailed t-test (since want pos and neg)
    
    ##### WILKS 2016 METHOD
    
    N = np.sum(np.isfinite(pvals))
    
    p_FDR = np.ravel(pvals)[np.ravel(np.isfinite(pvals))] #vectorize the field if it was a matrix
    p_FDR.sort();p_FDR[p_FDR==1.0]=np.nan #sort the vector of p-values and reset nans (in my field p=1 was nan)

    for i in range(N):
        if (p_FDR[i]>((i+1)/np.float(N)*pVal*2)): #find the cut-point
            p_FDR[i]=0  #set above the cut point to 0
            
    p_FDR = np.nanmax(p_FDR) #find the max value ahead of the cut point, where all values were set to 0
    
    msk = pvals < p_FDR
    
    np.savetxt("MJ_TW_phase" + np.str(phase) + "_sigmask.csv", msk, delimiter=",")
    
    print('MJ TW phase '+ np.str(phase) + ' saved.')


# #### JA

# In[ ]:


# Calculate student t-test for each grid point within each phase

pVal = 0.05

# Phase 0
mean_0 = JA_TW_mean.sel(OMIphase = 0).Tw
std_0 = JA_TW_std.sel(OMIphase = 0).Tw
n0 = len(JA_TWdata.sel(time = JA_TWdata['OMIphase'] == 0).time)

# Iterate through each phase
for phase in np.arange(1,9):
    
    ##### STUDENT T-TEST
    
    mean_1 = JA_TW_mean.sel(OMIphase = phase).Tw
    std_1 = JA_TW_std.sel(OMIphase = phase).Tw
    
    # Calculate number of days in each phase, seasonally
    n1 = len(JA_TWdata.sel(time = JA_TWdata['OMIphase'] == phase).time)
    
    # Calculate t
    sdelta = np.sqrt((std_1**2)/n1 + (std_0**2)/n0)
    num = mean_1 - mean_0
    t_stat = num/sdelta

    df_num = ((std_1**2)/n1 + (std_0**2)/n0)**2
    df_denom = (((std_1**2)/n1)**2)/(n1-1) + (((std_0**2)/n0)**2)/(n0-1)
    df = df_num/df_denom
    
    pvals = t.sf(np.abs(t_stat), df)*2 ## Two tailed t-test (since want pos and neg)
    
    ##### WILKS 2016 METHOD
    
    N = np.sum(np.isfinite(pvals))
    
    p_FDR = np.ravel(pvals)[np.ravel(np.isfinite(pvals))] #vectorize the field if it was a matrix
    p_FDR.sort();p_FDR[p_FDR==1.0]=np.nan #sort the vector of p-values and reset nans (in my field p=1 was nan)

    for i in range(N):
        if (p_FDR[i]>((i+1)/np.float(N)*pVal*2)): #find the cut-point
            p_FDR[i]=0  #set above the cut point to 0
            
    p_FDR = np.nanmax(p_FDR) #find the max value ahead of the cut point, where all values were set to 0
    
    msk = pvals < p_FDR
    
    np.savetxt("JA_TW_phase" + np.str(phase) + "_sigmask.csv", msk, delimiter=",")
    
    print('JA TW phase '+ np.str(phase) + ' saved.')

