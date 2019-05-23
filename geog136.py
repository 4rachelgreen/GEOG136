#!/usr/bin/env python3
# -*- coding: utf-8 -*-


__author__ = 'Rachel Green'
__contact__ = 'rachel.green@geog.ucsb.edu'
__copyright__ = '(c) Rachel Green 2019'

__license__ = 'MIT'
__date__ = 'Tue May 14 11:55:03 2019'
__version__ = '1.0'
__status__ = 'initial release'
__url__ = ''


"""

Name:           filename.py
Compatibility:  Python 3.7.0
Description:    Description of what program does

URL:            https://

Requires:       list of libraries required

Dev ToDo:       None

AUTHOR:         Rachel Green
ORGANIZATION:   University of California, Santa Barbara
Contact:        rachel.green@geog.ucsb.edu
Copyright:      (c) Rachel Green 2019


"""
#%% IMPORTS
import os
import numpy as np
import pandas as pd
from pandas import DataFrame as df
from pandas import Series as s
import math as m

import statsmodels.api as sm
from scipy import stats

import matplotlib.pyplot as plt
import matplotlib as mpl

import matplotlib.ticker as plticker
from datetime import datetime, timedelta
import pytz

# t-test for independent samples
from math import sqrt
from numpy.random import seed
from numpy.random import randn
from numpy import mean
from scipy.stats import sem
from scipy.stats import t

%matplotlib qt5 

#%% FUNCTIONS
# function for calculating the t-test for two independent samples
def independent_ttest(data1, data2, alpha):
	# calculate means
	mean1, mean2 = mean(data1), mean(data2)
	# calculate standard errors
	se1, se2 = sem(data1), sem(data2)
	# standard error on the difference between the samples
	sed = sqrt(se1**2.0 + se2**2.0)
	# calculate the t statistic
	t_stat = (mean1 - mean2) / sed
	# degrees of freedom
	df = len(data1) + len(data2) - 2
	# calculate the critical value
	cv = t.ppf(1.0 - alpha, df)
	# calculate the p-value
	p = (1.0 - t.cdf(abs(t_stat), df)) * 2.0
	# return everything
	return t_stat, df, cv, p


#%% MAIN
def main():

    
#%%
    def dateparse(d,t):
        dt = d + " " + t
        return pd.datetime.strptime(dt, '%Y-%m-%d %H:%M')
    
#%%
os.chdir('/Users/rgreen/Documents/UCSB/Courses/GEOG 136/sedgwick-towers/')
#os.chdir('C:/Users/grad/Documents/Courses/Geog136')

pacific = pytz.timezone('US/Pacific')

coyote_full = pd.read_csv('coyote_full.csv', header = 1, skiprows = [2], parse_dates={'DATE': ['date', 'time']}, date_parser=dateparse)
coyote_full['DATE'] = pd.to_datetime(coyote_full['DATE'])
coyote_full.insert(1, 'DATEP', coyote_full['DATE'].dt.tz_localize("GMT").dt.tz_convert('America/Los_Angeles').dt.tz_localize(None))
coyote_full.insert(2,'time', coyote_full['DATEP'].dt.time)

#to set datetime index
#coyote_full = coyote_full.set_index(pd.DatetimeIndex(coyote_full['DATE']))
coyote_biomet = pd.read_csv('coyote_biomet.csv', skiprows = [1], parse_dates={'DATE': ['date', 'time']}, date_parser=dateparse)
coyote_biomet['DATE'] = pd.to_datetime(coyote_biomet['DATE'])
coyote_biomet.insert(1, 'DATEP', coyote_biomet['DATE'].dt.tz_localize("GMT").dt.tz_convert('America/Los_Angeles').dt.tz_localize(None))
coyote_biomet.insert(2,'time', coyote_biomet['DATEP'].dt.time)

grass_full = pd.read_csv('grass_full.csv', header = 1, skiprows = [2], parse_dates={'DATE': ['date', 'time']}, date_parser=dateparse)
grass_full['DATE'] = pd.to_datetime(grass_full['DATE'])
grass_full.insert(1, 'DATEP', grass_full['DATE'].dt.tz_localize("GMT").dt.tz_convert('America/Los_Angeles').dt.tz_localize(None))
grass_full.insert(2,'time', grass_full['DATEP'].dt.time)

grass_biomet = pd.read_csv('coyote_biomet.csv', skiprows = [1], parse_dates={'DATE': ['date', 'time']}, date_parser=dateparse)
grass_biomet['DATE'] = pd.to_datetime(coyote_biomet['DATE'])
grass_biomet.insert(1, 'DATEP', coyote_biomet['DATE'].dt.tz_localize("GMT").dt.tz_convert('America/Los_Angeles').dt.tz_localize(None))
grass_biomet.insert(2,'time', grass_biomet['DATEP'].dt.time)


# =============================================================================
# datetime.strptime("2019-03-25 2:00", "%Y-%m-%d %H:%M")
# x = np.array([datetime.datetime(2013, 3, 25, i, 0) for i in range(24)])
# y = coyote_full.H[24:49]
# plot.plot(x,y)
# plt.show()
# =============================================================================


#Just plotting 3/25/2019
#H = corrected sensible heat flux
#Tau = corrected momentum flux

fig, axs = plt.subplots(2, 1, sharex = True, sharey = True)
fig.suptitle('Corrected latent heat flux 3/25/19')

axs[0].plot(coyote_full.time[52:72], coyote_full.LE[52:72])
axs[0].set_title('Coyote Bush')
axs[1].plot(grass_full.time[52:72], grass_full.LE[52:72])
axs[1].set_title('Grass Site')
axs[1].set_xlabel('Time')
axs[1].set_ylabel('LE (W m-2)')
plt.setp(axs[1].get_xticklabels(), rotation=45)    

plt.show()

fig = plt.figure()
plt.plot(coyote_full.time[53:77], coyote_full.LE[53:77], marker =None, color = 'cadetblue', linestyle = '--', label = 'Coyote Bush')
plt.plot(grass_full.time[53:77], grass_full.LE[53:77], marker =None, color = 'indigo', linestyle = '--', label = 'Grass')
plt.legend(loc = 'upper right')
plt.xlabel('Time')
plt.ylabel('LE (W m-2]')
plt.title('Corrected Latent Heat Flux 3/26/2019')



# =============================================================================
# ax.xaxis.set_major_locator(loc)
# 
# locator = AutoDateLocator()
# locator.intervald[HOURLY] = [3] 
# 
# 
# =============================================================================

#%%

#mask out nans 
mask = np.isfinite([coyote_full, grass_full]).all(axis=0)
(coyote_full[mask], coyote_full[mask])

# calculate the t test
alpha = 0.05
t_stat, df, cv, p = independent_ttest(coyote_full.LE[24:49], grass_full.LE[24:49], alpha)
print('t=%.3f, df=%d, cv=%.3f, p=%.3f' % (t_stat, df, cv, p))
# interpret via critical value
if abs(t_stat) <= cv:
	print('Accept null hypothesis that the means are equal.')
else:
	print('Reject the null hypothesis that the means are equal.')
# interpret via p-value
if p > alpha:
	print('Accept null hypothesis that the means are equal.')
else:
	print('Reject the null hypothesis that the means are equal.')
    


#%%

fig, ax2 = plt.subplots()
ax2.plot(coyote_full.DATE[24:49:], coyote_full.Tau[24:49])
plt.setp(ax2.get_xticklabels(), rotation=45)
plt.xlabel('Time')
plt.ylabel('Tau (kg m-1 s-2)')
plt.title('Coyote Bush: Corrected momentum flux 3/25/19')
    

    
#ax1.set_xlim(pd.Timestamp('2019-03-25 02:30'), pd.Timestamp('2019-03-26 02:30:00'))
#ax1.set_ticks(np.arrange(pd.Timestamp('2019-03-25 02:30'), pd.Timestamp('2019-03-26 02:30:00')))

fig, ax4 = plt.subplots()
ax4.plot(grass_full.DATE[24:49:], grass_full.Tau[24:49])
plt.setp(ax4.get_xticklabels(), rotation=45)
plt.xlabel('Time')
plt.ylabel('Tau (kg m-1 s-2)')
plt.title('Grass site: Corrected momentum flux 3/25/19')
    
#%%

#convert time zone
fmt = "%Y-%m-%d %H:%M:%S %Z%z"

# Current time in UTC
now_utc = datetime.now(timezone('UTC'))
print now_utc.strftime(fmt)

# Convert to US/Pacific time zone
now_pacific = now_utc.astimezone(timezone('US/Pacific'))
print now_pacific.strftime(fmt)

# Convert to Europe/Berlin time zone
now_berlin = now_pacific.astimezone(timezone('Europe/Berlin'))
    print now_berlin.strftime(fmt)





#%%
if __name__ == "__main__":
    main()
#%%


#%

