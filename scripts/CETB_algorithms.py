### CETB_analysis.py
### last edited: 5/20/18
### by: Mitch Johnson
### functions for analysis of CETB data, creates histograms and time series plots of Tb and DAV

from netCDF4 import Dataset, num2date
import numpy as np
import pandas as pd
import warnings

# getting a runtimewarning when using operators on numpy arrays with lots of NaNs, functions still perform, but using this command to suppress the warning
warnings.filterwarnings("ignore",category =RuntimeWarning)

# calculate seasonal melt onset date with the DAV/Tb Threshold algorithm. The user chooses the DAV and Tb thresholds, number (count) of melt occurrences, 
# and window of days for the algorithm to calculate the MOD.  3 occurrences of tripping Tb/DAV thresholds (252K/18K) in a 5-day (10 observation) window was previously used in Literature (Apgar/Ramage)
# the current form gets the first day of the year where any pixel in the subset experiences melt
def DAV_MOD(DAV_threshold, Tb_threshold, count, window, DAV, CETB_data, Years, cal_year, cal_date, rows_cols):
	np.errstate(invalid='ignore')	
	y_s=list(range(rows_cols[0],rows_cols[1]))
	x_s=list(range(rows_cols[2],rows_cols[3]))
	y_dims_list=list(range(len(CETB_data[0,:,0])))
	x_dims_list=list(range(len(CETB_data[0,0,:])))
	
	melt_condition_met = (DAV>DAV_threshold) & (CETB_data[:,:,:]>Tb_threshold)  #the melt condition is met when both the DAV and the Tb thresholds are exceeded
	flag = melt_condition_met.astype(int)
	matrix=pd.DataFrame()
	for i in y_dims_list:
    		for j in x_dims_list:
       			column=pd.DataFrame(data=flag[:,i,j], columns=[str(y_s[i])+','+str(x_s[j])])
        		matrix=pd.concat([matrix,column],axis=1)
	matrix=matrix.set_index(cal_date)
	shift_period=int(window/2)  #shift to get the first MOD trigger
	matrix=matrix.rolling(window, min_periods=3, center=True).sum().shift(-shift_period)
	matrix=matrix[matrix>=count]  # convert cells that do not meet criteria to NaN
	matrix=matrix.dropna(axis=0, how='all')  # deletes all rows of the dataframe that contain all NaN values, switch how='all' to how='any' to delete all rows that contain at least one NaN
	MOD=matrix.groupby(pd.Grouper(freq='A')).head(1)  #group the dataframe by year, then get the first row for that year
	MOD=MOD.dropna(axis=1, how='all')	
	return MOD  # returns a dataframe, each column is a pixel in the specified subset, each row is the algorithm-estimated seasonal melt onset date for that year

# cross-polarized gradient ratio (XPGR) melt algorithm from Abdalati and Steffen, 1995.  Threshold for Greenland is -0.0158 for SSMI
def XPGR(CETB_data, CETB_data_2):
	ymean=np.nanmean(CETB_data_2, axis=1)
	xmean=np.nanmean(CETB_data_2, axis=2)
	xmean2=np.nanmean(xmean, axis=1)
	ymean2=np.nanmean(ymean, axis=1)
	CETB_37V=np.nanmean([xmean2,ymean2],axis=0)
	CETB_19H=np.squeeze(CETB_data)	
	XPGR=(CETB_19H-CETB_37V)/(CETB_19H+CETB_37V)
	return XPGR

# dynamic-DAV from Tedesco et al 2009. This function returns a dataframe with a DAV threshold for each pixel for each year.  The threshold is calculated by taking the average
# DAV value for Jan-Feb for each year for each pixel.  - IN PROGRESS
def D_DAV(CETB_data, cal_date, cal_year, Years, rows_cols):
	y_s=list(range(rows_cols[0],rows_cols[1]))
	x_s=list(range(rows_cols[2],rows_cols[3]))
	y_dims_list=list(range(len(CETB_data[0,:,0])))
	x_dims_list=list(range(len(CETB_data[0,0,:])))

	matrix=pd.DataFrame()
	# this for loop creates a dataframe with time series of Tb for each pixel	
	for i in y_dims_list:
		for j in x_dims_list:
        		column=pd.DataFrame(data=CETB_data[:,i,j], columns=[str(y_s[i])+','+str(x_s[j])])
        		matrix=pd.concat([matrix,column],axis=1)
	matrix=matrix.set_index(cal_date)
	DAVpd=matrix.diff()  #take running difference to get DAV
	DAVpd=DAVpd.abs()  #absolute value
	DAV_monthly=DAVpd.groupby(pd.Grouper(freq='M')).mean()  #group by month and get average for each month
	DAV_monthly=DAV_monthly.dropna(axis=0, how='all')  #drop rows with all NaN values
	DAV_monthly=DAV_monthly.groupby(pd.Grouper(freq='A')).head(2)  #group by year and take the first two rows of each year (Jan-Feb)
	DAV_monthly=DAV_monthly.groupby(pd.Grouper(freq='A')).mean()  #
	DAV_monthly=DAV_monthly.set_index([Years])
	
	DAV_monthly=DAV_monthly+10
	return DAV_monthly

# winter DAV (Jan-Feb)
def Winter_DAV(CETB_data, cal_date, cal_year, Years, rows_cols):
	y_s=list(range(rows_cols[0],rows_cols[1]))
	x_s=list(range(rows_cols[2],rows_cols[3]))
	y_dims_list=list(range(len(CETB_data[0,:,0])))
	x_dims_list=list(range(len(CETB_data[0,0,:])))

	matrix=pd.DataFrame()
	# this for loop creates a dataframe with time series of Tb for each pixel	
	for i in y_dims_list:
		for j in x_dims_list:
        		column=pd.DataFrame(data=CETB_data[:,i,j], columns=[str(y_s[i])+','+str(x_s[j])])
        		matrix=pd.concat([matrix,column],axis=1)
	matrix=matrix.set_index(cal_date)
	DAVpd=matrix.diff()  #take running difference to get DAV
	DAVpd=DAVpd.abs()  #absolute value
	DAV_monthly=DAVpd.groupby(pd.Grouper(freq='M')).mean()  #group by month and get average for each month
	DAV_monthly=DAV_monthly.dropna(axis=0, how='all')  #drop rows with all NaN values
	DAV_monthly=DAV_monthly.groupby(pd.Grouper(freq='A')).head(2)  #group by year and take the first two rows of each year (Jan-Feb)
	DAV_monthly=DAV_monthly.groupby(pd.Grouper(freq='A')).mean()  #
	DAV_monthly=DAV_monthly.set_index([Years])
	
	return DAV_monthly

#End of High DAV period, gets the last day where DAV threshold and Tb threshold are both exceeded - IN PROGRESS
def end_high_DAV(DAV_threshold, Tb_threshold, count, window, DAV, CETB_data, Years, cal_year, cal_date, rows_cols):
	y_s=list(range(rows_cols[0],rows_cols[1]))
	x_s=list(range(rows_cols[2],rows_cols[3]))
	y_dims_list=list(range(len(CETB_data[0,:,0])))
	x_dims_list=list(range(len(CETB_data[0,0,:])))
	
	no_exceedance = (DAV>DAV_threshold) & (CETB_data[:,:,:]>Tb_threshold)  #the melt condition is met when both the DAV and the Tb thresholds are exceeded
	flag = no_exceedance.astype(int)
	matrix=pd.DataFrame()
	for i in y_dims_list:
    		for j in x_dims_list:
       			column=pd.DataFrame(data=flag[:,i,j], columns=[str(y_s[i])+','+str(x_s[j])])
        		matrix=pd.concat([matrix,column],axis=1)
	matrix=matrix.set_index(cal_date)
	matrix=matrix.rolling(window).sum()
	matrix=matrix[matrix>=count]  # convert cells that do not meet criteria to NaN
	matrix=matrix.dropna(axis=0, how='all')  # deletes all rows of the dataframe that contain all NaN values, switch how='all' to how='any' to delete all rows that contain at least one NaN
	EHD=matrix.groupby(pd.TimeGrouper('A')).tail(1)  #group the dataframe by year, then get the first row for that year
	EHD=EHD.dropna(axis=1, how='all')	
	return EHD


