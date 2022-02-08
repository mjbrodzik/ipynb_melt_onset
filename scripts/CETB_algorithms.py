### CETB_analysis.py
### functions for analysis of CETB data, creates histograms and time series plots of Tb and DAV

from netCDF4 import Dataset, num2date
import numpy as np
import pandas as pd
import pdb; # insert at places for breakpoints: pdb.set_trace()
import warnings

# getting a runtimewarning when using operators on numpy arrays with lots of NaNs, functions still perform, but using this command to suppress the warning
warnings.filterwarnings("ignore",category =RuntimeWarning)

# findMOD: Finds the first date in each column of the input DataFrame that meets melt criteria
# Input:
# df : DataFrame with 1 column for each pixel
#      and 1 row for each date, with pd.datetime64 index
#      Values are NaN for no melt, non-NaN when melt criteria have been triggered
#      Assumes that date rows are in chronological order, so that the
#      first non-NaN entry indicates the melt onset date (MOD)
# Output:
# DataFrame with columns to match input df columns
#      and 1 row, with the datetime value of the first row in df[col]
#      with non-NaN value
#
# Assumes non-melting dates are NaNs
# WE MAY NEED TO CHECK FOR THIS:
# Also assumes:
#    each input column has at least 1 non-NaN value
#    first index is a datetime64 obj that we can get the year from 
def findMOD(df):

    # The index (row label) of the output frame will be the integer year
    # Make a blank DataFrame for this row
    # The row hold the MOD for this year for each pixel
    myYearMOD = pd.DataFrame(index = [df.index[0].year], 
                             columns=df.columns).rename_axis(index='Year')

    # Treat each column of data as a separate entity
    # and look for the the first row with non-NaN entry
    for col in df.columns:

        # Find the places where this column satisfied the
        # melt criteria in our current algorithm settings
        # print('column is: %s' % col)
        # print(d2003[col])
        isMelting = ~np.isnan(df[col])

        # Get the first date when there was (persistent) melt
        # print(isMelting)
        # print("melt onset date:")
        myMOD = isMelting[isMelting == True].index[0]
        # print(myMOD)
    
        # Save this MOD in the output array for this pixel's column
        myYearMOD[col] = myMOD
        
    return myYearMOD


# calculate seasonal melt onset date with the DAV/Tb Threshold algorithm. 
# The user chooses the DAV and Tb thresholds, 
# number (count) of melt occurrences, 
# and window of days for the algorithm to calculate the MOD.  
# 3 occurrences of tripping Tb/DAV thresholds (e.g. 252K/18K) 
# in a 5-day (10 observation) window 
# was previously used in Literature (Apgar/Ramage)
# the current form gets the first day of the year 
# where any pixel in the subset experiences melt
def DAV_MOD(DAV_threshold, Tb_threshold, count, window, 
            DAV, CETB_data, Years, cal_year, cal_date, rows_cols):
    
    #FIXME: this is not usually a good idea to ignore errors
    #       unless you are absolutely sure they are spurious
    np.errstate(invalid='ignore')
    
    #the melt condition is met when both the DAV and the Tb thresholds are exceeded
    melt_condition_met = (DAV > DAV_threshold) & (CETB_data[:, :, :] > Tb_threshold)  
    flag = melt_condition_met.astype(int)
    
    # convert the melt condition array to a data frame with
    # date indexing and 1 column for each pixel in the subset region
    matrix = pd.DataFrame()
    for i in np.arange(rows_cols[0], rows_cols[1]):
        for j in np.arange(rows_cols[2], rows_cols[3]):
            column = pd.DataFrame(
                data=flag[:, i, j], 
                columns=["%d,%d" % (i, j)])
            matrix = pd.concat([matrix,column], axis=1)
    
    matrix.rename_axis(columns="Row,Col", inplace=True)
    matrix['date'] = np.array(cal_date)
    matrix.set_index('date', inplace=True)
    
    # Save the original data frame for output
    melt_flag_df = matrix.copy()
    
    # shift to get the first MOD trigger
    shift_period = int(window / 2)  
    matrix = matrix.rolling(window, min_periods=3, center=True).sum().shift(-shift_period)
    
    # convert cells that do not meet criteria to NaN
    matrix = matrix[matrix >= count]  
    
    # deletes all rows of the dataframe that contain all NaN values, 
    # switch how='all' to how='any' to delete all rows that contain at least one NaN
    matrix = matrix.dropna(axis=0, how='all') 
    
    # group the dataframe by year, then get the MOD for each pixel by year
    grouped = matrix.groupby(pd.Grouper(freq='A'))
    MOD_df = grouped.apply(findMOD)
    
    MOD_df.index = MOD_df.index.droplevel('date')
    
    # returns two dataframes:
    # MOD_df:
    #    each column is a pixel in the specified subset, 
    #    each row is the algorithm-estimated seasonal melt onset date at that pixel for that year
    # melt_flag_df: is is the complete data frame for all dates and subset pixels
    #    each column is a pixel in the specified subset,
    #    each row is 1 if melt conditions are met on this date, 0 otherwis
    # return MOD_df, melt_flag_df
    return MOD_df, melt_flag_df


def calc_DAV(CETB_data):

    # function takes the CETB_data that was read in read_Tb()
    # and returns the absolute value of the DAV         
    DAV=np.diff(CETB_data,n=1,axis=0)
    DAV_abs=np.absolute(DAV)

    # insert a 0 at beginning of array so same length as CETB_data for plotting together    
    DAV_abs=np.insert(DAV_abs, [0],[0], axis=0)

    return DAV_abs


# cross-polarized gradient ratio (XPGR) melt algorithm from Abdalati and Steffen, 1995.  
# Threshold for Greenland is -0.0158 for SSMI
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
	# this for loop creates a dataframe with time series of Tb for each pixel	
	matrix=pd.DataFrame()
	for i in np.arange(rows_cols[0], rows_cols[1]):
		for j in np.arange(rows_cols[2], rows_cols[3]):
			column = pd.DataFrame(
				data=CETB_data[:, i-rows_cols[0], j-rows_cols[2]], 
				columns=["%d,%d" % (i, j)])
			matrix = pd.concat([matrix,column], axis=1)
	
	matrix.rename_axis(columns="Row,Col", inplace=True)
	matrix['date'] = np.array(cal_date)
	matrix.set_index('date', inplace=True)
	  
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


