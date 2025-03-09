### CETB_analysis.py
### functions for analysis of CETB data, creates histograms and time series plots of Tb and DAV

from netCDF4 import Dataset, num2date
import numpy as np
import pandas as pd
import pdb; # insert at places for breakpoints: pdb.set_trace()
import warnings


# Importing necessary libraries and modules for Tb Threshold Algorithm by MB 
from scipy.signal import find_peaks, gaussian
import os
import matplotlib.pyplot as plt
import csv
import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
# End of necessary libraries for Tb Threshold Algorithm by MB 


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
    grouped = matrix.groupby(pd.Grouper(freq='YE'))
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
	DAV_monthly=DAVpd.groupby(pd.Grouper(freq='ME')).mean()  #group by month and get average for each month
	DAV_monthly=DAV_monthly.dropna(axis=0, how='all')  #drop rows with all NaN values
	DAV_monthly=DAV_monthly.groupby(pd.Grouper(freq='YE')).head(2)  #group by year and take the first two rows of each year (Jan-Feb)
	DAV_monthly=DAV_monthly.groupby(pd.Grouper(freq='YE')).mean()  #
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
	DAV_monthly=DAVpd.groupby(pd.Grouper(freq='ME')).mean()  #group by month and get average for each month
	DAV_monthly=DAV_monthly.dropna(axis=0, how='all')  #drop rows with all NaN values
	DAV_monthly=DAV_monthly.groupby(pd.Grouper(freq='YE')).head(2)  #group by year and take the first two rows of each year (Jan-Feb)
	DAV_monthly=DAV_monthly.groupby(pd.Grouper(freq='YE')).mean()  #
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



"""
Brightness Temperature (Tb) Threshold Optimization Algorithm for Snowmelt Detection
Author: Mahboubeh Boueshagh

Purpose: 
- Analyzes satellite brightness temperature data to determine optimal snowmelt thresholds
- Improves upon legacy 246K threshold using site-specific histogram analysis
- Handles different snow classifications with specialized processing parameters
"""

def extract_relevant_data(data_SIR, year, cal_year, cal_month, snow_class, Site):
    """
    Extract relevant brightness temperature data based on the snow class.
    
    Parameters:
    - data_SIR: DataFrame containing SIR brightness temperature data.
    - year: The year for which the data is being analyzed.
    - cal_year: Array or list of years corresponding to each data point in data_SIR.
    - cal_month: Array or list of months corresponding to each data point in data_SIR.
    - snow_class: The classification of snow for the site.
    - Site: The name of the site being analyzed.

    Returns:
    - data: A numpy array containing the filtered and combined brightness temperature data for the specified year. 
    """

    # If the snow class belongs to one of the specified classes ['Boreal Forest', 'Montane Forest', 'Tundra', 'Prairie']: 
    # Through empirical trials, we just use Jan to Sep of the year. By excluding October to December, the data analysis concentrates on the months when snowmelt is actively occurring (late winter to early summer). This helps to reduce noise in the data from the snow accumulation phase, which is predominant in these excluded months. The removal sharpens the focus on the transition from peak snowpack conditions to melting phases, which generally start in spring. This can help in more accurately identifying the onset of snowmelt. 
    # Focusing on the critical snowmelt months allows for more specialized handling of data variations due to environmental factors such as temperature fluctuations and solar radiation, which are more relevant to the melting processes than to the freezing or snow accumulation processes.
    # January to September covers winter, spring, and part of summer. This period includes the tail end of the Arctic winter, the entire spring thaw, and the majority of the summer melt season. Therefore, the histogram will mainly reflect the brightness temperatures associated with these specific seasonal conditions.
    # The histogram is likely to show a range of brightness temperatures that includes very cold values at the beginning of the period (reflecting frozen conditions) and warmer values towards the end as surface melting increases. 
    # By excluding the late autumn and early winter months (October to December), the histogram will not capture the re-freezing period and the onset of the snow season, which typically show lower Tb values associated with fresh snow and freezing conditions.
    # Without data from October to December, the histogram misses the period when temperatures drop, and surfaces begin to refreeze, which would normally provide a counterbalance to the melt season data, showing lower Tb values.
    
    if snow_class in ['Boreal Forest', 'Montane Forest', 'Tundra', 'Prairie']:
        print(f"  - Analyzing {Site} based on January to September of the year")
        # Selecting January to September data from the current year
        mask_curr_year = (cal_year == year) & (cal_month <= 9)
        data_curr_year = data_SIR['TB'][mask_curr_year]
        data = data_curr_year

    # For other snow classes, use full calendar year
    elif snow_class in ['Maritime', 'Ephemeral', 'Ice', 'Ocean']:
        print(f"  - Analyzing {Site} based on Calendar Year")
        data = data_SIR['TB'][(data_SIR['cal_year'] == year)]
    else:
        return None
        
    # Remove physically impossible measurements
    data = data[data > 0]
    return data

def compute_smoothed_histogram(data, snow_class, Site):
    """
    Compute and return a smoothed histogram for the brightness temperature data.
    - opt_kernel_width: Optimal kernel width for smoothing, calculated based on Silverman's rule. 
                         This rule is a commonly used method to estimate the bandwidth of a kernel density estimate.
                         The 0.4 factor is an empirically determined scaling factor specific to this project tested in a range of SNOTEL sites in AK.
                         A minimum value of 2 is set based on trial and error to ensure adequate smoothing.
    
    Parameters:
    - data: A numpy array containing brightness temperature data.
    - snow_class: The classification of snow for the site.
    - Site: The name of the site being analyzed.

    Returns:
    - hist: The computed histogram of the data.
    - bin_edges: The edges of the bins used in the histogram.
    - hist_smooth: The smoothed histogram obtained by convolution with a Gaussian kernel.
    """
    
    # Define histogram bins based on data range
    min_bin = int(np.floor(np.min(data)))
    max_bin = int(np.ceil(np.max(data)))
    bins = range(min_bin, max_bin)
    
    # Calculate the histogram
    hist, bin_edges = np.histogram(data, bins)

    # Calculate optimal kernel width for smoothing, based on Silverman's rule
    std_data = np.std(data)
    n = len(data)
    opt_kernel_width = 0.4 * 1.06 * std_data * (n ** (-1 / 5))
    if opt_kernel_width <= 1:
        print(f"Warning: Kernel width is too small for site {Site}. Using a minimum kernel width of 1.")
        opt_kernel_width = 2  # Use a minimum kernel width
    
    # Compute the smoothed histogram using convolution with a Gaussian kernel
    kernel = gaussian(int(opt_kernel_width), opt_kernel_width)
    hist_smooth = np.convolve(hist, kernel, mode='same')

    return hist, bin_edges, hist_smooth

def analyze_histogram(hist, bin_edges, hist_smooth, Site, year, sensor_SIR, channel_SIR, snow_class, ThresholdDir):
    """
    Analyze the smoothed histogram to identify peaks, valleys, and determine temperature thresholds.
    - threshold = np.percentile(hist, 70): This threshold is set to the 70th percentile of the histogram values
                                            to focus on the higher end of the temperature distribution, as determined 
                                            through experimentation across various SNOTEL sites in Alaska.
    - distance=16&14: These values are used in peak detection. The 'distance' parameter specifies the required 
                      minimum horizontal distance (in number of bins) between neighboring peaks. The values 16 and 14 
                      have been determined through trial and error to effectively distinguish peaks in different 
                      snow classes for Alaskan sites.
    """
    
    # Initial peak detection using 70th percentile threshold
    threshold = np.percentile(hist, 70)

    # Apply snow-class specific peak detection parameters
    if snow_class in ['Boreal Forest', 'Montane Forest', 'Tundra', 'Prairie']:
        peaks_idx, _ = find_peaks(hist_smooth, distance=16, height=threshold)
    elif snow_class in ['Maritime', 'Ephemeral', 'Ice', 'Ocean']:
        peaks_idx, _ = find_peaks(hist_smooth, distance=14, height=threshold)

    num_peaks = len(peaks_idx)

    # Process bimodal case - ideal for threshold detection
    if num_peaks == 2:
        print(f"2 peaks (bimodal) found for {Site} in year {year[0]}.")
        print(f"X-axis values of peaks: {bin_edges[peaks_idx]}")
        
        # Find optimal threshold in valley between peaks
        peak1, peak2 = sorted(peaks_idx[np.argsort(hist_smooth[peaks_idx])[-2:]])
        valley = hist_smooth[peak1:peak2]
        valley_min_idx = np.argmin(valley)
        absolute_min_idx = peak1 + valley_min_idx
        min_temperature = bin_edges[absolute_min_idx]  # Optimized Tb threshold

        plot_histogram(hist, bin_edges, hist_smooth, min_temperature, Site, year, 
                      snow_class, sensor_SIR, channel_SIR, ThresholdDir, peaks_idx=peaks_idx)
        return {'Site': Site, 'snow_class': snow_class, 'threshold': min_temperature}

    # Handle single peak case - requires manual review
    elif num_peaks < 2:
        print(f"Less than two peaks (unimodal) found for {Site} in year {year[0]}.")
        min_temperature = None
        plot_histogram(hist, bin_edges, hist_smooth, min_temperature, Site, year,
                      snow_class, sensor_SIR, channel_SIR, ThresholdDir, peaks_idx=peaks_idx)
        return {'threshold': [Site, "N/A"]}

    # Handle multiple peaks case - attempt refined detection
    else:
        print(f"More than 2 peaks found for {Site} in year {year[0]}.")

        # Use mean-based threshold for secondary peak detection
        meanHist = np.mean(hist)
        StdHist = np.std(hist)
        threshold2 = meanHist + 0*StdHist  # Mean value works efficiently for vertical constraint
        peaks_idx, _ = find_peaks(hist_smooth, distance=14, height=threshold2)
        print(f"X-axis values of peaks: {bin_edges[peaks_idx]}")
        
        if len(peaks_idx) < 3:
            print(f"Complex histogram: algorithm could not find the exact 3 peaks for {Site} in year {year}.")
            return {'threshold': [Site, "N/A"]}

        else:
            # Analyze three most prominent peaks
            peak1, peak2, peak3 = peaks_idx[np.argsort(hist_smooth[peaks_idx])[-3:]]
            
            # Find optimal threshold between peaks
            valley1 = hist_smooth[peak1:peak2]
            valley2 = hist_smooth[peak2:peak3]
            combined_valley = np.concatenate((valley1, valley2))
            
            valley_min_idx = np.argmin(combined_valley)
            
            # Determine valley location and adjust index
            if valley_min_idx < len(valley1):
                absolute_min_idx = peak1 + valley_min_idx
            else:
                absolute_min_idx = peak2 + (valley_min_idx - len(valley1))
            min_temperature = bin_edges[absolute_min_idx]  # Optimized Tb threshold
            
            plot_histogram(hist, bin_edges, hist_smooth, min_temperature, Site, year,
                         snow_class, sensor_SIR, channel_SIR, ThresholdDir, peaks_idx=peaks_idx)
            return {'Site': Site, 'snow_class': snow_class, 'threshold': min_temperature}

def plot_histogram(hist, bin_edges, hist_smooth, min_temperature, Site, year, sensor_SIR, channel_SIR, snow_class, ThresholdDir, peaks_idx=None):
    """
    Plot the histogram and smoothed curve, highlight peaks and threshold, and save the plot.
    
    Parameters:
    - hist: The histogram of brightness temperature data.
    - bin_edges: The edges of the bins used in the histogram.
    - hist_smooth: The smoothed histogram.
    - min_temperature: The minimum temperature threshold determined by the analysis.
    - Site: The name of the site.
    - year: The year for which the data is being analyzed.
    - sensor_SIR: Sensor information for the data.
    - channel_SIR: Channel information for the data.
    - snow_class: The classification of snow for the site.
    - ThresholdDir: Directory where the plot will be saved.
    - peaks_idx: Indices of detected peaks in the smoothed histogram (optional).

    Saves:
    - A plot of the histogram and smoothed curve as a PNG file.
    """

    save_filename = f"{Site}_{year[0]}_{sensor_SIR}.png"
    save_path = os.path.join(ThresholdDir, save_filename)

    # Create and configure plot
    fig, ax = plt.subplots()
    ax.set_title('All data SIR Histogram (' + str(year[0]) + ') ' + Site + ' (' + snow_class + ')' + ' for ' + sensor_SIR + ' ' + channel_SIR)
    ax.set_xlabel('Brightness Temp (K)')
    ax.hist(bin_edges[:-1], bins=bin_edges, weights=hist, alpha=0.5, label='Histogram')
    
    # Add appropriate visualization based on peak detection results
    if len(peaks_idx) >= 2:
        # Show fitted curve and thresholds for multi-peak cases
        ax.plot(bin_edges[:-1], hist_smooth, label='Fitted curve')
        ax.axvline(x=min_temperature, color='black', linestyle='--', 
                  label=f'New optimized Tb threshold at {min_temperature}K')
        Tb_threshold = 246
        ax.axvline(x=Tb_threshold, color='red', linestyle='--', 
                  label=f'Legacy Tb threshold at {Tb_threshold}K')
    else:
        # Add warning text for single-peak cases
        min_temperature = None
        ax.text(0.5, 0.5, 'Histogram has only 1 peak, check the histogram visually & manually for thresholds',
                horizontalalignment='center', verticalalignment='center',
                transform=ax.transAxes, fontsize=12)
    
    # Adjust plot appearance
    max_hist = np.max(hist)
    ax.set_ylim([0, max_hist + 100])  # Add space above highest bar
    plt.legend(loc='upper left')
    plt.savefig(save_path)
    plt.show()
