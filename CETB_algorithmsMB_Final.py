### CETB_analysis.py
### functions for analysis of CETB data, creates histograms and time series plots of Tb and DAV

from netCDF4 import Dataset, num2date
import numpy as np
import pandas as pd
import pdb; # insert at places for breakpoints: pdb.set_trace()
import warnings


# Importing necessary libraries and modules for Tb Threshold Algorithm by MB 
from scipy.signal import find_peaks, gaussian
from scipy.ndimage import gaussian_filter1d
from scipy.stats import norm
import os
import matplotlib.pyplot as plt
import csv
import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
log = logging.getLogger(__name__)
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

    # If the snow class belongs to one of the specified classes:
# ['Boreal Forest', 'Montane Forest', 'Tundra', 'Prairie']
#
# For high-latitude regions such as Alaska, empirical testing shows that using
# data from January to September provides the most effective temporal window
# for snowmelt detection. This range captures the complete progression from
# frozen winter conditions through spring thaw and into the main melt season.
#
# By excluding October to December, the analysis avoids periods dominated by
# snow accumulation and re-freezing, which introduce noise unrelated to melt
# processes. This refines the histogram’s focus on the transition from peak
# snowpack to active melt, improving the detection of melt onset thresholds.
#
# Within the January–September window:
#   - Brightness temperatures transition from very cold values (frozen state)
#     to warmer values (active melt and surface thaw).
#   - This interval covers late winter, spring, and early-to-mid summer,
#     when diurnal temperature variation and solar radiation strongly
#     influence melt dynamics.
#
# For lower-latitude Western U.S. regions (e.g., Rockies, Sierra Nevada,
# and Cascades), where snowpack develops and melts earlier, the time
# window should be adjusted to February–May to better capture the local
# melt period.
#
# Note: Maritime and Ocean snow classes are excluded from this process,
# as their thermal and hydrological regimes differ significantly from
# inland or high-latitude snow environments.

    
    #### For Alaska
    if snow_class in ['Boreal Forest', 'Montane Forest', 'Tundra', 'Prairie']:
        #print(f"  - Analyzing {Site} based on January to September of the year")
        # Selecting January to September data from the current year
        # mask_curr_year = (cal_year == year) & (cal_month <= 9)
        mask_curr_year = np.isin(cal_year, year) & (cal_month <= 9)
        data = data_SIR['TB'][mask_curr_year]

    # For other snow classes, use full calendar year
    elif snow_class in ['Maritime', 'Ephemeral', 'Ice', 'Ocean']: 
        #print(f"  - Analyzing {Site} based on Calendar Year")
        mask = np.isin(cal_year, year)
        data = data_SIR['TB'][mask]
    else:
        return None

    #### For lower latitudes 
    
    #print(f"  - Analyzing {Site} based on Feb - May of the year") 
    #mask = np.isin(cal_year, year) & (cal_month >= 2) & (cal_month <= 5)
    #data = data_SIR['TB'][mask]
    #data = data_SIR['TB'][(data_SIR['cal_year'] == year)]


    
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
        #print(f"Warning: Kernel width is too small for site {Site}. Using a minimum kernel width of 2.")
        opt_kernel_width = 2  # Use a minimum kernel width
    
    # Compute the smoothed histogram using convolution with a Gaussian kernel
    kernel = gaussian(int(opt_kernel_width), opt_kernel_width)
    hist_smooth = np.convolve(hist, kernel, mode='same')

    return hist, bin_edges, hist_smooth

def analyze_histogram(hist, bin_edges, hist_smooth, data, Site, year, sensor_SIR, channel_SIR, snow_class, ThresholdDir):

    # ⛔ Skip ocean and other snow classes — no fallback needed, None is correct
    skip_classes = ["Ocean", "Maritime"]
    if snow_class in skip_classes:
        return {
            'Site': Site,
            'snow_class': snow_class,
            'threshold': None
        }

    # --- helper: single entry point for segmentation fallback ---
    def try_segmentation_fallback():
        # For Ephemeral, the histogram is often genuinely unimodal
        # at lower latitudes — don't force a threshold from noise
        SKIP_FALLBACK_CLASSES = ['Ephemeral']
        if snow_class in SKIP_FALLBACK_CLASSES:
            return "N/A"

        seg = analyze_and_segment_histogram(
            data=data,
            Site=Site,
            year=year,
            snow_class=snow_class,
            ThresholdDir=ThresholdDir,
            bin_range=(150, 300),
            plot=False
        )
        if not seg:
            return "N/A"

        t = seg.get('threshold_intersection', None)
        if t is not None and np.isfinite(float(t)):
            return int(round(float(t)))

        t = seg.get('threshold_used', None)
        if t is not None and np.isfinite(float(t)):
            return int(round(float(t)))

        return "N/A"


    # Initial peak detection using 70th percentile threshold
    threshold = np.percentile(hist, 70)

    # Apply snow-class specific peak detection parameters
    if snow_class in ['Boreal Forest', 'Montane Forest', 'Tundra', 'Prairie']:
        peaks_idx, _ = find_peaks(hist_smooth, distance=16, height=threshold)
    elif snow_class in ['Ephemeral', 'Ice']:
        peaks_idx, _ = find_peaks(hist_smooth, distance=14, height=threshold)
    else:
        # Unknown snow class — attempt segmentation as best effort
        return {
            'Site': Site,
            'snow_class': snow_class,
            'threshold': try_segmentation_fallback()
        }

    num_peaks = len(peaks_idx)

    # ----------------------------------------------------------------
    # CASE 1: Bimodal — ideal, direct valley detection
    # ----------------------------------------------------------------
    if num_peaks == 2:
        peak1, peak2 = sorted(peaks_idx[np.argsort(hist_smooth[peaks_idx])[-2:]])
        valley = hist_smooth[peak1:peak2]
        valley_min_idx = np.argmin(valley)
        absolute_min_idx = peak1 + valley_min_idx
        min_temperature = bin_edges[absolute_min_idx]
        return {'Site': Site, 'snow_class': snow_class, 'threshold': min_temperature}

    # ----------------------------------------------------------------
    # CASE 2: Unimodal — try segmentation fallback
    # ----------------------------------------------------------------
    elif num_peaks < 2:
        fallback = try_segmentation_fallback()
        return {'Site': Site, 'snow_class': snow_class, 'threshold': fallback}

    # ----------------------------------------------------------------
    # CASE 3: More than 2 peaks — try refined 3-peak detection,
    #         then segmentation fallback
    # ----------------------------------------------------------------
    else:
        meanHist = np.mean(hist)
        threshold2 = meanHist  # Mean value works efficiently for vertical constraint
        peaks_idx, _ = find_peaks(hist_smooth, distance=14, height=threshold2)

        if len(peaks_idx) == 3:
            three_peaks = peaks_idx[np.argsort(hist_smooth[peaks_idx])[-3:]]
            peak1, peak2, peak3 = np.sort(three_peaks)

            valley1 = hist_smooth[peak1:peak2]
            valley2 = hist_smooth[peak2:peak3]
            combined_valley = np.concatenate((valley1, valley2))
            valley_min_idx = np.argmin(combined_valley)

            if valley_min_idx < len(valley1):
                absolute_min_idx = peak1 + valley_min_idx
            else:
                absolute_min_idx = peak2 + (valley_min_idx - len(valley1))

            min_temperature = bin_edges[absolute_min_idx]
            return {'Site': Site, 'snow_class': snow_class, 'threshold': min_temperature}

        else:
            # Complex histogram — segmentation is the last resort
            fallback = try_segmentation_fallback()
            return {'Site': Site, 'snow_class': snow_class, 'threshold': fallback}


def analyze_and_segment_histogram(data, Site, year, snow_class, ThresholdDir, bin_range=(150, 300), plot=True):

    # ⛔ Skip processing for specified classes
    skip_classes = ["Ocean", "Maritime"]
    if snow_class in skip_classes:
        return {
            'Site': Site,
            'year': year,
            'snow_class': snow_class,
            'threshold_used': None,
            'threshold_intersection': None,
            'threshold_valley': None,
            'segments': None,
            'gaussian_fits': None,
            'valley_index': None,
            'peak_distance': None,
            'effective_end_peak1': None,
            'effective_start_peak2': None
        }

    # --- Guard: warn if data has values outside bin_range ---
    data_min = float(np.nanmin(data)) if len(data) > 0 else np.nan
    data_max = float(np.nanmax(data)) if len(data) > 0 else np.nan
    if np.isfinite(data_min) and (data_min < bin_range[0] or data_max > bin_range[1]):
        log.warning(
            "%s: data range [%.0f, %.0f] extends outside bin_range [%d, %d] — "
            "some data will be excluded from segmentation histogram.",
            Site, data_min, data_max, bin_range[0], bin_range[1]
        )

    # --- empty data guard ---
    if data is None or len(data) == 0:
        return {
            'Site': Site, 'year': year, 'snow_class': snow_class,
            'threshold_used': None, 'threshold_intersection': None,
            'threshold_valley': None, 'segments': None,
            'gaussian_fits': None, 'valley_index': None,
            'peak_distance': None, 'effective_end_peak1': None,
            'effective_start_peak2': None
        }

    # --- histogram & smoothing ---
    bins = np.arange(bin_range[0], bin_range[1] + 1)
    hist, bin_edges = np.histogram(data, bins=bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    hist_smooth = gaussian_filter1d(hist, sigma=2)

    # --- peak finding (80th percentile height, min separation 16 bins) ---
    height_thr = np.percentile(hist, 80)
    distance = 16
    peaks_idx, _ = find_peaks(hist_smooth, distance=distance, height=height_thr)
    peaks_idx = np.sort(peaks_idx)

    # --- RELAXED FALLBACK: if not exactly 2 peaks, try to find best 2 ---
    if len(peaks_idx) != 2:
        # Try progressively relaxed height thresholds
        for percentile in [70, 60, 50]:
            height_thr_relaxed = np.percentile(hist, percentile)
            peaks_relaxed, _ = find_peaks(hist_smooth, distance=distance, height=height_thr_relaxed)
            peaks_relaxed = np.sort(peaks_relaxed)
            if len(peaks_relaxed) >= 2:
                # Pick the 2 tallest among candidates
                top2 = peaks_relaxed[np.argsort(hist_smooth[peaks_relaxed])[-2:]]
                peaks_idx = np.sort(top2)
                log.warning(
                    "%s: used relaxed peak detection (percentile=%d) — found %d peaks, selected 2.",
                    Site, percentile, len(peaks_relaxed)
                )
                break

    # If still not 2 peaks after relaxation, try absolute minimum: 
    # just find the 2 highest points in the smoothed histogram
    if len(peaks_idx) != 2:
        all_peaks, _ = find_peaks(hist_smooth, distance=distance)
        if len(all_peaks) >= 2:
            top2 = all_peaks[np.argsort(hist_smooth[all_peaks])[-2:]]
            peaks_idx = np.sort(top2)
            log.warning("%s: used top-2 peaks as last resort in segmentation.", Site)
        else:
            # Cannot find 2 peaks even after all relaxation attempts.
            # Do NOT manufacture a threshold from a flat/noisy histogram —
            # returning None is more honest than a spurious value.
            log.warning(
                "%s: could not find 2 peaks even with relaxation — "
                "returning None (no reliable threshold).", Site
            )
            return {
                'Site': Site, 'year': year, 'snow_class': snow_class,
                'threshold_used': None,
                'threshold_intersection': None,
                'threshold_valley': None,
                'segments': None,
                'gaussian_fits': None,
                'valley_index': None,
                'peak_distance': None,
                'effective_end_peak1': None,
                'effective_start_peak2': None
            }

    # --- valley (histogram minimum) between the 2 selected peaks ---
    peak1, peak2 = peaks_idx
    valley_region = hist_smooth[peak1:peak2]
    valley_min_idx = np.argmin(valley_region)
    valley_idx = peak1 + valley_min_idx
    boundary_valley = float(bin_centers[valley_idx])

    # --- initial segmentation by valley for Gaussian fits ---
    seg_peak1_init = data[data < boundary_valley]
    seg_peak2_init = data[data >= boundary_valley]

    def fit_gaussian(segment):
        return norm.fit(segment) if len(segment) > 0 else (np.nan, np.nan)

    mu1_init, sigma1_init = fit_gaussian(seg_peak1_init)
    mu2_init, sigma2_init = fit_gaussian(seg_peak2_init)

    # --- Gaussian intersection ---
    w1 = max(len(seg_peak1_init), 1)
    w2 = max(len(seg_peak2_init), 1)

    def gaussian_intersection(mu1, s1, w1, mu2, s2, w2):
        if not np.isfinite(s1) or not np.isfinite(s2) or s1 <= 0 or s2 <= 0:
            return np.nan
        A = 0.5*(1.0/s2**2 - 1.0/s1**2)
        B = -(mu2/(s2**2)) + (mu1/(s1**2))
        C = (mu2**2)/(2*s2**2) - (mu1**2)/(2*s1**2) - np.log((w2*s1)/(w1*s2))
        if abs(A) < 1e-12:
            if abs(B) < 1e-12:
                return 0.5*(mu1 + mu2)
            return -C / B
        disc = B*B - 4*A*C
        if disc < 0:
            return np.nan
        r1 = (-B + np.sqrt(disc)) / (2*A)
        r2 = (-B - np.sqrt(disc)) / (2*A)
        lo, hi = (mu1, mu2) if mu1 < mu2 else (mu2, mu1)
        candidates = [r for r in (r1, r2) if lo <= r <= hi]
        if candidates:
            mid = 0.5*(mu1 + mu2)
            return min(candidates, key=lambda r: abs(r - mid))
        mid = 0.5*(mu1 + mu2)
        return r1 if abs(r1 - mid) < abs(r2 - mid) else r2

    boundary_intersection = gaussian_intersection(mu1_init, sigma1_init, w1,
                                                  mu2_init, sigma2_init, w2)

    def is_valid_intersection(x, m1, m2):
        if not np.isfinite(x):
            return False
        lo, hi = (m1, m2) if m1 < m2 else (m2, m1)
        return (lo <= x <= hi)

    use_intersection = is_valid_intersection(boundary_intersection, mu1_init, mu2_init)
    boundary_final = float(boundary_intersection) if use_intersection else float(boundary_valley)
    which = "intersection" if use_intersection else "valley"

    # --- final segmentation and Gaussian fits ---
    seg_peak1 = data[data < boundary_final]
    seg_peak2 = data[data >= boundary_final]
    mu1, sigma1 = fit_gaussian(seg_peak1)
    mu2, sigma2 = fit_gaussian(seg_peak2)
    gaussian_fits = {'peak1': (mu1, sigma1), 'peak2': (mu2, sigma2)}

    effective_end_peak1 = mu1 + 1.7 * sigma1 if np.isfinite(mu1) and np.isfinite(sigma1) else np.nan
    effective_start_peak2 = mu2 - 2.3 * sigma2 if np.isfinite(mu2) and np.isfinite(sigma2) else np.nan
    peak_distance = (effective_start_peak2 - effective_end_peak1
                     if np.isfinite(effective_start_peak2) and np.isfinite(effective_end_peak1)
                     else np.nan)

    seg_rest = (data[(data >= effective_end_peak1) & (data <= effective_start_peak2)]
                if np.all(np.isfinite([effective_end_peak1, effective_start_peak2]))
                else data[[]])
    segments = {
        'peak1': pd.DataFrame({'TB': seg_peak1}),
        'peak2': pd.DataFrame({'TB': seg_peak2}),
        'rest':  pd.DataFrame({'TB': seg_rest})
    }

    # --- plotting ---
    title_year = f"{year[0]}" if len(year) == 1 else f"{year[0]}–{year[-1]}"
    save_filename = f"{Site}_{title_year}_ThresholdRanges.png"
    save_path = os.path.join(ThresholdDir, save_filename)

    if plot:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.hist(data, bins=bins, alpha=0.3, color='gray', label="Data")
        x_vals = np.linspace(bin_range[0], bin_range[1], 1000)
        bin_width = bins[1] - bins[0]
        if np.isfinite(mu1) and np.isfinite(sigma1) and sigma1 > 0:
            pdf1 = norm.pdf(x_vals, mu1, sigma1) * len(seg_peak1) * bin_width
            ax.plot(x_vals, pdf1, color='blue', label=f"Peak 1 (μ={mu1:.0f})")
        if np.isfinite(mu2) and np.isfinite(sigma2) and sigma2 > 0:
            pdf2 = norm.pdf(x_vals, mu2, sigma2) * len(seg_peak2) * bin_width
            ax.plot(x_vals, pdf2, color='green', label=f"Peak 2 (μ={mu2:.0f})")
        if np.isfinite(effective_end_peak1):
            ax.axvline(x=effective_end_peak1, color='black', linestyle='--',
                       label=f"End Peak 1 (+1.7σ): {effective_end_peak1:.0f}K")
        if np.isfinite(effective_start_peak2):
            ax.axvline(x=effective_start_peak2, color='red', linestyle='--',
                       label=f"Start Peak 2 (-2.3σ): {effective_start_peak2:.0f}K")
        ax.axvline(x=boundary_final, color='purple', linestyle=':',
                   label=f"Threshold used ({which}): {boundary_final:.0f}K")
        if not np.isnan(boundary_intersection) and abs(boundary_intersection - boundary_final) > 1e-6:
            ax.axvline(x=boundary_intersection, color='orange', linestyle='-.',
                       label=f"Intersection: {boundary_intersection:.0f}K")
        ax.set_title(f"{Site} ({title_year}) Histogram & Gaussian Fits")
        ax.set_xlabel("Brightness Temperature (K)")
        ax.set_ylabel("Frequency")
        ax.legend()
        plt.tight_layout()
        plt.close(fig)

    return {
        'Site': Site,
        'year': year,
        'snow_class': snow_class,
        'threshold_used': boundary_final,
        'threshold_intersection': float(boundary_intersection) if np.isfinite(boundary_intersection) else None,
        'threshold_valley': float(boundary_valley),
        'segments': segments,
        'gaussian_fits': gaussian_fits,
        'valley_index': int(valley_idx),
        'peak_distance': float(peak_distance) if np.isfinite(peak_distance) else None,
        'effective_end_peak1': int(round(effective_end_peak1)) if np.isfinite(effective_end_peak1) else None,
        'effective_start_peak2': int(round(effective_start_peak2)) if np.isfinite(effective_start_peak2) else None
    }



def plot_histogram(hist, bin_edges, hist_smooth, min_temperature, Site, year, snow_class, sensor_SIR, channel_SIR, ThresholdDir, peaks_idx=None):
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

    ### Note: If the site's name is saved as a descriptive name (e.g., Monument Creek), it indicates that the output is from a single-site analysis. However, if the site is labeled using row and column coordinates (e.g., 344,248), it refers to a pixel in a region-based analysis.
    if len(year) == 1:
        if 'processed_GRD_pixels' in locals() and processed_GRD_pixels:
            save_filename = f"region_{grd_pixel_key[0]}_{grd_pixel_key[1]}_{year[0]}_{sensor_SIR}.png"
        else:
            save_filename = f"{Site}_{year[0]}_{sensor_SIR}.png"
    else:
        if 'processed_GRD_pixels' in locals() and processed_GRD_pixels:
            save_filename = f"region_{grd_pixel_key[0]}_{grd_pixel_key[1]}_{year[0]}–{year[-1]}_{sensor_SIR}.png"
        else:
            save_filename = f"{Site}_{year[0]}–{year[-1]}_{sensor_SIR}.png"

        
    save_path = os.path.join(ThresholdDir, save_filename)

    # Create and configure plot
    fig, ax = plt.subplots()
    if len(year) == 1:
        ax.set_title('All data SIR Histogram (' + str(year[0]) + ') ' + Site + ' (' + snow_class + ')' + ' for ' + sensor_SIR + ' ' + channel_SIR)
    else:
        ax.set_title('All data SIR Histogram (' + str(year[0]) + '-' + str(year[-1]) + ') ' + Site + ' (' + snow_class + ') for ' + sensor_SIR + ' ' + channel_SIR)

    ax.set_xlabel('Brightness Temp (K)')
    ax.set_ylabel('Frequency')
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
    #plt.savefig(save_path)
    #plt.show()
    plt.close(fig)
