### CETB_analysis.py
### last edited: 5/20/18
### by: Mitch Johnson
### functions for analysis of CETB data, creates histograms and time series plots of Tb and DAV
## 

from netCDF4 import Dataset, num2date
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pdb; # insert at places for breakpoints: pdb.set_trace()
import seaborn as sns
import warnings

# getting a runtimewarning when using operators on numpy arrays with lots of NaNs, functions still perform, but using this command to suppress the warning
warnings.filterwarnings("ignore",category =RuntimeWarning)


# this function accepts the time-series Tb data created by the read_Tb function in CETB_read_functions.py and will save a figure displaying an annual histogram of Tb for the specified area
# the function also accepts the prefix (beginning of file path) and lat_start and lon_start for titling and naming the figure
def Tb_hist_annual(CETB_data, prefix, lat_start, lon_start):
    y = CETB_data[:,:,:] # CETB_data for all pixels in subset
    y = y[y>=0]
    bins = range(150,300)
    fig,ax=plt.subplots()
    ax.hist(y, bins)
    ax.set_title(prefix)
    ax.set_xlabel('Brightness Temp (K)')
    return

# plots 12 histogram subplots of all the Tb data that was loaded, one for each month
def Tb_hist_monthly(CETB_data, prefix, lat_start, lon_start, cal_month):
    #plot histogram of Tb measurements by month, accepts all pixels passed
    months=[0,1,2,3,4,5,6,7,8,9,10,11]
    x=[None]*12
    n=1
    names=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    fig, axs = plt.subplots(3, 4, sharey=True, facecolor='#2c7bb6')
    n_bins = range(150,350)
    axs = axs.ravel()
    for month in months:
        x[month]=CETB_data[cal_month==n]
        x[month]=x[month][x[month]>=0]
        n=n+1
        axs[month].hist(x[month],bins=n_bins)
        axs[month].set_title(names[month])
    fig.autofmt_xdate()
    fig.tight_layout()
    return

# this function accepts the time-series DAV data created by the calc_DAV function in CETB_read_functions.py and will save a figure displaying an annual histogram of DAV for the specified area
# the function also accepts the prefix (beginning of file path) and lat_start and lon_start for titling and naming the figure
def DAV_hist_annual(DAV, prefix, lat_start, lon_start):
    y = DAV[:,:,:] # DAV for all pixels in subset
    y = y[y>=0]
    bins = range(0,60)
    fig,ax=plt.subplots()
    ax.hist(y, bins, color='#2c7bb6')
    ax.set_title(prefix)
    ax.set_xlabel('DAV (K)')
    #plt.savefig('/home/mij216/Documents/WesternUS/Figures/histograms/DAV_'+prefix+'_'+str(lat_start)+'_'+str(lon_start)+'.png')
    return

# plot histogram of DAV measurements by month, accepts all pixels passed, HAS NOT BEEN CHECKED IN NB YET, SAME CODE AS TB MONTHLY
def DAV_hist_monthly(DAV, prefix, lat_start, lon_start, cal_month):
    months=[0,1,2,3,4,5,6,7,8,9,10,11]  # initialize list of months to loop through
    x=[None]*12  #initialize a list
    n=1
    names=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    fig, axs = plt.subplots(3, 4, sharey=True, facecolor='#2c7bb6')
    n_bins = range(0,60)
    axs = axs.ravel()
    for month in months:
        x[month]=DAV[cal_month==n]
        x[month]=x[month][x[month]>=0]
        n=n+1
        axs[month].hist(x[month],bins=n_bins)
        axs[month].set_title(names[month])
    fig.autofmt_xdate()
    fig.tight_layout()
    return

def TbDAV_series_one_year(CETB_data, DAV, cal_date, cal_year, year, Tb_threshold, DAV_threshold):
    #plots Tb and DAV for all pixels that were read in, enter year of interest
    fig, (ax1, ax2) = plt.subplots(2,1, sharex=True)  #create two subplots with a shared x-axis
    y_dims_list=list(range(len(CETB_data[0,:,0])))  # creates a list of the y-dimension pixel indices, used for plotting
    x_dims_list=list(range(len(CETB_data[0,0,:])))    # creates a list of the x-dimension pixel indices, for plotting
    
    thisYr = cal_year==year
    
    for i in y_dims_list:  #plot the time-series
        for j in x_dims_list:
            CETB_df = pd.DataFrame(CETB_data[thisYr, i, j], index=cal_date[thisYr])
            CETB_df.plot(ax=ax1, legend=False)
            DAV_df = pd.DataFrame(DAV[thisYr, i, j], index=cal_date[thisYr])
            DAV_df.plot(ax=ax2, legend=False)

    ax1.set_ylabel('Tb (K)')
    ax2.set_ylabel('DAV (K)')
    ax1.axhline(y=Tb_threshold,c="red",linewidth=0.5,zorder=0)  # horizontal line at Tb=252K, proposed melt threshold
    ax2.axhline(y=DAV_threshold,c="red",linewidth=0.5,zorder=0)   # horizontal line at DAV=18K, proposed DAV melt threshold
    fig.autofmt_xdate()  #rotates the date labels so they aren't bunched up
    return

# calculate early season melt events, accepts one pixel or a rectangular set of adjacent pixels
def early_melt_events(CETB_data, DAV, DAV_threshold, Tb_threshold, cal_date, Years, rows_cols):
    # count the number of DAV/Tb threshold exceedances in Jan/Feb of each year, returns a dataframe where each row is a pixel and each column is a year

    y_dims_list=list(range(len(CETB_data[0,:,0])))  # creates a list of the y-dimension pixel indices, used for plotting
    x_dims_list=list(range(len(CETB_data[0,0,:])))    # creates a list of the x-dimension pixel indices, for plotting
    y_s=list(range(rows_cols[0],rows_cols[1]))  # makes a list of the y(row) numbers so the dataframe of early melt events can be indexed
    x_s=list(range(rows_cols[2],rows_cols[3]))  # makes a list of the x(col) numbers so the dataframe of early melt events can be indexed
    melt_condition_met = (DAV[:]>DAV_threshold) & (CETB_data[:]>Tb_threshold)

    flag = melt_condition_met.astype(int)
    matrix=pd.DataFrame()
    for i in y_dims_list:
        for j in x_dims_list:
            column=pd.DataFrame(data=flag[:,i,j], columns=[str(y_s[i])+','+str(x_s[j])])
            matrix=pd.concat([matrix,column],axis=1)
    matrix=matrix.set_index(cal_date)
    events = pd.DataFrame()
    for year in Years:
            events[year-Years[0]]=matrix[str(year)+'-01-01':str(year)+'-03-31'].sum()
    events.columns=[Years]
    events[events==0]=np.NaN
    events=events.dropna(axis=0, how='all')
    events=events.dropna(axis=1, how='all')

    return events

# plot the early melt events on a heatmap, not georeferenced
def earlymelt_map(CETB_data, DAV, DAV_threshold, Tb_threshold, cal_date, Years, rows_cols):
    y_dims_list=list(range(len(CETB_data[0,:,0])))  # creates a list of the y-dimension pixel indices, used for plotting
    x_dims_list=list(range(len(CETB_data[0,0,:])))    # creates a list of the x-dimension pixel indices, for plotting
    y_s=list(range(rows_cols[0],rows_cols[1]))  # makes a list of the y(row) numbers so the dataframe of early melt events can be indexed
    x_s=list(range(rows_cols[2],rows_cols[3]))  # makes a list of the x(col) numbers so the dataframe of early melt events can be indexed
    melt_condition_met = (DAV[:]>DAV_threshold) & (CETB_data[:]>Tb_threshold)
    flag = melt_condition_met.astype(int)
    #flag.sum(axis=1)
    matrix=pd.DataFrame()
    for i in y_dims_list:
        for j in x_dims_list:
            column=pd.DataFrame(data=flag[:,i,j], columns=[str(y_s[i])+','+str(x_s[j])])
            matrix=pd.concat([matrix,column],axis=1)
    matrix=matrix.set_index(cal_date)
    events = pd.DataFrame()
    for year in Years:
        events[year-Years[0]]=matrix[str(year)+'-01-01':str(year)+'-03-31'].sum()
    events.columns=[Years]
    events[events==0]=np.NaN
    sum_events=events.sum(axis=1)
    sum_events=sum_events.as_matrix()
    sum_events=sum_events.reshape(8,8)

    # plot early melt events on heatmap
    fig,ax=plt.subplots()
    ax=sns.heatmap(sum_events, cmap='Blues')
    return

# this function only works when 2 CETB datasets are passed (2channel/alg notebook), specifically designed to take a subset of 3 km pixels and one 25km GRD pixel
# creates three arrays for the min, max, and avg at each time step in the 64 pixel subset, also takes an array for the 25 km pixel
# plots the four arrays on a histogram and the four arrays on a time-series, function takes a 'year' argument for plotting a desired year on the time-series 
def min_max_series(CETB_GRD, CETB_SIR, cal_date, cal_year, year, title):
    ymins=np.nanmin(CETB_SIR, axis=1)
    xmins=np.nanmin(CETB_SIR, axis=2)
    xmins2=np.nanmin(xmins, axis=1)
    ymins2=np.nanmin(ymins, axis=1)
    totalmin=np.fmin(xmins2,ymins2)

    ymax=np.nanmax(CETB_SIR, axis=1)
    xmax=np.nanmax(CETB_SIR, axis=2)
    xmax2=np.nanmax(xmax, axis=1)
    ymax2=np.nanmax(ymax, axis=1)
    totalmax=np.fmax(xmax2,ymax2)

    ymean=np.nanmean(CETB_SIR, axis=1)
    xmean=np.nanmean(CETB_SIR, axis=2)
    xmean2=np.nanmean(xmean, axis=1)
    ymean2=np.nanmean(ymean, axis=1)
    totalmean=np.nanmean([xmean2,ymean2],axis=0)

    y_GRD=np.nanmean(CETB_GRD, axis=1)
    x_GRD=np.nanmean(CETB_GRD, axis=2)
    x_GRD_2=np.nanmean(x_GRD, axis=1)
    y_GRD_2=np.nanmean(y_GRD, axis=1)
    total_GRD=np.nanmean([x_GRD_2,y_GRD_2],axis=0)

    # hist of min/max/mean of the 3 km pixels in subset and the 25km (GRD) pixel that envelopes them
    frame=pd.DataFrame(data={'min':totalmin, 'max':totalmax, 'avg':totalmean, 'GRD':total_GRD}, index=cal_date)

    return frame


# Helper function to parse "row,col" strings
def parse_row_col(s):
    return [int(str) for str in s.split(',')]


# plot map of average MOD for period of record
# Inputs:
#   datadir: top directory location with cube data
#   prefix:
#   CETB: dictionary with data read using read_Tb_whole
#         dict fields are expected for TB, latitude, longitude, gpd
#   DAV: DAV as calculated by calc_DAV, dimensions must match data['TB']
#   rows_cols: 4-element tuple with (row_begin, row_end, col_begin, col_end)
#         of the subset to read from the cube
#   Years: list of years to process
#   window: window for algorithm, 10 would be 5 days (assumes 2 measurements/day)
#   count: number of TB/DAV exceedances to trigger melt onset
#   EHD_window: window for end of high DAV, similar to window
#   EHD_count: number of high TB/low DAV occurrence after melt onset
#   DAV_threshold:
#   Tb_threshold:
def MOD_array(datadir, prefix, CETB, DAV, rows_cols, Years, window, count,
              EHD_window, EHD_count, DAV_threshold, Tb_threshold):

    # Find times/places when melt conditions are satisfied
    melt_condition_met = (
        DAV > DAV_threshold) & (CETB['TB'][:, :, :] > Tb_threshold)
    flag = melt_condition_met.astype(int)
    
    # Find times/places when end of high DAV conditions are satisfied
    EHD_condition_met = (
        DAV <= DAV_threshold) & (CETB['TB'][:, :, :] > Tb_threshold)
    EHD_flag = EHD_condition_met.astype(int)
    
    # Prepare a DataFrame to do the heavy-lifting on the algorithm:
    # Define a list of column_names with one for each array row, col
    # Populate each row with the flattened data for that date
    col_names = ["%s,%s" % (str(y),str(x)) 
                 for y in np.arange(rows_cols[0], rows_cols[1])
                 for x in np.arange(rows_cols[2], rows_cols[3])]

    # newdata: rows are dates, columns are pixels
    newdata = np.zeros([flag.shape[0],
                        flag.shape[1] * flag.shape[2]], dtype=flag.dtype)
    
    print("newdata.shape %s" % str(newdata.shape))
    print("moving flag array to newdata...")
    print("number of days = %d" % flag.shape[0])
    
    # EHD_newdata: rows are dates, columns are pixels
    EHD_newdata = np.zeros([EHD_flag.shape[0],
                        EHD_flag.shape[1] * EHD_flag.shape[2]], dtype=EHD_flag.dtype)
    
    print("EHD_newdata.shape %s" % str(EHD_newdata.shape))
    print("moving EHD_flag array to EHD_newdata...")
    print("number of days = %d" % EHD_flag.shape[0])

    for d in np.arange(flag.shape[0]):
        if np.mod(d, 100) == 0:
            print("Next d = %d" % d)
        newdata[d,] = flag[d, :, :].flatten()
        EHD_newdata[d,] = EHD_flag[d, :, :].flatten()
    
    # Converting numpy arrays to dataframe
    matrix = pd.DataFrame(data=newdata, columns=col_names)
    matrix.set_index(pd.Index(CETB['cal_date']), inplace=True)

    meltflag_df=matrix.copy(deep=True)
    
    #set breakpoint
    #pdb.set_trace()
    
    print("dataFrame is ready with flag data")
    print("doing rolling sums...")
    # MOD algorithm - calculate sum on a rolling window
    # (window= no. of obs, 2 per day)
    matrix=matrix.rolling(window).sum()
    
    # count= no. times thresholds are tripped for algo
    # returns all values in the matrix that are >= count and
    # meet the melt criteria (i.e. melt thresholds)
    matrix=matrix[(matrix>=count) & (meltflag_df==1)]
    # drop rows (dates) with all pixels NaN
    matrix=matrix.dropna(axis=0, how='all') 

    # gets a dataframe with one row for each pixel,
    # one column for each year,
    # MOD in DOY in each cell for that pixel and year
    df = pd.DataFrame()
    num_pixels = len(matrix.columns)
    
    #
    EHD_matrix = pd.DataFrame(data=EHD_newdata, columns=col_names)
    EHD_matrix.set_index(pd.Index(CETB['cal_date']), inplace=True)

    EHDflag_df=EHD_matrix.copy(deep=True)
    
    print("dataFrame is ready with EHD flag data")
    print("doing rolling sums...")
    # EHD algorithm - calculate sum on a rolling window
    # (window= no. of obs, 2 per day)
    EHD_matrix=EHD_matrix.rolling(EHD_window).sum()
    
    # count= no. times thresholds are tripped for algo
    EHD_matrix=EHD_matrix[EHD_matrix>=EHD_count]
    # drop rows (dates) with all pixels NaN
    EHD_matrix=EHD_matrix.dropna(axis=0, how='all') 

    # gets a dataframe with one row for each pixel,
    # one column for each year,
    # EHD in DOY in each cell for that pixel and year
    EHD_df = pd.DataFrame()
    EHD_num_pixels = len(EHD_matrix.columns)
    
    # Find the first value date for each year in each column
    # It's possible that no melt condition is met for a given year/column
    # but this shows up in different ways:
    # When matrix contains data for some years, but no data for a given year,
    # first_valid_index returns a KeyError
    # When matrix[year][column] is empty, first_valid_index returns None
    for year in Years:
        print("Next year = %d..." % year)
        dates = np.full(num_pixels, pd.NaT)
        for column_index, column in enumerate(matrix.columns):
            try:
                first_date = matrix.loc[str(year)][column].first_valid_index()
            except KeyError:
                print("MOD_array: no melt found for pixel %s in year %d" % (
                    column, year))
                continue

            if first_date is not None:
                dates[column_index] = first_date
            else:
                print("MOD_array: no melt found for pixel %s in year %d" % (
                    column, year))
                
        dates_series = pd.Series(dates)
        dates_series = dates_series.dt.dayofyear
        df = pd.concat([df, dates_series], axis=1)
    
    df.columns = Years
    
    # Above section copied and modified for EHD (IN PROGRESS) - need to
    # feed in df from line 331
    # If there is a melt onset date for this year and column, then look
    # for the ensuing EHD
    for year in Years:
        print("Next year = %d..." % year)
        dates = np.full(EHD_num_pixels, pd.NaT)
        for column_index, column in enumerate(EHD_matrix.columns):
            try:
                first_EHD_date = EHD_matrix.loc[str(year)][column].first_valid_index()
            except KeyError:
                print("MOD_array: no EHD found for pixel %s in year %d" % (
                    column, year))
                continue

            if first_EHD_date is not None:
                dates[column_index] = first_EHD_date
            else:
                print("MOD_array: no EHD found for pixel %s in year %d" % (
                    column, year))
                
        dates_series = pd.Series(dates)
        dates_series = dates_series.dt.dayofyear
        EHD_df = pd.concat([EHD_df, dates_series], axis=1)
    
    EHD_df.columns = Years

    # get the average MOD for each pixel and make it an array for plotting
    print("Getting average MOD at each pixel...")
    MOD = df.mean(axis=1).values
    MOD = np.ma.array(MOD)   # make it masked array
    MOD[MOD < 0] = np.ma.masked   #convert any invalid MODs to masked
    
    # Store the Avg MOD for these years as the last column in the data frame
    df['Avg'] = MOD
    data_columns = df.columns

    print("Setting geolocation information...")
    df['pixel'] = matrix.columns
    
    # get the average EHD for each pixel and make it an array for plotting
    print("Getting average EHD at each pixel...")
    EHD = EHD_df.mean(axis=1).values
    EHD = np.ma.array(EHD)   # make it masked array
    EHD[EHD < 0] = np.ma.masked   #convert any invalid EHDs to masked
    
    # Store the Avg EHD for these years as the last column in the data frame
    EHD_df['Avg'] = EHD
    EHD_data_columns = EHD_df.columns

    print("Setting geolocation information...")
    EHD_df['pixel'] = EHD_matrix.columns

    # Insert columns with subset pixel geolocations x, y, lat, lon, row, col
    # IN PROGRESS (may do EHD and MOD together)
    num_rows = len(df.index)
    xx, yy = np.meshgrid(CETB['x'], CETB['y'])
    df['x'] = np.reshape(xx, num_rows)
    df['y'] = np.reshape(yy, num_rows)
    df['latitude'] = np.reshape(CETB['latitude'], num_rows)
    df['longitude'] = np.reshape(CETB['longitude'], num_rows)
    df["row"] = df["column"] = ""
    df[["row", "column"]] = list(df.pixel.apply(parse_row_col))

    EHD_df['x'] = df['x']
    EHD_df['y'] = df['y']
    EHD_df['latitude'] = df['latitude']
    EHD_df['longitude'] = df['longitude']
    EHD_df["row"] = df["row"]
    EHD_df[["row", "column"]] = df[["row", "column"]]

    # Put the geolocation fields at the beginning of the columns
    geo_columns = df.columns[-7:]
    df = df[geo_columns.append(data_columns)]
    
    geo_columns = EHD_df.columns[-7:]
    EHD_df = EHD_df[geo_columns.append(EHD_data_columns)]
    
    return MOD, df, meltflag_df, EHD, EHD_df, EHDflag_df

# plot map of average MOD for year of interest
#def MOD_array_year(datadir, prefix, CETB_data, DAV,
#                   rows_cols, cal_date, year, window, count,
#                   DAV_threshold, Tb_threshold)
# THIS FUNCTION IS OBSOLETE, do the call to MOD_array and
# and then get the year that you are interested in from the output
