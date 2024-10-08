from datetime import datetime
from netCDF4 import Dataset, num2date
import numpy as np
import pandas as pd
from cetbtools.ease2conv import Ease2Transform
import warnings
import glob
import pdb; # insert at places for breakpoints: pdb.set_trace()

# getting a runtimewarning when using operators
# on numpy arrays with lots of NaNs, functions still perform,
# but using this command to suppress the warning
#FIXME: this shouldn't be necessary
warnings.filterwarnings("ignore",category =RuntimeWarning)

# Original version of load CETB data for a subset and sensor of interest
# We should eventually delete this in favor of the new version (below)
# that doesn't force entry of the subset to read and that
# returns gridding metadata from the .nc file
def read_Tb(datadir, prefix, Years,y_start,y_end,x_start,x_end):
    for year in Years:

        # Create filename
        filename=datadir+prefix+'.'+str(year)+'.TB.nc'
        list=glob.glob(filename)
        
        # load the raw data in
        rawdata = Dataset(list[-1], "r", format="NETCDF4")
                
        # Compile the CETB data, the TB variable is saved as (time, y, x) - 
        subset = rawdata.variables['TB'][0:,y_start:y_end,x_start:x_end]
        if year==Years[0]:
            CETB_data = subset
        else:
            CETB_data = np.concatenate((CETB_data, subset), axis=0)

        # Compile the date information
        d=rawdata.variables['time']
        date=d[:]
        greg_date = num2date(date[:], units=d.units,calendar=d.calendar)
        if year==Years[0]:
            cal_date = greg_date
        else:
            cal_date = np.concatenate((cal_date, greg_date), axis=0)

        # TB value of 60000 is (unscaled) missing data:
        # This value is used when there are TBs but the the
        # reconstruction algorithm didn't converge on an answer
        # TB value of 0 is no data value:
        # This value is used when there were no TB measurements
        # available to grid to this cell location
        #FIXME: these values should not be hardcoded, they should
        #be read from the variable metadata
        CETB_data[CETB_data==60000] = np.NaN
        CETB_data[CETB_data==0] = np.NaN

    # get date info for plotting
    cal_year = np.empty([len(cal_date)], dtype=float)
    for n in range(0,len(cal_date)):
        cal_year[n]=cal_date[n].year
            
    # set up an array with the month data, if want to examine Tb data for a particular month
    cal_month = np.empty([len(cal_date)], dtype=float)
    for n in range(0,len(cal_date)):
        cal_month[n]=cal_date[n].month

    return CETB_data, cal_date, cal_year, cal_month

# load CETB data for a sensor of interest
# if no subset arguments are included, reads the whole cube
# this is a more versatile version of read_Tb and should eventually
# replace it
# prefix can be a filename glob pattern
# Returns:
#   dict with fields for subsetted TB and geolocation information
#   read from the cube
def read_Tb_whole(datadir, prefix, Years,
                  y_start=None, y_end=None, x_start=None, x_end=None):

    # Read the whole cube if no subset is specified
    get_full_cube = False
    if not y_start and not y_end and not x_start and not x_end:
        get_full_cube = True
        print("No subset specified, fetching complete cube...")
        
    for year in Years:
        # Create filename
        cubePattern = "%s/%s.%4d.TB.nc" % (datadir, prefix, year)
        list = glob.glob(cubePattern)
        if not list:
            raise IOError("No cubefile found to match the pattern = %s" % (cubePattern))
        else:
            filename = list[-1]
            print("Next filename=%s..." % filename)

        # load the raw data in
        rawdata = Dataset(filename, "r", format="NETCDF4")

        # Compile the CETB data, the TB variable is saved as (time, y, x) -
        if get_full_cube:
            subset = rawdata.variables['TB'][:,:,:]
        else:
            subset = rawdata.variables['TB'][0:,y_start:y_end,x_start:x_end]
            
        if year==Years[0]:
            CETB_data = subset
        else:
            CETB_data = np.concatenate((CETB_data, subset), axis=0)

        # Compile the date information
        d=rawdata.variables['time']
        date=d[:]
        greg_date = num2date(date[:], units=d.units,calendar=d.calendar)
        if year==Years[0]:
            cal_date = greg_date
        else:
            cal_date = np.concatenate((cal_date, greg_date), axis=0)

        # TB value of 60000 is (unscaled) missing data:
        # This value is used when there are TBs but the the
        # reconstruction algorithm didn't converge on an answer
        # TB value of 0 is no data value:
        # This value is used when there were no TB measurements
        # available to grid to this cell location
        #FIXME: these values should not be hardcoded, they should
        #be read from the variable metadata
        CETB_data[CETB_data==60000] = np.NaN
        CETB_data[CETB_data==0] = np.NaN

        # First time through,
        # Fetch the gpd, (x, y) and (lat, lon) coordinates
        # of the subset of the cube
        if year == Years[0]:

            gpd = rawdata.variables['crs'].long_name
            if get_full_cube:
                x = rawdata.variables['x'][:]
                y = rawdata.variables['y'][:]
                latitude = rawdata.variables['latitude'][:]
                longitude = rawdata.variables['longitude'][:]
            else:
                x = rawdata.variables['x'][x_start:x_end]
                y = rawdata.variables['y'][y_start:y_end]
                latitude = rawdata.variables['latitude'][
                    y_start:y_end,x_start:x_end]
                longitude = rawdata.variables['longitude'][
                    y_start:y_end,x_start:x_end]

        # Close the file
        rawdata.close()

    # get date info for plotting
    cal_year = np.empty([len(cal_date)], dtype=float)
    for n in range(0,len(cal_date)):
        cal_year[n]=cal_date[n].year

    # set up an array with the month data,
    # if want to examine Tb data for a particular month
    cal_month = np.empty([len(cal_date)], dtype=float)
    for n in range(0,len(cal_date)):
        cal_month[n]=cal_date[n].month
        
    # Convert cal_dates masked array to pd.datetime format
    cal_date = np.array([
        pd.to_datetime(i.strftime("%m/%d/%Y, %H:%M:%S")) for i in cal_date.data ])

    return {'TB': CETB_data,
            'cal_date': cal_date,
            'cal_year': cal_year,
            'cal_month': cal_month,
            'gpd': gpd,
            'x': x,
            'y': y,
            'latitude': latitude,
            'longitude': longitude}


# load ALL TB data for a cubefile
def read_Tb_all(datadir, prefix, Years):
    for year in Years:
        # Create filename
        filename=datadir+prefix+'.'+str(year)+'.TB.nc'
        list=glob.glob(filename)
        # load the raw data in
        rawdata = Dataset(list[-1], "r", format="NETCDF4")
                
        # Compile the CETB data, the TB variable is saved as (time, y, x) - 
        subset = rawdata.variables['TB'][0:,:,:]
        if year==Years[0]:
            CETB_data = subset
        else:
            CETB_data = np.concatenate((CETB_data, subset), axis=0)
                
        # Compile the date information
        d=rawdata.variables['time']
        date=d[:]
        greg_date = num2date(date[:], units=d.units,calendar=d.calendar) 
        if year==Years[0]:
            cal_date = greg_date
        else:
            cal_date = np.concatenate((cal_date, greg_date), axis=0)

        # TB value of 60000 is (unscaled) missing data:
        # This value is used when there are TBs but the the
        # reconstruction algorithm didn't converge on an answer
        # TB value of 0 is no data value:
        # This value is used when there were no TB measurements
        # available to grid to this cell location
        #FIXME: these values should not be hardcoded, they should
        #be read from the variable metadata
        CETB_data[CETB_data==60000] = np.NaN
        CETB_data[CETB_data==0] = np.NaN
        
    # get date info for plotting        
    cal_year = np.empty([len(cal_date)], dtype=float)
    for n in range(0,len(cal_date)):
        cal_year[n]=cal_date[n].year
    
    # set up an array with the month data, if want to examine Tb data for a particular month
    cal_month = np.empty([len(cal_date)], dtype=float)
    for n in range(0,len(cal_date)):
        cal_month[n]=cal_date[n].month

    return CETB_data, cal_date, cal_year, cal_month

# read in Tb standard deviation data
def read_Tb_std_dev(datadir, prefix, Years,y_start,y_end,x_start,x_end):
    for year in Years:
        # Create filename
        cubePattern = "%s/%s.%4d.TB_std_dev.nc" % (datadir, prefix, year)
        list = glob.glob(cubePattern)
        if not list:
            raise IOError("No cubefile found to match the pattern = %s" % (cubePattern))
        else:
            filename = list[-1]
            print("Next filename=%s..." % filename)
        
        # load the raw data in
        rawdata = Dataset(filename, "r", format="NETCDF4")
                
        # Compile the CETB data, the TB variable is saved as (time, y, x) - 
        subset = rawdata.variables['TB_std_dev'][0:,y_start:y_end,x_start:x_end]
        if year==Years[0]:
            CETB_data = subset
        else:
            CETB_data = np.concatenate((CETB_data, subset), axis=0)
                
        # stddev value of 2^16 - 2 missing data:
        # This value is used when there are TBs but the the
        # reconstruction algorithm didn't converge on an answer
        # stddev value of 2^16 - 1 is no data value:
        # This value is used when there were no TB measurements
        # available to grid to this cell location
        # This approach works because 60000 < these values
        #FIXME: these values should not be hardcoded, they should
        #be read from the variable metadata
        CETB_data[CETB_data>=60000] = np.NaN
        CETB_data[CETB_data==0] = np.NaN

    return CETB_data


# read in the Tb_time variable - time of the satellite overpass
def read_Tb_time(datadir, prefix, Years,y_start,y_end,x_start,x_end):
    for year in Years:

        # Create filename
        filename=datadir+prefix+'.'+str(year)+'.TB_time.nc'
        list=glob.glob(filename)
        # load the raw data in
        rawdata = Dataset(list[-1], "r", format="NETCDF4")
        
        # Compile the CETB data, the TB variable is saved as (time, y, x) - 
        subset = rawdata.variables['TB_time'][0:,y_start:y_end,x_start:x_end]
        if year==Years[0]:
            CETB_data = subset
        else:
            CETB_data = np.concatenate((CETB_data, subset), axis=0)
                
        # TB_time value has no missing data values
        # TB_time value of min(INT16) is no data value:
        # This value is used when there were no TB measurements
        # available to grid to this cell location
        #FIXME: these values should not be hardcoded, they should
        #be read from the variable metadata
        #CETB_data[CETB_data>60000] = np.NaN
        #CETB_data[CETB_data==0] = np.NaN

    return CETB_data

def coords(datadir, prefix, lat_start, lat_end, lon_start, lon_end):
    # this function takes the user inputs of latitude and longitude
    # and returns the row/col numbers that identify Tb grids in the
    # cubefile
    # NOTE: change made April 2021: returned x_end and y_end are
    # incremented by 1. This is different from prior versions,
    # see note below. This may now be the correct action, we
    # need to test it.    
    # NOTE: this is a function I wrote, but it does basically the
    # same thing as "find_cube_offset" and "grid_locationsof..."
    # which are listed below and were coded by Mary Jo
    # Just read any year since the structure is the same,
    filename=datadir+prefix+'.*.TB.nc'
    list=glob.glob(filename)
    data=Dataset(list[-1], "r", format="NETCDF4")

    lat=data.variables['latitude']
    lat=lat[:]
    lon=data.variables['longitude']
    lon=lon[:]
    lat_lon=np.dstack((lat,lon))

    # if latitude and longitude are the same, for reading one pixel
    if ((lat_start==lat_end) & (lon_start==lon_end)):
        y_start=np.min(np.argmin(
            ((lat_start-lat_lon[:,:,0])**2+(lon_start-lat_lon[:,:,1])**2),
                axis=0))
        y_end=y_start+1
        x_start=np.max(np.argmin(
            ((lat_start-lat_lon[:,:,0])**2+(lon_start-lat_lon[:,:,1])**2),
                axis=1))
        x_end=x_start+1
        return y_start,y_end,x_start,x_end
    else:
        row_col_set=np.where(
            (lat_lon[:,:,0]>lat_start)&(
                lat_lon[:,:,0]<lat_end)&(
                     lat_lon[:,:,1]>lon_start)&(
                        lat_lon[:,:,1]<lon_end))
        y_start=np.min(row_col_set[0])
        y_end=np.max(row_col_set[0])
        x_start=np.min(row_col_set[1])
        x_end=np.max(row_col_set[1])
        # if x_start==x_end:
        #     x_end=x_end+1
        # if y_start==y_end:
        #     y_end=y_end+1
        x_end=x_end+1
        y_end=y_end+1

        # The returned x_end and y_end are 1 pixel further than the
        # required area, because the subsequent subsetting is an
        # open interval the ending row/col
        # Probably this should be handled in the subsetters, but for
        # now I'm going to test it by extending the bounds by 1
        # row and column here...
        return y_start,y_end,x_start,x_end

def calc_DAV(CETB_data):

    # function takes the CETB_data that was read in read_Tb() and returns the absolute value of the DAV         
    DAV=np.diff(CETB_data,n=1,axis=0)
    DAV_abs=np.absolute(DAV)    
    DAV_abs=np.insert(DAV_abs, [0],[0], axis=0)         # insert a 0 at beginning of array so same length as CETB_data for plotting together

    return DAV_abs

# find the offset to be used for locating the grid row/col from lat/lon coordinates
# for the requested subsetName
# subsetNames are, for e.g. "GLaIL" or "WesternCA" etc.
def find_cube_offset(subsetName, cubeDir=None, cubeType=None, verbose=False):

    # when no cubeDir is specified, assume it's in a regular
    # location on fringe        
    if not cubeDir:
        cubeDir = "%s%s" % (
            "/mnt/data3/cetb/cubes/AQUA_AMSRE/",
            subsetName )

    if not cubeType:
        cubeType = '3*V-SIR'

    cubePattern = "%s/CETB.cubefile.%s.*-%s-*-v*.*.TB.nc" % (
        cubeDir,
        subsetName,
        cubeType)
    
    # Just use the last one found
    # Check for no files found and return a reasonable error message
    list = glob.glob(cubePattern)
    if not list:
        raise IOError("No cubefiles found to match the pattern = %s" % (cubePattern))
    else:
        cubeFile = list[-1]
        print("Reading offset information from cubeFile=%s..." % (cubeFile))
        
    f = Dataset(cubeFile, "r", "NETCDF4")   
    lats = f.variables['latitude'][:]
    lons = f.variables['longitude'][:]
    baseGPD = f.variables['crs'].long_name
    f.close()

    # find and return the baseGPD (row, col) offset for cubeUL(0, 0) location
    grid = Ease2Transform(baseGPD)
    row_offset, col_offset = grid.geographic_to_grid(lats[0,0], lons[0,0])

    if verbose:
        print("%10s offsets for subsetName=%s and cubeType=%s" % (baseGPD, subsetName, cubeType))
        print("(Add these offsets to cube (row, col) to get row, col in full hemisphere)")
        print("offset row = %f, offset col = %f" % (row_offset, col_offset))
        
    return row_offset, col_offset  


# pass subsetName, latitudes and longitudes,
# returns the rows and columns on the EASE grid for cubefiles,
# uses the 'find_cube_offset' function
# this function assumes that at least one 36 or 37V SIR file,
# one 18 or 19V SIR file and one 18 or 19V GRD file are stored locally
def grid_locations_of_subset(subsetName, lat, lon, cubeDir=None):

    gpds = ["EASE2_N3.125km", "EASE2_N6.25km", "EASE2_N25km"]
    types = ["3?V-SIR", "1?V-SIR", "1?V-GRD"]
    #possibly change to H or V for flexibility
    
    print("Input lat, lon = %.6f, %.6f" % (lat, lon))
    rows = np.zeros((3))
    cols = np.zeros((3))
    for i, (thisGpd, thisType) in enumerate(zip(gpds, types)):
        grid = Ease2Transform(thisGpd)
        row, col = grid.geographic_to_grid(lat, lon)
        print("%15s         : row, col = %.6f, %.6f" % (thisGpd, row, col))
        
        # Get the cube offsets
        offset_row, offset_col = find_cube_offset(
            subsetName,
            cubeDir=cubeDir,
            cubeType=thisType)
        subRow = row - offset_row
        subCol = col - offset_col
        print("%15s(%s): row, col = %.6f, %.6f" % (
            subsetName, thisType, subRow, subCol))
        rows[i] = subRow
        cols[i] = subCol

    rows = (rows + 0.5).astype('int32')
    cols = (cols + 0.5).astype('int32')

    #FIXME: MJB: validate that this routine is returning the correct coordinates
    #FIXME: This routine should really return a dict structure, not an array
    print(gpds)
    print("%d %d %d %d %d %d" % (
        cols[0], rows[0],
        cols[1], rows[1],
        cols[2], rows[2]))
    
    return (rows, cols)


# Function for reading in ascii files - IN PROGRESS
def read_esri_ascii(asc_file, grid=None, reshape=False, name=None, halo=0):
    """Read :py:class:`~landlab.RasterModelGrid` from an ESRI ASCII file.

    Read data from *asc_file*, an ESRI_ ASCII file, into a
    :py:class:`~landlab.RasterModelGrid`.  *asc_file* is either the name of
    the data file or is a file-like object.

    The grid and data read from the file are returned as a tuple
    (*grid*, *data*) where *grid* is an instance of
    :py:class:`~landlab.RasterModelGrid` and *data* is a numpy
    array of doubles with that has been reshaped to have the number of rows
    and columns given in the header.

    .. _ESRI: http://resources.esri.com/help/9.3/arcgisengine/java/GP_ToolRef/spatial_analyst_tools/esri_ascii_raster_format.htm

    Parameters
    ----------
    asc_file : str of file-like
        Data file to read.
    reshape : boolean, optional
        Reshape the returned array, otherwise return a flattened array.
    name : str, optional
        Add data to the grid as a named field.
    grid : *grid* , optional
        Adds data to an existing *grid* instead of creating a new one.
    halo : integer, optional
        Adds outer border of depth halo to the *grid*. 

    Returns
    -------
    (grid, data) : tuple
        A newly-created RasterModel grid and the associated node data.
        
    Raises
    ------
    DataSizeError
        Data are not the same size as indicated by the header file.
    MismatchGridDataSizeError
        If a grid is passed, the size of the grid does not agree with the
        size of the data.
        
    Examples
    --------
    Assume that fop is the name of a file that contains text below
    (make sure you have your path correct):
    ncols         3
    nrows         4
    xllcorner     1.
    yllcorner     2.
    cellsize      10.
    NODATA_value  -9999
    0. 1. 2.
    3. 4. 5.
    6. 7. 8.
    9. 10. 11.
    --------
    >>> from landlab.io import read_esri_ascii
    >>> (grid, data) = read_esri_ascii('fop') # doctest: +SKIP
    >>> #grid is an object of type RasterModelGrid with 4 rows and 3 cols
    >>> #data contains an array of length 4*3 that is equal to
    >>> # [9., 10., 11., 6., 7., 8., 3., 4., 5., 0., 1., 2.]
    >>> (grid, data) = read_esri_ascii('fop', halo=1) # doctest: +SKIP
    >>> #now the data has a nodata_value ring of -9999 around it. So array is
    >>> # [-9999, -9999, -9999, -9999, -9999, -9999,
    >>> #  -9999, 9., 10., 11., -9999, 
    >>> #  -9999, 6., 7., 8., -9999, 
    >>> #  -9999, 3., 4., 5., -9999,
    >>> #  -9999, 0., 1., 2. -9999,
    >>> #  -9999, -9999, -9999, -9999, -9999, -9999]
    """
    from ..grid import RasterModelGrid

    if isinstance(asc_file, six.string_types):
        file_name = asc_file
        with open(file_name, 'r') as asc_file:
            header = read_asc_header(asc_file)
            data = _read_asc_data(asc_file)
    else:
        header = read_asc_header(asc_file)
        data = _read_asc_data(asc_file)
    
    #There is no reason for halo to be negative.
    #Assume that if a negative value is given it should be 0.
    if halo <= 0:
        shape = (header['nrows'], header['ncols'])
        if data.size != shape[0] * shape[1]:
            raise DataSizeError(shape[0] * shape[1], data.size)
    else:
        shape = (header['nrows'] + 2 * halo, header['ncols'] + 2 * halo)
        #check to see if a nodata_value was given.  If not, assign -9999.
        if 'nodata_value' in header.keys():
            nodata_value = header['nodata_value']
        else:
            header['nodata_value'] = -9999.
            nodata_value = header['nodata_value']
        if data.size != (shape[0] - 2 * halo) * (shape[1] - 2 * halo):
            raise DataSizeError(shape[0] * shape[1], data.size)
    spacing = (header['cellsize'], header['cellsize'])
    #origin = (header['xllcorner'], header['yllcorner'])   
    
    data = np.flipud(data)

    #REMEMBER, shape contains the size with halo in place
    #header contains the shape of the original data
    #Add halo below
    if halo > 0:
        helper_row = np.ones(shape[1]) * nodata_value
        #for the first halo row(s), add num cols worth of nodata vals to data
        for i in range(0, halo):
            data = np.insert(data,0,helper_row)
        #then for header['nrows'] add halo number nodata vals, header['ncols'] 
        #of data, then halo number of nodata vals
        helper_row_ends = np.ones(halo) * nodata_value
        for i in range(halo, header['nrows']+halo):
            #this adds at the beginning of the row
            data = np.insert(data,i * shape[1],helper_row_ends)
            #this adds at the end of the row
            data = np.insert(data,(i + 1) * shape[1] - halo,helper_row_ends)
        #for the last halo row(s), add num cols worth of nodata vals to data
        for i in range(header['nrows']+halo,shape[0]):
            data = np.insert(data,data.size,helper_row)
        
    if not reshape:
        data = data.flatten()
        
    if grid is not None:
        if (grid.number_of_node_rows != shape[0]) or \
        (grid.number_of_node_columns != shape[1]):
            raise MismatchGridDataSizeError(shape[0] * shape[1], \
            grid.number_of_node_rows * grid.number_of_node_columns )

    if grid is None:
        grid = RasterModelGrid(shape, spacing=spacing)
    if name:
        grid.add_field('node', name, data)

    return (grid, data)

def years_for(platform):
    """Returns years (and partial years) of operation for passive microwave platform being processed

    years corresponding to platform are returned as an array of integers
    
    Parameters
    ----------
    platform : str of platform ID
       for example 'F11' or 'AQUA'
        
    Returns
    -------
    [year, year, year]
                
    Raises
    ------
    NA
        
    Examples
    --------
    
    """
    if platform=='F8': #F8 SSMI dates Sept. 7, 1987-Dec 30, 1991
        Years = [1987,1988,1989,1990,1991]
    #skip F10; relatively elliptical orbit
    elif platform=='F11': #F11 SSMI dates XXXXX
        Years = [1995,1996,1997,1998,1999,2000,2001,2002,2003, 2004,2005,2006,2007,2008,2009]
    elif platform=='F13' : #F13 SSMI dates May 3,1995-Nov 19,2009
        Years = [1995,1996,1997,1998,1999,2000,2001,2002,2003, 2004,2005,2006,2007,2008,2009]
    elif platform=='F14': #F14 SSMI dates May 7, 1997 - Aug 23, 2008
        Years=[1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008]
    elif platform=='F15': #F15 SSMI dates Feb 23, 2000 - Dec 31 2019
        Years=[2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019]
    elif platform=='F16': #F16 SSMIS dates Nov. 1 2005 - Dec 31 2019
        Years=[2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019] 
    elif platform=='F17': #F17 SSMIS dates
        Years=[2010,2011,2012,2013,2014,2015,2016,2017,2018,2019] 
    elif platform=='F18': #F18 SSMIS dates
        Years=[2010,2011,2012,2013,2014,2015,2016,2017,2018,2019]  
    elif platform=='F18': #F18 SSMIS dates
        Years=[2010,2011,2012,2013,2014,2015,2016,2017,2018,2019]  
    #not including F19 (short time series)
    elif platform=='AQUA':
        Years=[2003,2004,2005,2006,2007,2008,2009,2010,2011]
    else: 
        raise IOError ('unrecognized platform')

    return (Years)




