# CETB_IO contains input/output and formatting conversion utilities
# for CETB data cubes and derived melt-onset-day images
from cetbtools.ease2conv import Ease2Transform
from datetime import datetime
import glob
from netCDF4 import Dataset, num2date
import numpy as np
from osgeo import gdal, gdal_array, osr   # noqa
import pandas as pd
from pathlib import Path
import pdb # insert at places for breakpoints: pdb.set_trace()
import re 
import warnings


# getting a runtimewarning when using operators
# on numpy arrays with lots of NaNs, functions still perform,
# but using this command to suppress the warning
#FIXME: this shouldn't be necessary
warnings.filterwarnings("ignore",category =RuntimeWarning)


# ###################################################################
# ############################ readers ##############################
# ###################################################################

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
        cube_time, cube_y, cube_x=rawdata.variables['TB'].shape
        if ( y_end > cube_y ) or ( x_end > cube_x ):
            raise ValueError('Subset area not found in this cube')
            
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
        print("Next time filename=%s..." % list[-1])
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


def read_cetb_geotiff(filename, verbose=False):
    """Reads the geotiff image and projected extent from geotiff

    Parameters
    ----------
    filename : name of geotiff

    Returns
    -------
    out : dictionary with
        img : 2-D data image
        extent : dictionary extent in projected coordinates
                 values are corners of corner pixels
        
    Raises
    ------
    n/a
    """
        
    # Read the image data as a numeric array
    img = gdal_array.LoadFile(filename)

    f = gdal.Open(filename, gdal.GA_ReadOnly)
    geoTransform = f.GetGeoTransform()
    minX = geoTransform[0]
    maxY = geoTransform[3]
    maxX = minX + geoTransform[1] * f.RasterXSize
    minY = maxY + geoTransform[5] * f.RasterYSize
    extent = {'minX': minX,
              'maxX': maxX,
              'minY': minY,
              'maxY': maxY}
    f = None

    # FIXME: hardcoding this for now, since all the subsets
    # we are working with are EASE2_N.  If we switch to another
    # hemisphere or start reading generic geotiffs, this
    # should be read directly from the .tif file
    epsg = "EPSG:6931"

    if verbose:
        print('Read geotiff image/extent from %s' % filename)

    return {'img': img,
            'epsg': epsg,
            'extent': extent}

    

# ###################################################################
# ############################ writers ##############################
# ###################################################################

# write_df_column_to_geotiff
# Writes a column of the melt-onset-date or end-of-high-DAV dataframe 
# to a geotiff
# Input:
#   df : MOD (or EHD) dataframe, with data columns for lat/lon, row/col, year and Avg
#        and 1 dataframe row for each pixel in the subset
#        assumes data rows are ordered from UL to LR
#   column : column name to write to geotiff
#   grid : Ease2Transform object with projection/grid transformations
#   outbasename : basename for outfile to write, will be appended with
#        column name and '.tif' extension
#   dtype='int16' : data type to write (data will be cast to this type if needed)
#   verbose=False : verbose information during operation
#
def write_df_column_to_geotiff(df, column, grid, outbasename,
                               dtype='int16', verbose=False):

    outfilename = "%s.%s.tif" % (outbasename, column)
    nrows = int(df.iloc[-1].row - df.iloc[0].row + 1)
    ncols = int(df.iloc[-1].column - df.iloc[0].column + 1)
    
    # Coerce array data to requested dtype
    data = np.array(df[column].values.data).reshape(nrows, ncols)
    data = data.astype(dtype)
    if verbose:
        print("Writing %s data as %s (%d, %d) to %s..." % (
            column, data.dtype, nrows, ncols, outfilename))
    
    if ("float32" == data.dtype):
        gdal_data_type = gdal.GDT_Float32
    elif ("int8" == data.dtype):
        gdal_data_type = gdal.GDT_Byte
    elif ("int16" == data.dtype):
        gdal_data_type = gdal.GDT_Int16
    else:
        print("%s : Unrecognized type %s " %                                                 
              (my_name, str(data.dtype)))
        raise ValueError
    # Initialize the output driver
    # Documentation for raster GTiff driver here:
    # https://gdal.org/drivers/raster/gtiff.html#raster-gtiff
    driver = gdal.GetDriverByName("GTiff")
    
    # use this to control block sizes:
    # dest_ds_options = ['COMPRESS=LZW', 'TILED=YES', 'BLOCKXSIZE=256', 'BLOCKYSIZE=256']  
    dest_ds_options = ['COMPRESS=LZW']
    dest_ds = driver.Create(outfilename, ncols, nrows, 1, gdal_data_type, dest_ds_options)
    
    # Initialize the projection information
    # The crs.proj4text attribute can also be used here,
    # but the srid attribute provides more complete PROJCRS metadata
    proj = osr.SpatialReference()
    proj.SetFromUserInput(grid.epsg)
    dest_ds.SetProjection(proj.ExportToWkt())
    
    # Initialize the grid information (extent and scale)                                        
    # Thanks to web page at:                                                                    
    # http://geoexamples.blogspot.com/2012/01/                                                  
    # creating-files-in-ogr-and-gdal-with.html                                                  
    # The geotransform defines the relation between the                                         
    # raster coordinates x, y and the                                                           
    # geographic coordinates, using the following definition:                                   
    # Xgeo = geotransform[0] + Xpixel*geotransform[1] + Yline*geotransform[2]                   
    # Ygeo = geotransform[3] + Xpixel*geotransform[4] + Yline*geotransform[5]                   
    # The first and fourth parameters define the origin of the upper left pixel                 
    # The second and sixth parameters define the pixels size.                                   
    # The third and fifth parameters define the rotation of the raster.                         
    # Values are meters                                                                         
    # The UL information is the center of the UL corner pixel in projected
    # coordinates
    ULrow, ULcol = grid.geographic_to_grid(df.iloc[0].latitude, df.iloc[0].longitude)
    ULrow = int(ULrow + 0.5)
    ULcol = int(ULcol + 0.5)
    
    map_ULx, map_ULy = grid.grid_to_map(ULrow, ULcol)
    if verbose:
        print("UL: ", ULrow, ULcol,
              df.iloc[0].latitude, df.iloc[0].longitude, map_ULx, map_ULy)
    
    LRrow, LRcol = grid.geographic_to_grid(df.iloc[-1].latitude, df.iloc[-1].longitude)
    LRrow = int(LRrow + 0.5)
    LRcol = int(LRcol + 0.5)
    map_LRx, map_LRy = grid.grid_to_map(LRrow, LRcol)
    if verbose:
        print("LR: ", LRrow, LRcol,
              df.iloc[-1].latitude, df.iloc[-1].longitude, map_LRx, map_LRy)
    
    # Get the projection scales by dividing the projected extents from UL and LR pixels
    # by row/col dimensions
    scale_x = (map_LRx - map_ULx) / np.float32(ncols - 1)
    scale_y = -1. * (map_ULy - map_LRy) / np.float32(nrows - 1)
    if verbose:
        print('scale x, y = %f, %f' % (scale_x, scale_y))

    cornerULx = map_ULx - (scale_x / 2.)
    cornerULy = map_ULy - (scale_y / 2.)
    if verbose:
        print("cornerUL x,y: ", cornerULx, cornerULy)
    geotransform = (cornerULx, scale_x, 0., cornerULy, 0., scale_y)
    dest_ds.SetGeoTransform(geotransform)                                                       
    dest_ds.GetRasterBand(1).WriteArray(data) 
                                                     
    dest_ds = None

    if verbose:                                                                                 
        print("\n%s geotiff image saved to: %s" %                                          
              (str(column), outfilename))

    return {'filename': outfilename,
            'geotransform': geotransform,
            'data': data}


# write_df_to_geotiff
# Writes data columns of the melt-onset-date or end-of-high-DAV dataframe to a geotiff
# year columns assumed to contain doy (1-366 or NaN) will be written as int16
# Avg column assumed to contain avg doy will be written as float32
# Input:
#   df : MOD or EHD dataframe, with data columns for lat/lon, row/col, year and Avg
#        and 1 dataframe row for each pixel in the subset
#        assumes data rows are ordered from UL to LR
#   gpd : projection/grid name, e.g. 'EASE2_N25km', 'EASE2_N6.25km', 'EASE2_N3.125km'
#   outbasename : basename for outfile to write, will be appended with
#        column name and '.tif' extension
#   verbose=False : verbose information during operation
#
def write_df_to_geotiff(df, gpd, outbasename, verbose=False):

    # initialize the grid once
    grid = Ease2Transform(gpd)

    # Figure out the non-geolocation columns (years or 'Avg')
    columns = df.columns
    columns = columns.drop(['pixel',
                            'latitude', 'longitude',
                            'row', 'column',
                            'x', 'y'])

    out = {}
    for col in columns:

        dtype = 'int16'
        if 'Avg' == col:
            dtype = 'float32'
        thisOut = write_df_column_to_geotiff(
            df, col, grid, outbasename, dtype=dtype, verbose=verbose)
        out[col] = thisOut

    return out


# ###################################################################
# ####################### geolocation utilities #####################
# ###################################################################

def coords(datadir, prefix, lat_start, lat_end, lon_start, lon_end):
    # converts input latitude and longitude area to row/col coordinates
    # in the cubefiles described by datadir and prefix
    # datadir : path to cubefiles
    # prefix : string wildcard to use to look for 'datadir/prefix.*.TB.nc'
    #          cubefiles
    # lat_start, lat_end : begin/end of latitude range
    # lon_start, lon_end : begin/end of longitude range
    #
    # returns: row/col coordinates for smallest rectangle in the cube's grid
    #          that includes the input lat/lon box, as
    #          
    #          row_start, row_end, col_start, col_end
    #      
    #          where (row,col) = (0, 0) refers to the UL cell of the cube area
    #   

    # Look for TB.nc cubefiles in the specified location
    # Assumes all cubes found here are for the same cube coverage,
    # so looking at the last one is good enough
    filename = str(Path(datadir, "%s*TB.nc" % prefix))
    files = glob.glob(filename)
    if not files:
        raise RuntimeError(
            "%s : can't find TB.nc cubefiles to match %s" % (
                __file__, filename))
    
    f = Dataset(files[-1], "r", format="NETCDF4")

    # Read the last file found to get the gpd name to use
    gpd = f.variables["crs"].long_name
    grid = Ease2Transform(gpd)

    # ...and get the lat/lon of the UL grid cell
    cube_UL_lat = f.variables["latitude"][0]
    cube_UL_lon = f.variables["longitude"][0]

    f.close()

    # Convert the lat/lon of the UL cube cell to row/col offset in the full hemisphere grid
    offset_row, offset_col = grid.geographic_to_grid(cube_UL_lat[0], cube_UL_lon[0])
    offset_row = np.int32(offset_row + 0.5)
    offset_col = np.int32(offset_col + 0.5)

    # verify that lat_start <= lat_end and lon_start <= lon_end, swap if necessary
    if lat_end < lat_start:
        tmp = lat_start
        lat_start = lat_end
        lat_end = tmp

    if lon_end < lon_start:
        tmp = lon_start
        lon_start = lon_end
        lon_end = tmp

    # make arrays of area corner lat/lons, start with NW corner, go clockwise
    num = 4
    corner_lats = np.array([lat_end, lat_end, lat_start, lat_start], dtype=np.float64)
    corner_lons = np.array([lon_start, lon_end, lon_end, lon_start], dtype=np.float64)
    corner_cols = np.zeros(num)
    corner_rows = np.zeros(num)

    # convert the area corners to integer row/col in the full hemisphere grid
    for i in np.arange(num):
        corner_rows[i], corner_cols[i] = grid.geographic_to_grid(
            corner_lats[i], corner_lons[i])
        
    corner_rows = np.int32(corner_rows + 0.5)
    corner_cols = np.int32(corner_cols + 0.5)

    # Get the smallest rectangle in grid space that encloses
    # the requested lat/lon area
    # Convert to offsets relative to UL of cube
    row_start = np.min(corner_rows) - offset_row
    row_end = np.max(corner_rows) - offset_row
    col_start = np.min(corner_cols) - offset_col
    col_end = np.max(corner_cols) - offset_col

    # The returned row_end and col_end are 1 pixel further than the
    # required area, because the subsequent subsetting is open interval 
    row_end = row_end + 1
    col_end = col_end + 1

    # Returned values are row/col relative to the UL=(0,0) cell in the cube
    return row_start, row_end, col_start, col_end
    

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
# returns the rows and columns on the EASE grid for each (25, 6.25, 3.125) cubefiles,
# uses the 'find_cube_offset' function
# this function assumes that at least one 36 or 37V SIR file,
# one 18 or 19V SIR file and one 18 or 19V GRD file are stored locally
def grid_locations_of_subset(subsetName, lat, lon, cubeDir=None):

    gpds = ["EASE2_N3.125km", "EASE2_N6.25km", "EASE2_N25km"]
    types = ["3?V-SIR", "1?V-SIR", "1?V-GRD"]
    #possibly change to H or V for flexibility

    dims = dict.fromkeys(['row', 'col'])
    out = dict.fromkeys(gpds)
    for g in gpds:
        out[g] = dims.copy()
    
    print("Input lat, lon = %.6f, %.6f" % (lat, lon))
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
        out[thisGpd]['row'] = np.int32(subRow + 0.5)
        out[thisGpd]['col'] = np.int32(subCol + 0.5)
    
    return out


# get_nearest_ease2_coordinates(grid, lat, lon)
# Input:
#   grid: Ease2Transform for the EASE-Grid 2.0 grid
#   lat/lon: latitude/longitude of a point
# Output:
#   Dictionary with row/col and x/y of nearest grid pixel to
#   lat/lon
def get_nearest_ease2_coordinates(grid, lat, lon):

    # Get the row/col - these are real-valued (exact location)
    row, col = grid.geographic_to_grid(lat, lon)

    # Round the row/col coordinates to the nearest integer
    # this will be location of the center of the nearest cell
    row = int(row + 0.5)
    col = int(col + 0.5)

    # Get projected x, y of the center of the nearest cell
    x, y = grid.grid_to_map(row, col)

    return {'row': row,
            'col': col,
            'x': x,
            'y': y}


# get_extent_xy(x, y)
# Given a cube subset x, y (1-D) arrays,
# calculate the projected extent of the corners of the corner pixels
# Input:
#   x, y : 1-D array of x, y locations of centers of cells
# Returns:
#   extent = dictionary with full extent (corners of corner cells)
#            ULx, LRx, ULy, LRy
def get_extent_xy(x, y):
    # Assumes pixels with regular spacing
    scale_x = x[1] - x[0]
    scale_y = y[1] - y[0]
    print("scales: ", scale_x, scale_y)
    
    extent = {'ULx': x[0] - (scale_x/2.),
              'LRx': x[-1] + (scale_x/2.),
              'ULy': y[-1] + (scale_y/2.),
              'LRy': y[0] - (scale_y/2.)}
              
    return extent


# ###################################################################
# ###################### miscellaneous utilities  ###################
# ###################################################################

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
        Years = [1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,
                 2005,2006,2007,2008,2009]
    elif platform=='F13' : #F13 SSMI dates May 3,1995-Nov 19,2009
        Years = [1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,
                 2005,2006,2007,2008,2009]
    elif platform=='F14': #F14 SSMI dates May 7, 1997 - Aug 23, 2008
        Years = [1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,
                 2007,2008]
    elif platform=='F15': #F15 SSMI dates Feb 23, 2000 - Dec 31 2019
        Years = [2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,
                 2010,2011,2012,2013,2014,2015,2016,2017,2018,2019]
    elif platform=='F16': #F16 SSMIS dates Nov. 1 2005 - Dec 31 2019
        Years = [2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,
                 2015,2016,2017,2018,2019] 
    elif platform=='F17': #F17 SSMIS dates
        Years = [2010,2011,2012,2013,2014,2015,2016,2017,2018,2019] 
    elif platform=='F18': #F18 SSMIS dates
        Years = [2010,2011,2012,2013,2014,2015,2016,2017,2018,2019]  
    #not including F19 (short time series)
    elif platform=='AQUA':
        Years = [2003,2004,2005,2006,2007,2008,2009,2010,2011]
    elif platform=='GCOMW1':
        Years = [2021,2022]
    else: 
        raise IOError ('unrecognized platform')

    return (Years)


def get_sir_info(channel, hem='N'):
    """Returns the sir-to-grd pixel factor and sir gpd name for the specified channel
    
    Parameters
    ----------
    channel : str channel name, may include polarization, which will be ignored,
              e.g. '19V', '1.4H'
        
    Returns
    -------
    (factor, sir_gpd)
                
    Raises
    ------
    NA
        
    Examples
    --------
    
    """

    # set the sir to grd factor, depends on the channel
    chans3km = ['^1.4', '^36', '^37', '^85', '^89', '^91']
    chans6km = ['^18', '^19', '^21', '^22', '^23']
    chans12km = ['^6', '^10']
    if re.search('|'.join(chans3km), channel):
        factor = 8 # assume 3.125 km to 25 km
        sir_gpd = 'EASE2_%s3.125km' % hem
    elif re.search('|'.join(chans6km), channel):
        factor = 4 # assume 6.25 km to 25 km
        sir_gpd = 'EASE2_%s6.25km' % hem
    elif re.search('|'.join(chans12km), channel):
        factor = 2 # assume 12.5 km to 25 km
        sir_gpd = 'EASE2_%s12.5km' % hem        
    else:
        raise ValueError("Cannot determine sir-to-grd factor from channel %s\n" % (channel) )
    
    return (factor, sir_gpd)


def get_site_boundaries(SiteLabel):
    """ Given a site label short string returns the lat/lon boundaries for that area
            - could be a single pixel or an area
            - also returns the full Site name for graphs etc
        Parameters
        -----------
        SiteLabel - string - short site label with no punctuation or spaces

        Returns
        -------
        lat_start, lat_end, lon_start, lon_end, Site
        
    """
    
    dict = {
        'NEUkraine':{
            'Site':'NE, Ukraine',
            'lat_start':51.050,
            'lon_start':34.500,
            'lat_end':51.050,
            'lon_end':34.500},
        'ChernihivForested':{
            'Site':'Chernihiv Oblast Forested, Ukraine',
            'lat_start':51.907788,
            'lon_start':31.295555,
            'lat_end':51.907788,
            'lon_end':31.295555},
        'ChernihivAg':{
            'Site':'Chernihiv Oblast Agricultural, Ukraine',
            'lat_start':51.698747,   
            'lon_start':31.449069,
            'lat_end':51.698747,   
            'lon_end':31.449069},
        'KomsomoletsIs':{
            'Site':'Komsomolets Is., Severnaya Zemlya AWS, Russia',
            'lat_start':80.51666667   ,
            'lon_start':94.81666667,
            'lat_end':80.51666667   ,
            'lon_end':94.81666667},
        'Alberta':{
            'Site':'Low Relief Test, Alberta',
            'lat_start':49.73,
            'lon_start':-111.14,
            'lat_end':49.73,
            'lon_end':-111.14},
        'GreatLakes':{
            'Site':'Great Lakes Test Site, Lake Superior E of Duluth',
            'lat_start':47.353097   ,
            'lon_start':-89.861310,
            'lat_end':47.353097   ,
            'lon_end':-89.861310},
 #AKYukon Sites
 # 3 km SE of the airport pixel, GRD and SIRs should not be coastal
 #Note is within the same GRD as the airport       
        'Barrow3kmSEAirport':{
            'Site':'Barrow SE, AK', 
            'lat_start':71.2709,
            'lon_start':-156.694,
            'lat_end':71.2709,
            'lon_end':-156.694},
        'Barrow airport, AK':{
            'Site':'BarrowAirport',
            'lat_start':71.28181,
            'lon_start':-156.772,
            'lat_end':71.28181,
            'lon_end':-156.772},
        'Barrow45kmAirport':{
            'Site':'45 km S of Barrow airport, AK',
            'lat_start':70.908886,
            'lon_start':-156.56875,
            'lat_end':70.908886,
            'lon_end':-156.56875},
        'Barrow114kmS':{
            'Site':'114 km South of Barrow airport, AK', #114 km south of airport to get an inland site
            'lat_start':70.28181,
            'lon_start':-156.772,
            'lat_end':70.28181,
            'lon_end':-156.772},
        'UpperKuparuk':{
            'Site':'Upper Kuparuk Basin, AK', #selected by VJ
            'lat_start':68.62421667,
            'lon_start':-149.5250722,
            'lat_end':68.62421667,
            'lon_end':-149.5250722},
        'TechsekpukLake':{
            'Site':'Teshekpuk Lake 1',#near Barrow, AK (Permafrost installation NW margin of lake)'
            'lat_start':70.722902,
            'lon_start':-153.836329,
            'lat_end':70.722902,
            'lon_end':-153.836329},
        'TechsekpukLake2':{
            'Site':'Teshekpuk Lake 2; center', #, near Barrow, AK (middle of lake)',
            'lat_start':70.58333675326588,
            'lon_start':-153.452493001114,
            'lat_end':70.58333675326588,
            'lon_end':-153.452493001114},
        
        
     #Western_CA Sites
        'NWT':{
            'Site':'NWT C57 04242006 Spring Migration',
            'lat_start':62.61,
            'lon_start':-109.64,
            'lat_end':62.61,
            'lon_end':-109.64},
        'GreatSlaveLake':{
            'Site':'Great Slave Lake',
            'lat_start':61.87167,
            'lon_start':-114.05,
            'lat_end':61.87167,
            'lon_end':-114.05},
         #GLAIL Sites
        'IcelandCentralVatnajokull':{
            'Site':'Iceland, Central Vatnajokull',
            'lat_start':64.442   ,
            'lon_start':-16.730 ,
            'lat_end':64.442   ,
            'lon_end':-16.730 },
        'IcelandKatlaCaldera':{
            'Site':'Iceland, Katla Caldera Test',
            'lat_start':63.64361111   ,
            'lon_start':-19.20138889,
            'lat_end':63.64361111   ,
            'lon_end':-19.20138889},
        'IcelandVatnajokull1':{
            'Site':'Iceland, Vatnajokull Site 1',
            'lat_start':64.64903   ,
            'lon_start':-17.32128,
            'lat_end':64.64903   ,
            'lon_end':-17.32128},
        'IcelandVatnajokull2':{
            'Site': 'Iceland, Vatnajokull Site 2',
            'lat_start':64.25860278,
            'lon_start':-16.93241667,
            'lat_end':64.25860278,
            'lon_end':-16.93241667},
        'IcelandVatnajokull3':{
            'Site':'Iceland, Vatnajokull Site 3',
            'lat_start':64.423333,
            'lon_start':-16.77305556,
            'lat_end':64.423333,
            'lon_end':-16.77305556},
        'IcelandVatnajokull4':{
            'Site':'Iceland, Vatnajokull Site 4',
            'lat_start':64.66222,
            'lon_start':-17.399444,
            'lat_end':64.66222,
            'lon_end':-17.399444},
             #Western_US Site
        'SASP': {
            'Site':'Senator Beck Basin (SASP), Colorado',
            'lat_start':37.9069   ,
            'lon_start':-107.710,
            'lat_end':37.9069   ,
            'lon_end':-107.710},
        'RabbitEars':{
            'Site':'Rabbit Ears',
            'lat_start':40.5 , #southern boundary
            'lat_end':40.5  ,#northern boundary
            'lon_start':-106.7,  #western boundary
            'lon_end':-106.7}, #eastern boundary
        'Fraser':{
            'Site':'Fraser',
            'lat_start':39.85,  #southern boundary
            'lat_end':40.0  ,#northern boundary
            'lon_start':-106.0,  #western boundary
            'lon_end':-105.9}, #eastern boundary
        'NorthPark':{
            'Site':'North Park',
            'lat_start':40.7  ,#southern boundary
            'lat_end':40.7  ,#northern boundary
            'lon_start':-106.15 , #western boundary
            'lon_end':-106.15 },#eastern boundary
        'CLPX-LSA':{
            'Site':'CLPX-LSA',
            'lat_start':39.7 , #southern boundary
            'lat_end':41.0  ,#northern boundary
            'lon_start':-106.7 , #western boundary
            'lon_end':-105.4 },#eastern boundary
        'CUSTOM':{
            'Site':'Custom',
            'lat_start':43 , #southern boundary
            'lat_end':43  ,#northern boundary
            'lon_start':-110 , #western boundary
            'lon_end':-110 },
        #eastern boundary
        'SenatorBeck':{
            'Site':'Senator Beck',
             # enter latitutde and longitude in decimal degrees
            'lat_start':37.9069  , #southern boundary
            'lat_end':37.9069   ,#northern boundary
            'lon_start':-107.726 ,  #western boundary
            'lon_end':-107.726 },
        'vatna':{
            'lat_start':63.75  ,
            'lat_end':64.88    ,
            'lon_start':-20 ,
            'lon_end':-15  ,
            'Site':'Vatnajokull, Iceland'},
        'hunza':{
            'lat_start':35.9  ,
            'lat_end':37.1   ,
            'lon_start':74 ,
            'lon_end':76 ,
            'Site':'Hunza Basin'},
        'gsl':{
            'lat_start':59.00  ,
            'lat_end':67.00   ,
            'lon_start':-119.00, 
            'lon_end':-107.00,
            'Site':'Great Slave Lake, Canada'},
        'bathurst_range':{
            'lat_start':60.00  ,
            'lat_end':67.25   ,
            'lon_start':-119.00, 
            'lon_end':-107.50,
            'Site':'Bathurst Caribou Range, NWT'},
        'bathurst_range2':{
            'lat_start':63.00  ,
            'lat_end':65.500   ,
            'lon_start':-117.500, 
            'lon_end':-112.00,
            'Site':'Bathurst Caribou Range subset, NWT'},
        'barrow':{
            'lat_start':69.50  ,
            'lat_end':71.50   , 
            'lon_start':-158 ,
            'lon_end':-152  ,
            'Site':'Barrow/Utkiagvik, AK'  },
        'fairbanks':{
            'lat_start':63.0,  
            'lat_end':66.7 ,   
            'lon_start':-151.8,
            'lon_end':-143.4  ,
            'Site':'Fairbanks, AK'}
    }
            
    
    try:
        return (dict[SiteLabel]['lat_start'],
        dict[SiteLabel]['lat_end'],
        dict[SiteLabel]['lon_start'],
        dict[SiteLabel]['lon_end'],
        dict[SiteLabel]['Site'])
    except KeyError as e:
        print("Invalid SiteLabel - valid values are")
        for k in dict.keys():
            print("%30s lat: %9.4f,%9.4f lon: %9.4f,%9.4f  %s"%(
                k,dict[k]['lat_start'], dict[k]['lat_end'], dict[k]['lon_start'], dict[k]['lon_end'], dict[k]['Site']))
    else:
        raise
       
