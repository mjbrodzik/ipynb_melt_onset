from netCDF4 import Dataset, num2date
import numpy as np
from cetbtools.ease2conv import Ease2Transform
import warnings

# getting a runtimewarning when using operators on numpy arrays with lots of NaNs, functions still perform, but using this command to suppress the warning
warnings.filterwarnings("ignore",category =RuntimeWarning)

# load CETB data for a subset and sensor of interest
def read_Tb(datadir, prefix, Years,y_start,y_end,x_start,x_end):
	for year in Years:
		# Create filename
		filename=datadir+prefix+'.'+str(year)+'.TB.nc'
	    
		# load the raw data in
		rawdata = Dataset(filename, "r", format="NETCDF4")
		
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

		# Handle missing data - Hard coded!!
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

# load ALL TB data for a cubefile
def read_Tb_all(datadir, prefix, Years):
	for year in Years:
		# Create filename
		filename=datadir+prefix+'.'+str(year)+'.TB.nc'
	    
		# load the raw data in
		rawdata = Dataset(filename, "r", format="NETCDF4")
		
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

		# Handle missing data - Hard coded!!
		CETB_data[CETB_data>=60000] = np.NaN
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
		filename=datadir+prefix+'.'+str(year)+'.TB_std_dev.nc'
	    
		# load the raw data in
		rawdata = Dataset(filename, "r", format="NETCDF4")
		
		# Compile the CETB data, the TB variable is saved as (time, y, x) - 
		subset = rawdata.variables['TB_std_dev'][0:,y_start:y_end,x_start:x_end]
		if year==Years[0]:
			CETB_data = subset
		else:
			CETB_data = np.concatenate((CETB_data, subset), axis=0)
		
		# Handle missing data - Hard coded!!
		CETB_data[CETB_data>=60000] = np.NaN
		CETB_data[CETB_data==0] = np.NaN

	return CETB_data

# read in the Tb_time variable - time of the satellite overpass
def read_Tb_time(datadir, prefix, Years,y_start,y_end,x_start,x_end):
	for year in Years:
	# Create filename
		filename=datadir+prefix+'.'+str(year)+'.TB_time.nc'
	   
		# load the raw data in
		rawdata = Dataset(filename, "r", format="NETCDF4")
	
		# Compile the CETB data, the TB variable is saved as (time, y, x) - 
		subset = rawdata.variables['TB_time'][0:,y_start:y_end,x_start:x_end]
		if year==Years[0]:
			CETB_data = subset
		else:
			CETB_data = np.concatenate((CETB_data, subset), axis=0)
		
		# Handle missing data - Hard coded!!
		#CETB_data[CETB_data>60000] = np.NaN
		#CETB_data[CETB_data==0] = np.NaN

	return CETB_data

def coords(datadir, prefix, lat_start, lat_end, lon_start, lon_end):
	# this function takes the user inputs of latitude and longitude and returns the row/col numbers that identify Tb grids in the cubefile	
	# NOTE: this is a function I wrote, but it does basically the same thing as "find_cube_offset" and "grid_locationsof..." which are listed below and were coded by Mary Jo 
	year=2003    #could choose any year since the structure is the same, chose 2003 since this year is included for all sensors
	filename=datadir+prefix+'.'+str(year)+'.TB.nc'	
	data=Dataset(filename, "r", format="NETCDF4")
	
	lat=data.variables['latitude'] 
	lat=lat[:]
	lon=data.variables['longitude']
	lon=lon[:]
	lat_lon=np.dstack((lat,lon))

	# if latitude and longitude are the same, for reading one pixel
	if ((lat_start==lat_end) & (lon_start==lon_end)):	
		y_start=np.min(np.argmin(((lat_start-lat_lon[:,:,0])**2+(lon_start-lat_lon[:,:,1])**2), axis=0))
		y_end=y_start+1
		x_start=np.max(np.argmin(((lat_start-lat_lon[:,:,0])**2+(lon_start-lat_lon[:,:,1])**2), axis=1))
		x_end=x_start+1
		return y_start,y_end,x_start,x_end
	else:
		row_col_set=np.where((lat_lon[:,:,0]>lat_start)&(lat_lon[:,:,0]<lat_end)&(lat_lon[:,:,1]>lon_start)&(lat_lon[:,:,1]<lon_end))	
		y_start=np.min(row_col_set[0])
		y_end=np.max(row_col_set[0])
		x_start=np.min(row_col_set[1])
		x_end=np.max(row_col_set[1])
		if x_start==x_end:
			x_end=x_end+1
		if y_start==y_end:		
			y_end=y_end+1
		return y_start,y_end,x_start,x_end

def calc_DAV(CETB_data):
	# function takes the CETB_data that was read in read_Tb() and returns the absolute value of the DAV	
	DAV=np.diff(CETB_data,n=1,axis=0)
	DAV_abs=np.absolute(DAV)	
	DAV_abs=np.insert(DAV_abs, [0],[0], axis=0)	# insert a 0 at beginning of array so same length as CETB_data for plotting together
	return DAV_abs

# find the offset to be used for locating the grid row/col from lat/lon coordinates (GLaIL)
def find_GLaIL_cube_offset(cubeType=None, verbose=False):
    if not cubeType:
        cubeType = '36V-SIR'
        
    cubeFile = "%s%s%s%s" % (
        "/home/mij216/Desktop/data3/cetb/cubes/AQUA_AMSRE/GLaIL/",
        "CETB.cubefile.GLaIL.AQUA_AMSRE-",
        cubeType,
        "-RSS-v1.3.2003.TB.nc")
        
    f = Dataset(cubeFile, "r", "NETCDF4")   
    lats = f.variables['latitude'][:]
    lons = f.variables['longitude'][:]
    baseGPD = f.variables['crs'].long_name
    f.close()
    
    # find and return the baseGPD (row, col) offset for cubeUL(0, 0) location
    grid = Ease2Transform(baseGPD)
    row_offset, col_offset = grid.geographic_to_grid(lats[0,0], lons[0,0])

    if verbose:
        print("%10s offsets for cubeType=%s" % (baseGPD, cubeType))
        print("(Add these offsets to cube (row, col) to get row, col in full hemisphere)")
        print("offset row = %f, offset col = %f" % (row_offset, col_offset))
        
    return row_offset, col_offset  

# offset for Western US cubes
def find_WesternUS_cube_offset(cubeType=None, verbose=False):
    if not cubeType:
        cubeType = '36V-SIR'
        
    cubeFile = "%s%s%s%s" % (
        "/home/mij216/Desktop/data3/cetb/cubes/AQUA_AMSRE/WesternUS/",
        "CETB.cubefile.WesternUS.AQUA_AMSRE-",
        cubeType,
        "-RSS-v1.3.2003.TB.nc")
        
    f = Dataset(cubeFile, "r", "NETCDF4")   
    lats = f.variables['latitude'][:]
    lons = f.variables['longitude'][:]
    baseGPD = f.variables['crs'].long_name
    f.close()
	# find and return the baseGPD (row, col) offset for cubeUL(0, 0) location
    grid = Ease2Transform(baseGPD)
    row_offset, col_offset = grid.geographic_to_grid(lats[0,0], lons[0,0])

    if verbose:
        print("%10s offsets for cubeType=%s" % (baseGPD, cubeType))
        print("(Add these offsets to cube (row, col) to get row, col in full hemisphere)")
        print("offset row = %f, offset col = %f" % (row_offset, col_offset))
        
    return row_offset, col_offset 

# offset for Upper Indus Basin (UIB) cubes
def find_UIB_cube_offset(cubeType=None, verbose=False):
    if not cubeType:
        cubeType = '36V-SIR'
        
    cubeFile = "%s%s%s%s" % (
        "/home/mij216/Desktop/data3/cetb/cubes/AQUA_AMSRE/UIB/",
        "CETB.cubefile.UIB.AQUA_AMSRE-",
        cubeType,
        "-RSS-v1.3.2003.TB.nc")
        
    f = Dataset(cubeFile, "r", "NETCDF4")   
    lats = f.variables['latitude'][:]
    lons = f.variables['longitude'][:]
    baseGPD = f.variables['crs'].long_name
    f.close() 

# find and return the baseGPD (row, col) offset for cubeUL(0, 0) location
    grid = Ease2Transform(baseGPD)
    row_offset, col_offset = grid.geographic_to_grid(lats[0,0], lons[0,0])

    if verbose:
        print("%10s offsets for cubeType=%s" % (baseGPD, cubeType))
        print("(Add these offsets to cube (row, col) to get row, col in full hemisphere)")
        print("offset row = %f, offset col = %f" % (row_offset, col_offset))
        
    return row_offset, col_offset 

# GLaIL - pass latitudes and longitudes, returns the rows and columns on the EASE grid for cubefiles, uses the 'find_cube_offset' function 
def grid_locations_of_GLaIL(lat, lon):
    gpds = ["EASE2_N3.125km", "EASE2_N6.25km", "EASE2_N25km"]
    types = ["36H-SIR", "18H-SIR", "18H-GRD"]
    
    print("Input lat, lon = %.6f, %.6f" % (lat, lon))
    rows = np.zeros((3))
    cols = np.zeros((3))
    for i, (thisGpd, thisType) in enumerate(zip(gpds, types)):
        grid = Ease2Transform(thisGpd)
        row, col = grid.geographic_to_grid(lat, lon)
        print("%15s         : row, col = %.6f, %.6f" % (thisGpd, row, col))
        
        # Get the cube offsets
        offset_row, offset_col = find_GLaIL_cube_offset(
            cubeType=thisType)
        UIBrow = row - offset_row
        UIBcol = col - offset_col
        print("%15s(%s): row, col = %.6f, %.6f" % ("GlaIL", thisType, UIBrow, UIBcol))
        rows[i] = UIBrow
        cols[i] = UIBcol

    rows = (rows + 0.5).astype('int32')
    cols = (cols + 0.5).astype('int32')

    print(gpds)
    print("%d %d %d %d %d %d" % (
        cols[0], rows[0],
        cols[1], rows[1],
        cols[2], rows[2]))
    
    
    return (rows, cols)

# Western US - pass latitudes and longitudes, returns the rows and columns on the EASE grid for cubefiles, uses the 'find_cube_offset' function 
def grid_locations_of_WesternUS(lat, lon):
    gpds = ["EASE2_N3.125km", "EASE2_N6.25km", "EASE2_N25km"]
    types = ["36H-SIR", "18H-SIR", "18H-GRD"]
    
    print("Input lat, lon = %.6f, %.6f" % (lat, lon))
    rows = np.zeros((3))
    cols = np.zeros((3))
    for i, (thisGpd, thisType) in enumerate(zip(gpds, types)):
        grid = Ease2Transform(thisGpd)
        row, col = grid.geographic_to_grid(lat, lon)
        print("%15s         : row, col = %.6f, %.6f" % (thisGpd, row, col))
        
        # Get the cube offsets
        offset_row, offset_col = find_WesternUS_cube_offset(
            cubeType=thisType)
        UIBrow = row - offset_row
        UIBcol = col - offset_col
        print("%15s(%s): row, col = %.6f, %.6f" % ("WesternUS", thisType, UIBrow, UIBcol))
        rows[i] = UIBrow
        cols[i] = UIBcol

    rows = (rows + 0.5).astype('int32')
    cols = (cols + 0.5).astype('int32')

    print(gpds)
    print("%d %d %d %d %d %d" % (
        cols[0], rows[0],
        cols[1], rows[1],
        cols[2], rows[2]))
    
    
    return (rows, cols)

# UIB - pass latitudes and longitudes, returns the rows and columns on the EASE grid for cubefiles, uses the 'find_cube_offset' function 
def grid_locations_of_UIB(lat, lon):
    gpds = ["EASE2_N3.125km", "EASE2_N6.25km", "EASE2_N25km"]
    types = ["36H-SIR", "18H-SIR", "18H-GRD"]
    
    print("Input lat, lon = %.6f, %.6f" % (lat, lon))
    rows = np.zeros((3))
    cols = np.zeros((3))
    for i, (thisGpd, thisType) in enumerate(zip(gpds, types)):
        grid = Ease2Transform(thisGpd)
        row, col = grid.geographic_to_grid(lat, lon)
        print("%15s         : row, col = %.6f, %.6f" % (thisGpd, row, col))
        
        # Get the cube offsets
        offset_row, offset_col = find_UIB_cube_offset(
            cubeType=thisType)
        UIBrow = row - offset_row
        UIBcol = col - offset_col
        print("%15s(%s): row, col = %.6f, %.6f" % ("UIB", thisType, UIBrow, UIBcol))
        rows[i] = UIBrow
        cols[i] = UIBcol

    rows = (rows + 0.5).astype('int32')
    cols = (cols + 0.5).astype('int32')

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



