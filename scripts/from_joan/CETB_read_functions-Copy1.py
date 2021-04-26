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

# offset of Western CA cubes
def find_WesternCA_cube_offset(cubeType=None, verbose=False):
    if not cubeType:
        cubeType = '36V-SIR'
        
    cubeFile = "%s%s%s%s" % (
        "/home/mij216/Desktop/data3/cetb/cubes/AQUA_AMSRE/WesternCA/",
        "CETB.cubefile.WesternCA.AQUA_AMSRE-",
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

# Western CA - pass latitudes and longitudes, returns the rows and columns on the EASE grid for cubefiles, uses the 'find_cube_offset' function 
def grid_locations_of_WesternCA(lat, lon):
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
        offset_row, offset_col = find_WesternCA_cube_offset(
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


