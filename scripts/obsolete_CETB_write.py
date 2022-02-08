from datetime import datetime
from cetbtools.ease2conv import Ease2Transform
import glob
from netCDF4 import Dataset, num2date
import numpy as np
from osgeo import gdal, osr   # noqa
import pandas as pd
import pdb; # insert at places for breakpoints: pdb.set_trace()

# write_MOD_df_column_to_geotiff
# Writes a column of the melt-onset-date dataframe to a geotiff
# Input:
#   df : MOD dataframe, with data columns for lat/lon, row/col, year and Avg
#        and 1 dataframe row for each pixel in the subset
#        assumes data rows are ordered from UL to LR
#   column : column name to write to geotiff
#   grid : Ease2Transform object with projection/grid transformations
#   outbasename : basename for outfile to write, will be appended with
#        column name and '.tif' extension
#   dtype='int16' : data type to write (data will be cast to this type if needed)
#   verbose=False : verbose information during operation
#
def write_MOD_df_column_to_geotiff(df, column, grid, outbasename,
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
    # Documentation for raster GTiff driver here: https://gdal.org/drivers/raster/gtiff.html#raster-gtiff
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

    return outfilename


# write_MOD_df_to_geotiff
# Writes data columns of the melt-onset-date dataframe to a geotiff
# year columns assumed to contain doy (1-366 or NaN) will be written as int16
# Avg column assumed to contain avg doy will be written as float32
# Input:
#   df : MOD dataframe, with data columns for lat/lon, row/col, year and Avg
#        and 1 dataframe row for each pixel in the subset
#        assumes data rows are ordered from UL to LR
#   gpd : projection/grid name, e.g. 'EASE2_N25km', 'EASE2_N6.25km', 'EASE2_N3.125km'
#   outbasename : basename for outfile to write, will be appended with
#        column name and '.tif' extension
#   verbose=False : verbose information during operation
#
def write_MOD_df_to_geotiff(df, gpd, outbasename, verbose=False):

    # initialize the grid once
    grid = Ease2Transform(gpd)

    # Figure out the non-geolocation columns (years or 'Avg')
    columns = df.columns
    columns = columns.drop(['pixel', 'latitude', 'longitude', 'row', 'column'])

    files = []
    for col in columns:

        dtype = 'int16'
        if 'Avg' == col:
            dtype = 'float32'
        outfilename = write_MOD_df_column_to_geotiff(
            df, col, grid, outbasename, dtype=dtype, verbose=verbose)
        files.append(outfilename)

    return files
