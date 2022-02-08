from datetime import datetime
from cetbtools.ease2conv import Ease2Transform
import glob
from netCDF4 import Dataset, num2date
import numpy as np
from osgeo import gdal, osr   # noqa
import pandas as pd
import pdb; # insert at places for breakpoints: pdb.set_trace()
import sys;

# get_extent_xy(x, y)
# Given a cube subset x, y (1-D) arrays,
# calculate the projected extent of the corners of the corner pixels
# Input:
#   x, y : 1-D array of x, y locations of centers of cells
# Returns:
#   extent = 4 element array of [ULx, LRx, ULy, LRy]
def get_extent_xy(x, y):
    # Assumes pixels with regular spacing
    scale_x = x[1] - x[0]
    scale_y = y[1] - y[0]
    print("scales: ", scale_x, scale_y)
    
    extent = [x[0] - (scale_x/2.),
              x[-1] + (scale_x/2.),
              y[-1] + (scale_y/2.),
              y[0] - (scale_y/2.)]
              
    return extent
