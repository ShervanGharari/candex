# section 1 load all the necessary modules and packages
import pandas as pd
import numpy as np
import re as re
import datetime
import matplotlib.pyplot as plt
import glob
import os
import netCDF4 as nc4
import time
import xarray as xr
import fiona
import shapely
import geopandas
import shutil
from shapely.geometry import Polygon
from collections import OrderedDict

def lat_lon_2D (lat, lon):
    """
    @ author:                  Shervan Gharari
    @ Github:                  ./shervangharari/repository
    @ author's email id:       sh.gharari@gmail.com
    @ lisence:                 GPL3
    
    This function gets lat and lon in one dimension and returns a 2D matrix of that lat and lon
    input for creating shapefile
    input:
        lat: the lat values with dimension of [n,]
        lon: the lat values with dimension of [m,]
    output:
        lat_2D: the 2D matrix of lat_2D [n,m,]
        lon_2D: the 2D matrix of lon_2D [n,m,]
    """
    # flattening the lat and lon
    lat = lat.flatten()
    lon = lon.flatten()
    
    # having the size of the lat and lon
    idx = lon.shape
    idy = lat.shape
    
    # creating empty numpys for puting the lat and lon values in 2-D format
    lon_2D = np.empty([idx[0],idy[0]])
    lat_2D = np.empty([idx[0],idy[0]])
    
    # for loops for creating the lon_2D and lat_2D
    for i in range(idx[0]):
        for j in range(idy[0]):
            lon_2D[i,j] = lon [i]
            lat_2D[i,j] = lat [j]
    
    # return lat_2D and lon_2D
    return lat_2D, lon_2D

def lat_lon_SHP (lat, lon):
    """
    @ author:                  Shervan Gharari
    @ Github:                  ./shervangharari/repository
    @ author's email id:       sh.gharari@gmail.com
    @license:                  GPL3
    
    This function gets a 2-D lat and lon and return the shapefile given the lat and lon matrices
    input:
        lat: the 2D matrix of lat_2D [n,m,]
        lon: the 2D matrix of lon_2D [n,m,]
    output:
        result: a shapefile with (n-2)*(m-2) elements depicting the provided 2-D lat and lon values
    """
    # getting the shape of the lat and lon (assuming that they have the same shape [n,m,])
    idx = lat.shape
    # preparing an empty shapefile
    result = geopandas.GeoDataFrame()
    # preparing the m whcih is a couter for the shapefile arbitrary ID
    m = 0.00
    
    # itterating to creat the shapes of the result shapefile
    for i in range(1,idx[0]-1):
        
        for j in range(1,idx[1]-1):
            
            # Creating the lat of the shapefile
            Lat_Up = (lat[i-1,j]+lat[i,j])/2
            Lat_UpRright = (lat[i-1,j]+lat[i-1,j+1]+lat[i,j+1]+lat[i,j])/4
            Lat_Right = (lat[i,j+1]+lat[i,j])/2
            Lat_LowRight = (lat[i,j+1]+lat[i+1,j+1]+lat[i+1,j]+lat[i,j])/4;
            Lat_Low = (lat[i+1,j]+lat[i,j])/2;
            Lat_LowLeft = (lat[i,j-1]+lat[i+1,j-1]+lat[i+1,j]+lat[i,j])/4;
            Lat_Left = (lat[i,j-1]+lat[i,j])/2;
            Lat_UpLeft = (lat[i-1,j-1]+lat[i-1,j]+lat[i,j-1]+lat[i,j])/4;
            
            # Creating the lon of the shapefile
            Lon_Up = (lon[i-1,j]+lon[i,j])/2;
            Lon_UpRright = (lon[i-1,j]+lon[i-1,j+1]+lon[i,j+1]+lon[i,j])/4;
            Lon_Right = (lon[i,j+1]+lon[i,j])/2;
            Lon_LowRight = (lon[i,j+1]+lon[i+1,j+1]+lon[i+1,j]+lon[i,j])/4;
            Lon_Low = (lon[i+1,j]+lon[i,j])/2;
            Lon_LowLeft = (lon[i,j-1]+lon[i+1,j-1]+lon[i+1,j]+lon[i,j])/4;
            Lon_Left = (lon[i,j-1]+lon[i,j])/2;
            Lon_UpLeft = (lon[i-1,j-1]+lon[i-1,j]+lon[i,j-1]+lon[i,j])/4;
            
            # craeting the polygongiven the lat and lon
            polys = Polygon([(Lon_Up,Lat_Up),(Lon_UpRright,Lat_UpRright),(Lon_Right,Lat_Right),(Lon_LowRight,Lat_LowRight),\
                    (Lon_Low,Lat_Low),(Lon_LowLeft,Lat_LowLeft),(Lon_Left,Lat_Left),(Lon_UpLeft,Lat_UpLeft), (Lon_Up,Lat_Up)])
            
            # putting the polygone into the shape file
            result.loc[m, 'geometry'] = polys
            result.loc[m, 'ID'] = m + 1.00 # inserting the couter
            result.loc[m, 'lat'] = lat[i,j] # inserting the lat
            result.loc[m, 'lon'] = lon[i,j] # inserting the lon
            
            # adding one to the couter
            m = m + 1.00
    
    # returning the result
    return result

def NetCDF_SHP_lat_lon (name_of_nc):
    """
    @ author:                  Shervan Gharari
    @ Github:                  ./shervangharari/repository
    @ author's email id:       sh.gharari@gmail.com
    @lisence:                  GPL3
    
    This function reads a nc file with 1D or 2D lat and lon
        name_of_nc: name of the nc file and 
    output:
        result: a shapefile that correspond to the 2D lat and lon
    """
    # open the nc file to read
    ncid = nc4.Dataset(name_of_nc, 'r')
    
    # getting the name of the dimentions for lat and lon (they should be the similar)
    lat_dimensions = ncid.variables['lat'].dimensions
    lon_dimensions = ncid.variables['lon'].dimensions
    # putting the names into a np.array
    lat_dimensions = np.array(lat_dimensions)
    lon_dimensions = np.array(lon_dimensions)
    
    # reading the lat and lon and converting them to np.array
    lat = ncid.variables['lat'][:]
    lat = np.array(lat)
    lon = ncid.variables['lon'][:]
    lon = np.array(lon)
    
    # check if lat and lon are 1 D, if yes then they should be converted to 2D lat and lon
    if len(lat.shape) == 1 and len(lon.shape) == 1:
        lat, lon = lat_lon_2D (lat, lon)
    
    # creating the shapefile
    result = lat_lon_SHP (lat, lon)
    
    return result
        

def interesection_shp (shp_1, shp_2):
    """
    @ author:                  Shervan Gharari
    @ Github:                  ./shervangharari/repository
    @ author's email id:       sh.gharari@gmail.com
    @ license:                 GPL3
    
    This function finds the intersection of two shapefiles and calculates their intersection shapes and area and the
    shared percentages.
    input:
        shp_1: the 2D matrix of lat_2D [n,m,]
        shp_2: the 2D matrix of lon_2D [n,m,]
    output:
        result: a shapefile with (n-2)*(m-2) elements depicting the provided 2-D lat and lon values
    """
    # Calculating the area of every shapefile (both should be in degree or meters)
    column_names = shp_1.columns
    column_names = list(column_names)
    column_names.remove('geometry') # removing the geometry from the column names
    
    # renaming the column with S_1
    for i in range(len(column_names)):
        shp_1 = shp_1.rename(columns={column_names[i]: 'S_1_'+column_names[i]})
    
    column_names = shp_2.columns
    column_names = list(column_names)
    column_names.remove('geometry') # removing the geometry from the colomn names
    
    # renaming the column with S_2
    for i in range(len(column_names)):
        shp_2 = shp_2.rename(columns={column_names[i]: 'S_2_'+column_names[i]})
    
    # Caclulating the area for shp1
    for index, row in shp_1.iterrows():
        shp_1.loc[index, 'S1_AREA'] = shp_1.area[index]
        shp_1.loc[index, 'S1_ID'] = index + 1.00
    
    # Caclulating the area for shp2
    for index, row in shp_2.iterrows():
        shp_2.loc[index, 'S2_AREA'] = shp_2.area[index]
        shp_2.loc[index, 'S2_ID'] = index + 1.00
    
    # making intesection
    result = geopandas.overlay(shp_1, shp_2, how='intersection')
    
    # Caclulating the area for shp2
    for index, row in result.iterrows():
        result.loc[index, 'AREA_INT'] = result.area[index]
        result.loc[index, 'AREA_PER_1'] = result.area[index]/result.loc[index, 'S1_AREA']
        result.loc[index, 'AREA_PER_2'] = result.area[index]/result.loc[index, 'S2_AREA']
    
    # normalizing the intersected area
    # to be done
    
    #return shp_1, shp_2
    return result

def read_value_lat_lon_nc (lat, lon, name_of_nc,name_of_variable,name_of_time_dim):
    """
    @ author:                  Shervan Gharari
    @ Github:                  ./shervangharari/repository
    @ author's email id:       sh.gharari@gmail.com
    @ license:                 GPL3
    
    This function funcitons read different grids and sum them up based on the weight provided to aggregate them over a
    larger area
    IMPORTANT NOTICE: the function assume that the dimensions are lat and lon, it should be changed accordingly
    input:
        lon: lon value [1,]
        lat: lat value [1,]
        name_of_nc: full or part of nc file(s) name including nc, string, example 'XXX/*01*.nc'
        name_of_variable: name of the varibale, string
        name_of_time_dim: name of time dimension, string
    output:
        agg_data: a 1-D array of data
    """
    names_all = glob.glob(name_of_nc)
    names_all.sort
    print(names_all)
    data = []
    
    # for to read on lat and lon given

    for names in names_all:
        # open nc file to read
        ncid = nc4.Dataset(names, 'r')
        # getting the length of time dimnesion
        time_steps = ncid.dimensions[name_of_time_dim].size
        print(time_steps)
        # opening the data set by xr package
        da = xr.open_dataset(names,decode_times=False)
        # read the data given a lat, lon, and varibale names, lat and lon should be called lat and lon in the file
        data_temp = da[name_of_variable].sel(latitude=lat, longitude=lon, method='nearest')
        print(data_temp)
        # put the read data into the data_temp
        data_temp = np.array(data_temp)
        data_temp = data_temp.reshape((time_steps,))
        print(data_temp)
        # append the data_temp
        if 'data' in globals():
            data = np.append(data, data_temp)
        else:
            data = data_temp
    return data
            
def area_ave (lat, lon, w, name_of_nc,name_of_variable,name_of_time_dim):
    """
    @ author:                  Shervan Gharari
    @ Github:                  ./shervangharari/repository
    @ author's email id:       sh.gharari@gmail.com
    @ license:                 GPL3
    
    This function funcitons read different grids and sum them up based on the weight provided to aggregate them over a
    larger area
    IMPORTANT NOTICE: the function assume that the dimensions are lat and lon, it should be changed accordingly
    input:
        lon: lon value [n,]
        lat: lat value [n,]
        w: weight [n,]
        name_of_nc: full or part of nc file(s) name including nc, string, example 'XXX/*01*.nc'
        name_of_variable: name of the varibale, string
        name_of_time_dim: name of time dimension, string
    output:
        agg_data: a 1-D array of data
    """
    for i in range(0, len(lat)):
        data_temp = read_value_lat_lon_nc (lat[i], lon[i], name_of_nc,name_of_variable,name_of_time_dim)
        if i is 0:
            data = data_temp*w[i]
        else:
            data = data + data_temp*w[i]
    return data
    
    
