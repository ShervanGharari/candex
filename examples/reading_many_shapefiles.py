# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 02:17:44 2019

@author: shg096
"""
import glob
import geopandas
import numpy as np
import pandas as pd
import xarray as xr
from shapely.geometry import Polygon
import netCDF4 as nc4
import time
from candex import box, NetCDF_SHP_lat_lon, intersection_shp, area_ave, write_netcdf

name = 'C:/Users/shg096/Dropbox/Individual_Shape_Meuse/Individual_Shape_Meuse/**/*.shp'
names_all = glob.glob(name)
names_all.sort()
print(names_all)

nc_name = '//datastore/GLOBALWATER/giws_research_water_share/ClimateForcing_Data/ClimateForcing_WFDEI/WFDEI_05d24hr/Tair_daily_WFDEI/*199*.nc'

i = 1
    
for name in names_all:
    shp_1 = geopandas.read_file(name) # reading each shapefile one by one...
    box_values = box (name)
    #assuming the netcdf lat lon are not changing and equals to the first netcdf file
    shp_2 = NetCDF_SHP_lat_lon ('//datastore/GLOBALWATER/giws_research_water_share/ClimateForcing_Data/ClimateForcing_WFDEI/WFDEI_05d24hr/Rainf_daily_WFDEI_CRU/Rainf_daily_WFDEI_CRU_201612.nc',box_values,'lat','lon')
    shp_int = intersection_shp (shp_1, shp_2)
    data = area_ave (shp_int.S_2_lat, shp_int.S_2_lon, shp_int.AREA_PER_1, nc_name,'Tair','tstep','lat','lon')
    
    ## netcdf wrting...
    nc_file_name = 'C:/Users/shg096/Dropbox/test'+str(i)+'.nc'
    variable_data = data
    variable_name = 'T_02'
    varibale_unit = 'K'
    varibale_long_name = 'surface temprature at 2 meters'
    lon_data = shp_1.LONG.iloc[0]
    lat_data = shp_1.LAT.iloc[0]
    ID_data = shp_1.ID.iloc[0]
    variable_time = np.arange(0,data.size)
    starting_date_string_and_units = 'days since 1990-01-01 00:00:00'
    
    write_netcdf (nc_file_name, variable_data, variable_name, varibale_unit, varibale_long_name,
                  lon_data, lat_data, ID_data,
                  variable_time, starting_date_string_and_units)
    
    i = i + 1;