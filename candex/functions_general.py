# section 1 load all the necessary modules and packages
import glob
import time
import geopandas
import netCDF4 as nc4
import numpy as np
import pandas as pd
import xarray as xr
from shapely.geometry import Polygon


def lat_lon_2D(lat, lon):
    """
    @ author:                  Shervan Gharari
    @ Github:                  ./shervangharari/repository
    @ author's email id:       sh.gharari@gmail.com
    @ license:                 Apache2

    This function gets lat and lon in one dimension and returns a 2D matrix of that lat and lon
    input for creating shapefile

    Arguments
    ---------
    lat : numpy.ndarray
        the lat values with dimension of [n,]
    lon : numpy.ndarray
        the lat values with dimension of [m,]

    Returns
    -------
    tuple[numpy.ndarray]
        lat_2D: the 2D matrix of lat_2D [n,m,]
        lon_2D: the 2D matrix of lon_2D [n,m,]
    """
    # flattening the lat and lon
    lat = lat.flatten()
    lon = lon.flatten()

    # return lat_2D and lon_2D
    return np.meshgrid(lat, lon)


def lat_lon_SHP(lat, lon, box):
    """
    @ author:                  Shervan Gharari
    @ Github:                  ./shervangharari/repository
    @ author's email id:       sh.gharari@gmail.com
    @license:                  Apache2

    This function gets a 2-D lat and lon and return the shapefile given the lat and lon matrices
    input:
        lat: the 2D matrix of lat_2D [n,m,]
        lon: the 2D matrix of lon_2D [n,m,]
        box: a 1D array [minlat, maxlat, minlon, maxlon]
    output:
        result: a shapefile with (n-2)*(m-2) elements depicting the provided 2-D lat and lon values
    """
    # getting the shape of the lat and lon (assuming that they have the same shape [n,m,])
    idx = lat.shape
    # preparing an empty shapefile
    result = geopandas.GeoDataFrame()
    # preparing the m whcih is a couter for the shapefile arbitrary ID
    m = 0.00

    # itterating to create the shapes of the result shapefile
    for i in range(1, idx[0] - 1):

        for j in range(1, idx[1] - 1):

            if lat[i, j] > box[0] and lat[i, j] < box[1] and lon[i, j] > box[
                    2] and lon[i, j] < box[3]:

                # Creating the lat of the shapefile
                Lat_Up = (lat[i - 1, j] + lat[i, j]) / 2
                Lat_UpRright = (lat[i - 1, j] + lat[i - 1, j + 1] +
                                lat[i, j + 1] + lat[i, j]) / 4
                Lat_Right = (lat[i, j + 1] + lat[i, j]) / 2
                Lat_LowRight = (lat[i, j + 1] + lat[i + 1, j + 1] +
                                lat[i + 1, j] + lat[i, j]) / 4
                Lat_Low = (lat[i + 1, j] + lat[i, j]) / 2
                Lat_LowLeft = (lat[i, j - 1] + lat[i + 1, j - 1] +
                               lat[i + 1, j] + lat[i, j]) / 4
                Lat_Left = (lat[i, j - 1] + lat[i, j]) / 2
                Lat_UpLeft = (lat[i - 1, j - 1] + lat[i - 1, j] + lat[i, j - 1]
                              + lat[i, j]) / 4

                # Creating the lon of the shapefile
                Lon_Up = (lon[i - 1, j] + lon[i, j]) / 2
                Lon_UpRright = (lon[i - 1, j] + lon[i - 1, j + 1] +
                                lon[i, j + 1] + lon[i, j]) / 4
                Lon_Right = (lon[i, j + 1] + lon[i, j]) / 2
                Lon_LowRight = (lon[i, j + 1] + lon[i + 1, j + 1] +
                                lon[i + 1, j] + lon[i, j]) / 4
                Lon_Low = (lon[i + 1, j] + lon[i, j]) / 2
                Lon_LowLeft = (lon[i, j - 1] + lon[i + 1, j - 1] +
                               lon[i + 1, j] + lon[i, j]) / 4
                Lon_Left = (lon[i, j - 1] + lon[i, j]) / 2
                Lon_UpLeft = (lon[i - 1, j - 1] + lon[i - 1, j] + lon[i, j - 1]
                              + lon[i, j]) / 4

                # craeting the polygon given the lat and lon
                polys = Polygon([(Lon_Up,Lat_Up),(Lon_UpRright,Lat_UpRright),(Lon_Right,Lat_Right),(Lon_LowRight,Lat_LowRight),\
                        (Lon_Low,Lat_Low),(Lon_LowLeft,Lat_LowLeft),(Lon_Left,Lat_Left),(Lon_UpLeft,Lat_UpLeft), (Lon_Up,Lat_Up)])

                # putting the polygone into the shape file
                result.loc[m, 'geometry'] = polys
                result.loc[m, 'ID'] = m + 1.00  # inserting the couter
                result.loc[m, 'lat'] = lat[i, j]  # inserting the lat
                result.loc[m, 'lon'] = lon[i, j]  # inserting the lon

                # adding one to the couter
                m = m + 1.00

    # returning the result
    return result


def NetCDF_SHP_lat_lon(name_of_nc, box, name_of_lat, name_of_lon):
    """
    @ author:                  Shervan Gharari
    @ Github:                  ./shervangharari/repository
    @ author's email id:       sh.gharari@gmail.com
    @lisence:                  Apache2

    This function reads a nc file with 1D or 2D lat and lon
        name_of_nc: name of the nc file and
    output:
        result: a shapefile that correspond to the 2D lat and lon
    """
    # open the nc file to read
    dataset = xr.open_dataset(name_of_nc, decode_times=False)

    # reading the lat and lon and converting them to np.array
    lat = dataset[name_of_lat].data
    lon = dataset[name_of_lon].data

    # check if lat and lon are 1 D, if yes then they should be converted to 2D lat and lon
    if len(lat.shape) == 1 and len(lon.shape) == 1:
        lat, lon = lat_lon_2D(lat, lon)

    # creating the shapefile
    result = lat_lon_SHP(lat, lon, box)

    return result


def intersection_shp(shp_1, shp_2):
    """
    @ author:                  Shervan Gharari
    @ Github:                  ./shervangharari/repository
    @ author's email id:       sh.gharari@gmail.com
    @ license:                 Apache2

    This function finds the intersection of two shapefiles and calculates their intersection shapes and area and the
    shared percentages.
    input:
        shp_1: the first shapefile
        shp_2: the second shapefile
    output:
        result: a shapefile which is the intersect of shp_1 shp_2
    """
    # Calculating the area of every shapefile (both should be in degree or meters)
    column_names = shp_1.columns
    column_names = list(column_names)

    # removing the geometry from the column names
    column_names.remove('geometry')

    # renaming the column with S_1
    for i in range(len(column_names)):
        shp_1 = shp_1.rename(
            columns={column_names[i]: 'S_1_' + column_names[i]})

    column_names = shp_2.columns
    column_names = list(column_names)

    # removing the geometry from the colomn names
    column_names.remove('geometry')

    # renaming the column with S_2
    for i in range(len(column_names)):
        shp_2 = shp_2.rename(
            columns={column_names[i]: 'S_2_' + column_names[i]})

    # Caclulating the area for shp1
    for index, _ in shp_1.iterrows():
        shp_1.loc[index, 'S1_AREA'] = shp_1.area[index]
        shp_1.loc[index, 'S1_ID'] = index + 1.00

    # Caclulating the area for shp2
    for index, _ in shp_2.iterrows():
        shp_2.loc[index, 'S2_AREA'] = shp_2.area[index]
        shp_2.loc[index, 'S2_ID'] = index + 1.00

    # making intesection
    result = geopandas.overlay(shp_1, shp_2, how='intersection')

    # Caclulating the area for shp2
    for index, _ in result.iterrows():
        result.loc[index, 'AREA_INT'] = result.area[index]
        result.loc[index, 'AREA_PER_1'] = result.area[index] / result.loc[
            index, 'S1_AREA']
        result.loc[index, 'AREA_PER_2'] = result.area[index] / result.loc[
            index, 'S2_AREA']

    # normalizing the intersected area
    # to be done

    #return shp_1, shp_2
    return result


def read_value_lat_lon_nc(lat, lon, name_of_nc, name_of_variable,
                          name_of_time_dim, name_of_lat, name_of_lon):
    """
    @ author:                  Shervan Gharari
    @ Github:                  ./shervangharari/repository
    @ author's email id:       sh.gharari@gmail.com
    @ license:                 Apache2

    This function funcitons read different grids and sum them up based on the
    weight provided to aggregate them over a larger area

    IMPORTANT NOTICE: the function assume that the dimensions are lat and lon,
    it should be changed accordingly

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
    names_all.sort()
    data = None

    select_dictionary = {name_of_lat: lat, name_of_lon: lon}

    # for to read on lat and lon given
    for names in names_all:
        # open nc file to read
        # opening the data set by xr package
        da = xr.open_dataset(names, decode_times=False)

        # getting the length of time dimension
        time_steps = da.dims[name_of_time_dim]

        # read the data given a lat, lon, and variable names,
        data_temp = da[name_of_variable].sel(
            select_dictionary, method='nearest')

        # put the read data into the data_temp
        data_temp = np.array(data_temp)
        data_temp = data_temp.reshape((time_steps, ))

        # append the data_temp
        if data is not None:
            data = np.append(data, data_temp)
        else:
            data = data_temp
    return data


def area_ave(lat, lon, w, name_of_nc, name_of_variable, name_of_time_dim,
             name_of_lat, name_of_lon):
    """
    @ author:                  Shervan Gharari
    @ Github:                  ./shervangharari/repository
    @ author's email id:       sh.gharari@gmail.com
    @ license:                 Apache2

    This function funcitons read different grids and sum them up based on the
    weight provided to aggregate them over a larger area
    IMPORTANT NOTICE: the function assume that the dimensions are lat and lon,
    it should be changed accordingly

    input:
        lon: lon value [n,]
        lat: lat value [n,]
        w: weight [n,]
        name_of_nc: full or part of nc file(s) name including nc, string,
            example 'XXX/*01*.nc'
        name_of_variable: name of the varibale, string
        name_of_time_dim: name of time dimension, string
    output:
        agg_data: a 1-D array of data
    """
    for i in range(0, len(lat)):
        data_temp = read_value_lat_lon_nc(lat.iloc[i], lon.iloc[i], name_of_nc,
                                          name_of_variable, name_of_time_dim,
                                          name_of_lat, name_of_lon)
        if i is 0:
            data = data_temp * w.iloc[i]
        else:
            data = data + data_temp * w.iloc[i]
    return data


def box(name_of_shp):
    """
    @ author:                  Shervan Gharari
    @ Github:                  ./shervangharari/repository
    @ author's email id:       sh.gharari@gmail.com
    @ license:                 Apache2

    This function funcitons calculates the bounding box of a single shapefile
    input:
        name_of_shp: full or part of shp file(s) name including .shp, string,
        example 'XXX/*01*.shp'
    output:
        box: which return [minlat maxlat minlon maxlon]
    """
    shp = geopandas.read_file(name_of_shp)
    shp['geometry'] = shp.geometry.buffer(1)  # adding a buffer of 1 degree
    A = shp.bounds
    box = np.array([A.miny, A.maxy, A.minx, A.maxx])
    return box


def write_netcdf(nc_file_name, variable_data, variable_name, varibale_unit,
                 varibale_long_name, lon_data, lat_data, ID_data,
                 variable_time, starting_date_string):
    """
    @ author:                  Shervan Gharari
    @ Github:                  ./shervangharari/repository
    @ author's email id:       sh.gharari@gmail.com
    @ license:                 Apache2

    This function takes in a single array of data with an ID and it lat and lon value and save it as nc file
    input:
        nc_file_name: the name of the nc file, string
        variable_data: and array of varibale data [time,], with dimension time
        variable_name: name of the varibale, string
        varibale_unit: unit of the varibale, string
        varibale_long_name: long name of the varibale,string, e.g. 'm'
        lon_data: lon of the data
        lat_data: lat of the data
        ID_data: ID of the data such as basin ID or grid ID
        variable_time: data varibale time
        starting_date_string: the time reseference, string, 'days since 1900-01-01 00:00'
    """

    with nc4.Dataset(nc_file_name, "w", format="NETCDF4") as ncid:

        dimid_N = ncid.createDimension('n', 1)  # only write one variable
        dimid_T = ncid.createDimension('time', variable_time.size)

        # Variables
        time_varid = ncid.createVariable('time', 'i4', ('time', ))
        # Attributes
        time_varid.long_name = 'time'
        time_varid.units = starting_date_string  # e.g. 'days since 1900-01-01 00:00'
        time_varid.calendar = 'gregorian'
        time_varid.standard_name = 'time'
        time_varid.axis = 'T'
        # Write data
        time_varid[:] = variable_time

        # Variables
        lat_varid = ncid.createVariable('lat', 'f8', ('n', ))
        lon_varid = ncid.createVariable('lon', 'f8', ('n', ))
        ID_varid = ncid.createVariable('ID', 'f8', ('n', ))
        # Attributes
        lat_varid.long_name = 'latitude'
        lon_varid.long_name = 'longitude'
        ID_varid.long_name = 'ID'
        lat_varid.units = 'degrees_north'
        lon_varid.units = 'degrees_east'
        ID_varid.units = '1'
        lat_varid.standard_name = 'latitude'
        lon_varid.standard_name = 'longitude'
        # Write data
        lat_varid[:] = lat_data
        lon_varid[:] = lon_data
        ID_varid[:] = ID_data

        # Variable
        data_varid = ncid.createVariable(variable_name, 'f8', (
            'n',
            'time',
        ))
        # Attributes
        data_varid.long_name = varibale_long_name
        data_varid.units = varibale_unit
        # Write data
        data_varid[:] = variable_data

        ##
        ncid.Conventions = 'CF-1.6'
        ncid.License = 'The data were written by Shervan Gharari. Under Apache2.'
        ncid.history = 'Created ' + time.ctime(time.time())
        ncid.source = 'Written by script from library of Shervan Gharari (https://github.com/ShervanGharari).'
