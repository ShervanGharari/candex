# section 1 load all the necessary modules and packages
import glob
import time
import geopandas as gpd
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

    This function gets lat and lon in one-dimension and returns a two-dimensional matrix of that lat and lon
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


def lat_lon_SHP(lat, lon, box, correct_360):
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
    result = gpd.GeoDataFrame()
    # preparing the m whcih is a couter for the shapefile arbitrary ID
    m = 0.00
    # making sure that the lon is less than 180
    if correct_360 is True:
        IN = lon>180 # index of more than 180
        lon[IN] = lon[IN]-360 # index of point with higher than are reduced to -180 to 0 instead
        print(lon)

    # itterating to create the shapes of the result shapefile
    for i in range(1, idx[0] - 1):
        for j in range(1, idx[1] - 1):
            if lat[i, j] > box[0] and lat[i, j] < box[1] and lon[i, j] > box[
                    2] and lon[i, j] < box[3]: # checking is lat and lon is located inside the provided box
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
                polys = Polygon([(Lon_Up,Lat_Up),(Lon_UpRright,Lat_UpRright),\
                                 (Lon_Right,Lat_Right),(Lon_LowRight,Lat_LowRight),\
                                 (Lon_Low,Lat_Low),(Lon_LowLeft,Lat_LowLeft),\
                                 (Lon_Left,Lat_Left),(Lon_UpLeft,Lat_UpLeft),\
                                 (Lon_Up,Lat_Up)])

                # putting the polygone into the shape file
                result.loc[m, 'geometry'] = polys
                result.loc[m, 'ID'] = m + 1.00  # inserting the couter, ID
                result.loc[m, 'lat'] = lat[i, j]  # inserting the lat
                if correct_360 is True:
                    result.loc[m, 'lon'] = lon[i, j]+360  # inserting the lon
                if correct_360 is False:
                    result.loc[m, 'lon'] = lon[i, j]  # inserting the lon
                    

                # adding one to the couter
                m = m + 1.00

    # returning the result
    return result


def NetCDF_SHP_lat_lon(name_of_nc, box, name_of_lat_var, name_of_lon_var, correct_360):
    """
    @ author:                  Shervan Gharari
    @ Github:                  ./shervangharari/repository
    @ author's email id:       sh.gharari@gmail.com
    @lisence:                  Apache2

    This function reads a nc file with 1D or 2D lat and lon
        name_of_nc: name of the nc file
        box:the defined box for making the shp file
        name_of_lat_var:
        name_of_lon_var:
    output:
        result: a shapefile that correspond to the 2D lat and lon
    """
    # open the nc file to read
    dataset = xr.open_dataset(name_of_nc, decode_times=False)

    # reading the lat and lon and converting them to np.array
    lat = dataset[name_of_lat_var].data
    lon = dataset[name_of_lon_var].data
    
    lat = np.array(lat)
    lon = np.array(lon)

    # check if lat and lon are 1 D, if yes then they should be converted to 2D lat and lon WARNING only for case 1 and 2
    if len(lat.shape) == 1 and len(lon.shape) == 1:
        lat, lon = lat_lon_2D(lat, lon)

    # creating the shapefile
    result = lat_lon_SHP(lat, lon, box, correct_360)

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
    result = spatial_overlays (shp_1, shp_2, how='intersection')
    # result = geopandas.tools.overlay(shp_1, shp_2, how='intersection')
    # result = geopandas.overlay(shp_1, shp_2, how='intersection')

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

def read_value_lat_lon_nc(case,
                          lat_target, lon_target, name_of_nc,
                          name_of_variable, name_of_time_dim,
                          name_of_lat_var, name_of_lon_var,
                          name_of_lat_dim, name_of_lon_dim):
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
        case: value [1,]
            1 is for 3-dimensional variable with 1-dimentional lat and lon
            2 is for 3-dimensional varibale with 2-dimentional lat and lon
            3 is for 2-dimensional variable with 1-dimentional lat and lon (time series)
        lon_target: lon value [1,]
        lat_target: lat value [1,]
        name_of_nc: full or part of nc file(s) name including nc, string, example 'XXX/*01*.nc'
        name_of_variable: name of the varibale, string
        name_of_time_dim: name of time dimension, string
        name_of_lat_dim: name of lat dimension, string
        name_of_lon_dim: name of lon dimension, string
        name_of_lat_var: name of lat variable, string
        name_of_lon_var: name of lon variable, string

    output:
        agg_data: a 1-D array of data
    """
    names_all = glob.glob(name_of_nc)
    names_all.sort()
    data = None

    # for to read on variouse nc files for target lat lon
    for names in names_all:
        
        da = xr.open_dataset(names, decode_times=False)
        
        # case 1, the varibale is 3-dimensional and lat and lon are one-dimnesional
        if case ==1:
            # finding the index for the lat for target_lat
            da_lat = da[name_of_lat_var] # reading the lat variable
            temp = np.array(abs(da_lat-lat_target)) # finding the distance to target_lat
            index_target_lat = np.array([temp.argmin()]) # finding the closest index to target_lat
            
            # finding the index for the lon for target_lon
            da_lon = da[name_of_lon_var] # reading the lon variable
            temp = np.array(abs(da_lon-lon_target)) # finding the distnace to target_lon
            index_target_lon = np.array([temp.argmin()]) # finding the closest index to target_lon
            
            # making sure that the lat and lon are only one value and not two
            index_target_lon = index_target_lon[0]
            index_target_lat = index_target_lat[0]
            
            # porder of dimensions for the target variable
            dataset = da[name_of_variable]
            order_time_dim = dataset.dims.index(name_of_time_dim)
            order_lat_dim = dataset.dims.index(name_of_lat_dim)
            order_lon_dim = dataset.dims.index(name_of_lon_dim)
            
            if order_time_dim == 0: # such as the varibaele dimension is time, lat, lon
                data_temp = dataset [:,index_target_lat,index_target_lon]
            if order_time_dim == 2: # such as the varibaele dimension is lon, lat, time
                data_temp = dataset [index_target_lon,index_target_lat,:]

        # case 2, the varibale is 3-dimensional and lat and lon are 2-dimentional such as rotated lat lon
        if case ==2:
            # finding the index for the lat
            da_lat = da[name_of_lat_var]
            da_lon = da[name_of_lon_var]
            temp = np.array(abs(da_lat-lat_target)+abs(da_lon-lon_target))
            ind = np.unravel_index(np.argmin(temp, axis=None), temp.shape)
            ind = np.array(ind)
            
            # order of dimensions for the target variable
            dataset = da[name_of_variable]
            order_time_dim = dataset.dims.index(name_of_time_dim)
            order_lat_dim = dataset.dims.index(name_of_lat_dim)
            order_lon_dim = dataset.dims.index(name_of_lon_dim)
            
            if order_time_dim == 0: # such as the varibaele dimension is time, lat, lon
                index_target_lat = ind[0]
                index_target_lon = ind[1]
                data_temp = dataset [:,index_target_lat,index_target_lon]
            if order_time_dim == 2: # such as the varibaele dimension is lon, lat, time
                index_target_lat = ind[1]
                index_target_lon = ind[0]
                data_temp = dataset [index_target_lon,index_target_lat,:]

        # case 3, the varibale is 2-dimnesional and lat and lon are 1-dimensional such as n time or time n
        if case ==3:
            da_lat = da[name_of_lat_var]
            da_lon = da[name_of_lon_var]
            temp = np.array(abs(da_lat-lat_target)+abs(da_lon-lon_target))
            ind = np.unravel_index(np.argmin(temp, axis=None), temp.shape)
            ind = np.array(ind)
            ind = ind[0]
            
            # order of dimensions for the target variable
            dataset = da[name_of_variable]
            order_time_dim = dataset.dims.index(name_of_time_dim)
            
            if order_time_dim == 0: # such as the varibaele dimension is time, n
                index_target_n = ind
                data_temp = dataset [:,index_target_n]
            if order_time_dim == 1: # such as the varibaele dimension is n, time
                index_target_n = ind
                data_temp = dataset [index_target_n,:]
        
        # getting the length of time dimension
        time_steps = da.dims[name_of_time_dim]
        
        # put the read data into the data_temp
        data_temp = np.array(data_temp)
        data_temp = data_temp.reshape((time_steps, ))

        # append the data_temp
        if data is not None:
            data = np.append(data, data_temp)
        else:
            data = data_temp
            
    return data


def area_ave(case,
             lat, lon, w,
             name_of_nc, name_of_variable,
             name_of_time_dim,
             name_of_lat_dim, name_of_lon_dim,
             name_of_lat_var, name_of_lon_var):
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
        lon: lon value [n,] target value to read from the nc file
        lat: lat value [n,] target values to read from the nc file
        w: weight [n,] weigth of each target values to read from the nc file
        name_of_nc: full or part of nc file(s) name including nc, string, example 'XXX/*01*.nc'
        name_of_variable: name of the varibale, string
        name_of_time_dim: name of time dimension, string
        name_of_lat_dim: name of lat dimension, string
        name_of_lon_dim: name of lon dimension, string
        name_of_lat_var: name of lat variable, string
        name_of_lon_var: name of lon variable, string
    output:
        data: a 1-D array of data
    """
    for i in range(0, len(lat)): # itterate over target values
        data_temp = read_value_lat_lon_nc(case,
                                          lat.iloc[i], lon.iloc[i], name_of_nc,
                                          name_of_variable, name_of_time_dim,
                                          name_of_lat_dim, name_of_lon_dim,
                                          name_of_lat_var, name_of_lon_var)
        if i is 0: # multiply the read value with their weight and sum
            data = data_temp * w.iloc[i]
        else:
            data = data + data_temp * w.iloc[i]
    return data


def box(name_or_singleframe_shp):
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
    if type(name_or_singleframe_shp) is str:
        shp = gpd.read_file(name_or_singleframe_shp)
    else:
        shp = name_or_singleframe_shp
    # checking if dataframe has only one row otherwise stop
    if shp.shape[0] > 1:
    #    raise Exception('the dataframe in the box function has more than one row')
        print('WARNING: your shapefile has more than one value! in box funcitons')
    shp['geometry'] = shp.geometry.buffer(1)  # adding a buffer of 1 degree
    A = shp.bounds
    box = np.array([A.miny.min(), A.maxy.max(), A.minx.min(), A.maxx.max()])
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
        data_varid = ncid.createVariable(variable_name, 'f8', ('n','time', ))
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

        
        
        
def spatial_overlays(df1, df2, how='intersection', reproject=True):
    """Perform spatial overlay between two polygons.

    Currently only supports data GeoDataFrames with polygons.
    Implements several methods that are all effectively subsets of
    the union.
    
    Ömer Özak
    ozak
    https://github.com/ozak
    https://github.com/geopandas/geopandas/pull/338

    Parameters
    ----------
    df1 : GeoDataFrame with MultiPolygon or Polygon geometry column
    df2 : GeoDataFrame with MultiPolygon or Polygon geometry column
    how : string
        Method of spatial overlay: 'intersection', 'union',
        'identity', 'symmetric_difference' or 'difference'.
    use_sindex : boolean, default True
        Use the spatial index to speed up operation if available.

    Returns
    -------
    df : GeoDataFrame
        GeoDataFrame with new set of polygons and attributes
        resulting from the overlay

    """
    df1 = df1.copy()
    df2 = df2.copy()
    df1['geometry'] = df1.geometry.buffer(0)
    df2['geometry'] = df2.geometry.buffer(0)
    if df1.crs!=df2.crs and reproject:
        print('Data has different projections.')
        print('Converted data to projection of first GeoPandas DatFrame')
        df2.to_crs(crs=df1.crs, inplace=True)
    if how=='intersection':
        # Spatial Index to create intersections
        spatial_index = df2.sindex
        df1['bbox'] = df1.geometry.apply(lambda x: x.bounds)
        df1['sidx']=df1.bbox.apply(lambda x:list(spatial_index.intersection(x)))
        pairs = df1['sidx'].to_dict()
        nei = []
        for i,j in pairs.items():
            for k in j:
                nei.append([i,k])
        pairs = gpd.GeoDataFrame(nei, columns=['idx1','idx2'], crs=df1.crs)
        pairs = pairs.merge(df1, left_on='idx1', right_index=True)
        pairs = pairs.merge(df2, left_on='idx2', right_index=True, suffixes=['_1','_2'])
        pairs['Intersection'] = pairs.apply(lambda x: (x['geometry_1'].intersection(x['geometry_2'])).buffer(0), axis=1)
        pairs = gpd.GeoDataFrame(pairs, columns=pairs.columns, crs=df1.crs)
        cols = pairs.columns.tolist()
        cols.remove('geometry_1')
        cols.remove('geometry_2')
        cols.remove('sidx')
        cols.remove('bbox')
        cols.remove('Intersection')
        dfinter = pairs[cols+['Intersection']].copy()
        dfinter.rename(columns={'Intersection':'geometry'}, inplace=True)
        dfinter = gpd.GeoDataFrame(dfinter, columns=dfinter.columns, crs=pairs.crs)
        dfinter = dfinter.loc[dfinter.geometry.is_empty==False]
        dfinter.drop(['idx1','idx2'], inplace=True, axis=1)
        return dfinter
    elif how=='difference':
        spatial_index = df2.sindex
        df1['bbox'] = df1.geometry.apply(lambda x: x.bounds)
        df1['sidx']=df1.bbox.apply(lambda x:list(spatial_index.intersection(x)))
        df1['new_g'] = df1.apply(lambda x: reduce(lambda x, y: x.difference(y).buffer(0), 
                                 [x.geometry]+list(df2.iloc[x.sidx].geometry)) , axis=1)
        df1.geometry = df1.new_g
        df1 = df1.loc[df1.geometry.is_empty==False].copy()
        df1.drop(['bbox', 'sidx', 'new_g'], axis=1, inplace=True)
        return df1
    elif how=='symmetric_difference':
        df1['idx1'] = df1.index.tolist()
        df2['idx2'] = df2.index.tolist()
        df1['idx2'] = np.nan
        df2['idx1'] = np.nan
        dfsym = df1.merge(df2, on=['idx1','idx2'], how='outer', suffixes=['_1','_2'])
        dfsym['geometry'] = dfsym.geometry_1
        dfsym.loc[dfsym.geometry_2.isnull()==False, 'geometry'] = dfsym.loc[dfsym.geometry_2.isnull()==False, 'geometry_2']
        dfsym.drop(['geometry_1', 'geometry_2'], axis=1, inplace=True)
        dfsym = gpd.GeoDataFrame(dfsym, columns=dfsym.columns, crs=df1.crs)
        spatial_index = dfsym.sindex
        dfsym['bbox'] = dfsym.geometry.apply(lambda x: x.bounds)
        dfsym['sidx'] = dfsym.bbox.apply(lambda x:list(spatial_index.intersection(x)))
        dfsym['idx'] = dfsym.index.values
        dfsym.apply(lambda x: x.sidx.remove(x.idx), axis=1)
        dfsym['new_g'] = dfsym.apply(lambda x: reduce(lambda x, y: x.difference(y).buffer(0), 
                         [x.geometry]+list(dfsym.iloc[x.sidx].geometry)) , axis=1)
        dfsym.geometry = dfsym.new_g
        dfsym = dfsym.loc[dfsym.geometry.is_empty==False].copy()
        dfsym.drop(['bbox', 'sidx', 'idx', 'idx1','idx2', 'new_g'], axis=1, inplace=True)
        return dfsym
    elif how=='union':
        dfinter = spatial_overlays(df1, df2, how='intersection')
        dfsym = spatial_overlays(df1, df2, how='symmetric_difference')
        dfunion = dfinter.append(dfsym)
        dfunion.reset_index(inplace=True, drop=True)
        return dfunion
    elif how=='identity':
        dfunion = spatial_overlays(df1, df2, how='union')
        cols1 = df1.columns.tolist()
        cols2 = df2.columns.tolist()
        cols1.remove('geometry')
        cols2.remove('geometry')
        cols2 = set(cols2).intersection(set(cols1))
        cols1 = list(set(cols1).difference(set(cols2)))
        cols2 = [col+'_1' for col in cols2]
        dfunion = dfunion[(dfunion[cols1+cols2].isnull()==False).values]
        return dfunion
