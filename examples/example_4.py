"""
@ author:                  Shervan Gharari
@ Github:                  https://github.com/ShervanGharari/candex
@ author's email id:       sh.gharari@gmail.com
@ license:                 Apache2
"""

# example for many shapefile 4
# the name of the NetCDF file to be read, all the file for 2010
nc_name = 'local_dir/*2010*.nc' # the local directory should be copied and pasted here.


shp_name = 'local_dir/4_LEB_boundary_subwatershed.shp' # the name of the shapefile with many shape
shp_1 = gpd.read_file(shp_name) # reading the shapefile as a geopandas data frame
shp_1.crs = {'init': 'epsg:4326'} # setting the cordinate system in case it doesnt exists

box_value = box (shp_1,1) # buffer of 1 degree

shp_2 = NetCDF_SHP_lat_lon('local_dir/Rainf_daily_WFDEI_CRU_201612.nc',box_value,'lat','lon',False) # creating a mesh of NetCDF
shp_2.crs = {'init': 'epsg:4326'} # setting the cordinate system

# intersection of the shapefiles
shp_int = intersection_shp (shp_1, shp_2)

# finding the unique shapes that are in intersection from the shp_1, their field in the intersection is called S_1_*, in this example
# S_1_OBJECTID
IDs_from_int = np.unique(shp_int['S_1_OBJECTID'])

data_all = None # creating empty field of the data_all

for i in IDs_from_int: # looping over values in ID
    
    # extracting the part of intersection with unique S_1_OBJECTID values
    shp_temp = shp_int.loc[shp_int['S_1_OBJECTID'] == i]
    
    # getting the lat/lon, weight into the numpy array
    lat = np.array(shp_temp.S_2_lat)
    lon = np.array(shp_temp.S_2_lon)
    W   = np.array(shp_temp.AP1N)

    case = 1 # case 1 for regualre lat lon
    # reading the data
    data = area_ave (case, lat, lon, W,\
                     nc_name,'Tair','tstep','lat','lon','lat','lon')
    # storign the data for every shape from shp_1
    if data_all is not None:
        data_all = np.vstack((data_all, data))
    else:
        data_all = data

# visualziation
plt.plot(data_all.transpose()-273.3)
