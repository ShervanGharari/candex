"""
@ author:                  Shervan Gharari
@ Github:                  ./shervangharari/candex
@ author's email id:       sh.gharari@gmail.com
@ license:                 Apache2
This is an example to read values from a NetCDF file for a shapefile that constains many shapes
"""

# example for many shapefile 4
nc_name = 'local_dir/*2010*.nc' # the local directory should be copied and pasted here.


shp_name = 'local_dir/4_LEB_boundary_subwatershed.shp'
shp_1 = gpd.read_file(shp_name)
shp_1.crs = {'init': 'epsg:4326'} # setting the cordinate system in case it doesnt exists

box_value = box (shp_1)
shp_2 = NetCDF_SHP_lat_lon('local_dir/Rainf_daily_WFDEI_CRU_201612.nc',box_value,'lat','lon',False)
shp_2.crs = {'init': 'epsg:4326'} # setting the cordinate system in case it doesnt exists

shp_int = intersection_shp (shp_1, shp_2)

IDs_from_int = np.unique(shp_int['S_1_OBJECTID'])

data_all = None

for i in IDs_from_int:
    
    shp_temp = shp_int.loc[shp_int['S_1_OBJECTID'] == i]
    
    lat = np.array(shp_temp.S_2_lat)
    lon = np.array(shp_temp.S_2_lon)
    W   = np.array(shp_temp.AP1N)

    case = 1
    data = area_ave (case, lat, lon, W,\
                     nc_name,'Tair','tstep','lat','lon','lat','lon')
    if data_all is not None:
        data_all = np.vstack((data_all, data))
    else:
        data_all = data
