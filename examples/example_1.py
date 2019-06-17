"""
@ author:                  Shervan Gharari
@ Github:                  ./shervangharari/candex
@ author's email id:       sh.gharari@gmail.com
@ license:                 Apache2
"""


nc_name = 'local_dir/*199*.nc' # the local directory should be copied and pasted here.
shp_name = 'local_dir/Bow.shp'
shp_1 = gpd.read_file(shp_name)
shp_1.crs = {'init': 'epsg:4326'} # setting the cordinate system in case it doesn't exists

box_values = box (shp_name)

##
shp_2 = NetCDF_SHP_lat_lon('//datastore/GLOBALWATER/giws_research_water_share/ClimateForcing_Data/ClimateForcing_WFDEI/WFDEI_05d24hr/Rainf_daily_WFDEI_CRU/Rainf_daily_WFDEI_CRU_201612.nc',box_values,'lat','lon',False)
shp_2.crs = {'init': 'epsg:4326'} # setting the cordinate system in case it doesnt exists

shp_int = intersection_shp (shp_1, shp_2)

lat = np.array(shp_int.S_2_lat)
lon = np.array(shp_int.S_2_lon)
W   = np.array(shp_int.AP1N)

case = 1
data = area_ave (case, lat, lon, W,\
                 nc_name,'Tair','tstep','lat','lon','lat','lon')
