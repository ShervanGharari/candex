"""
@ author:                  Shervan Gharari
@ Github:                  https://github.com/ShervanGharari/candex
@ author's email id:       sh.gharari@gmail.com
@ license:                 Apache2
"""

# reading the nc name or names, this example gets all the yearly or monthly values in 1990's
nc_name = 'local_dir/*199*.nc' # the local directory should be copied and pasted here.
shp_name = 'local_dir/Bow.shp' # the shapefile name for the Bow river basin in Saskatchewan River Basin in Canada
shp_1 = gpd.read_file(shp_name) # reading the shapefile
shp_1.crs = {'init': 'epsg:4326'} # setting the cordinate system in case it doesn't exists

box_values = box (shp_name,1) # finding the box around the Bow River Basin with buffer of 1 egree

## reading on of the nc file and creating the shapefile from that the corresponde to the NetCDF lat/lon
shp_2 = NetCDF_SHP_lat_lon('//datastore/GLOBALWATER/giws_research_water_share/ClimateForcing_Data/ClimateForcing_WFDEI/WFDEI_05d24hr/Rainf_daily_WFDEI_CRU/Rainf_daily_WFDEI_CRU_201612.nc',box_values,'lat','lon',False)
shp_2.crs = {'init': 'epsg:4326'} # setting the cordinate system

# interseting shp_1 the basin and nf file mesh
shp_int = intersection_shp (shp_1, shp_2)

# changing the data frame into the lat, lon and weight numpy arrays
lat = np.array(shp_int.S_2_lat)
lon = np.array(shp_int.S_2_lon)
W   = np.array(shp_int.AP1N)

# reading the data from NetCDF file
case = 1 # case 1, regualr lat/lon
# getting back the data
data = area_ave (case, lat, lon, W,\
                 nc_name,'Tair','tstep','lat','lon','lat','lon')

# visualziaing the data
plt.plot(data-273)
