"""
@ author:                  Shervan Gharari
@ Github:                  https://github.com/ShervanGharari/candex
@ author's email id:       sh.gharari@gmail.com
@ license:                 Apache2
"""

# reading the name of the NetCDF file
nc_name = 'local_dir/TMP.nc' # the local directory should be copied and pasted here.
shp_name = 'local_dir/Bow_360.shp' # rading the shapefile which is 0 to 360 
# shp_name = 'local_dir/Bow.shp' # reading the shapefile which is -180 to 180

shp_1 = gpd.read_file(shp_name) # reading the shapefile
shp_1.crs = {'init': 'epsg:4326'} # setting the spatial reference system if not exists

box_values = box (shp_name,1) # box arounf the Bow river basin with 1 degree buffer

shp_2 = NetCDF_SHP_lat_lon(nc_name,box_values,'latitude','longitude', False) # creat the NetCDF mesh shapefile not correct for 0 to 360
# shp_2 = NetCDF_SHP_lat_lon(nc_name,box_values,'latitude','longitude', True) # creat the NetCDF mesh shapefile correct for -180 to 180
shp_2.crs = {'init': 'epsg:4326'} # set the cordinate system

# intersection
shp_int = intersection_shp (shp_1, shp_2)

# change of the data frame into numpy array
lat = np.array(shp_int.S_2_lat)
lon = np.array(shp_int.S_2_lon)
W   = np.array(shp_int.AP1N)


case = 2 # case 2 rotated lat and lon
data = area_ave (case, lat, lon, W,\
                 nc_name,'TMP_40maboveground','time','latitude','longitude','y','x')

# visualization of data
plt.plot(data-273.3)
