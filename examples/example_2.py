"""
@ author:                  Shervan Gharari
@ Github:                  ./shervangharari/candex
@ author's email id:       sh.gharari@gmail.com
@ license:                 Apache2
"""

nc_name = 'local_dir/TMP.nc' # the local directory should be copied and pasted here.
shp_name = 'local_dir/Bow_360.shp'
# shp_name = 'local_dir/Bow.shp'

shp_1 = gpd.read_file(shp_name)
print(shp_1)
print(shp_1.crs)
shp_1.crs = {'init': 'epsg:4326'} # setting the cordinate

box_values = box (shp_name)

shp_2 = NetCDF_SHP_lat_lon(nc_name,box_values,'latitude','longitude', False)
# shp_2 = NetCDF_SHP_lat_lon(nc_name,box_values,'latitude','longitude', True)
shp_2.crs = {'init': 'epsg:4326'}

shp_int = intersection_shp (shp_1, shp_2)

lat = np.array(shp_int.S_2_lat)
lon = np.array(shp_int.S_2_lon)
W   = np.array(shp_int.AP1N)

case = 2
data = area_ave (case, lat, lon, W,\
                 nc_name,'TMP_40maboveground','time','latitude','longitude','y','x')
