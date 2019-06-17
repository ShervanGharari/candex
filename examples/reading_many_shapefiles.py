# example for many shapefile 4
nc_name = '//datastore/GLOBALWATER/giws_research_water_share/ClimateForcing_Data/\
ClimateForcing_WFDEI/WFDEI_05d24hr/Tair_daily_WFDEI/*2010*.nc' # the local directory should be copied and pasted here.
shp_name = 'F:/Intercomparison/VIC_GRU/Erie_subwatershed/4_LEB_boundary_subwatershed.shp'
shp_1 = gpd.read_file(shp_name)
shp_1.crs = {'init': 'epsg:4326'} # setting the cordinate system in case it doesnt exists
shp_1.plot(facecolor='none', edgecolor='k')
plt.savefig('C:/Users/shg096/Dropbox/candex/example_4/basin.jpg')

box_value = box (shp_1)
shp_2 = NetCDF_SHP_lat_lon('//datastore/GLOBALWATER/giws_research_water_share/\
ClimateForcing_Data/ClimateForcing_WFDEI/WFDEI_05d24hr/Rainf_daily_WFDEI_CRU/Rainf_daily_WFDEI_CRU_201612.nc',\
                           box_value,'lat','lon',False)
shp_2.crs = {'init': 'epsg:4326'} # setting the cordinate system in case it doesnt exists
shp_2.plot(facecolor='none', edgecolor='k')
plt.savefig('C:/Users/shg096/Dropbox/candex/example_4/netcdf_shp.jpg')


shp_int = intersection_shp (shp_1, shp_2)
shp_int.plot(facecolor='none', edgecolor='k')
plt.savefig('C:/Users/shg096/Dropbox/candex/example_4/intersection.jpg')

print(shp_int)
    
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
