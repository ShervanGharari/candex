[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2628350.svg)](https://doi.org/10.5281/zenodo.2628350)

# CANDEX: CAtchment NetcDf EXtractor
This package allows you to extract and aggregate the relevant values from a
cfconventions compliant netcdf files given shapefiles.

CANDEX is a collection of functions that allows extraction of the data from a NetCDF file for a given shapefile such as a basin, catchment. It can map gridded data or model output to any given shapefile and provide area average for a target variable. CANDEX can:

1. [Remap a regular lat/lon gridded data or model to any shapefile](./candex/example_regular_lat_lon.ipynb)
2. [Remap a rotate lat/lon gridded data or model to any shapefile](./candex/example_rotated_lat_lon.ipynb)
3. Remap a non-regular shapefile data, such as Thiessen polygon for example, to any shapefile such as sub-basin.

## The code can be used for the following purposes:

1. Remapping the relevant forcing variables, such as precipitation or temperature and other variables for the effortless model set up. This transfer can be from Thiessen polygon or gridded data, for example, to computational units, hydrological model for example.
2. Remapping the output of a hydrological or land surface model to force another model, such as providing the gridded model output in sub-basin for routing.

## Example

The two figures show remapping of the gridded temperature from ERA5 data set to subbasin of South Saskatchewan River at Medicine Hat.

### Original gridded temperature field:

<img src="./fig/Gird.png" width="2000">

### Remapped temperature field to the subbasins

<img src="./fig/Remapped.png" width="2000"> 

## Publication that have used CANDEX so far:

[Gharari, S., Clark, M. P., Mizukami, N., Knoben, W. J. M., Wong, J. S., and Pietroniro, A.: Flexible vector-based spatial configurations in land models, Hydrol. Earth Syst. Sci., 24, 5953–5971, https://doi.org/10.5194/hess-24-5953-2020, 2020] (https://doi.org/10.5194/hess-24-5953-2020)
[Sheikholeslami, R., Gharari, S. Papalexiou, S. M., Clark, M. P.: VISCOUS: A Variance-Based Sensitivity Analysis Using Copulas for Efficient Identification of Dominant Hydrological Processes, submitted to Water Resources Research, https://doi.org/10.1002/essoar.10505333.1, 2020] (https://doi.org/10.1002/essoar.10505333.1)
