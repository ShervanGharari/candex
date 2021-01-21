[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2628350.svg)](https://doi.org/10.5281/zenodo.2628350)

# CANDEX: CAtchment NetcDf EXtractor
This package allows you to extract and aggregate the relevant values from a
cfconventions compliant netcdf files given shapefiles.

CANDEX is a collection of functions that allows extraction of the data from a NetCDF file for a given shapefile such as a basin, catchment. It can map gridded data or model output to any given shapefile and provide area average for a target variable. CANDEX can:

1. [Read a regular lat/lon gridded data or model to any shapefile](./candex/example_regular_lat_lon.ipynb)
2. [Read a rotate lat/lon gridded data or model to any shapefile](./candex/example_rotated_lat_lon.ipynb)
3. Map a non-regular shapefile data, such as Thiessen polygon for example, to any shapefile such as sub-basin.

## The code can be used for the following purposes:

1. Calculation of forcing variables, such as precipitation or temperature and other variables for the effortless model set up. This transfer can be from Thiessen polygon or gridded data, for example.
2. Remapping the output of a hydrological or land surface model to force another model, such as providing the gridded model output in sub-basin for routing.

## Example

The pictures show mapping of the gridded averaged temperature for period for daily temprature of 2010 from the WFDEI dataset. The temperature field is mapped to the sub-basins of the Bow and Oldman River Basin (using Matlab topotoolbox and hydrologically conditioned Hydroshed DEM).

Temperature field in degree celsius for the grids of WFDEI data set

<img src="./fig/Gird.png" width="500">

Temperature field in degree celsius remapped for the sub-basin of Bow and Oldman River Basin

<img src="./fig/Remapped.png.png" width="500">
