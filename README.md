[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2628350.svg)](https://doi.org/10.5281/zenodo.2628350)

# CANDEX: CAtchment NetcDf EXtractor
This package allows you to extract and aggregate the relevant values from a
cfconventions compliant netcdf files given shapefiles.

CANDEX is a collection of functions that allows extraction of the data from a NetCDF file for a given shapefile such as a basin, catchment. It can map gridded data or model output to any given shapefile and provide areas average for a target variable. CANDEX can:

1. [Read a regular lat/lon gridded data or model to any shapefile](https://github.com/ShervanGharari/candex/wiki/Example-1,-Case-1)
2. [Read a rotate lat/lon gridded data or model to any shapefile](https://github.com/ShervanGharari/candex/wiki/Example-2,-Case-2)
3. Map a non-regular shapefile data, such as Thiessen polygon for example, to any shapefile such as sub-basin.
4. [Looping over many shapefile and NetCDF files](https://github.com/ShervanGharari/candex/wiki/Example-4,-Case-4)

## The code can be used for the following purposes:

1. Calculation of forcing variables, such as precipitation or temperature and other variables for the effortless model set up. This transfer can be from Thiessen polygon or gridded data, for example.
2. Mapping the output of a model to force another model, such as providing the gridded model output in sub-basin for routing.

## Features of the functions:

1. No manual ID is needed. The code needs an ID for any target value but it works with lat/lon instead of ID. This reduces the chance of confusion in mapping the source to sink polygons.
2. The code can handle varying resolution and jump in the NetCDF file. in any itteration, the user can ask the grided shapefile to be generated and intersected with the target polygon.
3. It can redo the calculation for a shapefile many times over many files if the files are separated for yearly, monthly or daily time period.

## Example

The pictures show mapping of the gridded averaged temparture for period for 2010 daily from the WFDEI dataset. The temrapture field is mapped to the subbasins of the Bow and Oldman river basin (derived by Matlab topotoolbox and hydrologically conditioned Hydroshed DEM).

![image](https://github.com/ShervanGharari/candex/blob/master/figures/general/temprature_grid.jpg =24x48)

<img src="https://github.com/ShervanGharari/candex/blob/master/figures/general/temprature_grid.jpg" width="48">

![image](https://github.com/ShervanGharari/candex/blob/master/figures/general/temprature_subbasin.jpg =24x48)
