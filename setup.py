from setuptools import find_packages, setup

setup(
    name='candex',
    version='0.1.0',
    license='Apache Software License',
    author=('Shervan Gharari', 'Berend Weel'),
    install_requires=[
        'numpy',
        'geopandas',
        'xarray',
        'pandas',
        'shapely',
    ],
    author_email='sh.gharari@gmail.com',
    description=(
        'Extract catchment data from netcdf file based on a catchment shapefile'
    ),
    long_description=(
        'Extract catchment data from netcdf file based on a catchment shapefile'
    ),
    packages=find_packages(),
    include_package_data=True,
    platforms='any',
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Development Status :: 3 - Alpha',
        'Natural Language :: English',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: OS Independent',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
)
