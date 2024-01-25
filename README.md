# A collection of tools to perform Data Analysis on Landsat7 data

**remsen_landsat7** offers command line tools to performs some analysis on Landsat7 data.

It was written for the University of Innsbruck's Atmospheric Radiation and Remote Sensing
course for the semester project

## HowTo

Make sure you have all dependencies installed. These are:
- numpy
- pandas
- matplotlib
- pathlib
- os
- osgeo
- ephem
- pyproj
- pytest

## Usage
Download the package in your folder. You can find the script

    $ landsat.py


That allows to perform the entire analysis. 
Make sure to download also the csv files:
- spectral_irradiance.csv
- solar_exoatm.csv
- wavelenght.csv

These files are use to convert pixel in spectral radiance and to calculate planetary extraatmospheric spectral radiance. 

## License

With the exception of the ``setup.py`` file, which was adapted from the
[sampleproject](https://github.com/pypa/sampleproject) package, all the
code in this repository is dedicated to the public domain.
