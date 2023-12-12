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
Download the package and install it in development mode. From the root directory,
do:

    $ pip install -e .



## Command line interface

``setup.py`` defines an "entry point" for a script to be used as a
command line program. Currently, the only command installed is ``wrfvis_gridcell``.

After installation, just type

    $ wrfvis_gridcell --help

to see what the tool can do.

## Testing

To develop!

I recommend to use [pytest](https://docs.pytest.org) for testing. To test
the package, run

    $ pytest .

in the package root directory.


## License

With the exception of the ``setup.py`` file, which was adapted from the
[sampleproject](https://github.com/pypa/sampleproject) package, all the
code in this repository is dedicated to the public domain.
