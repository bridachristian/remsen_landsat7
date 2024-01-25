# -*- coding: utf-8 -*-
"""
Created on Sun Nov 26 17:41:16 2023

@author: Christian
"""
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path
import os
from osgeo import gdal
import ephem
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pyproj
from pathlib import Path
from scipy.signal import convolve2d


os.getcwd()

global wavelenght
wavelenght = pd.read_csv(
    'wavelenght.csv', sep=';', index_col=0)

global spectral_irradiance
spectral_irradiance = pd.read_csv(
    'spectral_irradiance.csv', sep=';', index_col=0)

global solar_exoatm
solar_exoatm = pd.read_csv(
    'solar_exoatm.csv', sep=';', index_col=0)


def import_file(filepath):
    '''
    Import the data using gdal package
    '''

    dataset = gdal.Open(file_path)
    return dataset


def convert_double(dataset, band_index=1):
    '''
    1. Estract from gdal dataset the data array
    2- Convert the integer values of pixel in float for postprocessing operations
    '''

    # Read data corresponding to the index
    band = dataset.GetRasterBand(band_index)
    data = band.ReadAsArray()

    # Convert data to fload precision
    data_double = data.astype(float)
    return data_double


def get_projection(dataset):
    '''
    Estract projection, coordinates of upper-left corner and pixel size
    '''

    projection = dataset.GetProjection()
    geotransform = dataset.GetGeoTransform()

    # x-coordinate of the upper-left corner of the upper-left pixel.
    x_ul = geotransform[0]
    # w-e pixel resolution / pixel width.
    x_size = geotransform[1]
    # y-coordinate of the upper-left corner of the upper-left pixel.
    y_ul = geotransform[3]
    # n-s pixel resolution / pixel height (negative value for a north-up image).
    y_size = geotransform[5]

    return projection, x_ul, y_ul, x_size, y_size


def plot_band(data_double, band_index=1, value_range=(1, 255)):
    '''
    Plot pixel value raster for a specific band
    The band can be selected by the user. Default is band 1.
    '''

    projection, x_ul, y_ul, x_size, y_size = get_projection(dataset)

    x = np.linspace(x_ul, x_ul+x_size*dataset.RasterXSize, dataset.RasterXSize)
    y = np.linspace(y_ul, y_ul+y_size*dataset.RasterYSize, dataset.RasterYSize)

    xx, yy = np.meshgrid(x, y)
    plt.figure(figsize=(10, 6))  # Set the figure size

    plt.imshow(data_double, cmap='rainbow', extent=(np.min(x), np.max(x), np.min(
        y), np.max(y)), vmin=value_range[0], vmax=value_range[1])
    plt.colorbar(label='Pixel value')
    plt.title(f'Band {band_index}', fontsize=16)
    plt.xlabel('Easting [m]', fontsize=12)
    plt.ylabel('Northing [m]', fontsize=12)
    plt.show()


def get_dimensions(dataset):
    '''
    Get number of pixels in x and in y directions
    '''

    # Get image dimensions (width x height) using the raster's size
    image_width = dataset.RasterXSize
    image_height = dataset.RasterYSize

    return image_width, image_height


def get_filesize(dataset, image_width, image_height):
    '''
    Get the file size in bytes based on data type of a pixel.
    To do that we calculate multiply the data type size
    for the number of pixels and for the number of bands
    '''

    # Get the number of bytes required to store the image
    band1 = dataset.GetRasterBand(1)  # Assuming it's a single band image

    num_bands = dataset.RasterCount

    data_type = gdal.GetDataTypeSize(band1.DataType)
    file_size_bytes = image_width * image_height * \
        (data_type // 8) * num_bands

    return file_size_bytes


def calculate_solar_param(date, time, lat, lon):
    '''
    Calculate the solar parameter based on the location and the date and time of acquisition.
    We use the functions of package 'ephem'.
    '''
    # Observer (representing the Earth-based location)
    observer = ephem.Observer()
    observer.date = date + ' ' + time
    observer.lat = lat
    observer.lon = lon

    # Sun object
    sun = ephem.Sun(observer)

    # Calculate solar zenith angle and Earth-Sun distance
    solar_zenith_angle = sun.alt  # Calculate the angle in radianss
    earth_sun_distance = sun.earth_distance  # Distance in astronomical units

    # Convert radians to degrees for solar zenith angle
    solar_zenith_angle_deg = solar_zenith_angle * (180.0 / ephem.pi)

    return solar_zenith_angle_deg, earth_sun_distance


def calc_radiance(D, band_index=1):
    '''
    Calculate the spectral radiance.
    We use the Radiometric correction proposed in Rees book (p.355)
    Linear interpolation of pixel using the spectral radiance table provided.
    '''

    numer = D-1
    denom = 255 - 1

    sp_rad_min = spectral_irradiance.iloc[band_index-1, 0]
    sp_rad_max = spectral_irradiance.iloc[band_index-1, 1]

    d_new = (numer/denom) * (sp_rad_max - sp_rad_min) + sp_rad_min
    return d_new


def calc_planetary_reflec(R, band_index=1):
    '''
    Calculate the planetary reflactance.
    We use the Planetary reflectance proposed in Rees book (p.356)
    Correction by solar parameter of measured radiance by the satellite
    '''
    exoatm_param = solar_exoatm.iloc[band_index-1, 0]
    numer = (R * np.pi * (earth_sun_distance)**2)
    denom = (exoatm_param * np.cos(np.radians(solar_zenith_angle_deg)))

    rho = numer/denom
    return rho


def calculate_size(dataset):
    '''
    Extract projection, pixel size and map extend from dataset
    '''
    projection, x_ul, y_ul, x_size, y_size = get_projection(dataset)

    pixel_size_x_meters = abs(x_size)
    pixel_size_y_meters = abs(y_size)

    width_km = dataset.RasterXSize * pixel_size_x_meters / 1000
    height_km = dataset.RasterYSize * pixel_size_y_meters / 1000

    return projection, pixel_size_x_meters, pixel_size_y_meters, width_km, height_km


def plot_geographic(raster_data, point_lat, point_lon, band_index=1):
    '''
    Tools used to plot in geographical coordinates instead of UTM
    We also add the point of the airport, find on Google maps.
    Airport location is confirmed.
    '''

    # raster_data = data_spectral[0]

    projection = dataset.GetProjection()
    geotransform = dataset.GetGeoTransform()

    # Get the UTM zone information from the projection
    # Extract UTM zone from projection string
    utm_zone = int(projection.split('UTM zone ')[-1][:2])

    # Define the UTM and lat-lon coordinate systems
    utm_crs = pyproj.Proj(proj='utm', zone=utm_zone, ellps='WGS84')
    latlon_crs = pyproj.Proj(proj='latlong', datum='WGS84')

    # Get the coordinates of the corners of the raster in UTM
    # Upper-left corner coordinates
    x_ul, y_ul = geotransform[0], geotransform[3]
    # Lower-right corner x-coordinate
    x_lr = x_ul + geotransform[1] * dataset.RasterXSize
    # Lower-right corner y-coordinate
    y_lr = y_ul + geotransform[5] * dataset.RasterYSize

    # Transform UTM coordinates to lat-lon
    x_ul_lon, y_ul_lat = pyproj.transform(utm_crs, latlon_crs, x_ul, y_ul)
    x_lr_lon, y_lr_lat = pyproj.transform(utm_crs, latlon_crs, x_lr, y_lr)

    # Create a meshgrid of lat-lon coordinates for the raster
    x_lon = np.linspace(x_ul_lon, x_lr_lon, dataset.RasterXSize)
    y_lat = np.linspace(y_lr_lat, y_ul_lat, dataset.RasterYSize)
    xx, yy = np.meshgrid(x_lon, y_lat)

    # Plot raster data in lat-lon coordinates
    plt.figure(figsize=(10, 6))
    # Extent in lat-lon
    plt.imshow(raster_data, extent=(x_ul_lon, x_lr_lon, y_lr_lat, y_ul_lat),
               cmap='rainbow')
    plt.colorbar(
        label='Spectral radiance [$W m^{-2} sr^{-1} \mu m^{-1}$]')
    plt.scatter(point_lon, point_lat, marker='+', s=500, color='red')
    plt.title(f'Band {band_index} - Spectral Radiance', fontsize=16)
    plt.xlabel('Longitude [deg]', fontsize=12)
    plt.ylabel('Latitude [deg]', fontsize=12)
    plt.annotate('+ Airport', (29.5, 60.5), color='red', ha='left',
                 fontsize=12, bbox=dict(facecolor='white', alpha=0.5))
    plt.show()


def calc_brightness_temp(R, band_index=6):
    '''
    Calculate the brightness temperature based on Landat7 documentation.
    Formuala and K1,K2 param can be obtaindef from:
         Landsat 7 (L7) Data Users Handbook
    '''

    K1 = 666.09  # Watts*m-2 * sr-1 * µm-1
    K2 = 1282.71  # K
    numer = K2
    denom = np.log((K1/R) + 1)
    T = (numer/denom)
    return T


def plot_brightess(raster_data, value_range=(250, 360), save_fig=False):
    '''
    Tool used to plot the brightess temperature for the band 6 - Thermal
    '''

    # raster_data = Tb
    projection, x_ul, y_ul, x_size, y_size = get_projection(dataset)

    x = np.linspace(x_ul, x_ul+x_size*dataset.RasterXSize, dataset.RasterXSize)
    y = np.linspace(y_ul, y_ul+y_size*dataset.RasterYSize, dataset.RasterYSize)

    xx, yy = np.meshgrid(x, y)

    # cm = 1/2.54
    plt.figure(figsize=(10, 6), dpi=300)  # Set the figure size
    plt.imshow(raster_data, extent=(np.min(x), np.max(x), np.min(
        y), np.max(y)), cmap='plasma', vmin=value_range[0], vmax=value_range[1])
    plt.colorbar(label='Brightness Temperture [K]')
    plt.title(f'Brightness Temperture', fontsize=12)
    plt.xlabel('Easting [m]', fontsize=12)
    plt.ylabel('Northing [m]', fontsize=12)
    plt.show()


def create_custom_colormap(min_val, max_val, midpoint, color_low, color_high):
    # Define the colors for low, middle, and high values
    colors = [color_low, (1, 1, 1), color_high]

    # Define the positions for each color
    positions = [0, (midpoint - min_val) / (max_val - min_val), 1]

    # Create the colormap
    custom_cmap = LinearSegmentedColormap.from_list(
        "custom_colormap", list(zip(positions, colors)), N=256)

    return custom_cmap


def plot_contrail_brightess(raster_data, value_range, save_fig=False):
    '''
    Tool used to plot the brightess temperature for the band 6 - Thermal
    '''
    projection, x_ul, y_ul, x_size, y_size = get_projection(dataset)

    x = np.linspace(x_ul, x_ul+x_size*dataset.RasterXSize, dataset.RasterXSize)
    y = np.linspace(y_ul, y_ul+y_size*dataset.RasterYSize, dataset.RasterYSize)

    custom_palette = create_custom_colormap(
        value_range[0], value_range[1], midpoint=0,
        color_low='blue', color_high='red')

    plt.figure(figsize=(10, 6), dpi=300)  # Set the figure size
    plt.imshow(raster_data, extent=(np.min(x), np.max(x), np.min(
        y), np.max(y)), cmap=custom_palette, vmin=value_range[0], vmax=value_range[1])

    plt.colorbar(label='Relative Brightness Temperture [K]')
    plt.title(f'Relative Brightness Temperture', fontsize=12)
    plt.xlabel('Easting [m]', fontsize=12)
    plt.ylabel('Northing [m]', fontsize=12)

    fig_name = f'relative_brightness_band6.png'
    if save_fig:
        plt.savefig(outdir / fig_name)
        print(f"Figure saved as '{fig_name}'")
        plt.show()

    else:
        plt.show()


def band_selection(raster, index_r, index_g, index_b):
    red = raster[index_r - 1]  # Example red band data
    green = raster[index_g - 1]  # Example green band data
    blue = raster[index_b - 1]  # Example blue band data

    return red, green, blue


def calc_rgb(red_band, green_band, blue_band):
    ''' Calculate the rgb image combining the 3 bands'''

    rgb_image = np.stack((red_band, green_band, blue_band), axis=-1)
    rgb_image = rgb_image/np.max(rgb_image)

    rgb_pixels = np.stack((red_band, green_band, blue_band), axis=-1)

    return rgb_image, rgb_pixels


def plot_rgb_map(raster_3d, title):
    ''' plot the map in rgb format '''
    projection, x_ul, y_ul, x_size, y_size = get_projection(dataset)

    x = np.linspace(x_ul, x_ul+x_size*dataset.RasterXSize, dataset.RasterXSize)
    y = np.linspace(y_ul, y_ul+y_size*dataset.RasterYSize, dataset.RasterYSize)

    plt.figure(figsize=(10, 6), dpi=300)  # Set the figure size
    plt.imshow(raster_3d, extent=(np.min(x), np.max(x), np.min(
        y), np.max(y)))

    plt.title(f'{title}', fontsize=12)
    plt.xlabel('Easting [m]', fontsize=12)
    plt.ylabel('Northing [m]', fontsize=12)


def plot_rgb_hist(raster_3d, title):
    ''' Plot histogram of pixel values for the different bands in rgb '''

    plt.figure(figsize=(6, 5), dpi=300)  # Set the figure size
    plt.hist(raster_3d[:, :, 0].flatten(), bins=255, color='red', alpha=0.7,)
    plt.hist(raster_3d[:, :, 1].flatten(), bins=255, color='blue', alpha=0.7)
    plt.hist(raster_3d[:, :, 2].flatten(), bins=255, color='green', alpha=0.7)
    plt.xlim(xmin=1, xmax=255)
    plt.title(f'{title}', fontsize=12)
    plt.ylabel('Counts', fontsize=12)
    plt.xlabel('Pixel values', fontsize=12)
    plt.grid(alpha=0.2)
    plt.show()


def contrast_stretch(raster, min_val, max_val):
    ''' contrast manipulation:
        pixel value is stretched from the minimum and maximum value'''
    raster = (raster - min_val) * 255.0 / (max_val - min_val)
    return raster


def spatial_avg(raster, window=5):
    ''' spatial filtering: moving average with a fixed window '''

    # Define a 3x3 moving average kernel
    kernel = np.ones((window, window)) / window**2

    # Perform convolution
    result = convolve2d(raster, kernel, mode='same', boundary='wrap')

    return result


def contrast_manipulation(red_band, green_band, blue_band):
    ''' constrast manipulation applied to 3 bands for RGB calculation'''

    red_star = contrast_stretch(red_band, np.percentile(
        red_band.flatten(), 1), np.percentile(red_band.flatten(), 99))
    green_star = contrast_stretch(green_band, np.percentile(
        green_band.flatten(), 1), np.percentile(green_band.flatten(), 99))
    blue_star = contrast_stretch(blue_band, np.percentile(
        blue_band.flatten(), 1), np.percentile(blue_band.flatten(), 99))

    return red_star, green_star, blue_star


def spatial_filtering(red_band, green_band, blue_band):
    ''' spatial filtering applied to 3 bands for RGB calculation'''

    red_spatial = spatial_avg(red_band)
    green_spatial = spatial_avg(green_band)
    blue_spatial = spatial_avg(blue_band)

    return red_spatial, green_spatial, blue_spatial


def plot_NDVI_map(raster, value_range=(-1, 1)):
    ''' plot the NDVI map  '''
    projection, x_ul, y_ul, x_size, y_size = get_projection(dataset)

    x = np.linspace(x_ul, x_ul+x_size*dataset.RasterXSize, dataset.RasterXSize)
    y = np.linspace(y_ul, y_ul+y_size*dataset.RasterYSize, dataset.RasterYSize)

    custom_palette = create_custom_colormap(
        value_range[0], value_range[1], midpoint=0,
        color_low='brown', color_high='green')

    plt.figure(figsize=(10, 6), dpi=300)  # Set the figure size
    plt.imshow(raster, extent=(np.min(x), np.max(x), np.min(
        y), np.max(y)), cmap=custom_palette, vmin=value_range[0], vmax=value_range[1])

    plt.colorbar(label='NDVI')
    plt.title(f'Normalize Difference Vegetation Index', fontsize=12)
    plt.xlabel('Easting [m]', fontsize=12)
    plt.ylabel('Northing [m]', fontsize=12)
    plt.show()


def plot_NDSI_map(raster, value_range=(-1, 1)):
    ''' plot the NDSI map  '''
    projection, x_ul, y_ul, x_size, y_size = get_projection(dataset)

    x = np.linspace(x_ul, x_ul+x_size*dataset.RasterXSize, dataset.RasterXSize)
    y = np.linspace(y_ul, y_ul+y_size*dataset.RasterYSize, dataset.RasterYSize)

    custom_palette = create_custom_colormap(
        value_range[0], value_range[1], midpoint=0,
        color_low='darkgreen', color_high='deepskyblue')

    plt.figure(figsize=(10, 6), dpi=300)  # Set the figure size
    plt.imshow(raster, extent=(np.min(x), np.max(x), np.min(
        y), np.max(y)), cmap=custom_palette, vmin=value_range[0], vmax=value_range[1])

    plt.colorbar(label='NDSI')
    plt.title(f'Normalize Difference Snow Index', fontsize=12)
    plt.xlabel('Easting [m]', fontsize=12)
    plt.ylabel('Northing [m]', fontsize=12)
    plt.show()


def image_processing(red_band, green_band, blue_band):

    figure, hist = calc_rgb(red_band, green_band, blue_band)

    # Contrast Manipulation
    red_star, green_star, blue_star = contrast_manipulation(
        red_band, green_band, blue_band)
    figure_contrast, hist_contrast = calc_rgb(
        red_star, green_star, blue_star)

    # Spatial filtering
    red_spatial, green_spatial, blue_spatial = spatial_filtering(
        red_star, green_star, blue_star)

    figure_spatial, hist_spatial = calc_rgb(
        red_spatial, green_spatial, blue_spatial)

    return figure, hist, figure_contrast, hist_contrast, figure_spatial, hist_spatial


if __name__ == '__main__':

    ''' Point 1: import the data, convert to double precision, plotting '''

    mydir = Path("C:/Users/Christian/OneDrive/Desktop/Family/Christian/MasterMeteoUnitn/Corsi/3_terzo_semestre/AtmosphericRadiationRemoteSensing/Project")
    file_path = 'SPBL7.tif'

    filedir = mydir / file_path

    global outdir
    outdir = mydir / 'output'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    global dataset
    dataset = import_file(filedir)

    bands = range(1, 8)
    data_double = []

    for band_index in bands:
        data_double.append(convert_double(dataset, band_index))

    for band_index in bands:
        plot_band(data_double[band_index-1],  band_index)

    ''' Point 2: Size in pixels and the number of bytes '''

    image_width, image_height = get_dimensions(dataset)
    print(f"Image dimensions: {image_width} x {image_height} pixels")

    file_size_bytes = get_filesize(dataset, image_width, image_height)
    print(f"File size: {file_size_bytes} bytes")

    ''' Point 3a: Calculate solar zenith angle and the Earth-Sun distance '''

    # Acquisition date and time
    date = '2004/02/12'
    time = '12:40:00'  # in UTC ---> local time 15:40

    # Image location (S.Petersburg --> online)
    lat = '59.80310883711776'  # London latitude
    lon = '30.26872228346058'   # London longitude

    global solar_zenith_angle_deg
    global earth_sun_distance
    solar_zenith_angle_deg, earth_sun_distance = calculate_solar_param(
        date, time, lat, lon)

    # Print the calculated values
    print(f"Solar Zenith Angle: {solar_zenith_angle_deg} degrees")
    print(f"Earth-Sun Distance: {earth_sun_distance} AU")

    ''' Point 3b: Calculate the planetary reflectance for bands 1-5 and 7 '''
    data_spectral = []
    data_plan_reflact = []

    for band_index in bands:
        data_spectral.append(calc_radiance(
            data_double[band_index-1], band_index))
        data_plan_reflact.append(
            calc_planetary_reflec(data_spectral[band_index-1], band_index))

    ''' Point 3c: image projection, pixel size ( in metres) and width and height ( in km) '''
    projection, pixel_size_x_meters, pixel_size_y_meters, width_km, height_km = calculate_size(
        dataset)
    print(f"Projection: {projection}")
    print(
        f"Pixel size in meters (X, Y): {pixel_size_y_meters}, {pixel_size_y_meters}")
    print(f"Image Width: {width_km} km")
    print(f"Image Height: {height_km} km")

    ''' Point 4: Find the geographical coord(lat lon) of Pulkovo airport '''
    # aeroport coords from Google Maps
    point_lat = 59.80310883711776
    point_lon = 30.26872228346058

    for band_index in bands:
        plot_geographic(data_spectral[band_index-1],
                        point_lat, point_lon, band_index)

    ''' Point 5a: Brightness temperature of band 6 '''

    band_index = 6
    band6 = data_spectral[5]
    Tb = calc_brightness_temp(band6, band_index=6)

    plot_brightess(Tb, save_fig=False)

    ''' Point 5b: Condense of air, compare with water and ice temperature'''
    ''' Water - Baltic Sea'''
    water_mean_pixel = 96.842  # std = 1.244 From Zonal Statistics poligon shapefile Qgis
    water_mean_radiance = calc_radiance(water_mean_pixel, band_index=6)
    water_mean_Tb = calc_brightness_temp(
        water_mean_radiance, band_index=6)

    water_shifted_Tb = Tb - water_mean_Tb
    value_range = (-10, 25)
    plot_contrail_brightess(water_shifted_Tb, value_range, save_fig=False)

    ''' Ice - Lake Ladoga'''
    ice_mean_pixel = 87.933  # std = 1.404 From Zonal Statistics poligon shapefile Qgis.
    ice_mean_radiance = calc_radiance(ice_mean_pixel, band_index=6)
    ice_mean_Tb = calc_brightness_temp(
        ice_mean_radiance, band_index=6)

    ice_shifted_Tb = Tb - ice_mean_Tb
    value_range = (-5, 30)
    plot_contrail_brightess(ice_shifted_Tb, value_range, save_fig=False)

    ''' Point 6: contrast manipulation and spatial filtering '''

    '''a) 321 – true-colour composite '''
    red_band, green_band, blue_band = band_selection(data_double, 3, 2, 1)

    figure, hist, figure_contrast, hist_contrast, figure_spatial, hist_spatial = image_processing(
        red_band, green_band, blue_band)

    title = 'True color image: 321'
    plot_rgb_map(figure, title)
    plot_rgb_hist(hist, title)

    title = 'True color image contrast: 321'
    plot_rgb_map(figure_contrast, title)
    plot_rgb_hist(hist_contrast, title)

    title = 'True color image spatial: 321'
    plot_rgb_map(figure_spatial, title)
    plot_rgb_hist(hist_spatial, title)

    '''b) 432 – false colour infrared '''
    red_band, green_band, blue_band = band_selection(data_double, 4, 3, 2)
    figure, hist, figure_contrast, hist_contrast, figure_spatial, hist_spatial = image_processing(
        red_band, green_band, blue_band)

    title = 'False colour infrared image: 432'
    plot_rgb_map(figure, title)
    plot_rgb_hist(hist, title)

    title = 'False colour infrared image contrast: 432'
    plot_rgb_map(figure_contrast, title)
    plot_rgb_hist(hist_contrast, title)

    title = 'False colour infrared image spatial: 432'
    plot_rgb_map(figure_spatial, title)
    plot_rgb_hist(hist_spatial, title)

    '''c) 543 – snowcover '''

    red_band, green_band, blue_band = band_selection(data_double, 5, 4, 3)
    figure, hist, figure_contrast, hist_contrast, figure_spatial, hist_spatial = image_processing(
        red_band, green_band, blue_band)

    title = 'Snow cover image: 543'
    plot_rgb_map(figure, title)
    plot_rgb_hist(hist, title)

    title = 'Snow cover infrared contrast: 543'
    plot_rgb_map(figure_contrast, title)
    plot_rgb_hist(hist_contrast, title)

    title = 'Snow cover infrared spatial: 543'
    plot_rgb_map(figure_spatial, title)
    plot_rgb_hist(hist_spatial, title)

    '''c) 666 – thermal infrared channel in greyscale '''

    red_band, green_band, blue_band = band_selection(data_double, 6, 6, 6)
    figure, hist, figure_contrast, hist_contrast, figure_spatial, hist_spatial = image_processing(
        red_band, green_band, blue_band)

    title = 'Thermal infrared image: 666'
    plot_rgb_map(figure, title)
    plot_rgb_hist(hist, title)

    title = 'Thermal infrared infrared contrast: 666'
    plot_rgb_map(figure_contrast, title)
    plot_rgb_hist(hist_contrast, title)

    title = 'Thermal infrared infrared spatial: 666'
    plot_rgb_map(figure_spatial, title)
    plot_rgb_hist(hist_spatial, title)

    ''' Point 6b: dark-pixel for bands 2,3,4,5'''
    band2 = data_double[1]
    data2_stack = np.stack([band2, band2, band2])
    plot_rgb_hist(np.stack([band2, band2, band2]), 'Band2 Histogram ')
    band2_star = contrast_stretch(
        band2, np.percentile(band2.flatten(), 1), np.percentile(band2.flatten(), 99))
    plot_rgb_hist(
        np.stack([band2_star, band2_star, band2_star]), 'Band2 stretched Histogram ')
    band2_spatial = spatial_avg(band2_star)

    band3 = data_double[2]
    data3_stack = np.stack([band3, band3, band3])
    plot_rgb_hist(np.stack([band3, band3, band3]), 'Band3 Histogram ')
    band3_star = contrast_stretch(
        band3, np.percentile(band3.flatten(), 1), np.percentile(band3.flatten(), 99))
    plot_rgb_hist(
        np.stack([band3_star, band3_star, band3_star]), 'Band3 stretched Histogram ')
    band3_spatial = spatial_avg(band3_star)

    band4 = data_double[3]
    data4_stack = np.stack([band4, band4, band4])
    plot_rgb_hist(np.stack([band4, band4, band4]), 'Band4 Histogram ')
    band4_star = contrast_stretch(
        band4, np.percentile(band4.flatten(), 0.3), np.percentile(band4.flatten(), 99))
    plot_rgb_hist(
        np.stack([band4_star, band4_star, band4_star]), 'Band4 stretched Histogram ')
    band4_spatial = spatial_avg(band4_star)

    band5 = data_double[4]
    data5_stack = np.stack([band5, band5, band5])
    plot_rgb_hist(np.stack([band5, band5, band5]), 'Band5 Histogram ')
    band5_star = contrast_stretch(
        band5, np.percentile(band5.flatten(), 1), np.percentile(band5.flatten(), 99))
    plot_rgb_hist(
        np.stack([band5_star, band5_star, band5_star]), 'Band5 stretched Histogram ')
    band5_spatial = spatial_avg(band5_star)

    ''' Point 6c: NDVI, NDSI'''
    NDVI = (band4 - band3) / (band4 + band3)
    plot_NDVI_map(NDVI)

    NDVI_dark_filt = (band4_spatial - band3_spatial) / \
        (band4_spatial + band3_spatial)
    plot_NDVI_map(NDVI_dark_filt)

    NDSI = (band2 - band5) / (band2 + band5)
    plot_NDSI_map(NDSI)

    NDSI_dark_filt = (band2_spatial - band5_spatial) / \
        (band2_spatial + band5_spatial)
    plot_NDSI_map(NDSI_dark_filt)
