# %%
## INDEX BAND_1 VALUES FROM GEOTIFF ARRAY FOR SINGLE PIXEL
### so far has been done for sensor locations
import matplotlib.pyplot as plt
import numpy as np
import geopandas as gpd
from shapely.geometry import Point
from pyproj import Proj, transform
import utm
from PIL import Image
import panel as pn

# Reading in tif file and convert to np array 
file_name = "Q:\\My Drive\\Code Repositories\\R\\CCRG\\SeaDAS files\\CyanoIndices\\L2022157.L3m_DAY_CYAN_CI_cyano_CYAN_CONUS_300m_8_3.tif"
tiff_data = Image.open(file_name)  # read in TIFF data
cyan_index = np.array(tiff_data)  # convert to numpy array
tiff_data.close()

# Sample data for demonstration (replace this with your actual data)
band1 = cyan_index  # Example array, replace with your actual band1 data

# Lat and Lon bounds of grid (NOT MAPPED MUST USE UTM -- LEVEL-2 DATA)
# metadata from R raster package: (lat and lons in UTM) 1634403, 2234403, 1514826, 2114826  (xmin, xmax, ymin, ymax)
x_min = 1634403
x_max = 2234403
y_min = 1514826
y_max = 2114826

# Number of pixels and pixel size
num_pixels_x, num_pixels_y = 2000,2000 # could be changed back to this: num_pixels_x, num_pixels_y = 2000, 2000
pixel_size = 300  # in meters

dlon = abs((x_max-x_min)/(num_pixels_x-1))
dlat = abs((y_max-y_min)/(num_pixels_y-1))
print('change in lat, lon per pixel', dlat, dlon)

# testing to see if correct band_1 values pull up from indexing
# y_test_utm = 1566576 #latitude TEST ONE PASSED
y_test_utm = 1573476 #latitude TEST TWO PASSED
# x_test_utm = 1766252 # longitude TEST ONE PASSED
x_test_utm = 1762952 # TEST TWO PASSED

lon = x_test_utm
lat = y_test_utm

# Finding the indices for each lat lon pairing -- test with known values (i.e. specific pixel grids)
icol = int(round((((lon-x_max))/dlon)) + (num_pixels_x-1))  # where lon is the longitude of interest 
irow = int(round((y_max-lat)/dlat)) # where lat is the latitude of interest

# testing to ensure this works for corner boxes
print('[irow, icol]',  icol, irow)

band1_value = band1[irow,icol]
print(band1_value)

# next assign coordinate array indices for each sensor location, based on UTM inputs
# sensor.1318 = 1747047.4536078246, 1653374.2879349752
# sensor.1334 = 1706948.9064928016, 1638159.456854635
# sensor.1344 = 1747251.7496186004, 1662336.6205101442 
# sensor.1348 = 1718998.4433066146, 1621145.2433478357
# sensor.1358 =	1691219.0881850973, 1607922.9297252244	
# sensor.1362 = 1745034.263892109, 1642821.4103902446
# sensor.1378 = 1733938.620915762, 1634980.4588427036
# sensor.1562 = 1668123.2874667693, 1654566.0260411282
# sensor.1680 = 1709476.2430031477, 1618371.542780486
# sensor.1806 = 1794125.8191582672, 1644563.7608671915
# sensor.5822 = 1706936.5884716138, 1638430.688035448
# sensor.5838 = 1719341.67516028, 1621446.0612671818
# sensor.9875 = 1747663.1887545092, 1654645.5407499915

# Define the sensor coordinates (UTM)
sensors = {
    "sensor_1318": (1747047.4536078246, 1653374.2879349752),
    "sensor_1334": (1706948.9064928016, 1638159.456854635),
    "sensor_1344": (1747251.7496186004, 1662336.6205101442),
    "sensor_1348": (1718998.4433066146, 1621145.2433478357),
    "sensor_1358": (1691219.0881850973, 1607922.9297252244),
    "sensor_1362": (1745034.263892109, 1642821.4103902446),
    "sensor_1378": (1733938.620915762, 1634980.4588427036),
    "sensor_1562": (1668123.2874667693, 1654566.0260411282),
    "sensor_1680": (1709476.2430031477, 1618371.542780486),
    "sensor_1806": (1794125.8191582672, 1644563.7608671915),
    "sensor_5822": (1706936.5884716138, 1638430.688035448),
    "sensor_5838": (1719341.67516028, 1621446.0612671818),
    "sensor_9875": (1747663.1887545092, 1654645.5407499915),
}

# Function to find the indices for each sensor
def get_indices_and_value(lon, lat):
    icol = int(round((((lon-x_max))/dlon)) + (num_pixels_x-1))  # where lon is the longitude of interest 
    irow = int(round((y_max-lat)/dlat))
    return irow, icol, band1[irow, icol]

# Loop through the sensors and extract the index and band_1 value
# Array to collect results
results = []

# Loop through the sensors and extract the index and band_1 value
for sensor_name, (lon, lat) in sensors.items():
    irow, icol, band1_value = get_indices_and_value(lon, lat)
    results.append([sensor_name, irow, icol, band1_value])

# Convert results to a NumPy array
results_array = np.array(results, dtype=object)

# Specify the output folder and file name
output_folder = "Q:\My Drive\Code Repositories\R\CCRG\SeaDAS files\Band_1_py"
output_file = os.path.join(output_folder, "sensor_data.txt")

# Save the array to a tab-delimited text file
np.savetxt(output_file, results_array, delimiter='\t', fmt='%s', header="Sensor\tRow_Index\tColumn_Index\tBand_1_Value", comments='')

print(f"Results saved to {output_file}")

# note the four outputs are Sensor_ID, Row_Index, Column_Index, and Band_1. I can sync Row_Index and Column_Index to their actual lat and lons eventually

# %%
## INDEX BAND_1 VALUES FROM GEOTIFF ARRAY FOR SINGLE PIXEL
### so far has been done for sensor locations
import matplotlib.pyplot as plt
import numpy as np
import geopandas as gpd
from shapely.geometry import Point
from pyproj import Proj, transform
import utm
from PIL import Image
import panel as pn

# Reading in tif file and convert to np array 
file_name = "Q:\\My Drive\\Code Repositories\\R\\CCRG\\SeaDAS files\\CyanoIndices\\L2022157.L3m_DAY_CYAN_CI_cyano_CYAN_CONUS_300m_8_3.tif"
tiff_data = Image.open(file_name)  # read in TIFF data
cyan_index = np.array(tiff_data)  # convert to numpy array
tiff_data.close()

# Sample data for demonstration (replace this with your actual data)
band1 = cyan_index  # Example array, replace with your actual band1 data

# Extract sensor name, year, and Julian date from the file name
base_name = os.path.basename(file_name)
sensor_name = base_name.split('.')[0].split('_')[0][:1]
year = base_name.split('.')[0].split('_')[0][1:5]
julian_date = base_name.split('.')[0].split('_')[0][5:]

# Lat and Lon bounds of grid (NOT MAPPED MUST USE UTM -- LEVEL-2 DATA)
# metadata from R raster package: (lat and lons in UTM) 1634403, 2234403, 1514826, 2114826  (xmin, xmax, ymin, ymax)
x_min = 1634403
x_max = 2234403
y_min = 1514826
y_max = 2114826

# Number of pixels and pixel size
num_pixels_x, num_pixels_y = 2000,2000 # could be changed back to this: num_pixels_x, num_pixels_y = 2000, 2000
pixel_size = 300  # in meters

dlon = abs((x_max-x_min)/(num_pixels_x-1))
dlat = abs((y_max-y_min)/(num_pixels_y-1))
print('change in lat, lon per pixel', dlat, dlon)

# testing to see if correct band_1 values pull up from indexing
# y_test_utm = 1566576 #latitude TEST ONE PASSED
y_test_utm = 1573476 #latitude TEST TWO PASSED
# x_test_utm = 1766252 # longitude TEST ONE PASSED
x_test_utm = 1762952 # TEST TWO PASSED

lon = x_test_utm
lat = y_test_utm

# Finding the indices for each lat lon pairing -- test with known values (i.e. specific pixel grids)
icol = abs(int(round((((lon-x_max))/dlon)))) #+ (num_pixels_x-1))  # where lon is the longitude of interest 
irow = int(round((y_max-lat)/dlat)) # where lat is the latitude of interest

# testing to ensure this works for corner boxes
print('[irow, icol]',  icol, irow)

band1_value = band1[irow,icol]
print(band1_value)

# Define the sensor coordinates (UTM)
sensors = {
    "sensor_1318": (1747047.4536078246, 1653374.2879349752),
    "sensor_1334": (1706948.9064928016, 1638159.456854635),
    "sensor_1344": (1747251.7496186004, 1662336.6205101442),
    "sensor_1348": (1718998.4433066146, 1621145.2433478357),
    "sensor_1358": (1691219.0881850973, 1607922.9297252244),
    "sensor_1362": (1745034.263892109, 1642821.4103902446),
    "sensor_1378": (1733938.620915762, 1634980.4588427036),
    "sensor_1562": (1668123.2874667693, 1654566.0260411282),
    "sensor_1680": (1709476.2430031477, 1618371.542780486),
    "sensor_1806": (1794125.8191582672, 1644563.7608671915),
    "sensor_5822": (1706936.5884716138, 1638430.688035448),
    "sensor_5838": (1719341.67516028, 1621446.0612671818),
    "sensor_9875": (1747663.1887545092, 1654645.5407499915),
}

# Function to find the indices for each sensor
def get_indices_and_value(lon, lat):
    icol = int(round((((lon-x_max))/dlon)) + (num_pixels_x-1))  # where lon is the longitude of interest 
    irow = int(round((y_max-lat)/dlat))
    return irow, icol, band1[irow, icol]

# Loop through the sensors and extract the index and band_1 value
# Array to collect results
results_S = []

# specify boundary of pixels to scrape in m 
#offsets = [-1350, -1050, -750, -450, -150, 150, 450, 750, 1050, 1350]
offsets_S = np.arange(-1350, 1351, 300)

# Loop through the sensors and extract the index and band_1 value
for sensor_name_key, (sensor_lon, sensor_lat) in sensors.items():
    for dx in offsets_S:
        for dy in offsets_S:
            lon = sensor_lon + dx
            lat = sensor_lat + dy
            irow, icol, band1_value = get_indices_and_value(lon, lat)
            if band1_value is not None:
                results_S.append([sensor_name_key, irow, icol, lon, lat, band1_value])

# Convert results to a NumPy array
results_array_S = np.array(results_S, dtype=object)

# Specify the output folder and file name
output_folder = "Q:\My Drive\Code Repositories\R\CCRG\SeaDAS files\Band_1_py"
# Create the output file name
output_file_name_S = f"{sensor_name}{year}{julian_date}_S.txt"
output_file_S = os.path.join(output_folder, output_file_name_S)

# Save the array to a tab-delimited text file
np.savetxt(output_file_S, results_array_S, delimiter='\t', fmt='%s', header="Sensor\tRow_Index\tColumn_Index\tLon\tLat\tBand_1_Value", comments='')

# ------------------------------------------------------- M 
# Loop through the sensors and extract the index and band_1 value
# Array to collect results
results_M = []

# specify boundary of pixels to scrape in m 
offsets_M = np.arange(-4650, 4651, 300)

# Loop through the sensors and extract the index and band_1 value
for sensor_name_key, (sensor_lon, sensor_lat) in sensors.items():
    for dx in offsets_M:
        for dy in offsets_M:
            lon = sensor_lon + dx
            lat = sensor_lat + dy
            irow, icol, band1_value = get_indices_and_value(lon, lat)
            if band1_value is not None:
                results_M.append([sensor_name_key, irow, icol, lon, lat, band1_value])

# Convert results to a NumPy array
results_array_M = np.array(results_M, dtype=object)

# Specify the output folder and file name
output_folder = "Q:\My Drive\Code Repositories\R\CCRG\SeaDAS files\Band_1_py"
# Create the output file name
output_file_name_M = f"{sensor_name}{year}{julian_date}_M.txt"
output_file_M = os.path.join(output_folder, output_file_name_M)

# Save the array to a tab-delimited text file
np.savetxt(output_file_M, results_array_M, delimiter='\t', fmt='%s', header="Sensor\tRow_Index\tColumn_Index\tLon\tLat\tBand_1_Value", comments='')

# ------------------------------------------------------- L 
# Loop through the sensors and extract the index and band_1 value
# Array to collect results
results_L = []

# specify boundary of pixels to scrape in m 
offsets_L = np.arange(-13650, 13651, 300)

# Loop through the sensors and extract the index and band_1 value
for sensor_name_key, (sensor_lon, sensor_lat) in sensors.items():
    for dx in offsets_L:
        for dy in offsets_L:
            lon = sensor_lon + dx
            lat = sensor_lat + dy
            irow, icol, band1_value = get_indices_and_value(lon, lat)
            if band1_value is not None:
                results_L.append([sensor_name_key, irow, icol, lon, lat, band1_value])

# Convert results to a NumPy array
results_array_L = np.array(results_L, dtype=object)

# Specify the output folder and file name
output_folder = "Q:\My Drive\Code Repositories\R\CCRG\SeaDAS files\Band_1_py"
# Create the output file name
output_file_name_L = f"{sensor_name}{year}{julian_date}_L.txt"
output_file_L = os.path.join(output_folder, output_file_name_L)

# Save the array to a tab-delimited text file
np.savetxt(output_file_L, results_array_L, delimiter='\t', fmt='%s', header="Sensor\tRow_Index\tColumn_Index\tLon\tLat\tBand_1_Value", comments='')

# %% 
## TURNING THIS INTO A FOR LOOP SO I CAN READ MULTIPLE FILES AT ONCE 
### ONCE I GET MY DATA STORAGE SET UP I CAN RUN ALL OF THE FILES
import numpy as np
import os
from PIL import Image
from shapely.geometry import Point
import glob

# Define the sensor coordinates (UTM)
sensors = {
    "sensor_1318": (1747047.4536078246, 1653374.2879349752),
    "sensor_1334": (1706948.9064928016, 1638159.456854635),
    "sensor_1344": (1747251.7496186004, 1662336.6205101442),
    "sensor_1348": (1718998.4433066146, 1621145.2433478357),
    "sensor_1358": (1691219.0881850973, 1607922.9297252244),
    "sensor_1362": (1745034.263892109, 1642821.4103902446),
    "sensor_1378": (1733938.620915762, 1634980.4588427036),
    "sensor_1562": (1668123.2874667693, 1654566.0260411282),
    "sensor_1680": (1709476.2430031477, 1618371.542780486),
    "sensor_1806": (1794125.8191582672, 1644563.7608671915),
    "sensor_5822": (1706936.5884716138, 1638430.688035448),
    "sensor_5838": (1719341.67516028, 1621446.0612671818),
    "sensor_9875": (1747663.1887545092, 1654645.5407499915)
}

# Define the offsets for each range
# updating these offsets as of 08-27-2024 per Zorbas et al 2023 and information on smallest accurate size from Coffer et al 2021
#offsets = {
   # 'S': np.arange(-1350, 1351, 300),
    #'M': np.arange(-4650, 4651, 300),
   # 'L': np.arange(-13650, 13651, 300)
#}
#offsets = [-750, -450, -150, 150, 450, 750]
offsets = {
    'S': np.arange(-750, 751, 300),
    'M': np.arange(-1500, 1501, 300),
    'L': np.arange(-3000, 3001, 300)
}

# Function to find the indices for each sensor
def get_indices_and_value(lon, lat, x_max, y_max, dlon, dlat, band1):
    icol = abs(int(round((((lon-x_max))/dlon))) + (num_pixels_x-1))  # where lon is the longitude of interest 
    irow = int(round((y_max-lat)/dlat))
    return irow, icol, band1[irow, icol]

# Loop to process multiple files
input_folder = "C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Code Repositories\\R\\CCRG\\SeaDAS files\\CyanoIndices"
output_folder = "C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Code Repositories\\R\\CCRG\\SeaDAS files\\Band_1_py_2"

for file_name in glob.glob(os.path.join(input_folder, "*.tif")):
    # Reading in tif file and convert to np array
    tiff_data = Image.open(file_name)
    cyan_index = np.array(tiff_data)
    tiff_data.close()

    band1 = cyan_index  # Example array, replace with your actual band1 data

    # Extract sensor name, year, and Julian date from the file name
    base_name = os.path.basename(file_name)
    sensor_name = base_name.split('.')[0].split('_')[0][:1]
    year = base_name.split('.')[0].split('_')[0][1:5]
    julian_date = base_name.split('.')[0].split('_')[0][5:]

    # Lat and Lon bounds of grid (NOT MAPPED MUST USE UTM -- LEVEL-2 DATA)
    # metadata from R raster package: (lat and lons in UTM) 1634403, 2234403, 1514826, 2114826  (xmin, xmax, ymin, ymax)
    x_min = 1634403
    x_max = 2234403
    y_min = 1514826
    y_max = 2114826

    # Number of pixels and pixel size
    num_pixels_x, num_pixels_y = 2000, 2000
    pixel_size = 300  # in meters

    dlon = abs((x_max-x_min)/(num_pixels_x-1))
    dlat = abs((y_max-y_min)/(num_pixels_y-1))

    # Loop through each range and extract data
    for range_key, offset_values in offsets.items():
        results = []
        for sensor_name_key, (sensor_lon, sensor_lat) in sensors.items():
            for dx in offset_values:
                for dy in offset_values:
                    lon = sensor_lon + dx
                    lat = sensor_lat + dy
                    irow, icol, band1_value = get_indices_and_value(lon, lat, x_max, y_max, dlon, dlat, band1)
                    if band1_value is not None:
                        results.append([sensor_name_key, irow, icol, lon, lat, band1_value])

        # Convert results to a NumPy array
        results_array = np.array(results, dtype=object)

        # Create the output file name
        output_file_name = f"{sensor_name}{year}{julian_date}_{range_key}.txt"
        output_file = os.path.join(output_folder, output_file_name)

        # Save the array to a tab-delimited text file
        np.savetxt(output_file, results_array, delimiter='\t', fmt='%s', header="Name\tPixelY\tPixelX\tUTM_lon\tUTM_lat\tband_1", comments='')

# PixelX and PixelY are flipped in the outputs... I'm just flipping the labels and then going to do some QA/QCing 

print ('done')

# %%
