# %% 
import matplotlib.pyplot as plt
import geopandas as gpd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd
import numpy as np
from matplotlib.patches import FancyArrowPatch

# the windsor point needs to be updated, I think the Lat and Lon is incorrect for that one

def draw_compass_rose(ax, center, size):
    x, y = center
    directions = ['E', 'N', 'W', 'S']
    for i, direction in enumerate(directions):
        angle = np.pi / 2 * i
        dx = size * np.cos(angle)
        dy = size * np.sin(angle)
        ax.add_patch(FancyArrowPatch((x, y), (x + dx, y + dy), mutation_scale=size*15,
                                     arrowstyle='-|>', color='black', transform=ccrs.PlateCarree()._as_mpl_transform(ax)))
        ax.text(x + 1.5 * dx, y + 1.5 * dy, direction, ha='center', va='center', fontsize=10, 
                transform=ccrs.PlateCarree()._as_mpl_transform(ax))
        
# Function to draw distance marker
def draw_distance_marker(ax, start, distance_km, text, color='black'):
    x_start, y_start = start
    x_end = x_start + (distance_km / 111)  # Approximate conversion from km to degrees
    ax.plot([x_start, x_end], [y_start, y_start], color=color, transform=ccrs.PlateCarree())
    ax.text((x_start + x_end) / 2, y_start, text, ha='center', va='bottom', 
            transform=ccrs.PlateCarree()._as_mpl_transform(ax), fontsize=10, color=color)
    
# Load the GeoJSON file
shapefile_path = "C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Code Repositories\\R\\CCRG\\SeaDAS files\\AS_shapefiles"

# Load the shapefile
albemarle_sound = gpd.read_file(shapefile_path)
# Read the sensor data from the provided text file
sensor_data_path = "C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Code Repositories\\R\\CCRG\\SeaDAS files\\shapefiles\\purpleair_pins_numericlabels.txt"
sensor_data = pd.read_csv(sensor_data_path, sep='\t', comment='#')

# Extract relevant columns (longitude and latitude)
sensor_data = sensor_data[['Lon', 'Lat']]

# Apply jitter to sensor points
jitter_strength = 0.02 # Adjust jitter strength as necessary
sensor_data['Lon'] += np.random.uniform(-jitter_strength, jitter_strength, sensor_data.shape[0])
sensor_data['Lat'] += np.random.uniform(-jitter_strength, jitter_strength, sensor_data.shape[0])

# Create the plot
fig, ax = plt.subplots(1, 1, figsize=(10, 10), subplot_kw={'projection': ccrs.PlateCarree()})

# Add features to the map
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.LAND)

# Customize the OCEAN color
ocean = cfeature.NaturalEarthFeature(
    'physical', 'ocean', '10m',
    edgecolor='face',
    facecolor=cfeature.COLORS['water']
)
ax.add_feature(ocean, facecolor='#6699CC', alpha = 0.8) 

# Plot the Albemarle Sound
albemarle_sound.plot(ax=ax, edgecolor='#6699CC')

# Set the extent to focus on the Albemarle Sound region
ax.set_extent([-77.25, -75.5, 35.7, 36.6], crs=ccrs.PlateCarree())

# Plot the sensor locations
ax.scatter(sensor_data['Lon'], sensor_data['Lat'], edgecolor='black', facecolor='none', s=50, zorder=5, alpha = 0.8, transform=ccrs.PlateCarree(), label='PurpleAir Sensors', linewidths=1.5)  

# Add compass rose
draw_compass_rose(ax, center=(-77.1, 35.87), size=0.05)

# Add distance marker
draw_distance_marker(ax, start=(-77.15, 35.73), distance_km=10, text='10 km')

# Add gridlines
ax.gridlines(draw_labels=True)

# Title
ax.set_title('Albemarle Sound with Sensor Locations', fontsize=15)

# Add a legend
ax.legend()

# Show the plot
plt.show()



# %% 
# same thing but for Stockton area
import matplotlib.pyplot as plt
import geopandas as gpd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd

# Load the GeoJSON file
shapefile_path = "C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Code Repositories\\R\\Stockton\\stockton_shp"

# Load the shapefile
Stockton = gpd.read_file(shapefile_path)

# Load additional shapefiles for small rivers and lakes
lakes_shapefile_path = "C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Code Repositories\\R\\Stockton\\ponds\\Ponds.shp"
rivers_shapefile_path = "C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Code Repositories\\R\\Stockton\\waterlines\\WaterLines.shp"
rivers = gpd.read_file(rivers_shapefile_path)
lakes = gpd.read_file(lakes_shapefile_path)

# Read the sensor data from the provided text file
sensor_data_path = "C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Code Repositories\\R\\Stockton\\points.xlsx"
sensor_data = pd.read_excel(sensor_data_path)

# Extract relevant columns (longitude and latitude)
sensor_data = sensor_data[['Lon', 'Lat']]

# Create the plot
fig, ax = plt.subplots(1, 1, figsize=(10, 10), subplot_kw={'projection': ccrs.PlateCarree()})

# Add features to the map
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.LAND)

# Customize the OCEAN color
ocean = cfeature.NaturalEarthFeature(
    'physical', 'ocean', '10m',
    edgecolor='face',
    facecolor=cfeature.COLORS['water']
)
ax.add_feature(ocean, facecolor='#6699CC', alpha = 0.8)  # Change 'lightblue' to your desired water color

# Plot the Albemarle Sound
Stockton.plot(ax=ax, edgecolor='#6699CC') # facecolor='lightblue'

# Plot the small rivers and lakes
rivers.plot(ax=ax, color='blue', linewidth=1, zorder=1)  # Adjust color and linewidth as needed
lakes.plot(ax=ax, color='lightgreen', edgecolor='black', zorder=2)  # Adjust color and edgecolor as needed

# Set the extent to focus on the Albemarle Sound region
ax.set_extent([-122.2, -121.2, 37.7, 38.7], crs=ccrs.PlateCarree())

# Plot the sensor locations
ax.scatter(sensor_data['Lon'], sensor_data['Lat'], color='black', s=100, zorder=5, transform=ccrs.PlateCarree(), label='Sampling Sites')

# Add gridlines
ax.gridlines(draw_labels=True)

# Title
ax.set_title('Stockton with Sample Sites', fontsize=15)

# Add a legend
ax.legend()

# Show the plot
plt.show()
# %%
