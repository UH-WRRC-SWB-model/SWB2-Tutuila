""" # Tutuila SWB2 model code
by Chris Shuler and Aly El-Kadi

This code was written in python (version 3.6.2, |Continuum Analytics, Inc.| (default, Jul 20 2017, 12:30:02) [MSC v.1900 64 bit (AMD64)]), and specifically in Jupyter notebook (version 5.6.0) 
running in an Anaconda environment. Additionally, use of ArcPy modules from ArcPro 
version 2.1.0 were used for geoprocessing operations.

To run this code, the user should have all of the modules specified in the first 
code block below installed. If any modules are not installed on your computer, 
simply google the module name followed by "anaconda forge" to obtain the syntax 
that will allow you to download the missing package(s) using the command prompt or 
other package installation method.

"""


#### import modules
import numpy as np
import arcpy
import os
import sys
from arcpy.sa import *
import pandas as pd
import gdal
from arcpy import env
import shutil
import numpy.ma as ma
import netCDF4 as nc
import subprocess
import re
%matplotlib notebook
# this is a list of additional functions to load up, as to not clutter the script
%run ../../Std_input/COMMON/plot_and_table_functions

# set properties
arcpy.env.overwriteOutput = True # make sure overwrite files is on
# projection definition 
sr_project = arcpy.SpatialReference(32702)   # Project dataset into WGS84
cel_size = 20     # in m 
Control_File_Name = 'Tutuila200_controlFile.ctl'


#### General coverages and paths. More, basic model setup.
GIS_FOLDER = os.path.join('..', '..', 'Raw_GIS_Data')
STD_INPUT_FOLDER = os.path.join('..', '..', 'Std_input')
# path to the grid bound
Grid_shp = os.path.join(GIS_FOLDER, 'grid_bound.shp')

if not os.path.exists(os.path.join('..', 'output//')):
    os.makedirs(os.path.join('..', 'output//'))
                        
if not os.path.exists(os.path.join('..', 'output//Figures//')):
    os.makedirs(os.path.join('..', 'output//Figures//'))
fig_path =  (os.path.join('..', 'output//Figures//'))

# set/create General workspace
workspace = os.path.join('..', 'input/General')
if not os.path.exists(workspace):
    os.makedirs(workspace)

# this makes a domain grid of zeros of with cell size
Grid_rast0 =os.path.join(workspace, 'Grid_rast0')
arcpy.PolygonToRaster_conversion(os.path.join(GIS_FOLDER, 'grid_bound.shp'), "zero", Grid_rast0, cell_assignment="MAXIMUM_AREA",  cellsize=cel_size )

# this makes a domain grid of ones of with cell size
Grid_rast1 =os.path.join(workspace, 'Grid_rast1')
arcpy.PolygonToRaster_conversion(os.path.join(GIS_FOLDER, 'grid_bound.shp'), "one", Grid_rast1, cell_assignment="MAXIMUM_AREA",  cellsize=cel_size )


# The next code blocks process shapefile based (vector) input data
### Land use data: incluing canopy coverage and % impervious surfaces 
# set/create workspace
workspace = os.path.join('..', 'input/Land_use_Soils_data')
if not os.path.exists(workspace):
    os.makedirs(workspace)

# project file 
arcpy.Project_management(os.path.join(GIS_FOLDER, "Land_use\\Land_use_wRO_codes2.shp"), os.path.join(workspace, 'LU_shp_projected.shp'), sr_project)

#  Merge in the grid bound into the shapefile to create accurate grid coverage 
arcpy.Erase_analysis(Grid_shp ,os.path.join(workspace, 'LU_shp_projected.shp'),  os.path.join(workspace, 'LU_shp_bound.shp'))
arcpy.Merge_management([os.path.join(workspace, 'LU_shp_bound.shp'), os.path.join(workspace, 'LU_shp_projected.shp')], os.path.join(workspace, 'LU_shp_ready.shp'))

# convert to raster.asc
arcpy.PolygonToRaster_conversion(os.path.join(workspace, 'LU_shp_ready.shp'), "LU2", os.path.join(workspace, "LU_raster"), cell_assignment="MAXIMUM_AREA",  cellsize=cel_size)
arcpy.RasterToASCII_conversion(os.path.join(workspace, "LU_raster"), os.path.join(workspace, "LU_grid.asc"))

# canopy cover %ages with landuse map
arcpy.PolygonToRaster_conversion(os.path.join(workspace, 'LU_shp_ready.shp'), "fracCanCov", os.path.join(workspace, "CanCovRas"), cell_assignment="MAXIMUM_AREA",  cellsize=cel_size)
arcpy.RasterToASCII_conversion(os.path.join(workspace, "CanCovRas"), os.path.join(workspace, "CanCovRas.asc"))

# % impervious raster from landuse map
arcpy.PolygonToRaster_conversion(os.path.join(workspace, 'LU_shp_ready.shp'), "pct_Impv", os.path.join(workspace, "pImpvRas"), cell_assignment="MAXIMUM_AREA",  cellsize=cel_size )
arcpy.RasterToASCII_conversion(os.path.join(workspace, "pImpvRas"), os.path.join(workspace, "pImpvRas.asc"))

# % pervious raster from landuse map
arcpy.PolygonToRaster_conversion(os.path.join(workspace, 'LU_shp_ready.shp'), "pct_Perv", os.path.join(workspace, "pPERVRas"), cell_assignment="MAXIMUM_COMBINED_AREA",  cellsize=cel_size )
arcpy.RasterToASCII_conversion(os.path.join(workspace, "pPERVRas"), os.path.join(workspace, "pPERVRas.asc"))

# clean up workspace
arcpy.Delete_management(os.path.join(workspace, 'LU_shp_projected.shp'))
arcpy.Delete_management(os.path.join(workspace, 'LU_shp_bound.shp'))
arcpy.Delete_management(os.path.join(workspace, 'LU_shp_ready.shp'))
arcpy.Delete_management(os.path.join(workspace, 'LU_raster'))
arcpy.Delete_management(os.path.join(workspace, 'CanCovRas'))
arcpy.Delete_management(os.path.join(workspace, 'pImpvRas'))
arcpy.Delete_management(os.path.join(workspace, 'pPERVRas'))

# plot output 
lu_data, ly_gt, lu_proj, lu_xy                = read_raster( os.path.join(workspace, "LU_grid.asc") )  # land use
lu_cmap = discrete_irreg_cmap(discrete_vals=np.unique( lu_data.flatten()), base_cmap='nipy_spectral')
make_plot( x=lu_xy[0], y=lu_xy[1], var=lu_data, discrete=True, cmap=lu_cmap, title='Land-use type codes' )
plt.savefig(os.path.join(fig_path, "Land_use_{}m.png".format(cel_size)), dpi = 300)

CC_data, CC_gt, CC_proj, CC_xy                = read_raster( os.path.join(workspace, "CanCovRas.asc"))            # canopy cover
make_plot( x=CC_xy[0], y=CC_xy[1], var=CC_data, title='Percent Canopy Cover', barlabel='% canopy cover')
plt.savefig(os.path.join(fig_path, "Canopy_cover_{}m.png".format(cel_size)), dpi = 300)

PI_data, PI_gt, PI_proj, PI_xy                = read_raster( os.path.join(workspace, "pImpvRas.asc"))
make_plot( x=PI_xy[0], y=PI_xy[1], var=PI_data, title='Percent Impervious surface', barlabel='% impervious') #, savepath = fig_path)
plt.savefig(os.path.join(fig_path, "Percent_impervious_{}m.png".format(cel_size)), dpi = 300)

### Same with soil shapefile into WGS84 and then resterize to make soil moisture and h2o group rasters
# set/create workspace
workspace = os.path.join('..', 'input/Land_use_Soils_data')
if not os.path.exists(workspace):
    os.makedirs(workspace)

# Project dataset into WGS84
arcpy.Project_management(os.path.join(GIS_FOLDER, 'Soils\\Tut_Soil_clip2.shp'), os.path.join(workspace, 'soils_shp_projected.shp'), sr_project) 

arcpy.Erase_analysis(Grid_shp, os.path.join(workspace, 'soils_shp_projected.shp'),  os.path.join(workspace, 'soils_shp_bound.shp'))
arcpy.Merge_management([os.path.join(workspace, 'soils_shp_bound.shp'), os.path.join(workspace, 'soils_shp_projected.shp')], os.path.join(workspace, 'soils_shp_ready.shp'))

# convert to raster
arcpy.PolygonToRaster_conversion(os.path.join(workspace, 'soils_shp_ready.shp'),  "H2O_grp", os.path.join(workspace, "H2Ogp_rst"), cell_assignment="MAXIMUM_AREA",  cellsize=cel_size)
arcpy.PolygonToRaster_conversion(os.path.join(workspace, 'soils_shp_ready.shp'),  "SM_INpFT", os.path.join(workspace, "SMC_rst"), cell_assignment="MAXIMUM_AREA",  cellsize=cel_size)

arcpy.RasterToASCII_conversion(os.path.join(workspace, "SMC_rst"), os.path.join(workspace, "soil_moist_cap_grid.asc"))
arcpy.RasterToASCII_conversion(os.path.join(workspace, "H2Ogp_rst"), os.path.join(workspace, "H2Ogp_grid.asc"))

# clean up workspace
arcpy.Delete_management(os.path.join(workspace, 'soils_shp_projected.shp'))
arcpy.Delete_management(os.path.join(workspace, 'soils_shp_bound.shp'))
arcpy.Delete_management(os.path.join(workspace, 'soils_shp_ready.shp'))
arcpy.Delete_management(os.path.join(workspace, 'SMC_rst'))
arcpy.Delete_management(os.path.join(workspace, 'H2Ogp_rst'))

# plot output 
soils_data, soils_gt, soils_proj, soils_xy    = read_raster( os.path.join(workspace, "H2Ogp_grid.asc"))    
make_plot( x=soils_xy[0], y=soils_xy[1], var=soils_data, discrete=True, title='Hydrologic soil group', barlabel='Soil group number')
plt.savefig(os.path.join(fig_path, "soil_hydro_group{}m.png".format(cel_size)), dpi = 300)

ss_data, ss_gt, ss_proj, ss_xy                = read_raster( os.path.join(workspace, "soil_moist_cap_grid.asc") )  # soil moinsture capacity
make_plot( x=ss_xy[0], y=ss_xy[1], var=ss_data, title='Maximum Soil Storage', barlabel='Soil Storage, in inches')
plt.savefig(os.path.join(fig_path, "soil_moisture_capac_{}m.png".format(cel_size)), dpi = 300)

### Rasterize the rainfall station zones
workspace = os.path.join('..', 'input/Rain_stations')
if not os.path.exists(workspace):
    os.makedirs(workspace)

# project file 
arcpy.Project_management(os.path.join(GIS_FOLDER, 'Stations\\Thissen_poly_rain_clip_modified2.shp'), os.path.join(workspace, 'Thissen_poly_projected.shp'), sr_project)

#  Merge in the grid bound into the shapefile to create accurate grid coverage 
arcpy.Erase_analysis(Grid_shp ,os.path.join(workspace, 'Thissen_poly_projected.shp'),  os.path.join(workspace, 'Thissen_poly_bound.shp'))
arcpy.Merge_management([os.path.join(workspace, 'Thissen_poly_bound.shp'), os.path.join(workspace, 'Thissen_poly_projected.shp')], os.path.join(workspace, 'Thissen_poly_ready.shp'))

# convert to raster.asc
arcpy.PolygonToRaster_conversion(os.path.join(workspace, 'Thissen_poly_ready.shp'), "Gage_ID", os.path.join(workspace, "TP_raster"), cell_assignment="MAXIMUM_AREA",  cellsize=cel_size)
arcpy.RasterToASCII_conversion(os.path.join(workspace, "TP_raster"), os.path.join(workspace, "TP_grid.asc"))

# clean up workspace
arcpy.Delete_management(os.path.join(workspace, 'Thissen_poly_projected.shp'))
arcpy.Delete_management(os.path.join(workspace, 'Thissen_poly_bound.shp'))
arcpy.Delete_management(os.path.join(workspace, 'Thissen_poly_ready.shp'))
arcpy.Delete_management(os.path.join(workspace, 'TP_raster'))

rfs_data, rfs_gt, rfs_proj, rfs_xy    = read_raster( os.path.join(workspace, "TP_grid.asc"))    
make_plot( x=rfs_xy[0], y=rfs_xy[1], var=rfs_data, discrete=True, title='Rainfall station zones', barlabel='Zone ID')
plt.savefig(os.path.join(fig_path, "Rainfall_stations{}m.png".format(cel_size)), dpi = 300)

### Rasterise the runoff zones (runoff to rainfall watersheds)
workspace = os.path.join('..', 'input/Runoff_zones_WS')
if not os.path.exists(workspace):
    os.makedirs(workspace)
    
# project file 
arcpy.Project_management(os.path.join(GIS_FOLDER, 'Runofftorainfall2\\All_major_WS_modified3.shp'), os.path.join(workspace, 'RO_projected.shp'), sr_project)
arcpy.Erase_analysis(Grid_shp ,os.path.join(workspace, 'RO_projected.shp'),  os.path.join(workspace, 'RO_bound.shp'))
arcpy.Merge_management([os.path.join(workspace, 'RO_bound.shp'), os.path.join(workspace, 'RO_projected.shp')], os.path.join(workspace, 'RO_ready.shp'))

arcpy.PolygonToRaster_conversion(os.path.join(workspace, 'RO_ready.shp'), "Zone_ID", os.path.join(workspace, "Ro_Rast"), cell_assignment="MAXIMUM_AREA",  cellsize=cel_size)
arcpy.RasterToASCII_conversion(os.path.join(workspace, "Ro_Rast"), os.path.join(workspace, "Ro_Rast.asc"))

# clean up workspace
arcpy.Delete_management(os.path.join(workspace, 'RO_projected.shp'))
arcpy.Delete_management(os.path.join(workspace, 'RO_bound.shp'))
arcpy.Delete_management(os.path.join(workspace, 'RO_ready.shp'))
arcpy.Delete_management(os.path.join(workspace, 'Ro_Rast'))

rfs_data, rfs_gt, rfs_proj, rfs_xy    = read_raster( os.path.join(workspace, "Ro_Rast.asc"))    
make_plot( x=rfs_xy[0], y=rfs_xy[1], var=rfs_data, discrete=True, title='Runoff zones', barlabel='Zone ID')
plt.savefig(os.path.join(fig_path, "Runoff_zones{}m.png".format(cel_size)), dpi = 300)


"""  note that the UH-ASPA data value of 0.4115 for Fagaalu was replaced by the Olkeba-Alex value of 0.285015 because the Falu stream gauge was out of the water a good bit"""

### Rasterize the water lines to calculate direct infiltration from leaking lines
workspace = os.path.join('..', 'input/Direct_infiltration')
if not os.path.exists(workspace):
    os.makedirs(workspace)
    
# Calculate in the amount of NRW water equally distributed over the area of the water lines 
# note total ASPA 5 year production average = 12,740,000 and the Non-Rev water ave = 7,940,000 Gal/day   (264.172 gal per m3 of water) 
# conversions 5 year production average in m3/day = 48226.1557 and the Non-Rev water ave = 30056.175522 m3/day    
    
# create appropriately sized buffers of 1/2 the cell size
arcpy.Buffer_analysis(os.path.join(GIS_FOLDER, 'Direct_infiltration\\Transmission_water_mains.shp'), os.path.join(workspace, 'Water_lines_Buffer.shp'), "{} meters".format(cel_size/2), "FULL", "ROUND", "ALL")
# calculate area of water line polygons in m2
arcpy.AddField_management(os.path.join(workspace, 'Water_lines_Buffer.shp'), "Infl_inch", "DOUBLE")    # add Active cell unit field
arcpy.AddGeometryAttributes_management(os.path.join(workspace, 'Water_lines_Buffer.shp'), "AREA")    # calculate the area of the buffer zone
Area_list= [row[0] for row in arcpy.da.SearchCursor(os.path.join(workspace, 'Water_lines_Buffer.shp'), "POLY_AREA")]  # make a list of the area number
arcpy.CalculateField_management(os.path.join(workspace, 'Water_lines_Buffer.shp'), "Infl_inch", "!Id! + 30056.1755/Area_list[0]*39.3701", "PYTHON3") # calculate the appropriate amount of infitration in inches to equal 7.94 MGD as spread over the whole water line influence area

    
# Project dataset into WGS84
arcpy.Project_management(os.path.join(workspace, 'Water_lines_Buffer.shp'), os.path.join(workspace, 'Water_lines_Buffer_shp_projected.shp'), sr_project) 
arcpy.Erase_analysis(Grid_shp, os.path.join(workspace, 'Water_lines_Buffer_shp_projected.shp'),  os.path.join(workspace, 'Water_lines_Buffer_shp_bound.shp'))
arcpy.Merge_management([os.path.join(workspace, 'Water_lines_Buffer_shp_bound.shp'), os.path.join(workspace, 'Water_lines_Buffer_shp_projected.shp')], os.path.join(workspace, 'Water_lines_Buffer_shp_ready.shp'))

arcpy.PolygonToRaster_conversion(os.path.join(workspace, 'Water_lines_Buffer_shp_ready.shp'), "Infl_inch", os.path.join(workspace, "A_WL_Rast"), cell_assignment="MAXIMUM_AREA",  cellsize=cel_size)
arcpy.RasterToASCII_conversion(os.path.join(workspace, "A_WL_Rast"), os.path.join(workspace, "A_WL_Rast.asc"))

# clean up workspace
arcpy.Delete_management(os.path.join(workspace, 'Water_lines_Buffer.shp'))
arcpy.Delete_management(os.path.join(workspace, 'Water_lines_Buffer_shp_projected.shp'))
arcpy.Delete_management(os.path.join(workspace, 'Water_lines_Buffer_shp_bound.shp'))
arcpy.Delete_management(os.path.join(workspace, 'Water_lines_Buffer_shp_ready.shp'))
arcpy.Delete_management(os.path.join(workspace, 'A_WL_Rast'))

Pcmap = discrete_irreg_cmap(discrete_vals=np.unique( lu_data.flatten()), base_cmap='Greys')   # base_cmap='Greys' gist_gray
rfs_data, rfs_gt, rfs_proj, rfs_xy    = read_raster( os.path.join(workspace, "A_WL_Rast.asc")) 
make_plot( x=rfs_xy[0], y=rfs_xy[1], var=rfs_data, discrete=True, cmap=Pcmap, title='Water line direct infiltration', barlabel='direct infiltration [in]')
plt.savefig(os.path.join(fig_path, "Water_line_infiltration{}m.png".format(cel_size)), dpi = 300)

### Rasterize OSDS location infiltration
# FROM Shuler et al 2017   OSDS flow = 1.454203 m3/unit/Day
workspace = os.path.join('..', 'input/Direct_infiltration')
if not os.path.exists(workspace):
    os.makedirs(workspace)
    
arcpy.env.extent = os.path.join('..', 'input/General/Grid_rast1')  # set processing extent to our desired model boundary
    
pdensOut = arcpy.sa.PointDensity(os.path.join(GIS_FOLDER, 'Direct_infiltration\\OSDS_units_pts.shp'),"NONE", cel_size, NbrCircle(1, "CELL"), "SQUARE_METERS" )
pdensOut.save(os.path.join(workspace, 'temp_dens'))   # calculate the OSDS density per m2

arcpy.Times_3d(os.path.join(workspace, 'temp_dens'), cel_size**2, os.path.join(workspace, 'OSDS_dens'))  # make the OSDS density in OSDS units, meaning number of units per cell

infiltration_in_inches_per_day_per_OSDS = (1.454203/cel_size**2)*39.3701

arcpy.Times_3d(os.path.join(workspace, 'OSDS_dens'), infiltration_in_inches_per_day_per_OSDS, os.path.join(workspace, 'OSDS_inlf_in'))  # calculate amount of water infiltrated in each cell in inches

arcpy.RasterToASCII_conversion(os.path.join(workspace, 'OSDS_inlf_in'), os.path.join(workspace, "OSDS_inlf_in.asc"))

arcpy.Delete_management(os.path.join(workspace, 'temp_dens'))   # clean up workspace
arcpy.Delete_management(os.path.join(workspace, 'OSDS_inlf_in'))   # clean up workspace
arcpy.Delete_management(os.path.join(workspace, 'OSDS_inlf_in'))   # clean up workspace
arcpy.Delete_management(os.path.join(workspace, 'OSDS_dens'))   # clean up workspace

# Plot
Pcmap = discrete_irreg_cmap(discrete_vals=np.unique( lu_data.flatten()), base_cmap='Greys')   # base_cmap='Greys'   gist_gray
rfs_data, rfs_gt, rfs_proj, rfs_xy    = read_raster( os.path.join(workspace, "OSDS_inlf_in.asc")) 
make_plot( x=rfs_xy[0], y=rfs_xy[1], var=rfs_data, discrete=True, cmap=Pcmap, title='OSDS direct infiltration', barlabel='direct infiltration [in]')
plt.savefig(os.path.join(fig_path, "OSDS_infiltration{}m.png".format(cel_size)), dpi = 300)

#### add all the direct infiltration into one raster because apparently the control file wont take multiples???
arcpy.Plus_3d(os.path.join('..', 'input/Direct_infiltration',  "A_WL_Rast.asc"), os.path.join('..', 'input/Direct_infiltration', "OSDS_inlf_in.asc"), os.path.join('..', 'input/Direct_infiltration', "WLOSDrast"))
arcpy.RasterToASCII_conversion(os.path.join('..', 'input/Direct_infiltration', "WLOSDrast"), os.path.join('..', 'input/Direct_infiltration', "Total_inlf_in.asc"))


### Create Flow direction raster
####resampling the 10 m DEM to whatever cell size we is using 
workspace = os.path.join('..', 'input/DEM_Process')
if not os.path.exists(workspace):
    os.makedirs(workspace)
    
# this resamples the 10 , dem to whatever cell size we are working in
arcpy.ProjectRaster_management(os.path.join( GIS_FOLDER, 'DEM', '10M_DEM.tif'), os.path.join(workspace, 'DEM_projected.tif'), sr_project, cell_size = cel_size )
arcpy.Clip_management(os.path.join(workspace, 'DEM_projected.tif'), "515000 8429000 550000 8409000", os.path.join(workspace, 'DEM_clip.tif'), Grid_shp, -9999, 'ClippingGeometry', 'MAINTAIN_EXTENT')
arcpy.Plus_3d(os.path.join(workspace, 'DEM_clip.tif'), Grid_rast0,  os.path.join(workspace, 'DEM_ready.tif'))

outFlowDirection = FlowDirection(os.path.join(workspace, 'DEM_ready.tif'), "NORMAL")
outFlowDirection.save(os.path.join(workspace,'Flow_direction.tif'))
arcpy.RasterToASCII_conversion(os.path.join(workspace,'Flow_direction.tif'), os.path.join(workspace,'Flow_direction.asc'))

arcpy.Delete_management(os.path.join(workspace, 'DEM_projected.tif'))
arcpy.Delete_management(os.path.join(workspace, 'DEM_clip.tif'))
arcpy.Delete_management(os.path.join(workspace, 'DEM_ready.tif'))
arcpy.Delete_management(os.path.join(workspace, 'Flow_direction.tif'))


rfs_data, rfs_gt, rfs_proj, rfs_xy    = read_raster( os.path.join(workspace, "Flow_direction.asc")) 
make_plot( x=rfs_xy[0], y=rfs_xy[1], var=rfs_data, discrete=True,  title='D8 flow direction', barlabel='flow direction code')
plt.savefig(os.path.join(fig_path, "Flow_direction{}m.png".format(cel_size)), dpi = 300)



### Create mean annual wet canopy evaporation layer from wind and rainfall data
workspace = os.path.join('..', 'input/Evaporation')
if not os.path.exists(workspace):
    os.makedirs(workspace)
#files to use    
anual_rain = os.path.join(GIS_FOLDER, 'Gridded_rain', "year", "An_pcip_in.tif")
anual_wind =  os.path.join(GIS_FOLDER, 'Wind_map', 'wndsp_30m.tif')

# divide them to get w
outDivide = Divide(anual_wind, anual_rain)                     
outDivide.save(os.path.join(workspace, "w_WndOvrRain.tif")) # create variable w = quotient of mean annual wind speed, in meters per second, and mean annual rainfall, in inches

# now apply equation 18  V = 2.677 × (w) – 0.014, from "Spatially Distributed Groundwater Recharge for 2010 Land Cover Estimated Using a Water-Budget Model for the Island of O’ahu, Hawai’i"
outTimes = Times(os.path.join(workspace, "w_WndOvrRain.tif"), 2.677)                     
outTimes.save(os.path.join(workspace, "w_Wndtmp.tif")) # create variable w = quotient of mean annual wind speed, in meters per second, and mean annual rainfall, in inches
arcpy.Plus_3d(os.path.join(workspace, "w_Wndtmp.tif"),  -0.014, os.path.join(workspace, "V_evp2pcip.tif"))

# synchronize geometry with the model
arcpy.ProjectRaster_management(os.path.join(workspace, "V_evp2pcip.tif"), os.path.join(workspace, 'V_evp_projected.tif'), sr_project, cell_size = cel_size )
arcpy.Clip_management(os.path.join(workspace, 'V_evp_projected.tif'), "515000 8429000 550000 8409000", os.path.join(workspace, 'V_evp_clip.tif'), Grid_shp, -9999, 'ClippingGeometry', 'MAINTAIN_EXTENT')
arcpy.Plus_3d(os.path.join(workspace, 'V_evp_clip.tif'), Grid_rast0,  os.path.join(workspace, 'V_evp_ready.tif'))
# print .asc file
arcpy.RasterToASCII_conversion(os.path.join(workspace,'V_evp_ready.tif'), os.path.join(workspace,'V_evp2pcp.asc'))   

#clean
arcpy.Delete_management(os.path.join(workspace, 'w_WndOvrRain.tif'))
arcpy.Delete_management(os.path.join(workspace, 'w_Wndtmp.tif'))
arcpy.Delete_management(os.path.join(workspace, 'V_evp2pcip.tif'))
arcpy.Delete_management(os.path.join(workspace, 'V_evp_projected.tif'))
arcpy.Delete_management(os.path.join(workspace, 'V_evp_clip.tif'))
arcpy.Delete_management(os.path.join(workspace, 'V_evp_ready.tif'))

rfs_data, rfs_gt, rfs_proj, rfs_xy    = read_raster( os.path.join(workspace, "V_evp2pcp.asc")) 
make_plot( x=rfs_xy[0], y=rfs_xy[1], var=rfs_data, tick_factor=0.1, discrete=True,  title='Canopy Evaporation', barlabel='Evaporation [in.]')
plt.savefig(os.path.join(fig_path, "Canopy_evaporation{}m.png".format(cel_size)), dpi = 300)
    

### format raingall grids
workspace = os.path.join('..', 'input/Gridded_rain')
if not os.path.exists(workspace):
    os.makedirs(workspace)

mo_list = ["01", "02", "03","04", "05", "06","07", "08", "09","10", "11", "12"]

for i in mo_list:
    inpt =  os.path.join(GIS_FOLDER, 'Gridded_rain', "PRISM_ppt_tutuila_30yr_normal_80mM1_{}_asc.asc".format(i))
    outpy = os.path.join(workspace, 'PRISM_ppt_clip.tif')

    arcpy.ProjectRaster_management(inpt, outpy, sr_project, "BILINEAR", cel_size, "NAD_1983_To_WGS_1984_1", "#", "#") 
    outpy2 = Times(outpy, 0.000393701)   # needed to convert 100* mm to inches
    arcpy.Clip_management(outpy2, "515000 8429000 550000 8409000", os.path.join(workspace, 'file_clip.tif'), Grid_shp, -9999, 'ClippingGeometry', 'MAINTAIN_EXTENT')
    arcpy.Plus_3d(os.path.join(workspace, 'file_clip.tif'), Grid_rast0,  os.path.join(workspace, 'file_plus.tif'))    
    
    arcpy.RasterToASCII_conversion(os.path.join(workspace, 'file_plus.tif'), os.path.join(workspace, 'PRISM_ppt_tutuila_30yr_normal_{}.asc'.format(i)))
    
arcpy.Delete_management(os.path.join(workspace, 'PRISM_ppt_clip.tif'))
arcpy.Delete_management(os.path.join(workspace, 'file_clip.tif'))
arcpy.Delete_management(os.path.join(workspace, 'file_plus.tif'))

rfs_data, rfs_gt, rfs_proj, rfs_xy    = read_raster( os.path.join(workspace, 'PRISM_ppt_tutuila_30yr_normal_{}.asc'.format('12'))) 
make_plot( x=rfs_xy[0], y=rfs_xy[1], var=rfs_data, discrete=True,  title='Rainfall month_{}'.format("12"), barlabel='precipitation [in.]')
plt.savefig(os.path.join(fig_path, "Rainfall{}m.png".format(cel_size)), dpi = 300)


### format Air Temp grids
workspace = os.path.join('..', 'input/Gridded_temps')
if not os.path.exists(workspace):
    os.makedirs(workspace)

mo_list = ["01", "02", "03","04", "05", "06","07", "08", "09","10", "11", "12"]

#max temps
for i in mo_list:
    inpt =  os.path.join(GIS_FOLDER, 'Gridded_temps/Temp_max/Monthly', "PRISM_tmax_tutuila_30yr_normal_80mM1_{}_asc.asc".format(i))
    outpy = os.path.join(workspace,'PRISM_clip.tif')

    arcpy.ProjectRaster_management(inpt, outpy, sr_project, "BILINEAR", cel_size, "NAD_1983_To_WGS_1984_1", "#", "#") 
    outpy2 = Times(outpy, 0.01)   # needed to convert to deg celsius
    arcpy.Clip_management(outpy2, "515000 8429000 550000 8409000", os.path.join(workspace, 'file_clip.tif'), Grid_shp, -9999, 'ClippingGeometry', 'MAINTAIN_EXTENT')
    arcpy.Plus_3d(os.path.join(workspace, 'file_clip.tif'), Grid_rast0,  os.path.join(workspace, 'file_plus.tif'))    
    
    arcpy.RasterToASCII_conversion(os.path.join(workspace, 'file_plus.tif'), os.path.join(workspace, 'PRISM_tmax_tutuila_30yr_normal_{}.asc'.format(i)))
    
# min temps
for i in mo_list:
    inpt =  os.path.join(GIS_FOLDER, 'Gridded_temps/Temp_min/Monthly', "PRISM_tmin_tutuila_30yr_normal_80mM1_{}_asc.asc".format(i))
    outpy = os.path.join(workspace, 'PRISM_clip.tif')

    arcpy.ProjectRaster_management(inpt, outpy, sr_project, "BILINEAR", cel_size, "NAD_1983_To_WGS_1984_1", "#", "#") 
    outpy2 = Times(outpy, 0.01)   # needed to convert to deg celsius
    arcpy.Clip_management(outpy2, "515000 8429000 550000 8409000", os.path.join(workspace,'file_clip.tif'), Grid_shp, -9999, 'ClippingGeometry', 'MAINTAIN_EXTENT')
    arcpy.Plus_3d(os.path.join(workspace, 'file_clip.tif'), Grid_rast0,  os.path.join(workspace, 'file_plus.tif'))    
    
    arcpy.RasterToASCII_conversion(os.path.join(workspace, 'file_plus.tif'), os.path.join(workspace, 'PRISM_tmin_tutuila_30yr_normal_{}.asc'.format(i)))    
    
arcpy.Delete_management(os.path.join(workspace, 'PRISM_clip.tif'))
arcpy.Delete_management(os.path.join(workspace, 'file_clip.tif'))
arcpy.Delete_management(os.path.join(workspace, 'file_plus.tif'))

rfs_data, rfs_gt, rfs_proj, rfs_xy    = read_raster( os.path.join(workspace, 'prism_tmin_tutuila_30yr_normal_12.asc'))
make_plot( x=rfs_xy[0], y=rfs_xy[1], var=rfs_data, discrete=True,  title='Air, min Temp month_{}'.format("12"), barlabel='Temperature [F]')
plt.savefig(os.path.join(fig_path, "Min_temp{}m.png".format(cel_size)), dpi = 300)


### Create Monthly Evapotranspiration grids from the original Izuka 2005 paper maps
#Spline method 
workspace = os.path.join('..', 'input//ET_Process')
if not os.path.exists(workspace):
    os.makedirs(workspace)

files = ['Jan_ET.shp', 'Feb_ET.shp', 'March_ET.shp', 'April_ET.shp', 'May_ET.shp', 'June_ET.shp', 'July_ET.shp', 'Aug_ET.shp', 'Sep_ET.shp', 'Oct_ET.shp', 'Nov-Dec_ET.shp', 'Nov-Dec_ET.shp']
mos = ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec",]

for idx, i in enumerate(files):
    key_nam = mos[idx]
    outSPL = Spline(os.path.join(GIS_FOLDER, 'ET', i), "ET_inches", cel_size, "REGULARIZED")   # spline the points
    outSPL.save(os.path.join(workspace, 'ET.tif'))
    
    outExtractByMask = ExtractByMask(outSPL, os.path.join(GIS_FOLDER, 'Precip_bound_poly.shp'))                  # clip the raster
    outExtractByMask.save(os.path.join(workspace, 'ET_clipped.tif'))
    
# already projected   arcpy.ProjectRaster_management(os.path.join(workspace, 'junk' 'ET_clipped.tif'), os.path.join(workspace, 'junk', 'ET_proj.tif'), sr_project, cell_size = cel_size)
    arcpy.Clip_management(os.path.join(workspace, 'ET_clipped.tif'), "515000 8429000 550000 8409000", os.path.join(workspace, 'file_clip.tif'), Grid_shp, -9999, 'ClippingGeometry', 'MAINTAIN_EXTENT')   # make sure raster extent is good
    arcpy.Plus_3d(os.path.join(workspace, 'file_clip.tif'), Grid_rast0,  os.path.join(workspace, 'file_plus.tif'))
    
    arcpy.RasterToASCII_conversion(os.path.join(workspace, 'file_plus.tif'), os.path.join(workspace, key_nam+'_ET_clipped.asc'))    

    arcpy.Delete_management(os.path.join(workspace, 'ET.tif'))                                                     # delete extranious file
    arcpy.Delete_management(os.path.join(workspace, 'ET_clipped.tif'))    
    arcpy.Delete_management(os.path.join(workspace, 'file_clip.tif'))    
    arcpy.Delete_management(os.path.join(workspace, 'file_plus.tif'))  
    
rfs_data, rfs_gt, rfs_proj, rfs_xy    = read_raster( os.path.join(workspace, 'dec_et_clipped.asc'))
make_plot( x=rfs_xy[0], y=rfs_xy[1], var=rfs_data, discrete=True,  title='Potental ET, December', barlabel='PET [in.]')
plt.savefig(os.path.join(fig_path, "ET_dec{}m.png".format(cel_size)), dpi = 300)



### Silly rainfall adjustment grids,  not actual data, just a filler grid
workspace = os.path.join('..', 'input//RF_adj_grids')
if not os.path.exists(workspace):
    os.makedirs(workspace)

Field_to_rasterize = "one"

mos = ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec",]

for i in mos:
    Gen_name = "RFadj_{}".format(i)

    arcpy.PolygonToRaster_conversion(Grid_shp, Field_to_rasterize, os.path.join(workspace, '{}'.format(Gen_name)), cell_assignment="MAXIMUM_AREA",  cellsize=cel_size)
    arcpy.RasterToASCII_conversion( os.path.join(workspace, '{}'.format(Gen_name)), os.path.join(workspace, "{}.asc".format(Gen_name)))


# Move in other standard input files and modify control file to the shape of the current run
# modify the control file grid for the given run   (note this uses dimensions from the rainfall adjustment grid in april)
with open(os.path.join('..', 'input//RF_adj_grids', 'rfadj_apr.asc'), 'r') as dims_file:   # open an ASC file and get the dimenstions out of it 
    dimsfile1 = dims_file.read().splitlines(True)
    x_dim = float(re.findall('\d+', dimsfile1[0])[-1])    
    y_dim = float(re.findall('\d+', dimsfile1[1])[-1])

with open(os.path.join('.', Control_File_Name), 'r') as fin:   # open file 
    data = fin.read().splitlines(True)
with open(os.path.join('.', Control_File_Name), 'w') as fout:     # delete first line
    fout.writelines(data[1:])
new_first = 'GRID {} {} 515000. 8409000. {} '.format(x_dim, y_dim, cel_size)  # new first line 
with open(os.path.join('.', Control_File_Name), 'r+') as file:                # add in new first line and save file  
    file_data = file.read()
    file. seek(0, 0)
    file. write(new_first + '\n' + file_data)

# land use lookup file
shutil.copy2(os.path.join(STD_INPUT_FOLDER, 'Landuse_lookup_maui_mod5.txt') ,os.path.join('..', 'input'))    

# Simple RO : RF ratios file
shutil.copy2(os.path.join(GIS_FOLDER, 'Runofftorainfall2\\RO_Rf_ratios_real_monthly3_2000_2010.txt') ,os.path.join('..', 'input')) # note this is from the simplified version with zone IDs starting at 1

# Rain Fragments file
shutil.copy2(os.path.join(STD_INPUT_FOLDER, "Fragments", 'Rainfall_fragments_2001.prn') ,os.path.join('..', 'input'))  

#  Fragments sequence file
shutil.copy2(os.path.join(STD_INPUT_FOLDER, "Fragments", 'Sequence_file_2002.prn') ,os.path.join('..', 'input')) 

# need to run this before the 1st model run to re-fresh the direct net infiltration coverage to not include the MFR. 
arcpy.Plus_3d(os.path.join('..', 'input/Direct_infiltration',  "A_WL_Rast.asc"), os.path.join('..', 'input/Direct_infiltration', "OSDS_inlf_in.asc"), os.path.join('..', 'input/Direct_infiltration', "temprast"))
arcpy.RasterToASCII_conversion(os.path.join('..', 'input/Direct_infiltration', "temprast"), os.path.join('..', 'input/Direct_infiltration', "Total_inlf_in.asc"))


# RUN Da MODEL (with no MFR) 
os.chdir(os.path.join("..", "run"))
# Executable and control file copies
shutil.copy2(os.path.join("." , 'swb2.exe') ,os.path.join('..', 'output')) 
shutil.copy2(os.path.join("." , Control_File_Name) ,os.path.join('..', 'output')) 

os.chdir(os.path.join("..", "output"))
subprocess.call('swb2.exe {}'.format(Control_File_Name), shell=True)
os.chdir(os.path.join("..", "run"))

### Post process da files
outspace = os.path.join('..', "output", 'post_prcessed_no_MFR')
if not os.path.exists(outspace):
    os.makedirs(outspace)
    
# Parameters
Desired_files = ['actual_et',  'direct_net_infiltation', 'direct_soil_moisture',
             'interception', 'net_infiltration', 'rainfall', 'runoff'] # 'delta_soil_storage',  'irrigation', 
XLLCORNER =      515000.000
YLLCORNER =      8409000.000
CELLSIZE  =      cel_size

# functions
def create_file_reference( component_name ):
    '''
    This is a simple convenience function that will form a path and filename to a
    given water budget component
    '''
    # specify the prefix, path to SWB2 output, timeframe, and resolution
    #output_path = os.path.join(os.getcwd(), "output")
    #prefix      = '\\'
    start_year  = '2000-01-01'
    end_year    = '2009-12-31'
    ncol        = str(int(x_dim))
    nrow        = str(int(y_dim))
    return(  component_name + '__' + start_year + '_' 
          + end_year + '__' + nrow + '_by_' + ncol + '.nc' )

# some other functions to post process stuff

def ncdump(nc_fid):
    '''ncdump outputs dimensions, variables and their attribute information of netCDF4 files'''
    nc_attrs = nc_fid.ncattrs()
    nc_dims = [dim for dim in nc_fid.dimensions]  
    nc_vars = [var for var in nc_fid.variables] 
    return nc_attrs, nc_dims, nc_vars

def writeArrayToArcGrid(arr,filename,xll,yll,cellsize,no_data_val):
    """ this takes a 2d numpy array and turns it into an .asc file """
    arr                = np.copy(arr)
    arr[np.isnan(arr)] = no_data_val
    headerstring       = bytes('NCOLS %d\nNROWS %d\nXLLCENTER %f\nYLLCENTER %f\nCELLSIZE %f\nNODATA_value %f\n' % 
        (arr.shape[1], arr.shape[0], xll, yll, cellsize, no_data_val), 'UTF-8')

    with open(filename,'wb') as fout:
        fout.write(headerstring)
        np.savetxt(fout,arr,'%5.2f')
        

# post process the whole model domain
os.chdir(os.path.join("..", 'output'))  # difficulty in making the path to the file so need to change into the output directory
var = []; tot = []; nclist = []

# Step 1: make list of files that you wish to process 
for i in Desired_files:
    Da_file = create_file_reference(i)
    nclist.append(Da_file)
    
# Step 2 average the daily dimension (len(t) is # of days in the run) to annual 
for i, f in enumerate(nclist):
    nc_data = nc.Dataset(nclist[i])
    nc_attrs, nc_dims, nc_vars = ncdump(nc_data)
    nc_var = nc_vars[3]
    t = nc_data.variables['time'][:]
    y = nc_data.variables['y'][:]
    x = nc_data.variables['x'][:]
    nt = len(t)
    nrow = len(y)
    ncol = len(x)
    rd = np.zeros((nrow, ncol))  # create 0 array of the proper shape
    for day in range(nt):
        r_temp = nc_data.variables[str(nc_var)][day, :, :]
        r_filled = np.ma.filled(r_temp, fill_value=0)    # fills in missing values with 0s  
        rd = rd+r_filled                                 # sequentially add each day's value in each cell to the empty frame  
    r = rd/nt*365 # to create a one year average from all the years in model.  if want to add leap years add 0.25 
    
    # step 3: write each yearly average array to a .asc file
    keyname = Desired_files[i] 
    writeArrayToArcGrid(r, os.path.join(outspace, "{}_annual.asc".format(keyname)), XLLCORNER, YLLCORNER, CELLSIZE, -999)
    
    # Step 4: calculate total amounts of water in cubic meters per day and create statistics dataframe
    m3pd = ((cel_size**2)*r.sum()*.0254)/365 
    print("{} total  {} [m3/d]".format(keyname, '%.1f' % m3pd))
    var.append(keyname) ; tot.append(m3pd)     # make lists to populate pandas dataframe
    
    nc_data.close()          # make sure to close the nc file so it doesnt stay open

stat_frame = pd.DataFrame({'Variable' : var, 'total_[m3pd]': tot})    #in case you want the max and min#, "Max_[in]": mx, "Min_[in]":mn})
stat_frame["total_[MGD]"] = stat_frame["total_[m3pd]"]/3785.41178       # put things in MGD if interested, 3785.41178 is number of gal in m3      
Precip = list(stat_frame[stat_frame['Variable'] == 'rainfall']['total_[m3pd]'])[0]   # define the amount of calculated Precip
Dir_net_inf = list(stat_frame[stat_frame['Variable'] == 'direct_net_infiltation']['total_[m3pd]'])[0]   # define the amount of calculated infiltration
WB_ins = Precip + Dir_net_inf
stat_frame['pct_of_pcip'] = stat_frame["total_[m3pd]"]/WB_ins
stat_frame.to_csv(os.path.join(outspace, "stats_run7_{}m_cells.csv".format(cel_size)))

# how does the model balance? 
print("WATER BALANCE ratio: outs over ins water budget balanece =  {} % ".format(stat_frame['pct_of_pcip'].sum()-1))   # check water balance

os.chdir(os.path.join("..", 'run'))  # then back out to the home directory

# calculate statistics for individual watersheds
# note, for some reason will not overwrite csvs need to clear them out or recode to make this issue not an issue
#create workspace
outspace_table = os.path.join('..', 'output', 'post_prcessed_no_MFR', "tables")
if not os.path.exists(outspace_table):
    os.makedirs(outspace_table)
sheds = (os.path.join(GIS_FOLDER, 'Watersheds\\Runoff_zones_sheds_WGS2S_clip.shp'))

# process each raster layer into a table
for i in (os.listdir(outspace)):
    if i.endswith('.asc'):
        outZSaT = ZonalStatisticsAsTable(sheds, "SHED_NAME", os.path.join(outspace, i), os.path.join(outspace_table, "temptab.dbf"))  # in arc format
        arcpy.TableToTable_conversion(outZSaT, outspace_table, "Table_{}_1.csv".format(i))                                            # take table out of stupid arc format and put into csv format 
        
# this block takes each of the csvs, reads them and calculates water volumnes (m3/d) for each watershed
templist = []
for c in (os.listdir(os.path.join(outspace, "tables"))):
    if c.endswith('.csv'):
        data = pd.read_csv(os.path.join(outspace, "tables", c))
        keyname = c.split("Table_")[1].split("_annual")[0]                   # parameter being worked on
        data[keyname] = (data['MEAN']*.0254/365) * data['AREA'] 
        temp_frame = data[["SHED_NAME", keyname]]
        templist.append(temp_frame)
        
summarry_frame1 = data[['SHED_NAME']]                                        # this is just sticking them all together into one dataframe
for i in templist:
    summarry_frame1 = summarry_frame1.merge(i, on ='SHED_NAME', how='outer')
                          

# that was in actual volumns, now to convert each component into a fraction of the rainfall...
templist2 = []
summarry_frame2 = data[['SHED_NAME']]
for i in summarry_frame1.columns[1:]:
    temp_frame = data[['SHED_NAME']] ; temp_frame[i.split("-")[0]] = summarry_frame1[i]/summarry_frame1['rainfall']
    templist2.append(temp_frame)
    
summarry_frame3 = data[['SHED_NAME']]
for i in templist2:
    summarry_frame3 = summarry_frame3.merge(i, on ='SHED_NAME', how='outer')
                          
summarry_frame_4000 = summarry_frame1.set_index('SHED_NAME')
summarry_frame_4 = summarry_frame_4000.select_dtypes(exclude=['object'])*264.172/1000000   # convert to million gallons per day
    
summarry_frame3.to_csv(os.path.join(outspace, "watershed_summary_stats_percentages.csv"))
summarry_frame1.to_csv(os.path.join(outspace, "watershed_summary_stats_volume_m3pd.csv"))
summarry_frame_4.to_csv(os.path.join(outspace, "watershed_summary_stats_volumes_MGD.csv"))


### MFR calculations      
outspace = os.path.join('..', "output", 'post_prcessed_no_MFR')
if not os.path.exists(outspace):
    os.makedirs(outspace)

# caclulate how much runoff to dump into the MFR area
outspace_table = os.path.join('..', 'output', 'MFR_calcs', "tables")
if not os.path.exists(outspace_table):
    os.makedirs(outspace_table)
    
Contributing_area_leo = (os.path.join(GIS_FOLDER, 'MFR\\Contributing_MRF_Areas_leone.shp'))
Contributing_area_taf = (os.path.join(GIS_FOLDER, 'MFR\\Contributing_MRF_Areas_tafuna.shp'))

outZSaT_leo = ZonalStatisticsAsTable(Contributing_area_leo, "SHED_NAME", os.path.join(outspace, "runoff_annual.asc"), os.path.join(outspace_table, "temptab_leo.dbf"))  # in arc format
arcpy.TableToTable_conversion(outZSaT_leo, outspace_table, "runoff_MFR_leo.csv")                                            # take table out of stupid arc format and put into csv format 
outZSaT_leo = ZonalStatisticsAsTable(Contributing_area_taf, "SHED_NAME", os.path.join(outspace, "runoff_annual.asc"), os.path.join(outspace_table, "temptab_taf.dbf"))  # in arc format
arcpy.TableToTable_conversion(outZSaT_leo, outspace_table, "runoff_MFR_taf.csv") 

data_leo = pd.read_csv(os.path.join(outspace_table, "runoff_MFR_leo.csv"))
data_taf = pd.read_csv(os.path.join(outspace_table, "runoff_MFR_taf.csv"))

data_leo["AreaRunoff_m3pd"] = (data_leo['MEAN']*.0254/365) * data_leo['AREA']    # this is how much runoff is in each MFR contributionzone
data_taf["AreaRunoff_m3pd"] = (data_taf['MEAN']*.0254/365) * data_taf['AREA']    # this is how much runoff is in each MFR contributionzone
tot_MFR_leo = sum(data_leo['AreaRunoff_m3pd'])
tot_MFR_taf = sum(data_taf['AreaRunoff_m3pd'])

# calculate the MFR area and prepare input files
workspace = os.path.join('..', 'input/MFR')
if not os.path.exists(workspace):
    os.makedirs(workspace)

arcpy.Project_management(os.path.join(GIS_FOLDER, 'MFR\\MFR_infiltration_area_leone.shp'),  os.path.join(workspace, 'MFR_infiltration_area_leone_projected.shp'), sr_project) 
arcpy.AddField_management(os.path.join(workspace, 'MFR_infiltration_area_leone_projected.shp'), "MFR_inch", "DOUBLE")    # add Active cell unit field
arcpy.AddGeometryAttributes_management(os.path.join(workspace, 'MFR_infiltration_area_leone_projected.shp'), "AREA")
Total_MFR_area_leo = 0                                                                                                        # stupid block just to calculate the total area
with arcpy.da.SearchCursor(os.path.join(workspace, 'MFR_infiltration_area_leone_projected.shp'), "POLY_AREA") as cursor:
    for row in cursor:
        Total_MFR_area_leo = Total_MFR_area_leo + row[0]

arcpy.Project_management(os.path.join(GIS_FOLDER, 'MFR\\MFR_infiltration_area_tafuna.shp'),  os.path.join(workspace, 'MFR_infiltration_area_tafuna_projected.shp'), sr_project) 
arcpy.AddField_management(os.path.join(workspace, 'MFR_infiltration_area_tafuna_projected.shp'), "MFR_inch", "DOUBLE")    # add Active cell unit field
arcpy.AddGeometryAttributes_management(os.path.join(workspace, 'MFR_infiltration_area_tafuna_projected.shp'), "AREA")
Total_MFR_area_taf = 0                                                                                                        # stupid block just to calculate the total area
with arcpy.da.SearchCursor(os.path.join(workspace, 'MFR_infiltration_area_tafuna_projected.shp'), "POLY_AREA") as cursor:
    for row in cursor:
        Total_MFR_area_taf = Total_MFR_area_taf + row[0]
        
Inches_of_MFR_across_leo = (tot_MFR_leo/Total_MFR_area_leo/0.0254) * 0.75   # note this 75% number if directly from Izuka 2007
Inches_of_MFR_across_taf = (tot_MFR_taf/Total_MFR_area_taf/0.0254) * 0.75   # note this 75% number if directly from Izuka 2007

arcpy.CalculateField_management(os.path.join(workspace, 'MFR_infiltration_area_leone_projected.shp'), "MFR_inch", "!MFR_inch! + {}".format(Inches_of_MFR_across_leo), "PYTHON3") # calculate the appropriate amount of infitration in inches spread over all MFR zone
arcpy.CalculateField_management(os.path.join(workspace, 'MFR_infiltration_area_tafuna_projected.shp'), "MFR_inch", "!MFR_inch! + {}".format(Inches_of_MFR_across_taf), "PYTHON3") # calculate the appropriate amount of infitration in inches spread over all MFR zone

arcpy.Erase_analysis(Grid_shp, os.path.join(workspace, 'MFR_infiltration_area_leone_projected.shp'),  os.path.join(workspace, 'MFR_infiltration_area_leone_bound.shp'))
arcpy.Erase_analysis(Grid_shp, os.path.join(workspace, 'MFR_infiltration_area_tafuna_projected.shp'),  os.path.join(workspace, 'MFR_infiltration_area_tafuna_bound.shp'))

arcpy.Merge_management([os.path.join(workspace, 'MFR_infiltration_area_leone_bound.shp'), os.path.join(workspace, 'MFR_infiltration_area_leone_projected.shp')], os.path.join(workspace, 'MFR_infiltration_area_leone_ready.shp'))
arcpy.Merge_management([os.path.join(workspace, 'MFR_infiltration_area_tafuna_bound.shp'), os.path.join(workspace, 'MFR_infiltration_area_tafuna_projected.shp')], os.path.join(workspace, 'MFR_infiltration_area_tafuna_ready.shp'))

arcpy.PolygonToRaster_conversion(os.path.join(workspace, 'MFR_infiltration_area_leone_ready.shp'), "MFR_inch", os.path.join(workspace, "MFR_Rast_L"), cell_assignment="MAXIMUM_AREA",  cellsize=cel_size)
arcpy.PolygonToRaster_conversion(os.path.join(workspace, 'MFR_infiltration_area_tafuna_ready.shp'), "MFR_inch", os.path.join(workspace, "MFR_Rast_T"), cell_assignment="MAXIMUM_AREA",  cellsize=cel_size)

arcpy.RasterToASCII_conversion(os.path.join(workspace, "MFR_Rast_L"), os.path.join(workspace, "MFR_Rast_L.asc"))
arcpy.RasterToASCII_conversion(os.path.join(workspace, "MFR_Rast_T"), os.path.join(workspace, "MFR_Rast_T.asc"))

arcpy.Delete_management(os.path.join(workspace, 'MFR_infiltration_area_leone_projected.shp'))
arcpy.Delete_management(os.path.join(workspace, 'MFR_infiltration_area_tafuna_projected.shp'))
arcpy.Delete_management(os.path.join(workspace, 'MFR_infiltration_area_leone_bound.shp'))
arcpy.Delete_management(os.path.join(workspace, 'MFR_infiltration_area_tafuna_bound.shp'))
arcpy.Delete_management(os.path.join(workspace, 'MFR_infiltration_area_leone_ready.shp'))
arcpy.Delete_management(os.path.join(workspace, 'MFR_infiltration_area_tafuna_ready.shp'))
arcpy.Delete_management(os.path.join(workspace, 'MFR_Rast_L'))
arcpy.Delete_management(os.path.join(workspace, 'MFR_Rast_T'))

# now combine the MFR raster into the other direct infiltration rasters
arcpy.Plus_3d(os.path.join('..', 'input/MFR', "MFR_Rast_L.asc"), os.path.join('..', 'input/Direct_infiltration', "WLOSDrast"), os.path.join('..', 'input/Direct_infiltration', "temprast2"))
arcpy.Plus_3d(os.path.join('..', 'input/MFR', "MFR_Rast_T.asc"), os.path.join('..', 'input/Direct_infiltration', "temprast2"), os.path.join('..', 'input/Direct_infiltration', "temprast3"))
arcpy.RasterToASCII_conversion(os.path.join('..', 'input/Direct_infiltration', "temprast3"), os.path.join('..', 'input/Direct_infiltration', "Total_inlf_in.asc"))

print('MFR leo in MGD is {}'.format(tot_MFR_leo*264.172/1000000))
print('MFR taf in MGD is {}'.format(tot_MFR_taf*264.172/1000000))
print('MFR total in MGD is {}'.format((tot_MFR_leo+tot_MFR_taf)*264.172/1000000))

# Run da Model again, this time including the MFR
# Executable and control file copies
shutil.copy2(os.path.join("." , 'swb2.exe') ,os.path.join('..', 'output')) 
shutil.copy2(os.path.join("." , Control_File_Name) ,os.path.join('..', 'output')) 

os.chdir(os.path.join("..", "output"))
subprocess.call('swb2.exe {}'.format(Control_File_Name), shell=True)
os.chdir(os.path.join("..", "run"))

# Post process the files again, this time with the MFR added 
outspace = os.path.join('..', "output", 'post_prcessed_with_MFR')
if not os.path.exists(outspace):
    os.makedirs(outspace)

# post process the whole model domain
os.chdir(os.path.join("..", 'output'))  # difficulty in making the path to the file so need to change into the output directory
var = []; tot = []; nclist = []

# Step 1: make list of files that you wish to process 
for i in Desired_files:
    Da_file = create_file_reference(i)
    nclist.append(Da_file)
    
# Step 2 average the daily dimension (len(t) is # of days in the run) to annual 
for i, f in enumerate(nclist):
    nc_data = nc.Dataset(nclist[i])
    nc_attrs, nc_dims, nc_vars = ncdump(nc_data)
    nc_var = nc_vars[3]
    t = nc_data.variables['time'][:]
    y = nc_data.variables['y'][:]
    x = nc_data.variables['x'][:]
    nt = len(t)
    nrow = len(y)
    ncol = len(x)
    rd = np.zeros((nrow, ncol))  # create 0 array of the proper shape
    for day in range(nt):
        r_temp = nc_data.variables[str(nc_var)][day, :, :]
        r_filled = np.ma.filled(r_temp, fill_value=0)    # fills in missing values with 0s (i think) 
        rd = rd+r_filled                                 # sequentially add each day's value in each cell to the empty frame  
    r = rd/nt*365 # to create a one year average from all the years in model.  if want to add leap years add 0.25 
    
    # step 3: write each yearly average array to a .asc file
    keyname = Desired_files[i] 
    writeArrayToArcGrid(r, os.path.join(outspace, "{}_annual.asc".format(keyname)), XLLCORNER, YLLCORNER, CELLSIZE, -999)
    
    # Step 4: calculate total amounts of water in cubic meters per day and create statistics dataframe
    m3pd = ((cel_size**2)*r.sum()*.0254)/365 
    print("{} total  {} [m3/d]".format(keyname, '%.1f' % m3pd))
    var.append(keyname) ; tot.append(m3pd)     # make lists to populate pandas dataframe
    
    nc_data.close()          # make sure to close the nc file so it doesnt stay open

stat_frame = pd.DataFrame({'Variable' : var, 'total_[m3pd]': tot})    #in case you want the max and min#, "Max_[in]": mx, "Min_[in]":mn})
stat_frame["total_[MGD]"] = stat_frame["total_[m3pd]"]/3785.41178       # put things in MGD if interested, 3785.41178 is number of gal in m3      
Precip = list(stat_frame[stat_frame['Variable'] == 'rainfall']['total_[m3pd]'])[0]   # define the amount of calculated Precip
Dir_net_inf = list(stat_frame[stat_frame['Variable'] == 'direct_net_infiltation']['total_[m3pd]'])[0]   # define the amount of calculated infiltration
WB_ins = Precip + Dir_net_inf
stat_frame['pct_of_pcip'] = stat_frame["total_[m3pd]"]/WB_ins
stat_frame.to_csv(os.path.join(outspace, "stats_run7_{}m_cells.csv".format(cel_size)))

# how does the model balance? 
print("WATER BALANCE ratio: outs over ins water budget balanece =  {} % ".format(stat_frame['pct_of_pcip'].sum()-1))   # check water balance

os.chdir(os.path.join("..", 'run'))  # then back out to the home directory

# calculate statistics for individual watersheds
# note, for some reason will not overwrite csvs need to clear them out or recode to make this issue not an issue
#create workspace
outspace_table = os.path.join('..', 'output', 'post_prcessed_with_MFR', "tables")
if not os.path.exists(outspace_table):
    os.makedirs(outspace_table)
sheds = (os.path.join(GIS_FOLDER, 'Watersheds\\Runoff_zones_sheds_WGS2S_clip.shp'))

# process each raster layer into a table
for i in (os.listdir(outspace)):
    if i.endswith('.asc'):
        outZSaT = ZonalStatisticsAsTable(sheds, "SHED_NAME", os.path.join(outspace, i), os.path.join(outspace_table, "temptab.dbf"))  # in arc format
        arcpy.TableToTable_conversion(outZSaT, outspace_table, "Table_{}_1.csv".format(i))                                            # take table out of stupid arc format and put into csv format 
        
# this block takes each of the csvs, reads them and calculates water volumnes (m3/d) for each watershed
templist = []
for c in (os.listdir(os.path.join(outspace, "tables"))):
    if c.endswith('.csv'):
        data = pd.read_csv(os.path.join(outspace, "tables", c))
        keyname = c.split("Table_")[1].split("_annual")[0]                   # parameter being worked on
        data[keyname] = (data['MEAN']*.0254/365) * data['AREA'] 
        temp_frame = data[["SHED_NAME", keyname]]
        templist.append(temp_frame)
        
summarry_frame1 = data[['SHED_NAME']]                                        # this is just sticking them all together into one dataframe
for i in templist:
    summarry_frame1 = summarry_frame1.merge(i, on ='SHED_NAME', how='outer')
                          

# that was in actual volumns, now to convert each component into a fraction of the rainfall...
templist2 = []
summarry_frame2 = data[['SHED_NAME']]
for i in summarry_frame1.columns[1:]:
    temp_frame = data[['SHED_NAME']] ; temp_frame[i.split("-")[0]] = summarry_frame1[i]/summarry_frame1['rainfall']
    templist2.append(temp_frame)
    
summarry_frame3 = data[['SHED_NAME']]
for i in templist2:
    summarry_frame3 = summarry_frame3.merge(i, on ='SHED_NAME', how='outer')
                          
summarry_frame_4000 = summarry_frame1.set_index('SHED_NAME')
summarry_frame_4 = summarry_frame_4000.select_dtypes(exclude=['object'])*264.172/1000000   # convert to million gallons per day
    
summarry_frame3.to_csv(os.path.join(outspace, "watershed_summary_stats_percentages.csv"))
summarry_frame1.to_csv(os.path.join(outspace, "watershed_summary_stats_volume_m3pd.csv"))
summarry_frame_4.to_csv(os.path.join(outspace, "watershed_summary_stats_volumes_MGD.csv"))

## Plot results
for a in (os.listdir(outspace)):
    if a.endswith('.asc'):
        rfs_data, rfs_gt, rfs_proj, rfs_xy    = read_raster(os.path.join(outspace, a))
        make_plot( x=rfs_xy[0], y=rfs_xy[1], var=rfs_data, discrete=True,  title=a, barlabel='[in.]')
        plt.savefig(os.path.join(fig_path, "Map_of_{},_at{}m.png".format(a, cel_size)), dpi = 300)
       

       
       
       
       
       
       