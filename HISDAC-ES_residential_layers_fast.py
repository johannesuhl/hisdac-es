# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 16:47:27 2023

@author: Johannes H. Uhl, University of Colorado Boulder, USA.
"""

import os,sys
import subprocess
import geopandas as gp
import numpy as np
from osgeo import gdal
import scipy.stats
import time

#################################################################################################################

ES_building_centroids_merged = 'ES_building_centroids_merged.shp' #path to shp with country-level building centroids. to be created with HISDAC-ES_muni_stats.py
surface_folder=r'H:\SPAIN_DATA\ES_buildings_raster_laea' #folder to store output tifs.
gdal_edit = r'gdal_edit.py' #path to gdal_edit
gdalwarp = r'C:\OSGeo4W\bin\gdalwarp.exe'#path to gdalwarp
stepsize=10 #tine interval between epochs
years=np.arange(1900,2021,stepsize) #time range
resample_factor=100 #set to resolution of template_raster
template_raster = 'template_epsg3035_100m.tif'
crs_grid = 3035 #epsg of template_raster
    
#################################################################################################################

def gdalNumpy2floatRaster_compressed(array,outname,template_georef_raster,x_pixels,y_pixels,px_type):
    dst_filename = outname
    driver = gdal.GetDriverByName('GTiff')
    dataset = driver.Create(dst_filename,x_pixels, y_pixels, 1, px_type)   
    dataset.GetRasterBand(1).WriteArray(array)                
    mapraster = gdal.Open(template_georef_raster, gdal.GA_ReadOnly)
    proj=mapraster.GetProjection() #you can get from a existing tif or import 
    dataset.SetProjection(proj)
    dataset.FlushCache()
    dataset=None                
    #set bounding coords
    ulx, xres, xskew, uly, yskew, yres  = mapraster.GetGeoTransform()
    lrx = ulx + (mapraster.RasterXSize * xres)
    lry = uly + (mapraster.RasterYSize * yres)            
    mapraster = None                    
    gdal_cmd = gdal_edit+' -a_ullr %s %s %s %s "%s"' % (ulx,uly,lrx,lry,outname)
    print(gdal_cmd)
    response=subprocess.check_output(gdal_cmd, shell=True)
    print(response)    
    outname_lzw=outname.replace('.tif','_lzw.tif')
    gdal_translate = r'gdal_translate %s %s -co COMPRESS=LZW' %(outname,outname_lzw)
    print(gdal_translate)
    response=subprocess.check_output(gdal_translate, shell=True)
    print(response)
    os.remove(outname)
    os.rename(outname_lzw,outname)
    
#################################################################################################################      

xcoo_col,ycoo_col = 'x','y'  

raster = gdal.Open(template_raster)
cols = raster.RasterXSize
rows = raster.RasterYSize
geotransform = raster.GetGeoTransform()
topleftX = geotransform[0]
topleftY = geotransform[3]
pixelWidth = int(abs(geotransform[1]))
pixelHeight = int(abs(geotransform[5]))
rasterrange=[[topleftX,topleftX+pixelWidth*cols],[topleftY-pixelHeight*rows,topleftY]]    
del raster
       
statistic=np.nansum
bitdepth=gdal.GDT_Int32
yearlist=list(years)
      
starttime=time.time()

indf = gp.read_file(ES_building_centroids_merged)

indf = indf[indf.lu_harm=='residential']
if not indf.crs.to_epsg()==crs_grid:
    indf.geometry = indf.geometry.to_crs(epsg=crs_grid)      
indf[xcoo_col]=indf.geometry.x
indf[ycoo_col]=indf.geometry.y

### delete irrelevant columns to reduce file size
del indf['num_floors']
del indf['num_dwel']
del indf['num_bunits']

### rasterize the increments (ie newly built-up per time period defined by stepsize)

for year in years:
    if year==1900:            
        lower_year=1
    else:
        lower_year=year-stepsize           

    yrdf = indf[np.logical_and(indf['yearbuilt']<=year,indf.yearbuilt>lower_year)]

    if len(yrdf)==0:
        continue
    statistic_str='resonlysum_%s_%s' %(year,lower_year)
    
    out_surface_bia =np.zeros((cols,rows),dtype=np.uint16)
    out_surface_bufa =np.zeros((cols,rows),dtype=np.uint16)

    target_variable='area'
    yrdf[target_variable]=yrdf[target_variable].map(float).map(np.int32)
    yrdf=yrdf[yrdf[target_variable]>0]
    yrdf = yrdf.dropna(subset=[target_variable])
    #yrdf=yrdf[[xcoo_col,ycoo_col,target_variable]]                   
    statsvals = yrdf[target_variable].values.astype(np.int32)  
    curr_surface = scipy.stats.binned_statistic_2d(yrdf[xcoo_col].values,yrdf[ycoo_col].values,statsvals,statistic,bins=[cols,rows],range=rasterrange).statistic     
    curr_surface = np.nan_to_num(curr_surface).astype(np.int32)         
    out_surface_bufa = np.add(out_surface_bufa,curr_surface)  
    del curr_surface

    target_variable='offi_area'
    yrdf[target_variable]=yrdf[target_variable].str.replace(',','.')
    yrdf[target_variable]=yrdf[target_variable].map(float).map(np.int32)
    yrdf=yrdf[yrdf[target_variable]>0]
    yrdf = yrdf.dropna(subset=[target_variable])
    #yrdf=yrdf[[xcoo_col,ycoo_col,target_variable]]                   
    statsvals = yrdf[target_variable].values.astype(np.int32)  
    curr_surface = scipy.stats.binned_statistic_2d(yrdf[xcoo_col].values,yrdf[ycoo_col].values,statsvals,statistic,bins=[cols,rows],range=rasterrange).statistic     
    curr_surface = np.nan_to_num(curr_surface).astype(np.int32)         
    out_surface_bia = np.add(out_surface_bia,curr_surface)  
    del curr_surface
    
    print(year,'resbia,resbufa')

    gdalNumpy2floatRaster_compressed(np.rot90(out_surface_bia),surface_folder+os.sep+'ES_buildings_%s_%s_%s_increment.tif' %('resbia',statistic_str,resample_factor),template_raster,cols,rows,bitdepth)
    gdalNumpy2floatRaster_compressed(np.rot90(out_surface_bufa),surface_folder+os.sep+'ES_buildings_%s_%s_%s_increment.tif' %('resbufa',statistic_str,resample_factor),template_raster,cols,rows,bitdepth)
    del out_surface_bia,out_surface_bufa
    
##### now merge the increments to cumulative counts:

for year in years:        
    
    if year ==1900:
        lower_year=1
    else:
        lower_year=year-stepsize
    
    statistic_str='resonlysum_%s_%s' %(year,lower_year)                
    curr_bia_incr = surface_folder+os.sep+'ES_buildings_%s_%s_%s_increment.tif' %('resbia',statistic_str,resample_factor)
    statistic_str='resonlysum_%s' %(year)                        
    curr_bia_cum = surface_folder+os.sep+'ES_buildings_%s_%s_%s_cum.tif' %('resbia',statistic_str,resample_factor)
    
    statistic_str='resonlysum_%s_%s' %(year,lower_year)        
    curr_area_incr = surface_folder+os.sep+'ES_buildings_%s_%s_%s_increment.tif' %('resbufa',statistic_str,resample_factor)
    statistic_str='resonlysum_%s' %(year)                               
    curr_area_cum = surface_folder+os.sep+'ES_buildings_%s_%s_%s_cum.tif' %('resbufa',statistic_str,resample_factor)

    curr_bia_incr_arr=gdal.Open(curr_bia_incr).ReadAsArray()
    curr_area_incr_arr=gdal.Open(curr_area_incr).ReadAsArray()
    
    if year==1900:
        total_bia_surface=curr_bia_incr_arr.copy()
        total_area_surface=curr_area_incr_arr.copy()
    else:
        total_bia_surface=total_bia_surface+curr_bia_incr_arr
        total_area_surface=total_area_surface+curr_area_incr_arr
        
    gdalNumpy2floatRaster_compressed(total_bia_surface,curr_bia_cum,template_raster,cols,rows,bitdepth)
    gdalNumpy2floatRaster_compressed(total_area_surface,curr_area_cum,template_raster,cols,rows,bitdepth)

    print('cumulative %s' %year)
