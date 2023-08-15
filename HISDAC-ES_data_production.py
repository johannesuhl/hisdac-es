# -*- coding: utf-8 -*-
"""
Created on Sun Jun 13 15:05:27 2021

@author: Johannes H. Uhl, University of Colorado Boulder, USA.
"""

import os,sys
import pandas as pd
import requests
from xml.etree import ElementTree as ET   
import subprocess
from zipfile import ZipFile
import geopandas as gp
import numpy as np
import gdal
import scipy.stats
import time

###### input data ###########################
incsv='inspire_buildings_xml_links.csv'
gml_navarra='BuildingsNavarra.gml' #needs to be aquired from https://inspire.navarra.es/services/BU/wfs
gml_alava='alava.gml' #needs to be aquired from https://geo.araba.eus/WFS_Katastroa?SERVICE=WFS&VERSION=1.1.0&REQUEST=GetCapabilities
gml_gipuzkoa='Gipuzkoa.ES.GFA.BU.gml' #needs to be aquired from https://b5m.gipuzkoa.eus/web5000/es/utilidades/inspire/edificios/
path_bizkaia='./gml_bizkaia' #needs to contain gml files from https://web.bizkaia.eus/es/inspirebizkaia

## folders, to be created by user:
metadata_xml_folder='./xml'
data_gml_folder='./gml'
data_unzip_folder='./unzip'
shp_dir='./ES_buildings_shp'
shp_pt_dir='./ES_buildings_shp_pt'
shp_pt_dir_lu_harm='./ES_buildings_shp_pt_lu_harm'
rasterdir_utm='./ES_buildings_raster_utm'
rasterdir_laea='./ES_buildings_raster_laea'
rasterdir_can_regcan='./ES_buildings_raster_can_regcan'
bakdir='./backups'

## template raster layers. Rasters will be aligned to these grids.
template_raster_regcan = 'template_epsg4083_100m.tif'
template_raster_laea = 'template_epsg3035_100m.tif'
template_raster_utm = 'template_epsg25830_all_100m.tif'
lu_mapping_csv='landuse_mapping.csv'

## control variables for individual code blocks:
download=False ## downloads building footprint data per municipality in gml format for regions except Navarra and Basque Country.
convert_geopandas=False ## converts from GML to SHP for each municipality-level dataset.
harmonize_data=False ## takes the previously created shapefiles and the manually downloaded data for Basque country and Navarra, and harmonizes them.
mine_landuse=False ## analyzes the land use / building function types used in each data model
harmonize_landuse=False ##harmonizes the land use / building function classes used across data models.
rasterize_age=False ## create age-related gridded surfaces
rasterize_mutemp=False ## create evolutionary layers
rasterize_landuse_mutemp=False ##creates multitemporal land use layers
rasterize_physical_characteristics=False ##creates physical properties layers

## spatial and temporal dimensions:
resample_factor=100 #spatial resolution in meters
years=np.arange(1900,2026,5) #temporal resolution / coverage

## GDAL binaries:
gdal_edit = r'python C:\OSGeo4W\bin\gdal_edit.py' #path to gdal_edit.py
gdalwarp = r'C:\OSGeo4W\bin\gdalwarp.exe' #path to gdalwarp.exe

## some functions ##############################################################################################

def variety(x):
    return np.unique(x).shape[0]

def mode(x):
    vals,counts = np.unique(x, return_counts=True)
    index = np.argmax(counts)
    return vals[index]   
    
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
    
def getyear1(x):
    try:
        if '-' in x[:4]:
            return 0            
        else:
            return int(x[:4])
    except:
        return 0

xcoo_col,ycoo_col = 'x','y' ## don't change
################################################################################################
    
if download:
    xmlurldf=pd.read_csv(incsv,header=None)
    xmlurldf.columns=['url']
    for xmlurl in xmlurldf.url.values:
        province=xmlurl.split('_')[-1].replace('.xml','')    
        remotexmlfile = requests.get(xmlurl) 
        localxmlfile = metadata_xml_folder+os.sep+xmlurl.split('/')[-1]
        open(localxmlfile, 'wb').write(remotexmlfile.content)  
        try:                       
            root = ET.parse(localxmlfile).getroot()  
        except:
            continue
        for child in root:
            for el in child:
                if el.tag.split('}')[-1]=='id':
                    gmlzipurl = el.text                
                    ###download gml zip
                    remotefile = requests.get(gmlzipurl) 
                    localfile = data_gml_folder+os.sep+gmlzipurl.split('/')[-1]
                    open(localfile, 'wb').write(remotefile.content)
                    print(province,gmlzipurl)

################################################################################################

if convert_geopandas:

    count=0
    for zipfile in os.listdir(data_gml_folder):
        count+=1
        try:
            muni=zipfile.split('.')[-2] 
            province=muni[:2]
            with ZipFile(data_gml_folder+os.sep+zipfile, 'r') as zipObj:
               zipObj.extractall(data_unzip_folder)            
            ingml=data_unzip_folder+os.sep+zipfile.replace('.zip','.building.gml')
            #### get_epsg:
            inxml=data_unzip_folder+os.sep+zipfile.replace('.zip','.xml').replace('BU','BU.MD')
            root = ET.parse(inxml).getroot()  
            for child in root:
                for el in child:
                    if 'MD_ReferenceSystem' in el.tag:
                        for child2 in el:
                            if 'referenceSystemIdentifier' in child2.tag:  
                                for child3 in child2:
                                    if 'RS_Identifier' in child3.tag:  
                                        for child4 in child3:
                                            if 'code' in child4.tag:  
                                                for child5 in child4:
                                                    if 'CharacterString' in child5.tag:  
                                                        epsg=child5.text.split('/')[-1]
       
            ingdf=gp.read_file(ingml)
            ingdf.set_crs(epsg=int(epsg),inplace=True)
            ingdf['begin_int']=ingdf.apply(lambda row : getyear1(row['beginning']),axis=1).map(int)
            ingdf['end_int']=ingdf.apply(lambda row : getyear1(row['end']),axis=1).map(int)
            
            ingdf.to_file(filename=shp_dir+os.sep+zipfile.replace('.zip','.shp'))
            print('shp written',count,province,muni)
            for file in os.listdir(data_unzip_folder):
                os.remove(data_unzip_folder+os.sep+file)
        except:
            print('ERROR',count,province,muni)
            print('error %s %s' %(province,muni), file=open('ES_buildings_errors.txt','a'))
            #sys.exit(0)

################################################################################################

if harmonize_data:
    
    yearcol_lookup={'ES_INSPIRE':'beginning',
                    'BIZKAIA':'end',
                    'NAVARRA':'BU_dateconsbegin',
                    'ALAVA':'dateOfConstruction',
                    'GIPUZKOA':'dateOfConstruction|DateOfEvent|end'}
    
    usage_lookup={'ES_INSPIRE':'currentUse',
                    'BIZKAIA':'currentUse',
                    'NAVARRA':'BU_currentUseValue',
                    'ALAVA':'currentUse',
                    'GIPUZKOA':'currentUse'}
    
    idcol_lookup={'ES_INSPIRE':'gml_id',
                    'BIZKAIA':'gml_id',
                    'NAVARRA':'gml_id',
                    'ALAVA':'inspireId',
                    'GIPUZKOA':'gml_id'} 
    
    numfloors_lookup={'ES_INSPIRE':'numberOfFl',
                    'BIZKAIA':'numberOfFloorsAboveGround',
                    'NAVARRA':'BU_numberOfFloorsAboveGround',
                    'ALAVA':'numberOfFloorsAboveGround',
                    'GIPUZKOA':'numberOfFloorsAboveGround'}

    numdwel_lookup={'ES_INSPIRE':'numberOfDw',
                    'BIZKAIA':'numberOfDwellings',
                    'NAVARRA':'BU_numberOfDwellings',
                    'ALAVA':'',
                    'GIPUZKOA':'numberOfDwellings'}

    numbunits_lookup={'ES_INSPIRE':'numberOfBu',
                    'BIZKAIA':'numberOfBuildingUnits',
                    'NAVARRA':'BU_numberOfBuildingUnits',
                    'ALAVA':'numberOfBuildingUnits',
                    'GIPUZKOA':''}

    area_lookup=   {'ES_INSPIRE':'value',
                    'BIZKAIA':'',
                    'NAVARRA':'BU_officialArea',
                    'ALAVA':'',
                    'GIPUZKOA':''}

    areatype_lookup={'ES_INSPIRE':'officialAr',
                    'BIZKAIA':'',
                    'NAVARRA':'BU_officialAreaReference',
                    'ALAVA':'',
                    'GIPUZKOA':''}
    
    total_files=[]
    for file in os.listdir(shp_dir):
        if file.split('.')[-1]=='shp':
            total_files.append(['ES_INSPIRE',shp_dir+os.sep+file])
    
    for file in os.listdir(path_bizkaia):        
        total_files.append(['BIZKAIA',path_bizkaia+os.sep+file+os.sep+'ES.BFA.BU.gml']) 

    total_files.append(['NAVARRA',gml_navarra])
    total_files.append(['ALAVA',gml_alava])
    total_files.append(['GIPUZKOA',gml_gipuzkoa])

    total_files_df=pd.DataFrame(total_files,columns=['region','fullpath'])   
    year_stats=[]
    for i,row in total_files_df.iterrows():

        datadf=gp.read_file(row.fullpath)
        yearcol=yearcol_lookup[row.region]
        usagecol=usage_lookup[row.region]
        idcol=idcol_lookup[row.region]        
        numfloors_col=numfloors_lookup[row.region]
        numdwel_col=numdwel_lookup[row.region]
        numbunits_col=numbunits_lookup[row.region]       
        area_col=area_lookup[row.region]
        areatype_col=areatype_lookup[row.region]
        
        datadf['yearbuilt']=datadf.apply(lambda row : getyear1(row[yearcol]),axis=1).map(int)            
        datadf['yearbuilt']= datadf['yearbuilt'].replace(0,np.nan)
        datadf['area']=datadf.geometry.to_crs(epsg=3035).area # calculate footprint area in LAEA projection

        datadf['idx']=datadf.index
        datadf['source']=row.fullpath
        datadf['region']=row.region  
        try:
            datadf['usage']=datadf[usagecol]
        except:
            datadf['usage']=np.nan
            
        datadf['identifier']=datadf[idcol]
        datadf['num_floors']=datadf[numfloors_col]
        
        if numdwel_col=='':
            datadf['num_dwel']=0
        else:
            datadf['num_dwel']=datadf[numdwel_col]

        if numbunits_col=='':
            datadf['num_bunits']=0
        else:
            datadf['num_bunits']=datadf[numbunits_col]            

        if area_col=='':
            datadf['offi_area']=0
        else:
            datadf['offi_area']=datadf[area_col]          

        if areatype_col=='':
            datadf['areatype']=''
        else:
            datadf['areatype']=datadf[areatype_col] 
                
        datadf=datadf[['geometry','yearbuilt','area','idx','source','region','usage','identifier',
                       'num_floors','num_dwel','num_bunits','offi_area','areatype']]
        
        datadf.geometry=datadf.geometry.centroid

        if row.region=='GIPUZKOA':
            datadf['x']=datadf.geometry.centroid.x
            datadf['y']=datadf.geometry.centroid.y
            datadf = gp.GeoDataFrame(datadf, geometry=gp.points_from_xy(datadf['y'], datadf['x']))        
        if row.region=='ALAVA' or row.region=='NAVARRA' or row.region=='GIPUZKOA':
            datadf.set_crs(epsg=25830,inplace=True) 
        if row.region=='BIZKAIA':
            datadf['x']=datadf.geometry.centroid.x
            datadf['y']=datadf.geometry.centroid.y
            datadf = gp.GeoDataFrame(datadf, geometry=gp.points_from_xy(datadf['y'], datadf['x']))                        
            datadf=datadf.set_crs("+init=EPSG:4326",inplace=False)
            datadf=datadf.to_crs(epsg=25830,inplace=False)

        datadf['x']=datadf.geometry.centroid.x
        datadf['y']=datadf.geometry.centroid.y
                    
        year_compl=100*((len(datadf.dropna(subset=['yearbuilt'])))/float(len(datadf)))
        minyr,maxyr=min(datadf.yearbuilt.dropna().values),max(datadf.yearbuilt.dropna().values)
        print(i,row.region,row.fullpath,year_compl,minyr,maxyr)
        year_stats.append([row.region,row.fullpath,year_compl,minyr,maxyr])
        
        outfile=shp_pt_dir+os.sep+'es_bldg_harm_pt_%s_%s.shp' %(row.region,i)
        datadf.to_file(outfile)

################################################################################################
            
if rasterize_age:

    params=[]
    params.append([True,False,False])
    params.append([False,True,False])
    params.append([False,False,True])
    
    for param in params:
        do_can_regcan=param[0]
        do_all_excpt_can_utm=param[1]        
        do_all_laea=param[2]
            
        ################################################# canaries regcan
        if do_can_regcan:
            template_raster = template_raster_regcan
            crs_grid = 4083 #epsg of template_raster
            surface_folder = rasterdir_can_regcan
            do_canaries_only=True
        ################################################# iberic pen and canaries laea,
        if do_all_laea:
            template_raster = template_raster_laea 
            crs_grid = 3035 #epsg of template_raster
            surface_folder = rasterdir_laea
            do_canaries_only=False
        ################################################# iberic pen and bal, mel, ceu, UTM
        if do_all_excpt_can_utm:
            template_raster = template_raster_utm 
            crs_grid = 25830 #epsg of template_raster
            surface_folder = rasterdir_utm
            do_canaries_only=False
        
        if do_canaries_only:    
            can_files=[]
            for file in os.listdir(shp_dir):
                if file.split('.')[-1]=='shp':
                    if '.BU.35' in file or '.BU.38' in file:
                        can_files.append(['CAN',shp_dir+os.sep+file])
                    else:
                        can_files.append(['IB',shp_dir+os.sep+file])                    
            can_files_df=pd.DataFrame(can_files,columns=['region','fullpath'])
            can_idx=can_files_df[can_files_df.region=='CAN'].index.values
        
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
        
        stats=[]
        stats.append(['yearbuilt',np.nanmin,'min',gdal.GDT_Int16]) #gdal.GDT_Float32
        stats.append(['yearbuilt',variety,'variety',gdal.GDT_Int16]) #gdal.GDT_Float32
        stats.append(['yearbuilt',mode,'mode',gdal.GDT_Int16]) #gdal.GDT_Float32
        stats.append(['yearbuilt',np.nanmax,'max',gdal.GDT_Int16]) #gdal.GDT_Float32
        stats.append(['yearbuilt',np.nanmean,'mean',gdal.GDT_Float32]) #gdal.GDT_Float32
        stats.append(['yearbuilt',np.nanmedian,'median',gdal.GDT_Float32]) #gdal.GDT_Float32    
        stats.append(['count',np.nansum,'sum',gdal.GDT_Int16]) #gdal.GDT_Float32
        stats.append(['count_yearbuilt_missing',np.nansum,'sum',gdal.GDT_Int16]) #gdal.GDT_Float32
        stats.append(['area',np.nansum,'sum',gdal.GDT_Float32]) #gdal.GDT_Float32
                  
        for stat in stats:
            target_variable=stat[0]
            statistic=stat[1]
            statistic_str=stat[2]
            bitdepth=stat[3]
            #sys.exit(0)
            
            out_surface =np.zeros((cols,rows)).astype(np.float32)
                    
            if statistic_str=='min':
                out_surface = out_surface+2022
    
            counter=0
            for file in os.listdir(shp_pt_dir):
                
                if not file.split('.')[-1]=='shp':
                    continue
                
                if do_canaries_only:
                    curridx=file.replace('.shp','').split('_')[-1]
                    if not int(curridx) in can_idx:
                        continue
    
                starttime=time.time()            
                
                #if 'BIZKAIA' in file:
                #    continue
    
                indf = gp.read_file(shp_pt_dir+os.sep+file)
    
                if len(indf)==0:
                    continue
        
                if target_variable == 'count_yearbuilt_missing':
                    indf=indf[np.logical_not(indf.yearbuilt>0)]
    
                if target_variable == 'count_area_missing':
                    indf=indf[np.logical_not(indf['area']>0)]
                    
                if len(indf)==0:
                    continue
                
                indf.geometry = indf.geometry.to_crs(epsg=crs_grid)      
                indf[xcoo_col]=indf.geometry.x
                indf[ycoo_col]=indf.geometry.y
                                
                if target_variable == 'count' or target_variable == 'count_yearbuilt_missing' or target_variable == 'count_area_missing':
                    indf[target_variable]=1
                        
                #################################################################
             
                starttime=time.time()
                counter+=1
    
                indf = indf.dropna(subset=[target_variable])
                indf=indf[[xcoo_col,ycoo_col,target_variable]]                   
                statsvals = indf[target_variable].values.astype(np.int32) 
                
                curr_surface = scipy.stats.binned_statistic_2d(indf[xcoo_col].values,indf[ycoo_col].values,statsvals,statistic,bins=[cols,rows],range=rasterrange).statistic     
                
                if statistic_str=='min':   
                    curr_surface = np.nan_to_num(curr_surface)
                    curr_surface[curr_surface==0]=2022
                    out_surface = np.minimum(out_surface,curr_surface)    
                else:
                    curr_surface = np.nan_to_num(curr_surface)
                    out_surface = np.maximum(out_surface,curr_surface)                 
                
                print (counter,' of 7725',target_variable,statistic_str,file,time.time()-starttime)
                
                #if counter==10:
                #    break
    
            if statistic_str=='min':            
                out_surface[out_surface==2022]=0     
                    
            gdalNumpy2floatRaster_compressed(np.rot90(out_surface),surface_folder+os.sep+'ES_buildings_%s_%s_%s.tif' %(target_variable,statistic_str,resample_factor),template_raster,cols,rows,bitdepth)
            #sys.exit(0)

################################################################################################
        
if rasterize_mutemp:

    params=[]
    params.append([True,False,False])
    params.append([False,True,False])
    params.append([False,False,True])
    
    for param in params:
        do_can_regcan=param[0]
        do_all_excpt_can_utm=param[1]        
        do_all_laea=param[2]
            
        ################################################# canaries regcan
        if do_can_regcan:
            template_raster = template_raster_regcan
            crs_grid = 4083 #epsg of template_raster
            surface_folder = rasterdir_can_regcan
            do_canaries_only=True
        ################################################# iberic pen and canaries laea,
        if do_all_laea:
            template_raster = template_raster_laea 
            crs_grid = 3035 #epsg of template_raster
            surface_folder = rasterdir_laea
            do_canaries_only=False
        ################################################# iberic pen and bal, mel, ceu, UTM
        if do_all_excpt_can_utm:
            template_raster = template_raster_utm 
            crs_grid = 25830 #epsg of template_raster
            surface_folder = rasterdir_utm
            do_canaries_only=False
                   
        if do_canaries_only:    
            can_files=[]
            for file in os.listdir(shp_dir):
                if file.split('.')[-1]=='shp':
                    if '.BU.35' in file or '.BU.38' in file:
                        can_files.append(['CAN',shp_dir+os.sep+file])
                    else:
                        can_files.append(['IB',shp_dir+os.sep+file])                    
            can_files_df=pd.DataFrame(can_files,columns=['region','fullpath'])
            can_idx=can_files_df[can_files_df.region=='CAN'].index.values
    
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
     
        for year in years:
            
            if year==1900:            
                lower_year=1
            else:
                lower_year=year-5
                
            target_variable='count'
            statistic=np.sum
            statistic_str='sum_%s_%s' %(year,lower_year)
            bitdepth=gdal.GDT_Int16
            
            out_surface =np.zeros((cols,rows)).astype(np.float32)
            out_surface_area =np.zeros((cols,rows)).astype(np.float32)
                    
            counter=0
            for file in os.listdir(shp_pt_dir):
                
                if not file.split('.')[-1]=='shp':
                    continue            
    
                if do_canaries_only:
                    curridx=file.replace('.shp','').split('_')[-1]
                    if not int(curridx) in can_idx:
                        continue
                    
                #if not 'GIPUZKOA' in file:
                #    continue
                starttime=time.time()
    
                indf = gp.read_file(shp_pt_dir+os.sep+file)
                
                if len(indf)==0:
                    continue
                
                indf = indf[np.logical_and(indf['yearbuilt']<=year,indf.yearbuilt>lower_year)]
    
                if len(indf)==0:
                    continue
                                                    
                indf.geometry = indf.geometry.to_crs(epsg=crs_grid)      
                indf[xcoo_col]=indf.geometry.x
                indf[ycoo_col]=indf.geometry.y
                                
                if target_variable == 'count':
                    indf[target_variable]=1
                        
                #################################################################
             
                counter+=1
    
                indf = indf.dropna(subset=[target_variable])
                indf=indf[[xcoo_col,ycoo_col,target_variable,'area']]                   
                statsvals = indf[target_variable].values.astype(np.int32)            
                curr_surface = scipy.stats.binned_statistic_2d(indf[xcoo_col].values,indf[ycoo_col].values,statsvals,statistic,bins=[cols,rows],range=rasterrange).statistic     
                curr_surface = np.nan_to_num(curr_surface)
                out_surface = np.maximum(out_surface,curr_surface)  
    
                ### bldg area:
                statsvals = indf['area'].values.astype(np.int32)          
                curr_surface = scipy.stats.binned_statistic_2d(indf[xcoo_col].values,indf[ycoo_col].values,statsvals,statistic,bins=[cols,rows],range=rasterrange).statistic     
                curr_surface = np.nan_to_num(curr_surface)
                out_surface_area = np.maximum(out_surface_area,curr_surface)  
    
                print (counter,' of 7725',target_variable,statistic_str,file,time.time()-starttime)
    
            gdalNumpy2floatRaster_compressed(np.rot90(out_surface),surface_folder+os.sep+'ES_buildings_%s_%s_%s_increment.tif' %(target_variable,statistic_str,resample_factor),template_raster,cols,rows,bitdepth)
            gdalNumpy2floatRaster_compressed(np.rot90(out_surface_area),surface_folder+os.sep+'ES_buildings_%s_%s_%s_increment.tif' %('area',statistic_str,resample_factor),template_raster,cols,rows,bitdepth)
             
        ##### now merge the increments to cumulative counts:
                
        for year in years:        
            
            if year ==1900:
                lower_year=1
            else:
                lower_year=year-5
            
            target_variable='count'
            statistic_str='sum_%s_%s' %(year,lower_year)                
            curr_count_incr = surface_folder+os.sep+'ES_buildings_%s_%s_%s_increment.tif' %(target_variable,statistic_str,resample_factor)
            statistic_str='sum_%s' %(year)                        
            curr_count_cum = surface_folder+os.sep+'ES_buildings_%s_%s_%s_cum.tif' %(target_variable,statistic_str,resample_factor)
            
            target_variable='area'
            statistic_str='sum_%s_%s' %(year,lower_year)        
            curr_area_incr = surface_folder+os.sep+'ES_buildings_%s_%s_%s_increment.tif' %('area',statistic_str,resample_factor)
            statistic_str='sum_%s' %(year)                               
            curr_area_cum = surface_folder+os.sep+'ES_buildings_%s_%s_%s_cum.tif' %('area',statistic_str,resample_factor)
            statistic_str='%s' %(year)                                       
            curr_bua_cum = surface_folder+os.sep+'ES_buildings_%s_%s_%s_cum.tif' %('bua',statistic_str,resample_factor)
            
            curr_count_incr_arr=gdal.Open(curr_count_incr).ReadAsArray()
            curr_area_incr_arr=gdal.Open(curr_area_incr).ReadAsArray()
            
            if year==1900:
                total_count_surface=curr_count_incr_arr.copy()
                total_area_surface=curr_area_incr_arr.copy()
            else:
                total_count_surface=total_count_surface+curr_count_incr_arr
                total_area_surface=total_area_surface+curr_area_incr_arr
                
            total_area_surface[total_area_surface<0]=0
            total_area_surface[total_area_surface>resample_factor*resample_factor]=resample_factor*resample_factor # max bufa should not exceed surface area (will be chopped off at 10k sqm). this was set to 0 in the first version (copy-paste mistake).
            
            bua_surface = total_count_surface.copy()
            bua_surface[bua_surface>1]=1
                
            gdalNumpy2floatRaster_compressed(total_count_surface,curr_count_cum,template_raster,cols,rows,bitdepth)
            gdalNumpy2floatRaster_compressed(total_area_surface,curr_area_cum,template_raster,cols,rows,bitdepth)
            gdalNumpy2floatRaster_compressed(bua_surface,curr_bua_cum,template_raster,cols,rows,bitdepth)
    
            print('cumulative %s' %year)

################################################################################################

if mine_landuse:

    regions=[]
    for file in os.listdir(shp_pt_dir):
        regions.append(file.replace('es_bldg_harm_pt_','').replace(file.split('_')[-1],'')[:-1])
    regions=list(set(regions))
    
    for region in regions:
        lutypes_reg=[]
        sede_cad_count=0
        for file in os.listdir(shp_pt_dir):    
            if not file.split('.')[-1]=='shp':
                continue 
            
            if not region in file:
                continue
            
            if 'ES_INSPIRE' in file and sede_cad_count>30: ## we only need a small sample of sede cadastral data
                continue
            
            if 'ES_INSPIRE' in file:
                sede_cad_count+=1
            
            indf=gp.read_file(shp_pt_dir+os.sep+file)
            unique_lus=indf.usage.unique()
            if len(lutypes_reg)==0:
                lutypes_reg=list(unique_lus)
            else:
                lutypes_reg==list(set(lutypes_reg+list(unique_lus)))
        print(region,lutypes_reg)

################################################################################################
        
if harmonize_landuse:

    regions=[]
    for file in os.listdir(shp_pt_dir):
        regions.append(file.replace('es_bldg_harm_pt_','').replace(file.split('_')[-1],'')[:-1])
    regions=list(set(regions))
    
    lu_mapping_df=pd.read_csv(lu_mapping_csv)
    
    for region in regions:
        no_lu=False
        curr_lu_mapping_df=lu_mapping_df[lu_mapping_df.region==region]
        if len(curr_lu_mapping_df)==0:
            no_lu=True
        else:
            curr_lu_mapping_dict=dict(curr_lu_mapping_df[['source_lu','target_lu']].values)       
        
        for file in os.listdir(shp_pt_dir):    
            if not file.split('.')[-1]=='shp':
                continue            
            if not region in file:
                continue
            
            indf=gp.read_file(shp_pt_dir+os.sep+file)
            if no_lu:
                indf['lu_harm']=' '
            else:
                indf['lu_harm']=indf.usage.map(curr_lu_mapping_dict)
            indf.to_file(shp_pt_dir_lu_harm+os.sep+file)
            print('lu harmonized',file)
                      
################################################################################################
            
if rasterize_landuse_mutemp:

    target_lus=['agriculture', 'publicservices', 'other', 'residential',
       'industrial', 'office', 'commercial', ' ']
  
    params=[]
    params.append([True,False,False])
    params.append([False,True,False])
    params.append([False,False,True])
    
    for param in params:
        do_can_regcan=param[0]
        do_all_excpt_can_utm=param[1]        
        do_all_laea=param[2]
            
        ################################################# canaries regcan
        if do_can_regcan:
            template_raster = template_raster_regcan
            crs_grid = 4083 #epsg of template_raster
            surface_folder = rasterdir_can_regcan
            do_canaries_only=True
        ################################################# iberic pen and canaries laea,
        if do_all_laea:
            template_raster = template_raster_laea 
            crs_grid = 3035 #epsg of template_raster
            surface_folder = rasterdir_laea
            do_canaries_only=False
        ################################################# iberic pen and bal, mel, ceu, UTM
        if do_all_excpt_can_utm:
            template_raster = template_raster_utm 
            crs_grid = 25830 #epsg of template_raster
            surface_folder = rasterdir_utm
            do_canaries_only=False
        
        if do_canaries_only:    
            can_files=[]
            for file in os.listdir(shp_dir):
                if file.split('.')[-1]=='shp':
                    if '.BU.35' in file or '.BU.38' in file:
                        can_files.append(['CAN',shp_dir+os.sep+file])
                    else:
                        can_files.append(['IB',shp_dir+os.sep+file])                    
            can_files_df=pd.DataFrame(can_files,columns=['region','fullpath'])
            can_idx=can_files_df[can_files_df.region=='CAN'].index.values
            
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
        
        for year in years:#####################
        
            if year==1900:            
                lower_year=1
            else:
                lower_year=year-5
            
            out_dfs=[]
            for target_lu in target_lus:
                out_df =pd.DataFrame()
                out_dfs.append([out_df.copy(),target_lu])
                      
            counter=0
            for file in os.listdir(shp_pt_dir_lu_harm):
                
                if not file.split('.')[-1]=='shp':
                    continue
                
                if do_canaries_only:
                    curridx=file.replace('.shp','').split('_')[-1]
                    if not int(curridx) in can_idx:
                        continue
              
                counter+=1            
    
                indf = gp.read_file(shp_pt_dir_lu_harm+os.sep+file)
                if len(indf)==0:
                    continue
                
                indf = indf[np.logical_and(indf['yearbuilt']<=year,indf.yearbuilt>lower_year)]
    
                if len(indf)==0:
                    continue               
                
                #indf = indf[['geometry','lu_harm']]
                indf.lu_harm=indf.lu_harm.replace(np.nan,' ').replace('',' ').replace('None',' ')
                
                for target_lu in target_lus:
                    starttime=time.time()  
                    target_variable='lu_%s_count' %target_lu
                                 
                    currindf=indf[indf.lu_harm==target_lu]
                    currindf=gp.GeoDataFrame(currindf)
    
                    if not len(currindf)==0:
                                
                        currindf['geometry'] = currindf.geometry.to_crs(epsg=crs_grid)      
                        currindf[xcoo_col]=currindf.geometry.x
                        currindf[ycoo_col]=currindf.geometry.y
                                    
                        currindf[target_variable]=1
        
                        currindf = currindf.dropna(subset=[target_variable])
                        currindf=currindf[[xcoo_col,ycoo_col,target_variable]] 
                        totaldf=out_dfs[target_lus.index(target_lu)][0]
                        totaldf=totaldf.append(currindf)
                        out_dfs[target_lus.index(target_lu)][0]=totaldf.copy()  

                    print (year, counter,' of 7725',param,target_variable,file,time.time()-starttime)

            for target_lu in target_lus:#####################
                print('binning', target_lu, year)
                starttime=time.time() 
                statistic=np.nansum
                statistic_str=''
                bitdepth=gdal.GDT_Int16
                target_variable='lu_%s_count' %target_lu                
                curr_lu_df=out_dfs[target_lus.index(target_lu)][0]  
                if len(curr_lu_df)==0:
                    continue
                statsvals = curr_lu_df[target_variable].values.astype(np.int32)                 
                curr_surface = scipy.stats.binned_statistic_2d(curr_lu_df[xcoo_col].values,curr_lu_df[ycoo_col].values,statsvals,statistic,bins=[cols,rows],range=rasterrange).statistic                     
                curr_surface = np.nan_to_num(curr_surface)
                curr_lu=target_lu.replace(' ','')
                gdalNumpy2floatRaster_compressed(np.rot90(curr_surface),surface_folder+os.sep+'ES_buildings_%s_%s_%s_%s_%s_%s.tif' %('landuse','count_mutemp',curr_lu,resample_factor,lower_year,year),template_raster,cols,rows,bitdepth)
                print('done binning', target_lu, year, time.time()-starttime)
                            
        ##### now merge the increments to cumulative counts:
        for target_lu in target_lus:#####################
            curr_lu=target_lu.replace(' ','')
                    
            for year in years:        
                
                if year ==1900:
                    lower_year=1
                else:
                    lower_year=year-5
                    
                inraster = surface_folder+os.sep+'ES_buildings_%s_%s_%s_%s_%s_%s.tif' %('landuse','count_mutemp',curr_lu,resample_factor,lower_year,year)
                try:
                    curr_incr_arr=gdal.Open(inraster).ReadAsArray()
                except:
                    ## if not exist, no records in current time slice. we add 0
                    curr_incr_arr=np.zeros(total_count_surface.shape)
                
                if year==1900:
                    total_count_surface=curr_incr_arr.copy()
                else:
                    total_count_surface=total_count_surface+curr_incr_arr
                    
                total_count_surface[total_count_surface<0]=0
                
                curr_count_cum = surface_folder+os.sep+'ES_buildings_%s_%s_%s_%s_%s.tif' %('landuse','count_mutemp',curr_lu,resample_factor,year)
                gdalNumpy2floatRaster_compressed(total_count_surface,curr_count_cum,template_raster,cols,rows,bitdepth)
                print(target_lu, 'cumulative %s' %year)

################################################################################################
    
if rasterize_physical_characteristics:

    target_vars=['area','offi_area','num_dwel','num_bunits']
    stats=[np.nansum,np.nanmean]
    stats_str=['sum','mean']
    
    params=[]
    params.append([True,False,False])
    params.append([False,True,False])
    params.append([False,False,True])
    
    for param in params:
        do_can_regcan=param[0]
        do_all_excpt_can_utm=param[1]        
        do_all_laea=param[2]
            
        ################################################# canaries regcan
        if do_can_regcan:
            template_raster = template_raster_regcan
            crs_grid = 4083 #epsg of template_raster
            surface_folder = rasterdir_can_regcan
            do_canaries_only=True
        ################################################# iberic pen and canaries laea,
        if do_all_laea:
            template_raster = template_raster_laea 
            crs_grid = 3035 #epsg of template_raster
            surface_folder = rasterdir_laea
            do_canaries_only=False
        ################################################# iberic pen and bal, mel, ceu, UTM
        if do_all_excpt_can_utm:
            template_raster = template_raster_utm 
            crs_grid = 25830 #epsg of template_raster
            surface_folder = rasterdir_utm
            do_canaries_only=False
        
        if do_canaries_only:    
            can_files=[]
            for file in os.listdir(shp_dir):
                if file.split('.')[-1]=='shp':
                    if '.BU.35' in file or '.BU.38' in file:
                        can_files.append(['CAN',shp_dir+os.sep+file])
                    else:
                        can_files.append(['IB',shp_dir+os.sep+file])                    
            can_files_df=pd.DataFrame(can_files,columns=['region','fullpath'])
            can_idx=can_files_df[can_files_df.region=='CAN'].index.values
        
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
        
        target_var_count=-1
        for target_var in target_vars:
            target_var_count+=1
            
            out_surfaces=[]
            for stat_str in stats_str:
                out_surface =np.zeros((cols,rows)).astype(np.float32)
                out_surfaces.append([out_surface.copy(),'%s_%s' %(target_var,stat_str)])
                      
            counter=0
            for file in os.listdir(shp_pt_dir_lu_harm):
                
                if not file.split('.')[-1]=='shp':
                    continue

                
                if do_canaries_only:
                    curridx=file.replace('.shp','').split('_')[-1]
                    if not int(curridx) in can_idx:
                        continue
              
                counter+=1            
    
                indf = gp.read_file(shp_pt_dir_lu_harm+os.sep+file)
                if len(indf)==0:
                    continue

                indf[target_var]=indf[target_var].replace(np.nan,' ').replace('',' ').replace('None',' ').replace(0,' ').replace('0',' ')
                
                for stat_str in stats_str:
                    starttime=time.time()  
                    target_variable=target_var
                    statistic=stats[stats_str.index(stat_str)]
                    statistic_str=stat_str
                    bitdepth=gdal.GDT_Float32
                    curr__surface=out_surfaces[stats_str.index(stat_str)][0]
                    out_str=out_surfaces[stats_str.index(stat_str)][1]
                    
                    if stat_str=='missing':
                        currindf=indf[indf[target_var]==' ']
                        target_variable='count'
                        currindf['count']=1
                    else:
                        currindf=indf[indf[target_var]!=' ']
                    
                    try:
                        currindf=currindf[currindf[target_var]>=0]
                    except:
                        pass
                    
                    currindf=gp.GeoDataFrame(currindf)
    
                    if not len(currindf)==0:
                                
                        currindf['geometry'] = currindf.geometry.to_crs(epsg=crs_grid)      
                        currindf[xcoo_col]=currindf.geometry.x
                        currindf[ycoo_col]=currindf.geometry.y
                                    
                        #currindf[target_variable]=1
                        
                        try:
                            currindf[target_variable]=currindf[target_variable].str.replace(',','.')
                            currindf[target_variable]=currindf[target_variable].map(float)                        
                        except:
                            pass                        
                                
                        currindf = currindf.dropna(subset=[target_variable])
                        currindf=currindf[[xcoo_col,ycoo_col,target_variable]]                   
                        statsvals = currindf[target_variable].values.astype(np.int32)                 
                        curr_surface = scipy.stats.binned_statistic_2d(currindf[xcoo_col].values,currindf[ycoo_col].values,statsvals,statistic,bins=[cols,rows],range=rasterrange).statistic     
                        
                        curr_surface = np.nan_to_num(curr_surface)
                        curr__surface = np.maximum(curr__surface,curr_surface) 
                        out_surfaces[stats_str.index(stat_str)][0]=curr__surface.copy()           
    
                    print (counter,' of 7725',param,out_str,file,time.time()-starttime)
                    #if counter==25:
                    #    break
                
                    if counter%20==0:
                        ### backup surface
                        np.savez_compressed(bakdir+os.sep+'backup_arr_physprop_%s_%s.npz' %(out_str,counter), a=curr__surface)
                        ###delete previous versions
                        for x in os.listdir(bakdir):
                            if out_str in x:
                                if not int(x.split('_')[-1].split('.')[0])==counter:
                                    try:
                                        os.remove(bakdir+os.sep+x)
                                    except:
                                        pass
    
                        
            ###########################################################
            for out_surface in out_surfaces:
                curr_outsurface=out_surface[0]
                curr_str=out_surface[1]
                gdalNumpy2floatRaster_compressed(np.rot90(curr_outsurface),surface_folder+os.sep+'ES_buildings_%s_%s_%s_%s.tif' %('physprop','count',curr_str,resample_factor),template_raster,cols,rows,bitdepth)
                            
