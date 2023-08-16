# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 09:14:41 2021

@author: Johannes H. Uhl, University of Colorado Boulder, USA.
"""

import os,sys
import geopandas as gp
import pandas as pd
import numpy as np

#municipality boundaries need to be obtained from https://centrodedescargas.cnig.es/CentroDescargas/index.jsp:
muni_shp1 = './lineas_limite2023/SHP_ETRS89/recintos_municipales_inspire_peninbal_etrs89/recintos_municipales_inspire_peninbal_etrs89.shp'
muni_shp2 = './lineas_limite2023/SHP_REGCAN95/recintos_municipales_inspire_canarias_regcan95/recintos_municipales_inspire_canarias_regcan95.shp'
harm_pt_shp_dir='ES_buildings_shp_pt_lu_harm' # folder containing output from
outdir='MUNI_STATS' # folder for the outputs

###################################################################

prepare_muni_data=True #prepare municipality polygon data
merge_shps=True #create a country-wide shapefile of building centroids
spatial_join=True #append municipality ID to the building centroids, export as gpkg
muni_stats=True #create statistics for each municipality, export as csv
exp_muni_stats_gpkg=True #attache csv to gpkg
miss_stats=True #create completeness statistics per municipality, export as csv

###################################################################

years=np.arange(1900,2021,1)

if prepare_muni_data:
    muni_gdf1=gp.read_file(muni_shp1, encoding='ISO-8859-1')
    muni_gdf2=gp.read_file(muni_shp2, encoding='ISO-8859-1')
    muni_gdf1=muni_gdf1.to_crs(epsg=3035)
    muni_gdf2=muni_gdf2.to_crs(epsg=3035)
    gdf=muni_gdf1.append(muni_gdf2)
    gdf=gdf.reset_index()
    munigdf=gdf[['geometry','NATCODE','CODNUT1','CODNUT2','CODNUT3']]
    munigdf['CA_CODE']=munigdf.NATCODE.str.slice(2,4)
    munigdf['PROV_CODE']=munigdf.NATCODE.str.slice(4,6)   
    munigdf['LAU_CODE']=munigdf.NATCODE.str.slice(6,11)   
    munigdf.columns=['geometry', 'NATCODE', 'NUTS1', 'NUTS2', 'NUTS3', 'CA_CODE','PROV_CODE', 'LAU_CODE']
    munigdf.to_file(outdir+os.sep+'ES_municipalities_merged.shp')

if merge_shps: #####################################################

    count=0
    allgdf=pd.DataFrame()
    for shp in os.listdir(harm_pt_shp_dir):
        if not shp.split('.')[-1]=='shp':
            continue
        count+=1
        ingdf=gp.read_file(harm_pt_shp_dir+os.sep+shp)
        if not ingdf.crs.to_epsg()==3035:
            ingdf.geometry=ingdf.geometry.to_crs(epsg=3035)
            print('projected')
        ingdf.x=ingdf.geometry.centroid.x
        ingdf.y=ingdf.geometry.centroid.y
        allgdf=allgdf.append(ingdf[['yearbuilt','area','num_floors','num_dwel','num_bunits','offi_area','x','y','lu_harm']])
        print(count,shp)
    allgdf=gp.GeoDataFrame(allgdf,geometry=gp.points_from_xy(allgdf.x.values, allgdf.y.values))    
    allgdf.set_crs(epsg=3035)
    allgdf.to_file(outdir+os.sep+'ES_building_centroids_merged.shp')    

if spatial_join: #####################################################

    munigdf = gp.read_file(outdir+os.sep+'ES_municipalities_merged.shp')
    allgdf = gp.read_file(outdir+os.sep+'ES_building_centroids_merged.shp')
    source_crs =munigdf.geometry.crs
    target_crs =allgdf.geometry.crs
    if not source_crs==target_crs:
        munigdf.geometry=munigdf.geometry.to_crs(target_crs)
    print('spatial join...')
    allgdf = gp.sjoin(allgdf,munigdf, how="left")
    print('exporting...')
    allgdf.to_file(outdir+os.sep+'ES_building_centroids_merged_spatjoin.gpkg')    
    
if muni_stats: #####################################################
    
    allgdf=gp.read_file(outdir+os.sep+'ES_building_centroids_merged_spatjoin.gpkg') 
    munistats=[]
    counter=0
    for muni,munidf in allgdf.groupby('NATCODE'):
        counter+=1
        
        munidf=munidf.replace('None',np.nan)
        munidf=munidf.fillna(0)
        
        munidf.yearbuilt=munidf.yearbuilt.fillna(0)
        munidf.yearbuilt=munidf.yearbuilt.map(int)
        munidf['area']=munidf['area'].map(float)
        munidf['num_dwel']=munidf['num_dwel'].map(str).str.replace(',','.').map(float)
        munidf['num_bunits']=munidf['num_bunits'].map(str).str.replace(',','.').map(float)
        munidf['offi_area']=munidf['offi_area'].map(str).str.replace(',','.').map(float)
        
        for year in years:            
            muniyrdf=munidf[np.logical_and(munidf.yearbuilt>0,munidf.yearbuilt<=year)]
            numbldgs=len(muniyrdf)
            if not numbldgs==0:
                areasum=np.nansum(muniyrdf['area'].values)
                dwelsum=np.nansum(muniyrdf.num_dwel.values)
                bunitsum=np.nansum(muniyrdf.num_bunits.values)
                indoorareasum=np.nansum(muniyrdf.offi_area.values)
                resdf = muniyrdf[muniyrdf.lu_harm=='residential']
                areasum_residential = np.nansum(resdf['area'].values)
                indoorareasum_residential = np.nansum(resdf['offi_area'].values)
                munistats.append([muni,year,numbldgs,areasum,dwelsum,bunitsum,indoorareasum,areasum_residential,indoorareasum_residential])
                print(counter,muni,year,numbldgs,areasum,dwelsum,bunitsum,indoorareasum,areasum_residential,indoorareasum_residential)
                
        munimissdf=munidf[munidf.yearbuilt==0]
        numbldgs=len(munimissdf)
        if not numbldgs==0:
            areasum=np.nansum(munimissdf['area'].values)
            dwelsum=np.nansum(munimissdf.num_dwel.values)
            bunitsum=np.nansum(munimissdf.num_bunits.values)
            indoorareasum=np.nansum(munimissdf.offi_area.values)
            resdf = munimissdf[munimissdf.lu_harm=='residential']
            areasum_residential = np.nansum(resdf['area'].values)
            indoorareasum_residential = np.nansum(resdf['offi_area'].values)            
            munistats.append([muni,0,numbldgs,areasum,dwelsum,bunitsum,indoorareasum,areasum_residential,indoorareasum_residential])
            print(counter,muni,0,numbldgs,areasum,dwelsum,bunitsum,indoorareasum,areasum_residential,indoorareasum_residential,areasum_residential,indoorareasum_residential)
        else:    
            areasum=0
            dwelsum=0
            bunitsum=0
            indoorareasum=0
            areasum_residential=0 
            indoorareasum_residential=0
            munistats.append([muni,0,numbldgs,areasum,dwelsum,bunitsum,indoorareasum,areasum_residential,indoorareasum_residential])
            print(counter,muni,0,numbldgs,areasum,dwelsum,bunitsum,indoorareasum,areasum_residential,indoorareasum_residential)       

    munistatsdf=pd.DataFrame(munistats)
    munistatsdf.columns=['NATCODE','year','numbldgs','areasum','dwelsum','bunitsum','indoorareasum','areasum_residential','indoorareasum_residential']
    munistatsdf['CA_CODE']=munistatsdf.NATCODE.str.slice(2,4)
    munistatsdf['PROV_CODE']=munistatsdf.NATCODE.str.slice(4,6)   
    munistatsdf['LAU_CODE']=munistatsdf.NATCODE.str.slice(6,11)      
    munistatsdf.to_csv(outdir+os.sep+'ES_building_muni_stats_v2.csv',index=False)


if exp_muni_stats_gpkg: #####################################################

    munistatsdf=pd.read_csv(outdir+os.sep+'ES_building_muni_stats_v2.csv')
    gdf = gp.read_file(outdir+os.sep+'ES_municipalities_merged.shp')
    years_red=np.arange(1900,2021,10)
    gdf=gdf.merge(munistatsdf,left_on='NATCODE',right_on='NATCODE',how='right')
    gdf=gdf[gdf.year.isin(years_red)]
    gdf['building_dens']=np.divide(gdf['numbldgs'],gdf.geometry.area.values)
    gdf['area_dens']=np.divide(gdf['areasum'],gdf.geometry.area.values)
    gdf['offi_area_dens']=np.divide(gdf['indoorareasum'],gdf.geometry.area.values)
    gdf.to_file(outdir+os.sep+'ES_building_muni_stats_exp.gpkg',driver='GPKG')
            
if miss_stats: #####################################################
    
    allgdf=gp.read_file(outdir+os.sep+'ES_building_centroids_merged_spatjoin.gpkg')
    munistats=[]
    counter=0
    for muni,munidf in allgdf.groupby('NATCODE'):
        counter+=1
        munidf.yearbuilt=munidf.yearbuilt.fillna(0)
        munidf.yearbuilt=munidf.yearbuilt.map(int)
        munidf['area']=munidf['area'].map(float)
        munidf['num_floors']=munidf['num_floors'].map(str).str.replace(',','.').replace('None','0').map(float)
        munidf['num_dwel']=munidf['num_dwel'].map(str).str.replace(',','.').replace('None','0').map(float)
        munidf['num_bunits']=munidf['num_bunits'].map(str).str.replace(',','.').replace('None','0').map(float)
        munidf['offi_area']=munidf['offi_area'].map(str).str.replace(',','.').replace('None','0').map(float)
        munidf=munidf.replace('None',np.nan)
        munidf=munidf.replace('',np.nan)
        munidf=munidf.fillna(0)
        num_total=len(munidf)
        res_munidf = munidf[munidf.lu_harm=='residential']
        nres_munidf = munidf[munidf.lu_harm!='residential']
        num_total_residential = len(res_munidf)
        num_total_nonresidential = num_total - num_total_residential
        
        if num_total==0:
            continue
        
        munidf_valby = munidf[munidf.yearbuilt>0]
        try:
            minby=np.nanmin(munidf_valby.yearbuilt.values)
        except:
            minby=np.nan
        try:
            maxby=np.nanmax(munidf_valby.yearbuilt.values)
        except:
            maxby=np.nan               
        
        prop_bymiss=100*len(munidf[np.logical_not(munidf.yearbuilt>0)])/float(num_total)
        prop_lumiss=100*len(munidf[munidf.lu_harm==0])/float(num_total)
        prop_luother=100*len(munidf[munidf.lu_harm=='other'])/float(num_total)
        prop_num_floors_miss=100*len(munidf[np.logical_not(munidf.num_floors>0)])/float(num_total)        
        prop_num_dwel_miss=100*len(munidf[np.logical_not(munidf.num_dwel>0)])/float(num_total)
        if num_total_residential>0:
            prop_num_dwel_miss_res=100*len(res_munidf[np.logical_not(res_munidf.num_dwel>0)])/float(num_total_residential)         
        else:
            prop_num_dwel_miss_res=np.nan
        prop_num_bunits_miss=100*len(munidf[np.logical_not(munidf.num_bunits>0)])/float(num_total)
        if num_total_nonresidential>0:
            prop_num_bunits_miss_nores=100*len(nres_munidf[np.logical_not(nres_munidf.num_bunits>0)])/float(num_total_nonresidential)           
        else:
            prop_num_bunits_miss_nores=np.nan
        prop_offi_area_miss=100*len(munidf[np.logical_not(munidf.offi_area>0)]) /float(num_total)       
        prop_num_dwel_and_num_bunits_miss=100*len(munidf[np.logical_not(np.logical_and(munidf.num_dwel>0,munidf.num_bunits>0))])/float(num_total)  
        print([muni,num_total,prop_bymiss,prop_lumiss,prop_luother,prop_num_floors_miss,prop_num_dwel_miss,prop_num_bunits_miss,prop_offi_area_miss,prop_num_dwel_and_num_bunits_miss,prop_num_dwel_miss_res,prop_num_bunits_miss_nores,minby,maxby])
        munistats.append([muni,num_total,prop_bymiss,prop_lumiss,prop_luother,prop_num_floors_miss,prop_num_dwel_miss,prop_num_bunits_miss,prop_offi_area_miss,prop_num_dwel_and_num_bunits_miss,prop_num_dwel_miss_res,prop_num_bunits_miss_nores,minby,maxby])
    munistatsdf=pd.DataFrame(munistats)
    munistatsdf.columns=['NATCODE','num_total','perc_bymiss','perc_lumiss','perc_luother','perc_num_floors_miss','perc_num_dwel_miss','perc_num_bunits_miss','perc_offi_area_miss','perc_num_dwel_and_num_bunits_miss','prop_num_dwel_miss_res','prop_num_bunits_miss_nores','minby','maxby']
    munistatsdf['CA_CODE']=munistatsdf.NATCODE.str.slice(2,4)
    munistatsdf['PROV_CODE']=munistatsdf.NATCODE.str.slice(4,6)   
    munistatsdf['LAU_CODE']=munistatsdf.NATCODE.str.slice(6,11)          
    munistatsdf.to_csv(outdir+os.sep+'ES_building_muni_miss_stats.csv',index=False)  
    
