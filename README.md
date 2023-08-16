# Historical Settlement Data Compilation for Spain (HISDAC-ES) from 1900 to 2020

The Historical Settlement Data Compilation for Spain (HISDAC-ES) contains a range of gridded surfaces in GeoTIFF format, describing physical, temporal, landuse, and multitemporal characteristics of the built environment in Spain from 1900 to 2020. HISDAC-ES has been created based on publicly available cadastral data (i.e., building footprint data + thematic attributes). The spatial resolution of the raster data is 100m, the temporal resolution of raster time series is 5 years. The dataset is available at http://doi.org/10.6084/m9.figshare.22009643. The Python code in HISDAC-ES_data_production.py allows for recreating the dataset, and provides a numpy/scipy/gdal/geopandas based strategy for efficient, large-scale geospatial data gridding (raster aggregation).

#### The 240+ gridded datasets in HISDAC-ES include the DEVA layer series, highlighting historical extents of developed areas:
<img width="450" src="https://github.com/johannesuhl/hisdac-es/blob/main/hisdac_es_deva.gif">

More visualizations are available at https://doi.org/10.6084/m9.figshare.22064798.

Dataset:

Uhl, Johannes H.; Royé, Dominic; Burghardt, Keith; Aldrey Vázquez, José Antonio; Borobio Sanchiz, Manuel; Leyk, Stefan (2023): HISDAC-ES: Historical Settlement Data Compilation for Spain (1900-2020). figshare. Dataset. https://doi.org/10.6084/m9.figshare.22009643

Visualizations:

Uhl, Johannes H.; Royé, Dominic; Burghardt, Keith; Aldrey Vázquez, José Antonio; Borobio Sanchiz, Manuel; Leyk, Stefan (2023): Visualizing long-term urbanization and land development in Spain (1900-2020). figshare. Media. https://doi.org/10.6084/m9.figshare.22064798 

Code usage notes:

# UPDATE 08-2023: We published two new scripts, reflecting changes and additions to the data, accelerating the rasterization process significantly.
1) Script HISDAC-ES_data_production.py : We fixed two bugs that caused an incorrect rasterization of the multitemporal BUFA layers and of the PHYS layers (BUFA, BIA, DWEL, BUNITS sum and mean).
2) Script: HISDAC-ES_data_aggregation_municipalities.py :
     - Creates a harmonized building centroid vector dataset (output of block "spatial_join", which is provided in the figshare repository)
     - performs a spatial join with municipality polygons
     - Calculates building and completeness statistics per municipality
     - Outputs them as CSV and GPKG format
3) Script: HISDAC-ES_residential_layers_fast.py : Contains a revised, very fast method for rasterization, exemplified for the new residential building footprint and building indoor area layers (RESBIA, RESBUFA). While in HISDAC-ES_data_production.py, we loop through the buildings per municipality and add the raster cells per municipality to a "master" grid, in this new version, we use the country-level building centroid GeoPackage to perform the rasterization much faster, but also much more memory-intensive.







