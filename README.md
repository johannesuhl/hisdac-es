# Historical Settlement Data Compilation for Spain (HISDAC-ES) from 1900 to 2020

The Historical Settlement Data Compilation for Spain (HISDAC-ES) contains a range of gridded surfaces in GeoTIFF format, describing physical, temporal, landuse, and multitemporal characteristics of the built environment in Spain from 1900 to 2020. HISDAC-ES has been created based on publicly available cadastral data (i.e., building footprint data + thematic attributes). The spatial resolution of the raster data is 100m, the temporal resolution of raster time series is 5 years. The dataset is available at http://doi.org/10.6084/m9.figshare.22009643. Code in this repository allows for recreating the dataset, and provides a numpy/scipy/gdal/geopandas based strategy for efficient, large-scale geospatial data gridding (raster aggregation).

#### One of the 240 gridded datasets is the DEVA layer, showing historical extents of developed areas:
<img width="450" src="https://github.com/johannesuhl/hisdac-es/blob/main/hisdac_es_deva.gif">
