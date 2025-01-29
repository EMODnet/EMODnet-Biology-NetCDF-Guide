# NetCDF files for biodiversity products

| Date | Version | Authors | History |
| --- | --- | --- | --- |
| 2025-01-29 | V 0.1 | Joana Beja, Frederic Leclerq, Salva Fernandez | Document created |

# Contents

[Introduction 2](#_Toc119918628)

[Why Netcdf? 2](#_Toc119918629)

[ERDDAP constraints 2](#_Toc119918630)

[Python 3](#_Toc119918631)

[R 4](#_Toc119918632)

[CF Guidance 5](#_Toc119918633)

[Creating NetCDF files for EMODnet Biology 7](#_Toc119918634)

## Introduction

This document is meant as an internal guide for EMODnet Biology WP3 partners relating to NetCDF file creation. It was drafted in the context of the EMODnet centralisation which caused all thematic lots adjust the way their products are being created and subsequently integrated in the EMODnet viewer.

In EMODnet Biology we follow the FAIR principles (Wilkinson, M _et al._ The FAIR Guiding Principles for scientific data management and stewardship. _Sci Data_ **3** , 160018 (2016) [https://doi.org/10.1038/sdata.2016.18](https://doi.org/10.1038/sdata.2016.18)) for all our activities/services, be it relating to metadata, data, software or products.

### Why Netcdf?

During the centralisation procedure, a choice was made on the type of systems that were going to be implemented for data/products sub-setting and download. All thematic lots are responsible for implementing and/or maintaining OGC web services (e.g. WMS, WFS, WCS) for their data/products. As some thematic lots used ERDDAP ([https://github.com/BobSimons/erddap](https://github.com/BobSimons/erddap)), and it is a tool suited for gridded and tabular datasets, a decision was made to setup an EMODnet ERDDAP instance to allow for the sub-setting/merging for some products (e.g. essentially from the Chemistry and Physics thematic lots). It was also clear that the way Biology was making its products available via a map viewer was not ideal, as we weren't providing access to all the data but only to selected taxa. To bypass this, we found that creating NetCDFs ready for ERDDAP upload was the best option. It aligns with the Central Portal requirements and there will not be a need for a dedicated/customised solution for our products.

### ERDDAP constraints

During the Biology products' centralisation, we found that ERDDAP has a few technical constraints that we cannot bypass, namely, it doesn't accept files that contain more points that exceed the limits of the array size JAVA. (I.e. 2^31 – 1) nor does it deal well with strings and integers in the same file, something we require as we need to have the taxa LSID and the taxa name in the products. It also does not allow duplicates (e.g. duplicate LSIDs) nor does it allow the time dimension to be non-monotonic. Some of these constraints are features, others are unresolved issues.

As the Central Portal's requirements for the NetCDF files are that they need to comply with CF guidance, we have found a way that works for our type of data.

This document provides examples on how to create generic NetCDFs with Python or R, but note that this is illustrative, as all products are different, you will need to adjust the methodologies so that you can create your file.

## Python

More info: [https://unidata.github.io/netcdf4-python/#creatingopeningclosing-a-netcdf-file](https://unidata.github.io/netcdf4-python/#creatingopeningclosing-a-netcdf-file) . The library netCDF4 is widely supported by Unidata. It works well with numpy arrays.

```python
import netCDF4
import numpy as np
from numpy.random import uniform
from datetime import datetime, timedelta
from cftime import num2date, date2num
 
# Define length dimensions
nlon = 10
nlat = 10
ntime = 2
 
# Create nc file and add dimensions
ncfile = netCDF4.Dataset('foo.nc', mode = 'w', format = 'NETCDF4') 
dim_lon = ncfile.createDimension('lon', nlon)
dim_lat = ncfile.createDimension('lat', nlon)
dim_time = ncfile.createDimension('time', None)
 
# Create variables and add attributes
var_lon = ncfile.createVariable('lon', "f4", ("lon",))
var_lon.units = "degrees_east"
var_lon.long_name = "Longitude"
var_lon.reference_datum = "geographical coordinates, WGS84 projection"
 
var_lat = ncfile.createVariable('lat', "f4", ("lat",))
var_lat.units = "degrees_north"
var_lat.long_name = "Latitude"
var_lat.reference_datum = "geographical coordinates, WGS84 projection"
 
var_time = ncfile.createVariable('time', "i8", ("time",))
var_time.units = "days since 1970-01-01 00:00:00.0"
var_time.calendar = "gregorian"
 
var_temp = ncfile.createVariable('temp', "f8", ("lon", "lat", "time",))
 
# Add dummy data to variables
var_lon[:] = np.linspace(2.2383, 3.3704, nlon) 
var_lat[:] = np.linspace(51.0893, 51.8761, nlat)
var_time[:] = date2num([datetime(2020, 3, 13), datetime(2020, 3, 14)], 
                       units = var_time.units, calendar = var_time.calendar)
var_temp[0:nlon, 0:nlat, 0:ntime] = uniform(size=(nlon, nlat, ntime))
 
# Close and check
ncfile.close()
!ncdump -c -t foo.nc
```

## R

More info: [https://journal.r-project.org/archive/2013-2/michna-woods.pdf](https://journal.r-project.org/archive/2013-2/michna-woods.pdf)

This example uses the package RNetCDF as it is the **lowest level interface** to deal with NetCDF in R. Two other packages with higher level interfaces are ncdf4 and tidync.

- RNetCDF: [https://cran.r-project.org/web/packages/RNetCDF/index.html](https://cran.r-project.org/web/packages/RNetCDF/index.html)
- ncdf4: [https://cran.r-project.org/web/packages/ncdf4/index.html](https://cran.r-project.org/web/packages/ncdf4/index.html)
- tidync: [https://docs.ropensci.org/tidync/](https://docs.ropensci.org/tidync/)

```R
library(RNetCDF)
 
# Define length dimensions
nlon = 10
nlat = 10
ntime = 2
 
# Create nc file and add dimensions
nc <- create.nc("foo.nc")
dim.def.nc(nc, dimname = "lon", 10)
dim.def.nc(nc, dimname = "lat", 10)
dim.def.nc(nc, dimname = "time", unlim = TRUE)
 
# Create variables
var.def.nc(nc, varname = "lon", vartype = "NC_DOUBLE", dimensions = "lon")
var.def.nc(nc, varname = "lat", vartype = "NC_DOUBLE", dimensions = "lat")
var.def.nc(nc, varname = "time", vartype = "NC_INT", dimensions = "time")
var.def.nc(nc, varname = "temperature", vartype = "NC_DOUBLE", dimensions = c("lon", "lat", "time"))
 
# Add attributes. those starting by _* are special values used by the netcdf C library
att.put.nc(nc, variable = "temperature", name = "_FillValue", type = "NC_DOUBLE", value = -99999.9)
att.put.nc(nc, variable = "time", name = "units", type = "NC_CHAR", value = "days since 1970-01-01 00:00:00")
 
# Create dummy values
mylon <- seq(2.2383, 3.3704, length.out = nlon) # one dimension
mylat <- seq(51.0893, 51.8761, length.out = nlat) # one dimension
 
mytime <- matrix(nrow = ntime, ncol = 6) # one dimension
mytime[1,] <- c(2020, 03, 13, 12, 00, 00)
mytime[2,] <- c(2020, 03, 14, 16, 00, 00)
mytime_ut <- utinvcal.nc("days since 1970-01-01 00:00:00", mytime) # this is a special type of time units used by the netcdf C library
 
mytemperature_data <- runif(nlon*nlat*ntime) # three dimensions
mytemperature <- array(mytemperature_data, dim = c(nlon, nlat, ntime))
 
# Sync changes and check
sync.nc(nc)
print.nc(nc)
 
# Add values to netcdf file
var.put.nc(nc, variable = "lon", data = mylon)
var.put.nc(nc, variable = "lat", data = mylat)
var.put.nc(nc, variable = "time", data = mytime_ut)
var.put.nc(nc, variable = "temperature", data = mytemperature)
 
# Close and check
sync.nc(nc)
print.nc(nc)
close.nc(nc)

```

In [https://emodnet.github.io/EMODnet-Biology-products-erddap-demo/](https://emodnet.github.io/EMODnet-Biology-products-erddap-demo/) you can find an real life example and specific guidance on how to create a NetCDF file for Biology products.

## CF Guidance

Climate and Forecast conventions are widely used in some science disciplines. The most recent version can be found via [https://cfconventions.org/cf-conventions/cf-conventions.html](https://cfconventions.org/cf-conventions/cf-conventions.html). Chapter 6.1.2 describes the guidance for taxon names and identifiers and how they should be structured in a NetCDF file.

Below is an example of how to create a CF compliant NetCDF file. The aim is that the file can be manipulated without users having to go back to the metadata record. The header, which contains the metadata pertinent to each file, should be quite extensive and clear to any user (i.e no acronyms, no overly technical language, …).

•	Create a CF compliant netcdf file (R)

```R
library(RNetCDF)
 
# Create file
nc <- create.nc("./data/cf_comp", format = "netcdf4")
 
# Global attributes
att.put.nc(nc, variable = "NC_GLOBAL", name = "Conventions", type = "NC_CHAR", value = "CF-1.12")
att.put.nc(nc, variable = "NC_GLOBAL", name = "title", type = "NC_CHAR", value = "Example of CF compliant NetCDF file")
att.put.nc(nc, variable = "NC_GLOBAL", name = "institution", type = "NC_CHAR", value = "Flanders Marine Institute (VLIZ)")
att.put.nc(nc, variable = "NC_GLOBAL", name = "source", type = "NC_CHAR", value = "Bio-Oracle_2.1")
att.put.nc(nc, variable = "NC_GLOBAL", name = "history", type = "NC_CHAR", value = paste(Sys.time(), "File created"))
att.put.nc(nc, variable = "NC_GLOBAL", name = "comment", type = "NC_CHAR", value = "Uses attributes recommended by http://cfconventions.org")
att.put.nc(nc, variable = "NC_GLOBAL", name = "references", type = "NC_CHAR", value = "https://bio-oracle.org")
 
# Longitude dimension
dim.def.nc(nc, dimname = "lon", dimlength = 4320)
var.def.nc(nc, varname = "lon", vartype = "NC_FLOAT", dimensions = "lon")
att.put.nc(nc, variable = "lon", name = "standard_name", type = "NC_CHAR", value = "longitude")
att.put.nc(nc, variable = "lon", name = "long_name", type = "NC_CHAR", value = "longitude")
att.put.nc(nc, variable = "lon", name = "units", type = "NC_CHAR", value = "degrees_east")
 
# Latitude dimension
dim.def.nc(nc, dimname = "lat", dimlength = 2160)
var.def.nc(nc, varname = "lat", vartype = "NC_FLOAT", dimensions = "lat")
att.put.nc(nc, variable = "lat", name = "standard_name", type = "NC_CHAR", value = "latitude")
att.put.nc(nc, variable = "lat", name = "long_name", type = "NC_CHAR", value = "latitude")
att.put.nc(nc, variable = "lat", name = "units", type = "NC_CHAR", value = "degrees_north")
 
# CRS variable
var.def.nc(nc, varname = "crs", vartype = "NC_CHAR", dimensions = NA)
att.put.nc(nc, variable = "crs", name = "grid_mapping_name", type = "NC_CHAR", value = "latitude_longitude")
att.put.nc(nc, variable = "crs", name = "long_name", type = "NC_CHAR", value = "CRS definition")
att.put.nc(nc, variable = "crs", name = "longitude_of_prime_meridian", type = "NC_DOUBLE", value = 0.)
att.put.nc(nc, variable = "crs", name = "semi_major_axis", type = "NC_DOUBLE", value = 6378137.)
att.put.nc(nc, variable = "crs", name = "inverse_flattening", type = "NC_DOUBLE", value = 298.257223563)
att.put.nc(nc, variable = "crs", name = "spatial_ref", type = "NC_CHAR", value = 'GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AXIS[\"Latitude\",NORTH],AXIS[\"Longitude\",EAST],AUTHORITY[\"EPSG\",\"4326\"]]')
att.put.nc(nc, variable = "crs", name = "GeoTransform", type = "NC_CHAR", value = '-180 0.08333333333333333 0 90 0 -0.08333333333333333 ')
 
# Any random variable
var.def.nc(nc, varname = "data", vartype = "NC_FLOAT", dimensions = c("lon", "lat"))
att.put.nc(nc, variable = "data", name = "long_name", type = "NC_CHAR", value = "CRS definition")
att.put.nc(nc, variable = "data", name = "units", type = "NC_DOUBLE", value = "your_standard_cf_unit")
att.put.nc(nc, variable = "data", name = "standard_name", type = "NC_CHAR", value = "standardized_var_name") # Not needed for time coordinate variable
 
# Other recommended attributes
att.put.nc(nc, variable = "data", name = "grid_mapping", type = "NC_CHAR", value = "crs")
att.put.nc(nc, variable = "data", name = "_FillValue", type = "NC_FLOAT", value = -9999.9)

```



## Creating NetCDF files for EMODnet Biology

NetCDF files for the purposes required by EMODnet Biology will have slight modifications to what is presented to the CF compliance as we need to extend the guidance so that the file is also COARDS compliant. Those modifications include, but are not restricted to:

- Including extensive metadata in the Global attributes. VLIZ can provide you with the metadata record ID (in yellow) so you can incorporate it in your file.

```
// global attributes:
                :title = "Presence/Absence maps of phytoplankton in the Greater Baltic Sea" ;
                :summary = "The project aims to produce comprehensive data product of the occurence and absence of (phyto)plankton species. As a basis, data from EMODnet Biology are used. The selection of relevant datasets is optimized in order to find all planktonic species, and exclude all species that are not planktonic. The occurences from EMODnet Biology were complemenented with absence data assuming fixed species lists within each dataset and year. The products are presented as maps of the distribution of the 20 most common species of (phyto)plankton in the Greater Baltic Sea." ;
                :keywords = "Marine/Coastal, Baltic sea, Marine, Phytoplankton" ;
                :Conventions = "CF-1.8" ;
                :naming_authority = "https://emodnet.ec.europa.eu/en/biology" ;
                :history = "https://github.com/EMODnet/EMODnet-Biology-Phytoplankton-Greater-BalticSea" ;
                :source = "https://github.com/EMODnet/EMODnet-Biology-Phytoplankton-Greater-BalticSea" ;
                :license = "CC-BY" ;
                :standard_name_vocabulary = "CF Standard Name Table vNN" ;
                :date_created = "2022-08-10" ;
                :creator_name = "Figueroa Daniela" ;
                :creator_email = "daniela.figueroa@smhi.se" ;
                :creator_url = "www.smhi.se" ;
                :institution = "Swedish Meteorological and Hydrological Institute (SMHI)" ;
                :project = "EMODnet-Biology" ;
                :publisher_name = "EMODnet Biology Data Management Team" ;
                :publisher_email = "bio@emodnet.eu" ;
                :publisher_url = "https://emodnet.ec.europa.eu/en/biology" ;
                :geospatial_lat_min = 52.0285714285714 ;
                :geospatial_lat_max = 67.9714285714286 ;
                :geospatial_lon_min = 9.09333333333333 ;
                :geospatial_lon_max = 36.9066666666667 ;
                :sea_name = "Baltic Sea" ;
                :creator_institution = "Swedish Meteorological and Hydrological Institute (SMHI)" ;
                :publisher_institution = "Flanders Marine Institute (VLIZ)" ;
                :geospatial_lat_units = "degrees_north" ;
                :geospatial_lon_units = "degrees_east" ;
                :date_modified = "2021-03-08" ;
                :date_issued = "2020-12-24" ;
                :date_metadata_modified = "2021-03-08" ;
                :product_version = "1" ;
                :metadata_link = "https://marineinfo.org/imis?module=dataset&dasid=6618" ;
                :comment = "Uses attributes recommended by http://cfconventions.org" ;
                :citation = "Daniela Figueroa, Markus Lindh, Luuk van der Heijden, Willem Stolte & Lisa Sundqvist (2020). Presence/Absence maps of phytoplankton in the Greater Baltic Sea. Integrated data products created under the European Marine Observation Data Network (EMODnet) Biology project CINEA/EMFAF/2022/3.5.2/SI2.895681, funded by the by the European Union under Regulation (EU) No 508/2014 of the European Parliament and of the Council of 15 May 2014 on the European Maritime and Fisheries Fund" ;
                :acknowledgement = "European Marine Observation Data Network (EMODnet) Biology project CINEA/EMFAF/2022/3.5.2/SI2.895681, funded by the European Union under Regulation (EU) No 508/2014 of the European Parliament and of the Council of 15 May 2014 on the European Maritime and Fisheries Fund" ;
```



- Adding units to all variables, including those that do not have them

```
variables:
        double lon(lon) ;
                lon:units = "degrees_east" ;
                lon:standard_name = "longitude" ;
                lon:long_name = "Longitude" ;
        double lat(lat) ;
                lat:units = "degrees_north" ;
                lat:standard_name = "latitude" ;
                lat:long_name = "Latitude" ;
        int aphiaid(aphiaid) ;
                aphiaid:long_name = "Life Science Identifier - World Register of Marine Species" ;
                aphiaid:units = "level" ;
        char crs ;
                crs:long_name = "Coordinate Reference System" ;
                crs:geographic_crs_name = "WGS 84" ;
                crs:grid_mapping_name = "latitude_longitude" ;
                crs:reference_ellipsoid_name = "WGS 84" ;
                crs:prime_meridian_name = "Greenwich" ;
                crs:longitude_of_prime_meridian = 0. ;
                crs:semi_major_axis = 6378137. ;
                crs:semi_minor_axis = 6356752.31424518 ;
                crs:inverse_flattening = 298.257223563 ;
                crs:spatial_ref = "GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AXIS[\"Latitude\",NORTH],AXIS[\"Longitude\",EAST],AUTHORITY[\"EPSG\",\"4326\"]]" ;
                crs:GeoTransform = "-180 0.08333333333333333 0 90 0 -0.08333333333333333 " ;
        double probability_occurrence(aphiaid, lat, lon) ;
                probability_occurrence:_FillValue = -99999.9 ;
                probability_occurrence:long_name = "Probability of occurrence of biological entity" ;
                probability_occurrence:units = "1" ;
```

- Add aphiaid as a dimension (differs from CF guidance). The size will refer to the # of species in your file
- Add string as a dimension. The size will reflect the length needed for the lsid URL (80 should be sufficient)

```
dimensions:
        lon = 150 ;
        lat = 280 ;
        aphiaid = 331 ;
        string80 = 80 ;
```

- aphiaid should be a variable (differs from CF guidance). Also bear in mind that you cannot repeat taxa in your netCDF files and that taxa should be organised in a monotonic way (e.g. ever increasing or ever decreasing)

```
int aphiaid(aphiaid) ;
         aphiaid:long_name = "Life Science Identifier - World Register of Marine Species" ;
```

- taxon\_name should be a variable and defined as (same as in CF guidance)

```
char taxon_name(aphiaid, string80) ;
                taxon_name:standard_name = "biological_taxon_name" ;
                taxon_name:long_name = "Scientific name of the taxa" ;
```

- taxon\_lsid should be a variable and defined as (same as in CF guidance). Note that this is the full URL for each taxa, not just the ID number, e.g "urn:lsid:marinespecies.org:taxname:802"

```
char taxon_lsid(aphiaid, string80) ;
                taxon_lsid:standard_name = "biological_taxon_lsid" ;
                taxon_lsid:long_name = "Life Science Identifier - World Register of Marine Species" ;
```

- Your variable(s) will need to be dependent on aphiaid (as well as other dimensions, you can also add time, if it makes sense to your data) (differs from CF guidance), e.g

```
double probability_occurrence(aphiaid, lat, lon) ;
                probability_occurrence:_FillValue = -99999.9 ;
                probability_occurrence:long_name = "Probability of occurrence of biological entity" ;
```

These changes are needed as we can then ensure that by uploading the data to ERDDAP and adding the products in the map viewer as one layer per product, users will be able to filter for taxa by name, subset and download the data in various formats, including ascci and binary.
