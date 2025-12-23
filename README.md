Alpine Treeline Ecotone Detection Using Google Earth Engine

This Google Earth Engine (GEE) script implements a complete, reproducible workflow for detecting and mapping the Alpine Treeline Ecotone (ATE) using multi-temporal Landsat surface reflectance data and topographic information in the Eastern Slopes of the Canadian Rockies. The code is designed for large-area, high-resolution analysis of treeline dynamics in mountainous environments and follows a physically and ecologically meaningful index-based approach.

The workflow begins with study-area definition and map initialization, followed by robust preprocessing of Landsat imagery (Landsat 5 and 8 Collection 2 Level-2 products). Cloud, cloud shadow, snow, and water pixels are masked using QA_PIXEL bitwise operations, and surface reflectance scaling is applied to ensure radiometric consistency across sensors and time. Vegetation dynamics are characterized using the Normalized Difference Vegetation Index (NDVI), calculated for each image and used to generate a greenest-pixel composite during the peak growing season (July–September).

To reduce noise and emphasize ecotonal patterns, NDVI and elevation layers are smoothed using circular focal median kernels. Three ecologically meaningful components are then derived:
C1 – spatial magnitude of NDVI gradients (abrupt vegetation transitions),
C2 – intermediate NDVI values modeled with a Gaussian function, and
C3 – spatial covariation between NDVI and elevation gradients, capturing topographic control on vegetation structure.

Each component is standardized using z-score normalization over the study area. These standardized components are combined through a logistic regression formulation to compute the Alpine Treeline Ecotone Index (ATEI), producing a continuous probability surface (0–1) representing the likelihood of treeline ecotone presence.

The script includes extensive visualization options for intermediate products (NDVI, gradients, components, and ATEI), as well as post-processing steps such as spatial smoothing and export to Google Drive as cloud-optimized GeoTIFFs. Overall, this code provides a scalable and transparent framework for long-term treeline monitoring, comparative watershed analysis, and climate–vegetation interaction studies in alpine regions.
If you need any help for using this code and processing your ATE detection,  please reach out to me at Hooshyarkhah@uleth.ca
Please reference this code in this paper " Mapping Four Decades of Treeline Ecotone Migration: Remote Sensing of Alpine Ecotone Shifts on the Eastern Slopes of the Canadian Rocky Mountains"  
https://doi.org/10.3390/rs17244004



