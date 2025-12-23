Alpine Treeline Ecotone Detection Using Google Earth Engine

This Google Earth Engine (GEE) script implements a complete, reproducible workflow for detecting and mapping the Alpine Treeline Ecotone (ATE) using multi-temporal Landsat surface reflectance data and topographic information in the Eastern Slopes of the Canadian Rockies. The code is designed for large-area, high-resolution analysis of treeline dynamics in mountainous environments and follows a physically and ecologically meaningful index-based approach.

The workflow begins with study-area definition and map initialization, followed by robust preprocessing of Landsat imagery (Landsat 5 and 8 Collection 2 Level-2 products). Cloud, cloud shadow, snow, and water pixels are masked using QA_PIXEL bitwise operations, and surface reflectance scaling is applied to ensure radiometric consistency across sensors and time. Vegetation dynamics are characterized using the Normalized Difference Vegetation Index (NDVI), calculated for each image and used to generate a greenest-pixel composite during the peak growing season (July–September).

To reduce noise and emphasize ecotonal patterns, NDVI and elevation layers are smoothed using circular focal median kernels. Three ecologically meaningful components are then derived:
C1 – spatial magnitude of NDVI gradients (abrupt vegetation transitions),
C2 – intermediate NDVI values modelled with a Gaussian function, and
C3 – spatial covariation between NDVI and elevation gradients, capturing topographic control on vegetation structure.

Each component is standardized using z-score normalization over the study area. These standardized components are combined through a logistic regression formulation to compute the Alpine Treeline Ecotone Index (ATEI), producing a continuous probability surface (0–1) representing the likelihood of treeline ecotone presence.

The script includes extensive visualization options for intermediate products (NDVI, gradients, components, and ATEI), as well as post-processing steps such as spatial smoothing and export to Google Drive as cloud-optimized GeoTIFFs. Overall, this code provides a scalable and transparent framework for long-term treeline monitoring, comparative watershed analysis, and climate–vegetation interaction studies in alpine regions.
If you need any help using this code and processing your ATE detection,  please reach out to me at Hooshyarkhah@uleth.ca
If you use this code, please cite the following paper: " Mapping Four Decades of Treeline Ecotone Migration: Remote Sensing of Alpine Ecotone Shifts on the Eastern Slopes of the Canadian Rocky Mountains".  
https://doi.org/10.3390/rs17244004


/*******************************************************
 * Alpine Treeline Ecotone (ATE) Detection using GEE
 *
 * This script detects Alpine Treeline Ecotones using
 * Landsat surface reflectance imagery and topographic
 * data. The workflow includes cloud masking, NDVI
 * computation, spatial smoothing, gradient analysis,
 * and derivation of the Alpine Treeline Ecotone Index
 * (ATEI) following a logistic model.
 *
 * Author: [Behnia Hooshyarkhah]
 * Platform: Google Earth Engine(GEE)
 * Spatial Resolution: 30 m
 * Temporal Focus: Peak growing season (Jul–Sep)
 *******************************************************/
// Center the interactive map on the study area (ROI) at zoom level 6
Map.centerObject(ROI, 6);

// Add the study area boundary to the map for reference (set to hidden by default)
Map.addLayer(ROI, false);

// ------------------------------------------------------------
// Function: add_ndvi_l8
// Purpose : Calculate and append NDVI to Landsat 8 imagery
// Input   : ee.Image (Landsat 8 surface reflectance image)
// Output  : ee.Image with an added NDVI band
// ------------------------------------------------------------
function add_ndvi_l8(image) {

  // Compute NDVI using Near-Infrared (Band 4) and Red (Band 3)
  var ndvi = image
    .normalizedDifference(['B4', 'B3'])
    .rename('NDVI');

  // Add the NDVI band to the original image
  return image.addBands(ndvi);
}
// ------------------------------------------------------------
// Function: clip_rmnp
// Purpose : Spatially clip imagery to the ESCR study boundary
// Input   : ee.Image (full-extent image)
// Output  : ee.Image clipped to the RMNP / ESCR boundary
// ------------------------------------------------------------
function clip_rmnp(image) {

  // Restrict all pixel processing to the study area geometry
  return image.clip(rmnp_boundary);
}

// ------------------------------------------------------------
// Function: mask_landsat8
// Purpose : Remove clouds, cloud shadows, snow, and water pixels
//           from Landsat 8 imagery using the QA_PIXEL band
// Input   : ee.Image (unmasked Landsat 8 image)
// Output  : ee.Image (quality-masked image)
// ------------------------------------------------------------
function mask_landsat8(image) {

  // Select the QA_PIXEL band containing pixel-level quality flags
  var qa_band = image.select('QA_PIXEL');
  
  // Extract water mask (bit 2): identifies water-covered pixels
  var water_bitmask = extract_qa_bits(qa_band, 2, 2, 'water_bitmask');

  // Extract cloud shadow mask (bit 3): identifies cloud-shadowed pixels
  var cloud_shadow_bitmask = extract_qa_bits(
    qa_band, 3, 3, 'cloud_shadow_bitmask'
  );

  // Extract snow/ice mask (bit 4): identifies snow-covered pixels
  var snow_bitmask = extract_qa_bits(qa_band, 4, 4, 'snow_bitmask');

  // Extract cloud mask (bit 5): identifies cloudy pixels
  var cloud_bitmask = extract_qa_bits(qa_band, 5, 5, 'cloud_bitmask');

  // Apply combined mask: retain only pixels free of all contaminants
  var image_masked = image.updateMask(
    water_bitmask.eq(0)              // Exclude water pixels
      .and(cloud_shadow_bitmask.eq(0)) // Exclude cloud shadows
      .and(snow_bitmask.eq(0))         // Exclude snow/ice
      .and(cloud_bitmask.eq(0))        // Exclude clouds
  );
  
  // Return the masked Landsat 8 image
  return image_masked;
}
// ------------------------------------------------------------
// Function: add_ndvi_l8
// Purpose : Compute NDVI and append it as a new band
// Note    : Uses band names B4 (NIR) and B3 (Red) from your renamed stack
// ------------------------------------------------------------
function add_ndvi_l8(image) {

  // Compute NDVI = (NIR - Red) / (NIR + Red), then name the output band "NDVI"
  var ndvi = image.normalizedDifference(['B4', 'B3']).rename('NDVI');

  // Add NDVI as an additional band to the input image
  return image.addBands(ndvi);
}


// ------------------------------------------------------------
// Create a Landsat 8 ImageCollection for the study area and time window
// Dataset : LANDSAT/LC08/C02/T1_L2 (Collection 2, Tier 1, Level-2 SR)
// الهدف   : Build a clean growing-season collection for composite/ATEI steps
// ------------------------------------------------------------
var rmnp_2021_collection = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")

  // Filter imagery to the analysis period (update as needed for your study)
  .filterDate('2021-03-01', '2023-12-31')

  // Restrict to peak growing season months (July–September) to reduce phenology noise
  .filter(ee.Filter.calendarRange(7, 9, 'month'))

  // Keep only scenes intersecting the study boundary (ROI)
  .filterBounds(ROI)

  // Sort scenes by metadata cloud cover (lowest cloud cover first)
  .sort("CLOUD_COVER")

  // Remove very cloudy scenes using metadata threshold (here: < 4%)
  .filterMetadata('CLOUD_COVER', 'less_than', 4)

  // Optional: apply pixel-level QA masking (cloud/shadow/snow/water) if desired
  // .map(mask_landsat8)

  // Convert SR bands to reflectance and standardize band names for later NDVI calls
  .map(function(img){

    // Select Landsat 8 SR bands and rename them to a consistent B1–B4 scheme
    // SR_B2=Blue, SR_B3=Green, SR_B4=Red, SR_B5=NIR
    var bands = img.select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5'],
                           ['B1',   'B2',   'B3',   'B4']);

    // Convert scaled integers to surface reflectance using Landsat C2 SR scale/offset
    var ref = bands.multiply(2.75e-05).add(-0.2);

    // Preserve original image metadata (e.g., time, cloud cover, IDs)
    return ref.copyProperties(img, img.propertyNames());
  })

  // Add NDVI band to every image in the collection
  .map(add_ndvi_l8)

  // Clip each image to the study boundary to reduce processing time and export size
  .map(clip_rmnp);


// Print the final filtered and processed collection to the Console for verification
print('rmnp_2021_collection', rmnp_2021_collection);

// ------------------------------------------------------------
// Create a multi-temporal composite from the filtered Landsat 8
// collection using the median reducer to suppress residual noise
// ------------------------------------------------------------
var Comp21 = rmnp_2021_collection.median();

// Optional: print composite for inspection
// print('Comp21', Comp21);


// ------------------------------------------------------------
// Evaluate average scene-level cloud cover across the collection
// (metadata-based, not pixel-level QA masking)
// ------------------------------------------------------------

// Extract CLOUD_COVER metadata values into a list
var cloudCoverList = rmnp_2021_collection.aggregate_array('CLOUD_COVER');

// Compute the mean cloud cover percentage for quality assessment
var cloudCoverMean = cloudCoverList.reduce(ee.Reducer.mean());

// Optional: display average cloud cover in the Console
// print('Average Cloud Cover Range (%):', cloudCoverMean);


// ------------------------------------------------------------
// Visualization parameters for Landsat 8 RGB composite
// Bands: B4 (NIR), B3 (Red), B2 (Green)
// ------------------------------------------------------------
var vis_params_landsat8_rgb = {
  bands: ['B4', 'B3', 'B2'],  // False-color composite highlighting vegetation
  min: 0,
  max: 0.25
};


// ------------------------------------------------------------
// Add the median composite to the map for visual inspection
// Layer is turned off by default to reduce map clutter
// ------------------------------------------------------------
Map.addLayer(
  Comp21,
  vis_params_landsat8_rgb,
  'Comp21 Median Composite (2021–2023)',
  false
);
// ------------------------------------------------------------
// Step 1) Build a "greenest-pixel" composite using maximum NDVI
// qualityMosaic('NDVI') selects, per pixel, the image with the highest NDVI
// ------------------------------------------------------------
var rmnp_2021_greenest = rmnp_2021_collection.qualityMosaic('NDVI');

// Print the selected greenest-pixel composite to verify bands/metadata
print('Greenest Pixel:', rmnp_2021_greenest);


// ------------------------------------------------------------
// Step 2) Define spatial smoothing kernel for NDVI and elevation
// Circular kernel radius is given in pixels (30 m pixels for Landsat)
// NOTE: Your comment mentions reducing radius; update radius value accordingly
// ------------------------------------------------------------
var ate_kernel = ee.Kernel.circle({
  radius: 10,          // Kernel radius in pixels (e.g., 10 px ≈ 300 m at 30 m)
  units: 'pixels'
});


// ------------------------------------------------------------
// Step 3) Smooth NDVI using a focal median filter to suppress speckle/noise
// ------------------------------------------------------------

// Extract NDVI from the greenest-pixel composite
var ndvi_greenest = rmnp_2021_greenest.select('NDVI');

// Apply median filter and clip to the study boundary
var ndvi_greenest_smoothed = ndvi_greenest
  .focal_median({kernel: ate_kernel, iterations: 1})  // Spatial smoothing
  .clip(rmnp_boundary);                                // Limit to study area


// ------------------------------------------------------------
// NDVI visualization parameters (0–1 typical vegetation range; -1–1 full NDVI)
// ------------------------------------------------------------
var vis_params_ndvi = {
  bands: ['NDVI'],          // Display NDVI band
  min: -1,                  // Lower NDVI bound
  max: 1,                   // Upper NDVI bound
  palette: ['ffffff', 'f7fcb9', 'addd8e', '31a354'] // Low → high vegetation
};


// Add NDVI layers for visual QA (both off by default)
Map.addLayer(ndvi_greenest, vis_params_ndvi,
  'ESCR - 2021 - Greenest Pixel - NDVI', false);

Map.addLayer(ndvi_greenest_smoothed, vis_params_ndvi,
  'ESCR - 2021 - Greenest Pixel - NDVI Smoothed', false);


// ------------------------------------------------------------
// Step 4) Prepare elevation layer (DEM) and apply the same smoothing approach
// ------------------------------------------------------------

// Assign the DEM input (must already exist in your script as "DEM")
var CDEM_ned = DEM;

// Clip DEM to study area to reduce processing
var rmnp_elevation = CDEM_ned.clip(rmnp_boundary);

// Smooth elevation to reduce local micro-topographic noise
var rmnp_elevation_smoothed = rmnp_elevation
  .focal_median({kernel: ate_kernel, iterations: 1})
  .clip(rmnp_boundary);


// ------------------------------------------------------------
// Elevation visualization parameters (adjust min/max to your region)
// ------------------------------------------------------------
var vis_params_elevation_rmnp = {
  min: 1000.0,
  max: 5000.0,
  palette: ['blue', 'green', 'yellow', 'orange', 'red', 'brown', 'white']
};


// ------------------------------------------------------------
// Gradient-direction visualization (0–360 degrees)
// Used later for NDVI and elevation gradient direction products
// ------------------------------------------------------------
var vis_params_gradient_direction = {
  min: 0,
  max: 360,
  palette: [
    '#666666', '#bf5b17', '#f0027f', '#386cb0',
    '#ffff99', '#fdc086', '#beaed4', '#7fc97f'
  ]
};


// Add DEM layers for visual QA (both off by default)
Map.addLayer(rmnp_elevation, vis_params_elevation_rmnp, 'ESCR - Elevation', false);
Map.addLayer(rmnp_elevation_smoothed, vis_params_elevation_rmnp, 'ESCR - Elevation Smoothed', false);


// ============================================================
// Step 5) Data Processing for ATEI Components
// ============================================================


// ------------------------------------------------------------
// C1: Abrupt spatial transitions in NDVI
// Compute gradient of the smoothed NDVI surface
// ------------------------------------------------------------
var ndvi_gradient = ndvi_greenest_smoothed.gradient();


// ------------------------------------------------------------
// Compute gradient magnitude: sqrt(dx^2 + dy^2)
// This represents the strength of spatial NDVI change (edge intensity)
// ------------------------------------------------------------
var ndvi_gradient_magnitude = ndvi_gradient
  .select('x').pow(2)                           // dx^2
  .add(ndvi_gradient.select('y').pow(2))        // + dy^2
  .sqrt()                                       // sqrt(dx^2 + dy^2)
  .select(['x'], ['ndvi_grad_mag']);            // rename output band


// ------------------------------------------------------------
// Compute gradient direction in degrees (0–360)
// atan2(dy, dx) gives angle; then convert to degrees and normalize to 0–360
// ------------------------------------------------------------
var ndvi_gradient_direction = ndvi_gradient
  .select('y').atan2(ndvi_gradient.select('x')) // angle in radians
  .multiply(180).divide(Math.PI)                // radians → degrees
  .add(360).mod(360)                            // normalize to 0–360
  .select(['y'], ['ndvi_grad_dir']);             // rename output band


// Add NDVI gradient products for visual QA (off by default)
Map.addLayer(ndvi_gradient, {}, 'NDVI Gradient', false);

Map.addLayer(
  ndvi_gradient_magnitude,
  {
    min: 0,
    max: 0.005956791228225612,
    palette: ['404040', 'bababa', 'ffffff', 'f4a582', 'ca0020']
  },
  'NDVI Gradient Magnitude', false
);

Map.addLayer(
  ndvi_gradient_direction,
  vis_params_gradient_direction,
  'NDVI Gradient Direction', false
);


// ------------------------------------------------------------
// Create the C1 component by renaming gradient magnitude to "c1"
// ------------------------------------------------------------
var c1 = ndvi_gradient_magnitude.select(['ndvi_grad_mag'], ['c1']);


// ------------------------------------------------------------
// C2: Intermediate NDVI values modeled as a Gaussian response
// Parameters b and c control the Gaussian center and spread (from your paper)
// ------------------------------------------------------------
var b = 0.42;   // Gaussian center (NDVI peak location)
var c = 0.06;   // Gaussian spread (controls width)


// Compute C2 using a Gaussian-like function of smoothed NDVI
var c2 = ee.Image(Math.E).clip(rmnp_boundary)
  .pow(
    ndvi_greenest_smoothed.subtract(b).pow(2)
      .divide(ee.Number(c).pow(2).multiply(2))
      .multiply(-1)
  )
  .select(['constant'], ['c2']);


// ------------------------------------------------------------
// C3: NDVI–elevation spatial covariation
// Start by computing the gradient of smoothed elevation
// ------------------------------------------------------------
var elevation_gradient = rmnp_elevation_smoothed.gradient();


// Compute elevation gradient magnitude: sqrt(dx^2 + dy^2)
var elevation_gradient_magnitude = elevation_gradient
  .select('x').pow(2)
  .add(elevation_gradient.select('y').pow(2))
  .sqrt()
  .select(['x'], ['elevation_grad_mag']);

// ------------------------------------------------------------
// Step 6) Compute elevation gradient direction (0–360 degrees)
// Method: atan2(dy, dx) → degrees → normalize to 0–360
// Output: "elevation_grad_dir" (direction of steepest increase)
// ------------------------------------------------------------
var elevation_gradient_direction = elevation_gradient
  .select('y')                                 // Use dy component of elevation gradient
  .atan2(elevation_gradient.select('x'))       // Compute angle (radians) using atan2(dy, dx)
  .multiply(180).divide(Math.PI)               // Convert radians to degrees
  .add(360).mod(360)                           // Normalize to 0–360 (0=N, 90=E, 180=S, 270=W)
  .select(['y'], ['elevation_grad_dir']);      // Rename band for clarity

// Optional: print to Console for debugging
// print('Elevation Gradient Direction:', elevation_gradient_direction);


// ------------------------------------------------------------
// Step 7) Compute angular difference between NDVI and elevation
// gradients (theta). This captures alignment/misalignment of
// vegetation change vs topographic change.
// ------------------------------------------------------------
var theta = ndvi_gradient_direction
  .subtract(elevation_gradient_direction)      // Difference in direction (degrees)
  .add(360).mod(360);                          // Normalize to 0–360 for interpretable angles
// Interpretation:
// - theta < 90 or theta > 270 → gradients are broadly aligned
// - 90 < theta < 270           → gradients are broadly opposite


// ------------------------------------------------------------
// Step 8) Define parameter for the C3 response curve (from paper)
// n controls how strongly C3 penalizes misalignment (sharper response for larger n)
// ------------------------------------------------------------
var n = 10;


// ------------------------------------------------------------
// Step 9) Compute C3: NDVI–elevation directional covariation term
// Formula: C3 = [ (1 - cos(theta))^n ] / (2^n)
// - When theta ≈ 0 (same direction): cos(theta)=1 → C3≈0
// - When theta ≈ 180 (opposite):     cos(theta)=-1 → C3≈1
// ------------------------------------------------------------
var c3 = theta
  .cos()                                       // cos(theta) (theta must be in degrees in GEE trigonometry)
  .multiply(-1).add(1)                         // (1 - cos(theta))
  .pow(n)                                      // raise to power n to increase contrast
  .divide(ee.Number(2).pow(n))                 // normalize by 2^n to keep range ~0–1
  .select(['ndvi_grad_dir'], ['c3']);          // Rename output band to "c3"

print('C3 Component:', c3);


// ------------------------------------------------------------
// Optional visualization: show theta (direction difference) on map
// Colors indicate aligned vs opposite gradient directions
// ------------------------------------------------------------
Map.addLayer(
  theta,
  {
    min: 0,
    max: 360,
    palette: [
      '0571b0', // 0–90   : mostly aligned
      'ca0020', // 90–180 : mostly opposite
      'ca0020', // 180–270: mostly opposite
      '0571b0'  // 270–360: mostly aligned
    ]
  },
  'Difference in Gradient Direction (theta)', false
);


// ------------------------------------------------------------
// Step 10) Standardize C1, C2, C3 (z-score) for comparability
// Standardization removes scale differences between components
// ------------------------------------------------------------
var c1_standardized = standardize_component(c1, rmnp_boundary, 'c1');
var c2_standardized = standardize_component(c2, rmnp_boundary, 'c2');
var c3_standardized = standardize_component(c3, rmnp_boundary, 'c3');

// Optional: print standardized layers
// print('C1 Standardized:', c1_standardized);
// print('C2 Standardized:', c2_standardized);
// print('C3 Standardized:', c3_standardized);


// ------------------------------------------------------------
// Optional visualization: elevation gradient products for QA/QC
// ------------------------------------------------------------
Map.addLayer(elevation_gradient, {}, 'Elevation Gradient', false);

Map.addLayer(
  elevation_gradient_magnitude,
  {
    min: 0,
    max: 3.137629831784505,
    palette: ['404040', 'bababa', 'ffffff', 'f4a582', 'ca0020']
  },
  'Elevation Gradient Magnitude', false
);

Map.addLayer(
  elevation_gradient_direction,
  vis_params_gradient_direction,
  'Elevation Gradient Direction', false
);


// ============================================================
// Step 11) Compute Alpine Treeline Ecotone Index (ATEI)
// Model form: ATEI = e^x / (e^x + 1)
// Where: x = b0 + b1*C1z + b2*C2z + b3*C3z  (z = standardized)
// Coefficients taken from the referenced paper
// ============================================================

// Intercept term of the logistic model
var b0 = -1.47;

// Coefficients for standardized components (from paper)
var b1 = 0.42;
var b2 = 0.58;
var b3 = 0.56;


// ------------------------------------------------------------
// Compute the linear predictor (x) of the logistic model
// ------------------------------------------------------------
var atei_sum = c1_standardized.multiply(b1)     // b1*C1z
  .add(c2_standardized.multiply(b2))            // + b2*C2z
  .add(c3_standardized.multiply(b3))            // + b3*C3z
  .add(b0)                                      // + intercept
  .select(['c1_standardized'], ['atei_sum_comp']); // Rename for clarity


// ------------------------------------------------------------
// Compute exp(x) term of the logistic equation
// ------------------------------------------------------------
var atei_exponent = ee.Image(Math.E)
  .clip(rmnp_boundary)
  .pow(atei_sum)
  .select(['constant'], ['atei_exp_comp']);     // Rename for clarity


// ------------------------------------------------------------
// Compute ATEI probability surface: exp(x) / (exp(x) + 1)
// Output range is approximately 0–1
// ------------------------------------------------------------
var atei = atei_exponent
  .divide(atei_exponent.add(1))
  .select(['atei_exp_comp'], ['atei']);

print('ATEI Image:', atei);


// ------------------------------------------------------------
// Optional: compute min/max of ATEI (useful for QA and reporting)
// ------------------------------------------------------------
var atei_min = atei.reduceRegion({
  reducer: ee.Reducer.min(),
  geometry: rmnp_boundary.geometry(),
  scale: 100,
  maxPixels: 1e12
});

var atei_max = atei.reduceRegion({
  reducer: ee.Reducer.max(),
  geometry: rmnp_boundary.geometry(),
  scale: 100,
  maxPixels: 1e12
});

// Optional: print min/max
// print('ATEI Min:', atei_min);
// print('ATEI Max:', atei_max);


// ------------------------------------------------------------
// Visualization: add ATEI to the map (off by default)
// ------------------------------------------------------------
Map.addLayer(
  atei,
  {
    min: 0,
    max: 1,
    palette: [
      'fff7f3', 'fde0dd', 'fcc5c0', 'fa9fb5',
      'f768a1', 'dd3497', 'ae017e', '7a0177', '49006a'
    ]
  },
  'Alpine Treeline Ecotone Index (ATEI)', false
);


// ============================================================
// Optional visualization of components and standardized layers
// ============================================================

// Add raw components (off by default)
Map.addLayer(c1, {min: 0, max: 0.005956791228225612}, 'C1 Component', false);
Map.addLayer(c2, {min: 0, max: 1}, 'C2 Component', false);
Map.addLayer(c3, {}, 'C3 Component', false);

// Add standardized components (off by default)
Map.addLayer(
  c1_standardized,
  {palette: ['0571b0', '92c5de', 'f7f7f7', 'f4a582', 'ca0020']},
  'C1 Standardized', false
);

Map.addLayer(
  c2_standardized,
  {palette: ['0571b0', '92c5de', 'f7f7f7', 'f4a582', 'ca0020']},
  'C2 Standardized', false
);

Map.addLayer(
  c3_standardized,
  {palette: ['fef0d9', 'fdcc8a', 'fc8d59', 'e34a33', 'b30000']},
  'C3 Standardized', false
);


// ============================================================
// Post-processing: smooth ATEI before export to reduce speckle
// ============================================================

// Define a smaller kernel for final ATEI smoothing (3-pixel radius)
var atei_kernel = ee.Kernel.circle({radius: 3, units: 'pixels'});

// Apply focal median smoothing and clip to study boundary
var atei_smoothed = atei
  .focal_median({kernel: atei_kernel, iterations: 1})
  .clip(rmnp_boundary);


// Visualize smoothed ATEI (off by default)
Map.addLayer(
  atei_smoothed,
  {
    min: 0,
    max: 1,
    palette: [
      'fff7f3', 'fde0dd', 'fcc5c0', 'fa9fb5',
      'f768a1', 'dd3497', 'ae017e', '7a0177', '49006a'
    ]
  },
  'ATEI Smoothed (3-pixel Circular Median)', false
);


// ============================================================
// Export: send smoothed ATEI to Google Drive as a Cloud-Optimized GeoTIFF
// ============================================================
Export.image.toDrive({
  image: atei_smoothed,                         // Export target raster
  description: 'ATEI_Smoothed_2021',            // Task name in the Tasks tab
  folder: 'Your_Google_Drive_Folder_Name',      // Optional Drive folder name
  region: rmnp_boundary,                        // Export region = study area boundary
  scale: 30,                                    // Landsat spatial resolution (meters)
  maxPixels: 1e12,                              // Large limit for mountainous regions
  formatOptions: { cloudOptimized: true }       // Create a Cloud Optimized GeoTIFF (COG)
});
