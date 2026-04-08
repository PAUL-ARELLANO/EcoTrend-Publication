/**
 * =========================================================================================
 * ECOTREND FRAMEWORK: INTRA-ANNUAL FOREST MONITORING SCRIPT
 * =========================================================================================
 * * DESCRIPTION:
 * This script implements the Standardized Composite Trend Index (SCTI) as described in
 * [Your Paper Title/Link]. It integrates multi-sensor data to quantify forest structural 
 * degradation and recovery at a 10-meter spatial resolution.
 * * METHODOLOGY:
 * 1. DATA HARMONIZATION: 
 * - Sentinel-2 (Optical): Computes NDVI using 95th percentile monthly composites to 
 * capture peak physiological greenness while mitigating cloud interference.
 * - Sentinel-1 (SAR): Computes Radar Vegetation Index (RVI) from VV/VH backscatter 
 * using monthly median reducers to capture structural canopy density.
 * * 2. TREND ANALYSIS:
 * - Performs a pixel-wise Ordinary Least Squares (OLS) linear regression across the 
 * 12-month time series. The resulting slope represents the intra-annual trajectory 
 * of forest health.
 * * 3. SENSORY FUSION (SCTI):
 * - Standardizes NDVI and RVI slopes using Z-score normalization (Mean = 0, SD = 1).
 * - Computes the SCTI as the mean of the standardized slopes, ensuring equal weight
 * between optical (photosynthetic) and radar (structural) signals.
 * * 4. OUTPUT:
 * - A 10m GeoTIFF where negative values indicate degradation/biomass loss and 
 * positive values indicate vegetation gain/recovery.
 * * CONTACT & CITATION:
 * Dr. Paul Arellano/School of Informatics, Computing and Cuber Systems, Northern Arizona University, Flagstaff, Arizona, USA
 * Reference: GitHub Repo: https://github.com/PAUL-ARELLANO/EcoTrend-Publication.git]
 * =========================================================================================
 */

// --- 1. SETUP: Study Area and Timeframe ---
var aoi_collection = ee.FeatureCollection('projects/paul-gee/assets/Priority_Areas_Merged');
var studyArea = aoi_collection.filter(ee.Filter.eq('OBJECTID', 332)); // Set your ID here
Map.centerObject(studyArea, 12);
Map.addLayer(studyArea, {color: 'FF0000'}, 'Study Area');

var targetYear = 2022; // Set your Year here
var months = ee.List.sequence(1, 12);

// --- 2. HELPER FUNCTIONS ---

// Function to compute NDVI
function computeNDVI(image) {
  var ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI');
  return image.addBands(ndvi);
}

// Function to compute RVI (Radar Vegetation Index)
function addRVI(image) {
  var vv_lin = ee.Image(10).pow(image.select('VV').divide(10));
  var vh_lin = ee.Image(10).pow(image.select('VH').divide(10));
  var rvi = vh_lin.multiply(4).divide(vv_lin.add(vh_lin)).rename('RVI');
  return image.addBands(rvi);
}

// --- 3. MONTHLY NDVI PROCESSING ---

var s2_collection = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
  .filterBounds(studyArea)
  .filterDate(ee.Date.fromYMD(targetYear, 1, 1), ee.Date.fromYMD(targetYear, 12, 31));

var monthlyNDVI = months.map(function(m) {
  var startDate = ee.Date.fromYMD(targetYear, m, 1);
  var endDate = startDate.advance(1, 'month');
  
  var filtered = s2_collection.filterDate(startDate, endDate)
    .map(computeNDVI)
    .select('NDVI');
  
  // ERROR HANDLING: If month is empty, return a masked blank image to avoid crash
  var composite = ee.Image(ee.Algorithms.If(
    filtered.size().gt(0),
    filtered.reduce(ee.Reducer.percentile([95])).rename('NDVI'),
    ee.Image.constant(0).rename('NDVI').mask(ee.Image.constant(0))
  ));
    
  return composite.clip(studyArea)
    .set('month', m)
    .set('system:time_start', startDate.millis());
});

var ndviCol = ee.ImageCollection.fromImages(monthlyNDVI);

// --- 4. MONTHLY RVI PROCESSING ---

var s1_collection = ee.ImageCollection("COPERNICUS/S1_GRD")
  .filterBounds(studyArea)
  .filterDate(ee.Date.fromYMD(targetYear, 1, 1), ee.Date.fromYMD(targetYear, 12, 31))
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  // Consider removing the orbit filter if data gaps persist
  .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING')); 

var monthlyRVI = months.map(function(m) {
  var startDate = ee.Date.fromYMD(targetYear, m, 1);
  var endDate = startDate.advance(1, 'month');
  
  var filtered = s1_collection.filterDate(startDate, endDate)
    .map(function(img) { return img.focal_median(3, 'square', 'meters', 1); }) // Speckle filter
    .map(addRVI)
    .select('RVI');
  
  // ERROR HANDLING: If month is empty, return a masked blank image to avoid crash
  var composite = ee.Image(ee.Algorithms.If(
    filtered.size().gt(0),
    filtered.reduce(ee.Reducer.median()).rename('RVI'),
    ee.Image.constant(0).rename('RVI').mask(ee.Image.constant(0))
  ));
    
  return composite.clip(studyArea)
    .set('month', m)
    .set('system:time_start', startDate.millis());
});

var rviCol = ee.ImageCollection.fromImages(monthlyRVI);

// --- 5. COMPUTE TRENDS (Intrayear Slopes) ---

// NDVI monthly slope
var ndviTrend = ndviCol.map(function(image) {
  var xBand = ee.Image.constant(image.get('month')).rename('x').toFloat();
  var yBand = image.select('NDVI').rename('y');
  return xBand.addBands(yBand);
}).reduce(ee.Reducer.linearFit()).rename('NDVI_scale', 'NDVI_offset');

// RVI monthly slope
var rviTrend = rviCol.map(function(image) {
  var xBand = ee.Image.constant(image.get('month')).rename('x').toFloat();
  var yBand = image.select('RVI').rename('y');
  return xBand.addBands(yBand);
}).reduce(ee.Reducer.linearFit()).rename('RVI_scale', 'RVI_offset');

var combinedTrend = ndviTrend.addBands(rviTrend);

// --- 6. FUSE INTO STANDARDIZED COMPOSITE TREND INDEX (SCTI) ---

// Statistics for Z-score normalization
var trendStats = combinedTrend.reduceRegion({
  reducer: ee.Reducer.mean().combine({reducer2: ee.Reducer.stdDev(), sharedInputs: true}),
  geometry: studyArea.bounds(),
  scale: 10,
  maxPixels: 1e13
});

var ndvi_mu = ee.Number(trendStats.get('NDVI_scale_mean'));
var ndvi_sigma = ee.Number(trendStats.get('NDVI_scale_stdDev'));
var rvi_mu = ee.Number(trendStats.get('RVI_scale_mean'));
var rvi_sigma = ee.Number(trendStats.get('RVI_scale_stdDev'));

// Standardize and Average
var SCTI = combinedTrend.expression(
    '(((NDVI_SLOPE - NDVI_MU) / NDVI_SIGMA) + ((RVI_SLOPE - RVI_MU) / RVI_SIGMA)) / 2',
    {
        'NDVI_SLOPE': combinedTrend.select('NDVI_scale'),
        'RVI_SLOPE': combinedTrend.select('RVI_scale'),
        'NDVI_MU': ndvi_mu,
        'NDVI_SIGMA': ndvi_sigma,
        'RVI_MU': rvi_mu,
        'RVI_SIGMA': rvi_sigma
    }
).rename('SCTI');

// Force 10m resolution and clip to precise AOI
SCTI = SCTI.reproject({
    crs: 'EPSG:4326', 
    scale: 10
}).clip(studyArea);

// --- 7. VISUALIZATION AND EXPORT ---

var sctiVizParams = {
    min: -1.0, 
    max: 1.0, 
    palette: ['#d7191c', '#fdae61', 'white', '#a6d96a', '#1a9641'] 
};

Map.addLayer(SCTI, sctiVizParams, 'Intrayear SCTI Trend @ 10m');

// Export to Google Drive
Export.image.toDrive({
    image: SCTI.multiply(1000).toInt16().unmask(-9999),
    description: 'Intrayear_SCTI_' + targetYear + '_10m',
    fileNamePrefix: 'SCTI_Intrayear_' + targetYear + '_OBJ' + 334,
    folder: 'QAQC',
    region: studyArea.bounds(),
    scale: 10,
    fileFormat: 'GEO_TIFF',
    formatOptions: { cloudOptimized: true, noData: -9999 },
    maxPixels: 1e13
});

print('Processing for OBJECTID 334 complete.');
