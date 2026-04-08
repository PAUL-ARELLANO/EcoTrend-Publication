/**
 * ==========================================================================================
 * SCRIPT: MULTI-YEAR SCTI (PERIOD SELECTOR)
 * USE: Switch between 2022-2023 and 2023-2024
 * ==========================================================================================
 * SCRIPT: MULTI-YEAR ECOTREND FRAMEWORK (PERIOD SELECTOR)
 * Dr. Paul Arellano/School of Informatics, Computing and Cuber Systems, Northern Arizona University, Flagstaff, Arizona, USA
 * Reference: GitHub Repo: https://github.com/PAUL-ARELLANO/EcoTrend-Publication.git]
 * VERSION: 2026.03
 * * DESCRIPTION:
 * 
 * This script implements the EcoTrend geospatial framework to monitor forest structural 
 * degradation and physiological stress. It calculates the Standardized Composite Trend 
 * Index (SCTI) by fusing optical (NDVI) and radar (RVI) signals across a user-defined 
 * multi-year period (e.g., 2022-2023 or 2023-2024).
 * 
 * * METHODOLOGICAL WORKFLOW:
 * 
 * 1. HIERARCHICAL DATA REDUCTION:
 * - Sentinel-2 (Optical): Uses a 95th percentile reducer to capture peak annual 
 * photosynthetic vigor while filtering atmospheric haze and cloud shadows.
 * - Sentinel-1 (C-band SAR): Applies a 3m focal median spatial filter to suppress 
 * speckle noise, followed by a median temporal reducer to stabilize the structural 
 * baseline against moisture fluctuations.
 * 
 * * 2. INTER-ANNUAL TREND ANALYSIS:
 * - Performs a pixel-wise least-squares linear regression across the selected years.
 * - Generates slope coefficients (S) representing the trajectory of forest health.
 * 
 * * 3. Z-SCORE STANDARDIZATION & FUSION:
 * - Normalizes NDVI and RVI slopes independently using regional mean (μ) and 
 * standard deviation (σ).
 * - This standardization effectively isolates relative local anomalies from 
 * broad-scale inter-annual climatic fluctuations or sensor calibration shifts.
 * - Fuses the results into the SCTI: (Z_NDVI + Z_RVI) / 2.
 * 
 * * OUTPUT:
 * 
 * - A 10m resolution continuous metric where:
 * Negative Values (Red): Indicate forest decline or structural loss.
 * Positive Values (Green): Indicate vegetation recovery or reforestation.
 * Values near Zero (White): Indicate forest stability.
 * * ==========================================================================================
 */

// --- 1. CONTROL PANEL: SET YOUR TARGET PERIOD HERE ---
var startYear = 2023; // Set to 2022 for Period 1, or 2023 for Period 2
var endYear = 2024;   // Set to 2023 for Period 1, or 2024 for Period 2

// Standard JavaScript concatenation (Fixes the ee.String error)
var dateStart = startYear + '-01-01';
var dateEnd = endYear + '-12-31';

var years = ee.List.sequence(startYear, endYear);

// --- 2. SETUP: Study Area ---
var aoi_collection = ee.FeatureCollection('projects/paul-gee/assets/Priority_Areas_Merged');
var studyArea = aoi_collection.filter(ee.Filter.eq('OBJECTID', 334));
Map.centerObject(studyArea, 11);
Map.addLayer(studyArea, {color: 'FF0000'}, 'Study Area');

// ======================================================================
// 1. NDVI PROCESSING
// ======================================================================
var sentinel2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
  .filterBounds(studyArea)
  .filterDate(dateStart, dateEnd);

function computeNDVI(image) {
  var ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI');
  return image.addBands(ndvi);
}

var yearlyNDVI = years.map(function(y) {
  var year = ee.Number(y);
  var startDate = ee.Date.fromYMD(year, 1, 1);
  var composite = sentinel2.map(computeNDVI)
    .filterDate(startDate, ee.Date.fromYMD(year, 12, 31))
    .select('NDVI')
    .reduce(ee.Reducer.percentile([95]))
    .rename('NDVI');
  
  return composite.clip(studyArea)
    .set('year', year)
    .set('system:time_start', startDate.millis());
});

var ndviCol = ee.ImageCollection.fromImages(yearlyNDVI);

// ======================================================================
// 2. RVI PROCESSING
// ======================================================================
var s1_collection = ee.ImageCollection("COPERNICUS/S1_GRD")
  .filterBounds(studyArea)
  .filterDate(dateStart, dateEnd)
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'));

function addRVI(image) {
  var vv_lin = ee.Image(10).pow(image.select('VV').divide(10));
  var vh_lin = ee.Image(10).pow(image.select('VH').divide(10));
  var rvi = vh_lin.multiply(4).divide(vv_lin.add(vh_lin)).rename('RVI');
  return image.addBands(rvi);
}

var yearlyRVI = years.map(function(y) {
  var year = ee.Number(y);
  var startDate = ee.Date.fromYMD(year, 1, 1);
  var filtered = s1_collection.filterDate(startDate, ee.Date.fromYMD(year, 12, 31));
  
  var composite = ee.Image(ee.Algorithms.If(
    filtered.size().gt(0),
    filtered.map(function(img) { return img.focal_median(3, 'square', 'meters', 1); })
      .map(addRVI).select('RVI').reduce(ee.Reducer.median()).rename('RVI'),
    ee.Image.constant(0).rename('RVI').mask(ee.Image.constant(0))
  ));
  
  return composite.clip(studyArea)
    .set('year', year)
    .set('system:time_start', startDate.millis());
});

var rviCol = ee.ImageCollection.fromImages(yearlyRVI);

// ======================================================================
// 3. TRENDS & SCTI FUSION
// ======================================================================
var ndviTrend = ndviCol.map(function(img) {
  return ee.Image.constant(img.get('year')).rename('x').toFloat().addBands(img.select('NDVI').rename('y'));
}).reduce(ee.Reducer.linearFit()).select('scale').rename('NDVI_slope');

var rviTrend = rviCol.map(function(img) {
  return ee.Image.constant(img.get('year')).rename('x').toFloat().addBands(img.select('RVI').rename('y'));
}).reduce(ee.Reducer.linearFit()).select('scale').rename('RVI_slope');

var combinedTrend = ndviTrend.addBands(rviTrend);

var stats = combinedTrend.reduceRegion({
  reducer: ee.Reducer.mean().combine({reducer2: ee.Reducer.stdDev(), sharedInputs: true}),
  geometry: studyArea.bounds(),
  scale: 10,
  maxPixels: 1e13
});

var SCTI = combinedTrend.expression(
  '(((NDVI - MU_N) / SIG_N) + ((RVI - MU_R) / SIG_R)) / 2', {
    'NDVI': combinedTrend.select('NDVI_slope'),
    'MU_N': ee.Number(stats.get('NDVI_slope_mean')),
    'SIG_N': ee.Number(stats.get('NDVI_slope_stdDev')),
    'RVI': combinedTrend.select('RVI_slope'),
    'MU_R': ee.Number(stats.get('RVI_slope_mean')),
    'SIG_R': ee.Number(stats.get('RVI_slope_stdDev'))
}).rename('SCTI').reproject('EPSG:4326', null, 10).clip(studyArea);

// ======================================================================
// 4. VISUALIZATION & EXPORT
// ======================================================================
var viz = {min: -1.5, max: 1.5, palette: ['#d7191c', '#fdae61', 'white', '#a6d96a', '#1a9641']};
Map.addLayer(SCTI, viz, 'SCTI Trend ' + startYear + '-' + endYear);

Export.image.toDrive({
  image: SCTI.multiply(1000).toInt16().unmask(-9999),
  description: 'SCTI_Trend_' + startYear + '_' + endYear,
  fileNamePrefix: 'SCTI_' + startYear + '_' + endYear + '_OBJ334',
  folder: 'QAQC',
  region: studyArea.bounds(),
  scale: 10,
  fileFormat: 'GEO_TIFF',
  maxPixels: 1e13
});
