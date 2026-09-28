/* =====================================================================
 * PLUME RECOVERY — minimal standalone script
 *
 * Paste this whole file into a blank Code Editor script and run it.
 * Nothing else is needed: it does not import the burn inventory, the
 * terrain fingerprint or the index library, which is what exhausted the
 * interactive memory on the full script. Those parts compute 6820
 * perimeters and 5450 candidate events before reaching the export, and
 * none of it is relevant here.
 *
 * What it needs to reproduce is narrow. The raw B2..B12 columns in
 * oil_burn_samples.csv come from featureStack's `single` term, which is
 * the POST scene — so only pickScene(..., false) has to be reproduced,
 * and it is copied verbatim from narli_burn_confusion_gee.js along with
 * every constant it reads.
 *
 * Expected console output:
 *   post scene date  2018-06-20      <- must match the archive
 *   pixels in window ~22380
 *
 * Then start the task, download oil_window_pixels.csv next to
 * oil_burn_samples.csv, and run:
 *   python3 recover_plume.py
 * ===================================================================== */

// --- copied verbatim from narli_burn_confusion_gee.js ---------------
var OIL_LON = 37.130123;
var OIL_LAT = 37.350183;
var OIL_POINT = ee.Geometry.Point([OIL_LON, OIL_LAT]);
var OIL_WINDOW = OIL_POINT.buffer(1500).bounds();

var OIL_EVENT_DATE = '2018-06-19';   // rupture

var CS_PLUS_BAND   = 'cs_cdf';
var CS_PLUS_THRESH = 0.60;
var CLEAR_MIN      = 0.90;
var PAIR_WINDOW_D  = 45;

var SAMPLE_SCALE   = 20;
var EXPORT_FOLDER  = 'oil_burn_separability';

var S2  = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED');
var CSP = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');

var S2_BANDS = ['B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12'];

function s2Clear(geom, start, end) {
  return S2.filterBounds(geom)
           .filterDate(start, end)
           .linkCollection(CSP, [CS_PLUS_BAND])
           .map(function (img) {
             var clear = img.select(CS_PLUS_BAND).gte(CS_PLUS_THRESH);
             return img.select(S2_BANDS)
                       .divide(10000)
                       .updateMask(clear)
                       .copyProperties(img, ['system:time_start']);
           });
}

function pickScene(geom, dateStr, before) {
  var d = ee.Date(dateStr);
  var start = before ? d.advance(-PAIR_WINDOW_D, 'day') : d.advance(1, 'day');
  var end   = before ? d.advance(-1, 'day') : d.advance(PAIR_WINDOW_D, 'day');

  var col = s2Clear(geom, start, end).map(function (img) {
    var frac = img.select('B8').mask().reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: geom,
      scale: 60,
      maxPixels: 1e8
    }).get('B8');
    return img.set('clear_frac', ee.Number(frac));
  }).filter(ee.Filter.gte('clear_frac', CLEAR_MIN));

  var sorted = col.sort('system:time_start', !before);
  return ee.Image(sorted.first());
}

// --- recovery export ------------------------------------------------

var post = pickScene(OIL_WINDOW, OIL_EVENT_DATE, false);

print('post scene date (must read 2018-06-20):',
      post.date().format('YYYY-MM-dd'));

var withXY = post.select(S2_BANDS).addBands(ee.Image.pixelLonLat());

// No numPixels and no factor, so sample() returns every pixel in the
// region rather than a random subset.
var pixels = withXY.sample({
  region: OIL_WINDOW,
  scale: SAMPLE_SCALE,
  dropNulls: true,
  geometries: false
});

print('pixels in window (expect about 22380):', pixels.size());

Export.table.toDrive({
  collection: pixels,
  description: 'oil_window_pixels',
  folder: EXPORT_FOLDER,
  fileFormat: 'CSV',
  selectors: ['longitude', 'latitude'].concat(S2_BANDS)
});

Map.centerObject(OIL_POINT, 14);
Map.addLayer(post.select(['B12','B11','B8A']), {min: 0, max: 0.4},
             'post SWIR composite');
Map.addLayer(OIL_WINDOW, {color: 'yellow'}, 'Recovery window', false);
