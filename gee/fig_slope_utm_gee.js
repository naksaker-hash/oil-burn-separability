/* =====================================================================
 * Recompute slope in a conformal local projection.
 *
 * The analysis derived slope from the GLO-30 DEM reprojected to
 * EPSG:3857. Web Mercator carries a scale factor of 1/cos(phi), so at
 * 37.35 N a nominal 30 m cell spans 23.9 m on the ground and slope is
 * understated by about 26%. The bias is identical for the spill site,
 * every candidate burn event and every sampled pixel, so it does not
 * affect any comparison between classes -- but the absolute values
 * printed in the methods deserve a projection that does not distort
 * distance.
 *
 * This prints both, so the manuscript can quote a measured figure
 * instead of the analytic correction.
 * ===================================================================== */

var OIL_POINT = ee.Geometry.Point([37.130123, 37.350183]);
var demCol = ee.ImageCollection('COPERNICUS/DEM/GLO30').select('DEM');

function meanSlope(crs, label) {
  var dem = demCol.mosaic().setDefaultProjection(crs, null, 30).rename('DEM');
  var slope = ee.Terrain.products(dem).select('slope');
  print(label, slope.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: OIL_POINT.buffer(200),
    scale: 30,
    maxPixels: 1e9
  }).get('slope'));
}

// As used in the analysis.
meanSlope('EPSG:3857', 'mean slope, EPSG:3857 (as analysed):');

// UTM zone 37N, conformal and true to scale over this site.
meanSlope('EPSG:32637', 'mean slope, EPSG:32637 (UTM 37N):');
