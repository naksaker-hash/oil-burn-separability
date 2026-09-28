/* =====================================================================
 * Türkiye outline for Figure 1.
 *
 * fig05_turkey_outline.geojson, used to draw Figure 5, is not on disk
 * any more. This re-exports it from the same GAUL source the main
 * script already uses, so Figure 1 and Figure 5 share one boundary.
 *
 * Paste into a blank Code Editor script, run, start the task, and put
 * the CSV next to make_fig01_site.py. The .geo column holds the
 * geometry; make_fig01_site.py reads it directly.
 * ===================================================================== */

var gaul = ee.FeatureCollection('FAO/GAUL_SIMPLIFIED_500m/2015/level0');
var tr = gaul.filter(ee.Filter.eq('ADM0_NAME', 'Turkey'));

print('features (expect 1):', tr.size());

Export.table.toDrive({
  collection: tr,
  description: 'fig01_turkey_outline',
  folder: 'oil_burn_separability',
  fileFormat: 'GeoJSON'
});

Map.addLayer(tr, {color: 'grey'}, 'Türkiye');
Map.centerObject(tr, 5);
