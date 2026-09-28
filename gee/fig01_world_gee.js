/* =====================================================================
 * World country outlines for Figure 1, panel (c).
 *
 * Same source as the Türkiye outline and the burn inventory, so all
 * three panels of Figure 1 and Figure 5 share one set of boundaries.
 * Using a different world dataset here would put two coastlines in one
 * figure.
 *
 * Geometries are simplified to 20 km before export. At the size this
 * panel is printed that is far below a pixel, and it keeps the file
 * small enough to handle comfortably.
 *
 * Paste into a blank Code Editor script, run, start the task, and put
 * the GeoJSON next to make_fig01_site.py.
 * ===================================================================== */

var SIMPLIFY_M = 20000;

var world = ee.FeatureCollection('FAO/GAUL_SIMPLIFIED_500m/2015/level0');

print('countries:', world.size());

var simplified = world.map(function (f) {
  return ee.Feature(f.geometry().simplify(SIMPLIFY_M))
           .set('ADM0_NAME', f.get('ADM0_NAME'));
});

Export.table.toDrive({
  collection: simplified,
  description: 'fig01_world_outline',
  folder: 'oil_burn_separability',
  fileFormat: 'GeoJSON'
});

Map.addLayer(simplified, {color: 'grey'}, 'world');
