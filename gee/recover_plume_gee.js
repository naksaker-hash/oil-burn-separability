/* =====================================================================
 * PLUME RECOVERY  —  append to narli_burn_confusion_gee.js
 *
 * The hand-digitised plume polygon was never saved: OIL_PLUME_GEOM and
 * OIL_PLUME_ASSET are both null in the archived script, and the export
 * dropped pixel geometries, so neither the label nor the locations of
 * the 218 sampled oil pixels survive anywhere.
 *
 * They can be recovered from the spectra. The oil rows in
 * oil_burn_samples.csv hold reflectance for ten bands, and that
 * ten-band tuple is unique across all 2311 pixels sampled at the Narlı
 * site — no two share one. So exporting every pixel in the 3 x 3 km
 * window with its coordinates, then matching on the ten values,
 * relocates each sampled pixel exactly.
 *
 * The window is 150 x 150 pixels at 20 m, so this is a small export.
 *
 * TO RUN
 *   1. Paste this whole block at the END of narli_burn_confusion_gee.js.
 *   2. Set RUN_STAGE = 'RECOVER' at line 41.
 *   3. Run, then start the task. It writes one CSV to Drive.
 *   4. Download it next to oil_burn_samples.csv and run:
 *        python3 recover_plume.py oil_window_pixels.csv
 *
 * Nothing above is modified and nothing is overwritten. The stack is
 * rebuilt by the same featureStack() call the original sampling used,
 * over the same window and the same event date, so the reflectance
 * values are the ones already in the archive.
 * ===================================================================== */

if (RUN_STAGE === 'RECOVER') {
  print('--- PLUME RECOVERY  window pixel export ---');

  var recStack = featureStack(OIL_WINDOW, OIL_EVENT_DATE);

  print('pre_date  (must match the archive):', recStack.get('pre_date'));
  print('post_date (must match the archive):', recStack.get('post_date'));

  // Pixel centres come through as ordinary properties, which keeps the
  // exported table flat and avoids a geometry column.
  var recWithXY = recStack.select(S2_BANDS).addBands(ee.Image.pixelLonLat());

  // No numPixels and no factor: sample() then returns every pixel in
  // the region rather than a random subset.
  var recPixels = recWithXY.sample({
    region: OIL_WINDOW,
    scale: SAMPLE_SCALE,
    dropNulls: true,
    geometries: false
  });

  print('pixels in window (expect roughly 22500):', recPixels.size());

  Export.table.toDrive({
    collection: recPixels,
    description: 'oil_window_pixels',
    folder: EXPORT_FOLDER,
    fileFormat: 'CSV',
    selectors: ['longitude', 'latitude'].concat(S2_BANDS)
  });

  Map.centerObject(OIL_POINT, 14);
  Map.addLayer(OIL_WINDOW, {color: 'yellow'}, 'Recovery window');
}
