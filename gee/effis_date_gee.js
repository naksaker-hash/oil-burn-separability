/* Recover the EFFIS ingest date for Sect. 2.2.
 *
 * The manuscript needs the date the EFFIS burnt area database was
 * downloaded. The upload date of the Earth Engine asset it was ingested
 * into is a defensible proxy and is recorded by Earth Engine itself.
 *
 * Paste into a blank Code Editor script and run. Read startTime /
 * updateTime from the printed object. */

print(ee.data.getAsset('projects/propane-forge-479411-k2/assets/effis_layer'));
