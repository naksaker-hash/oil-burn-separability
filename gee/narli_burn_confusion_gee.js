/**
 * =====================================================================
 * OIL vs BURN SEPARABILITY IN SENTINEL-2
 * Burn-scar inventory, substrate matching, and matched sample export
 * Study region: Kahramanmaras / Adana / Osmaniye, Turkiye
 * Oil reference: Kirkuk-Ceyhan pipeline rupture, Narli, 19 June 2018
 * =====================================================================
 *
 * PURPOSE
 *   Low et al. (2021, J. Environ. Manag. 298:113424) reported that
 *   vegetation fires confused their "burnt area" and "slightly oil
 *   polluted" classes, with class F1 falling to 40.0 % at one site.
 *   Burn indices (NBR, BAI) are simultaneously the strongest oil
 *   discriminators on bare soil. This script assembles the data needed
 *   to test whether those indices separate oil from burn scars, or
 *   only oil from clean soil.
 *
 * STAGE A (PARTS 1-5)  Inventory + substrate matching. Cheap. Run first.
 *   Answers: how many burn scars exist on Narli-like substrate, how
 *   large are they, and how many Sentinel-2 pixels do they carry?
 *   Prints a candidate event table and diagnostic charts.
 *
 * STAGE B (PARTS 6-9)  Index computation and sample export.
 *   Run only after Stage A confirms sufficient matched burns, and after
 *   pasting the selected event IDs into SELECTED_EVENTS below.
 *
 * Set RUN_STAGE = 'A' or 'B'.
 *
 * NOTE ON PROVENANCE
 *   The index library in PART 6 is applied identically to the oil site
 *   and to every burn event. Feature comparability across classes is
 *   the whole point; do not add a per-class preprocessing step.
 * =====================================================================
 */

// Stage A complete on EFFIS perimeters: 58 matched events across 23
// provinces, 7301 ha. Use this setting ONLY after the
// burn_inventory_asset export has finished, or Stage B will sample the
// previous MODIS inventory still sitting at that asset path.
// 'A' inventory, 'B' sample export, 'OIL' plume digitising view (PART 10).
var RUN_STAGE = 'OIL';

/* =====================================================================
 * PART 1  CONFIGURATION
 * ===================================================================== */

// --- Oil reference site (Narli) ------------------------------------
// Coordinates as corrected: latitude 37.350183 N, longitude 37.130123 E.
// Note the order: ee.Geometry.Point takes [lon, lat].
var OIL_LON = 37.130123;
var OIL_LAT = 37.350183;
var OIL_POINT = ee.Geometry.Point([OIL_LON, OIL_LAT]);

// 3 x 3 km analysis window, matching the earlier study's footprint.
var OIL_WINDOW = OIL_POINT.buffer(1500).bounds();

var OIL_DATE_PRE  = '2018-06-15';   // last clear pre-event scene
var OIL_DATE_POST = '2018-06-20';   // first clear post-event scene
var OIL_EVENT_DATE = '2018-06-19';  // rupture

// Plume mask from the earlier analysis. Export it from that project as
// an asset and paste the path here. Leave null to skip the oil rows in
// Stage B and export burn rows only.
// Oil plume label. Two ways to supply it, in order of preference:
//
// 1. OIL_PLUME_GEOM - a polygon digitised by hand in the Code Editor.
//    Set RUN_STAGE = 'OIL' to display the pre and post composites at
//    the site, draw around the plume with the polygon tool, rename the
//    drawn layer to plumeGeom, then set OIL_PLUME_GEOM = plumeGeom.
//    This is the preferred route: the burn labels are human-delineated
//    perimeters from EFFIS, so a hand-digitised plume keeps both
//    classes defined the same way. An algorithmic plume mask derived
//    from a burn index would reintroduce the circularity the whole
//    design is built to avoid.
//
// 2. OIL_PLUME_ASSET - a binary raster mask exported from elsewhere.
//
// If both are set, the geometry wins.
var OIL_PLUME_GEOM = null;
var OIL_PLUME_ASSET = null;

// --- Study region --------------------------------------------------
// REGION_SCOPE controls how wide the search goes.
//   'provinces' - Narli's own province plus Adana and Osmaniye.
//                 160 EFFIS perimeters, of which 18 survived the area
//                 and season filters. Too few to support a strict
//                 substrate envelope.
//   'country'   - every GAUL level-1 unit in the countries listed in
//                 REGION_COUNTRIES. Hundreds of perimeters, which is
//                 what makes strict matching affordable again.
//   'bbox'      - REGION_BBOX, for ad hoc tests.
//
// Widening the region is not a loosening of the design. The substrate
// envelope is what controls comparability, and a larger candidate pool
// lets that envelope be tightened rather than relaxed. Consider
// restoring SLOPE_MAX_DEG to 5, MATCH_SOIL_TEXTURE to true and the
// NDVI percentiles to p25-p75 once the matched count allows it, and
// report which envelope produced the final sample.
var REGION_SCOPE = 'country';

// EFFIS covers Europe, the Middle East and North Africa. Adding
// neighbours widens the pool but weakens any claim that the burns are
// climatically comparable to the spill site, so start with Turkiye
// alone and only extend if the matched count is still too low.
var REGION_COUNTRIES = ['Turkey'];

var PROVINCE_NAMES = ['Kahramanmaras', 'Adana', 'Osmaniye'];

// Fallback bounding box, used for REGION_SCOPE 'bbox' and if the
// province match returns nothing.
var REGION_BBOX = ee.Geometry.Rectangle([34.8, 36.5, 37.8, 38.4]);

// Cap on events entering the expensive annotation phase. The pre-fire
// NDVI composite is by far the costliest step, so terrain and soil are
// evaluated first and only the survivors, up to this many, get an NDVI
// composite. Raise it if the console copes; lower it if Stage A times
// out. Selection is by a seeded random column, so it is reproducible.
var CANDIDATE_LIMIT = 400;
var CANDIDATE_SEED = 17;

// Cap on MATCHED events written to the inventory and sampled in
// Stage B. The Turkiye-wide run matched 223 events, which at
// PIX_PER_EVENT would be roughly 67,000 burn rows and a very long
// Stage B. Capping keeps the sample large enough for stable estimates
// while staying tractable, and the excluded events remain in the
// exported inventory CSV so the selection is transparent.
// Set to 0 for no cap.
var MATCHED_EVENT_LIMIT = 80;

// Stage B row-count printing. Off by default: at 58 events it exceeds
// the interactive memory limit, because each print forces every
// event's scene pair and feature stack to evaluate.
var PRINT_STAGE_B_DIAGNOSTICS = false;

// Per-pixel substrate masking. OFF by default, for two reasons found
// on the first EFFIS run.
//
// 1. Asymmetry. The mask was applied when sampling burn and control
//    pixels but never in sampleOilSite, so burn pixels were pre-
//    selected to a narrow pre-fire NDVI band and oil pixels were not.
//    That biases the very comparison this study makes. Worse, pre-event
//    NDVI is one of the features under test, so selecting pixels on it
//    is the same circularity trap as thresholding on dNBR.
// 2. Attrition. Events pass on their MEAN substrate, then lose most of
//    their pixels to the per-pixel test: median burn rows per event
//    fell to 13, and only 20 of 57 events had 30 or more rows in both
//    classes. The sparse events were not smaller fires (median 47 ha
//    against 41 ha for the well-populated ones), so this was the mask,
//    not the geometry.
//
// Substrate is now controlled at the event level and recorded per
// pixel, so it can be handled as a covariate in the analysis, applied
// identically to every class.
var PIXEL_SUBSTRATE_MASK = false;

// --- Temporal extent -----------------------------------------------
// Sentinel-2 L2A over Turkiye is usable from 2017 onward.
var YEAR_START = 2017;
var YEAR_END   = 2026;

// --- Substrate matching tolerances ---------------------------------
// A burn scar in dense Mediterranean pine on a 20-degree slope is not a
// legitimate comparison for oil on a near-flat basin floor. These
// tolerances define the envelope. Widen them only with justification,
// and report whatever you use.
//
// Elevation matching is OFF by default. On the first full run it
// rejected 865 of 985 events, by far the harshest criterion, and it is
// the weakest justified: elevation does not itself alter the spectral
// response of char or crude. It only proxies for cover and soil, and
// both of those are matched directly by the NDVI envelope and the soil
// texture class. Keeping it would discard most of the sample to
// control for something already controlled. Set true to reinstate, and
// report the choice either way.
var MATCH_ELEVATION    = false;
var ELEV_TOLERANCE_M   = 300;   // used only if MATCH_ELEVATION is true
// Slope: relaxed from 5 to 15 degrees. It rejected 15 of 18 EFFIS
// candidates, because wildfire in this region burns hill country while
// the MODIS agricultural patches were on the plain. Slope acts on the
// spectra through illumination, and illumination largely cancels in the
// differenced features that carry this analysis, so 5 degrees was
// stricter than the physics requires. Recorded per event either way.
// Tightened to 5 degrees now that the candidate pool supports it.
// The Turkiye-wide distribution has median slope 12.4 deg (p25 6.6,
// p75 18.4) while Narli is 0.8 deg, so a 15 deg cutoff was comparing
// hillside fires against a flat basin floor. Roughly a fifth of
// candidates fall below 5 deg, which still leaves about a thousand.
var SLOPE_MAX_DEG      = 5.0;
// Soil texture: no longer a filter. OpenLandMap at 250 m cannot
// resolve substrate within this basin, and requiring an exact class
// match rejected 11 of 18 candidates on the strength of a product
// coarser than the events themselves. The class is still recorded on
// every event and every sample row, so it can be used as a covariate
// or a stratifier in the analysis. This is the decision flagged
// earlier: drop it rather than loosen it, and state that it was
// dropped.
// Reinstated. The Turkiye-wide soil_mode histogram is 4188 events on
// class 7 and 1159 on class 4, which is Narli's class, so requiring a
// match is affordable here in a way it was not with 18 candidates.
// Matching is done on the ROUNDED class: ee.Reducer.mode on a
// reprojected image returns near-integer floats (6.9999999999999964),
// and an exact equality test would silently reject genuine matches.
var MATCH_SOIL_TEXTURE = true;
// Pre-event NDVI envelope percentiles.
// The first run gave [0.090, 0.699] over the full 3 x 3 km window,
// which is so wide it excludes almost nothing: the window contains both
// bare basin floor and irrigated fields. Two options, both defensible,
// pick one and report it:
//   (a) tighten the percentiles (p25-p75), or
//   (b) measure the envelope from the plume's immediate surroundings
//       rather than the whole window (NDVI_ENV_RADIUS below).
// Widened from p25-p75 to p10-p90. The narrow envelope rejected 13 of
// 18 candidates. This remains the one substrate criterion doing real
// work, since pre-fire cover is the direct spectral confounder.
// Back to the interquartile envelope. On the capped subset the
// candidate NDVI spread was p25 0.21, p50 0.33, p75 0.49 against a
// Narli envelope of [0.17, 0.47], so p10-p90 was admitting cover well
// above the spill site.
var NDVI_PCT_LOW       = 25;
var NDVI_PCT_HIGH      = 75;

// Radius in metres over which the NDVI envelope is measured. Set to
// null to use the full OIL_WINDOW as before.
var NDVI_ENV_RADIUS    = 400;

// --- Burn perimeter source -----------------------------------------
// MCD64A1 at 500 m proved inadequate as a label source: on the first
// full Stage B run only 5 of 25 events showed a burn signal against
// their own control ring (median per-event AUC on d_NBR was 0.624), and
// 5 events scored below 0.5, meaning the char sat in the control
// annulus rather than in the labelled polygon. That is geolocation
// error, not weak fire, and eroding the polygon made those cases worse.
//
// Set PERIMETER_ASSET to a FeatureCollection of burn perimeters from an
// independent source: EFFIS burnt-area polygons downloaded from the
// portal and ingested as an asset, or perimeters digitised by hand on a
// Sentinel-2 SWIR composite. Both are independent of Sentinel-2
// reflectance indices, so neither contaminates the features under test.
//
// Do NOT screen or delineate with dNBR, dBAI or any other burn index.
// Those are the features being evaluated; selecting the burn class with
// them would guarantee the result.
//
// Each feature needs a date property named by PERIMETER_DATE_PROP,
// formatted 'YYYY-MM-DD'. Leave PERIMETER_ASSET null to fall back to
// MCD64A1, which is adequate for scoping but not for the analysis.
var PERIMETER_ASSET = 'projects/propane-forge-479411-k2/assets/effis_layer';

// EFFIS field names, read off the uploaded attribute table.
var PERIMETER_DATE_PROP     = 'FIREDATE';   // fire start
var PERIMETER_END_DATE_PROP = 'FINALDATE';  // fire end; null to reuse start
var PERIMETER_AREA_PROP     = 'AREA_HA';    // EFFIS published area; null to compute

// EFFIS date fields carry a time component with fractional seconds,
// e.g. '2024-09-08 08:40:23.496', and the precision is not consistent
// across records. Rather than chase format strings, the parser takes
// the leading 10 characters ('YYYY-MM-DD') and ignores the time
// entirely. Scene selection works at daily resolution anyway, so the
// time of day is irrelevant here.
//
// Set to a format string only if a perimeter source uses some layout
// where the first 10 characters are not an ISO date.
var PERIMETER_DATE_FORMAT = null;

// EFFIS reports the land cover composition of each burn as percentages
// (AGRIAREAS, BROADLEA, CONIFER, SCLEROPH, TRANSIT, MIXED, ARTIFSURF,
// OTHERNATLC, OTHERLC). AGRIAREAS is the useful one here: a burn that
// was mostly agricultural is char on cultivated soil, which is far
// closer to crude on bare ground than a canopy fire is. This is a land
// cover attribute from an external product, so filtering on it does
// not touch the spectral indices under test.
//
// Set to 0 to disable. 50 keeps burns at least half agricultural.
// Recorded, not filtered. Requiring 50 % agricultural cover cut 160
// perimeters to 21 and, combined with the substrate envelope, left
// nothing. agri_pct stays on every row so agricultural and canopy
// burns can be compared directly, which is a more useful result than
// excluding one of them. Set above 0 to reinstate the hard filter.
var MIN_AGRIAREAS_PCT = 0;

// Land cover fields carried through to the sample rows so composition
// can be reported and used for stratification in the analysis.
var LC_PROPS = ['AGRIAREAS', 'BROADLEA', 'CONIFER', 'SCLEROPH',
                'TRANSIT', 'MIXED', 'ARTIFSURF', 'OTHERNATLC', 'OTHERLC'];

// --- Burn event filtering ------------------------------------------
// MCD64A1 pixels are 500 m, i.e. 25 ha each, so the smallest possible
// patch is one pixel and this 2 ha floor never actually binds. Kept
// only as a guard against sliver polygons at the region boundary.
var BURN_MIN_AREA_HA = 2.0;
var BURN_MAX_AREA_HA = 5000.0;  // exclude megafires: substrate too heterogeneous
var BURN_MONTH_MIN   = 4;       // restrict to the fire season so that
var BURN_MONTH_MAX   = 10;      // pre/post pairs are cloud-viable

// --- Burn label refinement -----------------------------------------
// A vectorised 500 m MODIS pixel is a square, not the shape of the burn
// inside it, so its edge pixels are mixed. Eroding the polygon inward
// keeps only the interior, where the burn label is most likely correct.
//
// CRITICAL: do NOT refine these polygons with dNBR, dBAI or any other
// burn index. Those indices are the features under test. Defining the
// burn class by a large NBR drop, then reporting that oil confuses with
// burn because oil also drops NBR, is circular and biased toward the
// hypothesis. Refinement must use geometry or an independent sensor
// only. Erosion at 100 m leaves a 300 m core from a single MODIS pixel,
// about 225 Sentinel-2 pixels at 20 m, which still clears the sampling
// cap for most events.
// Erosion is appropriate for well-delineated perimeters (EFFIS, manual)
// where the edge is real and only mixed pixels need trimming. It is
// actively harmful for MODIS squares that are offset from the fire, as
// it shrinks the core further away from the char. Reduced from 100 m;
// raise it again once PERIMETER_ASSET is in use.
var BURN_ERODE_M = 40;

// Control ring: annulus around the original (uneroded) polygon. The
// inner gap must exceed the MODIS geolocation uncertainty, or controls
// will contain unlabelled burn.
var CTRL_INNER_GAP_M = 300;
var CTRL_OUTER_M     = 1200;

// --- Sentinel-2 scene selection ------------------------------------
var CS_PLUS_BAND   = 'cs_cdf';  // Cloud Score+ band
var CS_PLUS_THRESH = 0.60;      // per-pixel clear threshold
var CLEAR_MIN      = 0.90;      // min clear fraction over the AOI
var PAIR_WINDOW_D  = 45;        // search radius for pre and post scenes

// --- Sampling ------------------------------------------------------
var PIX_PER_EVENT   = 300;      // cap so one large fire cannot dominate
var CONTROL_PER_EVT = 300;      // unburned matched-substrate controls
var SAMPLE_SCALE    = 20;       // metres; 20 m is the coarsest band used
var EXPORT_FOLDER   = 'oil_burn_separability';

// Event IDs for Stage B. Leave empty to use every event in the
// inventory asset, which is what you want when all matched events are
// going into the analysis. Populate it only to run a subset, for
// example to test one event or to exclude the large outlier.
var SELECTED_EVENTS = [];

// Maximum events Stage B will process. toList needs a client-side
// number, so this is a ceiling rather than a filter; 25 matched events
// leaves plenty of headroom.
var MAX_EVENTS = 200;

// Inventory asset exported by Stage A.
// Earth Engine refuses to overwrite an existing asset: re-exporting to
// burn_inventory_matched failed with "An asset with the ID already
// exists", the export never ran, and Stage B silently kept reading the
// old MODIS inventory. Version the asset name per perimeter source
// instead of overwriting, which also keeps the MODIS and EFFIS
// inventories side by side for comparison.
var INVENTORY_ASSET =
  'projects/propane-forge-479411-k2/assets/burn_inventory_effis';

/* =====================================================================
 * PART 2  SUBSTRATE FINGERPRINT OF THE OIL SITE
 * Everything the burn scars have to match is measured here, not assumed.
 * ===================================================================== */

// GLO30 is delivered as an ImageCollection of tiles, not a single
// Image. It must be mosaicked, and the mosaic must be given an explicit
// projection before ee.Terrain touches it: on an unprojected mosaic
// slope silently returns 0 rather than failing.
var demCol = ee.ImageCollection('COPERNICUS/DEM/GLO30');
var dem = demCol.select('DEM').mosaic()
            .setDefaultProjection('EPSG:3857', null, 30)
            .rename('DEM');
var terrain = ee.Terrain.products(dem);
var slope = terrain.select('slope');

// OpenLandMap USDA texture class, 0-5 cm.
var soilTexture = ee.Image('OpenLandMap/SOL/SOL_TEXTURE-CLASS_USDA-TT_M/v02')
                    .select('b0');

// Narli elevation and slope.
var oilElev = dem.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: OIL_POINT.buffer(200),
  scale: 30,
  maxPixels: 1e9
}).get('DEM');

var oilSlope = slope.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: OIL_POINT.buffer(200),
  scale: 30,
  maxPixels: 1e9
}).get('slope');

// Modal soil texture class over the analysis window.
var oilSoilClass = soilTexture.reduceRegion({
  reducer: ee.Reducer.mode(),
  geometry: OIL_WINDOW,
  scale: 250,
  maxPixels: 1e9
}).get('b0');

// Rounded: the mode reducer on a reprojected image returns floats a
// hair off the integer class value.
oilSoilClass = ee.Number(oilSoilClass).round();

print('--- PART 2  Narli substrate fingerprint ---');
print('Elevation (m), 200 m radius:', oilElev);
print('Mean slope (deg), 200 m radius:', oilSlope);
print('Modal soil texture class (USDA):', oilSoilClass);

/* ---- Sentinel-2 access with Cloud Score+ masking ------------------ */

var S2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED');
var CSP = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');

var S2_BANDS = ['B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12'];

/**
 * Cloud-masked, scaled Sentinel-2 collection over a geometry and window.
 * Reflectance is returned as a fraction (0-1), not the DN 0-10000.
 */
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

/**
 * Pick the single scene nearest to a date with adequate clear coverage.
 * before = true  -> latest clear scene strictly before the date
 * before = false -> earliest clear scene strictly after the date
 *
 * Returns an ee.Image, or an image with no bands if nothing qualifies.
 * Single scenes rather than composites: compositing alters the pixel
 * values and changes the noise characteristics of a differenced pair.
 */
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

/* ---- Pre-event NDVI envelope at the oil site --------------------- */

var oilPre = pickScene(OIL_WINDOW, OIL_EVENT_DATE, true);

var oilNdvi = oilPre.normalizedDifference(['B8', 'B4']).rename('NDVI');

var NDVI_ENV_GEOM = NDVI_ENV_RADIUS === null
  ? OIL_WINDOW
  : OIL_POINT.buffer(NDVI_ENV_RADIUS);

var ndviEnv = oilNdvi.reduceRegion({
  reducer: ee.Reducer.percentile([NDVI_PCT_LOW, NDVI_PCT_HIGH]),
  geometry: NDVI_ENV_GEOM,
  scale: 20,
  maxPixels: 1e9
});

var NDVI_LO = ee.Number(ndviEnv.get('NDVI_p' + NDVI_PCT_LOW));
var NDVI_HI = ee.Number(ndviEnv.get('NDVI_p' + NDVI_PCT_HIGH));

print('Pre-event scene date:', oilPre.date().format('YYYY-MM-dd'));
print('Pre-event NDVI envelope [p' + NDVI_PCT_LOW + ', p' + NDVI_PCT_HIGH + ']:',
      NDVI_LO, NDVI_HI);

/* =====================================================================
 * PART 3  BURN INVENTORY
 * Two burned-area products plus active-fire points. MCD64A1 and
 * FireCCI51 are both 250-500 m and will miss small stubble burns, which
 * is precisely the class most likely to resemble oil on bare soil. The
 * FIRMS points are included to catch those.
 * ===================================================================== */

// Study region.
// Name matching for Kahramanmaras failed twice: exact match, then the
// diacritic-free substring 'Kahramanmara'. GAUL 2015 evidently stores
// it under some other form. Rather than guess again, select the
// province that spatially contains the rupture point. That is
// definitionally the right polygon whatever it is called.
var gaul = ee.FeatureCollection('FAO/GAUL_SIMPLIFIED_500m/2015/level1');
var turkey = gaul.filter(ee.Filter.eq('ADM0_NAME', 'Turkey'));

var homeProvince = turkey.filterBounds(OIL_POINT);

var namedProvinces = turkey.filter(ee.Filter.or(
  ee.Filter.stringContains('ADM1_NAME', 'Adana'),
  ee.Filter.stringContains('ADM1_NAME', 'Osmaniye')
));

var provinces = homeProvince.merge(namedProvinces).distinct('ADM1_NAME');

var countryUnits = gaul.filter(
  ee.Filter.inList('ADM0_NAME', REGION_COUNTRIES));

// Diagnostic. Click the List to expand it; the name GAUL actually uses
// for Kahramanmaras will be visible, and can go into PROVINCE_NAMES for
// the record.
print('All Turkiye ADM1_NAME values:',
      turkey.aggregate_array('ADM1_NAME').sort());
print('Province containing the rupture point:',
      homeProvince.aggregate_array('ADM1_NAME'));
print('Matched province names:', provinces.aggregate_array('ADM1_NAME'));

// Resolve the region from REGION_SCOPE, with the bbox as a fallback if
// a name-based selection comes back empty.
var REGION;
if (REGION_SCOPE === 'bbox') {
  REGION = REGION_BBOX;
} else if (REGION_SCOPE === 'country') {
  REGION = ee.Geometry(ee.Algorithms.If(
    countryUnits.size().gt(0), countryUnits.geometry(), REGION_BBOX));
} else {
  REGION = ee.Geometry(ee.Algorithms.If(
    provinces.size().gt(0), provinces.geometry(), REGION_BBOX));
}

print('--- PART 3  Burn inventory ---');
print('REGION_SCOPE:', REGION_SCOPE);
print('GAUL province matches (0 means bbox fallback in use):', provinces.size());
print('GAUL units in REGION_COUNTRIES:', countryUnits.size());
print('Region area (ha):', REGION.area(1000).divide(1e4));

var DATE_START = ee.Date.fromYMD(YEAR_START, 1, 1);
var DATE_END   = ee.Date.fromYMD(YEAR_END, 9, 1);

// MODIS MCD64A1: BurnDate band gives day-of-year of burn.
var mcd64 = ee.ImageCollection('MODIS/061/MCD64A1')
              .filterDate(DATE_START, DATE_END)
              .select('BurnDate');

// ESA FireCCI51 ends earlier than MCD64A1; filterDate handles that.
var firecci = ee.ImageCollection('ESA/CCI/FireCCI/5_1')
                .filterDate(DATE_START, DATE_END)
                .select('BurnDate');

// FIRMS active fire detections.
var firms = ee.ImageCollection('FIRMS')
              .filterDate(DATE_START, DATE_END)
              .select('T21');

/**
 * Vectorise burned pixels for one monthly burned-area image into
 * discrete events, each carrying its burn date.
 */
function burnEventsFromImage(img) {
  var t = ee.Date(img.get('system:time_start'));
  var year = t.get('year');
  var burned = img.gt(0).selfMask().rename('burned');

  // Median burn day-of-year per connected patch, converted to a date.
  // setOutputs is not optional here. reduceToVectors names the reduced
  // property after the reducer's own output name, not after the input
  // band, so ee.Reducer.median() produced a property called 'median'
  // and 'doy' was null on every patch. Before the notNull guard that
  // surfaced as a subtract-on-null crash; after it, it silently threw
  // away all 825 patches. Naming the output explicitly fixes both.
  var vectors = burned.addBands(img.rename('doy')).reduceToVectors({
    geometry: REGION,
    scale: 500,
    geometryType: 'polygon',
    eightConnected: true,
    labelProperty: 'patch',
    reducer: ee.Reducer.median().setOutputs(['doy']),
    maxPixels: 1e10
  });

  // Patches clipped by the region boundary can still end up with no
  // unmasked day-of-year pixel. Drop them rather than defaulting the
  // date, which would invent a burn date.
  vectors = vectors.filter(ee.Filter.notNull(['doy']));

  return vectors.map(function (f) {
    var doy = ee.Number(f.get('doy'));
    var burnDate = ee.Date.fromYMD(year, 1, 1).advance(doy.subtract(1), 'day');
    var areaHa = f.geometry().area(30).divide(1e4);
    return f.set({
      burn_date: burnDate.format('YYYY-MM-dd'),
      burn_year: year,
      burn_month: burnDate.get('month'),
      burn_doy: doy,
      area_ha: areaHa,
      source: 'MCD64A1'
    });
  });
}

var rawEvents;
if (PERIMETER_ASSET === null) {
  rawEvents = ee.FeatureCollection(mcd64.map(burnEventsFromImage)).flatten();
} else {
  var perim = ee.FeatureCollection(PERIMETER_ASSET);

  // Inspect the upload before trusting it. Wrong field name or an
  // unparseable date shows up here rather than as an empty result
  // twenty lines later.
  print('Perimeter features (all):', perim.size());
  print('Perimeter property names:',
        ee.Feature(perim.first()).propertyNames());
  print('Perimeter sample feature:', perim.first());

  // EFFIS records carry empty strings in numeric and date fields for
  // some features, and ee.Number.parse('') throws, which kills the
  // whole collection over one bad record. Every attribute read below
  // goes through a defaulted lookup.
  //
  // Feature.toDictionary().get(key, default) supplies the default when
  // the property is absent; the length check then catches present but
  // empty values. A stored 0 survives, because the conversion to string
  // happens before any truthiness test.
  var safeNum = function (f, prop, fallback) {
    var raw = ee.Dictionary(f.toDictionary()).get(prop, '');
    var s = ee.Algorithms.String(raw);
    return ee.Number(ee.Algorithms.If(
      s.length().gt(0), ee.Number.parse(s), fallback));
  };

  var safeDate = function (f, prop, fallback) {
    var raw = ee.Dictionary(f.toDictionary()).get(prop, '');
    var s = ee.Algorithms.String(raw);
    return ee.Date(ee.Algorithms.If(
      s.length().gte(10),
      PERIMETER_DATE_FORMAT === null
        ? ee.Date(s.slice(0, 10))
        : ee.Date.parse(PERIMETER_DATE_FORMAT, s),
      fallback));
  };

  // Features with no usable fire date cannot be paired at all, so they
  // are dropped up front rather than given a fabricated date.
  var inRegion = perim.filterBounds(REGION);
  var dated = inRegion
    .filter(ee.Filter.notNull([PERIMETER_DATE_PROP]))
    .filter(ee.Filter.neq(PERIMETER_DATE_PROP, ''));

  print('Perimeter features in region:', inRegion.size());
  print('  with a usable fire date:', dated.size());

  rawEvents = dated.map(function (f) {
    var d = safeDate(f, PERIMETER_DATE_PROP, ee.Date('1970-01-01'));

    // Fire end date. A multi-day fire needs its post scene after the
    // fire stopped, not after it started, or the post image can catch
    // the burn mid-progress. Falls back to the start date if absent.
    var dEnd = PERIMETER_END_DATE_PROP === null
      ? d
      : safeDate(f, PERIMETER_END_DATE_PROP, d);

    // If FINALDATE precedes FIREDATE in a malformed record, use the
    // start date rather than a negative duration.
    dEnd = ee.Date(ee.Algorithms.If(dEnd.millis().gte(d.millis()), dEnd, d));

    // area_ha falls back to the polygon's own area when AREA_HA is
    // blank, so a missing attribute never removes an event.
    var areaGeom = f.geometry().area(30).divide(1e4);
    var area = PERIMETER_AREA_PROP === null
      ? areaGeom
      : safeNum(f, PERIMETER_AREA_PROP, areaGeom);

    // agri_pct uses -1 for "not reported", which is distinguishable
    // from a genuine 0. Note that MIN_AGRIAREAS_PCT > 0 therefore also
    // excludes unreported events; that is intended, but it means the
    // filter count includes missing data, not only canopy fires.
    return f.set({
      burn_date: d.format('YYYY-MM-dd'),
      burn_end_date: dEnd.format('YYYY-MM-dd'),
      burn_duration_d: dEnd.difference(d, 'day'),
      burn_year: d.get('year'),
      burn_month: d.get('month'),
      burn_doy: d.getRelative('day', 'year').add(1),
      area_ha: area,
      agri_pct: safeNum(f, 'AGRIAREAS', -1),
      source: 'PERIMETER_ASSET'
    });
  });

  // AGRIAREAS distribution. On the first EFFIS run the filter removed
  // nothing (160 of 160 passed), which is worth verifying rather than
  // assuming: a string-typed property would compare unexpectedly.
  // A spike at -1 means the field is mostly unreported.
  print('AGRIAREAS distribution (-1 = not reported):',
        rawEvents.aggregate_histogram('agri_pct'));
  print('Fire duration (days) distribution:',
        rawEvents.aggregate_histogram('burn_duration_d'));

  if (MIN_AGRIAREAS_PCT > 0) {
    var beforeAgri = rawEvents.size();
    rawEvents = rawEvents.filter(
      ee.Filter.gte('agri_pct', MIN_AGRIAREAS_PCT));
    print('Before AGRIAREAS filter:', beforeAgri);
    print('After AGRIAREAS >= ' + MIN_AGRIAREAS_PCT + '%:', rawEvents.size());
  }
}

print('Perimeter source:',
      PERIMETER_ASSET === null ? 'MCD64A1 (scoping only)' : PERIMETER_ASSET);

var sizedEvents = rawEvents
  .filter(ee.Filter.gte('area_ha', BURN_MIN_AREA_HA))
  .filter(ee.Filter.lte('area_ha', BURN_MAX_AREA_HA))
  .filter(ee.Filter.gte('burn_month', BURN_MONTH_MIN))
  .filter(ee.Filter.lte('burn_month', BURN_MONTH_MAX));

print('Raw burn patches:', rawEvents.size());

// Property-name check. A silently renamed property is what turned 825
// patches into 0. If 'doy' is absent from this list, read off the real
// name and put it in setOutputs above.
print('Property names on a raw patch:',
      ee.Feature(rawEvents.first()).propertyNames());
print('Sample raw patch:', rawEvents.first());

print('After area and season filters:', sizedEvents.size());

/* =====================================================================
 * PART 4  SUBSTRATE MATCHING
 * A burn event survives only if its own substrate falls inside the
 * Narli envelope. Matching is done on the event mean, then again
 * per-pixel at sampling time in Stage B.
 * ===================================================================== */

var elevLo = ee.Number(oilElev).subtract(ELEV_TOLERANCE_M);
var elevHi = ee.Number(oilElev).add(ELEV_TOLERANCE_M);

/**
 * Null-safe numeric coalesce. reduceRegion returns null for a band when
 * no valid pixel centre falls inside the geometry, which happens for
 * small polygons and over DEM or soil no-data (coastal fires, water).
 * One null then killed the whole collection on Number.lte.
 *
 * ee.Algorithms.If cannot be used here because it treats a genuine 0 as
 * false; firstNonNull distinguishes null from zero correctly.
 *
 * Fallbacks should be chosen so a missing value FAILS the match rather
 * than passing it silently.
 */
function orElse(v, fallback) {
  return ee.Number(ee.List([v, fallback]).reduce(ee.Reducer.firstNonNull()));
}

// Sentinels that fail their respective tests.
var NO_ELEV  = -9999;   // fails the elevation floor
var NO_SLOPE = 9999;    // fails slope <= max
var NO_SOIL  = -1;      // fails the class equality
var NO_NDVI  = -9;      // fails the NDVI floor

/**
 * Phase 1, cheap. Terrain and soil only: raster lookups, no imagery.
 * Runs over every candidate however many there are.
 */
function annotateTerrain(f) {
  var geom = f.geometry();

  // bestEffort is deliberately NOT set. It silently coarsens the scale
  // to fit the pixel budget, and the scale it picks can differ between
  // runs, which made pre_ndvi shift between bins on repeat runs of an
  // otherwise identical script. Fixed scales with an ample maxPixels
  // are reproducible; if a reduction now fails on pixel count, raise
  // maxPixels rather than reinstating bestEffort.
  var elevMean = orElse(dem.reduceRegion({
    reducer: ee.Reducer.mean(), geometry: geom, scale: 30,
    maxPixels: 1e10}).get('DEM'), NO_ELEV);
  var slopeMean = orElse(slope.reduceRegion({
    reducer: ee.Reducer.mean(), geometry: geom, scale: 30,
    maxPixels: 1e10}).get('slope'), NO_SLOPE);
  var soilMode = orElse(soilTexture.reduceRegion({
    reducer: ee.Reducer.mode(), geometry: geom, scale: 250,
    maxPixels: 1e10}).get('b0'), NO_SOIL).round();

  // Flag rather than hide missing substrate data, so the count of
  // events dropped for no-data is reportable and not confused with
  // events dropped for genuine mismatch.
  var dataOk = ee.Number(slopeMean).neq(NO_SLOPE)
                 .and(ee.Number(elevMean).neq(NO_ELEV))
                 .and(ee.Number(soilMode).neq(NO_SOIL));

  var elevOk  = MATCH_ELEVATION
                  ? ee.Number(elevMean).gte(elevLo)
                      .and(ee.Number(elevMean).lte(elevHi))
                  : ee.Number(1);
  var slopeOk = ee.Number(slopeMean).lte(SLOPE_MAX_DEG);
  var soilOk  = MATCH_SOIL_TEXTURE
                  ? ee.Number(soilMode).eq(ee.Number(oilSoilClass))
                  : ee.Number(1);

  return f.set({
    elev_mean: elevMean,
    slope_mean: slopeMean,
    soil_mode: soilMode,
    substrate_data_ok: dataOk,
    ok_elev: elevOk,
    ok_slope: slopeOk,
    ok_soil: soilOk,
    terrain_ok: elevOk.and(slopeOk).and(soilOk).and(dataOk)
  });
}

/**
 * Phase 2, expensive. Pre-fire NDVI from a median composite of the 90
 * days before the burn. The year before the burn is deliberately not
 * used: interannual substitution would confound land use change with
 * the substrate match.
 *
 * This is the costly step, so it runs only on candidates that already
 * passed terrain and soil, capped at CANDIDATE_LIMIT.
 */
function annotateNdvi(f) {
  var geom = f.geometry();
  var burnDate = ee.Date(f.get('burn_date'));

  var preCol = s2Clear(geom, burnDate.advance(-90, 'day'),
                             burnDate.advance(-2, 'day'));
  var preNdviImg = preCol.median().normalizedDifference(['B8', 'B4'])
                         .rename('NDVI');
  var preNdvi = orElse(preNdviImg.reduceRegion({
    reducer: ee.Reducer.mean(), geometry: geom, scale: 20,
    maxPixels: 1e10}).get('NDVI'), NO_NDVI);

  var nPre = preCol.size();
  var ndviOk = ee.Number(preNdvi).gte(NDVI_LO).and(ee.Number(preNdvi).lte(NDVI_HI));

  return f.set({
    pre_ndvi: preNdvi,
    n_pre_scenes: nPre,
    ok_ndvi: ndviOk,
    matched: ee.Number(f.get('terrain_ok')).and(ndviOk).and(nPre.gt(0))
  });
}

// Phase 1 over everything, then cap, then phase 2.
var terrainAnnotated = sizedEvents.map(annotateTerrain);
var terrainPassed = terrainAnnotated.filter(ee.Filter.eq('terrain_ok', 1));

// Seeded random subsample so the cap is reproducible rather than
// dependent on collection order.
var capped = terrainPassed
  .randomColumn('rnd', CANDIDATE_SEED)
  .sort('rnd')
  .limit(CANDIDATE_LIMIT);

var annotated = capped.map(annotateNdvi);
var matchedAll = annotated.filter(ee.Filter.eq('matched', 1));

// Seeded subsample so the cap is reproducible rather than dependent on
// collection order.
var matchedEvents = MATCHED_EVENT_LIMIT > 0
  ? matchedAll.randomColumn('rnd2', CANDIDATE_SEED + 1)
              .sort('rnd2').limit(MATCHED_EVENT_LIMIT)
  : matchedAll;

print('Candidates after area and season filters:', sizedEvents.size());
print('Passed terrain and soil (cheap phase):', terrainPassed.size());
print('Entering NDVI phase (capped at ' + CANDIDATE_LIMIT + '):',
      capped.size());

// Stable event IDs, derived from each event's own properties.
// Earlier versions sorted the collection, pulled it to a List, and
// indexed into it. Both the indexOf and the sequence form force the
// whole annotated chain to evaluate at once, which kept exceeding the
// interactive memory limit even for .limit(5). Encoding year, day of
// year and rounded centroid coordinates needs no cross-feature
// dependency at all, is deterministic, and carries the location in the
// identifier, which is better provenance than a running number.
function addEventId(f) {
  f = ee.Feature(f);
  var c = f.geometry().centroid(10).coordinates();
  var lon = ee.Number(c.get(0)).multiply(1000).round().format('%d');
  var lat = ee.Number(c.get(1)).multiply(1000).round().format('%d');
  var y = ee.Number(f.get('burn_year')).format('%d');
  var d = ee.Number(f.get('burn_doy')).round().format('%03d');
  return f.set('event_id',
    ee.String('BURN_').cat(y).cat('_').cat(d).cat('_').cat(lon).cat('_').cat(lat));
}

var matchedIded = matchedEvents.map(addEventId);

/* =====================================================================
 * PART 5  STAGE A DIAGNOSTICS
 * The feasibility question: are there enough matched burns to train on?
 * ===================================================================== */

if (RUN_STAGE === 'A') {
  // Substrate distributions across ALL candidates, before matching.
  // Printed so the tolerances above can be set from the data rather
  // than guessed and then revised after each failed run. Compare these
  // against the Narli fingerprint in PART 2.
  print('--- candidate substrate distributions (pre-match) ---');
  var withData = terrainAnnotated.filter(ee.Filter.eq('substrate_data_ok', 1));
  print('candidates with usable substrate data:', withData.size());
  print('slope_mean (deg) quartiles:',
        withData.reduceColumns(
          ee.Reducer.percentile([0, 25, 50, 75, 100]), ['slope_mean']));
  print('soil_mode classes:',
        withData.aggregate_histogram('soil_mode'));
  print('area_ha quartiles:',
        terrainAnnotated.reduceColumns(
          ee.Reducer.percentile([0, 25, 50, 75, 100]), ['area_ha']));
  print('agri_pct quartiles:',
        terrainAnnotated.reduceColumns(
          ee.Reducer.percentile([0, 25, 50, 75, 100]), ['agri_pct']));
  print('pre_ndvi quartiles (capped subset):',
        annotated.reduceColumns(
          ee.Reducer.percentile([0, 25, 50, 75, 100]), ['pre_ndvi']));

  print('--- PART 5  Stage A diagnostics ---');
  print('Matched events (all):', matchedAll.size());
  print('Matched events (after cap of ' + MATCHED_EVENT_LIMIT + '):',
        matchedEvents.size());
  print('Provinces represented:',
        matchedEvents.aggregate_histogram('PROVINCE'));
  print('Total matched burn area (ha):',
        matchedEvents.aggregate_sum('area_ha'));

  // Approximate Sentinel-2 pixel count at 20 m: 25 pixels per hectare.
  print('Approx. matched burn pixels at 20 m:',
        ee.Number(matchedEvents.aggregate_sum('area_ha')).multiply(25).round());

  // Why events were rejected. If one criterion dominates, reconsider it
  // rather than quietly widening every tolerance.
  // Terrain rejections are counted over all candidates; NDVI
  // rejections only over the capped subset that reached phase 2.
  print('Rejected on elevation (of all):',
        terrainAnnotated.filter(ee.Filter.eq('ok_elev', 0)).size());
  print('Rejected on slope (of all):',
        terrainAnnotated.filter(ee.Filter.eq('ok_slope', 0)).size());
  print('Rejected on soil (of all):',
        terrainAnnotated.filter(ee.Filter.eq('ok_soil', 0)).size());
  print('Dropped for missing substrate data (of all):',
        terrainAnnotated.filter(ee.Filter.eq('substrate_data_ok', 0)).size());
  print('Rejected on NDVI (of capped):',
        annotated.filter(ee.Filter.eq('ok_ndvi', 0)).size());
  print('No pre-fire imagery (of capped):',
        annotated.filter(ee.Filter.eq('n_pre_scenes', 0)).size());

  print(ui.Chart.feature.histogram(matchedEvents, 'area_ha', 30)
          .setOptions({title: 'Matched burn event area (ha)'}));

  print(ui.Chart.feature.histogram(matchedEvents, 'pre_ndvi', 30)
          .setOptions({title: 'Pre-fire NDVI of matched events'}));

  // The by-year chart and the 100-row candidate table both exceeded the
  // interactive memory limit: annotateEvent runs several reduceRegions
  // and a Sentinel-2 median composite per event, and the console
  // recomputes the whole chain for every widget. Batch exports have no
  // such limit, so the inventory goes to Drive and to an asset instead.
  //
  // Run the asset export first, then set INVENTORY_ASSET below and
  // re-run. Stage B will read the finished asset rather than rebuilding
  // the inventory, which is what makes Stage B tractable at all.
  // assetId must match INVENTORY_ASSET. Bump the suffix whenever the
  // perimeter source or the substrate envelope changes, rather than
  // trying to overwrite: the overwrite fails and the stale asset stays
  // in place.
  Export.table.toAsset({
    collection: matchedIded,
    description: 'burn_inventory_asset_effis',
    assetId: INVENTORY_ASSET
  });

  Export.table.toDrive({
    collection: matchedIded,
    description: 'burn_inventory_effis',
    folder: EXPORT_FOLDER,
    fileFormat: 'CSV'
  });

  // No candidate table is printed. Even .limit(5) forced the full
  // annotated chain to evaluate and exceeded the memory limit. Take the
  // event_id list from the exported CSV instead, which is the same data
  // without the interactive budget.
  print('Inventory will be exported. Read event_id values from the CSV.');

  Map.centerObject(OIL_POINT, 11);
  Map.addLayer(REGION, {color: 'grey'}, 'Study region', false);
  Map.addLayer(matchedEvents, {color: 'red'}, 'Matched burn events');
  Map.addLayer(OIL_WINDOW, {color: 'blue'}, 'Narli analysis window');
}

/* =====================================================================
 * PART 6  INDEX LIBRARY
 * Applied identically to oil and burn. Monotone-equivalent formulations
 * are noted: they cannot differ in AUC, and reporting them as separate
 * discriminators inflates the apparent size of the feature space.
 * ===================================================================== */

function addIndices(img) {
  var B2  = img.select('B2'),  B3  = img.select('B3');
  var B4  = img.select('B4'),  B5  = img.select('B5');
  var B6  = img.select('B6'),  B7  = img.select('B7');
  var B8  = img.select('B8'),  B8A = img.select('B8A');
  var B11 = img.select('B11'), B12 = img.select('B12');

  // --- Burn family ---
  var nbr  = B8.subtract(B12).divide(B8.add(B12)).rename('NBR');
  var nbr2 = B11.subtract(B12).divide(B11.add(B12)).rename('NBR2');

  // BAI is a red/NIR index, not a SWIR index. Distance in the
  // red-NIR plane from a reference charcoal point.
  var bai = ee.Image(1).divide(
      B4.subtract(0.1).pow(2).add(B8.subtract(0.06).pow(2))
    ).rename('BAI');
  var logbai = bai.log10().rename('logBAI');

  // MIRBI: designed for burn mapping in sparsely vegetated terrain,
  // which is the relevant regime here. Coefficients assume reflectance
  // as a fraction.
  var mirbi = B12.multiply(10).subtract(B11.multiply(9.8)).add(2).rename('MIRBI');

  // BAIS2: Sentinel-2 specific burned area index.
  var bais2 = ee.Image(1)
      .subtract(B6.multiply(B7).multiply(B8A).divide(B4).sqrt())
      .multiply(
        B12.subtract(B8A).divide(B12.add(B8A).sqrt()).add(1)
      ).rename('BAIS2');

  // --- Mineral / moisture family ---
  // FMR = B11/B8. MSI is the same ratio, and NDBI is a monotone
  // transform of it, so all three rank pixels identically. Kept for
  // reporting continuity with the literature, flagged as equivalent.
  var fmr  = B11.divide(B8).rename('FMR');
  var ndbi = B11.subtract(B8).divide(B11.add(B8)).rename('NDBI');

  // --- Broadband albedo (Liang 2001) ---
  var albedo = B2.multiply(0.356)
      .add(B4.multiply(0.130))
      .add(B8.multiply(0.373))
      .add(B11.multiply(0.085))
      .add(B12.multiply(0.072))
      .subtract(0.0018)
      .divide(1.016).rename('ALBEDO');

  // --- Vegetation / pigment / red edge ---
  var ndvi  = B8.subtract(B4).divide(B8.add(B4)).rename('NDVI');
  var dvi   = B8.subtract(B4).rename('DVI');
  var arvi  = B8.subtract(B4.multiply(2).subtract(B2))
                .divide(B8.add(B4.multiply(2).subtract(B2))).rename('ARVI');
  var ari   = ee.Image(1).divide(B3).subtract(ee.Image(1).divide(B5)).rename('ARI');
  var ireci = B7.subtract(B4).divide(B5.divide(B6)).rename('IRECI');
  var ndwi  = B8.subtract(B11).divide(B8.add(B11)).rename('NDWI');

  return img.addBands([
    nbr, nbr2, bai, logbai, mirbi, bais2,
    fmr, ndbi, albedo,
    ndvi, dvi, arvi, ari, ireci, ndwi
  ]);
}

var FEATURE_BANDS = [
  'B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12',
  'NBR','NBR2','BAI','logBAI','MIRBI','BAIS2',
  'FMR','NDBI','ALBEDO',
  'NDVI','DVI','ARVI','ARI','IRECI','NDWI'
];

/**
 * Build a single-date and differenced feature stack for one event.
 * Differenced bands are prefixed 'd_' and are post minus pre, matching
 * the sign convention of the earlier study.
 */
function featureStack(geom, preAnchor, postAnchor) {
  postAnchor = postAnchor === undefined ? preAnchor : postAnchor;
  var pre  = addIndices(pickScene(geom, preAnchor, true));
  var post = addIndices(pickScene(geom, postAnchor, false));

  var single = post.select(FEATURE_BANDS);
  var diff = post.select(FEATURE_BANDS)
                 .subtract(pre.select(FEATURE_BANDS))
                 .rename(FEATURE_BANDS.map(function (b) { return 'd_' + b; }));

  // Per-pixel substrate covariates, attached to the same stack that
  // every class is sampled from, so oil, burn and control all carry
  // them on identical terms.
  var covars = pre.select('NDVI').rename('pre_NDVI_px')
    .addBands(slope.rename('slope_px'))
    .addBands(soilTexture.rename('soil_px'));

  return single.addBands(diff).addBands(covars)
    .set({
      pre_date:  pre.date().format('YYYY-MM-dd'),
      post_date: post.date().format('YYYY-MM-dd'),
      pre_millis: pre.date().millis(),
      post_millis: post.date().millis()
    });
}

/* =====================================================================
 * PART 7  PER-PIXEL SUBSTRATE MASK
 * Event-level matching in PART 4 is coarse. This applies the same
 * envelope pixel by pixel, so a matched event contributes only its
 * matched pixels.
 * ===================================================================== */

function substrateMask(geom, refDate) {
  var elevOk  = MATCH_ELEVATION
                  ? dem.gte(elevLo).and(dem.lte(elevHi))
                  : ee.Image.constant(1);
  var slopeOk = slope.lte(SLOPE_MAX_DEG);
  var soilOk  = MATCH_SOIL_TEXTURE
                  ? soilTexture.eq(ee.Image.constant(oilSoilClass))
                  : ee.Image.constant(1);

  var preNdvi = addIndices(pickScene(geom, refDate, true)).select('NDVI');
  var ndviOk = preNdvi.gte(ee.Image.constant(NDVI_LO))
                      .and(preNdvi.lte(ee.Image.constant(NDVI_HI)));

  return elevOk.and(slopeOk).and(soilOk).and(ndviOk).rename('substrate_ok');
}

/* =====================================================================
 * PART 8  SAMPLE EXTRACTION
 * Classes: oil, burn, and unburned matched-substrate control.
 * Per-event pixel caps prevent a single large fire from dominating.
 * ===================================================================== */

/**
 * Sample one burn event: burn pixels plus surrounding unburned controls
 * on the same substrate, from the same image pair.
 */
function sampleBurnEvent(feat) {
  feat = ee.Feature(feat);
  var geom = feat.geometry();
  var eventDate = ee.String(feat.get('burn_date'));
  var endDate = ee.String(ee.Algorithms.If(
    feat.get('burn_end_date'), feat.get('burn_end_date'), feat.get('burn_date')));
  var eventId = ee.String(feat.get('event_id'));

  // Burn core: polygon eroded inward, so mixed edge pixels of the
  // 500 m MODIS square are excluded. Purely geometric, so it introduces
  // no dependence on the indices being tested. If erosion empties a
  // single-pixel polygon the event yields no burn rows and will show up
  // as missing in the per-event row counts.
  var core = geom.buffer(-BURN_ERODE_M);

  // Control ring: annulus around the ORIGINAL polygon, not the eroded
  // core, with a gap wider than the MODIS geolocation uncertainty.
  var ring = geom.buffer(CTRL_OUTER_M)
                 .difference(geom.buffer(CTRL_INNER_GAP_M), 1);
  var aoi = geom.buffer(CTRL_OUTER_M);

  var stack = featureStack(aoi, eventDate, endDate);
  var stackMasked = PIXEL_SUBSTRATE_MASK
    ? stack.updateMask(substrateMask(aoi, eventDate))
    : stack;

  var burnSample = stackMasked.sample({
    region: core,
    scale: SAMPLE_SCALE,
    numPixels: PIX_PER_EVENT,
    seed: 42,
    dropNulls: true,
    geometries: true
  }).map(function (f) {
    return f.set({
      class_name: 'burn',
      event_id: eventId,
      event_date: eventDate,
      event_end_date: endDate,
      area_ha: feat.get('area_ha'),
      pre_ndvi: feat.get('pre_ndvi'),
      agri_pct: feat.get('agri_pct'),
      slope_mean: feat.get('slope_mean'),
      soil_mode: feat.get('soil_mode'),
      province: feat.get('PROVINCE'),
      pre_date: stack.get('pre_date'),
      post_date: stack.get('post_date')
    });
  });

  var ctrlSample = stackMasked.sample({
    region: ring,
    scale: SAMPLE_SCALE,
    numPixels: CONTROL_PER_EVT,
    seed: 43,
    dropNulls: true,
    geometries: true
  }).map(function (f) {
    return f.set({
      class_name: 'control_unburned',
      event_id: eventId,
      event_date: eventDate,
      event_end_date: endDate,
      area_ha: feat.get('area_ha'),
      pre_ndvi: feat.get('pre_ndvi'),
      agri_pct: feat.get('agri_pct'),
      slope_mean: feat.get('slope_mean'),
      soil_mode: feat.get('soil_mode'),
      province: feat.get('PROVINCE'),
      pre_date: stack.get('pre_date'),
      post_date: stack.get('post_date')
    });
  });

  return burnSample.merge(ctrlSample);
}

/**
 * Sample the oil site. Requires the plume mask asset. Includes a plume
 * edge ring, because a negative class drawn only from distant clean
 * soil returns perfect separability as an artefact of the design.
 */
function sampleOilSite() {
  if (OIL_PLUME_GEOM === null && OIL_PLUME_ASSET === null) {
    print('WARNING: no plume label set. Skipping oil samples. ' +
          'Run RUN_STAGE = \'OIL\' to digitise the plume.');
    return ee.FeatureCollection([]);
  }

  var stack = featureStack(OIL_WINDOW, OIL_EVENT_DATE);

  // A hand-drawn polygon takes precedence over a raster mask.
  var plumeGeom;
  if (OIL_PLUME_GEOM !== null) {
    plumeGeom = ee.Geometry(OIL_PLUME_GEOM);
  } else {
    plumeGeom = ee.Image(OIL_PLUME_ASSET).gt(0).selfMask()
      .reduceToVectors({
        geometry: OIL_WINDOW, scale: 10, geometryType: 'polygon',
        eightConnected: true, maxPixels: 1e10
      }).geometry();
  }

  print('Plume area used for sampling (ha):',
        plumeGeom.area(1).divide(1e4));

  var edgeRing = plumeGeom.buffer(30).difference(plumeGeom, 1);
  var farClean = OIL_WINDOW.difference(plumeGeom.buffer(100), 1);

  function tag(fc, name) {
    return fc.map(function (f) {
      return f.set({
        class_name: name,
        event_id: 'OIL_NARLI_2018',
        event_date: OIL_EVENT_DATE,
        pre_date: stack.get('pre_date'),
        post_date: stack.get('post_date')
      });
    });
  }

  var oil = tag(stack.sample({
    region: plumeGeom, scale: SAMPLE_SCALE, numPixels: 2000,
    seed: 42, dropNulls: true, geometries: true}), 'oil');

  var edge = tag(stack.sample({
    region: edgeRing, scale: SAMPLE_SCALE, numPixels: 500,
    seed: 43, dropNulls: true, geometries: true}), 'oil_edge');

  var clean = tag(stack.sample({
    region: farClean, scale: SAMPLE_SCALE, numPixels: 2000,
    seed: 44, dropNulls: true, geometries: true}), 'control_clean');

  return oil.merge(edge).merge(clean);
}

/* =====================================================================
 * PART 9  STAGE B EXPORT
 * ===================================================================== */

if (RUN_STAGE === 'B') {
  print('--- PART 9  Stage B export ---');

  if (INVENTORY_ASSET === null) {
    print('INVENTORY_ASSET is null. Run Stage A and export the asset first.');
  } else {
    var inventory = ee.FeatureCollection(INVENTORY_ASSET);

    // Empty SELECTED_EVENTS means every event in the asset.
    var selected = SELECTED_EVENTS.length === 0
      ? inventory
      : inventory.filter(ee.Filter.inList('event_id', SELECTED_EVENTS));

    print('Events in inventory asset:', inventory.size());
    // Stale-asset guard. PERIMETER_ASSET only affects Stage A, so an
    // inventory exported before the perimeters were switched will keep
    // feeding Stage B the old events. This shows which source the
    // asset actually came from.
    print('Inventory source (should match PERIMETER_ASSET setting):',
          inventory.aggregate_histogram('source'));
    print('Events selected for sampling:', selected.size());
    print('Event IDs:', selected.aggregate_array('event_id'));

    var burnSamples = ee.FeatureCollection(
      selected.toList(MAX_EVENTS).map(sampleBurnEvent)
    ).flatten();

    var oilSamples = sampleOilSite();
    var allSamples = burnSamples.merge(oilSamples);

    // Row counts per class and per event. These force every event's
    // pre/post pair and 50-band stack to evaluate at once, which
    // exceeded the interactive memory limit at 58 events. They are off
    // by default; the same counts are trivially recoverable from the
    // exported CSV, which is where they should be read from anyway.
    //
    // Set PRINT_STAGE_B_DIAGNOSTICS true only with a small
    // MATCHED_EVENT_LIMIT, for debugging a handful of events.
    if (PRINT_STAGE_B_DIAGNOSTICS) {
      print('Rows per class:',
            allSamples.aggregate_histogram('class_name'));
      print('Rows per event:',
            allSamples.aggregate_histogram('event_id'));
    } else {
      print('Row counts not printed (memory limit). ' +
            'Start the oil_burn_samples export and read them from the CSV.');
    }

    Export.table.toDrive({
      collection: allSamples,
      description: 'oil_burn_samples',
      folder: EXPORT_FOLDER,
      fileFormat: 'CSV',
      selectors: ['class_name', 'event_id', 'event_date', 'event_end_date',
                  'pre_date', 'post_date', 'area_ha', 'pre_ndvi', 'agri_pct',
                  'slope_mean', 'soil_mode', 'province',
                  'pre_NDVI_px', 'slope_px', 'soil_px']
                 .concat(FEATURE_BANDS)
                 .concat(FEATURE_BANDS.map(function (b) { return 'd_' + b; }))
    });

    Map.centerObject(OIL_POINT, 11);
    Map.addLayer(selected, {color: 'red'}, 'Selected burn events');
    Map.addLayer(OIL_WINDOW, {color: 'blue'}, 'Narli window');
  }
}

/* =====================================================================
 * PART 10  PLUME DIGITISING VIEW
 * RUN_STAGE = 'OIL' displays the Narli site before and after the
 * rupture so the plume can be drawn by hand.
 *
 * How to use:
 *   1. Set RUN_STAGE = 'OIL' and run.
 *   2. Toggle between the 19 June (pre) and 20 June (post) SWIR
 *      composites. The plume is the dark feature that is present in the
 *      post image and absent in the pre image, next to the rupture
 *      point marker.
 *   3. Use the polygon tool in the top-left of the map to draw around
 *      it. Draw the visible extent only; do not extend into ground that
 *      looks unchanged between the two dates.
 *   4. In the Geometry Imports panel, rename the drawn layer to
 *      plumeGeom, then set OIL_PLUME_GEOM = plumeGeom above.
 *   5. Re-run with RUN_STAGE = 'OIL' to see the reported area. The
 *      earlier bi-temporal analysis of this event gave 5.64 ha, so a
 *      hand-drawn polygon in that neighbourhood is a good sign. Do not
 *      force it to match: the drawn extent is its own measurement.
 *   6. Set RUN_STAGE = 'B' and run.
 * ===================================================================== */

if (RUN_STAGE === 'OIL') {
  var oilPreImg  = addIndices(pickScene(OIL_WINDOW, OIL_EVENT_DATE, true));
  var oilPostImg = addIndices(pickScene(OIL_WINDOW, OIL_EVENT_DATE, false));

  print('--- PART 10  plume digitising view ---');
  print('pre scene:', oilPreImg.date().format('YYYY-MM-dd'));
  print('post scene:', oilPostImg.date().format('YYYY-MM-dd'));

  // SWIR composite: crude on soil is darkest in the SWIR, so B12/B8/B4
  // separates it from dry soil better than true colour does.
  var swirVis = {bands: ['B12', 'B8', 'B4'], min: 0.02, max: 0.45};
  var rgbVis  = {bands: ['B4', 'B3', 'B2'], min: 0.02, max: 0.30};

  Map.centerObject(OIL_POINT, 15);
  Map.addLayer(oilPreImg,  swirVis, 'PRE  SWIR (B12/B8/B4)');
  Map.addLayer(oilPostImg, swirVis, 'POST SWIR (B12/B8/B4)');
  Map.addLayer(oilPreImg,  rgbVis,  'PRE  true colour', false);
  Map.addLayer(oilPostImg, rgbVis,  'POST true colour', false);
  Map.addLayer(OIL_POINT, {color: 'red'}, 'Rupture point');

  if (OIL_PLUME_GEOM !== null) {
    var drawnArea = ee.Geometry(OIL_PLUME_GEOM).area(1).divide(1e4);
    print('Drawn plume area (ha):', drawnArea);
    Map.addLayer(OIL_PLUME_GEOM, {color: 'yellow'}, 'Drawn plume');
  } else {
    print('No OIL_PLUME_GEOM set yet. Draw the polygon, rename the ' +
          'geometry import to plumeGeom, then set OIL_PLUME_GEOM.');
  }
}

/* =====================================================================
 * FIRST-RUN CHECKS
 * This script has not been executed. Verify, in order:
 *
 * 1. PART 2 prints an elevation near 566.6 m and a slope near 0.69 deg.
 *    If not, the DEM projection default is wrong or the point is
 *    transposed. This is the same trap as before: Point takes [lon, lat].
 * 2. GAUL province matches is 3, not 0. Turkish diacritics in
 *    ADM1_NAME may not match; if 0, the bbox fallback is in use and the
 *    region is approximate. Print distinct ADM1_NAME values to fix.
 * 3. Pre-event scene date resolves to 2018-06-15 or near it. If empty,
 *    lower CLEAR_MIN or widen PAIR_WINDOW_D.
 * 4. Cloud Score+ coverage: it does not extend to the earliest S2
 *    archive everywhere. If oilPre comes back empty, switch to the
 *    QA60 or SCL mask for 2018 and state the inconsistency in methods.
 * 5. MIRBI and BAIS2 assume reflectance as a fraction. Reflectance is
 *    divided by 10000 in s2Clear; confirm MIRBI values are order 1-10,
 *    not 1e4.
 * 6. Rejection counts in PART 5. If nearly everything fails on soil,
 *    OpenLandMap at 250 m is too coarse for this basin and soil
 *    matching should be dropped rather than loosened, with that
 *    decision stated.
 * 7. If matched burns are dominated by cropland stubble burning rather
 *    than wildfire, that is a good outcome, not a problem: char on
 *    agricultural soil is the closest analogue to oil on bare ground.
 *    Record the land use of matched events so it can be reported.
 * ===================================================================== */
