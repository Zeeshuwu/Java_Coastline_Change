var ROI = ROI8;

function maskL457(image) {
  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 3).eq(0).and(qa.bitwiseAnd(1 << 4).eq(0));
  return image.updateMask(mask);
}

function maskL89(image) {
  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 3).eq(0).and(qa.bitwiseAnd(1 << 4).eq(0));
  return image.updateMask(mask);
}

function scaleImg(img) {
  return img.multiply(0.0000275).add(-0.2);
}

function extractCoastline(image, label) {
  var mndwi = image.normalizedDifference(['green', 'swir']).rename('MNDWI');
  if (label === '2025' || label === '2020' || label === '2015') {
    mndwi = mndwi.focal_median({ radius: 2, iterations: 1 });
  }
  var landMask = mndwi.lt(0.0);
  var minSize = 50;
  var cleanLand = landMask.updateMask(landMask.connectedPixelCount(minSize, true).gte(minSize));
  var kernel = ee.Kernel.circle({ radius: 1.5 });
  var smoothMask = cleanLand
    .focal_min({ kernel: kernel, iterations: 1 })
    .focal_max({ kernel: kernel, iterations: 1 })
    .focal_mode({ radius: 2, iterations: 1 });
  return ee.Algorithms.CannyEdgeDetector({
    image: smoothMask,
    threshold: 0.5,
    sigma: 1
  }).selfMask();
}

function getVectorLine(image, region) {
  var imageInt = image.gt(0).selfMask().toByte();
  var vectors = imageInt.reduceToVectors({
    geometry: region,
    scale: 30,
    geometryType: 'polygon',
    eightConnected: true,
    maxPixels: 1e13,
    labelProperty: 'b1'
  });
  return vectors.map(function(f) {
    var polygonCoords = f.geometry().coordinates();
    return f.setGeometry(ee.Geometry.LineString(polygonCoords.get(0)));
  });
}

var l4 = ee.ImageCollection("LANDSAT/LT04/C02/T1_L2");
var l5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2");
var l7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2");
var l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2");
var l9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2");

var col2002 = l5.filterDate('2002-01-01', '2002-12-31').filterBounds(ROI).map(maskL457)
  .merge(l7.filterDate('2002-01-01', '2002-12-31').filterBounds(ROI).map(maskL457));
var img2002 = col2002.map(function(img) {
  return img.select(['SR_B2', 'SR_B5'], ['green', 'swir']);
}).map(scaleImg).median().clip(ROI);
var coast2002 = extractCoastline(img2002, '2002');

var col2025 = l8.filterDate('2024-01-01', '2025-12-31').filterBounds(ROI).map(maskL89)
  .merge(l9.filterDate('2024-01-01', '2025-12-31').filterBounds(ROI).map(maskL89));
var img2025 = col2025.map(function(img) {
  return img.select(['SR_B3', 'SR_B6'], ['green', 'swir']);
}).map(scaleImg).median().clip(ROI);
var coast2025 = extractCoastline(img2025, '2025');

var mndwi2002 = img2002.normalizedDifference(['green', 'swir']).rename('MNDWI');
var mndwi2025 = img2025.normalizedDifference(['green', 'swir']).rename('MNDWI');

Map.centerObject(ROI, 13);
Map.addLayer(img2025.normalizedDifference(['green', 'swir']), { min: -0.5, max: 0.5, palette: ['red', 'white', 'blue'] }, 'Check MNDWI 2025', false);
Map.addLayer(coast2002, { palette: ['FF0000'], min: 0, max: 1, strokeWidth: 2 }, 'Garis Pantai 2002');
Map.addLayer(coast2025, { palette: ['FFFF00'], min: 0, max: 1, strokeWidth: 2 }, 'Garis Pantai 2025');

var panel = ui.Panel({ style: { position: 'bottom-left', padding: '8px' } });
panel.add(ui.Label('Garis Pantai', { fontWeight: 'bold' }));
panel.add(ui.Label('Merah: 2002', { color: 'red' }));
panel.add(ui.Label('Kuning: 2025', { color: '#d4af37' }));
Map.add(panel);

var shp2002_line = getVectorLine(coast2002, ROI);
var shp2025_line = getVectorLine(coast2025, ROI);

Export.table.toDrive({ collection: shp2002_line, description: 'Shoreline_2002_Line', folder: 'NEW_SHORELINES', fileFormat: 'SHP' });
Export.table.toDrive({ collection: shp2025_line, description: 'Shoreline_2025_Line', folder: 'NEW_SHORELINES', fileFormat: 'SHP' });
Export.image.toDrive({ image: mndwi2002, description: 'MNDWI_2002', folder: 'NEW_SHORELINES', scale: 30, region: ROI, maxPixels: 1e13 });
Export.image.toDrive({ image: mndwi2025, description: 'MNDWI_2025', folder: 'NEW_SHORELINES', scale: 30, region: ROI, maxPixels: 1e13 });

var years = [1990, 1995, 2005, 2010, 2015, 2020];
var colors = ['00FFFF', 'FF00FF', '00FF00', 'FFA500', '800080', 'FFC0CB'];

years.forEach(function(year, index) {
  var col, img;
  var start = year + '-01-01';
  var end = year + '-12-31';
  var label = year.toString();

  if (year < 1999) {
    col = l4.filterDate(start, end).filterBounds(ROI).map(maskL457)
      .merge(l5.filterDate(start, end).filterBounds(ROI).map(maskL457));
    img = col.map(function(i) { return i.select(['SR_B2', 'SR_B5'], ['green', 'swir']); });
  } else if (year < 2013) {
    col = l5.filterDate(start, end).filterBounds(ROI).map(maskL457)
      .merge(l7.filterDate(start, end).filterBounds(ROI).map(maskL457));
    img = col.map(function(i) { return i.select(['SR_B2', 'SR_B5'], ['green', 'swir']); });
  } else {
    col = l8.filterDate(start, end).filterBounds(ROI).map(maskL89);
    img = col.map(function(i) { return i.select(['SR_B3', 'SR_B6'], ['green', 'swir']); });
  }

  var composite = img.map(scaleImg).median().clip(ROI);
  var mndwi = composite.normalizedDifference(['green', 'swir']).rename('MNDWI');
  var coast = extractCoastline(composite, label);

  Map.addLayer(coast, { palette: [colors[index]], min: 0, max: 1, strokeWidth: 2 }, 'Garis Pantai ' + year, false);

  var shp_line = getVectorLine(coast, ROI);

  Export.table.toDrive({ collection: shp_line, description: 'Shoreline_' + year + '_Line', folder: 'NEW_SHORELINES', fileFormat: 'SHP' });
  Export.image.toDrive({ image: mndwi, description: 'MNDWI_' + year, folder: 'NEW_MNDWI', scale: 30, region: ROI, maxPixels: 1e13 });
});
