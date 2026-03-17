wgs84 = wgs84Ellipsoid("km");

lat0 = 37.5642135;
lon0 = 127.0016985;

r = 25;
[lat,lon] = scircle1(lat0,lon0,r,[],wgs84);

h = geoplot(lat, lon, ".", "LineWidth", 2)
geobasemap streets
