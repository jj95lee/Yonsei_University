 wgs84 = wgs84Ellipsoid("km");

lat0 = 37.5642135;
lon0 = 126.0016985;

r = 25;
[lat,lon] = scircle(lat0,lon0,r,[],wgs84);

geoplot(lat, lon, "k", "LineWidth", 2)
geobasemap streets
