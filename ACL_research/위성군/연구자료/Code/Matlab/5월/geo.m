wgs84 = wgs84Ellipsoid("km");

lat0 = 37.57;
lon0 = 127
r = 50;
[lat, lon] = scircle1(lat0, lon0, r, [], wgs84);

geoplot(lat,lon,"k", "Linewidth", 2)
geobasemap streets