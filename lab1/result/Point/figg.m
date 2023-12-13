% 坐标
coords = [
    2251420.790, 862817.219, 5885476.705; % KIRU
    2587384.310, -1043033.522, 5716564.039; % REYK
    3645667.836, -107277.235, 5215053.530  % MORP
];

% 名称
names = {'KIRU', 'REYK', 'MORP'};

% ECEF坐标转换为地理坐标（纬度、经度）
wgs84 = wgs84Ellipsoid('meter');
[lat, lon, h] = ecef2geodetic(wgs84, coords(:,1), coords(:,2), coords(:,3));

% 创建一个地图并使用正射投影
figure;
axesm('MapProjection', 'ortho', 'Origin', [90 0 0], 'MapLatLimit', [0 90]);

% 显示陆地
land = shaperead('landareas.shp', 'UseGeoCoords', true);
geoshow([land.Lat], [land.Lon], 'DisplayType', 'polygon', 'FaceColor', [0.5 0.7 0.5]);

% 标注站点位置
plotm(lat, lon, 'r.', 'MarkerSize', 20);

% 在每个站点位置添加文本
for i = 1:length(lat)
    textm(lat(i), lon(i), names{i}, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
end

% 添加网格和标签
setm(gca, 'MLabelParallel', 0, 'MLabelLocation', 60);
setm(gca, 'PLabelMeridian', 0, 'PLabelLocation', 30);
plabel on;
mlabel on;
gridm on;
framem on;

% 隐藏坐标轴
axis off;

% 将图形保存为图片
saveas(gcf, 'north_pole_stations_with_land.png');
