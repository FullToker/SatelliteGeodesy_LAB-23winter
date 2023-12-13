clear all
close all
clc 


%Tasks 1/3
%1

fid=fopen('MORP_ig1.xyz');%Read the file
data=textscan(fid,'%s %s %f %f %f %f %f %f %f %f %f %f');%Acording to the file
fclose(fid);

data=cell2mat(data(:,3:end));%Convert third to last columns to numeric matrix
ref=[3645667.836,-107277.235,5215053.530];%Read the Geocentric cartesian coordinates

ref_lam=atan2d(ref(2),ref(1));
ref_phi=atan2d(ref(3),sqrt(ref(1)^2+ref(2)^2));

R2=[cosd(-ref_phi),0,-sind(-ref_phi);
    0,1,0;
    sind(-ref_phi),0,cosd(-ref_phi)];
R3=[cosd(ref_lam),sind(ref_lam),0;
    -sind(ref_lam),cosd(ref_lam),0;
    0,0,1];

data_uen=zeros(size(data,1),3);

for i=1:size(data,1)
    data_uen(i,:)=R2*R3*(data(i,5:7)'-ref');
end

ref_epoch=[1858,11,17,00,00,01];
t=decyear(data(:,2)+datenum(ref_epoch));

figure('Name','Series in the local horizontal system')
subplot '311'
plot(t,data_uen(:,1)')
title 'KIRU-up'
grid minor

subplot '312'
plot(t,data_uen(:,2))
title 'KIRU-east'
grid minor

subplot '313'
plot(t,data_uen(:,3))
title 'KIRU-north'
grid minor

%2

A=[ones(size(data,1),1),t,cos(2*pi*t),sin(2*pi*t)];

para_u=(A'*A)\(A'*data_uen(:,1));
para_e=(A'*A)\(A'*data_uen(:,2));
para_n=(A'*A)\(A'*data_uen(:,3));

up_linear=para_u(2)*1000;
east_linear=para_e(2)*1000;
north_linear=para_n(2)*1000;

up_annual=sqrt(para_u(3)^2+para_u(4)^2)*1000;
east_annual=sqrt(para_e(3)^2+para_e(4)^2)*1000;
north_annual=sqrt(para_n(3)^2+para_n(4)^2)*1000;

%3

figure('Name','Question 3')
subplot '311'
plot(t,data_uen(:,1)-para_u(1)-para_u(2)*t)
title 'KIRU-up'
grid minor

subplot '312'
plot(t,data_uen(:,2)-para_e(1)-para_e(2)*t)
title 'KIRU-east'
grid minor

subplot '313'
plot(t,data_uen(:,3)-para_n(1)-para_n(2)*t)
title 'KIRU-north'
grid minor

%4
ending=2.010037899100203e+03;
t1=(t<=ending);
A1=[ones(size(data(t1),1),1),t(t1),cos(2*pi*t(t1)),sin(2*pi*t(t1))];

para_u1=(A1'*A1)\(A1'*data_uen(t1,1));
para_e1=(A1'*A1)\(A1'*data_uen(t1,2));
para_n1=(A1'*A1)\(A1'*data_uen(t1,3));


figure
subplot '311'
hold on
plot(t(t1),data_uen(t1,1)-para_u1(1)-para_u1(2)*t(t1) ...
    -para_u1(3)*cos(2*pi*t(t1))-para_u1(4)*sin(2*pi*t(t1)))

title 'up'

subplot '312'
hold on
plot(t(t1),data_uen(t1,2)-para_e1(1)-para_e1(2)*t(t1) ...
    -para_e1(3)*cos(2*pi*t(t1))-para_e1(4)*sin(2*pi*t(t1)))

title 'east'

subplot '313'
hold on
plot(t(t1),data_uen(t1,3)-para_n1(1)-para_n1(2)*t(t1) ...
    -para_n1(3)*cos(2*pi*t(t1))-para_n1(4)*sin(2*pi*t(t1)))

title 'north'

%TASK2

%5

fid=fopen('NNR_NUVEL1A.txt');
nuvel=textscan(fid, '%s %f %f %f', 'headerlines',2);
fclose(fid);
nuvel=cell2mat(nuvel(:,2:end));

omega=nuvel(11,:);
velocity=cross(omega,ref);
vel_nuvel=R2*R3*velocity';

%6

vel_gps=[para_u1(2);para_e1(2);para_n1(2)];
load 'coast30.mat'

figure('Name','Question6')
hold on
plot(lam,phi,'k')
plot(wrapTo360(ref_lam),ref_phi, 'Marker','o', 'MarkerFaceColor','r')
q1=quiver(wrapTo360(ref_lam),ref_phi,vel_gps(2)*1e3,vel_gps(3)*1e3,'r');
q2=quiver(wrapTo360(ref_lam),ref_phi,vel_nuvel(2)/1e3,vel_nuvel(3)/1e3,'b');
legend([q1;q2],{'GPS','NUEVL'},'Location','southeast')
axis equal
axis tight

% 8
load('crust_ICE4G.mat');
load('crust_ICE5G.mat');

x_orig = 0.5:359.5;
y_orig = -89.5:89.5;
[X_orig, Y_orig] = meshgrid(x_orig, y_orig); 

x_new = 0.5:0.1:359.5;
y_new = 89.5:-0.1:-89.5;
[X_new, Y_new] = meshgrid(x_new, y_new); 

crust_ICE4G_high_res = interp2(X_orig, Y_orig, flipud(crust_ICE4G), X_new, Y_new, 'linear');
crust_ICE5G_high_res = interp2(X_orig, Y_orig, flipud(crust_ICE5G), X_new, Y_new, 'linear');

x = linspace(0.5, 360 - 0.5 / 3591, 3591);
y = linspace(89.5, -89.5 + 0.1 / 1791, 1791);

figure('Name', 'High-Resolution ICE 4G and 5G')
subplot '211'
hold on
imagesc(x, y, crust_ICE4G_high_res); 
plot(lam, phi, 'k');
plot(wrapTo360(ref_lam), ref_phi, 'Marker', 'o', 'MarkerFaceColor', 'w')
title 'High-Res ICE 4G'
xlabel('Longitude')
ylabel('Latitude')
colormap 'jet'
colorbar
caxis([-5 15]) 
axis equal
axis tight

subplot '212'
hold on
imagesc(x, y, crust_ICE5G_high_res) 
plot(lam, phi, 'k');
plot(wrapTo360(ref_lam), ref_phi, 'Marker', 'o', 'MarkerFaceColor', 'w')
title 'High-Res ICE 5G'
xlabel('Longitude')
ylabel('Latitude')
colormap 'jet'
colorbar
caxis([-5 15]) 
axis equal
axis tight

% 9
up_ice4_high_res = interp2(X_new, Y_new, crust_ICE4G_high_res, wrapTo360(ref_lam), ref_phi);
up_ice5_high_res = interp2(X_new, Y_new, crust_ICE5G_high_res, wrapTo360(ref_lam), ref_phi);



%10
up_gps=para_u1(2)*1000;

%11 Site description
coords = [
    2251420.790, 862817.219, 5885476.705; 
    2587384.310, -1043033.522, 5716564.039; 
    3645667.836, -107277.235, 5215053.530  
];

names = {'KIRU', 'REYK', 'MORP'};

wgs84 = wgs84Ellipsoid('meter');
[lat, lon, h] = ecef2geodetic(wgs84, coords(:,1), coords(:,2), coords(:,3));

figure;
axesm('MapProjection', 'ortho', 'Origin', [90 0 0], 'MapLatLimit', [0 90]);

land = shaperead('landareas.shp', 'UseGeoCoords', true);
geoshow([land.Lat], [land.Lon], 'DisplayType', 'polygon', 'FaceColor', [0.5 0.7 0.5]);

plotm(lat, lon, 'r.', 'MarkerSize', 20);

for i = 1:length(lat)
    textm(lat(i), lon(i), names{i}, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
end

setm(gca, 'MLabelParallel', 0, 'MLabelLocation', 60);
setm(gca, 'PLabelMeridian', 0, 'PLabelLocation', 30);
plabel on;
mlabel on;
gridm on;
framem on;

axis off;

saveas(gcf, 'north_pole_stations_with_land.png');
