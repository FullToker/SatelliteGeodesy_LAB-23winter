clear all
close all
clc 


%Tasks 1/3
%1

fid=fopen('REYK_ig1.xyz');%Read the file
data=textscan(fid,'%s %s %f %f %f %f %f %f %f %f %f %f');%Acording to the file
fclose(fid);

data=cell2mat(data(:,3:end));%Convert third to last columns to numeric matrix
ref=[2587384.328,-1043033.510,5716564.045];%Read the Geocentric cartesian coordinates

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
xlabel('Time (year)')
ylabel(" (m)")
title 'REYK-up'
grid minor

subplot '312'
plot(t,data_uen(:,2))
xlabel('Time (year)')
ylabel(" (m)")
title 'REYK-east'
grid minor

subplot '313'
plot(t,data_uen(:,3))
xlabel('Time (year)')
ylabel(" (m)")
title 'REYK-north'
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
xlabel('Time (year)')
ylabel(" (m)")
title 'REYK-up'
grid minor

subplot '312'
plot(t,data_uen(:,2)-para_e(1)-para_e(2)*t)
xlabel('Time (year)')
ylabel(" (m)")
title 'REYK-east'
grid minor

subplot '313'
plot(t,data_uen(:,3)-para_n(1)-para_n(2)*t)
xlabel('Time (year)')
ylabel(" (m)")
title 'REYK-north'
grid minor

%4

jump1=2000.463;
jump2=2003.320;
jump3=2008.200;

t1=(t<=jump1);
A1=[ones(size(data(t1),1),1),t(t1),cos(2*pi*t(t1)),sin(2*pi*t(t1))];

para_u1=(A1'*A1)\(A1'*data_uen(t1,1));
para_e1=(A1'*A1)\(A1'*data_uen(t1,2));
para_n1=(A1'*A1)\(A1'*data_uen(t1,3));

t2=((jump1<t)&(t<=jump2));
A2=[ones(size(data(t2),1),1),t(t2),cos(2*pi*t(t2)),sin(2*pi*t(t2))];

para_u2=(A2'*A2)\(A2'*data_uen(t2,1));
para_e2=(A2'*A2)\(A2'*data_uen(t2,2));
para_n2=(A2'*A2)\(A2'*data_uen(t2,3));

t3=((jump2<t)&(t<=jump3));
A3=[ones(size(data(t3),1),1),t(t3),cos(2*pi*t(t3)),sin(2*pi*t(t3))];

para_u3=(A3'*A3)\(A3'*data_uen(t3,1));
para_e3=(A3'*A3)\(A3'*data_uen(t3,2));
para_n3=(A3'*A3)\(A3'*data_uen(t3,3));

t4=(t>jump3);
A4=[ones(size(data(t4),1),1),t(t4),cos(2*pi*t(t4)),sin(2*pi*t(t4))];

para_u4=(A4'*A4)\(A4'*data_uen(t4,1));
para_e4=(A4'*A4)\(A4'*data_uen(t4,2));
para_n4=(A4'*A4)\(A4'*data_uen(t4,3));


figure
subplot '311'
hold on
plot(t(t1),data_uen(t1,1)-para_u1(1)-para_u1(2)*t(t1) ...
    -para_u1(3)*cos(2*pi*t(t1))-para_u1(4)*sin(2*pi*t(t1)))
plot(t(t2),data_uen(t2,1)-para_u2(1)-para_u2(2)*t(t2) ...
    -para_u2(3)*cos(2*pi*t(t2))-para_u2(4)*sin(2*pi*t(t2)))
plot(t(t3),data_uen(t3,1)-para_u3(1)-para_u3(2)*t(t3) ...
    -para_u3(3)*cos(2*pi*t(t3))-para_u3(4)*sin(2*pi*t(t3)))
plot(t(t4),data_uen(t4,1)-para_u4(1)-para_u4(2)*t(t4) ...
    -para_u4(3)*cos(2*pi*t(t4))-para_u4(4)*sin(2*pi*t(t4)))
xlabel('Time (year)')
ylabel(" (m)")
title 'REYK-up'

subplot '312'
hold on
plot(t(t1),data_uen(t1,2)-para_e1(1)-para_e1(2)*t(t1) ...
    -para_e1(3)*cos(2*pi*t(t1))-para_e1(4)*sin(2*pi*t(t1)))
plot(t(t2),data_uen(t2,2)-para_e2(1)-para_e2(2)*t(t2) ...
    -para_e2(3)*cos(2*pi*t(t2))-para_e2(4)*sin(2*pi*t(t2)))
plot(t(t3),data_uen(t3,2)-para_e3(1)-para_e3(2)*t(t3) ...
    -para_e3(3)*cos(2*pi*t(t3))-para_e3(4)*sin(2*pi*t(t3)))
plot(t(t4),data_uen(t4,2)-para_e4(1)-para_e4(2)*t(t4) ...
    -para_e4(3)*cos(2*pi*t(t4))-para_e4(4)*sin(2*pi*t(t4)))
xlabel('Time (year)')
ylabel(" (m)")

title 'REYK-east'

subplot '313'
hold on
plot(t(t1),data_uen(t1,3)-para_n1(1)-para_n1(2)*t(t1) ...
    -para_n1(3)*cos(2*pi*t(t1))-para_n1(4)*sin(2*pi*t(t1)))
plot(t(t2),data_uen(t2,3)-para_n2(1)-para_n2(2)*t(t2) ...
    -para_n2(3)*cos(2*pi*t(t2))-para_n2(4)*sin(2*pi*t(t2)))
plot(t(t3),data_uen(t3,3)-para_n3(1)-para_n3(2)*t(t3) ...
    -para_n3(3)*cos(2*pi*t(t3))-para_n3(4)*sin(2*pi*t(t3)))
plot(t(t4),data_uen(t4,3)-para_n4(1)-para_n4(2)*t(t4) ...
    -para_n4(3)*cos(2*pi*t(t4))-para_n4(4)*sin(2*pi*t(t4)))
xlabel('Time (year)')
ylabel(" (m)")

title 'REYK-north'

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
vel_gps=[(para_u1(2)+para_u2(2)+para_u3(2)+para_u4(2))/4;(para_e1(2) ...
    +para_e2(2)+para_e3(2)+para_e4(2))/4;(para_n1(2)+para_n2(2)+para_n3(2)+para_n4(2))/4];

load 'coast30.mat'

figure('Name','Question6')
hold on
plot(lam,phi,'k')
plot(wrapTo360(ref_lam),ref_phi, 'Marker','o', 'MarkerFaceColor','r')
q1=quiver(wrapTo360(ref_lam),ref_phi,vel_gps(2)*1e3,vel_gps(3)*1e3,'r');
q2=quiver(wrapTo360(ref_lam),ref_phi,vel_nuvel(2)/1e3,vel_nuvel(3)/1e3,'b');
legend([q1;q2],{'GPS','NUEVL'},'Location','southeast')
xlabel('Latitude')
ylabel( 'Longitude')
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
up_gps=(para_u1(2)+para_u2(2)+para_u3(2)+para_u4(2))/4*1000;

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
