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

ref_lam=atan2d(ref(2),ref(1));%Calculating Longitude
ref_phi=atan2d(ref(3),sqrt(ref(1)^2+ref(2)^2));%Calculating Latitude

R2=[cosd(-ref_phi),0,-sind(-ref_phi);%Calculate the rotation matrix
    0,1,0;
    sind(-ref_phi),0,cosd(-ref_phi)];
R3=[cosd(ref_lam),sind(ref_lam),0;
    -sind(ref_lam),cosd(ref_lam),0;
    0,0,1];

data_uen=zeros(size(data,1),3);

for i=1:size(data,1)
    data_uen(i,:)=R2*R3*(data(i,5:7)'-ref');%Transform to the local horizontal system
end

ref_epoch=[1858,11,17,00,00,01];
t=decyear(data(:,2)+datenum(ref_epoch));

figure('Name','Series in the local horizontal system')
subplot '311'
plot(t,data_uen(:,1)')
title 'MORP-up'
grid minor

subplot '312'
plot(t,data_uen(:,2))
title 'MORP-east'
grid minor

subplot '313'
plot(t,data_uen(:,3))
title 'MORP-north'
grid minor

%2
%In order to get the linear trend and the amplitude of a seasonal signal：
A=[ones(size(data,1),1),t,cos(2*pi*t),sin(2*pi*t)];
%The first column is all 1, used to calculate the intercept
%The second column is time t, used to calculate the linear trend
%The third and fourth columns are trigonometric functions of time that capture cyclical changes

para_u=(A'*A)\(A'*data_uen(:,1));%least squares calculation
para_e=(A'*A)\(A'*data_uen(:,2));
para_n=(A'*A)\(A'*data_uen(:,3));

up_linear=para_u(2)*1000;%Calculate linear trend
east_linear=para_e(2)*1000;
north_linear=para_n(2)*1000;

up_annual=sqrt(para_u(3)^2+para_u(4)^2)*1000;%Calculate amplitude
east_annual=sqrt(para_e(3)^2+para_e(4)^2)*1000;
north_annual=sqrt(para_n(3)^2+para_n(4)^2)*1000;

%3
%Display geographic data with linear trends removed
figure('Name','Question 3')
subplot '311'
plot(t,data_uen(:,1)-para_u(1)-para_u(2)*t)
title 'MORP-up'
grid minor

subplot '312'
plot(t,data_uen(:,2)-para_e(1)-para_e(2)*t)
title 'MORP-east'
grid minor

subplot '313'
plot(t,data_uen(:,3)-para_n(1)-para_n(2)*t)
title 'MORP-north'
grid minor

%4

ending=2.010037899100203e+03;
t1=(t<=ending);
A1=[ones(size(data(t1),1),1),t(t1),cos(2*pi*t(t1)),sin(2*pi*t(t1))];
%Calculate parameters for each direction (up, east, north) using least squares method
para_u1=(A1'*A1)\(A1'*data_uen(t1,1));
para_e1=(A1'*A1)\(A1'*data_uen(t1,2));
para_n1=(A1'*A1)\(A1'*data_uen(t1,3));

%Plot data with linear and cyclical trends removed
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
%Convert the read data (all columns except the first column of strings) into a numeric matrix.

omega=nuvel(11,:);%Compute the horizontal movement rates
velocity=cross(omega,ref);
vel_nuvel=R2*R3*velocity';

%6

vel_gps=[para_u1(2);para_e1(2);para_n1(2)];%Calculate GPS speed

load 'coast30.mat'%Load map data

figure('Name','Question6')
hold on
plot(lam,phi,'k')
plot(wrapTo360(ref_lam),ref_phi, 'Marker','o', 'MarkerFaceColor','r')
q1=quiver(wrapTo360(ref_lam),ref_phi,vel_gps(2)*1e3,vel_gps(3)*1e3,'r');%Add velocity vector using 'quiver'，Red indicates GPS model
q2=quiver(wrapTo360(ref_lam),ref_phi,vel_nuvel(2)/1e3,vel_nuvel(3)/1e3,'b');%Add velocity vector using 'quiver'，Blue indicates nuvel model
legend([q1;q2],{'GPS','NUEVL'},'Location','southeast')
axis equal
axis tight

% 8
load('crust_ICE4G.mat');
load('crust_ICE5G.mat');

% Perform interpolation
% Assuming the original data has a resolution of 1 degree, we need to define a new grid to increase the resolution to 0.1 degrees
x_orig = 0.5:359.5;
y_orig = -89.5:89.5;
[X_orig, Y_orig] = meshgrid(x_orig, y_orig); % Original grid

x_new = 0.5:0.1:359.5;
y_new = 89.5:-0.1:-89.5;
[X_new, Y_new] = meshgrid(x_new, y_new); % New grid

crust_ICE4G_high_res = interp2(X_orig, Y_orig, flipud(crust_ICE4G), X_new, Y_new, 'linear');
crust_ICE5G_high_res = interp2(X_orig, Y_orig, flipud(crust_ICE5G), X_new, Y_new, 'linear');

% Define x and y for plotting
x = linspace(0.5, 360 - 0.5 / 3591, 3591);
y = linspace(89.5, -89.5 + 0.1 / 1791, 1791);

% Plot the interpolated data
figure('Name', 'High-Resolution ICE 4G and 5G')
subplot '211'
hold on
imagesc(x, y, crust_ICE4G_high_res); % Use defined x and y axis data
plot(lam, phi, 'k');
plot(wrapTo360(ref_lam), ref_phi, 'Marker', 'o', 'MarkerFaceColor', 'w')
title 'High-Res ICE 4G'
xlabel('Longitude')
ylabel('Latitude')
colormap 'jet'
colorbar
caxis([-5 15]) % Use caxis instead of clim
axis equal
axis tight

subplot '212'
hold on
imagesc(x, y, crust_ICE5G_high_res) % Use defined x and y axis data
plot(lam, phi, 'k');
plot(wrapTo360(ref_lam), ref_phi, 'Marker', 'o', 'MarkerFaceColor', 'w')
title 'High-Res ICE 5G'
xlabel('Longitude')
ylabel('Latitude')
colormap 'jet'
colorbar
caxis([-5 15]) % Use caxis instead of clim
axis equal
axis tight

% 9
% Interpolate latitude and longitude to obtain data values at locations
up_ice4_high_res = interp2(X_new, Y_new, crust_ICE4G_high_res, wrapTo360(ref_lam), ref_phi);
up_ice5_high_res = interp2(X_new, Y_new, crust_ICE5G_high_res, wrapTo360(ref_lam), ref_phi);


%10
up_gps=para_u1(2)*1000;