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
%Define unlimited time points and split the data
jump1=2000.463;
jump2=2003.320;
jump3=2008.200;

t1=(t<=jump1);
%For the four sets of data specified by t1,t2,t3,t4, construct design matrices A1,A2,A3,A4 respectively.
A1=[ones(size(data(t1),1),1),t(t1),cos(2*pi*t(t1)),sin(2*pi*t(t1))];
%Calculate parameters for each direction (up, east, north) using least squares method
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

%Plot data with linear and cyclical trends removed
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
%Convert the read data (all columns except the first column of strings) into a numeric matrix.

omega=nuvel(11,:);%Compute the horizontal movement rates
velocity=cross(omega,ref);
vel_nuvel=R2*R3*velocity';

%6
vel_gps=[(para_u1(2)+para_u2(2)+para_u3(2)+para_u4(2))/4;(para_e1(2) ...
    +para_e2(2)+para_e3(2)+para_e4(2))/4;(para_n1(2)+para_n2(2)+para_n3(2)+para_n4(2))/4];%Calculate GPS speed

load 'coast30.mat'%Load map data

figure('Name','Question6')
hold on
plot(lam,phi,'k')
plot(wrapTo360(ref_lam),ref_phi, 'Marker','o', 'MarkerFaceColor','r')
q1=quiver(wrapTo360(ref_lam),ref_phi,vel_gps(2)*1e3,vel_gps(3)*1e3,'r');%Add velocity vector using 'quiver'，Red indicates GPS model
q2=quiver(wrapTo360(ref_lam),ref_phi,vel_nuvel(2)/1e3,vel_nuvel(3)/1e3,'b');%Add velocity vector using 'quiver'，Blue indicates nuvel model
legend([q1;q2],{'GPS','NUEVL'},'Location','southeast')
xlabel('Latitude')
ylabel( 'Longitude')
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
up_gps=(para_u1(2)+para_u2(2)+para_u3(2)+para_u4(2))/4*1000;