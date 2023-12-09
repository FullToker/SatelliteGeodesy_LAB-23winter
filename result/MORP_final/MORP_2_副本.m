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

%8
load('crust_ICE4G.mat');
load('crust_ICE5G.mat');

figure('Name','Question8')
subplot '211'
hold on
imagesc(0.5:359.5,-89.5:89.5,flipud(crust_ICE4G))
plot(lam,phi,'k');
plot(wrapTo360(ref_lam),ref_phi,'Marker','o','MarkerFaceColor','w')
title 'ICE 4G'
colormap 'jet'
colorbar
clim([-5 15])
axis equal
axis tight

subplot '212'
hold on
imagesc(0.5:359.5,-89.5:89.5,flipud(crust_ICE5G))
plot(lam,phi,'k');
plot(wrapTo360(ref_lam),ref_phi,'Marker','o','MarkerFaceColor','w')
title 'ICE 5G'
colormap 'jet'
colorbar
clim([-5 15])
axis equal
axis tight

%9
x=0.5:359.5;
y=89.5:-1:-89.5;
[X,Y]=meshgrid(x,y);

up_ice4=interp2(X,Y,crust_ICE4G,wrapTo360(ref_lam),ref_phi);
up_ice5=interp2(X,Y,crust_ICE5G,wrapTo360(ref_lam),ref_phi);

%10
up_gps=para_u1(2)*1000;