clear all
clc

load ./data/lab2_data.mat

aohi_ewh=aohi_ewh; %mm
aohi_def=aohi_def *1000; % mm

lat=size(aohi_ewh,1);
long=size(aohi_ewh,2);

%%% task1
%ewh linear
para_ewh=zeros(lat,long,4);
for lam=1:lat
    for phi=1:long
        data=squeeze(aohi_ewh(lam,phi,:));
        A=[ones(size(t)), t/365,cos(2*pi*t/365), sin(2*pi*t/365)];
        para_ewh(lam, phi, :)=(A' * A)\(A'*data);
    end
end

linear_ewh=para_ewh(:, :, 2);
annual_ewh=sqrt(para_ewh(:, :, 3).^2 + para_ewh(:, :, 4).^2);

%def linear
para_def=zeros(lat,long,4);
for lam=1:lat
    for phi=1:long
        data=squeeze(aohi_def(lam,phi, : ));
        A=[ones(size(t)), t/365,cos(2*pi*t/365), sin(2*pi*t/365)];
        para_def(lam, phi, :)=(A' * A)\(A'*data);
    end
end

linear_def=para_def(:, :, 2);
annual_def=sqrt(para_def(:, :, 3).^2 + para_def(:, :, 4).^2);

%%% figure
load ./data/coast30.mat

%aohi_ewh
figure
hold on
imagesc(0.5:359.5, 89.5:-1:-89.5, linear_ewh);
colormap 'jet'
colorbar
plot(lam,phi,'k')
title 'Linear trend of global surface densities (mm ewh/year)'
xlabel 'Longitude'
ylabel 'Latitude'
axis equal
axis tight

figure
hold on
imagesc(0.5:359.5, 89.5:-1:-89.5, annual_ewh);
colormap 'jet'
colorbar
plot(lam,phi,'k')
title 'Amplitude of annual global surface densities (mm ewh)'
xlabel 'Longitude'
ylabel 'Latitude'
axis equal
axis tight

%aohi_def
figure
hold on
imagesc(0.5:359.5, 89.5:-1:-89.5, linear_def);
colormap 'jet'
colorbar
plot(lam,phi,'k')
title 'Linear trend of global vertical crustal deformations (mm/year)'
xlabel 'Longitude'
ylabel 'Latitude'
axis equal
axis tight

figure
hold on
imagesc(0.5:359.5, 89.5:-1:-89.5, annual_def);
colormap 'jet'
colorbar
plot(lam,phi,'k')
title 'Amplitude of annual global vertical crustal deformations(mm)'
xlabel 'Longitude'
ylabel 'Latitude'
axis equal
axis tight


%%% task2 station time series

ref_kiru=[2251420.790, 862817.219, 5885476.705];
ref_morp=[3645667.836,-107277.235,5215053.530];
ref_reyk=[2587384.328,-1043033.510,5716564.045];

[ref_lam, ref_phi]=ref2ll(ref_kiru);

x=0.5: 359.5;
y= 89.5 : -1 : -89.5;
[X,Y]=meshgrid(x,y);

ref_ewh=zeros(size(t));
for i = 1: size(t,1)
    ref_ewh(i)=interp2(X, Y, aohi_ewh(:,:,i), wrapTo360(ref_lam), ref_phi);
end
ref_def=zeros(size(t));
for i = 1: size(t,1)
    ref_def(i)=interp2(X, Y, aohi_def(:,:,i), wrapTo360(ref_lam), ref_phi);
end
%visualize
figure
yyaxis left
plot(2003 + (t-t(1))/365, ref_ewh);
ylabel 'mm ewh';

yyaxis right
plot(2003 + (t-t(1))/365, ref_def);
ylabel 'mm';

xlabel 'year'
title 'EWH and DEF at station'
legend('EWH', 'DEF');
grid minor
axis tight

%%%Task3



