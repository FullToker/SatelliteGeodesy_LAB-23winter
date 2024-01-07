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
        A=[ones(size(t)), t/365,cos(2 * pi * t/365), sin(2 * pi * t/365)];
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
xlabel 'Longitude(°)'
ylabel 'Latitude(°)'
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

[ewh_kiru, def_kiru]=cal_draw_sta(ref_kiru,  'KIRU',  aohi_ewh,aohi_def,t);
[ewh_morp, def_morp]=cal_draw_sta(ref_morp, 'MORP',  aohi_ewh,aohi_def,t);
[ewh_reyk, def_reyk]=cal_draw_sta(ref_reyk, 'REYK',   aohi_ewh,aohi_def,t);


%%%Task3

%from lab1 : gnss observations
% kiru =[up_linear, up_anaual, up_ice4, up_ice5] unit: mm/y, mm, mm/y mm/y
up_kiru=[7.3089, 5.1123, 5.5568, 6.1335]; 
up_morp=[0.4201, 3.9678, 0.0140, -0.0408]; 
up_reyk=[1.2973,4.1558, -0.0464, 0.7236];
% linear and annual def_ref
def_kiru_rate=lin_ann_def(def_kiru, t);
def_morp_rate=lin_ann_def(def_morp, t);
def_reyk_rate=lin_ann_def(def_reyk, t);

%error
% linear trend: AOHI + GIA = GNSS;
% annual amplitude: AOHI = GNSS
% err_station=[error_linear, error_ice4, error_ice5, error_annual]
err_kiru=cal_err(up_kiru,def_kiru_rate);
err_morp=cal_err(up_morp,def_morp_rate);
err_reyk=cal_err(up_reyk,def_reyk_rate);

%%% visualize
% AOHI +GIA vs GNSS




