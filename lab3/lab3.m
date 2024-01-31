%lab 3

clear all
clc

load 'lab3_data.mat'

% Tasks 1/2: Analysis of global gravity variations based on monthly GRACE solutions

% 1. Chose a suitable time interval (or points in time) where both GRACE %
% and HI data are available. Use this period to estimate a linear term 
% as well as an annual signal for every grid point (cf. Lab 1 and Lab 2).

% GRACE
start_GRACE=1;

end_GRACE=44;
new_t_GRACE = t_GRACE(start_GRACE:end_GRACE)/365;

para_GRACE=zeros(180,360,4);

for lat = 1:100
    for long = 1:360
        data=squeeze(ewh_GRACE(lat,long,start_GRACE:end_GRACE));
        A=[ones(size(new_t_GRACE)),new_t_GRACE,cos(2*pi * new_t_GRACE),sin(2*pi*new_t_GRACE)];
        para_GRACE(lat,long, : )=(A'*A) \ ( A'* data);
    end
end

linear_GRACE=para_GRACE(:,:,2);
annual_GRACE=sqrt(para_GRACE(:,:,3).^2 + para_GRACE(:,:,4).^2);



%HI
start_HI=3;
end_HI=48;
new_t_HI=t_HI(start_HI:end_HI)/365;

para_HI=zeros(180,360,4);

for lat = 1:100
    for long = 1:360
        data=squeeze(ewh_HI(lat,long,start_HI:end_HI));
        A=[ones(size(new_t_HI)),new_t_HI,cos(2*pi * new_t_HI),sin(2*pi*new_t_HI)];
        para_HI(lat,long, : )=(A'*A)\(A'*data);
    end
end

linear_HI=para_HI(:,:,2);
annual_HI=sqrt(para_HI(:,:,3).^2 + para_HI(:,:,4).^2);



% 2. Visualize the linear trends and total amplitudes of the GRACE and
% HI data (per grid point) and give a brief interpretation of what you see.

load 'coast30.mat'

% Histogram of Amplitude of annual data
figure
histogram(annual_GRACE)
title 'Histogram of GRACE annual trend'

figure
histogram(annual_HI)
title 'Histogram of HI annual trend'


% Amplitude of annual GRACE data
figure

hold on
imagesc(0.5:359.5,89.5:-1:-89.5,annual_GRACE);
colormap 'jet'
colorbar
caxis([0 200])
plot(lam, phi,'k')

title 'Amplitude of annual GRACE data(mm ewh)'
xlabel 'Longitude(º)'
ylabel 'Latitude(º)'

axis equal
axis tight


% Amplitude of annual HI data
figure

hold on
imagesc(0.5:359.5,89.5:-1:-89.5,annual_HI);
colormap 'jet'
colorbar
caxis([0 200])

plot(lam, phi,'k')

title 'Amplitude of annual HI data(mm ewh)'
xlabel 'Longitude(º)'
ylabel 'Latitude(º)'

axis equal
axis tight


% Histogram of linear trands
figure
histogram(linear_GRACE)
title 'Histogram of GRACE linear trend'

figure
histogram(linear_HI)
title 'Histogram of HI linear trend'

% Linear trend of GRACE data
figure

hold on
imagesc(0.5:359.5,89.5:-1:-89.5, linear_GRACE);
colormap 'jet'
colorbar
caxis([-60 60])

plot(lam, phi,'k')

title 'Linear trend of GRACE data(mm ewh/year)'
xlabel 'Longitude(º)'
ylabel 'Latitude(º)'

axis equal
axis tight


%Linear trend of HI data
figure

hold on
imagesc(0.5:359.5,89.5:-1:-89.5, linear_HI);
colormap 'jet'
colorbar
caxis([-60 60])

plot(lam, phi,'k')

title 'Linear trend of HI data(mm ewh/year)'
xlabel 'Longitude(º)'
ylabel 'Latitude(º)'

axis equal
axis tight





% 3. Carry out a suitable validation of this statement and discuss your 
% results. The comparison can be done by an interpretation of the global 
% maps (trends + annual amplitudes) and their differences.
% Amplitude of annual GRACE data

% Linear trend of GIA data
figure

hold on
imagesc(0.5:359.5,89.5:-1:-89.5, ewh_ICE5G);
colormap 'jet'
colorbar
caxis([-60 60])

plot(lam,phi,'k')

title 'Linear trend GIA data(mm ewh/year)'
xlabel 'Longitude(º)'
ylabel 'Latitude'

axis equal
axis tight


% Linear trend of HI + GIA data
figure

hold on
imagesc(0.5:359.5,89.5:-1:-89.5, linear_HI + ewh_ICE5G);
colormap 'jet'
colorbar
caxis([-60 60])

plot(lam,phi,'k')

title 'Linear trend HI + GIA data(mm ewh/year)'
xlabel 'Longitude(º)'
ylabel 'Latitude'

axis equal
axis tight


% Tasks 2/2: Calculation of mass changes within the Greenland ice sheet
figure

hold on
imagesc(0.5:359.5,89.5:-1:-89.5, greenlandmask);
plot(lam, phi,'k')

title 'Greenland Mask'
xlabel 'Longitude(º)'
ylabel 'Latitude'

axis equal
axis tight

% 1. Use the GRACE data to calculate integrated mass variations over 
% Greenland with the eiven mask for all 108 months. For this purpose, 
% adjust the GIA effect with the ICE-56 model. Repeat this task for the HI data

% Integrated mass variation of GRACE data
sum_GRACE=zeros(size(t_GRACE));
for lat=1:180
    for long =1:360
        sigma=squeeze(ewh_GRACE(lat,long,:));% kg/m2
        sigma=sigma-ewh_ICE5G(lat,long)*(t_GRACE-t_GRACE(1))/365;% kg/m2
        sigma=greenlandmask(lat,long)*sigma*111195*111195*cosd(90-lat);% kg
        sigma=sigma/1e12;
        sum_GRACE=sum_GRACE + sigma;
    end
end

% Integrated mass variation of HI data
sum_HI=zeros(size(t_HI));
for lat=1:180
    for long =1:360
        sigma=squeeze(ewh_HI(lat,long,:));% kg/m2
        sigma=greenlandmask(lat,long)*sigma*111195*111195*cosd(90-lat);% kg
        sigma=sigma/1e12;
        sum_HI=sum_HI + sigma;
    end
end

% 2. Calculate the 11ncar trend [Gt/ycarJ and the yearly amplitude [Gt]
% of the estimated ice mass changes (both for GRACE and HI).

% Linear trend and yearly amplitude of integrated GRACE mass variation
t = (t_GRACE - t_GRACE(1))/365;
A = [ ones(size(t)),t,cos(2*pi * t),sin(2*pi * t)];
para_sum_GRACE = (A' * A)\(A' * sum_GRACE);

linear_sum_GRACE = para_sum_GRACE(2)
annual_sum_GRACE = sqrt(para_sum_GRACE(3)^2 + para_sum_GRACE(4)^2)

adjusted_sum_GRACE = A * para_sum_GRACE;
reduced_sum_GRACE = sum_GRACE - para_sum_GRACE(1)-para_sum_GRACE(2)*t;



% Linear trend and yearly amplitude of integrated HI mass variation
t = (t_HI - t_HI(1))/365;
A = [ ones(size(t)),t,cos(2*pi * t),sin(2*pi * t)];
para_sum_HI = (A' * A)\(A' * sum_HI);

linear_sum_HI = para_sum_HI(2)
annual_sum_HI = sqrt(para_sum_HI(3)^2 + para_sum_HI(4)^2)

adjusted_sum_HI = A * para_sum_HI;
reduced_sum_HI = sum_HI - para_sum_HI(1)-para_sum_HI(2)*t;



% 3, Visualize the results of the mass integration over Greenland for GRACE 
% and the model data in one figure. Discuss your results with respect to 
% - amongst others - similarities and differences and why they occur.

date_GRACE = datetime(t_GRACE, 'ConvertFrom','datenum');
date_HI = datetime(t_HI, 'ConvertFrom','datenum');

% Visualization of integrated mass variations of GRACE and HI data
figure

hold on
plot(date_GRACE,sum_GRACE)
plot(date_HI,sum_HI)

title 'Integrated mass variations of GRACE and HI data'
xlabel 'year'
ylabel 'Gt'
legend 'GRACE' 'HI'
axis tight

% visualization of adjusted mass variations of GRACE and HI data
figure

hold on
plot(date_GRACE,adjusted_sum_GRACE)
plot(date_HI,adjusted_sum_HI)

title 'Adjusted mass variations of GRACE and HI data'
xlabel 'year'
ylabel 'Gt'
legend 'GRACE' 'HI'
axis tight

% visualization of reduced mass variations ef GRACE and HI data
figure

hold on
plot(date_GRACE(start_GRACE:end_GRACE),reduced_sum_GRACE(start_GRACE:end_GRACE))
plot(date_HI(start_HI:end_HI),reduced_sum_HI(start_HI:end_HI))

title 'Reduced mass variations of GRACE and HI data'
xlabel 'year'
ylabel 'Gt'
legend 'GRACE' 'HI'
axis tight












