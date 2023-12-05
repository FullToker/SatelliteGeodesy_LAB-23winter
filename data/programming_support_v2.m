% GGOS Matlab tips (exercise 1)
% Hints für reading, transforming and visualizing the data in Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import coordinate time series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% There are numerous possibilities for the data import to MATLAB
% example: 

filename = 'FAIR_ig1.xyz';  % example: station in Fairbanks
[col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12] ...
      = textread(filename, '%s%s%n%n%n%n%n%n%n%n%n%n\n');
XYZ = [col7, col8, col9];
GPSWEEK = col5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time transformation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% It is recommended to do all calculations with serial date numbers (SDN)
% which is MATLAB's time format

% 1. BEST: Use MJD to calculate the SDN with "datenum"
% MJD refers to 17 Nov 1858
% SDN refers to 01.01.0000

% 2. MEDIOCRE alternative: Use the GPS weeks to calculate the SDN

% Beginning of GPS-week 873  
%                =  29.9.1996 0:00 GMT  
%                =  matlab serial date number 729297.
% Mid of GPS-week 873
%                =  matlab serial date number 729300.5
% Which results in following transformation between GPS-week and the
% "serial data number"
% serialdatenumber = 729300.5 + (gpsweek - 873) * 7;


% Following Matlab functions could be used for further steps.
% - datevec
% - datestr
% - datetick

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coordinate transformation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A useful set of MATLAB codes for geodetic coordinate transformations is called
% 'Geodetic transformation toolbox' and was distributed  and used in the GuG Bachelor. 
% You can ask your colleagues about it, 

% use MATLAB intern routines (see documentation for ecef2geodetic)

% or include your own transformation routines. 

%%%%%%%%%%%%%%%%%%%%%%
% Visualize ICE4G
%%%%%%%%%%%%%%%%%%%%%%

load('ICE4G')

imagesc(0.5:359.5,89.5:-1:-89.5, crust_ICE4G); 
axis xy
% "axis xy" changes the direction of the y-axis (Increasing values from
% bottom to top).
caxis([-5,20])  % define colorbar
colorbar
hold on

% In case you want to present ICE4G on a map from -180° to 180°
imagesc(-179.5:1:179.5,89.5:-1:-89.5, crust_ICE4G(:,[181:360,1:180]));


%%%%%%%%%%%%%%%%%%%%%%
% Coastlines
%%%%%%%%%%%%%%%%%%%%%%

load('coast30'); %-> lam& phi 

% Visualize coastlines
% longitude -180° to 180°
plot( mod(lam+180,360)-180, phi, 'Color', 0.6*[1,1,1]); 
hold on
% longitude 0° to 360°
plot( lam, phi, 'Color', 0.6*[1,1,1]); hold on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple map with arrows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example for an arrow presentation
% first two input values for "origin", with (East, North) = (0.02, 0.0)
quiver(0, 0, 1000* 0.02, 1000* 0, 'Color', 'black')
% first two input values for "origin", (East, North) = (0.0, 0.02)
quiver(0, 0, 1000* 0,1000*  0.02, 'Color', 'red')
% first two input values for "origin", (East, North) = (0.014, 0.014)
quiver(0, 0, 1000* 0.014, 1000* 0.014, 'Color', 'red')

% Multiplication with 1000 is chosen randomly but identical for all
% arrows. Changing it you can choose a useful size for your arrows.
% Combining these arrows with a legend shows their real values (e.g. plot
% an arrow with 0.01m/year as part of the legend).

hold off
