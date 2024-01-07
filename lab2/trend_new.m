% from GNSS
%REYK
fid=fopen('./data/REYK_ig1.xyz');
% ref_reyk=[2587384.328,-1043033.510,5716564.045];

data=textscan(fid, '%s %s %f %f %f %f %f %f %f %f %f %f');
fclose(fid);
data=cell2mat(data(:, 3:end));

% calculate rotation matrix
[ref_lam, ref_phi]=ref2ll(ref_reyk);

R2=[cosd(-ref_phi),0,-sind(-ref_phi);
    0,1,0;
    sind(-ref_phi),0,cosd(-ref_phi)];
R3=[cosd(ref_lam),sind(ref_lam),0;
    -sind(ref_lam),cosd(ref_lam),0;
    0,0,1];
data_uen=zeros(size(data,1),3);
for i=1:size(data,1)
    data_uen(i,:)=R2*R3*(data(i,5:7)'-ref_reyk') *1e3;
end

%compute the linear trend and amplitude of annual signal
ref_epoch=[1858,11,17,00,00,01];
t2=decyear(data(:,2) + datenum(ref_epoch));
data_uen=data_uen(t2>2003 & t2<2007, :);
t2=t2(t2>2003 & t2<2007);

A=[ones(size(t2)),t2,cos(2*pi * t2), sin(2*pi * t2)];
para_u=(A' * A)\(A' *  data_uen(:,1));
up_linear_reyk=para_u(2);
up_annual_reyk= sqrt(para_u(3)^2 + para_u(4)^2);

%%%AOHI + GIA vs GNSS
offset=0;
figure
hold on
plot(t2, data_uen(:,1));
plot(2003 + (t-t(1))/365, def_reyk + up_reyk(3)* (t-t(1))/365 + offset);
plot(2003 + (t-t(1))/365, def_reyk + up_reyk(4)* (t-t(1))/365 + offset);

title 'GNSS and DEF + GIA at station REYK'
legend ('GNSS', 'DEF+ICE4G', "DEF+ICE5G")
ylabel 'mm'
xlabel 'year'
grid minor


%KIRU
fid=fopen('./data/KIRU_ig1.xyz');

data=textscan(fid, '%s %s %f %f %f %f %f %f %f %f %f %f');
fclose(fid);
data=cell2mat(data(:, 3:end));

% calculate rotation matrix
[ref_lam, ref_phi]=ref2ll(ref_kiru);

R2=[cosd(-ref_phi),0,-sind(-ref_phi);
    0,1,0;
    sind(-ref_phi),0,cosd(-ref_phi)];
R3=[cosd(ref_lam),sind(ref_lam),0;
    -sind(ref_lam),cosd(ref_lam),0;
    0,0,1];
data_uen=zeros(size(data,1),3);
for i=1:size(data,1)
    data_uen(i,:)=R2*R3*(data(i,5:7)'-ref_kiru') *1e3;
end


t2=decyear(data(:,2) + datenum(ref_epoch));
data_uen=data_uen(t2>2003 & t2<2007, :);
t2=t2(t2>2003 & t2<2007);

A=[ones(size(t2)),t2,cos(2*pi * t2), sin(2*pi * t2)];
para_u=(A' * A)\(A' *  data_uen(:,1));
up_linear_kiru=para_u(2);
up_annual_kiru= sqrt(para_u(3)^2 + para_u(4)^2);

%%%AOHI + GIA vs GNSS
offset=0;
figure
hold on
plot(t2, data_uen(:,1));
plot(2003 + (t-t(1))/365, def_kiru + up_kiru(3)* (t-t(1))/365 + offset);
plot(2003 + (t-t(1))/365, def_kiru + up_kiru(4)* (t-t(1))/365 + offset);

title 'GNSS and DEF + GIA at station KIRU'
legend ('GNSS', 'DEF+ICE4G', "DEF+ICE5G")
ylabel 'mm'
xlabel 'year'
grid minor

%MORP
fid=fopen('./data/MORP_ig1.xyz');

data=textscan(fid, '%s %s %f %f %f %f %f %f %f %f %f %f');
fclose(fid);
data=cell2mat(data(:, 3:end));

% calculate rotation matrix
[ref_lam, ref_phi]=ref2ll(ref_morp);

R2=[cosd(-ref_phi),0,-sind(-ref_phi);
    0,1,0;
    sind(-ref_phi),0,cosd(-ref_phi)];
R3=[cosd(ref_lam),sind(ref_lam),0;
    -sind(ref_lam),cosd(ref_lam),0;
    0,0,1];
data_uen=zeros(size(data,1),3);
for i=1:size(data,1)
    data_uen(i,:)=R2*R3*(data(i,5:7)'-ref_morp') *1e3;
end


t2=decyear(data(:,2) + datenum(ref_epoch));
data_uen=data_uen(t2>2003 & t2<2007, :);
t2=t2(t2>2003 & t2<2007);

A=[ones(size(t2)),t2,cos(2*pi * t2), sin(2*pi * t2)];
para_u=(A' * A)\(A' *  data_uen(:,1));
up_linear_morp=para_u(2);
up_annual_morp= sqrt(para_u(3)^2 + para_u(4)^2);

%%%AOHI + GIA vs GNSS
offset=0;
figure
hold on
plot(t2, data_uen(:,1));
plot(2003 + (t-t(1))/365, def_morp + up_morp(3)* (t-t(1))/365 + offset);
plot(2003 + (t-t(1))/365, def_morp + up_morp(4)* (t-t(1))/365 + offset);

title 'GNSS and DEF + GIA at station MORP'
legend ('GNSS', 'DEF+ICE4G', "DEF+ICE5G")
ylabel 'mm'
xlabel 'year'
grid minor
