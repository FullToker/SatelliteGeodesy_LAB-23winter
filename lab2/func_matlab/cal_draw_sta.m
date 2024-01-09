function [ref_ewh,ref_def] = cal_draw_sta(ref,name, aohi_ewh,aohi_def,t)
%   calculate the ref_ewh and ref_def according to the station ref
[ref_lam, ref_phi]=ref2ll(ref);

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
title(['EWH and DEF at station ' name]);
legend('EWH', 'DEF');
grid minor
axis tight

end