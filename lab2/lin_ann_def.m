function rate = lin_ann_def(def_ref, t)
%   calculate linear trend and annual amplitude of deformations

A=[ones(size(t)), t/365,cos(2*pi * t/365), sin(2*pi * t/365)];
para_def=(A' * A)\(A' * def_ref);
out1 = para_def(2);
out2 = sqrt(para_def(3)^2 + para_def(4)^2);
rate=[out1, out2];
end