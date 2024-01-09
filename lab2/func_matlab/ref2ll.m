function [ref_lam,ref_phi] = ref2ll(ref)
%   station coordiante to longitude and latitude
ref_lam = atan2d(ref(2),ref(1));
ref_phi = atan2d(ref(3), sqrt(ref(1)^2 + ref(2)^2));
end