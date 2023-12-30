function err= cal_err(up_sta,def_rate)
%
%  err_station=[error_linear, error_ice4, error_ice5, error_annual]
error_linear = abs(def_rate(1) - up_sta(1) )/ up_sta(1)  * 100;
error_ice4 = abs(def_rate(1) - up_sta(1) + up_sta(3) )/ up_sta(1)  * 100;
error_ice5 = abs(def_rate(1) - up_sta(1) + up_sta(4) )/ up_sta(1)  * 100;
error_annual = abs(def_rate(2) - up_sta(2) )/ up_sta(2)  * 100;

err=[error_linear, error_ice4, error_ice5, error_annual];
end