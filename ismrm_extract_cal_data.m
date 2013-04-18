function cal_data = ismrm_extract_cal_data(data, sp)
%
%   cal_data = ismrm_extract_cal_data(data, sp)
%   
%   Extract region of calibration data from an internally calibrated data
%   set.
%
%   INPUT:
%       data [kx,ky,coil]           : accelerated data set with zeros in
%                                     the unacquired frames.
%       sp [kx, ky]                 : sampling pattern mask.
%                                     0 = not acquired
%                                     >0 = acquired
%
%   OUTPUT:
%       cal_data [kx, ky, coil]     : calibration data sub matrix
%
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%   Philip J. Beatty (philip.beatty@sri.utoronto.ca)
%
    ky_projection = max(sp, [], 1);
    kx_projection = max(sp, [], 2);
    cal_indx = find(kx_projection > 1);
    cal_indy = find(ky_projection > 1);

    kx_cal_bounds = [min(cal_indx):max(cal_indx)];
    ky_cal_bounds = [min(cal_indy):max(cal_indy)];
    cal_data = data(kx_cal_bounds, ky_cal_bounds, :);