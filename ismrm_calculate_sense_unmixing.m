function [unmix, gmap] = ismrm_calculate_sense_unmixing(acc_factor, csm)
%
% function [unmix, gmap] = ismrm_calculate_sense_unmixing(acc_factor, csm)
%
% Calculates the unmixing coefficients for a 2D image
%
% INPUT:
%       acc_factor  scalar       : Acceleration factor, e.g. 2
%       csm1d       [x, y, coil] : Coil sensitivity map 
%
% OUTPUT:
%       unmix       [x, y, coil] : Image unmixing coefficients for a single x location 
%       gmap        [x, y]       : Noise enhancement map 
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%   Philip Beatty (philip.beatty@sri.utoronto.ca)
%

unmix = zeros(size(csm));

for x=1:size(csm,1), 
    unmix(x,:,:) = ismrm_calculate_sense_unmixing_1d(acc_factor, squeeze(csm(x,:,:))); 
end

if (nargout > 1),
   gmap = sqrt(sum(abs(unmix).^2,3)).*sqrt(sum(abs(csm).^2,3));
end
