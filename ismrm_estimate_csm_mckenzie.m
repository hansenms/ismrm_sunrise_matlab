function [csm] = ismrm_estimate_csm_mckenzie(img)
%
%   [csm] = ismrm_estimate_csm_mckenzie(img)
%
%   Estimates relative coil sensitivity maps from a set of
%   channel-by-channel images.
%
%   McKenzie et al. (Magn Reson Med
%   2002;47:529-538.)
%
%   INPUT:
%     - img     [x,y,coil]   : Coil images
%
%   OUTPUT:
%
%     - csm     [x,y,coil    : Relative coil sensitivity maps
%
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%   Philip J. Beatty (philip.beatty@sri.utoronto.ca)
%

csm = img ./ repmat( sqrt(sum(abs(img).^2,3)), [1 1 size(img,3)]);