function dmtx = ismrm_calculate_noise_decorrelation_mtx(noise_samples, scale_factor)
%
%  dmtx = ismrm_calculate_noise_decorrelation_mtx(noise_samples, scale_factor)
%
%  Calculates the noise decorrelation matrix based on noise samples
%
%  INPUT:
%    - noise_samples   [samples, coils]   : Input noise samples
%                                           In general it is assumed that
%                                           last dimension is coil, [x,y,c]
%                                           is rehsaped to [x*y,c]
%    - scale_factor    scalar             : optional scale factor to take
%                                           noise sample BW, etc. into consideration
%
%  OUTPUT:
%    - dmtx  [coils, coils]               : Noise de-correlation
%                                           (prewhitening) matrix
%
%  If the goal of the noise decorrelation is to create data with noise SD = 1.0 in
%  real and imaginary channels, the scale_factor should be:
%
%     scale_factor = (T_acq_dwell/T_noise_dwell)*NoiseReceiverBandwidthRatio
%
%  Where T_acq_dwell is the sampling dwell time in the acquisition data being
%  decorrelated and T_noise_dwell is the dwell time in the noise sample
%  acquistion (they may not be the same). The NoiseReceiverBandwithRatio is
%  a number (typically smaller than 1) that takes into account that the
%  noise acquisition may not have flat freqeuency response. 
%
%  see Kellman P, McVeigh ER. Magn Reson Med. 2005 Dec;54(6):1439-47. Erratum in: Magn Reson Med. 2007 Jul;58(1):211-2.
%
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%

if nargin < 2,
    scale_factor = 1;
end

noise = reshape(noise_samples,numel(noise_samples)/size(noise_samples,length(size(noise_samples))), size(noise_samples,length(size(noise_samples))));
noise = permute(noise,[2 1]);
M = size(noise,2);
dmtx = (1/(M-1))*(noise*noise');
dmtx = inv(chol(dmtx,'lower'));
dmtx = dmtx*sqrt(2)*sqrt(scale_factor);
    
return
