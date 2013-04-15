%
% A simple script to test the creation of noise decorrelation functions
%
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%   Philip J. Beatty (philip.beatty@sri.utoronto.ca)
%
nsamples = 100000;
load noise_covariances.mat

Rn = eye(2);
Rn(2,1) = 0.1 + .2i;
Rn(1,2) = 0.1 - .2i;
noise = ismrm_generate_correlated_noise(nsamples, Rn_broken_8);
noise = ismrm_generate_correlated_noise(nsamples, Rn);

dmtx1 = ismrm_calculate_noise_decorrelation_mtx(noise);

Rn_estimate = ismrm_estimate_covariance_matrix(noise);

dmtx2 = ismrm_calculate_noise_decorrelation_mtx_from_covariance_mtx(Rn_estimate);