function noise_covariance_matrix = ismrm_estimate_covariance_matrix(noise_data)
%
%   Estimate covariance matrix from noise samples
%
%    INPUT:
%        TODO
%    OUTPUT:
%        TODO
% 
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%   Philip J. Beatty (philip.beatty@sri.utoronto.ca)
%

nc = size(noise_data, ndims(noise_data));

num_samples = prod(size(noise_data)) / nc

noise_data_matrix = reshape(noise_data, [num_samples, nc]);

noise_covariance_matrix = noise_data_matrix.' * conj(noise_data_matrix);    