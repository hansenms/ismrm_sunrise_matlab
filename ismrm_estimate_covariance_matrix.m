function n_matrix = ismrm_estimate_covariance_matrix(noise_data)
%
%   Estimate covariance matrix from noise samples
%
%    INPUT:
%        noise_data [samples, coils] : zero mean gaussian complex noise samples
%    OUTPUT:
%        n_matrix [coils, coils] : estimate of the noise covariance matrix
% 
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%   Philip J. Beatty (philip.beatty@sri.utoronto.ca)
%

nc = size(noise_data, ndims(noise_data));

num_samples = numel(noise_data) / nc;

noise_data_matrix = reshape(noise_data, [num_samples, nc]);

n_matrix = noise_data_matrix.' * conj(noise_data_matrix);
n_matrix = n_matrix / (num_samples * 2);
n_matrix = 0.5 * (n_matrix + n_matrix');
