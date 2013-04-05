function noise = ismrm_generate_correlated_noise(im_shape, noise_covariance_matrix)
%
%   Generates noise that is correlated between channels
%
%    INPUT:
%        im_shape: matrix shape ([nx ny])
%        noise_covariance_matrix: nc x nc noise covariance matrix.  Make sure it is valid (positive definite)
%
%    OUTPUT:
%        nx x ny x nc matrix with correlated gaussian noise
% 
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%   Philip J. Beatty (philip.beatty@sri.utoronto.ca)
%

%%
% Validate input

assert(nargin == 2, 'valid number of arguments: 2')


%%
% Compute noise
nc = size(noise_covariance_matrix,1);
uncorrelated_noise = complex(randn(nc, prod(im_shape)), randn(nc, prod(im_shape)));

correlation_transform = chol(noise_covariance_matrix);

correlated_noise = correlation_transform * uncorrelated_noise;

noise = reshape( permute(correlated_noise, [2 1]), [im_shape nc]);

