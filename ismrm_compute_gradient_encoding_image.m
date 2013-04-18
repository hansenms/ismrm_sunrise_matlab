function g = ismrm_compute_gradient_encoding_image(kx, ky, nx, ny)
%
% g = ismrm_compute_gradient_encoding_image(kx, ky, nx, ny)
%
% Computes nx x ny complex exponential image:
%
% g(x, y) = exp(i2pi * kx*x + ky*y)
%
% x samples are spaced 1/nx in the range [-0.5,0.5)
% y samples are spaced 1/ny in the range [-0.5, 0.5)
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%   Philip J. Beatty (philip.beatty@sri.utoronto.ca)
%

x = repmat((((1:nx) - 1)/nx - 0.5).', [1 ny]);
y = repmat((((1:ny) - 1)/ny - 0.5)  , [nx 1]);

g = exp(2i*pi * (x * kx + y * ky));
