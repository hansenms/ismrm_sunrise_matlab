function im_kernel = ismrm_transform_kernel_to_image_space(kernel, out_size)
%
%  im_kernel = ismrm_transform_kernel_to_image_space(kspace_kernel, out_size)
%
%  Transforms a k-space convolution kernel to image space (for
%  multiplication)
%
%  INPUT:
%    - kernel    [kx,ky,..., source_nc, target_nc]  :  k-space convolution kernel
%    - out_size  [size_x,size_y]     :  Image size, e.g. 128,128
%
%  OUTPUT:
%    - im_kernel  [nx,ny,nc]    : Image multiplication kernel
%                            sum(w) == numel(w)
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%   Philip Beatty (philip.beatty@sri.utoronto.ca)
%

spatial_dimensions = length(out_size);
source_coils = size(kernel,spatial_dimensions+1);
target_coils = size(kernel,spatial_dimensions+2);

%Flip kernel since we will be FFTing to image space for pixel wise mult
for d = 1:spatial_dimensions,
    kernel = flipdim(kernel,d);
end

out_dimensions = [out_size(:);source_coils;target_coils].';

im_kernel = ismrm_transform_kspace_to_image(kernel, [spatial_dimensions:-1:1], out_dimensions);
im_kernel = im_kernel .* sqrt(size(im_kernel,1)*size(im_kernel,2)); % trying not to change MH scaling.

return 