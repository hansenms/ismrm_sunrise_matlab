function [img] = ismrm_transform_kspace_to_image(k, dim, img_shape)
%
%  [img] = ismrm_transform_kspace_to_image(k, dim)
%
%  Fourier transform from k-space to image space along a given or all 
%  dimensions
%
%  INPUT:
%    - k       [kx,ky,..]    : k-space data
%    - dim     vector        : Vector with dimensions to transform
%    - img_shape vector      : Set shape of output image
%
%  OUPUT:
%    - img    [x,y,...]      : Data in image space (along transformed
%                                                   dimensions)
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%   Philip J. Beatty (philip.beatty@sri.utoronto.ca)
%

if nargin < 2,
    dim = [];
end    

if nargin < 3,
    img_shape = [];
end
   
if isempty(dim),
    dim = ndims(k):-1:1;
end

if isempty(img_shape)
    img_shape = size(k);
end
    
img = k;
for d=1:length(dim),
    img = transform_one_dim(img, d, img_shape(d));
end

return

function img = transform_one_dim(k, dim, img_extent)
    img_shape = size(k);
    img_shape(dim) = img_extent;
    
    k_indices = repmat({':'},1, ndims(k));
    k_indices{dim} = (1:size(k,dim))+bitshift(img_extent-size(k,dim)+1,-1);
    
    img = zeros(img_shape);
    img(k_indices{:}) = k;
    
    img = fftshift(ifft(ifftshift(img, dim), [], dim), dim) .* sqrt(size(img,dim)); 
 return   
    
    
    
    
    
    
    