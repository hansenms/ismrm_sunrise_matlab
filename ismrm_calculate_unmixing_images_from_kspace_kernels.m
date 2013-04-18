function unmix = ismrm_calculate_unmixing_images_from_kspace_kernels(kernel, ccm)
%
%  unmix = ismrm_calculate_unmixing_images_from_kspace_kernels(kernel, ccm)
%
%  Compute unmixing images from k-space unaliasing kernels and channel
%  combination maps.
%
%  INPUT:
%        kernels [kx, ky, num_source_coils, num_target_coils] : k-space
%                 unaliasing kernels (for uniform undersampling pattern)
%
%        ccm [x, y, num_target_coils] :  channel combination maps
%
nx = size(ccm,1);
ny = size(ccm,2);
num_source_coils = size(kernel,3);
num_target_coils = size(kernel,4);
assert( num_target_coils == size(ccm,3), 'num_target_coils in kernels does not match ccm');

unmix = zeros([nx ny num_source_coils]);

im_kernel = ismrm_transform_kernel_to_image_space(kernel, [nx ny]);

for c = 1:num_source_coils,
    unmix(:,:,c) = sum(squeeze(im_kernel(:,:,c,:)) .* ccm, 3);
end

return
        
