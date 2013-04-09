function unmix = ismrm_calculate_unmixing_images_from_kspace_kernels(kernel, ccm)
%
% kernels: nx x ny x num_source_coils x num_target_coils
%
% ccm: channel combination maps
%
nx = size(ccm,1);
ny = size(ccm,2);
num_source_coils = size(kernel,3);
num_target_coils = size(kernel,4);
assert( num_target_coils == size(ccm,3), 'num_target_coils in kernels does not match ccm');

unmix = zeros([nx ny num_source_coils]);

im_kernel = ismrm_transform_kernel_to_image_space(kernel, [nx ny]);

for c = 1:num_source_coils,
  %  ismrm_imshow(abs(squeeze(im_kernel(:,:,c,:))));
    unmix(:,:,c) = sum(squeeze(im_kernel(:,:,c,:)) .* ccm, 3);
end

return
        
