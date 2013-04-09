function unmix = ismrm_calculate_jer_unmixing(jer_lookup, acc_factor, ccm, verbose)
%
%   [unmix, gmap] = ismrm_calculate_grappa_unmixing(source_data, kernel_size, acc_factor, csm, target_data, data_mask, verbose)
%   
%   Calculates b1-weighted image space GRAPPA unmixing coefficients.
%
%   INPUT:
%       source_data [kx,ky,coil]   : Source data for grappa kernel estimation (k-space)
%       kernel_size [kx,ky]        : e.g. [4 5]
%       acc_factor  scalar         : Acceleration factor, e.g. 2
%       data_mask   [kx,ky]        : '1' = calibration data, '0' = ignore
%       csm         [x,y,c]        : Coil sensitivity map, if empty, it
%                                    will be estimated from the reference lines.
%       target_data [kx,ky,coil]   : Target coil data, defaults to source data
%       verbose     bool           : Set true for verbose output
%
%   OUTPUT:
%       unmix [x,y,coil]           : Image unmixing coefficients
%       gmap  [x, y]               : Noise enhancement map 
%
%   Typical usage:
%       [unmix] = calculate_grappa_unmixing(source_data, [5 4], 4);
%
%
%   Notes:
%     - The unmixing coefficients produced by this routine produce uniform 
%       noise distribution images when there is no acceleration, i.e. the
%       noise in each pixel will be input noise * g-factor, where g-factor
%       is sqrt(sum(abs(unmix).^2,3)).
%
%       If you have coil sensitivities where the RSS of the coil
%       sensitivites is not 1 in each pixel, e.g. as obtained with a
%       seperate calibration scan using a body coil, and you would like a
%       uniform sensitivity image. You must apply that weighting after the
%       parallel imaging reconstruction by dividin with the RSS of the coil
%       sensitivites. 
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%


%%
% Compute unaliasing kernels
    
if (verbose),
    fprintf('Calculating unaliasing kernels...\n');
end
tic
kernel_size = [size(jer_lookup,1) size(jer_lookup,2)];
target_location = bitshift(kernel_size, -1)+1;
num_channels = size(ccm, 3);
kernel = zeros(kernel_size(1),kernel_size(2), num_channels, num_channels); %*acc_factor,coils,target_coils);

for ic = 1:num_channels,
    kernel(target_location(1), target_location(2), ic, ic) = 1;
end

%[kx_cal,ky_cal] = ind2sub(size(data_mask),[find(data_mask == 1,1,'first') find(data_mask == 1,1,'last')]);
for s=1:(acc_factor-1),
   kernel_mask = zeros(kernel_size);
   kernel_mask(:,s:acc_factor:end) = 1;
   k = ismrm_compute_kspace_unaliasing_coefficients(jer_lookup, kernel_mask);
   kernel = kernel + k;
end
toc

%%
% Form unmixing images from channel combination maps and kernels
if (verbose),
    fprintf('Merging unaliasing and channel combination images...\n');
end
tic
unmix = ismrm_calculate_unmixing_images_from_kspace_kernels(kernel, ccm);
toc

if (verbose),
    fprintf('done.\n');
end

return
        
