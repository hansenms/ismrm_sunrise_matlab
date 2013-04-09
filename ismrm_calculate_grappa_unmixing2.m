function unmix = ismrm_calculate_grappa_unmixing2(source_data, kernel_size, acc_factor, data_mask, csm, target_data, im_shape, verbose)
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
% Validate Input
if nargin < 3,
   error('At least 4 arguments needed'); 
end

if nargin < 4,
    data_mask = [];
end

if nargin < 5,
    csm = [];
end

if nargin < 6,
    target_data = [];
end

if nargin < 7,
    im_shape = []
end
if nargin < 8,
    verbose = false;
end

if (isempty(target_data)),
        target_data = source_data;
end

if( isempty(im_shape) ),
    im_shape = [size(source_data,1) size(source_data,2)];
end

if (isempty(data_mask)),
    data_mask = ones(size(source_data,1),size(source_data,2));
end

if (length(size(source_data)) == 2),
    coils = 1;
else
    coils = size(source_data,length(size(source_data)));
end

if (length(size(target_data)) == 2),
    target_coils = 1;
else
    target_coils = size(target_data,length(size(target_data)));
end

%%
% Compute grappa kernels
    
if (verbose),
    fprintf('Calculating grappa kernels...\n');
end
tic
kernel = zeros(kernel_size(1),kernel_size(2)*acc_factor,coils,target_coils);

%[kx_cal,ky_cal] = ind2sub(size(data_mask),[find(data_mask == 1,1,'first') find(data_mask == 1,1,'last')]);
for s=1:acc_factor,
   kernel_mask = zeros(kernel_size(1),kernel_size(2)*acc_factor);
   kernel_mask(:,s:acc_factor:end) = 1;
   k = ismrm_estimate_convolution_kernel2(source_data,kernel_mask,target_data);
   kernel = kernel + k;
end
toc

%%
% coil sensitivity maps, if not provided

%If csm is not provided, we will estimate it.
if (isempty(csm)),
    if (verbose),
        fprintf('Estimating coil sensitivity...');
    end
    
    tic
    %Apply some filtering to avoid ringing
    f = hamming(max(sum(data_mask,1))) * hamming(max(sum(data_mask,2)))';
    fmask = zeros(size(source_data));
    fmask((1:size(f,1))+bitshift(size(source_data,1),-1)-bitshift(size(f,1),-1), ...
          (1:size(f,2))+bitshift(size(source_data,2),-1)-bitshift(size(f,2),-1), :) = ...
          repmat(f, [1 1 size(source_data,3)]);
    csm = ismrm_transform_kspace_to_image(source_data .* fmask, [1 2], [im_shape size(source_data,3)]);
    csm = ismrm_estimate_csm_walsh(csm); %Estimate coil sensitivity maps.
    toc
    if (verbose),
        fprintf('done.\n');
    end
end

%%
% Compute channel combination maps
tic
ccm = ismrm_compute_ccm(csm) .* (1.0/acc_factor);
toc

%%
% Form unmixing images from channel combination maps and kernels
tic
unmix = ismrm_calculate_unmixing_images_from_kspace_kernels(kernel, ccm);
toc

if (verbose),
    fprintf('done.\n');
end

return
        
