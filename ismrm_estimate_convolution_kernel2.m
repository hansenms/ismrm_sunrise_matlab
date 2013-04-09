function kernel = ismrm_estimate_convolution_kernel2(source_data, kernel_mask, target_data, verbose)
%
%   kernel = ismrm_estimate_convolution_kernel(source_data, kernel_mask, target_data, verbose)
%   
%   Estimates a convolution kernel based on the supplied calibration data
%   (source and target) and the kernel mask. 
%
%   INPUT:
%       source_data [kx,ky,coil]   : Source data kernel estimation (k-space)
%       kernel_mask [kx,ky]        : e.g [1 1 1; 0 0 0; 1 1 1] for a 3x3
%                                    kernel for a factor 2 GRAPPA.
%       target_data [kx,ky,coil]   : Target coil data, defaults to source data
%       verbose     bool           : Set true for verbose output
%
%   OUTPUT:
%       kernel [kx,ky,coil,coil]   : Kernel coefficients according to mask
%
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%

if nargin < 2,
    error('At least 3 arguments needed');
end

if nargin < 3,
    target_data = [];
end

if nargin < 4,
    verbose = 0;
end

if isempty(target_data),
    target_data = source_data;
end

if ndims(source_data) ~= 3,
    error('Source data should have 3 dimensions kx, ky, coil');
end

if ~ismatrix(kernel_mask),
    error('kernel mask should have 2 dimensions kx, ky');
end

if ndims(target_data) ~= 3,
    error('target data should have 3 dimensions kx, ky, coil');
end

[kx_offsets, ky_offsets] = ind2sub(size(kernel_mask), find(kernel_mask == 1));
kx_offsets = kx_offsets - bitshift(size(kernel_mask,1),-1)-1;
ky_offsets = ky_offsets - bitshift(size(kernel_mask,2),-1)-1;

kx_range = (1-min(kx_offsets)):size(source_data,1)-max(kx_offsets);
ky_range = (1-min(ky_offsets)):size(source_data,2)-max(ky_offsets);

A = zeros(length(kx_range)*length(ky_range),length(kx_offsets)*size(source_data,3));
b = zeros(length(kx_range)*length(ky_range),size(source_data,3));

for c=1:size(source_data,3),
    for p=1:length(kx_offsets),
        A(:,(c-1)*length(kx_offsets)+p) = reshape(source_data(kx_range+kx_offsets(p),ky_range+ky_offsets(p),c),length(kx_range)*length(ky_range),1);
    end
    b(:,c) = reshape(target_data(kx_range,ky_range,c),length(kx_range)*length(ky_range),1);
end

if (verbose == 1),
   if (size(A,1) < 3*size(A,2)),
       message('WARNING: The number of calibration points smaller than 3 times the number of kernel coefficients.');
       message('Calibration data potentially insufficient');
   end
end

%tic
%S = svd(A,0);
%A_inv = pinv(A'*A + eye(size(A'*A)).*(1e-3*max(abs(S(:)))).^2)*A';
%x = A_inv*b;
%toc

%tic
Rss = A'*A;
num_basis = size(Rss,1);
Rst = A'*b;
Rss = Rss + eye(num_basis)*0.01 * trace(Rss) / num_basis;
x = Rss \ Rst;
%toc

kernel = repmat(kernel_mask,[1 1 size(source_data,3) size(source_data,3)]);
kernel(kernel == 1) = x(:);

return
