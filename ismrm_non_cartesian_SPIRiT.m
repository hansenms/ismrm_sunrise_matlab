function [img,snr,g,noise_psf] = ismrm_non_cartesian_SPIRiT(inp,k,w,mat_size,cal,csm,replicas)
%
%   [img,snr,g,noise_psf] = ismrm_non_cartesian_SPIRiT(inp,k,w,mat_size,cal,csm,replicas)
%
%   Non-Cartesian SPIRiT reconstruction. Uses Matlab LSQR to solve.
%
%   It is recommended to input data with pre-whitened noise scale to sd=1.
%
%   Gridding weights should be scaled such that sum(w(:)) is equal to the
%   fraction of k-space that the samples cover, i.e. pi*0.5^2 for a circle.
%
%   INPUT:
%     - inp         [nsamples,coils]     : Input k-space data (vector)
%     - k           [nsamples,2]         : k-space coordinates, range-0.5:0.5
%     - w           [nsamples]           : vector of gridding weights
%     - mat_size    size_x,size_y        : Reconstruction matrix size
%     - cal         [kx,ky,coil]         : Fully sampled k-space calibration data 
%     - csm         [x,y,coil]           : Optional coil sensitivity map (for coil combination)
%     - replicas    scalar (dafault 100) : Number of replicas to run if SNR
%                                          is requested
%
%   OUTPUT:
%     - img         [x,y]                : Output image
%     - snr                              : An image in SNR units.
%    -  g                                : A g-map (assuming image_formation_func doesn't scale)
%    -  noise_psf                        : Point spread function of the noise
%
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%

max_iterations = 100;
limit = 1e-3;

if nargin<7,
    replicas = 100;
end

if nargin<6,
    csm = [];
end
samples = numel(w);
coils = numel(inp)/samples;
recon_elements = prod(mat_size)*coils;

%Estimate SPIRiT kernel
kernel_size = 7;
kernel_mask = ones(kernel_size,kernel_size);
kernel_mask(bitshift(kernel_size+1,-1),bitshift(kernel_size+1,-1)) = 0;
kernel = ismrm_estimate_convolution_kernel(cal, kernel_mask);

padded_kernel = ismrm_transform_kernel_to_image_space(kernel,mat_size);


%Prepare NUFFT
N = mat_size;
J = [5 5];
K = N*2;

w = w*prod(K);

nufft_st = nufft_init(k*2*pi,N,J,K,N/2,'minmax:kb');


E = @(x,tr) ismrm_system_non_cartesian_SPIRiT(x,nufft_st,w(:),padded_kernel,tr);
img_spirit = lsqr(E, [inp(:) .* repmat(sqrt(w(:)),[coils,1]);zeros(recon_elements,1)], limit,max_iterations);
img_spirit = reshape(img_spirit,mat_size(1),mat_size(2),coils);

if (isempty(csm)),
    %We have no coil sensitivity map, do RMS coil combination
    img = sqrt(sum(abs(img_spirit).^2,3));
else
    csm_sq = sum(csm .* conj(csm),3); csm_sq(csm_sq < eps) = 1;
    img = sum(conj(csm) .* img_spirit,3) ./ csm_sq;
end


if (nargout > 1),
    if (isempty(csm)),
        image_formation_func = @(x) sum(abs(reshape(lsqr(E,[x .* repmat(sqrt(w),[coils,1]);zeros(recon_elements,1)],limit,max_iterations),[mat_size(:); coils]')).^2,3);     
    else
        image_formation_func = @(x) sum(conj(csm) .* reshape(lsqr(E,[x .* repmat(sqrt(w),[size(csm,3),1]);zeros(recon_elements,1)],limit,max_iterations),[mat_size(:); coils]'),3) ./ csm_sq;
    end
    [snr,g,noise_psf] = ismrm_pseudo_replica(inp(:), image_formation_func,replicas);
    if (~isempty(csm)),
        g = g .* sqrt(csm_sq);
    end
end

