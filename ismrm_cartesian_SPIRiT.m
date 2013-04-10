function [img,snr,g,noise_psf] = ismrm_cartesian_SPIRiT(inp,samp_mat,cal,csm,replicas)
%
%   [img,snr,g,noise_psf] = ismrm_cartesian_SPIRiT(inp,samp_mat,cal,csm,replicas)
%
%   Cartesian SPIRiT reconstruction. Uses Matlab LSQR to solve.
%
%   It is recommended to input data with pre-whitened noise scale to sd=1.
%
%
%   INPUT:
%     - inp         [nsamples,coils]     : Input k-space data (vector)
%     - samp_mat    [kx,ky]              : mask 1 or 0
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

max_iterations = 50;
limit = 1e-3;

if nargin<5,
    replicas = 100;
end

if nargin<4,
    csm = [];
end

samples = sum(samp_mat(:) > 0);
coils = numel(inp)/samples;
recon_elements = numel(samp_mat)*coils;
mat_size = size(samp_mat);

%Estimate SPIRiT kernel
kernel_size = 7;
kernel_mask = ones(kernel_size,kernel_size);
kernel_mask(bitshift(kernel_size+1,-1),bitshift(kernel_size+1,-1)) = 0;
%kernel_mask(4,4) = 0;
kernel = ismrm_estimate_convolution_kernel(cal, kernel_mask);

padded_kernel = ismrm_transform_kernel_to_image_space(kernel,mat_size);

E = @(x,tr) ismrm_system_cartesian_SPIRiT(x,samp_mat,padded_kernel,tr);
img = lsqr(E, [inp(:);zeros(recon_elements,1)], limit,max_iterations);
img = reshape(img,mat_size(1),mat_size(2),coils);

if (isempty(csm)),
    %We have no coil sensitivity map, do RMS coil combination
    img = sqrt(sum(abs(img_spirit).^2,3));
else
    csm_sq = sum(csm .* conj(csm),3); csm_sq(csm_sq < eps) = 1;
    img = sum(conj(csm) .* img,3) ./ csm_sq;
end


if (nargout > 1),
    if (isempty(csm)),
        image_formation_func = @(x) sum(abs(reshape(lsqr(E,[x;zeros(recon_elements,1)],limit,max_iterations),[mat_size(:); coils]')).^2,3);     
    else
        image_formation_func = @(x) sum(conj(csm) .* reshape(lsqr(E,[x;zeros(recon_elements,1)],limit,max_iterations),[mat_size(:); coils]'),3) ./ csm_sq;
    end
    [snr,g,noise_psf] = ismrm_pseudo_replica(inp(:), image_formation_func,replicas);
    if (~isempty(csm)),
        g = g .* sqrt(csm_sq);
    end
end

