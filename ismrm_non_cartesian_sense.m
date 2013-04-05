function [img,snr,g,noise_psf] = ismrm_non_cartesian_sense(inp,k,w,csm,reg,replicas)
%
%   [img,snr,g,noise_psf] = ismrm_non_cartesian_sense(inp,k,w,csm,replicas)
%
%   Non-Cartesian SENSE reconstruction. Uses Matlab LSQR to solve.
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
%     - csm         [x,y,coil]           : Coil sensitivities 
%     - reg         [x,y]                : Image space regularization mask
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

if nargin<5,
    reg = [];
end

if nargin<6,
    replicas = 100;
end

%Prepare NUFFT
N = [size(csm,1) size(csm,2)];
J = [5 5];
K = N*2;

w = w*prod(K);

nufft_st = nufft_init(k*2*pi,N,J,K,N/2,'minmax:kb');

E = @(x,tr) ismrm_encoding_non_cartesian_SENSE(x,csm,nufft_st,w,tr);

img = lsqr(E, inp(:) .* repmat(sqrt(w),[size(csm,3),1]), 1e-3,30);
img = reshape(img,size(csm,1),size(csm,2));

if (nargout > 1),
    image_formation_func = @(x) reshape(lsqr(E,x .* repmat(sqrt(w),[size(csm,3),1]),1e-3,30),[size(csm,1) size(csm,2)]);
    [snr,g,noise_psf] = ismrm_pseudo_replica(inp(:), image_formation_func,replicas);
    csm_sq = sum(csm .* conj(csm),3); csm_sq(csm_sq < eps) = 1;
    g = g .* sqrt(csm_sq);
end

return