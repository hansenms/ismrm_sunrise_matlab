function [img,snr,g,noise_psf] = ismrm_non_cartesian_sense(inp,k,w,csm,reg,lambda,replicas)
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
%     - lambda      scalar               : regularization factor
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


max_iterations = 20;
limit = 1e-3;

if nargin<5,
    reg = [];
end

if nargin<6,
    lambda = 1;
end

if nargin<7,
    replicas = 100;
end

%Prepare NUFFT
N = [size(csm,1) size(csm,2)];
J = [5 5];
K = N*2;

w = w*prod(K);

nufft_st = nufft_init(k*2*pi,N,J,K,N/2,'minmax:kb');

csm_sq = sum(csm .* conj(csm),3); csm_sq(csm_sq < eps) = 1;
M = spdiag(sqrt(csm_sq)); %Preconditioner

if (isempty(reg)),
    E = @(x,tr) ismrm_encoding_non_cartesian_SENSE(x,csm,nufft_st,w,tr);
    reg_out = [];
else
    reg_mask = 1./reg;
    reg_mask = reg_mask * (numel(reg_mask)/sum(reg_mask(:)));
    E_base = @(x,tr) ismrm_encoding_non_cartesian_SENSE(x,csm,nufft_st,w,tr);
    E = @(x,tr) regularized_E(x,E_base,reg_mask,lambda,tr);
    reg_out = zeros(numel(reg_mask),1);
    M = M + spdiag(reg_mask(:));
end

img = lsqr(E, [inp(:) .* repmat(sqrt(w),[size(csm,3),1]); reg_out], limit,max_iterations,M);
img = reshape(img,size(csm,1),size(csm,2));

if (nargout > 1),
    image_formation_func = @(x) reshape(lsqr(E,[x .* repmat(sqrt(w),[size(csm,3),1]); reg_out],limit,max_iterations,M),[size(csm,1) size(csm,2)]);
    [snr,g,noise_psf] = ismrm_pseudo_replica(inp(:), image_formation_func,replicas);
    csm_sq = sum(csm .* conj(csm),3); csm_sq(csm_sq < eps) = 1;
    g = g .* sqrt(csm_sq);
end

return

function out = regularized_E(x,E,reg_mask,lambda,transpose_indicator)
    numimgel = length(reg_mask(:));
    if (strcmp(transpose_indicator,'transp')),
        numkel = length(x(:))-numimgel;
        out = E(x(1:numkel),transpose_indicator) + lambda*reg_mask(:).*x((numkel+1):end);
    elseif (strcmp(transpose_indicator, 'notransp')),
        out = [E(x,transpose_indicator);lambda*reg_mask(:).*x];
    else
        error('Transpose flag not appropriately defined');
    end
return