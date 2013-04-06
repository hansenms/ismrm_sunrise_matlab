function [img,snr,g,noise_psf] = ismrm_cartesian_iterative_SENSE(inp,samp_mat,csm,reg,replicas)
%
%   ismrm_cartesian_iterative_SENSE(inp,samp_pat,csm,reg,replicas)
%
%   Non-Cartesian SENSE reconstruction. Uses Matlab LSQR to solve.
%
%   It is recommended to input data with pre-whitened noise scale to sd=1.
%
%   Gridding weights should be scaled such that sum(w(:)) is equal to the
%   fraction of k-space that the samples cover, i.e. pi*0.5^2 for a circle.
%
%   INPUT:
%     - inp         [kx,ky,coils]        : Input k-space data (vector)
%     - samp_mat    [kx,ky]              : mask 1 or 0
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

if nargin<4,
    reg = [];
end

if nargin<4,
    replicas = 100;
end

csm_sq = sum(csm .* conj(csm),3); csm_sq(csm_sq < eps) = 1;
M = spdiag(sqrt(csm_sq)); %Preconditioner

if (isempty(reg)),
    E = @(x,tr) ismrm_encoding_cartesian_SENSE(x,csm,samp_mat,tr);
    reg_out = [];
else
    reg_mask = 1./reg;
    reg_mask = reg_mask * (numel(reg_mask)/sum(reg_mask(:)));
    E_base = @(x,tr) ismrm_encoding_cartesian_SENSE(x,csm,samp_mat,tr);
    E = @(x,tr) regularized_E(x,E_base,reg_mask,tr);
    reg_out = zeros(numel(reg_mask),1);
    M = M + spdiag(reg_mask(:));
end

img = lsqr(E, [inp(:); reg_out], 1e-3,30,M);
img = reshape(img,size(csm,1),size(csm,2));

if (nargout > 1),
    image_formation_func = @(x) reshape(lsqr(E,[x; reg_out],1e-3,30,M),[size(csm,1) size(csm,2)]);
    [snr,g,noise_psf] = ismrm_pseudo_replica(inp(:), image_formation_func,replicas);
    csm_sq = sum(csm .* conj(csm),3); csm_sq(csm_sq < eps) = 1;
    g = g .* sqrt(csm_sq);
end

return

function out = regularized_E(x,E,reg_mask,transpose_indicator)
    numimgel = length(reg_mask(:));
    if (strcmp(transpose_indicator,'transp')),
        numkel = length(x(:))-numimgel;
        out = E(x(1:numkel),transpose_indicator) + reg_mask(:).*x((numkel+1):end);
    elseif (strcmp(transpose_indicator, 'notransp')),
        out = [E(x,transpose_indicator);reg_mask(:).*x];
    else
        error('Transpose flag not appropriately defined');
    end
return