function [snr,g,noise_psf] = ismrm_pseudo_replica(in, image_formation_func, reps)
%
%  [snr,g] = ismrm_pseudo_replica(in, image_formation_func, reps)
%
%  Performs multiple reconstructions with the supplied image formation
%  function while adding noise with stdev=1.0 to each replica.
%
%  The function 'image_formation_func' can be used to generate an image
%  with:
%      image = image_formation_func(in)
%
%  This function assumes:
%    a) The noise in the original input data is white (after decorrelation)
%    b) The image_formation_function does has overall scale factor of 1.0
%       (Any additional scaling with scale the g-map).
%
%
%  INPUT:
%    - in        : Input data, any format that the image formation function
%                  expects.
%  OUTPUT:
%    - snr       : An image in SNR units.
%    - g         : A g-map (assuming image_formation_func doesn't scale)
%    - noise_psf : Point spread function of the noise
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%

baseline = image_formation_func(in);

for r=1:reps,
    fprintf('Running pseudo replica %d/%d\n',r,reps);
    n = complex(randn(size(in)),randn(size(in)));
    s = in + n;
    tmp = image_formation_func(s);
    img_noise_rep(:,r) = tmp(:);
end

img_noise_rep = reshape(img_noise_rep,[size(tmp),reps]);
rep_dim = length(size(img_noise_rep));

g = std(abs(img_noise_rep + max(abs(img_noise_rep(:)))),[],rep_dim); %Measure variation, but add offset to create "high snr condition"
g(g < eps) = 1;
snr = mean(img_noise_rep,3)./g;

img_noise_rep = img_noise_rep - repmat(baseline,[ones(1,length(size(baseline))) reps]);

ftdims = 1:(rep_dim-1);
noise_psf = ismrm_transform_kspace_to_image(mean(abs(ismrm_transform_image_to_kspace(img_noise_rep,ftdims)).^2,3));

return