function [k_out,w] = ismrm_calculate_spiral_trajectory(acceleration)
%
%   [k,w] = ismrm_calculate_spiral_trajectory()
%
%   Simple wrapper for Brian Hargreaves vdspiral code:
%
%   http://mrsrl.stanford.edu/~brian/vdspiral/
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%


smax = 15000;	 % 150 T/m/s
gmax = 4;	 % G/cm
T = .000004;	 % Seconds
N = 32;		 % Interleaves
Fcoeff = [30]; 	% FOV decreases linearly from 24 to 12cm.
res = 1;
rmax = 5/res;		% cm^(-1), corresponds to 1mm resolution.

[k,g] = vds(smax,gmax,T,N,Fcoeff,rmax);

k_out = zeros(numel(k),2);
N = floor(N/acceleration);
for p=1:N,
    rotation = ((p-1) * 2 * pi)/N;
    k_out((1:numel(k))+(p-1)*numel(k),1) = (real(k) * cos(rotation)) + (imag(k) * sin(rotation));
    k_out((1:numel(k))+(p-1)*numel(k),2) = -(real(k) * sin(rotation)) + (imag(k) * cos(rotation));  
end

w = abs(repmat(g(:),[N 1])) .* abs(sin(angle(repmat(g(:),[N 1]))-angle(repmat(k(:),[N 1]))));
area_weights = 1;
w = w .* (area_weights/sum(w(:)));        

k_out = k_out ./ (2*max(abs(k_out(:))));
