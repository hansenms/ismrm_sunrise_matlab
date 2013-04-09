function y = ismrm_rss(x,dim)
%
%   [mag] = ismrm_rss(samples, dim)
%
%   Computes root-sum-of-squares along a single dimension.
%
%
%   INPUT:
%     - x   : multi-dimensional array of samples
%     - dim : dimension of reduction; defaults to last dimension
%
%   OUTPUT:
%
%     - y       : root sum of squares result
%
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%   Philip J. Beatty (philip.beatty@sri.utoronto.ca)
%


if nargin==1
    dim=ndims(x);
else
    if isempty(dim); dim=ndims(x); end
end

y = squeeze(sqrt(sum(real(x).^2 + imag(x).^2,dim)));
