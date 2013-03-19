 function [xf, yf] = mri_object_gauss(param)
%function [xf, yf] = mri_object_gauss(param)
%	generate inline 2D gaussian functions (and spectrum thereof)
%	for crude MRI simulations with analytical Fourier samples
%	in:
%		param		[N,5]	[xcent ycent xwidth ywidth value]
%	out:
%		xf	inline of \sum_n a * gauss((x-xc)/xw) * gauss((y-yx)/yw)
%		yf	inline of 2D Fourier of xf (gauss * gauss)
%
%	Copyright 2003-7-23	Jeff Fessler	The University of Michigan

if nargin < 1, help(mfilename), error args, end
warning 'obsolete: use mri_objects instead'

xf = '0';
yf = '0';
for ii=1:nrow(param)
	wx = param(ii,3) / sqrt(log(256)) * sqrt(2*pi);
	wy = param(ii,4) / sqrt(log(256)) * sqrt(2*pi);
	tmp = '+ %g * exp(-pi*((x-%g)/%g).^2) .* exp(-pi*((y-%g)/%g).^2)';
	tmp = sprintf(tmp, param(ii,5), param(ii,1), wx, param(ii,2), wy);
	xf = [xf tmp];

	tmp = '+ %g * exp(-pi*(u*%g).^2) .* exp(-pi*(v*%g).^2)';
	tmp = [tmp ' .* exp(-1i*2*pi*(u*%g+v*%g))'];
	tmp = sprintf(tmp, param(ii,5) * wx * wy, wx, wy, ...
		param(ii,1), param(ii,2));
	yf = [yf tmp];
end

xf = inline(xf, 'x', 'y');
yf = inline(yf, 'u', 'v');
