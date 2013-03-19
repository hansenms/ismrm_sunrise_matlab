 function [xf, yf] = mri_object_rect(param)
%function [xf, yf] = mri_object_rect(param)
%	generate inline 2D rect functions and 2D sinc functions
%	for crude MRI simulations with analytical Fourier samples
%	in:
%		param		[N,5]	[xcent ycent xwidth ywidth value]
%	out:
%		xf	inline of \sum_n a * rect((x-xc)/xw) * rect((y-yx)/yw)
%		yf	inline of 2D Fourier of xf (sinc * sinc)
%
%	Copyright 2003-7-23	Jeff Fessler	The University of Michigan

if nargin < 1, help(mfilename), error args, end
warning 'obsolete: use mri_objects instead'

xf = '0';
yf = '0';
for ii=1:nrow(param)
	tmp = '+ %g * rect((x-%g)/%g) .* rect((y-%g)/%g)';
	tmp = sprintf(tmp, param(ii,5), param(ii,1), param(ii,3), ...
		param(ii,2), param(ii,4));
	xf = [xf tmp];

	tmp = '+ %g * sinc(u*%g) .* sinc(v*%g)';
	tmp = [tmp ' .* exp(-1i*2*pi*(u*%g+v*%g))'];
	tmp = sprintf(tmp, ...
		param(ii,5) * param(ii,3) * param(ii,4), ...
		param(ii,3), param(ii,4), ...
		param(ii,1), param(ii,2));
	yf = [yf tmp];
end

xf = inline(xf, 'x', 'y');
yf = inline(yf, 'u', 'v');
