 function [xf, yf] = mri_objects_2d(varargin)
%function [xf, yf] = mri_objects_2d('type1', params1, 'type', params2, ...)
% Generate inline functions for describing image-domain objects
% and Fourier domain spectra of simple structures such as rectangles.
% These functions are useful for simple "idealized" MRI simulations
% where the data is modeled as analytical Fourier samples,
% i.e., no field inhomogeneity and no relaxation effects.
% in:
%	(type, params)	e.g. 'rect', params, 'gauss', params, ...
%
% params:
%	'rect'		[N,5]	[xcent ycent xwidth ywidth value]
%	'gauss'		[N,5]	[xcent ycent xwidth ywidth value]
%
% out:
%	xf	inline(x,y) that returns 2D image-domain picture
%	yf	inline(u,v) that returns 2D Fourier-domain values
%
% Copyright 2004-4-20, Jeff Fessler, The University of Michigan

if nargin < 1, help(mfilename), error args, end
if nargin == 1 && streq(varargin{1}, 'test')
	[xf yf] = mri_objects_2d_test;
return
end

xf = '';
yf = '';
while length(varargin)
	type = varargin{1};
	params = varargin{2};
	varargin = {varargin{3:end}};
	if streq(type, 'rect')
		[xx yy] = mri_objects_2d_rect(params);
	elseif streq(type, 'gauss')
		[xx yy] = mri_objects_2d_gauss(params);
	else
		error(sprintf('unknown object type %s', type))
	end

	if streq(xf, '')
		xf = xx;
		yf = yy;
	else
		xf = [xf ' + ' xx];
		yf = [yf ' + ' yy];
	end
end

xf = inline(xf, 'x', 'y');
yf = inline(yf, 'u', 'v');


% rect / sinc
% in:
%	param		[N,5]	[xcent ycent xwidth ywidth value]
% out:
%	xf	string: \sum_n a * rect((x-xc)/xw) * rect((y-yx)/yw)
%	yf	string: 2D Fourier of xf (sinc * sinc)
function [xf, yf] = mri_objects_2d_rect(param)

xf = '';
yf = '';
for ii=1:nrow(param)
	if ii > 1
		xf = [xf ' + '];
		yf = [yf ' + '];
	end
	tmp = '%g * rect((x-(%g))/%g) .* rect((y-(%g))/%g)';
	tmp = sprintf(tmp, param(ii,5), param(ii,1), param(ii,3), ...
		param(ii,2), param(ii,4));
	xf = [xf tmp];

	tmp = '%g * sinc(u*%g) .* sinc(v*%g)';
	tmp = [tmp ' .* exp(-1i*2*pi*(u*(%g)+v*(%g)))'];
	tmp = sprintf(tmp, ...
		param(ii,5) * param(ii,3) * param(ii,4), ...
		param(ii,3), param(ii,4), ...
		param(ii,1), param(ii,2));
	yf = [yf tmp];
end


% gauss
% generate inline 2D gaussian functions (and spectrum thereof)
% in:
%	param		[N,5]	[xcent ycent xwidth ywidth value]
% out:
%	xf	string: \sum_n a * gauss((x-xc)/xw) * gauss((y-yx)/yw)
%	yf	string:  2D Fourier of xf (gauss * gauss)
function [xf, yf] = mri_objects_2d_gauss(param)

xf = '';
yf = '';
for ii=1:nrow(param)
	if ii > 1
		xf = [xf ' + '];
		yf = [yf ' + '];
	end
	wx = param(ii,3) / sqrt(log(256)) * sqrt(2*pi);
	wy = param(ii,4) / sqrt(log(256)) * sqrt(2*pi);
	tmp = '%g * exp(-pi*((x-(%g))/(%g)).^2) .* exp(-pi*((y-(%g))/(%g)).^2)';
	tmp = sprintf(tmp, param(ii,5), param(ii,1), wx, param(ii,2), wy);
	xf = [xf tmp];

	tmp = '%g * exp(-pi*(u*%g).^2) .* exp(-pi*(v*%g).^2)';
	tmp = [tmp ' .* exp(-1i*2*pi*(u*%g+v*%g))'];
	tmp = sprintf(tmp, param(ii,5) * wx * wy, wx, wy, ...
		param(ii,1), param(ii,2));
	yf = [yf tmp];
end


function [xf, yf] = mri_objects_2d_test

[xf yf] = mri_objects_2d(...
	'rect', ...
	[	0 0	200 200	1;
		-50 -50	40 40	1;
		50 -50	20 20	1;
		0 50	50 50	1;
	], ...
	'gauss', ...	% gaussian bumps
	[	-70 0	1 1	1;
		-60 0	2 2	1;
		-50 0	3 3	1;
		-40 0	4 4	1;
		-20 0	5 5	1;
		00 0	6 6	1;
		20 0	7 7	1;
		50 0	8 8	1;
	]);
