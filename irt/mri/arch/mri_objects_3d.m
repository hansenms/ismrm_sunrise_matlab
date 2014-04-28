 function [xf, yf] = mri_objects_3d(varargin)
%function [xf, yf] = mri_objects_3d('type1', params1, 'type', params2, ...)
% Generate inline functions for describing image-domain objects
% and Fourier domain spectra of simple structures such as rectangles.
% These functions are useful for simple "idealized" MRI simulations
% where the data is modeled as analytical Fourier samples,
% i.e., no field inhomogeneity and no relaxation effects.
% in:
%	(type, params)	e.g. 'rect', params, 'gauss', params, ...
%
% params:
%	'rect'		[N,7]	[xcent ycent zcent xwidth ywidth zwidth value]
%	'gauss'		[N,7]	[xcent ycent zcent xwidth ywidth zwidth value]
%
% out:
%	xf	inline(x,y,z) that returns 3D image-domain picture
%	yf	inline(u,v,w) that returns 3D Fourier-domain values
%
% Copyright 2004-4-20, Jeff Fessler, The University of Michigan

if nargin < 1, help(mfilename), error args, end
if nargin == 1 && streq(varargin{1}, 'test'), mri_objects_3d_test; return, end

warn 'todo: superceded by mri_objects - need to update code to not use this!'

fov = [];
if streq(varargin{1}, 'fov')
	fov = varargin{2};
	varargin = {varargin{3:end}};
end

if streq(varargin{1}, 'rect2')
	[xf yf] = mri_objects_3d_rect2(fov, varargin{2:end});
return
end
if streq(varargin{1}, 'test2')
	[xf yf] = mri_objects_3d_test2(varargin{2:end});
return
end
if streq(varargin{1}, 'test3')
	[xf yf] = mri_objects_3d_test3(varargin{2:end});
return
end

xf = '';
yf = '';
while length(varargin)
	type = varargin{1};
	params = varargin{2};
	varargin = {varargin{3:end}};
	switch type
	case 'rect'
		[xs ys] = mri_objects_3d_rect(params);
	case 'gauss'
		[xs ys] = mri_objects_3d_gauss(params);
	otherwise
		error('unknown object type %s', type)
	end

	if streq(xf, '')
		xf = xs;
		yf = ys;
	else
		xf = [xf ' + ' xs];
		yf = [yf ' + ' ys];
	end
end

xf = inline(xf, 'x', 'y', 'z');
yf = inline(yf, 'u', 'v', 'w');


%
% mri_objects_3d_rect()
% rect / sinc
% in:
%	param		[N,7]	[xcent ycent zcent xwidth ywidth zwidth value]
% out:
%	xf	string: \sum_n a * rect((x-xc)/xw) * rect((y-yc)/yw) * rect((z-zc)/zw)
%	yf	string: 3D Fourier of xf (sinc * sinc * sinc)
%
function [xf, yf] = mri_objects_3d_rect(param)

xf = '';
yf = '';
for ii=1:nrow(param)
	if ii > 1
		xf = [xf ' + '];
		yf = [yf ' + '];
	end
	xs = sprintf('%g', param(ii,7));
	ys = '';
	yfac = param(ii,7);
	xl = {'x', 'y', 'z'};
	yl = {'u', 'v', 'w'};
	for id=1:3
		cc = param(ii,id);
		ww = param(ii,id+3);
		if ~isinf(ww)
			xs = [xs sprintf(' .* rect((%s-(%g))/%g)', ...
				xl{id}, cc, ww)];
			yfac = yfac * ww;
			ys = [ys sprintf(' .* sinc(%s*(%g))', ...
				yl{id}, ww)];
			ys = [ys sprintf(' .* exp(-2i*pi*(%s*(%g)))', ...
				yl{id}, cc)];
		end

	end
	ys = [sprintf('%g', yfac) ys];
	xf = [xf xs];
	yf = [yf ys];
end


%
% mri_objects_3d_gauss()
% generate inline 2D gaussian functions (and spectrum thereof)
% in:
%	param	[N,7]	[xcent ycent zcent xwidth ywidth zwidth value]
% out:
%	xf	string: \sum_n a * gauss((x-xc)/xw) * gauss((y-yx)/yw)
%	yf	string: 3D Fourier of xf (gauss * gauss * gauss)
%
function [xf, yf] = mri_objects_3d_gauss(param)

xf = '';
yf = '';
for ii=1:nrow(param)
	if ii > 1
		xf = [xf ' + '];
		yf = [yf ' + '];
	end

	xs = sprintf('%g', param(ii,7));
	ys = '';
	yfac = param(ii,7);
	xl = {'x', 'y', 'z'};
	yl = {'u', 'v', 'w'};
	for id=1:3
		cc = param(ii,id);
		ww = param(ii,id+3) / sqrt(log(256)) * sqrt(2*pi);
		if ~isinf(ww)
			xs = [xs sprintf(' .* exp(-pi*((%s-(%g))/(%g)).^2)', ...
				xl{id}, cc, ww)];
			yfac = yfac * ww;

			ys = [ys sprintf(' .* exp(-pi*(%s*(%g)).^2)', ...
				yl{id}, ww)];
			ys = [ys sprintf(' .* exp(-2i*pi*%s*(%g))', ...
				yl{id}, cc)];
		end
	end
	ys = [sprintf('%g', yfac) ys];
	xf = [xf xs];
	yf = [yf ys];
end


% mri_objects_3d_rect2()
% test case: rect, size is half of fov
function [xf, yf] = mri_objects_3d_rect2(fov, arg)

if isvar('arg') && ~isempty(arg), error 'no arg for rect2', end
if isempty(fov)
	error 'fov required for rect2'
elseif length(fov) == 1
	fov = [fov fov inf];
elseif length(fov) == 2
	fov = [fov inf];
elseif length(fov) ~= 3
	error 'bad fov for rect2'
end
[xf yf] = mri_objects_3d('rect', [0 0 0 fov/2 1]);


% mri_objects_3d_test2()
function [xf, yf] = mri_objects_3d_test2(arg)

rp = [ ...	% rect
	0 0 0		200 200 inf	1;
	-50 -50 0	40 40 inf	1;
	50 -50 0	20 20 inf	1;
	0 50 0		50 50 inf	1;
];

gp = [ ...	% gaussian bumps
	-70 0 0		1 1 inf	1;
	-60 0 0		2 2 inf	1;
	-50 0 0		3 3 inf	1;
	-40 0 0		4 4 inf	1;
	-20 0 0		5 5 inf	1;
	00 0 0		6 6 inf	1;
	20 0 0		7 7 inf	1;
	50 0 0		8 8 inf	1;
];

if isvar('arg') && streq(arg, 'cm')
	rp(:,1:6) = rp(:,1:6) / 10;
	gp(:,1:6) = gp(:,1:6) / 10;
end

[xf yf] = mri_objects_3d('rect', rp, 'gauss', gp);


% mri_objects_3d_test3()
function [xf, yf] = mri_objects_3d_test3(arg)

rp = [ ...	% rect
	0 0 0		200 200 10	1;
	-50 -50 0	40 40 10	1;
	50 -50 10	20 20 10	1;
	0 50 -10	50 50 10	1;
];

gp = [ ...	% gaussian bumps
	-70 0 0		1 1 1	1;
	-60 0 0		2 2 2	1;
	-50 0 0		3 3 3	1;
	-40 0 0		4 4 4	1;
	-20 0 0		5 5 5	1;
	00 0 0		6 6 6	1;
	20 0 0		7 7 7	1;
	50 0 0		8 8 8	1;
];

if isvar('arg') && streq(arg, 'cm')
	rp(:,1:6) = rp(:,1:6) / 10;
	gp(:,1:6) = gp(:,1:6) / 10;
end

[xf yf] = mri_objects_3d('rect', rp, 'gauss', gp);


function mri_objects_3d_test
if 0
	[xo yo] = mri_objects_2d('gauss', [1 2 3 4 5]);
	[xf yf] = mri_objects_3d('gauss', [1 2 0 3 4 inf 5]);
end
N = 128;
fov = 256; % mm
[xf yf] = mri_objects_3d('test2');
%[xf yf] = mri_objects_3d('fov', fov, 'rect2');
x1 = [-N/2:N/2-1] / N * fov;
[x1d x2d] = ndgrid(x1, x1);
xtrue = xf(x1d, x2d, 0);
k1 = [-N/2:N/2-1]/fov; % cartesian
[kk1 kk2] = ndgrid(k1, k1);
ytrue = yf(kk1, kk2, 0);
im(x1, x1, xtrue), cbar
if 0 % compare against old 2d version
	[xfo yfo] = mri_objects_2d('test');
	xold = xfo(x1d, x2d);
	yold = yfo(kk1, kk2);
	max_percent_diff(xtrue, xold)
	max_percent_diff(ytrue, yold)
end
