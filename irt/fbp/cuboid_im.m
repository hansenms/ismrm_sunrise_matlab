 function [phantom, params] = cuboid_im(ig, params, varargin)
%function [phantom, params] = cuboid_im(ig, params, varargin)
%
% generate ellipsoids phantom image from parameters:
%	[x_center y_center z_center x_diameter y_diameter z_diameter
%		xy_angle_degrees z_angle_degrees amplitude]
% in
%	ig		image_geom()
%	params		cuboid parameters.  if empty, use default
% Note!! diameter not radius
% option
%	'oversample'	over-sampling factor
%   'type'	char	'' (default);
%|          'lowmen1' one slice per time with default method
% out
%	phantom		[nx,ny,nz] image
%
% Yong Long, 2008-08-28, adapted from 
% ellipsoid_im()
% Copyright 2004-8-13, Patty Laskowsky, Nicole Caparanis, Taka Masuda,
% and Jeff Fessler, The University of Michigan

if nargin == 1 && streq(ig, 'test'), cuboid_im_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

arg.oversample = 1;
arg.type = '';
arg.showmem = false;
arg = vararg_pair(arg, varargin);

if ~isvar('params') || isempty(params)
	params = default_parameters(ig.fov/2, ig.fov/2, ig.zfov/2);
end

switch arg.type
case 'lowmem1'
	fun = @cuboid_im_do_lowmem1;   
case ''
	fun = @cuboid_im_do;
otherwise
	fail('bad type %s', arg.type)
end

[phantom params] = fun(ig.nx, ig.ny, ig.nz, params, ...
	ig.dx, ig.dy, ig.dz, ig.offset_x, ig.offset_y, ig.offset_z, ...
	arg.oversample, arg.showmem);

end % cuboid_im


%
% cuboid_im_do()
%
function [phantom params] = cuboid_im_do(nx, ny, nz, params, ...
	dx, dy, dz, offset_x, offset_y, offset_z, ...
	over, showmem);

if size(params,2) ~= 9
	error 'bad cuboid parameter vector size'
end

phantom = zeros(nx*over, ny*over, nz*over);

wx = (nx*over-1)/2 + offset_x*over;
wy = (ny*over-1)/2 + offset_y*over;
wz = (nz*over-1)/2 + offset_z*over;
xx = ([0:nx*over-1] - wx) * dx / over;
yy = ([0:ny*over-1] - wy) * dy / over;
zz = ([0:nz*over-1] - wz) * dz / over;
[xx yy zz] = ndgrid(xx, yy, zz);

ticker reset
ne = nrow(params);
for ie = 1:ne;
	ticker(mfilename, ie, ne)

	ell = params(ie, :);
	cx = ell(1);	rx = ell(4);
	cy = ell(2);	ry = ell(5);
	cz = ell(3);	rz = ell(6);

	theta = deg2rad(ell(7));
	phi = deg2rad(ell(8));
	if phi, error 'z rotation not done', end
	x = cos(theta) * (xx-cx) + sin(theta) * (yy-cy);
	y = -sin(theta) * (xx-cx) + cos(theta) * (yy-cy);
	z = zz - cz;
    % rx, ry and rz are diameters not radius
    tmp = abs(x/rx) <= 1/2 & abs(y/ry) <= 1/2 & abs(z/rz) <= 1/2;

	phantom = phantom + ell(9) * tmp;
end

if showmem, jf whos, end
phantom = downsample3(phantom, over);
end % ellipsoid_im_do()


%
% cuboid_im_do_lowmem1()
% this version does "one slice at a time" to reduce memory
%
function [phantom params] = cuboid_im_do_lowmem1(nx, ny, nz, params, ...
	dx, dy, dz, offset_x, offset_y, offset_z, over, showmem)

phantom = zeros(nx, ny, nz, 'single');
for iz=1:nz
	offset_z_new = (nz-1)/2 + offset_z - (iz-1);
	phantom(:,:,iz) = cuboid_im_do(nx, ny, 1, params, ...
		dx, dy, dz, offset_x, offset_y, offset_z_new, over, ...
		showmem && iz == 1);
end

end % ellipsoid_im_lowmem1()


%
% default_parameters_parameters()
% for test
%
function params = default_parameters(xfov, yfov, zfov)
params = [...
	0	0	0	0.69	0.92	0.9	0	2.0;
%	0	0	0	0.6624	0.874	0.88	0	-0.98;
%	-0.22	0	-0.25	0.41	0.16	0.21	108	-0.02;
%	0.22	0	-0.25	0.31	0.11	0.22	72	-0.02;
%	0	0.1	-0.25	0.046	0.046	0.046	0	0.02;
%	0	0.1	-0.25	0.046	0.046	0.046	0	0.02;
%	-0.8	-0.65	-0.25	0.046	0.023	0.02	0	0.01;
%	0.06	-0.065	-0.25	0.046	0.023	0.02	90	0.01;
%	0.06	-0.105	0.625	0.56	0.04	0.1	90	0.02;
%	0	0.1	-0.625	0.056	0.056	0.1	0	-0.02];
];
params(:,[1 4]) = params(:,[1 4]) * xfov;
params(:,[2 5]) = params(:,[2 5]) * yfov;
params(:,[3 6]) = params(:,[3 6]) * zfov;
params(:,9) = params(:,8);
params(:,8) = 0; % z rotation
end 


%
% ellipsoid_im_test()
%
function cuboid_im_test
ig = image_geom('nx', 2^5, 'ny', 2^5, 'nz', 2^4, 'fov', 240); % negative dz to match aspire
im pl 2 2
cub = [0 0 0 2*ig.dx 2*ig.dy 2*ig.dz 0 0 1]; 
if 1
	phantom = cuboid_im(ig, cub, 'oversample', 2);
	im(phantom, 'default cuboid'), cbar
end

end % ellipsoid_im_test()
