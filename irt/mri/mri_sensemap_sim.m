  function smap = mri_sensemap_sim(varargin)
%|function smap = mri_sensemap_sim(varargin)
%|
%| Simulate sensitivity maps for sensitivity-encoded MRI
%| based grivich:00:tmf doi:10.1119/1.19461
%|
%| option
%|	nx, ny, dx, dy, ncoil, rcoil, orbit (see below)
%| out
%|	smap	[nx ny ncoil]	simulated sensitivity maps (complex!)
%|
%| Copyright 2005-6-20, Jeff Fessler and Amanda Funai, University of Michigan

if nargin == 1 && streq(varargin{1}, 'test'), mri_sensemap_sim_test, return, end

arg.nx = 64;
arg.ny = [];
arg.dx = 3; % pixel size in mm
arg.dy = [];
arg.ncoil = 4; % # of coils
arg.rcoil = 100; % coil radius
arg.orbit = 360;
arg.coil_distance = 1.2; % multiplies fov/2
arg.flag_old = false; % old way based on wang:00:dop
arg.chat = nargout == 0;

arg = vararg_pair(arg, varargin);

if isempty(arg.dy), arg.dy = arg.dx; end
if isempty(arg.ny), arg.ny = arg.nx; end
if isempty(arg.rcoil), arg.rcoil = arg.dx * arg.nx / 2 * 0.50; end

smap = mri_sensemap_sim_do(arg.nx, arg.ny, arg.dx, arg.dy, ...
	arg.ncoil, arg.rcoil, arg.orbit, arg.coil_distance, ...
	arg.flag_old, arg.chat);

if ~nargout
	clear smap
end

%
% mri_sensemap_sim_do()
%
function smap = mri_sensemap_sim_do(nx, ny, dx, dy, ncoil, rcoil, orbit, ...
		coil_distance, flag_old, chat)

rlist = rcoil * ones(ncoil, 1); % coil radii

plist = zeros(ncoil,3); % position of coil center
nlist = zeros(ncoil,3); % normal vector (inward) from coil center

% circular coil configuration, like head coils
alist = deg2rad(orbit)/ncoil * [0:ncoil-1]; % list of angles in radians
for ii=1:ncoil
	phi = alist(ii);
	Rad = max(nx/2 * dx, ny/2 * dy) * coil_distance;
	plist(ii,:) = Rad * [cos(phi) sin(phi) 0];
	nlist(ii,:) = -[cos(phi) sin(phi) 0];
	olist(ii,:) = [-sin(phi) cos(phi) 0]; % unit vector orthogonal to nlist
end

% object coordinates for slice z=0
x = ([1:nx]-(nx+1)/2)*dx;
y = ([1:ny]-(ny+1)/2)*dy;
z = 0;
[xx yy zz] = ndgrid(x,y,z);

smap = zeros(nx, ny, ncoil);
for ii=1:ncoil
	% rotate coordinates to correspond to coil orientation
	zr =	(xx - plist(ii,1)) .* nlist(ii,1) + ...
		(yy - plist(ii,2)) .* nlist(ii,2) + ...
		(zz - plist(ii,3)) .* nlist(ii,3);
	xr =	xx .* nlist(ii,2) - yy .* nlist(ii,1);
	[sx sy sz] = mri_smap1(xr, 0, zr, rlist(ii)); % in coil coordinates

	if flag_old
		smap(:,:,ii) = sz; % old way (wrong!) based on smap_z

	else
		if nlist(ii,3) || olist(ii,3)
			fail 'unsupported'
		end
		% assume z component of plist and nlist are 0
		bx = sz * nlist(ii,1) + sx * olist(ii,1);
		by = sz * nlist(ii,2) + sx * olist(ii,2);
		smap(:,:,ii) = bx + 1i * by;
	end
end
smap = smap * rlist(1) / (2*pi); % trick: scale so near unity maximum

% plot array geometry in z=0 plane
if chat && im
	nshow = min(max(ncoil,2),4);
	im('plc', 3, nshow)
	for ii=1:min(ncoil,nshow)
		clim = [0 max(abs(smap(:)))];
		tmp = smap(:,:,ii);
		im(ii, x, y, abs(tmp), clim, 'Magnitude'), cbar
%		xmax = max(max(abs(x)), max(plist(:,1)));
%		ymax = max(max(abs(y)), max(plist(:,2)));
		xmax = max([max(abs(x)) max(abs(y)) max(col(plist(:,[1 2])))]);
%		axis([-xmax xmax -ymax ymax]*1.05)
		axis(xmax*[-1 1 -1 1]*1.1)

		hold on
		plot(0,0,'.', plist(:,1), plist(:,2), 'o')
		xdir = nlist(ii,2);
		ydir = nlist(ii,1);
		r = rlist(ii);
		plot(plist(ii,1)+r*xdir*[-1 1], plist(ii,2)+r*ydir*[1 -1], '-')
		hold off

		% trick: unwrap phase for pretty display
		ph = unwrap(angle(tmp).').';
		ph = unwrap(ph);
		im(ii+nshow, x, y, ph, 'Phase'), cbar
		axis(xmax*[-1 1 -1 1]*1.1)
	end
	sos = sum(abs(smap), 3);
	sos = sos / sos(end/2,end/2);
	im(2*nshow+1, x, y, sos, 'SoS (normalized)'), cbar

	if ncoil == 1
		subplot(212)
%		clf
		bx = real(smap);
		by = imag(smap);
		quiver(xx, yy, bx, by), title 'Field pattern in x-y plane'
		axis equal, axis tight
	end

%	clf, im(angle(smap(:,:,2))), cbar
end

if ~nargout, clear smap, end


%
% mri_smap_r(r, z)
% function for testing near 0
%
function out = mri_smap_r(r, z)
M = 4 * r ./ ((1 + r).^2 + z.^2); % = k^2, see ellipke
[K E] = ellipke(M);
out = 2 * z ./ r .* ((1+r).^2 + z.^2).^(-0.5) .* ...
	((1 + r.^2 + z.^2) ./ ((1-r).^2 + z.^2) .* E - K);


%
% based on grivich:00:tmf
% for a circular coil in "x-y plane" of radius a
% note that coil x-y plane is not same as object x-y plane!
%
function [smap_x smap_y smap_z] = mri_smap1(x, y, z, a)
x = x ./ a;
y = y ./ a;
z = z ./ a;
r = sqrt(x.^2 + y.^2);
M = 4 * r ./ ((1 + r).^2 + z.^2); % = k^2, see ellipke
[K E] = ellipke(M);

% the following is B_z in eqn (18) in grivich:00:tmf
smap_z = 2 * ((1+r).^2 + z.^2).^(-0.5) .* ...
	(K + (1 - r.^2 - z.^2) ./ ((1-r).^2 + z.^2) .* E);
smap_z = smap_z / a;

if 0 && any(r(:) == 0) % test code to explore when r is near 0
	r0 = linspace(0,5e-7,101);
	z0 = 0.4;
	t0 = mri_smap_r(r0, z0);
	slope = 3*pi * z0 / ((1+z0^2)^2.5);
	clf, plot(r0, t0, '-', r0, slope * r0, '--'); grid, prompt
end

% the following is B_r in eqn (17) in grivich:00:tmf
smap_r = 2 * z ./ r .* ((1+r).^2 + z.^2).^(-0.5) .* ...
	((1 + r.^2 + z.^2) ./ ((1-r).^2 + z.^2) .* E - K);
bad = abs(r) < 1e-6;
smap_r(bad) = 3 * pi * z(bad) ./ ((1 + z(bad).^2).^2.5) .* r(bad);
smap_r = smap_r / a;

if any(isnan(smap_r(:))) || any(isnan(smap_z(:)))
	keyboard
end

phi = atan2(y, x);
%smap_x = smap_r .* x ./ r; % divide by 0 problems
%smap_y = smap_r .* y ./ r;
smap_x = smap_r .* cos(phi);
smap_y = smap_r .* sin(phi);


%
% mri_sensemap_sim_test
%
function mri_sensemap_sim_test
smap = mri_sensemap_sim('chat', 1, ...
	'rcoil', [], 'ncoil', 4, 'coil_distance', 1.2);
%, 'nx', 32, 'dx', 6);

if 0
	m = linspace(0,1,201);
	m = 0;
	[k e] = ellipke(m)
	plot(m, k, '-', m, e, '--')
end

if 0
	a = 1;
	x = linspace(-2,2,102);
	y = linspace(-2,2,104);
	z = [0.1 0.2 0.5 1.0];
	[xx yy zz] = ndgrid(x,y,z);
	[smap_x smap_y smap_z] = mri_smap1(xx, yy, zz, a);
	im pl 4 1
	im row 1
	im(1, smap_x, 'x'), cbar
	im(2, smap_y, 'y'), cbar
	im(3, smap_z, 'z'), cbar
	im(4, sqrt(smap_x.^2+smap_y.^2), 'r'), cbar
end
