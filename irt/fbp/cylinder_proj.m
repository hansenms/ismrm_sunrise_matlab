  function proj = cylinder_proj(cg, params, varargin)
%|function proj = cylinder_proj(cg, params, [options])
%|
%| Create set of 2d projection views of one or more (elliptical) cylinders.
%| Works for parallel beam geometry or for flat-detector cone-beam geometry.
%|
%| in
%|	cg			ct_geom()
%|	params [ne 8]		elliptical cylinder parameters:
%|		[centx centy centz  radx rady zlength  angle_degrees  amplitude]
%|
%| options
%|	oversample		oversampling factor for emulating "strips"
%|				(to account for finite detector size)
%|
%| out
%|	proj	[ns nt na]	projections
%|
%| Copyright 2003-10-22, Patty Laskowsky, Nicole Caparanis, Taka Masuda,
%| Ajay Paidi and Jeff Fessler, University of Michigan

if ~nargin, help(mfilename), error(mfilename), end
if nargin == 1 && streq(cg, 'test'), cylinder_proj_test, return, end

% defaults
arg.oversample = 1;
arg = vararg_pair(arg, varargin);

proj = cylinder_proj_do(params, cg.s, cg.t, cg.ar, ...
		cg.dso, cg.dod, cg.dfs, arg.oversample);


%
% cylinder_proj_do()
%
function proj = cylinder_proj_do(params, ss, tt, ...
		beta, ... % radians
		dso, dod, dfs, ...
		oversample)

if size(params, 2) ~= 8, error '8 parameters per cylinder', end 
if oversample ~= 1, error 'only oversample=1 done', end

dsd = dso + dod;

ns = length(ss);
nt = length(tt);

[sss ttt] = ndgrid(ss, tt);

% pvar = sss * dso / dsd; % kak eq 155
% zvar = ttt * dso / dsd;

na = length(beta); % source angles

%
% from cone beam system to an equivalent parallel projection.
%
if isinf(dfs) % flat detector
%	uu = pvar * dso ./ sqrt(dso^2 + pvar.^2); % kak eq 156
%	vv = zvar * dso ./ sqrt(dso^2 + zvar.^2); % kak eq 158

	uu = dso * sss ./ sqrt(dsd^2 + sss.^2); % calculating ray origin (uu,vv,0)
	vv = dso * ttt ./ sqrt(dsd^2 + sss.^2 + ttt.^2) ...
		.* dsd ./ sqrt(dsd^2 + sss.^2);

else
	error 'not done'
end

proj = zeros(ns,nt,na);

% loop over cylinders
for ie = 1:size(params,1)
	ell = params(ie,:);

	cx = ell(1);	rx = ell(4);
	cy = ell(2);	ry = ell(5);
	cz = ell(3);	z_length = ell(6);
	eang = ell(7) / 180 * pi;
	val = ell(8);
	rz = inf; % initializing rz to infinity so that ellipsoid is transformed to infinite cylinder along z

%	gam = atan(zvar/dso); % kak eq 159
	gam = -atan(ttt ./ sqrt(dsd^2 + sss.^2));

	for ib = 1:length(beta)
%		theta = beta(ib) + atan(pvar/dso); % kak eq 158
		theta = beta(ib) + atan(sss / dsd) - eang;

		% shift property of 3D transform:
		ushift = cx * cos(theta) + cy * sin(theta);
		vshift = (cx * sin(theta) - cy * cos(theta)) .* sin(gam) ...
			+ cz * cos(gam);

		% substitute parametric equation of ray (p = p0 + t*p1) in the
		% equation for infinite cylinder (ellipsoid with rz = infinity) 

		p1 = (uu-ushift) .* cos(theta) + (vv-vshift) .* sin(theta) .* sin(gam);
		p2 = (uu-ushift) .* sin(theta) - (vv-vshift) .* cos(theta) .* sin(gam);
		p3 = (vv-vshift) .* cos(gam);

		e1 = -sin(theta) .* cos(gam); % calculating direction cosines
		e2 = cos(theta) .* cos(gam);
		e3 = sin(gam);

		A = e1.^2 / rx^2 + e2.^2 / ry^2;
		B = 2*(p1 .* e1 / rx^2 + p2 .* e2 / ry^2);
		C = p1.^2 / rx^2 + p2.^2 / ry^2 - 1;

	% calculate t for points of intersection of ray along infinite cylinder
	t0 = (-B + (sqrt(B.^2 - 4*A.*C)))./ (2*A);
	t0(find(imag(t0)~=0))= 0; % complex roots mean no intersection with cylinder and hence setting to zero
	t1 = (-B - (sqrt(B.^2 - 4*A.*C)))./ (2*A);
	t1(find(imag(t1) ~= 0)) = 0;

	% finding z values at the points of intersection
	z_alpha = (vv-vshift) .* cos(gam) + t0 .* e3;
	z_beta = (vv-vshift) .* cos(gam) + t1 .* e3;
	% re-arranging z_alpha and z_beta so that z_alpha contains the minimum
	% z values and z_beta contains the maximus z values
	temp = z_alpha;
	z_alpha = min(z_alpha, z_beta);
	z_beta = max(temp, z_beta);

	% re-arranging t0 and t1 according to new z values
	t0_new = (z_alpha - (vv-vshift) .* cos(gam)) ./ e3;
	t1_new = (z_beta - (vv-vshift) .* cos(gam)) ./ e3;
	t0_new(find(imag(t0_new) ~= 0)) = 0;
	t1_new(find(imag(t1_new) ~= 0)) = 0;

	% setting the upper and lower z-bounds on the infinite z cylinder
	t_up = ((z_length/2) - (vv-vshift) .* cos(gam)) ./ e3;
	t_down = ((-z_length/2) - (vv-vshift) .* cos(gam)) ./ e3;
	t_up(find(imag(t_up) ~= 0)) = 0;
	t_down(find(imag(t_down) ~= 0)) = 0;

	% calculating all t-values within this bound
	index_middle = find(z_beta < (z_length/2) & z_beta > (-z_length/2) ...
			& z_alpha < (z_length/2) & z_alpha > (-z_length/2));
	t_alpha = zeros(size(t0)); t_alpha(index_middle) = t0_new(index_middle);
	t_beta = zeros(size(t1)); t_beta(index_middle) = t1_new(index_middle);

	% if z values exceed this bound, then set those z values to the
	% limiting value
	index_up = find(z_beta > (z_length/2) & z_alpha > (-z_length/2) ...
			& z_alpha < (z_length/2));
	t_alpha(index_up) = t0_new(index_up);
	t_beta(index_up) = t_up(index_up);

	index_down = find(z_alpha < (-z_length/2) & z_beta < (z_length/2) ...
			& z_beta > (-z_length/2));
	t_alpha(index_down) = t_down(index_down);
	t_beta(index_down) = t1_new(index_down);

	index_both = find(z_alpha < (-z_length/2) & z_beta > (z_length/2)); % case corresponding to ray entering and leaving
	t_alpha(index_both) = t_down(index_both); % cylinder at diametrically opposite ends
	t_beta(index_both) = t_up(index_both);

	proj(:,:,ib) = proj(:,:,ib) + val * abs(t_alpha - t_beta); % .*(sqrt(A1));
	proj = real(proj);
    end
end


%
% cylinder_proj_test()
% internal test routine
%
function cylinder_proj_test
down = 4;
cg = ct_geom('fan', 'ns', 272, 'nt', 256, 'na', 112, ...
	'ds', 4, 'dsd', 949, 'dod', 408, 'dfs', inf, ...
	'offset_s', 0.25, 'down', down);

param = [ ... % defrise phantom type disks
	[20 0 -82	80 80 20 0 10];
	[20 0 -42	80 80 20 0 10];
	[20 0 +0	80 80 20 0 10];
	[20 0 42	80 80 20 0 10];
	[20 0 82	80 80 20 0 10];
];

oversample = 1;
proj = cylinder_proj(cg, param, 'oversample', oversample);
im clf, im(proj, 'matlab cone-beam projections'), cbar
%im(proj(:,:,24))
