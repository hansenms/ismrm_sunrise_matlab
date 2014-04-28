  function proj = cylinder_proj(cg, params, varargin)
%|function proj = cylinder_proj(cg, params, [options])
%|
%| Compute set of 2d line-integral projection views of (elliptical) cylinder(s).
%| Works for these 3D geometries:
%|	parallel beam
%|	flat-detector cone-beam
%|	arc-detector cone-beam (3rd generation CT)
%|	todo: would be nice to have tent cone-beam too!
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
%|	proj	[ns nt na]	projection views
%|
%| Copyright 2003-10-22, Patty Laskowsky, Nicole Caparanis, Taka Masuda,
%| Rebecca Malinas, Ajay Paidi and Jeff Fessler, University of Michigan

if nargin == 1 && streq(cg, 'test'), cylinder_proj_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

% defaults
arg.oversample = 1;
arg = vararg_pair(arg, varargin);

proj = cylinder_proj_do(params, cg.s, cg.t, cg.ar, cg.source_zs, ...
		cg.dso, cg.dod, cg.dfs, arg.oversample);

end % cylinder_proj()


% cylinder_proj_do()
function proj = cylinder_proj_do(params, ss, tt, ...
		beta, ... % [radians]
		source_zs, dso, dod, dfs, oversample)

if size(params, 2) ~= 8, error '8 parameters per cylinder', end 

if oversample > 1
	ds = ss(2) - ss(1);
	dt = tt(2) - tt(1);
	if any(abs(diff(ss) / ds - 1) > 1e-6) ...
	|| any(abs(diff(tt) / dt - 1) > 1e-6)
		error 'uniform spacing required for oversampling'
	end
	No = oversample;
	% determine new finer sampling positions
	ss = outer_sum([-(No-1):2:(No-1)]'/(2*No)*ds, ss(:)'); % [No ns]
	tt = outer_sum([-(No-1):2:(No-1)]'/(2*No)*dt, tt(:)'); % [No nt]
	proj = cylinder_proj_do(params, ss(:), tt(:), beta, source_zs, dso, dod, dfs, 1);
	proj = downsample3(proj, [No No 1]);
return
end


% determine equivalent parallel-beam projection coordinates, at beta=0
ns = length(ss);
nt = length(tt);
[sss ttt] = ndgrid(ss, tt);

if isinf(dso) % parallel beam
	uu = sss;
	vv = ttt;
	azim0 = zeros(size(uu));
	polar = zeros(size(uu));

elseif isinf(dfs) % cone-beam with flat detector
	[uu vv azim0 polar] = ir_coord_cb_flat_to_par(sss, ttt, dso, dod);

elseif dfs == 0 % cone-beam with arc detector
	[uu vv azim0 polar] = ir_coord_cb_arc_to_par(sss, ttt, dso, dod);

else
	fail 'not done'
end

clear sss ttt

cpolar = cos(polar);
spolar = sin(polar);
proj = zeros(ns, nt, numel(beta));

% loop over cylinders
for ip = 1:size(params,1)
	par = params(ip,:);

	cx = par(1);	rx = par(4);
	cy = par(2);	ry = par(5);
	cz = par(3);	zh = par(6) / 2; % half of z_length of cylinder
	eang = deg2rad(par(7)); % xy-plane rotation
	val = par(8);

	for ib = 1:length(beta)
		az = beta(ib) + azim0;

		% shift property of 3D transform:
		cz_eff = cz - source_zs(ib); % center relative to source
		ushift = cx * cos(az) + cy * sin(az);
		vshift = (cx * sin(az) - cy * cos(az)) .* spolar + cz_eff * cpolar;
		az = az - eang;

		% substitute parametric equation of ray (p = p0 + l*p1)
		% in equation for infinite elliptical cylinder

		p1 = (uu-ushift) .* cos(az) + (vv-vshift) .* sin(az) .* spolar;
		p2 = (uu-ushift) .* sin(az) - (vv-vshift) .* cos(az) .* spolar;
		p3 = (vv-vshift) .* cpolar;

		e1 = -sin(az) .* cpolar; % direction cosines of ray
		e2 = cos(az) .* cpolar;
		e3 = spolar;

		A = e1.^2 / rx^2 + e2.^2 / ry^2;
		B = p1 .* e1 / rx^2 + p2 .* e2 / ry^2;
		C = p1.^2 / rx^2 + p2.^2 / ry^2 - 1;

%		bad = @(x) any(isnan(x(:)));

		% calculate l at intersection points of ray with inf. cylinder
		det = B.^2 - A.*C;
		good = det >= 0; % real roots => ray intersects inf. cylinder
		tmp = sqrt(det(good));
		l0 = zeros(size(det));
		l1 = zeros(size(det));
		A = A(good);
		B = B(good);
		l0(good) = (-B - tmp) ./ A;
		l1(good) = (-B + tmp) ./ A;

		% z values at the points of intersection
		z0 = p3 + l0 .* e3;
		z1 = p3 + l1 .* e3;

		% re-arrange so that z0 <= z1
		zswap = z0 > z1;
		[z0 z1] = ir_swap(z0, z1, zswap);

		% re-arrange l0 and l1 according to new z ordering
		[l0 l1] = ir_swap(l0, l1, zswap);

%		if any(l0 > l1), warn('l0 > l1'), end

		t0 = zeros(size(l0));
		t1 = zeros(size(l1));

		% store l-values for cases entirely within finite cylinder
		tmp = (-zh < z0) & (z0 < zh) & (-zh < z1) & (z1 < zh);
pr sum(tmp(:))
		t0(tmp) = l0(tmp);
		t1(tmp) = l1(tmp);

		% set upper and lower z-bounds on the infinite z cylinder
		t_up = (zh - p3) ./ e3; % todo: divide by 0
		t_down = (-zh - p3) ./ e3;

		% if z values exceed this bound, then set those z values to the
		% limiting value
		tmp = (zh < z1) & (-zh < z0) & (z0 < zh); % "index_up"
pr sum(tmp(:))
		t0(tmp) = l0(tmp);
		t1(tmp) = t_up(tmp);

		tmp = (z0 < -zh) & (-zh < z1) & (z1 < zh); % "index_down"
pr sum(tmp(:))
		t0(tmp) = t_down(tmp);
		t1(tmp) = l1(tmp);

 		% ray entering and leaving cylinder at opposite ends:
		tmp = (z0 < -zh) & (zh < z1); % "index_both"
pr sum(tmp(:))
		t0(tmp) = t_down(tmp);
		t1(tmp) = t_up(tmp);

		proj(:,:,ib) = proj(:,:,ib) + val * abs(t1 - t0);
	end % ib
end % ip

end % cylinder_proj()


% ir_swap()
% swap specified elements
function [o0 o1] = ir_swap(i0, i1, swap)
o0 = i0;
o1 = i1;
o0(swap) = i1(swap);
o1(swap) = i0(swap);

end % ir_swap()


% cylinder_proj_test()
% internal test routine
function cylinder_proj_test

param = [ ... % defrise phantom type disks
	[20 0 -82	80 80 20 0 10];
	[20 0 -42	80 80 20 0 10];
	[20 0 +0	80 80 20 0 10];
	[20 0 42	80 80 20 0 10];
	[20 0 82	80 80 20 0 10];
];

param = [20 10 -10	80 70 160 0 10]; % todo: quick simple cylinder test
param = [20 10 -90	180 70 80 0 10]; % todo: quick simple cylinder test

fun_proj = @(cg) cylinder_proj(cg, param, 'oversample', 1); % analytical
fun_im = @(ig) cylinder_im(ig, param, 'oversample', 2, 'checkfov', true);

ir_proj3_compare1(fun_proj, fun_im, 'chat', 1, 'dsd', 949, 'dfs', 0, ...
	'arg_ct_geom', {'na', 1}, ...
	'nt', 65, ... % todo stress odd sin(polar) = 0
	'source_z0', 7, 'pitch', 0);

return

ir_proj3_compare1(fun_proj, fun_im, 'downi', 4, ...
	'chat', 1, 'dsd', inf, 'dfs', 0, 'pitch', 0);
ir_proj3_compare1(fun_proj, fun_im, 'chat', 1);

if 0
end

down = 4;
oversample = 1;

im plc 2 2

dfs_list = [0 inf];
for ii = 1:numel(dfs_list)
	dfs = dfs_list(ii);
	cg = ct_geom('fan', 'ns', 272, 'nt', 256, 'na', 112, ...
		'ds', 4, 'dsd', 949, 'dod', 408, 'dfs', dfs, ...
		'offset_s', 0.25, 'down', down);

	proj = cylinder_proj(cg, param, 'oversample', oversample);
	im(ii, proj, 'matlab cone-beam projections'), cbar
	titlef('matlab cone-beam projections, Dfs=%g', dfs)
	im(2+ii, proj(:,:,end/2))
end
%im(proj(:,:,24))

end % cylinder_proj_test()
