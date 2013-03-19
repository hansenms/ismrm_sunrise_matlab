  function proj = ellipsoid_proj(cg, ells, varargin)
%|function proj = ellipsoid_proj(cg, ells, varargin)
%|
%| Compute a set of 2d line-integral projection views of one or more ellipsoids.
%| Works for both parallel-beam and cone-beam geometry.
%|
%| in
%|	cg		ct_geom()
%|	ells [ne 9]	ellipsoid parameters:
%|			[x_center y_center z_center  x_radius y_radius z_radius
%|				xy_angle_degrees z_angle_degrees  amplitude]
%| options
%|	oversample	over-sampling factor (approximates finite detector size)
%|
%| out
%|	proj	[ns nt na]	projection views
%|
%| Copyright 2003-10-22, Patty Laskowsky, Nicole Caparanis, Taka Masuda,
%| and Jeff Fessler, University of Michigan

if nargin == 1 && streq(cg, 'test'), ellipsoid_proj_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

arg.oversample = 1;
arg = vararg_pair(arg, varargin);

proj = ellipsoid_proj_do(ells, cg.s, cg.t, cg.ar, cg.source_zs, ...
		cg.dso, cg.dod, cg.dfs, arg.oversample);

end % ellipsoid_proj()


%
% ellipsoid_proj_do()
%
function proj = ellipsoid_proj_do(ells, ss, tt, ...
		beta, ... % [radians]
		source_zs, dso, dod, dfs, oversample)

if size(ells, 2) ~= 9, error '9 parameters per ellipsoid', end

if oversample > 1
	ds = ss(2) - ss(1);
	dt = tt(2) - tt(1);
	if any(abs(diff(ss) / ds - 1) > 1e-6) ...
	|| any(abs(diff(tt) / dt - 1) > 1e-6) ...
		error 'uniform spacing required for oversampling'
	end
	No = oversample;
	% determine new finer sampling positions
	ss = outer_sum([-(No-1):2:(No-1)]'/(2*No)*ds, ss(:)'); % [No ns]
	tt = outer_sum([-(No-1):2:(No-1)]'/(2*No)*dt, tt(:)'); % [No nt]
	proj = ellipsoid_proj_do(ells, ss(:), tt(:), beta, source_zs, dso, dod, dfs, 1);
	proj = downsample3(proj, [No No 1]);
return
end

Ds = dso;
Dd = dod;
Dc = Ds + Dd;

%
% determine equivalent (u,v; phi) parallel projection coordinates, at beta=0.
%

ns = length(ss);
nt = length(tt);
[sss ttt] = ndgrid(ss, tt);

if isinf(dso) % parallel beam
	uu = sss;
	vv = ttt;
	phi0 = 0;
	theta = 0;
pr dso


elseif isinf(dfs) % cone-beam with flat detector

%	pvar = sss * Ds / Dc; % kak eq 155
%	zvar = ttt * Ds / Dc;
%	uu = pvar * Ds ./ sqrt(Ds^2 + pvar.^2); % kak eq 156
%	vv = zvar * Ds ./ sqrt(Ds^2 + zvar.^2); % kak eq 158

	uu = Ds * sss ./ sqrt(Dc^2 + sss.^2);
	vv = Ds * ttt ./ sqrt(Dc^2 + sss.^2 + ttt.^2) ...
		.* Dc ./ sqrt(Dc^2 + sss.^2);

%	theta = atan(zvar/Ds); % kak eq 159
	theta = -atan(ttt ./ sqrt(Dc^2 + sss.^2)); % trick: empirical negative

%	phi0 = atan(pvar/Ds); % kak eq 158
	phi0 = atan(sss / Dc);


else % cone-beam with arc detector
	pos_src = [0 dso 0];
	Rf = dfs + Dc; % focal radius
	pos_det(:,:,1) = Rf * sin(sss / Rf);
	pos_det(:,:,2) = Ds + dfs - Rf * cos(sss / Rf);
	pos_det(:,:,3) = ttt;

	% ray between source and detector element center
	ee1 = pos_det(:,:,1) - pos_src(1);
	ee2 = pos_det(:,:,2) - pos_src(2);
	ee3 = pos_det(:,:,3) - pos_src(3);
	enorm = sqrt(ee1.^2 + ee2.^2 + ee3.^2);
	ee1 = ee1 ./ enorm;
	ee2 = ee2 ./ enorm;
	ee3 = ee3 ./ enorm;

	theta = -asin(ee3); % trick: empirical negative
	phi0 = -atan(ee1 ./ ee2);
	uu = cos(phi0) * pos_src(1) + sin(phi0) * pos_src(2);
	vv = (sin(phi0) * pos_src(1) - cos(phi0) * pos_src(2)) .* sin(theta) ...
		+ pos_src(3) * cos(theta);
end

clear sss ttt

cthet = cos(theta);
sthet = sin(theta);
proj = zeros(ns, nt, length(beta));

%
% loop over ellipsoids
%
for ie = 1:size(ells,1)
	ell = ells(ie,:);

	cx = ell(1);	rx = ell(4);
	cy = ell(2);	ry = ell(5);
	cz = ell(3);	rz = ell(6);
	eang = deg2rad(ell(7)); % xy-plane rotation of ellipsoid
	zang = deg2rad(ell(8)); % z-plane rotation of ellipsoid
	if zang, error 'z rotation not done', end
	val = ell(9);

	for ib = 1:length(beta)
%		phi = beta(ib) + phi0 - eang; % assume source rotate in xy plane
		phi = beta(ib) + phi0; % correction due to Lei Zhu of Stanford

		% shift property of 3D transform:
		cz_eff = cz - source_zs(ib); % center relative to source
		ushift = cx*cos(phi) + cy*sin(phi);
		vshift = (cx*sin(phi) - cy*cos(phi)) .* sthet + cz_eff * cthet;
		phi = phi - eang; % correction due to Lei Zhu of Stanford

		p1 = (uu-ushift) .* cos(phi) + (vv-vshift) .* sin(phi) .* sthet;
		p2 = (uu-ushift) .* sin(phi) - (vv-vshift) .* cos(phi) .* sthet;
		p3 = (vv-vshift) .* cthet;

		e1 = -sin(phi) .* cthet;
		e2 = cos(phi) .* cthet;
		e3 = sthet;

		A = e1.^2 / rx^2 + e2.^2 / ry^2 + e3.^2 / rz^2;
		B = 2 * (p1.*e1 / rx^2 + p2.*e2 / ry^2 + p3.*e3 / rz^2);
		C = p1.^2 / rx^2 + p2.^2 / ry^2 + p3.^2 / rz^2 - 1;

		proj(:,:,ib) = proj(:,:,ib) + val * sqrt(B.^2 - 4 * A.*C) ./ A;
	end
end

% trick: anywhere proj of a single ellipsoid is imaginary, the real part is 0.
proj = real(proj);
end % ellipsoid_proj_do()


%
% ellipsoid_proj_test_pitch()
% internal test routine
%
function ellipsoid_proj_test_pitch(pitch, varargin)
down = 30;
cg = ct_geom('fan', 'ns', round(888/down), 'nt', 64, ...
	'na', 12, ...
	'ds', 1.0*down, 'dt', 1.1*down, ...
	'down', 1, ... % only downsample s,t, beta
	'offset_s', 0.25, ... % quarter detector
	'offset_t', 0.0, ...
	'pitch', pitch, 'source_z0', -200, ... % test helix too
	varargin{:});

ell = [20 0*50 -40 200 100 50 90 0 10;
	0 50 100 80 80 20 0 0 10];

pe = ellipsoid_proj(cg, ell, 'oversample', 2); % analytical

if im && 0
	t = sprintf('matlab cone-beam projections, dfs=%g', cg.dfs);
	clf, im(cg.s, cg.t, pe, t), cbar
	xlabel s, ylabel t
prompt
end

if 1
	downi = 4;
	ig = image_geom('nx', 512/downi, 'nz', 256/downi, ...
		'offset_x', -9, 'offset_y', 7, 'offset_z', -5, ... % stress test
		'dx', 1*downi, 'dz', 1*downi, 'down', 1);
	x = ellipsoid_im(ig, ell, 'oversample', 2, 'checkfov', true);
	im plc 2 2
	im(1, x)
	im(2, pe), cbar
end

% compare with Gcone
if 1 && has_mex_jf, printm 'compare to Gcone discrete projector'
	A = Gcone(cg, ig, 'type', 'sf1', 'use_hct2', 0);

	pc = A * x;
	max_percent_diff(pe, pc)

	if im
		im(3, pc), cbar, xlabel(A.arg.type)
		im(4, pc - pe, 'pc - pe'), cbar
	prompt
	end
end

% compare with aspire
if has_mex_jf && isinf(cg.dfs) && ~isinf(cg.dsd) % only flat cone-beam
	printm 'compare to 3l@ discrete projector'
	systype = aspire_pair(cg, ig, 'system', '3l');
	A = Gtomo3(systype, ig.mask, ig.nx, ig.ny, ig.nz, ...
		'chat', 0, 'checkmask', false, 'permute213', true);

	p3 = A * x;
	max_percent_diff(pe, p3)

	if im
		im plc 2 2
		im(1, x)
		im(2, pe), cbar
		im(3, p3), cbar
		im(4, p3 - pe, 'p3 - pe'), cbar
	prompt
	end

	if im && 0
		im_toggle(pe, p3, [0 max(pe(:))])
		clf, im(p3-pe)
		prompt
%		im clf, movie2(pe)
%	prompt
	end
end

end % ellipsoid_proj_test_pitch()


%
% ellipsoid_proj_test()
%
function ellipsoid_proj_test
dsd_list = [949 949 inf]; % cone-arc, cone-flat, parallel-beam
dfs_list = [0 inf 0];
for ii=1:length(dsd_list)
	dsd = dsd_list(ii);
	dfs = dfs_list(ii);
	for pitch = [0.5 0]
		pr pitch
		ellipsoid_proj_test_pitch(pitch, ...
			'dsd', dsd, 'dod', 408, 'dfs', dfs);
	end
end

end % ellipsoid_proj_test()
