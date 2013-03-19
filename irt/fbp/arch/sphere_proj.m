 function proj = sphere_proj(ells, ss, tt, beta, ...
			dis_src_iso, dis_iso_det, dis_foc_src)
% create set of 2d projections of one or more ellipsoid objects.
% works for parallel beam geometry or for flat-detector cone-beam geometry.
%
% this is now obsolete due to ellipsoid_proj.m
%
% in:
%	ells [ne,8]		ellipsoid parameters:
%		[centx centy centz radx rady radz angle_degrees	amplitude]
%	ss [nh,1]		projection sample locations: horizontal
%	tt [nv,1]		projection sample locations: vertical
%	beta [na,1]		source angles, ccw from vertical. [radians]
%
% optional:
%	dis_src_iso		distance from source to isocenter
%				default is Inf which gives parallel projections
%	dis_iso_det		distance from isocenter to detector
%	dis_foc_src		distance from arc focal point to source
%				default is 0, i.e., 3rd generation x-ray CT
%				use Inf for flat fan-beam detector
%
% out:
%	proj	[nh,nv,na]	projections
%
% Copyright 2003-10-22, Patty Laskowsky, Nicole Caparanis, Taka Masuda
% and Jeff Fessler, The University of Michigan

warning 'sphere_proj is obsolete.  use ellipsoid_proj'

if ~nargin, help(mfilename), error(mfilename), end
if streq(ells, 'test'), sphere_proj_test, return, end

if size(ells, 2) ~= 8, error '8 parameters per ellipsoid', end 

% default is parallel beam
if ~isvar('dis_src_iso') | isempty(dis_src_iso),
	dis_src_iso = inf; dis_iso_det = 1;
end
if ~isvar('dis_iso_det') | isempty(dis_iso_det)
	error 'dis_iso_det is required'
end
if ~isvar('dis_foc_src') | isempty(dis_foc_src)
	dis_foc_src = 0;
end

Dso = dis_src_iso;
Dod = dis_iso_det;
Dsd = Dso + Dod;

nh = length(ss);
nv = length(tt);

[sss ttt] = ndgrid(ss, tt);

%pvar = sss * Dso / Dsd; % kak eq 155
%zvar = ttt * Dso / Dsd;

% source angles
if length(beta) == 1
	na = beta;
	beta = [0:(na-1)]'/na * 2 * pi;	% default 360 degree rotation
else
	na = length(beta);
end

%
% from cone beam system to an equivalent parallel projection.
% (must be a flat detector to work)
%
if isinf(dis_foc_src)	% for a flat detector
%	uu = pvar * Dso ./ sqrt(Dso^2 + pvar.^2); % kak eq 156
%	vv = zvar * Dso ./ sqrt(Dso^2 + zvar.^2); % kak eq 158

	uu = Dso * sss ./ sqrt(Dsd^2 + sss.^2);
	vv = Dso * ttt ./ sqrt(Dsd^2 + sss.^2 + ttt.^2) ...
		.* Dsd ./ sqrt(Dsd^2 + sss.^2);

else
	error 'not done'
end
	
proj = zeros(nh,nv,na);

% loop over spheres
for ie = 1:size(ells,1)
	ell = ells(ie,:);

	cx = ell(1);	rx = ell(4);
	cy = ell(2);	ry = ell(5);
	cz = ell(3);	rz = ell(6);
	if (ry ~= rx | rz ~= rx), error 'only spheres done', end
%	eang = ell(7) / 180 * pi; % needed only for an ellipsoid
	val = ell(8);

%	theta = atan(zvar/Dso); % kak eq 159
	theta = -atan(ttt ./ sqrt(Dsd^2 + sss.^2));

	for ib = 1:length(beta)
%		phi = beta(ib) + atan(pvar/Dso); % kak eq 158, p. 101
		phi = beta(ib) + atan(sss / Dsd);

		% shift property of 3D transform:
		ushift = cx*cos(phi) + cy*sin(phi);
		vshift = (cx*sin(phi) - cy*cos(phi)) .* sin(theta) + cz*cos(theta);
		proj(:,:,ib) = proj(:,:,ib) + 2 * val * ...
			sqrt(rx^2 - (uu - ushift).^2 - (vv - vshift).^2);
	end

	% trick: anywhere proj is imaginary, the real part is 0.
	proj = real(proj);
end


%
% internal test routine
%
function sphere_proj_test
down = 1/8;
nh = 272*down;
nv = 256*down;
na = 200*down;
ell = [1*50 0*50 1*40 200 200 200 0 10];
dis_src_iso = 949.075 - 408.075;
dis_iso_det = 408.075;
dis_foc_src = inf;
ds = 4 / down;
dt = ds; % square detector pixels for now
bin_offset = 0.25;	% quarter detector
ss = ([-(nh-1)/2:(nh-1)/2]-bin_offset) * ds;
tt = ([-(nv-1)/2:(nv-1)/2]) * dt * (-1); % trick: -1 so z goes up

proj = sphere_proj(ell, ss, tt, na, ...
	dis_src_iso, dis_iso_det, dis_foc_src);

im(proj, 'matlab cone-beam projections'), cbar
