 
function proj = cuboid_proj(cg, cubs, varargin)
%function proj = cuboid_proj(cg, cubs, varargin)
%
% Compute a set of 2d line-integral projection views of one or more cuboids.
% Works for both parallel-beam and cone-beam geometry.
%
% in
%	cg		ct_geom()
%	ells [ne,9]	cuboid parameters:
%			[x_center y_center z_center  x_radius y_radius z_radius
%				xy_angle_degrees z_angle_degrees  amplitude]
% options
%	oversample	over-sampling factor (approximates finite detector size)
%
% out
%	proj	[ns,nt,na]	projection views
%
% Yong Long, 2008-08-28, adapted from
% ellipseoid_proj()
% Copyright 2003-10-22, Patty Laskowsky, Nicole Caparanis, Taka Masuda,
% and Jeff Fessler, The University of Michigan

if nargin == 1 && streq(cg, 'test'), cuboid_proj_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

arg.oversample = 1;
arg = vararg_pair(arg, varargin);

%if ~isempty(cg.zshifts), error 'zshifts not done', end
proj = cuboid_proj_do(cubs, cg.s, cg.t, cg.ar, ...
		cg.dso, cg.dod, cg.dfs, arg.oversample);

end % cuboid_proj()


%
% cuboid_proj_do()
%
function proj = cuboid_proj_do(cubs, ss, tt, ...
		beta, ... % [radians]
		dso, dod, dfs, oversample)

if size(cubs, 2) ~= 9, error '9 parameters per cuboid', end

if oversample > 1
	ds = ss(2) - ss(1);
	dt = tt(2) - tt(1);
	if any(abs(diff(ss) / ds - 1) > 1e-6) ...
	|| any(abs(diff(tt) / dt - 1) > 1e-6) ...
		error 'uniform spacing required for oversampling'
	end
	No1 = oversample;
    No2 = oversample;
	% determine new finer sampling positions
	ss = outer_sum([-(No1-1):2:(No1-1)]'/(2*No1)*ds, ss(:)'); % [No,ns]
	tt = outer_sum([-(No2-1):2:(No2-1)]'/(2*No2)*dt, tt(:)'); % [No,ns]
	proj = cuboid_proj_do(cubs, ss(:), tt(:), beta, dso, dod, dfs, 1);
	if (ndims(proj) == 3)
        proj = downsample3(proj, [No1 No2 1]);
    elseif (ndims(proj) == 2)
        proj = downsample2(proj, [No1 No2]);
    end
return
end

%Ds = dso;
%Dd = dod;
%Dc = Ds + Dd;
dsd = dso + dod;
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


elseif isinf(dfs) % cone-beam with flat detector

	uu = dso * sss ./ sqrt(dsd^2 + sss.^2);
	vv = dso * ttt ./ sqrt(dsd^2 + sss.^2 + ttt.^2) ...
		.* dsd ./ sqrt(dsd^2 + sss.^2);

	theta = -atan(ttt ./ sqrt(dsd^2 + sss.^2)); % trick: empirical negative

	phi0 = atan(sss / dsd);


else % cone-beam with arc detector
	pos_src = [0 dso 0];
	Rf = dfs + dsd; % focal radius
	pos_det(:,:,1) = Rf * sin(sss / Rf);
	pos_det(:,:,2) = dso + dfs - Rf * cos(sss / Rf);
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
clear pos_det ee1 ee2 ee3 enorm

cthet = cos(theta);
sthet = sin(theta);
proj = zeros(ns, nt, length(beta));

%
% loop over cuboids
%
for ie = 1:size(cubs,1)
	cub = cubs(ie,:);

	cx = cub(1);	rx = cub(4);
	cy = cub(2);	ry = cub(5);
	cz = cub(3);	rz = cub(6);
	eang = deg2rad(cub(7)); % xy-plane rotation of cuboid
	zang = deg2rad(cub(8)); % z-plane rotation of cuboid
	if zang, error 'z rotation not done', end
	val = cub(9);

	for ib = 1:length(beta)
%		phi = beta(ib) + phi0 - eang; % assume source rotate in xy plane
		phi = beta(ib) + phi0; % correction due to Lei Zhu of Stanford

		% shift property of 3D transform:
		ushift = cx*cos(phi) + cy*sin(phi);
		vshift = (cx*sin(phi) - cy*cos(phi)) .* sthet + cz * cthet;
		phi = phi - eang; % correction due to Lei Zhu of Stanford

		p1 = (uu-ushift) .* cos(phi) + (vv-vshift) .* sin(phi) .* sthet;
		p2 = (uu-ushift) .* sin(phi) - (vv-vshift) .* cos(phi) .* sthet;
		p3 = (vv-vshift) .* cthet;

		e1 = -sin(phi) .* cthet; %  x = p1 + l*e1
		e2 = cos(phi) .* cthet; % y = p2 + l*e2
		e3 = sthet; % z = p3 + l*e3
        
        clear phi
        
        % bounds of l corresponding to rect(x/rx)
        lxmin = (- rx/2 - p1) ./ e1;
        lxmax = (rx/2 - p1) ./ e1;
        % re-arrange the bounds so that lxmin contains the minimum l values
        % and lxmax contains the maximum l values
        temp = lxmin;
        lxmin = min(lxmin, lxmax);
        lxmax = max(temp, lxmax);
        % exclude the points where e1=0 by setting lxmin = -Inf
        % and lxmax = Inf
        lxmin(find(lxmin == inf)) = -inf;
        lxmax(find(lxmax == -inf)) = inf;
        clear p1
        % bounds of l corresponding to rect(y/ry)
        lymin = (- ry/2 - p2) ./ e2;
        lymax = (ry/2 - p2) ./ e2;
        % re-arrange the bounds so that lymin contains the minimum l values
        % and lymax contains the maximum l values
        temp = lymin;
        lymin = min(lymin, lymax);
        lymax = max(temp, lymax);
        lymin(find(lymin == inf)) = -inf;
        lymax(find(lymax == -inf)) = inf;
        clear p2
        % bounds of l corresponding to rect(z/rz)
        lzmin = (- rz/2 - p3) ./ e3;
        lzmax = (rz/2 - p3) ./ e3;
        % re-arrange the bounds so that lzmin contains the minimum l values
        % and lzmax contains the maximum l values
        temp = lzmin;
        lzmin = min(lzmin, lzmax);
        lzmax = max(temp, lzmax);
        lzmin(find(lzmin == inf)) = -inf;
        lzmax(find(lzmax == -inf)) = inf;
        clear p3
        % lower bound for l
        lmin = max(lxmin, lymin);
        lmin = max(lmin, lzmin);
        
         % upper bound for l
        lmax = min(lxmax, lymax);
        lmax = min(lmax, lzmax);
        
        l = max(lmax - lmin, 0); % intersection only when lmax > lmin
        
        % three cases where e(k) = 0.
        if (e3 == 0)% sin(theta) = 0; 
            zero_e = (-rz/2 <= (vv-vshift)) & ((vv-vshift) <= rz/2); 
            proj_i = val * l .* zero_e; 
        else
            proj_i = val * l;
        end
              
        if (e1 == 0)% sin(phi) = 0; theta samll, cos(theta) will not be zero
            zero_e = (-rx/2 <= (uu-ushift)) & ((uu-ushift) <= rx/2); 
            proj_i = proj_i .* zero_e; 
        end
       
        if (e2 == 0)% cos(phi) = 0; 
            zero_e = (-ry/2 <= (uu-ushift)) & ((uu-ushift) <= ry/2); 
            proj_i = proj_i .* zero_e; 
        end
        proj(:,:,ib) = proj(:,:,ib) + proj_i;
        clear proj_i
	end
end

% trick: anywhere proj of a single cuboid is imaginary, the real part is 0.
%proj = real(proj);
end % cuboid_proj_do()


%
% cuboid_proj_test()
% internal test routine
%
function cuboid_proj_test
down = 30;
cg = ct_geom('fan', 'ns', 888/down, 'nt', 64, 'na', 984/down, ...
 'ds', 1.0*down, 'down', 1, ... % only downsample s and beta
 'offset_s', 0.25, ... % quarter detector
 'offset_t', 0.0, ...
 'dsd', 949, 'dod', 408, 'dfs', inf); % flat detector
% 'dsd', inf); % parallel-beam
% 'dsd', 949, 'dod', 408, 'dfs', 0); % 3rd gen CT
%cg.rmax

cubs = [0*50 0*50 0*50 200 100 100 90 0 10];

proj = cuboid_proj(cg, cubs, 'oversample', 2);

t = sprintf('matlab cone-beam projections, dfs=%g', cg.dfs);
im clf, im(cg.s, cg.t, proj, t), cbar
xlabel s, ylabel t
%im clf, im(cg.s, cg.ad, permute(proj, [1 3 2])), cbar
%xlabel s, ylabel '\beta'
end % cuboid_proj_test()
