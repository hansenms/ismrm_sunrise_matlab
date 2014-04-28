% mri_pixel_size_example.m
% Show the effects of small pixels on iterative reconstruction for MRI.
% Copyright 2003-7-23, Jeff Fessler, The University of Michigan

%
% create k-space samples, in units of 1/mm
%
if ~isvar('kspace'), printm 'setup kspace'
	fov = 256;	% [mm] typical brain FOV
	N0 = 64;	% nominal image size
	kmax = N0/2*(1/fov);	% for display axes

	ktype = 'cartesian';
	ktype = 'random';
	ktype = 'spiral';

	if streq(ktype, 'random')
		rng(0)
		kspace = (rand(N0*N0, 2)-0.5)*N0/fov;	% random k-space
	elseif streq(ktype, 'spiral')
		t = linspace(0, N0/2*2*pi, N0^2)';	% crude spiral:
		kspace = N0/2*(1/fov)*[cos(t) sin(t)] .* (t(:,[1 1]) / max(t));
	elseif streq(ktype, 'cartesian')
		k1 = [-N0/2:N0/2-1]/fov;		% cartesian
		[kk1, kk2] = ndgrid(k1, k1);
		kspace = [kk1(:), kk2(:)];
		1/diff(kspace(1:2,1))
	end, clear t k1 kk1 kk2

	im clf, im pl 2 3
	if im
		im subplot 1
		if 0
			plot(kspace(:,1), kspace(:,2), '.')
			axis(1.1*[-1 1 -1 1]*kmax), axis square
			xlabel 'k_1 [mm^{-1}]', ylabel 'k_2 [mm^{-1}]'
		else
			text(0.5, 0.5, '(not shown)', 'horiz', 'center')
		end
		title(sprintf('%d %s k-space samples', size(kspace,1), ktype))
	end
prompt
end

%
% true object and analytical k-space data
%
if 0 || ~isvar('xtrue'), printm 'setup object'

	% display images with many pixels...
	Ndisp = 512;
	x1d = [-Ndisp/2:Ndisp/2-1] / Ndisp * fov;
	[x1dd, x2dd] = ndgrid(x1d, x1d);

	% parameter units all in [mm]
	obj = mri_objects('case1');
	xtrue = obj.image(x1dd, x2dd);
	ytrue = obj.kspace(kspace(:,1), kspace(:,2));
	clear x1dd x2dd obj

	clim = [0 2];
	im(2, x1d, x1d, xtrue, 'x true', clim), cbar

	% add noise
	rng(0)
	yd = ytrue + 0 * randn(size(ytrue));

	% lazy kspace data gridding for display
	Ng = 1*N0;
	[xhatg yd_g xg kg] = mri_grid_linear(kspace, yd, Ng, fov);

	disp(imax(yd_g, 2))
	im(3, kg{1}, kg{2}, abs(yd_g), '|y_d|'), cbar

	if streq(ktype, 'cartesian') && im
		im subplot 6
		k1 = kspace(1:N0, 1);
		im(6, k1, k1, abs(reshape(ytrue, 64, 64)), 'reshape'), cbar
		im(5, k1, k1, abs(reshape(ytrue, 64, 64)-yd_g), 'diff'), cbar
	end, clear k1

	im(4, xg{1}, xg{2}, abs(xhatg), '|x| "gridding"', clim), cbar
prompt
end, clear Ng yd_g xg kg

%
% now loop over several pixel sizes
%
if ~isvar('xpcg')

Nlist = [2.^[5:9]];
niter = 20;
beta = 2^-10 * size(kspace,1);	% good for quadratic
cost = zeros(niter, length(Nlist));
xpcg = {};

% loop over image sizes
for in=1:length(Nlist)
	N = Nlist(in);

	% gridding estimate to initialize iterations
	[xhatg, yhatg, xg] = mri_grid_linear(kspace, yd, N, fov);
	tmp = sprintf('|x| gridding %d', N);
	im(5, xg{1}, xg{2}, abs(xhatg), tmp, clim), cbar

	% create Gnufft class object
	if 1 || ~isvar('G'), printf('setup Gnufft object for N=%d', N)
		omega = 2*pi*kspace*fov/N;
%		minmax(omega(:,1))
		G = Gnufft({omega, [N N], [6 6], 2*[N N], [N/2 N/2], ...
			'table', 2^11, 'kaiser'});
	end

	% reconstruct by PCG
	if 1 || ~isvar('xpcg'), disp 'PCG with quadratic penalty'
		mask = true(N);
		R = Robject(mask, 'beta', beta);
		ytmp = yd(:) * (1/fov)^2 * N^2;	% scaling!
		xiter = qpwls_pcg(xhatg(:), G, 1, ytmp, 0, R.C, ...
				1, niter);
		cost(:,in) = pwls_cost(xiter, G, 1, ytmp, R);
		xpcg{in} = embed(xiter(:,end), mask);
		im(6, xg{1}, xg{2}, abs(xpcg{in}), '|x| pcg quad'), cbar
	end
end % for
end % if

if im
	im clf, im pl 2 3
	im(1, x1d, x1d, xtrue, clim, 'x true'), cbar
	xtick([-fov/2 0 fov/2-1])
	ytick([-fov/2 0 fov/2-1])
	for in=1:length(Nlist)
		N = Nlist(in);
		x1g = [-N/2:N/2-1]'/N * fov;
		tmp = sprintf('N=%d', Nlist(in));
		im(in+1, x1g, x1g, abs(xpcg{in}), tmp, clim), cbar
		xtick([-fov/2 0 fov/2-1])
		ytick([-fov/2 0 fov/2-1])
	end
%	savefig fig_mri_pixel_size_image
prompt
end

% cost figure: not much to learn here.
if 0 && im
	clf
	tmp = 1 ./ cost(1,:);
	tmp = cost * diag(tmp);
	plot(0:niter-1, tmp)
	arg = {};
	for in=1:length(Nlist)
		arg = {arg{:}, sprintf('N=%d', Nlist(in))};
	end
	legend(arg)
%	savefig fig_mri_pixel_size_cost
return
end

if im
	clf
	h(1) = plot(x1d, xtrue(:,Ndisp/2+1), 'c-');
	hold on
	llist = {'r:', 'y--', 'g-.', 'r.', 'm-'};
	arg = {'true'};
	for in=1:length(Nlist)
		N = Nlist(in);
		x1g = [-N/2:N/2-1]'/N * fov;
		tmp = xpcg{in};
		h(in+1) = plot(x1g, abs(tmp(:,N/2+1)), llist{in});
		arg = {arg{:}, sprintf('N=%d', N)};
	end
	hold off
	axis([-fov/2 fov/2 0 2.1])
	xlabel 'horizontal position [mm]'
	ylabel '|f(x,0)|'
	set(0, 'DefaultTextFontSize', 12)
%	set(0, 'DefaultAxesFontSize', 14)
	legend(h, arg{:})
%	savefig fig_mri_pixel_size_profile
end

printm(['Except for N=32, which is very under-sampled, the profiles and ', ...
'images of the other N are nearly indistinguishable.  So using "too many" ', ...
'"too small" pixels simply wastes compute time rather than degrades image', ...
' quality, thanks to regularization.'])
