% mri_example_3d.m
% Example illustrating 3D regularized iterative reconstruction for MRI
% from nonuniform k-space samples.
% (This example does not include field inhomogeneity or relaxation.)
% Copyright 2004-12-16, Jeff Fessler, The University of Michigan

% this is an "honest" simulation where the Fourier data is calculated
% analytically from a continuous-space object model, even though the
% reconstructions are all discrete-space.
if ~isvar('xtrue'), printm 'setup object'

	ig = image_geom('nx', 64, 'nz', 32, ...
		'offsets', 'dsp', ... % [-n/2:n/2-1] for mri
		'fov', [20 20 6]); % 20 cm transaxial FOV

	xs = mri_objects('fov', ig.fovs, 'test4'); % object model
	xtrue = xs.image(ig.xg, ig.yg, ig.zg); % samples of continuous-space

	im clf, im pl 2 3
	clim = [0 1] * max(xtrue(:));
	im(1, xtrue, 'x true', clim), cbar
prompt
end


%
% k-space trajectory, k-space data
%
if ~isvar('yi'), printm 'trajectory'
	f.traj = 'cartesian';
	f.traj = 'spiral3';
	f.traj = 'radial';
	f.dens = {'voronoi'};
	f.dens = {};
	cpu etic
	[kspace omega wi_traj] = mri_trajectory(f.traj, {}, ig.dim, ig.fovs, f.dens);
	cpu etoc 'trajectory setup time'
	if 1, minmax(wi_traj), end
	if 0, im(2, reshape(wi_traj, ig.dim)), cbar, return, end

	if 0 && im
		im clf; plot3(omega(:,1), omega(:,2), omega(:,3), '-.')
		title(sprintf('%s: %d', f.traj, size(omega,1)))
		axis(pi*[-1 1 -1 1 -1 1]), axis_pipi
	return
	end

	printm 'setup data, based on continuous-space object model'
	ytrue = xs.kspace(kspace(:,1), kspace(:,2), kspace(:,3));

	% add noise, if desired
	rng(0)
	yi = ytrue + 0 * (randn(size(ytrue)) + 1i * randn(size(ytrue)));

	if streq(f.traj, 'cartesian')
		im(4, reshape(abs(yi), ig.dim), '|yi|'), cbar
	end
prompt
end

if 0 % sanity check in cartesian case
	ytmp = fftshift(fftn(fftshift(xtrue))); % cartesian
	ytmp = ytmp * abs(ig.dx * ig.dy * ig.dz);
	im(5, abs(ytmp), 'FFT(xtrue)'), cbar

	y3 = reshape(yi, ig.dim);
	im(6, abs(y3-ytmp), 'fft err'), cbar
	max_percent_diff(y3, ytmp)

	xtmp = fftshift(ifftn(fftshift(y3)));
	xtmp = xtmp / abs(ig.dx * ig.dy * ig.dz);
	im(3, xtmp), cbar
	max_percent_diff(xtrue, xtmp)

	im subplot 6
	plot(	ig.z, squeeze(xtrue(end/2+1,end/2+1,:)), '-o', ...
		ig.z, squeeze(xtmp(end/2+1,end/2+1,:)), '-x')
return
end

% create Gmri or Gnufft class object
if ~isvar('Gm'), printm 'system model'
	printm 'setup G objects'
	N = ig.dim;
	nufft_args = {N, 6*ones(size(N)), 2*N, N/2, 'table', 2^10, 'minmax:kb'};
	mask = true(N);
	clear N
%	f.basis = {'rect'}; % fix: using this has scaling issues
	f.basis = {'dirac'};
	Gm = Gmri(kspace, mask, ...
		'fov', ig.fovs, 'basis', f.basis, 'nufft', nufft_args);
end


if ~isvar('xcp'), printm 'conj. phase reconstruction'
	minmax(Gm.arg.basis.transform)
	% trick! adjust wi's to undo the basis effect
	wi_basis = wi_traj ./ Gm.arg.basis.transform;
%	wi_basis = ig.dx*ig.dy*ig.dz * wi_basis;
	minmax(wi_basis)

	xcp = Gm' * (wi_basis .* yi);
	xcp = embed(xcp, mask);
	im(2, xcp, 'Conj. Phase Recon'), cbar
	minmax(xtrue)
	minmax(xcp)

	if im
	im subplot 3
	plot(	ig.x, squeeze(xtrue(:,end/2+1,end/2+1)), '-o', ...
		ig.x, squeeze(xcp(:,end/2+1,end/2+1)), '-x')
	end
prompt
end

if ~isvar('R'), printm 'regularizer'
	beta = 2^-7 * size(omega,1);	% good for quadratic
	R = Reg1(ig.mask, 'beta', beta);
	if 0 % check resolution: [1.1 1.1 1]
		qpwls_psf(Gm, R, 1, ig.mask, 1, 'fwhmtype', 'profile');
	return
	end
end

if ~isvar('xpcg'), printm 'PCG with quadratic penalty'
	niter = 10;
	ytmp = yi / abs(ig.dx*ig.dy*ig.dz); % trick: analytical data vs DSFT
	xpcg = qpwls_pcg1(1*xcp(ig.mask), Gm, 1, ytmp(:), R.C, 'niter', niter);
	xpcg = ig.embed(xpcg);
	im(5, xpcg, '|x| pcg quad', clim), cbar

	if im
	im subplot 6
	plot(	ig.x, squeeze(xtrue(:,end/2+1,end/2+1)), '-o', ...
		ig.x, squeeze(xpcg(:,end/2+1,end/2+1)), '-+')
	end
end

nrms(xcp, xtrue)
nrms(xpcg, xtrue)
