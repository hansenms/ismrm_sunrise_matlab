% mri_example.m
% A small example illustrating iterative reconstruction for MRI
% from a few nonuniform k-space samples using edge-preserving regularization.
% (This example does not include field inhomogeneity or relaxation.)
% Or more generally, this shows how to go from nonuniform samples
% in the frequency domain back to uniform samples in the space domain
% by an interative algorithm.
% Copyright 2003, Jeff Fessler, The University of Michigan

%
% create Gnufft class object
%
if ~isvar('Gm'), printm 'setup Gnufft object'
	ig = image_geom('nx', 32, 'ny', 28, 'dx', 1, 'offsets', 'dsp');
	ig.mask = ig.circ(1+ig.nx/2, 1+ig.ny/2) > 0;
	N = ig.dim;
	J = [6 6];

	[kspace omega wi] = mri_trajectory('spiral0', {}, N, ig.fov, {'voronoi'});
	warn 'todo: crudely fixing dcf at edge of spiral'
	wi(350:end) = wi(350);

	im pl 3 3
	if im
		clf, im subplot 1, plot(omega(:,1), omega(:,2), '.')
		axis([-1 1 -1 1]*pi), axis square
		titlef('%d k-space samples', size(omega,1))
	end
	nufft_args = {N, J, 2*N, N/2, 'table', 2^10, 'minmax:kb'};

	cpu etic
	Gn = Gnufft(ig.mask, {omega, nufft_args{:}});
	cpu etoc 'Gnufft build time'
	cpu etic
	Gm = Gmri(kspace, ig.mask, ...
		'fov', ig.fov, 'basis', {'rect'}, 'nufft', nufft_args);
	cpu etoc 'Gmri build time'
prompt, clear kspace omega nufft_args
end

if ~isvar('Tm'), printm 'Tm'
	tic
	try
		Tm = build_gram(Gm, 1);
		printm('build gram time: %g', toc)
	catch
		Tm = []; % fake empty matrix
	end
	if 0
%		x0 = ig.unitv;
		x1 = ig.embed(Gm' * (Gm * x0(ig.mask)));
		x2 = ig.embed(Tm * x0(ig.mask));
		im clf, im(stackup(x1,x2))
	return
	end
end

%
% test data
%
clim = [0 2];
if 0 || ~isvar('x'), printm 'setup object'
	x = ig.zeros;
	if 1
		x(5:25,5:25) = 1;
		x(10:20,10:15) = 2;
		x(20:22,18:22) = 0;
		x(7:11,18:22) = 0;
		x(15,20) = 2;
		angtrue = zeros(size(x));
		[tx ty] = ndgrid(ig.x/ig.nx, ig.y/ig.ny);
%		ndgrid([-N(1)/2:N(1)/2-1]/N(1), [-N(2)/2:N(2)/2-1]/N(2));
%		angtrue(15:25,15:25) = 0.3;
		angtrue = 0.25 * (1 + cos(min(sqrt(tx.^2+ty.^2)*2*pi,pi)));
		x = x .* exp(1i * angtrue);
	else
%		x = ones(N);
		x(N(1)/2+1,10) = 1;
	end, clear tx ty
	im(2, abs(x), '|x| true', clim), cbar

	pmask = ones(size(x));	% phase-map mask
	plim = [-0.5 0.5];
	im(3, angtrue .* (abs(x) > 0), '\angle x true', plim), cbar
prompt
end

if ~isvar('yi'), printm 'setup data'
	if 1	% generate "true" data using exact DTFT (to be fair)
		yd = dtft2(x, Gn.arg.st.om, Gn.arg.st.n_shift);
	else	% "cheat" by using NUFFT to generate data
		yd = Gn * x;
	end
	yi = yd; % noiseless for now
end

if 0
	[oo1 oo2] = ndgrid(	2*pi*([0:N(1)-1]/N(1) - 0.5), ...
				2*pi*([0:N(2)-1]/N(2) - 0.5));
	yd_g = griddata(omega(:,1), omega(:,2), yi, oo1, oo2, 'cubic');
	yd_g(isnan(yd_g)) = 0;

	disp(imax(yd_g, 2))
%	im(3, abs(yd_g), '|y_d|'), cbar
end

% horribly lazy attempt at gridding
if 0 && ~isvar('xg'), printm 'crude gridding reconstruction'

	xg = ifft2(fftshift(yd_g));
%	xg = fftshift(xg);
%	xg = xg(:, [N(2)/2+[1:N(2)/2] 1:N(2)/2]);
	im(4, abs(xg), '|x| "gridding"', clim), cbar
end

if ~isvar('xcp'), printm 'conj phase recon'
	xcp = Gn' * (wi .* yi);
	xcp = ig.embed(xcp);
	im(4, abs(xcp), '|x| "conj phase"', clim), cbar
	im(7, angle(xcp), '\angle x "conj phase"', plim), cbar
prompt
end

if ~isvar('Rq'), printm 'Rq quadratic'
	beta = 2^-7 * size(yi,1);	% good for quadratic
	Rq = Robject(ig.mask, 'beta', beta);
	psf = qpwls_psf(Gm, Rq.C, 1, ig.mask, 1);
	im(8, psf, 'psf')
prompt, clear psf
end


%
if ~isvar('xpcg'), printm 'PCG with quadratic penalty'
	niter = 40;
%	xinit = zeros(N);
	xinit = xcp(ig.mask(:));
	tic
	xpcg = qpwls_pcg(xinit, Gm, 1, yi(:), 0, Rq.C, 1, niter, ig.mask);
	tim.v1 = toc;
	xpcg = ig.embed(xpcg(:,end));
	[magq angq] = mag_angle_real(xpcg);
	im(5, magq, '|x| pcg quad', clim), cbar
	im(8, angq.*pmask, '\angle x pcg quad', plim), cbar
prompt
end


% compare speed of NUFFT and Toeplitz approach
if ~isvar('xpcg2') && ~isempty(Tm), printm 'Toeplitz approach'
	bb = Gm' * yi(:);
	xpcg2 = qpwls_pcg2(xinit, Tm, bb, Rq.C, 'niter', niter);
	tim.v2 = toc;
	printm('Times: nufft=%g toep=%g', tim.v1, tim.v2)
	xpcg2 = ig.embed(xpcg2(:,end));
	printm('nufft vs toeplitz mpd=%g%%', max_percent_diff(xpcg, xpcg2))
	if 0
		im clf, im(131, xpcg), im(132, xpcg2), im(133, xpcg-xpcg2)
	end
	[mag2 ang2] = mag_angle_real(xpcg);
	im(6, mag2, '|x| pcg quad', clim), cbar
	im(9, ang2.*pmask, '\angle x pcg quad', plim), cbar
prompt, clear mag2 ang2 bb
end

if ~isvar('xh'), printm 'PCG with edge-preserving penalty'
	Rn = Robject(ig.mask, 'type_denom', 'matlab', ...
		'potential', 'hyper3', 'beta', 2^2*beta, 'delta', 0.3);
	xh = pwls_pcg1(xpcg(ig.mask), Gm, 1, yi(:), Rn, 'niter', 2*niter);
	xh = ig.embed(xh);
	[magn angn] = mag_angle_real(xh);
	im(6, magn, '|x| pcg edge', clim), cbar
	im(9, angn.*pmask, '\angle x pcg edge', plim), cbar

printm(['So with only %d k-space samples the iterative method has' ...
 ' reconstructed %d of %d image pixels.  Edge-preserving regularization' ...
 ' worked particularly well, albeit without noise.'], ...
		size(yi,1), ig.np, prod(N))
	if exist('mri_recon2.m', 'file'), prompt, end
prompt
end, clear magn angn

if ~exist('mri_recon2.m', 'file'), return, end

% work in progress
if ~isvar('x2'), printm 'trying new'
	beta_phase = beta * 2^4;
	beta_mag = beta / 2^2;
	arg = {'edge_type', 'leak', 'type_denom', 'matlab'};
	Rphase = Robject(ig.mask, arg{:}, 'beta', beta_phase);
	Rmag = Robject(ig.mask, arg{:}, 'beta', beta_mag);
	if 0	% correct magnitude, noisy phase estimate
		finit = (abs(x)+eps) .* exp(1i * angle(xpcg));
%		finit = x;
	elseif 0 % pcg magnitude, with correct phase
		finit = (abs(xpcg)+eps) .* exp(1i * angle(x));
	else	% pcg, with cleaned up phase
		finit = xpcg;
%		finit = abs(xpcg);	% no phase info.!?
		ttt = mri_phase_denoise(xpcg, -5, 50, 0);
		finit = (abs(xpcg)+eps) .* exp(1i * ttt);
	end
	if ~isvar('x2')
		x2 = mri_recon2(finit(ig.mask), Gm, yd, Rmag, Rphase);
		x2 = ig.embed(x2);
	end
	[mag2 ang2] = mag_angle_real(x2);
	im(4, mag2, '|x| new', clim), cbar
	im(7, ang2.*pmask, '\angle x new', plim), cbar
	printm 'ang2 range:', disp(minmax(ang2.*pmask)')
prompt
end

if 1, % figure for ISBI
	im clf, im pl 3 2
	plim = [-0.1 0.6];
	im(1, abs(x), '|x| true', clim), cbar
	im(2, angtrue, '\angle x true', plim), cbar
	t = sprintf('|x| old, NRMS=%2.f%%', 100*nrms(magq(:), abs(x(:))));
	im(3, magq, t, clim), cbar
	im(4, angq.*pmask, '\angle x old', plim), cbar
	t = sprintf('|x| new, NRMS=%2.f%%', 100*nrms(mag2(:), abs(x(:))));
	im(5, mag2, t, clim), cbar
	im(6, ang2.*pmask, '\angle x new', plim), cbar
	printm 'mag nrms:', disp(nrms([magq(:) mag2(:)], abs(x(:)))')
end
