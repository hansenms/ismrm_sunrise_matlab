% mri_sensemap_demo1.m
% example of regularized estimation of MR sensitivity maps
% from body/surface coils
%
% Copyright 2005-9-12, Jeff Fessler, The University of Michigan

if ~isvar('ftrue')
	f.dir = [path_find_dir('mri') '/../data/mri/'];
	f.xtrue = [f.dir 'brainweb_t1.jpg'];
	ftrue = double6(imread(f.xtrue)');
	ftrue = ftrue(2:end-1,2:end-1); % make it 256^2
	ftrue = downsample2(ftrue, 4); % now 64^2
	im(ftrue)
	[nx ny] = size(ftrue);
	if 0 % simple synthetic rho map
		ftrue = zeros(size(ftrue));
		ftrue(nx/4:3*nx/4,ny/4:3*ny/4) = 2^7;
		ftrue(5+nx/2+[-2:2],ny/2+[-10:10]) = 2^2;
%		ftrue = 128 * ones(size(ftrue));
	end
end

if ~isvar('zlj'), printm 'zlj'
	im row 2
	ncoil = 4;
	smap = mri_sensemap_sim('nx', nx, 'ny', ny, 'dx', 192/nx, ...
		'ncoil', ncoil, 'orbit', 90*ncoil, 'rcoil', 100);
	if ncoil == 4, smap = smap(:,:,[4 1 3 2]); end
%	smap(nx/2+[-1:1]*10,ny/2+[-1:1]*10,1) = 1.5;

	f.sig = 1.0; % noise in image: bodycoil and array coils
	rng(0)
	yj = ftrue + f.sig * (randn(size(ftrue)) + 1i * randn(size(ftrue)));
	zlj = smap .* repmat(ftrue, [1 1 ncoil]) + f.sig * ...
		(randn(size(smap)) + 1i * randn(size(smap)));
%	f.snr_wrong = dB(norm(ftrue(:)) / norm(yj-ftrue));
%	pr f.snr_wrong
	f.snr = dB(norm(ftrue(:)) / norm(col(yj-ftrue))); % corrected 2012-02-21
	pr f.snr

	im pl 2 2
	im(1, smap, 'True sense maps'), cbar
	im(2, abs(zlj), '|zlj|'), cbar
end

% conventional "ratio" estimate of sensitivity map
if ~isvar('shat1')
	mask0 = conv2(double(ftrue > 0), ones(5), 'same') > 0;

	for imap = 1:ncoil
		zj = zlj(:,:,imap);
		tmp = zj ./ yj; % usual ratio
		if 1 % set all uncertain map values to median of good ones
			good = abs(zj) > 0.05 * max(abs(zj(:)));
			tmp(~good) = median(tmp(good));
		end
		shat1(:,:,imap) = tmp;
	end
	clear good zj tmp

	slim = [0 1.4];
	im(3, abs(shat1), 'Ratio sense map estimates', slim), cbar
prompt
end

%
% QPWLS approach
%
if ~isvar('args')
	f.niter = 2^8;
	f.l2b = -3;
f.l2b = 2;

	args = {zlj(:,:,:), 'order', 2, 'l2b', f.l2b, 'niter', f.niter, ...
		'bodycoil', yj};
end

if 0 & ~isvar('shat2b'), printm 'try precon'
	[shat2a dummy cost2a] = mri_sensemap_denoise(args{:}, ...
		'precon', 1, 'isave', 'all');
	[shat2b dummy cost2b] = mri_sensemap_denoise(args{:}, ...
		'precon', 'diag', 'isave', 'all');
%	plot cost function to see if smaller faster - very small difference?

	clf, plot(0:f.niter, cost2a(:,1), '-o', ...
		0:f.niter, cost2b(:,1), '-x')
	axisx(0,20)
	legend('normal', 'diag')
return
end

if 0 | ~isvar('shat2'), printm 'regularized method'
	[shat2 dummy cost2c] = mri_sensemap_denoise(args{:});
%prompt
end

elim = [0 0.25];
if 0, printm 'slow way for comparing' % to test convergence vs exact inverse
	shat3 = mri_sensemap_denoise(args{:}, 'slow', true);
	im clf, im pl 3 3
	im(1, ftrue, 'f')
	im(2, zlj(:,:,1), 'z1')
	im(3, smap(:,:,1), 'smap1', slim)
	im(4,  shat1(:,:,1), 'ratio', slim)
	im(5,  shat2, 'iter', slim)
	im(6,  shat3, 'slow', slim)
	im(8,  abs(shat2-smap(:,:,1)), 'err', elim)
	im(9,  abs(shat3-smap(:,:,1)), 'err', elim)
return
end


%
% old approach (obsolete due to mri_sensemap_denoise)
%
if 0 | ~isvar('shat2'), printm 'regularized method'
	mask1 = true(nx,ny);
	A = diag_sp(yj(mask1));
%	A = diag_sp(ftrue(mask1)); 'testing with ftrue'

%	R = Robject(mask1, 'beta', 2^4);
	R = Robject(mask1, 'beta', 2^f.l2b, 'order', 2, 'distance_power', 2);

	wj = abs(yj);
	wj = median(wj(wj > 0.05 * max(wj(:))));
	W = 1 ./ wj^2; % trick: these weights make the beta "universal"
	if 1 % examine psf
		qpwls_psf(1., R, 1., mask1, 1., 'offset', [8 0]);
		psf = qpwls_psf(A, R, 1., mask1, W, 'offset', [8 0]);
		im(224, psf, 'psf')
	prompt
	end

	for imap = 1:ncoil
		init = shat1(:,:,imap);

		ztmp = zlj(:,:,imap);
		if 0 % test precon (doesn't help!?)
			f.niter = 30;
			M = embed(feval(R.handle_diag, R), mask1);
			M = 1e0 + abs(ftrue).^2 + M;
			% im(M), return
			M = diag_sp(1 ./ M(mask1));
			tmp2 = qpwls_pcg1(init(mask1), A, W, ztmp(:), R.C, ...
				'isave', 'all', 'niter', f.niter, 'precon', M);
			tmp2 = embed(tmp2, mask1);
			tmp1 = qpwls_pcg1(init(mask1), A, W, ztmp(:), R.C, ...
				'isave', 'all', 'niter', f.niter);
			tmp1 = embed(tmp1, mask1);
			cost1 = pwls_cost(tmp1, A, W, ztmp(:), R, mask1);
			cost2 = pwls_cost(tmp2, A, W, ztmp(:), R, mask1);
			plot(0:f.niter, cost1, '-o', 0:f.niter, cost2, '-+')
		return
		end

		tmp = qpwls_pcg1(init(mask1), A, W, ...
			ztmp(:), R.C, 'niter', f.niter);
		shat2(:,:,imap) = embed(tmp, mask1);

		if 0 % cf unif init
			init = ones(nx,ny);
%			ej = zeros(nx,ny); ej(nx/2+1,ny/2+1) = 1;
%			ztmp = ztmp + reshape(A * ej(mask1), size(ztmp));
			tmp = qpwls_pcg1(init(mask1), A, W, ...
				ztmp(:), R.C, 'niter', f.niter);
			shat2u(:,:,imap) = embed(tmp, mask1);
		end
	end
	if 0
		max_percent_diff(shat2, shat2u)
		mask3 = repmat(mask0, [1 1 ncoil]);
		max_percent_diff(shat2 .* mask3, shat2u .* mask3)
		im(221, shat2u, 'shat2u', slim), cbar
	end
	im(224, shat2, 'Regularized sense map estimates', slim), cbar
prompt
end
%psf = reale(shat2u-shat2, 'warn');
%psf = psf(:,:,1);
%im(223, psf, 'psf'), cbar
%fwhm2(psf)

if 1 % figures for viewing
	im clf, im pl 4 1
	im row 1
	im(1, abs(zlj), 'Array coil images', [0 120]), cbar
	im(2, smap, 'True sensivity maps', slim), cbar
	im(3, abs(shat1), 'Ratio sensitivity maps', slim), cbar
	im(4, abs(shat2), 'Regularized sensitivity maps', slim), cbar
return
end

if 1 % figures for publication.  todo: look at real/imag parts
	im clf
	im(smap, 'True sensivity maps', slim), cbar
%	savefig 'fig_mr_sensemap1_map'
	im(abs(zlj), 'Array coil images')
%	savefig 'fig_mr_sensemap1_zl'
	im(abs(shat1), 'Ratio sensitivity maps', slim), cbar
%	savefig 'fig_mr_sensemap1_ratio'
	im(abs(shat2), 'Regularized sensitivity maps', slim), cbar
%	savefig 'fig_mr_sensemap1_reg'
return
end

if 1
	mask3 = 1;
	serr1 = (shat1 - smap) .* mask3;
	serr2 = (shat2 - smap) .* mask3;

	im(121, abs(serr1), 'Error: ratio estimates', elim)%, cbar h
	im(122, abs(serr2), 'Error: regularized estimates', elim)%, cbar h
%	savefig 'fig_mr_sensemap1_err'
return
end
