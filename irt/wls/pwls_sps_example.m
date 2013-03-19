% pwls_sps_example.m
%
% Complete example m-file illustrating algorithms for PWLS cost function
% (nonquadratically penalized weighted least squares).
% Algorithms: SPS (separable paraboloidal surrogates),
% SPS-OS (order subsets verions) and CG (conjugate gradents)
%
% Copyright 2001-8-16, Jeff Fessler, The University of Michigan

%
% generate data
%
if ~isvar('yi'), disp 'setup pwls_sps_example'
	f.count = 2e5;
	em_wls_test_setup
	W = diag_sp(wi(:));
prompt
end


%
% regularizer
%
if ~isvar('Rn'), disp 'make R'
	f.l2b_n = 12;		% log_2(beta)
	f.l2b_q = 9;

	% kappa terms from fessler:96:srp paper
	kappa = sqrt( (G' * wi(:)) ./ (G' * ones(size(wi(:)))) );
	kappa = ig.embed(kappa);
	im(8, kappa, 'kappa'), cbar

	Rq = Robject(kappa, 'potential', 'quad', 'beta', 2^f.l2b_q);

	%
	% Play with regularization parameter (beta) to look at PSF.
	% I manually adjusted f.l2b above until I got fwhm=1.2 pixels or so.
	%
	if 1
		psf = qpwls_psf(G, Rq, 1, ig.mask, W);
		printm('fwhm = %g', fwhm2(psf))
		im(7, psf, 'PSF'), cbar
		clear psf
	prompt
	end

	f.delta = 0.2;		% this depends on your units!
	f.pot = 'huber';	% nonquadratic, edge-preserving penalty function
	Rn = Robject(kappa, 'type_denom', 'matlab', ...
		'potential', f.pot, 'beta', 2^f.l2b_n, 'delta', f.delta);
end


%
% Unconstrained conjugate gradient with quadratic penalty
%
if ~isvar('xpcg'), disp 'PCG'
	f.niter = 16;
	%xinit = ones(size(ig.mask));	% uniform initial image
	xinit = max(xfbp,0);		% FBP initial image

	xpcg = qpwls_pcg(xinit(ig.mask), G, W, yi(:), 0, ...
		Rq.C, 'circ0', f.niter, ig.mask);
	xpcg = ig.embed(xpcg);
	im clf, im(xpcg, 'QPWLS-CG'), cbar
prompt
end


%
% SPS iterations
%
if ~isvar('xsps'), disp 'do SPS'
	f.sps_niter = 25;
	xsps = pwls_sps_os(xinit(ig.mask), yi, wi, G, Rn, ...
		f.sps_niter, inf, [], [], 1, 0);
	xsps = ig.embed(xsps);
	im clf, im(xsps, 'PWLS-SPS'), cbar
prompt
end


%
% SPS-OS iterations
% fix: replace with incremental version
%
if ~isvar('Gb'), disp 'make Gblock'
	f.nsubset = 8;
	Gb = Gblock(G, f.nsubset, 0);
	clear tmp
end
if ~isvar('xsos'), disp 'SPS-OS'
	xsos = pwls_sps_os(xinit(ig.mask), yi, wi, Gb, Rn, ...
		f.sps_niter, inf, [], [], 1, 0);
	xsos = ig.embed(xsos);
	im clf, im(xsos, 'PWLS-SPS-OS'), cbar
prompt
end


%
% look at final results
%
if 1
	xpcgn = max(xpcg(:,:,end),0);	% final nonnegativity for CG
	xfbpn = max(xfbp,0);
	im clf
	clim = [0 9];
	elim = [-3 3];
	im(251, xtrue, 'xtrue', clim), cbar horiz
	im(252, xfbpn, 'FBP', clim), cbar horiz
	im(253, xpcgn, 'QPWLS-PCG (+)', clim), cbar horiz
	im(254, xsps(:,:,end), 'PWLS-SPS', clim), cbar horiz
	im(255, xsos(:,:,end), 'PWLS-SPS-OS', clim), cbar horiz
	efbp = xfbpn - xtrue;
	epcg = xpcgn - xtrue;
	esps = xsps(:,:,end) - xtrue;
	esos = xsos(:,:,end) - xtrue;
	nfbp = norm(efbp(:)) / norm(xtrue(:)) * 100;
	npcg = norm(epcg(:)) / norm(xtrue(:)) * 100;
	nsps = norm(esps(:)) / norm(xtrue(:)) * 100;
	nsos = norm(esos(:)) / norm(xtrue(:)) * 100;
	im(257, efbp, sprintf('NRMSE %4.1f%%', nfbp), elim), cbar horiz
	im(258, epcg, sprintf('NRMSE %4.1f%%', npcg), elim), cbar horiz
	im(259, esps, sprintf('NRMSE %4.1f%%', nsps), elim), cbar horiz
	if im, subplot(2,5,10), end
	im(esos, sprintf('NRMSE %4.1f%%', nsos), elim), cbar horiz
end
