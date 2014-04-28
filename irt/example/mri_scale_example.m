% mri_scale_example.m
% Examine the effect of quadratic regularization on the "scale" of
% iteratively reconstructed MRI images.  (There is no effect.)
% This was a test to resolve an internal debate...
% (This example does not include field inhomogeneity or relaxation.)
% jf

%
% create Gnufft class object
%
if ~isvar('G'),	printm 'setup Gnufft object'
	N = [32 30];
	J = [6 6];
	K = 2*N;
	% spiral trajectory with reasonable sampling along k-space axes.
	omega = linspace(0, (max(N)-1)*2*pi/2, max(N)^2+1)';
	omega = pi*[cos(omega) sin(omega)].*omega(:,[1 1])/max(omega);

	if im
		im pl 2 3
		clf, im subplot 1, plot(omega(:,1), omega(:,2), '.')
		axis([-1 1 -1 1]*pi), axis square
		title(sprintf('%d k-space samples', size(omega,1)))
	end
	G = Gnufft({omega, N, J, K});
end

%
% test data
%
if 0 || ~isvar('x'), printm 'setup object'
	x = zeros(N);
	x(round(N(1)/4+1:3*N(1)/4), round(N(2)/4+1:3*N(2)/4)) = 1;
	clim = [0 2];
	im(2, x, 'x true', clim), cbar
end

if 0 || ~isvar('yd'), printm 'setup data (slow!)'
	% generate "true" data using exact DTFT (to be fair)
	yd = dtft2(x, omega);
end


% lazy attempt at gridding
if ~isvar('xhatg'), printm 'crude gridding reconstruction'
	[xhatg, yhatg, xg] = mri_grid_linear(omega/(2*pi), yd, N, N);
	xhatg = fftshift(xhatg);

	disp(imax(yhatg, 2))
	im(3, abs(yhatg), '|y_{grid}|'), cbar
	im(4, abs(xhatg), '|x| "gridding"', clim), cbar
end

%
if 1 || ~isvar('xpcg'), printm 'PCG with quadratic penalty'
	mask = true(N);
	niter = 40;
	% playing around with this beta does not change the summations
	% at the bottom very much, except if beta is way too large
	beta = 2^-3 * size(omega,1);
	R = Robject(mask, 'beta', beta);
	xinit = zeros(N);
	xpcg = qpwls_pcg(xinit(:), G, 1, yd(:), 0, R.C, 1, niter, mask);
	xpcg = embed(xpcg(:,end), mask);
	im(5, abs(xpcg), '|x| pcg quad', clim), cbar
end

if im
%	im subplot 6
	clf
	pro = xpcg(:,N(2)/2+1);
	plot(	...
		xg{1}, x(:,N(2)/2+1), 'y-', ...
		xg{1}, real(pro), 'c--', ...
		xg{1}, imag(pro), 'r--', ...
		xg{1}, abs(pro), 'm--')
end

4*[mean(x(:)) mean(abs(real(xpcg(:)))) mean(abs(imag(xpcg(:)))) mean(abs(xpcg(:)))]
