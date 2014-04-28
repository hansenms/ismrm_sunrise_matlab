 function [smap, sinit, cost] = mri_sensemap_denoise(ykj, varargin)
%function [smap, sinit, cost] = mri_sensemap_denoise(ykj, [options])
%|
%| Regularized estimation of (smooth) sensitivity maps for parallel MRI.
%| Signal model: y_kj = s_kj f_j + noise_kj
%|	s_kj	sensitivity map (relative to "bodycoil" reference image)
%|	f_j	unknown underlying object
%|	k: coil index
%|	j: voxel index
%| Signal model for "body coil" image: z_j = f_j + noise_j.
%| If no body coil image provided, then use sum-of-squares of y_kj,
%| multiplied by phase of the first coil image.
%|
%| This method avoids the problematic "ratio" of conventional methods.
%| It also smoothly interpolates over regions with signal voids.
%|
%| in
%|	ykj	[(N) ncoil]	noisy complex images (2D or 3D) for each coil
%|
%| options
%|	bodycoil [(N)]		reference image (optional)
%|				(default: sum-of-squares with phase of 1st coil)
%|	l2b			log_2(beta), regularization parameter (def: -5)
%|	psf 1|0			report fwhm of psf? (def: true)
%|	order			regularization order (default: 2)
%|	niter			# of iterations (default: 150)
%|	thresh			fraction of magnitude maximum used for median
%|				initial value in "background" (default: 0.05)
%|	init	[(N) ncoil]	initial sense maps for iterations
%|				(default: [] = standard ratio)
%|	isave			which iterations to save. (default: last)
%|	'chol'	0|1		use non-iterative cholesky decomposition method
%|				that works well for small images. default=0
%|
%| out
%|	smap	[(N) ncoil]	denoised sensitivity maps
%|	sinit	""		initial maps (if not provided)
%|	cost	[niter ncoil]	(optional) cost function vs iteration
%|
%| This method was explored in the 2008 ISMRM abstract by Kim etal, p. 1267,
%| "Smoothing effect of sensitivity map on fMRI data using a novel
%| regularized self-calibrated estimation method"	kim:08:seo
%| See also the Mar. 2013 IEEE T-MI paper by M J Allison et al
%| doi 10.1109/TMI.2012.2229711
%|
%| Copyright 2005, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(ykj, 'test')
	smap = mri_sensemap_denoise_test;
	if ~nargout, clear smap, end
return
end

% defaults
arg.l2b = -5;
arg.order = 2;
arg.niter = 150;
arg.isave = [];
arg.bodycoil = [];
arg.thresh = 0.05;
arg.init = [];
arg.psf = true;
arg.precon = 1; % default: no preconditioner
arg.chol = false; % set to true to use "exact" matrix inverse (slow if big)

arg = vararg_pair(arg, varargin, 'subs', {'slow', 'chol'});
if isempty(arg.isave), arg.isave = arg.niter; end

tmp = size(ykj);
NN = tmp(1:end-1);
ncoil = tmp(end);
coildim = numel(tmp);

% default reference image
if isempty(arg.bodycoil)
	arg.bodycoil = sqrt(sum(abs(ykj).^2, coildim)); % sum-of-squares
	tmp = angle(stackpick(ykj,1)); % crude estimate of phase of image f_j
	arg.bodycoil = arg.bodycoil .* exp(1i * tmp);
end

% initial maps
if isempty(arg.init)
	arg.init = zeros([prod(NN) ncoil]);
	good = abs(arg.bodycoil) > arg.thresh * max(abs(arg.bodycoil(:))); % 2010-06-24
	for ic = 1:ncoil
		zj = stackpick(ykj,ic);
		tmp = zj ./ arg.bodycoil; % usual ratio
		if 1 % set all uncertain map values to median of good ones
%			good = abs(zj) > arg.thresh * max(abs(zj(:))); % prior to 2010-06-24
			tmp(~good) = median(abs(tmp(good)));
		end
		arg.init(:,ic) = tmp(:);
	end
	sinit = reshape(arg.init, [NN ncoil]); % return to caller if wanted
end


% regularizer. todo: could penalize z different from xy for non-isotropic voxels
mask = true(NN); % estimate / extrapolate to *all* pixels
if 0 % experiment with excluding image boundary
	mask(:,[1 end]) = false;
	mask([1 end],:) = false;
end
args = {mask, 'beta', 2^arg.l2b, 'order', arg.order, 'distance_power', 2};
if arg.chol
	R = Reg1(args{:}, 'type_penal', 'mat', 'type_diff', 'spmat');
else
	R = Reg1(args{:});
	%R = Robject(mask, 'beta', 2^arg.l2b, 'order', arg.order, 'distance_power', 2);
end

% report expected blur (at image center)
if arg.psf
pr arg.l2b
	qpwls_psf(1., R, 1., mask);
end

% trick: normalize data by median of non-background value in bodycoil image
% so that the effective regularizer beta is "universal" (scale invariant)
tmp = abs(arg.bodycoil);
tmp = median(tmp(tmp > arg.thresh * max(tmp(:))));
%pr tmp
arg.bodycoil = arg.bodycoil / tmp;
ykj = ykj / tmp;

if arg.chol, printm 'make sparse hessian'
	if arg.niter ~= 1
		fail '"chol" version expects niter=1'
	end
	A = Gdiag(arg.bodycoil(mask(:)), 'mask', mask);
	A = sparse(A);
	C = R.C;
	H = A' * A + C' * C;
	if 0 % test equivalence of gradients
		Ro = Reg1(args{:});
		t1 = C' * (C * sinit(:));
		t2 = Ro.cgrad(Ro, sinit(:));
		equivs(t1, t2)
	end
else
	A = Gdiag(arg.bodycoil(mask(:)), 'mask', mask);
end

if streq(arg.isave, 'all')
	smap = zeros(prod(NN), ncoil, arg.niter+1);
else
	smap = zeros(prod(NN), ncoil, numel(arg.isave));
end

if streq(arg.precon, 'diag')
	if R.order == 2
		arg.precon = 6 * sum(R.beta); % (-1)^2 + 2^2 + (-1)^2
	else
		arg.precon = 2 * sum(R.beta); % (-1)^2 + (1)^2
	end
	arg.precon = abs(arg.bodycoil).^2 + arg.precon;
%	im(arg.precon), prompt
	arg.precon = 1 ./ arg.precon;
	arg.precon = Gdiag(arg.precon);
end

for ic = 1:ncoil
	ytmp = stackpick(ykj,ic);
	if arg.chol % numerical inverse for small problems
		ytmp = double(ytmp);
		tmp = A' * ytmp(:);
		tmp = H \ tmp;
	else % run qpwls algorithm for regularized fitting
		init = stackpick(arg.init,ic);
		tmp = qpwls_pcg1(init(mask), A, 1, ...
			ytmp(mask), R.C, ...
			'precon', arg.precon, ...
			'niter', arg.niter, 'isave', arg.isave);
		if nargout > 2 % cost
			cost(:,ic) = pwls_cost(tmp, A, 1, ytmp(mask(:)), R);
		end
	end
	smap(:,ic,:) = col(embed(tmp, mask));
end

if streq(arg.isave, 'all')
	smap = reshape(smap, [NN ncoil arg.niter+1]);
else
	smap = reshape(smap, [NN ncoil length(arg.isave)]);
end


% mri_sensemap_denoise_test()
% built-in test/example
function smap = mri_sensemap_denoise_test
mri_sensemap_demo1
