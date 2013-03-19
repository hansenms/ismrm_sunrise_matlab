  function [smap, sinit, cost] = mri_sensemap_denoise(ykj, varargin)
%|function [smap, sinit, cost] = mri_sensemap_denoise(ykj, [options])
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
%|	ykj	[nx ny ncoil]	noisy complex images for each coil
%| todo: handle 3D!
%|	ykj	[(*N) ncoil]	noisy complex images for each coil
%|
%| options
%|	bodycoil [nx ny]	reference image
%|				(default: sum-of-squares with phase of 1st coil)
%|	l2b			log_2(beta), regularization parameter (def: -5)
%|	psf 1|0			report fwhm of psf? (def: true)
%|	order			regularization order (default: 2)
%|	niter			# of iterations (default: 150)
%|	thresh			fraction of magnitude maximum used for median
%|				initial value in "background" (default: 0.05)
%|	init	[nx ny ncoil]	initial sense maps for iterations
%|				(default: standard ratio)
%|	isave			which iterations to save.  (default: last)
%|
%| out
%|	smap	[nx ny ncoil]	denoised sensitivity maps
%|	sinit	""		initial maps (if not provided)
%|	cost	[niter ncoil]	(optional) cost function vs iteration
%|
%| This method was explored in the 2008 ISMRM abstract by Kim etal, p. 1267,
%| "Smoothing effect of sensitivity map on fMRI data using a novel
%| regularized self-calibrated estimation method"	kim:08:seo
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
arg.slow = false; % set to true to use slow but "exact" matrix inverse

arg = vararg_pair(arg, varargin);
if isempty(arg.isave), arg.isave = arg.niter; end

[nx ny ncoil] = size(ykj);

% default reference image
if isempty(arg.bodycoil)
	arg.bodycoil = sqrt(sum(abs(ykj).^2, 3)); % sum-of-squares
	tmp = angle(ykj(:,:,1)); % crude estimate of phase of image f_j
	arg.bodycoil = arg.bodycoil .* exp(1i * angle(ykj(:,:,1)));
end

% initial maps
if isempty(arg.init)
	arg.init = zeros(nx,ny,ncoil);;
	good = abs(arg.bodycoil) > arg.thresh * max(abs(arg.bodycoil(:))); % 2010-06-24
	for ic = 1:ncoil
		zj = ykj(:,:,ic);
		tmp = zj ./ arg.bodycoil; % usual ratio
		if 1 % set all uncertain map values to median of good ones
%			good = abs(zj) > arg.thresh * max(abs(zj(:))); % prior to 2010-06-24
			tmp(~good) = median(tmp(good));
		end
		arg.init(:,:,ic) = tmp;
	end
	sinit = arg.init; % return to caller if wanted
end


% regularizer
mask = true(nx,ny); % estimate / extrapolate to *all* pixels
if 0 % experiment with excluding image boundary
	mask(:,[1 end]) = false;
	mask([1 end],:) = false;
end
args = {mask, 'beta', 2^arg.l2b, 'order', arg.order, 'distance_power', 2};
if arg.slow
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

if arg.slow, printm 'make sparse hessian'
	A = spdiag(arg.bodycoil(mask(:)), 'nowarn');
	C = R.C;
	H = A' * A + C' * C;
	if 0 % test equivalence of gradients
		Ro = Reg1(args{:});
		t1 = C' * (C *  sinit(:));
		t2 = Ro.cgrad(Ro, sinit(:));
		equivs(t1, t2)
	end
else
	A = diag_sp(arg.bodycoil(mask(:)));
end

if streq(arg.isave, 'all')
	smap = zeros(nx,ny,ncoil,arg.niter+1);
else
	smap = zeros(nx,ny,ncoil,length(arg.isave));
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
	arg.precon = diag_sp(arg.precon);
end

for ic = 1:ncoil
	ytmp = ykj(:,:,ic);
	if arg.slow % numerical inverse for testing
		ytmp = double(ytmp);
		tmp = A' * ytmp(:);
		tmp = H \ tmp;
	else % run qpwls algorithm for regularized fitting
		init = arg.init(:,:,ic);
%		init = arg.init(:,...,:,ic);
%		init = arg.init(:,ic);
		tmp = qpwls_pcg1(init(mask), A, 1, ...
			ytmp(mask(:)), R.C, ...
			'precon', arg.precon, ...
			'niter', arg.niter, 'isave', arg.isave);
		if nargout > 2 % cost
			cost(:,ic) = pwls_cost(tmp, A, 1, ytmp(mask(:)), R);
		end
	end
	smap(:,:,ic,:) = embed(tmp, mask);
end

%
% built-in test/example
%
function smap = mri_sensemap_denoise_test %(type, varargin)
mri_sensemap_demo1
