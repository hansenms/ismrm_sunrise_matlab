% mri_sample_density_example.m
% NOT DONE!
% Evaluation various sample density compensators
% Copyright 2003-7-29	Jeff Fessler	The University of Michigan

%
% create k-space samples, in units of 1/mm
%
if ~isvar('kspace'), printm 'setup kspace'
	fov = 256;	% [mm] typical brain FOV
	N0 = 64;	% nominal image size
	kmax = N0/2*(1/fov);	% for display axes

%	ktype = 'cartesian';
%	ktype = 'random';
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

	if im
		im clf, im pl 2 3
		im subplot 1
		plot(kspace(:,1), kspace(:,2), '.')
		axis(1.1*[-1 1 -1 1]*kmax), axis square
		xlabel 'k_1 [mm^{-1}]', ylabel 'k_2 [mm^{-1}]'
		title(sprintf('%d %s k-space samples', size(kspace,1), ktype))
	end
prompt
end


%
% true object and analytical k-space data
%
if 0 || ~isvar('xtrue'), printm 'setup object'
	% display images with many pixels...
	x1d = [-N0/2:N0/2-1] / N0 * fov;
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
	[xg yd_g x1g k1g] = mri_grid_linear(kspace, yd, N0, fov);

	disp(imax(yd_g, 2))
	im(3, k1g{1}, k1g{2}, abs(yd_g), '|y_d|'), cbar
	im(4, x1g{1}, x1g{2}, abs(xg), '|x| "gridding"', clim), cbar

	mask = true(N0);
prompt
end, clear yd_g k1g

%
% create Gnufft class object
%
if 0 || ~isvar('G'), printm 'setup Gnufft object'
	omega = 2*pi*kspace*fov/N0;
%	disp(minmax(omega(:,1)))
	G = Gnufft({omega, [N0 N0], [6 6], 2*[N0 N0], [N0/2 N0/2]});
%	G = Gnufft({omega, [N0 N0], [8 8], 4*[N0 N0], [N0/2 N0/2], 'kaiser'});
prompt
end

% gridding
if ~isvar('xgr.unif') && 0
	xgr.unif = embed(G' * yd, mask);
	im(5, x1g, x1g, abs(xgr.unif), 'no comp')
end

if ~isvar('wt.pipe')
	w = ones(size(yd));
	P = G.arg.st.p;
	for ii=1:20
		tmp = P * (P' * w);
		w = w ./ real(tmp);
	end
	minmax(tmp)
	wt.pipe = w;
	if im
		plot(wt.pipe)
	end
end

scale = G.arg.st.sn(end/2,end/2)^(-2) / fov^2 / prod(G.arg.st.Kd) * N0^2;

% gridding
if ~isvar('xgr.pipe')
	w = wt.pipe * scale;
	xgr.pipe = embed(G' * (w .* yd), mask);
	im(6, x1g{1}, x1g{2}, real(xgr.pipe), 'pipe', clim), cbar
end
sum(xgr.pipe) / sum(xtrue)
printf('nrms pipe %g%%', 100*nrms(xgr.pipe(:), xtrue(:)))

return %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% reconstruct by PCG
%
if 0 || ~isvar('xpcg'), printm 'PCG with quadratic penalty'

	niter = 20;
	beta = 2^-10 * size(omega,1);

	R = Robject(mask, 'beta', beta);
	ytmp = yd(:) * (1/fov)^2 * N0^2;	% scaling!
	xiter = qpwls_pcg(xg(:), G, 1, ytmp, 0, R.C, 1, niter);
	xpcg = embed(xiter(:,end), mask);
	im(3, x1g, x1g, abs(xpcg), '|x| pcg quad'), cbar
prompt
end

% ir_savefig fig_mri_sample_density
