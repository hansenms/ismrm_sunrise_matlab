% l1_regress_example.m
% Example of l_1 regression for robust estimation
%
% Copyright 2005-4-24, Jeff Fessler, University of Michigan

%% example data
ti = [0:10]'; % sample locations
rng(0)
yi = 3 + 2 * ti;
yi = yi + 0.2 * randn(size(yi));
yi(9) = yi(1); % outlier

%% run iterative reweighted LS algorithm (Huber's method)
if ~isvar('x1') || 0, printm 'l1 regression'
	f.niter = 10;
	% must choose 'delta' small enough to reject outliers, but not
	% so small that convergence is too slow.  play around with it...
	[x1s x2] = l1_regress_fun(ti, yi, 'niter', f.niter, 'delta', 0.01, ...
		'linear', false); % use affine
	x1 = x1s(:,end);
end

%% ADMM version
if ~isvar('x1_admm') || 1, printm 'l1 regression via admm'
	A = [ones(size(ti)) ti];
	x1s_admm = l1_regress_admm1(yi, A, 'niter', f.niter, ...
		'x0', x1, ...
		'x0', [0; 0], ...
		'x0', [], ...
		'rho', 2^-3, ... % trick: had to tune :(
		'isave', 'all');
	x1_admm = x1s_admm(:,end);
end


%% plot parameters vs iteration and robust regression lines
if im
	clf
	subplot(211)
	plot(0:f.niter, x1s', '-o', 0:f.niter, x1s_admm', '-+')
	xlabel 'iteration'
	axisy(0,4)
	ylabel 'estimate'

	tt = linspace(min(ti), max(ti), 101);

	subplot(212)
	plot(ti, yi, 'bo', ...
		tt, x1(1) + x1(2) * tt, 'g-', ...
		tt, x1_admm(1) + x1_admm(2) * tt, 'm-', ...
		tt, x2(1) + x2(2) * tt, 'r--')
	legend('data', 'l_1 regression', 'l_1 via admm', 'l_2 regression', ...
		'location', 'north')
end
