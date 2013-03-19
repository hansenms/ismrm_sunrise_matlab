% l1_regress_example.m
% Example of l_1 regression for robust estimation
%
% Copyright 2005-4-24, Jeff Fessler, The University of Michigan

ti = [0:10]'; % sample locations
rng(0)
yi = 3 + 2 * ti;
yi = yi + 0.2 * randn(size(yi));
yi(9) = yi(1); % outlier

if ~isvar('x1') | 1, printm 'l1 regression'
	% must choose 'delta' small enough to reject outliers, but not
	% so small that convergence is too slow.  play around with it...
	f.niter = 20;
	[x1 x2] = l1_regress_fun(ti, yi, 'niter', f.niter, 'delta', 0.01, ...
		'linear', false); % use affine

	if im
		clf
		subplot(211)
		plot(x1'), xlabel 'iteration'
		axisy(0,4)
		ylabel 'estimate'
	end
	x1 = x1(:,end);
end

tt = linspace(min(ti), max(ti), 101);
if im
	subplot(212)
	plot(ti, yi, 'o', ...
		tt, x1(1) + x1(2) * tt, '-', ...
		tt, x2(1) + x2(2) * tt, '--')
	legend('data', 'l_1 regression', 'l_2 regression', 2)
end
