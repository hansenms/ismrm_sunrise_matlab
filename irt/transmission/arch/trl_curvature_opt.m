 function ni = trl_curvature_opt(yi, bi, ri, li, yb)
%function ni = trl_curvature_opt(yi, bi, ri, li, yb)
%
%	Compute optimal surrogate parabola curvatures
%	for Poisson transmission model based on Erdogan's formula.
%	fix: align with C version of trpl2,3
%	The minimum returned curvature will be zero.
%	It is the user's responsibility to impose additional bounds
%	if desired for certain algorithms.
%	Copyright 2002-1-28	Jeff Fessler	The University of Michigan

	h = inline('y.*log(b.*exp(-l)+r)-(b.*exp(-l)+r)', 'y','b','r','l');
	dh = inline('(1 - y ./ (b.*exp(-l)+r)) .* b.*exp(-l)', 'y','b','r','l');

%
%	the default is to demonstrate an example surrogate parabola!
%
if nargin == 0
	help(mfilename)
	if 0
		l = linspace(0,2,101)';
		n = trl_curvature_opt(2+0*l, 3, 3, l);
		plot(l, n, '-o'), xlabel l, ylabel n
	else
		y = 4; b = 3; r = 1; ln = 3.8;
	%	y = 2; b = 3; r = 3; ln = 0.05;
		n = trl_curvature_opt(y, b, r, ln);
		l = linspace(-0.2,5,101);
		hn = h(y,b,r,ln);
		dhn = dh(y, b, r, ln);%
		ql = hn + dhn * (l-ln) - 0.5 * n * (l-ln).^2;
		hl = h(y,b,r,l);
		clf, plot(l, hl, '-',  l, ql, '--',  ln, hn, 'o'), grid
		axis([minmax(l)' minmax(hl)']), legend('h(l)', 'q(l)')
	end
return
end

if ~isvar('yb')
	yb = bi .* exp(-li) + ri;
end

	if any(bi <= 0), error 'bi=0 not done', end

	ni_max = bi .* (1 - yi .* ri ./ (bi + ri).^2);	% curvature at l=0
	ni_max = max(ni_max, 0);
	ni = ni_max;

	if 0
		il = li <= 0;
	else % trick in C program due to numerical precision issues
		il = li < 0.1;
	end

	tmp = h(yi,bi,ri,li) - h(yi,bi,ri,0) - li .* dh(yi,bi,ri,li);
	i = ~il;
	ni(i) = 2 ./ li(i).^2 .* max(tmp(i),0);

	if any(ni > ni_max)
	%	plot([ni_max(:) ni(:) ni(:)>ni_max(:)])
		warning 'large ni'
	end
%	minmax(ni)
