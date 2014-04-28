 function ci = trl_curvature_pre(yi, bi, ri)
%function ci = trl_curvature_pre(yi, bi, ri)
%
%	Precomputed approximating parabola curvatures
%	for Poisson transmission model.
%	The minimum returned curvature will be zero.
%	This is compatible with trpl/trp_init_der02_sino() in aspire.
%
%	Copyright 2002-1-28	Jeff Fessler	The University of Michigan


%
%	the default is to demonstrate an example parabola!
%
if nargin == 0
	help(mfilename)
	h = inline('y.*log(b.*exp(-l)+r)-(b.*exp(-l)+r)', 'y','b','r','l');
	dh = inline('(1 - y ./ (b.*exp(-l)+r)) .* b.*exp(-l)', 'y','b','r','l');
	y = 3; b = 5; r = 1;
	ln = 1.4 * log(b / (y - r));
	n = trl_curvature_pre(y, b, r);
	l = linspace(-0.2,5,101);
	hn = h(y,b,r,ln);
	dhn = dh(y, b, r, ln);%
	ql = hn + dhn * (l-ln) - 0.5 * n * (l-ln).^2;
	hl = h(y,b,r,l);
	clf, plot(l, hl, '-',  l, ql, '--',  ln, hn, 'o'), grid
	axis([minmax(l)' minmax(hl)']), legend('h(l)', 'q(l)')
return
end

trl_check(yi, bi, ri);

if ~isvar('bi') || isempty(bi)
	bi = ones(size(yi));
end
if ~isvar('ri') || isempty(ri)
	ri = zeros(size(yi));
end

% ci = (yi-ri)^2 / yi, if yi > ri >= 0 and bi > 0
ii = (yi > ri) & (ri >= 0) & (bi > 0); % good rays
ci = zeros(size(yi));
ci(ii) = (yi(ii) - ri(ii)).^2 ./ yi(ii);
