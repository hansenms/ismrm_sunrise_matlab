 function [basis, Hw] = deconv_z(psf, fun, t)
% psf: 3-tap symmetric psf: [a _b_ a]
% H(z) = 1 / (az + b + a/z)

if length(psf) == 1
	basis = fun(t) / psf;
	Hw = sprintf('%g * ones(size(t))', psf);
	Hw = inline(Hw, 't');
return
end

if length(psf) ~= 3, error 'not done', end
a = psf(1);
b = psf(2);
c = b / a / 2;
p = -c + sign(c) * sqrt(c^2 - 1); % pole

scale = 1/a * 1/(p - 1/p);
basis = fun(t);
for n=1:9
	basis = basis + p^n * (fun(t-n) + fun(t+n));
end
basis = scale * basis;

Hw = sprintf('1 ./ (%g + 2 * %g * cos(om))', b, a);
Hw = inline(Hw, 'om');
