  function inv1 = de_ftab_inv1(fit, s1, varargin)
%|function inv1 = de_ftab_inv1(fit, s1, [options])
%|
%| Build object that does 1D inverse of BH function for 1st material component,
%| (usually water), for conventional "water only" beam-hardening correction.
%|
%| in
%|	fit	strum	initialized by de_ftab_fit()
%|	s1	[N]	s_1 sample values for building table (ftab.sls.sl{1})
%|
%| option
%|	'll'	1	which component (default: 1)
%|	'show'	0|1	plot it?
%|
%| out
%|	inv1	strum
%|	methods:
%|		inv1.fun(hf1) 	map log values into corrected values
%|			(corresponding to line-integral of material density)
%|		inv1.plot(fit)	show fit
%|
%| Copyright 2008-09-28, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(fit, 'test'), fail('run de_ftab test'), return, end

arg.show = false;
arg.ll = 1; % which material, one of {1, ..., LL}
arg = vararg_pair(arg, varargin);

% although we need only the 1D array of s1 values,
% fit.fmfun requires a [n1 1 LL] array.  here the other components are zero.
LL = fit.LL;
sll = s1; sll(1,1,LL) = 0; % [n1 1 L]

% evaluate nonliner BH function for each of the M spectra
f1 = fit.fmfun(sll); % [n1 1 M]
f1 = squeeze(f1(:,1,:)); % [n1 M]

st.s1 = s1;
st.f1 = f1;

meth = {'fun', @de_ftab_inv1_fun, '(fh1 [N M]) -> [N M]';
	'plot', @de_ftab_inv1_plot, '(fit)'};

inv1 = strum(st, meth);

if arg.show && im
	inv1.plot(fit);
end


%
% de_ftab_inv1_fun()
% inverse based on interpolation
%
function s1hat = de_ftab_inv1_fun(st, fh1)

s1 = st.s1;
f1 = st.f1;

MM = size(f1,2);
dim = size(fh1);
fh1 = reshapee(fh1, [], MM);
s1hat = zeros(size(fh1), 'single');
for mm=1:MM
	s1hat(:,mm) = interp1(f1(:,mm), s1, fh1(:,mm), 'cubic', 'extrap');
end
s1hat = reshape(s1hat, dim);


%
% de_ftab_inv1_plot()
%
function out = de_ftab_inv1_plot(st, fit)

if nargin ~= 2, fail 'need "fit" argument', end

s1 = st.s1;

s4 = linspace(0, max(s1), 4*length(s1)+1)'; % [n4] fine sampling
LL = fit.LL;
MM = fit.MM;
sll = s4; sll(1,1,LL) = 0; % [n4 1 L]

f4 = fit.fmfun(sll); % [n4 1 M]
f4 = squeeze(f4(:,1,:)); % [n4 M]

s4i = st.fun(f4); % [n4 M] one correction for each spectrum
err = s4i - repmat(s4, [1 MM]); % [n4 M]

clf, pl = @(n,m) subplot(200 + MM*10 + (n-1)*MM + m);
for mm=1:MM
	pl(1,mm)
	plot(s4, f4(:,mm), '-', s4i(:,mm), f4(:,mm), '.', ...
		s4, s4 * fit.mac_eff(mm), '--') % monenergetic line
	m = num2str(mm);
	axis tight, xlabel 's1', ylabel(['f' m]), title(['inv1 fit m=' m])

	pl(2,mm)
	plot(s4, err(:,mm), '.-')
	axis tight, xlabel 's1', ylabel(['err' m]), title 'inv1 error'
end
prompt

if nargout, out = []; end
