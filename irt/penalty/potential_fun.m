  function pot = potential_fun(type, delta, param, varargin)
%|function pot = potential_fun(type, delta, param, varargin)
%|
%| Define roughness penalty potential functions (as strum).
%|
%| The penalty will have the form
%|	R(x) = sum_k w_k * potential_k([Cx]_k, delta_k)
%| where w_k is provided elsewhere, not here!
%|
%| in
%|	type		quad broken huber hyper2 hyper3 cauchy qgg2 genhub
%|			lange1 lange3 (fair) geman&mcclure gf1 li98cfs
%|			Recommended: 'hyper3'
%|			Use potential_fun('list') to list all choices.
%|	delta		scalar, or image-sized array;
%|			"cutoff" parameter for edge-preserving regularization
%|	param		optional additional parameter(s) for some choices:
%|				'gf1' (generalized Fair) : [a b]
%|				'qgg2' : q
%|				'genhub' & 'stevenson94dpr' : [p q]
%|				'table1' : [dt, dpot([1:K] * dt)]
%|
%| option
%|	'dummy'		include extra dummy argument for backward compatibility
%|			e.g., pot.potk([], C*x).  default: 0
%|
%| out
%|	pot		strum object, with data: delta and param
%|	methods:
%|		pot.potk(C*x)	potential function value
%|		pot.wpot(C*x)	potential 'weights' (aka half-quad. curvatures)
%|		pot.dpot(C*x)	potential derivative
%| 		pot.shrink(z, reg)	shrinkage operation:
%|					argmin_x 1/2 |z-x|^2 + reg * pot(x)
%|		pot.plot()	plot all of the above functions
%|
%| trick: for now these methods all require an extra dummy argument
%| for compatability with old version as follows:
%|	pot.potk(nan, C*x)
%|
%| trick: for type 'gf1-fit', the param argument should be:
%|	{'name of potential to fit', points, param}
%| and this function returns the [a b] parameters needed for a subsequent
%| call with type 'gf1'
%|
%| Copyright 2004-5-18, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if streq(type, 'test') potential_fun_test, return, end

if streq(type, 'list') % return list of all choices
	pot = {'quad', 'broken', 'huber', 'hyper2', 'hyper3', 'cauchy', ...
	'qgg2', 'genhub', 'lange1', 'lange3', 'geman&mcclure', 'gf1', ...
	'table0', 'table1', 'li98cfs'};
%	'stevenson94dpr' omit this one; use 'genhub' instead
return
end

if streq(type, 'list1') % return list of 'l1' type choices (mostly non-smooth)
	pot = {'l0', 'l1', 'tav', 'fair-l1'};
return
end

if ~isvar('delta'), delta = []; end
if ~isvar('param'), param = []; end

arg.dummy = false;
arg.scale = 1;
arg = vararg_pair(arg, varargin);

if streq(type, 'gf1-fit') % trick: just return parameters for this case
	pot = potential_fun_gf1_fit(param{:});
return
end

[pot.type pot.delta pot.param potk wpot dpot shrink] = ...
	potential_fun_parse(type, delta(:), param(:));

if arg.scale ~= 1
	potk = @(pot, t) arg.scale * potk(pot, t);
	dpot = @(pot, t) arg.scale * dpot(pot, t);
	wpot = @(pot, t) arg.scale * wpot(pot, t);
end

if arg.dummy
	potk = @(pot, dum, t) potk(pot, t);
	dpot = @(pot, dum, t) dpot(pot, t);
	wpot = @(pot, dum, t) wpot(pot, t);
end

meth = {'potk', potk, 'wpot', wpot, 'dpot', dpot, 'shrink', shrink, ...
	'plot', @ir_potential_fun_plot};
pot = strum(pot, meth);

end % potential_fun()


% potential_fun_gf1_fit()
% trick: gf1-fit means match gf1 to named potential at given points [0, s1, s2]
% where the given point values are s = t / delta, i.e., in "delta" units
% param: {'name of potential to fit', points, param}
function ab = potential_fun_gf1_fit(type, sv, param)
sv = sv(:);
pot = potential_fun(type, 1, param); % trick: delta=1 here
pt = pot.wpot(sv);

% ab = [1 1; -pt']' \ ((pt-1) ./ sv);

s1 = sv(1);
s2 = sv(2);
w1 = pt(1);
w2 = pt(2);
ab(1) = (w2 * s2 * (1 - w1) - w1 * s1 * (1 - w2)) ...
		/ ((w1 - w2) * s1 * s2);
ab(2) = (s2 * (1 - w1) - s1 * (1 - w2)) ...
		/ ((w1 - w2) * s1 * s2);
end % potential_fun_gf1_fit()


% potential_fun_parse()
function [type, delta, param, potk, wpot, dpot, shrink] = ...
	potential_fun_parse(type, delta, param)

dpot = [];
shrink = [];

% trick: huber2 is just a huber with delta / 2
% so that weighting function drops to 1/2 at delta, like hyper3 etc.
if streq(type, 'huber2')
	type = 'huber';
	delta = delta / 2;
end

% trick: hyper3 is just a hyperbola with delta scaled by sqrt(3)
% to approximately "match" the properties of 'cauchy' (old erroneous 'hyper')
if streq(type, 'hyper3')
	type = 'hyper2';
	delta = delta / sqrt(3);
end

if streq(type, 'gf1') && param(1) == 0
	if param(2) ~= 1
		fail 'only b=1 makes sense for gf1 with a=0'
	end
	type = 'lange3'; param = []; % trick: gf1 with a=0 and b=1 is lange3
end

switch type

% quadratic potential function
case 'quad'
	potk = @(pot, t) (abs(t).^2) / 2;
	wpot = @(pot, t) ones(size(t));
	dpot = @(pot, t) t;
	shrink = @(pot, z, reg) z ./ (1 + reg);

% broken parabola
case 'broken'
	potk = @(pot, t) min(t.^2, pot.delta.^2)/2;
	wpot = @(pot, t) abs(t) < pot.delta;
	dpot = @(pot, t) t .* (abs(t) < pot.delta);
	shrink = @(pot, z, reg) broken_shrink(z, reg, pot.delta);

% l0 (hard threshold)
case 'l0'
	potk = @(pot, t) tonumeric(t ~= 0, t);
	wpot = @(pot, t) nan(size(t), class(t));
	dpot = @(pot, t) nan(size(t), class(t));
	shrink = @(pot, z, reg) z .* (abs(z) > sqrt(2*reg));

% absolute value (i.e., l1, for soft threshold)
case {'abs', 'l1'}
	potk = @(pot, t) abs(t);
	wpot = @(pot, t) 1 ./ abs(t); % invalid at t=0
	dpot = @(pot, t) sign(t);
	shrink = @(pot, z, reg) sign(z) .* max(abs(z) - reg, 0);

% truncated absolute value
case 'tav'
	potk = @(pot, t) min(abs(t), pot.delta);
	wpot = @(pot, t) (1 ./ abs(t)) .* (abs(t) < pot.delta); % bad at t=0
	dpot = @(pot, t) sign(t) .*  (abs(t) < pot.delta);
	shrink = @(pot, z, reg) tav_shrink(z, reg, pot.delta);

% huber potential function
case 'huber'
	potk = @(pot, t) huber_pot(t, pot.delta);
	wpot = @(pot, t) huber_wpot(t, pot.delta);
	dpot = @(pot, t) huber_dpot(t, pot.delta);
	shrink = @(pot, z, reg) huber_shrink(z, reg, pot.delta);

% cauchy penalty: d^2 / 2 * log(1 + (t/d)^2) (not convex!)
case 'cauchy'
	potk = @(pot, t) pot.delta.^2 / 2 .* log(1 + abs(t ./ pot.delta).^2);
	wpot = @(pot, t) 1 ./ (1 + abs(t ./ pot.delta).^2);
	dpot = @(pot, t) t ./ (1 + abs(t ./ pot.delta).^2);
	shrink = @(pot, z, reg) cauchy_shrink(z, reg, pot.delta);

% Geman&McClure penalty: d^2 / 2 * (t/d)^2 / (1 + (t/d)^2)
% Not convex!
case 'geman&mcclure'
	potk = @(pot, t) pot.delta.^2 / 2 .* (t/pot.delta).^2 ./ (1 + abs(t ./ pot.delta).^2);
	wpot = @(pot, t) 1 ./ (1 + abs(t ./ pot.delta).^2).^2;
	dpot = @(pot, t) t ./ (1 + abs(t ./ pot.delta).^2).^2;

% gf1: Generalized Fair 1st-order
% wpot(t) = (1 + a * |t/d|) / (1 + b * |t/d|)
case 'gf1'
	potk = @(pot, t) gf1_potk(t, pot.delta, pot.param(1), pot.param(2));
	wpot = @(pot, t) (1 + pot.param(1) .* abs(t ./ pot.delta)) ...
		./ (1 + pot.param(2) .* abs(t ./ pot.delta));
	shrink = @(pot, z, reg) ...
		gf1_shrink(z, reg, pot.delta, pot.param(1), pot.param(2));


% hyperbola penalty: d^2 * [ sqrt(1 + (t/d)^2) - 1 ]
case 'hyper2'
	potk = @(pot, t) pot.delta.^2 .* (sqrt(1 + abs(t ./ pot.delta).^2) - 1);
	wpot = @(pot, t) 1 ./ sqrt(1 + abs(t ./ pot.delta).^2);
	dpot = @(pot, t) t ./ sqrt(1 + abs(t ./ pot.delta).^2);

case 'hyper'
	error 'use "cauchy" or "hyper3" not "hyper" now'

% Lange1 penalty
case 'lange1'
	potk = @(pot, t) t.^2 / 2 ./ (1+abs(t./pot.delta));
	wpot = @(pot, t) (1 + abs(t ./ pot.delta) / 2) ./ (1 + abs(t ./ pot.delta)).^2;

% Lange3 penalty
case {'lange3', 'fair'}
	potk = @(pot, t) pot.delta.^2 .* (abs(t./pot.delta) - log(1+abs(t./pot.delta)));
	wpot = @(pot, t) 1 ./ (1 + abs(t ./ pot.delta));
	dpot = @(pot, t) t ./ (1 + abs(t ./ pot.delta));

% Fair potential "rounded corner" approximation to l1
case 'fair-l1'
	potk = @(pot, t) abs(t) - pot.delta .* log(1+abs(t./pot.delta));
	wpot = @(pot, t) 1 ./ (pot.delta + abs(t));
	dpot = @(pot, t) t ./ (pot.delta + abs(t));
	shrink = @(pot, z, reg) fair_l1_shrink(z, reg, pot.delta);

% li98cfs
case 'li98cfs'
	% f = inline('atan(x) / x - 0.5'); fsolve(f, 2.3)
	delta = delta / 2.3311;
	potk = @(pot, t) li98cfs_potk(t, pot.delta);
	wpot = @(pot, t) li98cfs_wpot(t, pot.delta);

% qgg2: q-generalized gaussian for p=2, due to Thibault, Sauer, Bouman
% q = "param", same as lange1 when q=1
case 'qgg2'
	potk = @(pot, t) t.^2 / 2 ./ (1+abs(t./pot.delta).^(2-pot.param));
	wpot = @(pot, t) (1 + abs(t ./ pot.delta).^(2-pot.param) * pot.param ...
		 / 2) ./ (1 + abs(t ./ pot.delta).^(2-pot.param)).^2;

% genhub : generalized Huber (switch between two generalized gaussians)
% same as Huber when p=2 and q=1
% p is power near 0, q is asymptotic power
case 'genhub'
	potk = @(pot, t) genhub_potk(t, pot.delta, pot.param(1), pot.param(2));
	wpot = @(pot, t) genhub_wpot(t, pot.delta, pot.param(1), pot.param(2));

% from stevenson:94:dpr
% p = param(1), q = param(2), same as Huber when p=2 and q=1 ???
case 'stevenson94dpr'
	warn('Potential %s does not work well; use "genhub" instead', type)
	potk = @(pot, t) stevenson94dpr_potk(t, pot.delta, pot.param(1), pot.param(2));
	wpot = @(pot, t) ones(size(t)); % bogus attempt at upper bound

% tabulate derivative of pot at t_k = dt * k for k=1,...,K
% param is samples of derivative dpot(t_k)
% dpot(0) = 0, and use sample-and-hold interpolation of dpot()
% except use linear interpolation over [0 dt]
case 'table0'
	if ~isnumeric(param) || numel(param) < 10
%	if ~iscell(param) || numel(param) ~= 2
%		fail 'for table0 param must be {dt, dpot([1:K]*dt}'
		fail 'for table0 param must be [dt, dpot([1:K]*dt)]'
	end
%	param = table0_setup(param{1}, col(param{2}));
	param = table0_setup(param(1), col(param(2:end)));
	potk = @(pot, t) table0_potk(t, pot.param);
	wpot = @(pot, t) table0_wpot(t, pot.param);
	shrink = @(pot, z, reg) table0_shrink(z, reg, pot.param);

% tabulate derivative of pot at t_k = dt * k for k=1,...,K
% param is samples of derivative dpot(t_k)
% dpot(0 = 0, and use linear interpolation of dpot(t)
case 'table1'
	if ~isnumeric(param) || numel(param) < 10
%	if ~iscell(param) || numel(param) ~= 2
%		fail 'for table1 param must be {dt, dpot([1:K]*dt}'
		fail 'for table1 param must be [dt, dpot([1:K]*dt)]'
	end
%	param = table1_setup(param{1}, col(param{2}));
	param = table1_setup(param(1), col(param(2:end)));
	potk = @(pot, t) table1_potk(t, pot.param);
	wpot = @(pot, t) table1_wpot(t, pot.param);
	shrink = @(pot, z, reg) table1_shrink(z, reg, pot.param);

otherwise
	fail('Unknown potential "%s"', type)
end

if isempty(dpot) % default is t * wpot(t)
	dpot = @(pot, t) t .* wpot(pot, t);
end
if isempty(shrink)
	shrink = @potential_fun_shrink;
end

end % potential_fun_parse()


% potential_fun_shrink()
% default method uses fzero() which will be slow!
function out = potential_fun_shrink(pot, z, reg)
out = zeros(size(z));
for ii=1:numel(z)
	zi = z(ii);
	fun = @(t) t + reg * pot.dpot(t) - zi;
	out(ii) = fzero(fun, zi);
end
end % potential_fun_shrink


% table0_setup()
% dt	[1]	t spacing
% dk	[K]	dpot([1:K] * dt)
function out = table0_setup(dt, dk)
sk = dk(1)*dt/2 + dt * [0; cumsum(dk(1:end-1))]; % [K]
out.dt = dt;
out.dk = dk; % [1:K] samples of dpot
out.sk = sk; % [K] cumulative sums for pot(t)
K = numel(dk);
out.K = K;
end % table0_setup


% table0_potk()
function out = table0_potk(t, param)
dt = param.dt;
dk = param.dk;
sk = param.sk;
t = abs(t);
k = floor(t / dt);
k = min(k, param.K); % at most K
big = t > dt;
out(~big) = 0.5 * dk(1) / dt * t(~big).^2;
k = k(big);
out(big) = sk(k) + dk(k) .* (t(big) - k * dt);
end % table0_potk


% table0_wpot()
function out = table0_wpot(t, param)
dt = param.dt;
dk = param.dk; % [K]
t = abs(t);
k = floor(t / dt);
k = min(k, param.K); % at most K
out = repmat(dk(1) / dt, size(t));
big = k > 0;
k = k(big);
t = t(big);
out(big) = dk(k) ./ t;
end % table0_wpot


% table0_shrink()
function out = table0_shrink(z, reg, param)
K = param.K;
dk = param.dk;
dt = param.dt;
tk = (1:K)' * dt;
bk = tk + reg * dk; % must be done here because depends on reg
ck = bk + dt;
tmp = [0; col([bk ck]')];
out = [tk tk+dt];
out = [0; col([tk tk+dt]')];
out = sign(z) .* interp1(tmp, out, abs(z), 'linear', 'extrap');
end % table0_shrink


% table1_setup()
% dt	[1]	t spacing
% dk	[K]	dpot([1:K] * dt)
function out = table1_setup(dt, dk)
K = numel(dk);
k = [0:K]'; % [K+1]
tk = dt * k; % [K+1]
dk0 = [0; dk]; % [K+1] prepend sample at zero
ck = [diff(dk0); 0] / dt; % [K+1] curvatures 0:K
tmp = (dk0 - tk .* ck) * dt ...
	+ dt^2/2 * ck .* ((k+1).^2 - k.^2); % [K+1]
sk = cumsum( [0; tmp] );
out.dt = dt;
out.dk = dk; % [1:K] samples of dpot
out.ck = ck; % [K+1] curvatures 0:K
out.sk = sk; % [K] cumulative sums for pot(t)
out.K = K;
end % table1_setup


% table1_potk()
function out = table1_potk(t, param)
sk = param.sk;
dt = param.dt;
ck = param.ck;
dk0 = [0; param.dk]; % [K+1] prepend sample at zero
t = abs(t);
k = floor(t / dt);
k = min(k, param.K); % at most K
sk = sk(1+k); % matlab indexing
ck = ck(1+k); % matlab indexing
dk = dk0(1+k);
tk = k * dt;
ck = reshape(ck, size(t));
dk = reshape(dk, size(t));
sk = reshape(sk, size(t));
out = sk + (dk - tk .* ck) .* (t - tk) + ck / 2 .* (t.^2 - tk.^2);
end % table1_potk


% table1_wpot()
function out = table1_wpot(t, param)
dt = param.dt;
ck = param.ck;
dk = param.dk;
t = abs(t);
k = floor(t / dt);
k = min(k, param.K); % at most K
out = repmat(dk(1) / dt, size(t));
big = k > 0;
k = k(big);
t = t(big);
dk = reshape(dk(k), size(t));
ck = reshape(ck(k+1), size(t));
out(big) = (dk + (t - k * dt) .* ck) ./ t;
end % table1_wpot


% table1_shrink()
% unfortunately this routine works only for a single scalar reg value.
function out = table1_shrink(z, reg, param)
K = param.K;
tk = (1:K)' * param.dt;
bk = tk + reg * param.dk; % must be done here because depends on reg
out = sign(z) .* interp1([0; bk], [0; tk], abs(z), 'linear', 'extrap');
end % table1_shrink


% gf1_potk()
% gf1: generalized fair 1st-order potential
function pot = gf1_potk(t, delta, a, b)
atd = abs(t ./ delta);
pot = delta.^2 ./ (2 * b.^3) * (...
	2 * b.^2 .* atd + a .* b.^2 .* atd.^2 ...
	- 2 * a .* b .* atd + 2 * (a-b) .* log(1 + b .* atd));

if 0 % symbolic check
	syms x
	syms a positive
	syms b positive
	syms t positive
	int(x*(1+a*x) / (1+b*x), x, 0, t)
end
end % gf1_potk



% broken_shrink()
function out = broken_shrink(z, reg, delta)
out = z ./ (1 + reg);
big = delta * (1 + reg) < abs(z);
out(big) = z(big);
end % broken_shrink


% cauchy_shrink()
function out = cauchy_shrink(z, reg, delta)
z = z(:);
coef = ones(numel(z),1);
coef = [1/delta^2*coef -z/delta^2 (1+reg)*coef -z];
out = zeros(size(z));
for ii=1:numel(z)
	tmp = roots(coef(ii,:));
	pick = tmp == real(tmp); % empirically, 3rd root is often real
	if sum(pick) ~= 1 % if multiple real roots, empirically pick largest
		pick = imax(abs(tmp));
	end
	out(ii) = tmp(pick);
end
end % cauchy_shrink


% huber_shrink()
function out = huber_shrink(z, reg, delta)
out = z ./ (1 + reg);
big = delta * (1 + reg) < abs(z);
out(big) = z(big) .* (1 - reg .* delta ./ abs(z(big)));
end % huber_shrink


% fair_l1_shrink()
function out = fair_l1_shrink(z, reg, delta)
out = sign(z) .* (abs(z) - (delta + reg) ...
        + sqrt( (delta + reg - abs(z)).^2 + 4 * delta .* abs(z) )) ./ 2;
end % fair_l1_shrink


% gf1_shrink()
function out = gf1_shrink(z, reg, delta, a, b)
u = a / delta;
v = b / delta;
out = sign(z) .* (v .* abs(z) - (1+reg) ...
        + sqrt( (1 + reg - v.*abs(z)).^2 + 4 * (v + reg .* u) .* abs(z) )) ...
        ./ (2 * (v + reg .* u));
end % gf1_shrink

% tav_shrink()
function out = tav_shrink(z, reg, delta)
out = zeros(size(z), class(z));
big = reg < abs(z) & abs(z) < reg + delta;
out(big) = z(big) .* (1 - reg ./ abs(z(big)));
big = reg + delta <= abs(z);
out(big) = z(big);
end % tav_shrink


function pot = li98cfs_potk(t, d)
pot = d.^2 .* ((t ./ d) .* atan(t ./ d) - 0.5 * log(1 + (t ./ d).^2));
end

function pot = genhub_potk(t, d, p, q)
pot = 0.5 * abs(t) .^ p .* (abs(t) <= d) ...
 + 0.5 * (p ./ q .* d .^ (p-q) .* abs(t) .^ q ...
 + (1 - p ./ q) .* d .^ p) .* (abs(t) > d);
end

function pot = genhub_wpot(t, d, p, q)
pot = p / 2 .* (d .^ (p-q)) .* (abs(t) .^ (q-2)) .* (abs(t) > d);
ii = abs(t) <= d;
pot(ii) = p / 2 .* (abs(t(ii)) .^ (p-2));
%pot = p / 2 .* abs(t) .^ (p-2) .* (abs(t) <= d) ...
% + p / 2 .* d .^ (p-q) .* abs(t) .^ (q-2) .* (abs(t) > d);
end

function pot = stevenson94dpr_potk(t, d, p, q)
% pr [d p q]
tmp1 = (0.5 * abs(t) .^ p) .* (abs(t) <= d);
tmp2 = 0.5 * ( (p .* (d .^ (p-1)) .* abs(t) - p .* (d .^ p) ...
 + (1 ./ q) .^ (1 ./ (q-1)) ) .^ q ...
 + d .^ p - (1 ./ q) .^ (q ./ (q-1)) ) .* (abs(t) > d);
pot = tmp1 + tmp2;
end


% tonumeric(x, y)
% convert x to type of y
function out = tonumeric(x, y)
switch class(y)
case 'double'
	out = double(x);
case 'single'
	out = single(x);
otherwise
	fail('unknown type %s', class(y))
end
end % tonumeric


% ir_potential_fun_plot()
function dummy = ir_potential_fun_plot(pot)
t = linspace(-1,1,101)*2;
if isvar('pot.delta')
	t = t * pot.delta;
end
if ~im, return, end
im plc 2 2
im subplot 1
plot(t, pot.potk(t), '-', t, t.^2/2, ':', t, abs(t), ':')
axis([min(t) max(t) minmax(pot.potk(t))'])

im subplot 2
plot(t, pot.dpot(t), '-', t, t, ':')
axis([min(t) max(t) minmax(pot.dpot(t))'])

im subplot 3
plot(t, pot.wpot(t), '-', t, 1+0*t, ':')
axis([min(t) max(t) 0 1])

im subplot 4
reg = 1;
tmp = pot.shrink(t, reg);
plot(t, tmp, '-', t, t, ':')
axis([min(t) max(t) minmax(tmp)'])

dummy = [];
end % ir_potential_fun_plot


% potential_fun_test()
% test routine
% examine potential functions after rescaling.
function potential_fun_test

delta = 10; tmax = 4 * delta; reg = 1.5 * delta;
plist = potential_fun('list');
%plist = {'l1', 'fair-l1'}; delta = 0.5;
%plist = {'quad', 'li98cfs', 'hyper3', 'huber2'}; % show li98cfs roughly hyper3
%plist = {'quad', 'genhub', 'huber', 'stevenson94dpr'};
%plist = {'genhub'}
%plist = {'hyper3', 'qgg2', 'huber'}; delta = 10; tmax = 50;
%plist = {'lange3', 'qgg2'}; tmax = 200;
%plist = {'qgg2', 'gf1-fit'};
%plist = potential_fun('list1'); delta = 0.5;
%plist = {plist{:}, 'quad', 'gf1', 'huber', 'broken'}; % todo: more!
%plist = {'qgg2', 'table1', 'table0'};
%plist = {'qgg2', 'table1'}; fname = 'fig_reg_pot_table1_qgg2';
%plist = {'qgg2', 'table0'}; fname = 'fig_reg_pot_table0_qgg2';
%plist = {'qgg2', 'table1', 'gf1'}; fname = 'fig_reg_pot_table0_gf1_qgg2';
%plist = {'cauchy', 'l1'};
t = tmax * linspace(-1, 1, 2001)';
zs = 3 * max(delta + reg, delta * (1+reg)) * linspace(-1, 1, 301)';
ps = [];
lshrink = {};
for ii=1:length(plist)
	type = plist{ii}; % pr type
	if streq(type, 'quad')
		leg{ii} = type;
	else
		leg{ii} = [type ', \delta = ' num2str(delta)];
	end

	switch type
	case {'gf1', 'gf1-fit'}
		param = potential_fun('gf1-fit', nan, {'qgg2', [1 10], 1.2});
		type = 'gf1'; % trick
		leg{ii} = [leg{ii} sprintf(' %.3g %.4f', param(1), param(2))];
	case 'qgg2'
		param = 1.2;
		leg{ii} = [leg{ii} ', q = ' num2str(param)];
	case 'genhub'
		param = [2.0 1.2];
		leg{ii} = [leg{ii} sprintf('p=%g q=%g', param(1), param(2))];
	case 'stevenson94dpr'
		param = [2 1.2];
		leg{ii} = [leg{ii} sprintf('p=%g q=%g', param(1), param(2))];
	case {'table0', 'table1'}
		dt = delta / 20;
		tmp = potential_fun('qgg2', delta, 1.2);
%		param = {dt, tmp.dpot([1:1e4] * dt)};
		param = [dt, tmp.dpot([1:1e4] * dt)];
%		tmp = sprintf(' K = %d dt = %g', numel(param{2}), param{1});
		tmp = sprintf(', K = %d, dt = %g', numel(param)-1, param(1));
		leg{ii} = [leg{ii} tmp];
	otherwise
		param = [];
	end

	pot = potential_fun(type, delta, param);
	pp(:,ii) = pot.potk(t);
	pw(:,ii) = pot.wpot(t);
	pd(:,ii) = pot.dpot(t);

	% replace 0 with 1 to make (slow) figure showing fzero-based shrinkage
	if 0 || ~streq(func2str(pot.meth.shrink), 'potential_fun_shrink')
		tmp = pot.shrink(t, reg);
		if any(tmp ~= 0)
			ps(:,end+1) = pot.shrink(zs, reg);
			lshrink{end+1} = leg{ii};
		end
	end

	if 0 % test vs old
		opot = potential_func(type, delta, param);
		opp = opot.potk(opot, t);
		opw = opot.wpot(opot, t);
		opd = opot.dpot(opot, t);
		if any(opp ~= pp(:,ii)), 'p bug', type, keyboard, end
		if any(opw ~= pw(:,ii)), 'w bug', type, keyboard, end
		if any(opd ~= pd(:,ii)), 'd bug', type, keyboard, end
	end
end

if im
	set(0,'DefaultAxesLineStyleOrder', '-|:')
	clf
	subplot(411), plot(t, pp), title 'potk'
	axis tight, axisy([-0.0 2.5] * delta^2)
	legend(leg, 'location', 'north')
	subplot(412), plot(t, pw), title 'wpot'
	axis tight, axisy(0, 1.1), ytick([0 1])
	subplot(413), plot(t, pd), title 'dpot'
	axis tight, axisy([-1 1] * 1.3 * delta)
	xlabel 't'
%	savefig eps_c fname
	% check derivatives
	subplot(414)
	plot(t, pd)
	title 'dpot check'
	hold on
	d = diffc(pp) / (t(2)-t(1));
	plot(t(1:end-1), d(1:end-1,:), '--', ...
		t(1:end-1), d(1:end-1,:)-pd(1:end-1,:), ':')
	hold off
	axis tight, axisy([-1 1] * 1.3 * delta)
end

if 1 && im
	if ~isempty(lshrink)
		prompt
		clf
		plot(zs, ps, '-', zs, zs, '-')
		legend(lshrink{:}, 'location', 'southeast')
		xlabelf 'z'
%		ylabelf 'xhat(z)'
		ylabelf 't(z)'
		axis equal, axis square
		axis([-1 1 -1 1]*400)
		xtick([-1 0 1] * 400)
		ytick([-1 0 1] * 400)
%		savefig cw fig_reg_pot_table0_qgg2_shrink

		if 0 % figure for book
			plot(zs, ps(:,[3 2]) - repmat(ps(:,1), [1 2]), 'o')
			plot(	zs, ps(:,3) - ps(:,1), 'bo', ...
				zs, ps(:,2) - ps(:,1), 'gx')
			xlabel 'z', ylabel 'shrinkage error for QGG2'
			axis([0 500 -0.05 0.5]), ytick([0 0.05 0.1 0.5])
			legend(lshrink{[3 2]}, 1)
%			savefig eps_c fig_reg_pot_table01_qgg2_shrinker
			keyboard
		end
	end
end

% check 'scale' option
pot1 = potential_fun('lange3', 10, []);
pot2 = potential_fun('lange3', 10, [], 'scale', 2);
jf_equal(2 * pot1.potk(t), pot2.potk(t))
jf_equal(2 * pot1.dpot(t), pot2.dpot(t))
jf_equal(2 * pot1.wpot(t), pot2.wpot(t))
 
end % potential_fun_test()
