  function ob = Gdelaysum(varargin)
%|function ob = Gdelaysum([options])
%|
%| Construct Gdelaysum object, which performs weighted sums of a delayed signal.
%| This is useful for model-based reconstruction of THz images.
%| See Gdelaysum_test() at end for example usage.
%|
%| y[n;m] = sum_{k=0}^{Nx-1} x[k] h[n - d[k;m]], n=0,...,Ny-1, m=0,...,Nm-1
%| for 0 <= n - d[k] <= Nh - 1
%|
%| required
%|	'Ny'	[1]		output signal length
%|	'delay'	[Nx Nm]		delays
%|
%| options
%|	'h'	[Nh 1]		impulse response, default [1]
%|	'nthread' [1]		# threads, default 1
%|
%| out
%|	ob	[Ny*Nm Nx]	Fatrix object
%|
%| Copyright 2006-8-25, Jeff Fessler, University of Michigan

if nargin == 1 && streq(varargin{1}, 'test'), Gdelaysum_test, return, end
if nargin < 4, help(mfilename), error(mfilename), end

% defaults
arg.Ny = [];
arg.delay = [];
arg.h = 1;
arg.nthread = 1;
arg.class = 'fatrix2';
arg.chat = false;

% options
arg = vararg_pair(arg, varargin);

if isempty(arg.delay), fail 'delay required', end
if isempty(arg.Ny), fail 'Ny required', end

arg.h = single(arg.h);
arg.delay = single(arg.delay);
arg.nthread = int32(arg.nthread);
arg.chat = int32(arg.chat);

[arg.Nx arg.Nm] = size(arg.delay);
arg.Nh = length(arg.h);

arg.str_forw = 'delaysum1,forw';
arg.str_back = 'delaysum1,back';
if arg.nthread > 1 
	arg.str_forw = 'delaysum1,forw,thr';
	arg.str_back = 'delaysum1,back,thr';
end

switch arg.class
case 'Fatrix'
	arg.dim = [arg.Ny * arg.Nm, arg.Nx];
	ob = Fatrix(arg.dim, arg, ...
		'forw', @Gdelaysum_forw_Fatrix, 'back', @Gdelaysum_back_Fatrix);
case 'fatrix2'

	odim = arg.Ny*arg.Nm;
	forw = @(arg, x) delaysum1_mex(arg.str_forw, arg.h, arg.delay, ...
		arg.nthread, single(x), int32(arg.Ny), arg.chat);
	back = @(arg, y) delaysum1_mex(arg.str_back, arg.h, arg.delay, ...
		arg.nthread, single(y), arg.chat);
	ob = fatrix2('idim', arg.Nx, 'odim', odim, 'arg', arg, ...
		'forw', forw, 'back', back);
otherwise
	fail 'bug'
end


% Gdelaysum_forw_Fatrix(): y = A * x
% in
%	x	[Nx L]
% out
%	y	[Ny*Nm L]
%
function y = Gdelaysum_forw_Fatrix(arg, x)

LL = size(x, 2);
y = zeros(arg.Ny*arg.Nm, LL, 'single');
for ll=1:LL
	tmp = single(x(:,ll));
	tmp = delaysum1_mex(arg.str_forw, arg.h, arg.delay, ...
		arg.nthread, tmp, int32(arg.Ny), arg.chat);
	y(:,ll) = tmp(:);
end


% Gdelaysum_back_Fatrix(): x = A' * y
% in
%	y	[Ny*Nm L]
% out
%	x	[Nx L]
%
function x = Gdelaysum_back_Fatrix(arg, y)

LL = size(y, 2);
x = zeros(arg.Nx, LL, 'single');
for ll=1:LL
	tmp = single(y(:,ll));
	tmp = delaysum1_mex(art.str_back, arg.h, arg.delay, ...
		arg.nthread, tmp, arg.chat);
	x(:,ll) = tmp;
end


% Gdelaysum_test
function Gdelaysum_test

delay = [10 20 40; 50 51 52];
h = [4 3 2 1]';
Ny = 80;
ats = {'nthread', jf('ncore')};
A = Gdelaysum1('Ny', Ny, 'h', h, 'delay', delay, ats{:});

x = [1 0]';
y1 = A * x;
y1 = reshapee(y1, Ny, []);
if im
	clf, plot(y1, '-o')
end

%mask = true(A.arg.Nx,1);
%Fatrix_test_basic(A, mask)
fatrix2_tests(A)
test_adjoint(A);

if 1 % compare Fatrix and fatrix2
	Ny = 2^11;
	Nx = 2^10;
	Nm = 2^9;
	Nh = 2^8;
	rng(0)
	delay = rand(Nx, Nm);
	x = rand(Nx, 1);
	args = {'Ny', Ny, 'h', (1:Nh)', 'delay', delay, ats{:}};
	A1 = Gdelaysum1(args{:}, 'class', 'Fatrix');
	A2 = Gdelaysum1(args{:}, 'class', 'fatrix2');
	cpu etic
	y1 = A1 * x;
	cpu etoc 'Fatrix'
	cpu etic
	y2 = A2 * x;
	cpu etoc 'fatrix2'
	jf_equal(y1, y2)
end
