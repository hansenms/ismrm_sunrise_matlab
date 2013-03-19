  function M = qpwls_precon(type, sys, C, mask, varargin)
%|function M = qpwls_precon(type, sys, C, mask, varargin)
%| usage:
%| M = qpwls_precon('circ0', {T}, C, mask);
%| M = qpwls_precon('circ0', {A, W}, C, mask);
%|
%| build a preconditioner for QPWLS problems;
%| approximations to inv([A'WA + C'C])
%|
%| in
%|	type	string		'circ0' : circulant based on center of FOV
%|				'dcd0' : diagonal * circulant * diagonal
%|	sys	cell		{T} or {A, W}
%|	A	[nd np]		system matrix
%|	W	[nd nd]		data weighting matrix
%|	C	[nc np]		penalty 'derivatives' (R = \Half C'*C)
%|	mask	[nx ny]		which pixels are updated
%|
%| options
%|	kappa	[nx ny]		needed for dcd0
%|	apod	[nx ny]		optional apodizer for A'WAe
%|	'class'	''		'fatrix2' or 'Fatrix' (default)
%|
%| out
%|	M	[np np]		fatrix2 (todo) or Fatrix (default) object
%|
%| The 'dcd0' preconditioner is based on Fessler&Booth IEEE T-IP May 1999
%|
%| Copyright 2004-6-29, Jeff Fessler, University of Michigan

if nargin < 4, help(mfilename), error(mfilename), end
arg.chat = false;
arg.kappa = [];
arg.apod = [];
arg.class = 'Fatrix';
arg = vararg_pair(arg, varargin);


dim = [1 1] * numel(mask);
switch type
case 'circ0' % circulant preconditioner
	arg.mask = mask;
	arg.Mdft = qpwls_precon_circ_Mdft(sys, C, mask, arg.apod, arg.chat);
	switch arg.class
	case 'Fatrix'
		arg.st = strum(arg, {'plot', @qpwls_precon_plot, '()'});
		M = Fatrix(dim, arg, 'forw', @qpwls_precon_circ_mult);
	otherwise
		fail 'class'
	end

case 'dcd0' % diagonal / circulant / diagonal
	arg.mask = mask;
	arg = qpwls_precon_dcd0_init(sys, C, arg);
	switch arg.class
	case 'Fatrix'
		M = Fatrix(dim, arg, 'forw', @qpwls_precon_dcd0_mult);
	otherwise
		fail 'class'
	end

case 'diag' % diagonal
	fail 'todo: bug jeff'

otherwise
	fail('unknown preconditioner "%s"', type)
end


% qpwls_precon_plot()
% show spectrum and psf
function psf = qpwls_precon_plot(st)

[nx ny] = size(st.mask);
tmp = zeros(nx, ny);
tmp(nx/2+1,ny/2+1) = 1;
ix = [-nx/2:nx/2-1];
iy = [-ny/2:ny/2-1];

im plc 1 2
im(1, ix, iy, fftshift(st.Mdft)), cbar
xtick([-nx/2 0 nx/2-1])
ytick([-ny/2 0 ny/2-1])

psf = qpwls_precon_circ_mult(st, tmp(st.mask));
psf = embed(psf, st.mask);
%psf = reale(fftshift(ifftn(st.Mdft))); % same
im(2, ix, iy, psf), cbar
xtick([-nx/2 0 nx/2-1])
ytick([-ny/2 0 ny/2-1])


% qpwls_precon_dcd0_init()
% initialize dcd preconditioner based on center pixel
function arg = qpwls_precon_dcd0_init(sys, R, arg)
if length(sys) ~= 2, error 'need cell(2)', end
arg.diag = 1 ./ arg.kappa(arg.mask(:));

ej = qpwls_precon_e0(arg.mask);
ctc = R.cgrad(R, 1e-2 * ej) / 1e-2; % trick
ctc = embed(ctc, arg.mask);
ctc = ctc / arg.kappa(end/2+1,end/2+1)^2; % trick

sys = {sys{1}, 1}; % trick
arg.Mdft = qpwls_precon_circ_Mdft(sys, ctc, arg.mask, [], arg.chat);


% qpwls_precon_dcd0_mult()
% multiply using diag * circ * diag
function y = qpwls_precon_dcd0_mult(arg, x)

x = x .* arg.diag;
y = ifftn_fast(arg.Mdft .* fftn_fast(embed(x, arg.mask)));
y = y(arg.mask(:));
y = y .* arg.diag;
if isreal(x)
	y = reale(y, 'warn');
end



% setup circulant preconditioner based on center pixel
function Mdft = qpwls_precon_circ_Mdft(sys, C, mask, apod, chat)

ej = qpwls_precon_e0(mask);

% T * x or A'WA*x
if ~iscell(sys), error 'sys must be cell', end
if length(sys) == 2
	A = sys{1};
	W = sys{2};
	awa = A' * (W * (A * ej));
elseif length(sys) == 1
	awa = sys{1} * ej; % T * ej
else
	error 'unknown cell'
end

awa = embed(awa, mask);
if ~isempty(apod)
%	im(awa, 'awa'), cbar, prompt
	awa = awa .* apod;
%	im(awa, 'awa apodized'), cbar, prompt
end

if isnumeric(C) && isequal(size(C), size(awa)) % trick
	ccc = C;

elseif isstruct(C) % 'R'
	R = C;
	ccc = R.cgrad(R, 1e-2 * ej) / 1e-2; % trick
	ccc = embed(ccc, mask);

else
	ccc = C' * (C * ej);
	ccc = embed(ccc, mask);
end
%im(ccc, 'ccc'), cbar, prompt

f.awa = fftn_fast(fftshift(awa));
f.ccc = fftn_fast(fftshift(ccc));
f.ccc = reale(f.ccc, 'warn'); % these should be nearly real
if any(f.ccc(:) < - 1e-6 * max(f.ccc(:)))
	printm('ccc min = %g max = %g', min(f.ccc(:)), max(f.ccc(:)))
	clf, im(ccc), keyboard
	error 'bug: circulant penalty is not nonnegative definite!?'
end
f.ccc = max(f.ccc, 0);
f.awa = reale(f.awa, 'warn');
if min(f.awa(:)) < 0
	printm('setting %g%% to zero', ...
		min(f.awa(:)) / max(f.awa(:)) * 100)
	f.awa = max(f.awa, 0);
end
%minmax(f.awa), minmax(f.ccc)
f.h = f.awa + f.ccc;	% approximate hessian in freq. domain
if min(f.h(:)) <= 0
	error 'circulant preconditioner has zero!?  do you regularize?'
end
if chat
	printm('approximate condition number: %g', ...
		max(f.h(:)) / min(f.h(:)))
end

Mdft = 1 ./ f.h;



% qpwls_precon_circ_mult()
% multiply using fft
function y = qpwls_precon_circ_mult(arg, x)

y = ifftn_fast(arg.Mdft .* fftn_fast(embed(x, arg.mask)));
y = y(arg.mask(:));
if isreal(x)
	y = reale(y, 'warn');
end


% qpwls_precon_e0()
% unit impulse at "center" voxel
function e0 = qpwls_precon_e0(mask)

e0 = zeros(size(mask));
switch ndims(mask)
case 2
	if size(mask,2) == 1 % 1d
		e0(end/2+1) = 1;
	else % 2d
		e0(end/2+1,end/2+1) = 1;
	end
case 3
	e0(end/2+1,end/2+1,end/2+1) = 1;
otherwise
	fail 'not done'
end
e0 = e0(mask(:));
