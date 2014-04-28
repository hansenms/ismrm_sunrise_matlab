  function y = ir_conv(x, psf, varargin)
%|function y = ir_conv(x, psf, varargin)
%|
%| ND convolution of signal x with psf
%|
%| in
%|	x	[(N)]	signal
%|	psf	[(N)]	psf
%|
%| option
%|	'adj'	0|1	default: 0 (if 1 then adjoint)
%|	'per'	0|1|[]	default: 0 (if 1 then periodic boundary conditions)
%|			if logical array with numel(dims(x)), then selectively
%|			periodic along corresponding directions
%|
%| out
%|	y	[(N)]	x * psf
%|
%| 2012-05-22, Jeff Fessler, Univ. of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if streq(x, 'test'), ir_conv_test, return, end

arg.per = false;
arg.adj = false;
arg = vararg_pair(arg, varargin);

if arg.adj
	psf = conj(flipdims(psf, 'odd', 1));
end

ndim = ndims(x);

switch arg.per
case false
	y = convn(x, psf, 'same');

case true
	y = ir_conv_per(x, psf, true(1,ndim));

otherwise
	y = ir_conv_per(x, psf, arg.per);
end


% ir_conv_pad()
function x = ir_conv_pad(x, psf, per)
ndim = ndims(x);
if ~islogical(per) || numel(per) ~= ndim
	fail('per must be logical and length %d', ndim)
end

zeross = @(a) zeros(size(a), class(a));

spad = floor(size(psf) / 2); % pad at start
epad = size(psf) - spad - 1; % pad at end

switch ndim
case 2

	if spad(1) || epad(1)
		pad = spad(1);
		tmp1 = x((end-pad+1):end,:);
		pad = epad(1);
		tmp2 = x(1:pad,:);

		if ~per(1)
			tmp1 = zeross(tmp1);
			tmp2 = zeross(tmp2);
		end

		x = cat(1, tmp1, x, tmp2);
	end

	if spad(2) || epad(2)
		pad = spad(2);
		tmp1 = x(:,(end-pad+1):end);
		pad = epad(2);
		tmp2 = x(:,1:pad);
		if ~per(2)
			tmp1 = zeross(tmp1);
			tmp2 = zeross(tmp2);
		end
		x = cat(2, tmp1, x, tmp2);
	end

case 3
	if ndims(psf) == 2
		spad = [spad 0];
		epad = [epad 0];
	end

	if spad(1) || epad(1)
		pad = spad(1);
		tmp1 = x((end-pad+1):end,:,:);
		pad = epad(1);
		tmp2 = x(1:pad,:,:);

		if ~per(1)
			tmp1 = zeross(tmp1);
			tmp2 = zeross(tmp2);
		end

		x = cat(1, tmp1, x, tmp2);
	end

	if spad(2) || epad(2)
		pad = spad(2);
		tmp1 = x(:,(end-pad+1):end,:);
		pad = epad(2);
		tmp2 = x(:,1:pad,:);
		if ~per(2)
			tmp1 = zeross(tmp1);
			tmp2 = zeross(tmp2);
		end
		x = cat(2, tmp1, x, tmp2);
	end

	if spad(3) || epad(3)
		pad = spad(3);
		tmp1 = x(:,:,(end-pad+1):end);
		pad = epad(3);
		tmp2 = x(:,:,1:pad);
		if ~per(3)
			tmp1 = zeross(tmp1);
			tmp2 = zeross(tmp2);
		end
		x = cat(3, tmp1, x, tmp2);
	end

otherwise
	fail('todo: ndim=%d', ndim)
end


% y = ir_conv_per()
function y = ir_conv_per(x, psf, per)

x = ir_conv_pad(x, psf, per);
y = convn(x, psf, 'valid');


% ir_conv_test()
function ir_conv_test

psf = [1 -2]';
x = (1:6)';
y = ir_conv(x, psf, 'per', true);
z = ifft(fft(x) .* fft(psf, numel(x)));
equivs(z, y)

mask = true(6, 4); psf = [1 5 3]';
mask(7) = false; % stress
A = Gblur(mask, 'psf', psf, 'type', 'conv,per');
B = Gblur(mask, 'psf', psf, 'type', 'fft,same');

if 1
	im plc 2 3
	im(1, full(A)')
	im(2, full(A'))
	im(3, full(A') - full(A)')
	im(4, full(B)')
	im(5, full(B'))
	im(6, full(B') - full(B)')
end
equivs(full(A), full(B))
equivs(full(A'), full(B'))

test_adjoint(A, 'complex', 1);
