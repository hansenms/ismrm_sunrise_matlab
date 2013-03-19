  function xpad = pad_into_center(x, npad, varargin)
%|function xpad = pad_into_center(x, npad, varargin)
%|
%| Zero pad an input signal x symmetrically around "0" (image center).
%| Useful for DFT/FFT of PSF.
%| using this avoids matlab's padarray.m which is in image proc toolbox
%|
%| option
%|	'circ'	0|1	if 1, circular padding instead of 0.  (default 0)
%|
%| Originally by A. Yendiki, modified by Jeff Fessler, 2005-7-26.

if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(x, 'test'), pad_into_center_test, return, end

arg.circ = false;
arg = vararg_pair(arg, varargin);

if length(npad) == 1
	if min(size(x)) == 1 % 1d.  kludge: wrong if size(x) = [nx 1 nz]
		if size(x,1) == 1
			npad = [1 npad];
		else
			npad = [npad 1];
		end
	else % n-dimensional; pad all dimensions the same amount
		npad = repmat(npad, [1 ndims(x)]);
	end
end


ndim = ndims(x);
if ndim ~= length(npad)
	error 'Incorrect number of dimensions'
end

if any(size(x) > npad)
	pr size(x)
	pr npad
	fail('npad too small')
end
shift = ceil((npad - size(x))/2);

if arg.circ % circular boundary condition padding
	if any(round(shift) ~= shift)
		fail 'circular needs same padding on both sides'
	end
	args = cell(ndim,1);
	for id = 1:ndim
		nold = size(x,id);
		nnew = npad(id);
		args{id} = 1 + mod([0:nnew-1]-shift(id), size(x,id));
	end
	args = ndgrid_jf('cell', args{:});
	index = sub2ind(size(x), args{:});
	xpad = x(index);

else
	args = cell(ndim,1);
	for id = 1:ndim
		nold = size(x,id);
		nnew = npad(id);
		if nold > nnew
			fail('Padding[%d]=%d too small cf %d', id, nnew, nold)
		end
		args{id} = [1:nold] + ceil((nnew - nold)/2);
	end
	if islogical(x)
		xpad = false(npad);
	else
		xpad = zeros(npad, class(x));
	end
	xpad(args{:}) = x;
end


function pad_into_center_test
pad_into_center([1 2 1], [7])
pad_into_center([1 2 1]', [7])'
pad_into_center([1 2 1], [7 3])
pad_into_center([1 2 1]', [7 3])
pad_into_center(ones(3), [5 7])

x = reshape(1:3*4, 3, 4)
y = pad_into_center(x, [3+4 4+2], 'circ', 1)
if exist('padarray', 'file')
	z = padarray(x, [2 1], 'circular', 'both')
	jf_equal(y, z)
end
