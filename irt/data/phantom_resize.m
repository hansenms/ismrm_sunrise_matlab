 function x = phantom_resize(x, nx, ny)
%function x = phantom_resize(x, nx, ny)
% resize a [mx,my] phantom image to be [nx,ny]
% by combination of downsampling and cropping

mx = size(x,1);
my = size(x,2);
if nx == mx && ny == my
	return
end

if nx > mx || ny > my
	error 'requested size too large, implement interpolation!'
end

% if power of 2 size, then downsample
if 2^floor(log2(nx)) == nx && nx == ny
	x = downsample2(x, mx/nx);
return
end

nn = 2^max(ceil(log2([nx ny])));
x = downsample2(x, mx/nn);

% trim if non power of 2
if nx < nn
	x = x([1:nx]+round((nn-nx)/2),:);
end
if ny < nn
	x = x(:,[1:ny]+round((nn-ny)/2));
end
