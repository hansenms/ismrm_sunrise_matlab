function sinit = ir_mri_sensemap_init(init, ykj, bodycoil, sizeI, thresh, maskObj)

nx = sizeI(1);
ny = sizeI(2);
ncoil = sizeI(3);

if ~isempty(maskObj)
    warning('using the maskObj instead of thresholding to determine good pixels');
end

if isempty(init) || ischar(init)
        sinit = zeros(nx, ny, ncoil);
	for ic = 1:ncoil
		zj = ykj(:,:,ic);
		tmp = zj ./ bodycoil; % usual ratio
                
                % determine "good" pixels
                if isempty(maskObj)
                    good = abs(bodycoil) > thresh * ...
                           max(abs(bodycoil(:)));
                else
                    good = maskObj;
                end
                
                if isempty(init) || strcmp(init, 'ratio')
                    % set all uncertain map values to median of good ones
                    disp('Using zero background.');
                    tmp(~good) = 0; % median(abs(tmp(good)));; 
                    tmp = reshape(tmp, [nx,ny]);
                elseif isempty(init) || strcmp(init, 'median')
                    % set all uncertain map values to median of good ones
                    disp('Using median background.');
                    tmp(~good) = median(tmp(good));
                    tmp = reshape(tmp, [nx,ny]);
                elseif isempty(init) || strcmp(init, 'avg')
                    % set all uncertain map values to median of good ones
                    disp('Using mean background.');
                    tmp(~good) = mean(abs(tmp(good))).*exp(i*mean(angle(tmp(good))));
                    tmp = reshape(tmp, [nx,ny]);
                elseif strcmp(init, 'zeros')
                    tmp = zeros([nx ny]);
                elseif strcmp(init, 'order1')
                    disp('Using 1st order fit for background.');
                    expo = [0 0; 1 0; 0 1];
                    tmp = ortho_init(nx,ny,expo,bodycoil,good);
                elseif strcmp(init, 'order2')
                    disp('Using 2nd order fit for background.');
                    expo = [0 0; 1 0; 0 1; 1 1; 2 0; 0 2];
                    xx = ndgrid_jf('mat', 0:nx-1, 0:ny-1);
                    tmp = ortho_init(nx,ny,expo,bodycoil,good);
                elseif strcmp(init,'order3')
                    disp('Using 3rd order fit for background.');
                    expo = [0 0; 1 0; 0 1; 1 1; 2 0; 2 1; 0 2; 1 2; ...
                           3 0; 0 3];
                    tmp = ortho_init(nx,ny,expo,bodycoil,good);
                elseif strcmp(init,'order4')
                    disp('Using 4th order fit for background.');
                    expo = [0 0; 1 0; 0 1; 1 1; 2 0; 2 1; 0 2; 1 2; ...
                           3 0; 0 3; 2 2; 3 1; 1 3; 4 0; 0 4];
                    tmp = ortho_init(nx,ny,expo,bodycoil,good);
                else 
                    error('Do not recognize initialization type.')
                end
		sinit(:,:,ic) = tmp;
	end
else
    sinit = init;
end

end

function init = ortho_init(nx,ny,expo,bodycoil,good)
xx = ndgrid_jf('mat', 0:nx-1, 0:ny-1);
A = cheby2d(nx, ny, expo);
A = reshape(A, [nx*ny, length(expo)]);
maskBodyI = good .* bodycoil;
Aw = repmat(maskBodyI(:),[1 length(expo)]) .* A;
theta = pinv(Aw) * zj(:); clear Aw;
oInit = A * theta; %not Aw
oInit = reshape(oInit, [nx ny]);
end

function polys = cheby2d(nx,ny,expo)
nbases = length(expo);
x = linspace(-1,1,nx);
y = linspace(-1,1,ny);
polys = zeros(nx,ny,nbases);
for n = 1:nbases
    tmp = polyval(ChebyshevPoly(expo(n,1)),x);
    tmp = tmp';
    X = repmat(tmp,[1 ny]);
    tmp = polyval(ChebyshevPoly(expo(n,2)),y);
    Y = repmat(tmp,[nx 1]);
    polys(:,:,n) = X .* Y;
end
end
