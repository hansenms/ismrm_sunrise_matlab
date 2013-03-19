% Gtomo3_test_adj.m
% Test adjoint of Gtomo3 object(s)

if ~isvar('G'), printm 'setup Gtomo3_test_adj'
	nx = 8;
	ny = 8;
	nz = 10;

	if 0 % 3l
	%	3l@dis_src_obj,dis_obj_det,nu,nv,su,sv,sz,cx,cy,cz,cu,cv@-@-2d,nphi@nder2.file-
		f.sys_type = '3l@200,60,20,10,2,2,1,0,0,0,0,0@-@-2d,5@-';
	else % 3s
		if 1
			f.fwhm_collimator = 1;
		 	f.fwhm_iso = 2; % depth-dependent gaussian blur
			f.psfs = '-';
			f.blur = sprintf(',gauss,%g,%g', ...
				f.fwhm_collimator, f.fwhm_iso);
		elseif 0
			f.psfs = '-';
		else
			f.psfs = '/tmp/t,psfs.fld';
			psfs = make_3s_psfs(ny, 1, 1.2*nx, 0, 2/nx);
			fld_write(f.psfs, psfs)
		end
		f.mumap = '-';
		f.sfilter = 1;
		f.sys_type = sprintf('3s@1,1,1,360,0,%d%s@%s@%s@-%d,%d,%d', ...
			f.sfilter, f.blur, f.mumap, f.psfs, nx, nz, 5)
	%	3s@$sx,$sx,$sx,$orbit,$ostart,$sfilter@$mumap@$psfs@-$nx,$nz,$na
	end

	mask = true([nx ny nz]);
	G = Gtomo3(f.sys_type, mask, nx, ny, nz, ...
		'nthread', 1, 'chat', 0, 'checkmask', 0);
end

Fatrix_test_basic(G, mask)

[A Aa] = test_adjoint(G);
im clf, im pl 1 3, im(1, A, 'A'), im(1, Aa', 'Aa''')

x = reshape(G' * ones(nrow(G),1), [nx ny nz]);
im(1, x, 'G''1'), cbar horiz

% im clf, plot(squeeze(sum(sum(x,1),2))/nx/ny)
