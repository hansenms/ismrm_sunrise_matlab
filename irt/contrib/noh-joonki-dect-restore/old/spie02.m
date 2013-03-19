% xray simulation for spie '02
% 2001-04-25 started
% 2002-01-27 continued with new more linear tables
% 2004-05-04 updated to use new de_* routines
% Copyright 2001-04-25, Jeff Fessler, The University of Michigan

savefig off % disable figure saving

if ~isvar('ftab'), printm 'table'
	f.dir = '/home/fessler/l/dat/sim/spie,02/'; % output
%	f.stype = 'mono,60,100';
	f.stype = 'ps1';
	ftab = de_ftab_build({f.stype});
%	f.table = ['table/' f.stype];
	f.fig = ['../fig,' f.stype  ',128,new/'];
	Nrms = inline('100*nrms(col(x(:,:,i)), col(y(:,:,i)))', 'x', 'y', 'i');
prompt
end

%
% true density map (soft tissue and bone)
%
if ~isvar('xtrue'), printm 'xtrue'
	f.isdirect = logical(0);
	if 0	% simple squares
		n.x = 30; n.y = 32;
		n.b = n.x; n.a = n.y;
		f.dx = 0.45 * 128 / n.x;
		xtrue = zeros([n.x n.y 2]);
		xtrue(5:26,5:28,1) = 1;
		xtrue(11:20,11:22,1) = eps;
		xtrue(16:20,11:22,2) = 2;
		f.isdirect = logical(1);

	elseif 0	% ellipses
		f.tmp = [f.dir 't0'];
		n.x = 128; n.y = 104;	n.b = 140; n.a = 128;
%		n.x = 32; n.y = 26;	n.b = 35; n.a = 32;
		f.dx = 0.45 * 128 / n.x;	dx = f.dx;
		com = sprintf('op ellipse %s %d %d', f.tmp, n.x, n.y) + ...
			sprintf(' 0 0 %g %g 0 1 3', 45/dx/2, 35/dx/2) + ...
			sprintf(' %g 0 %g %g 0 -0.7 3', -11/dx, 12/dx/2, 16/dx/2) + ...
			sprintf(' %g 0 %g %g 0 -0.7 3', +11/dx, 12/dx/2, 16/dx/2);
		os_run(com)
		xtrue(:,:,1) = double(fld_read(f.tmp));
		com = sprintf('op ellipse %s %d %d', f.tmp, n.x, n.y) + ...
		sprintf(' 0 %g %g %g 0 2 3', -10/dx, 4/dx/2, 4/dx/2);
		os_run(com)
		xtrue(:,:,2) = double(fld_read(f.tmp));
		xtrue(:,:,1) = xtrue(:,:,1) - xtrue(:,:,2) / 2;
		clear dx com

	else	% ncat

	    f.xdir = strrep(path_find_dir([filesep 'ct']), 'ct', 'data');
		f.xfile = [f.xdir '/ncat,256,slice,140,ct,x100.fld'];
		x = fld_read(f.xfile) / 100;
		x = double(x(:,20:end-29));
		x = x([end 1:end-1],:);		% shift over 1 pixels
		f.dx = 0.16 * 2;
		x1 = x .* (x < 1.5);	% soft tissue
		x2 = x .* (x > 1.5);	% bone
        %x2 = (x2 > 1.5)*0.1;    % for only PET
        x1 = ( x1(1:2:end,:) + x1(2:2:end,:) ) / 2;	% downsample
		x1 = ( x1(:,1:2:end,:) + x1(:,2:2:end) ) / 2;
		x2 = ( x2(1:2:end,:) + x2(2:2:end,:) ) / 2;
		x2 = ( x2(:,1:2:end,:) + x2(:,2:2:end) ) / 2;
		xtrue(:,:,1) = x1; %[128 104 1]
		xtrue(:,:,2) = x2; %[128 104 1]
		[n.x n.y] = size(x1);
		n.b = 140; n.a = 128;
		clear x1 x2 x
	end

	%
	% mask based on dilating true object
	%
	mask = sum(xtrue, 3) > 0;
	if n.x > 100
%		mask = erode(mask, ones(3), 10);
%		mask = imerode(mask, strel('disk', 4));
%		mask = dilate(mask, ones(3), 15);
		mask = imdilate(mask, strel('disk', 5));
	else
		mask = dilate(mask, ones(3), 1);
	end
	mask = logical(mask);

	im clf
	c.soft = [0 1.2];
	c.bone = [0 2.2];
	c.dens = [0 2.2];
	im(221, xtrue(:,:,1), 'SoftTissue Density', c.soft), cbar([0 1])
	im(222, xtrue(:,:,2), 'Bone Density', c.bone), cbar([0 2])
	im(223, sum(xtrue,3), 'Density Map', c.dens), cbar([0 2])
%	xlabel x, ylabel y
	im(224, (sum(xtrue,3)>0.5) + 20*double(mask), 'Reconstruction Support')
	cbar hide
%	subplot(235), plot(sum(mask,2), 'o'), axis tight

	savefig(f.fig, 'fig_object')
prompt
	clear G ymi
end


%
% system matrix
%
if ~isvar('G'), printm 'G'
	ig = image_geom('nx', n.x, 'ny', n.y, 'dx', f.dx, 'mask', mask);
	sg = sino_geom('par', 'nb', n.b, 'na', n.a, 'dr', f.dx);
	G = Gtomo2_strip(sg, ig);
    %G = Gtomo2_table(sg, ig);
    %G = Gtomo2_wtmex(sg, ig);

if 0 % old G - delete?
	f.dsc = [f.dir 't.dsc'];
	f.wtf = strrep(f.dsc, '.dsc', '.wtf');
	f.wtr = strrep(f.wtf, '.wtf', ',row.wtf');
	f.mask = [f.dir 'mask.fld'];
	mat_write(f.mask, mask, '-nocheck', '-type', 'byte');

	if f.isdirect
		% direct observation case for fast debugging
		n.b = n.x;	n.a = n.y;
		delete(f.wtf)
		wtfmex(f.wtf, speye(n.x*n.y), n.x, n.y, n.b, n.a)
	else
		f.dim = sprintf('nx %d ny %d nb %d na %d', n.x, n.y, n.b, n.a);
		f.size = sprintf('pixel_size %g', f.dx);
		com = 'wt -chat 0 dsc -support "file %s" 2 %s scale 0 %s >! %s';
		com = sprintf(com, f.mask, f.dim, f.size, f.dsc);
		os_run(com)
		os_run(sprintf('echo y | wt gen %s', f.dsc))
		os_run(sprintf('echo y | wt col2row %s %s', f.wtr, f.wtf))
	end
	G = Gtomo2_wtfmex(f.wtr, 0);
%	G = Gtomo2_sparse(f.wtf, 0);

	if any(mask ~= G.mask), error mask, end
	G = G(:,mask);	% compact form

end
	im(mask + xtrue(:,:,1), 'mask + xtrue')

	gi = reshape(sum(G'), [n.b n.a]);
	im clf, im(121, gi, 'gi'), im(122, gi>0, 'gi>0')
prompt
end


%
% ideal and noisy projections
%
if ~isvar('ymi'), printm 'ymi'
	tmp = reshape(xtrue, n.x*n.y, 2);
	tic
	strue = reshape(G * tmp(mask,:), [n.b n.a 2]);
	printm('forward projection time %g', toc)

	im clf
	im(421, strue(:,:,1), 's1'), cbar
	im(422, strue(:,:,2), 's2'), cbar

%	load(f.table)

	s1max = max(col(strue(:,:,1)));
	s2max = max(col(strue(:,:,2)));
	printm('type1 headroom: %g%%', (max(ftab.sl{1})/s1max-1)*100)
	printm('type2 headroom: %g%%', (max(ftab.sl{2})/s2max-1)*100)
	if s1max > max(ftab.sl{1}) | s2max > max(ftab.sl{2})
		error 'table too small'
	end

	tic
%	ftrue = ftab.feval(ftab, strue(:,:,1), strue(:,:,2));
	ftrue = ftab.fm_fun(ftab, {strue(:,:,1), strue(:,:,2)});
	if any(isnan(ftrue(:)) | isinf(ftrue(:))), error nan-inf, end
%	printm('interpolation time %g', toc)
	im(423, ftrue(:,:,1), 'f1'), cbar
	im(424, ftrue(:,:,2), 'f2'), cbar
	ybi = zeros(n.b, n.a, 2);
	ybi(:,:,1) = ftab.xray.I(1) * exp(-ftrue(:,:,1));
	ybi(:,:,2) = ftab.xray.I(2) * exp(-ftrue(:,:,2));
	im(425, ybi(:,:,1), 'ybar1'), cbar
	im(426, ybi(:,:,2), 'ybar2'), cbar

	if f.isdirect | 0
		f.scale = inf;
	else
		f.scale = 10^6 / ybi(1,1,2);
	end
	printm('f.scale = %g', f.scale)
	if isinf(f.scale)
		ymi = ybi;	% noiseless
	else
		ymi = poisson(f.scale*ybi, 0) / f.scale;
	end

	f.title1 = sprintf('%2.0f kVp', ftab.xray.kvp(1));
	f.title2 = sprintf('%2.0f kVp', ftab.xray.kvp(2));
	if isinf(f.scale)
		f.title1 = [f.title1 ' (Noiseless)'];
		f.title2 = [f.title2 ' (Noiseless)'];
	else
		printm('total counts %g', sum(ymi(:)) * f.scale)
	end

	im clf, printm 'pretty figure of ymi'
	t.y1 = ymi(:,:,1) * f.scale;	% show counts!
	t.y2 = ymi(:,:,2) * f.scale;
	c.ymi1 = [0 floor(max(t.y1(:))/1e3)*1e3];
	c.ymi2 = [0 floor(max(t.y2(:))/1e4)*1e4];
	im(221, t.y1, f.title1, c.ymi1), cbar(c.ymi1)
	im(222, t.y2, f.title2, c.ymi2), cbar(c.ymi2)
	savefig(f.fig, 'fig_ymi')

	clear fhat s_water tmp xpl Rml Rwls denom wi s1max s2max t
prompt
end


%
% filter log sinograms to prepare for (FBP) reconstruction
%
if ~isvar('fhat'), printm 'fhat'
	if isinf(f.scale) | 1
		f.kernel = [1]';
	else
		f.kernel = [1 8 1]';
	end
	f.kernel = f.kernel / sum(f.kernel);	% filter for fbp
	ymi_filt = convn(ymi, f.kernel, 'same');

	% borrow neighbor(s) for any log(0) values
	off = 1;
	while any(ymi_filt(:) == 0)
		ii = find(ymi_filt(:) == 0);
		printm('fixing %d zeros in ymi_filt with %d', length(ii), off)
		ymi_filt(ii) = max(ymi_filt(ii+off), ymi_filt(ii-off));
		off = off + 1;
	end, clear off

%	fhat.fil(:,:,1) = -log(ymi_filt(:,:,1) / ftab.xray.I(1));
%	fhat.fil(:,:,2) = -log(ymi_filt(:,:,2) / ftab.xray.I(2));
	fhat.raw(:,:,1) = -log(ymi(:,:,1) / ftab.xray.I(1));
	fhat.raw(:,:,2) = -log(ymi(:,:,2) / ftab.xray.I(2));
	fhat.raw(isinf(fhat.raw)) = 0;

	im clf
	c.fhat = [0 ceil(max(fhat.raw(:)))];
	im(221, fhat.raw(:,:,1), f.title1, c.fhat), cbar
	im(222, fhat.raw(:,:,2), f.title2, c.fhat), cbar
	if 0
		if length(f.kernel) == 1 | 1
			tmp = 'Log-processed';
		else
			tmp = 'Smoothed and Log-processed';
		end
		text(-0.3*n.x, 1.7 * n.y, tmp, 'horizontalalign', 'center')
		text(-0.3*n.x, 1.9 * n.y, 'Sinogram Measurements', ...
			'horizontalalign', 'center')
	end
	savefig(f.fig, 'fig_fhat')
prompt
	clear ymi_filt shat
end

if ~isvar('fhat.err'), printm 'pretty figure of fhat errors'
	im clf
	if im
		plot(fhat.raw(:,:,1)-ftrue(:,:,1), fhat.raw(:,:,2)-ftrue(:,:,2), 'y.')
		axis([-1 1 -0.5 0.5])
%		plot(fhat.fil(:,:,1)-ftrue(:,:,1), fhat.fil(:,:,2)-ftrue(:,:,2), 'y.')
		ytick, xlabel 'f_1 error', ylabel 'f_2 error'
	end

%	fhat.err.fil(1) = max_percent_diff(fhat.fil(:,:,1), ftrue(:,:,1));
%	fhat.err.fil(2) = max_percent_diff(fhat.fil(:,:,2), ftrue(:,:,2));
%	printm('fhat.fil err %g%% %g%%', fhat.err.fil(1), fhat.err.fil(2))
	fhat.err.raw(1) = max_percent_diff(fhat.raw(:,:,1), ftrue(:,:,1));
	fhat.err.raw(2) = max_percent_diff(fhat.raw(:,:,2), ftrue(:,:,2));
	printm('fhat.raw err %g%% %g%%', fhat.err.raw(1), fhat.err.raw(2))
	savefig(f.fig, 'fig_fhat_err')
prompt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% todo: needs updating after here %
% could use de_ftab_s_iter        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% water correction only via 1d table
%
if ~isvar('s_water'), printm 's_water'
	s_water(:,:,1) = ftab.inv_water_eval(ftab, 1, fhat.raw(:,:,1));
	s_water(:,:,2) = ftab.inv_water_eval(ftab, 2, fhat.raw(:,:,2));
	if any(isnan(s_water(:)) | isinf(s_water(:))), error 'nan-inf', end
	im(221, s_water(:,:,1), 's water 1'), cbar
	im(222, s_water(:,:,2), 's water 2'), cbar
	subplot(223), plot([s_water(:,end/2,1) strue(:,end/2,1)])
	subplot(224), plot([s_water(:,end/2,2) strue(:,end/2,1)])
prompt
end


%
% estimate component ray integrals
%
if ~isvar('shat'), printm 'shat'

%	shat.raw = zeros(size(fhat.raw));
	shat.raw = ftab.inv.eval(ftab, fhat.raw);
%	shat.fil = ftab.inv.eval(ftab, fhat.fil);

	if 0
		err = abs(shat.fil - strue);
		err = err > 3;
		err(:,:,1) = err(:,:,1) | err(:,:,2);
		err(:,:,2) = err(:,:,1) | err(:,:,2);
		shat.fil = ftab.inv.eval(ftab, fhat.fil .* err);
	%	shat.fil = ftab.inv.eval(ftab, ftrue);
	return
	end

	% ad hoc constraints on s1 and s2
	shat.raw(:,:,1) = max(shat.raw(:,:,1), -1);
	shat.raw(:,:,1) = min(shat.raw(:,:,1), 60);
	shat.raw(:,:,2) = max(shat.raw(:,:,2), -1);
	shat.raw(:,:,2) = min(shat.raw(:,:,2), 20);
%	shat.fil(:,:,1) = max(shat.fil(:,:,1), -1);
%	shat.fil(:,:,1) = min(shat.fil(:,:,1), 60);
%	shat.fil(:,:,2) = max(shat.fil(:,:,2), -1);
%	shat.fil(:,:,2) = min(shat.fil(:,:,2), 20);

	err = shat.raw - strue;
	printm('shat.raw max error %g', max(abs(err(:))))
	err = abs(err);

	if f.isdirect
		c.s1 = [0 1.2];
		c.s2 = [0 2.4];
		c.s_err = [0 0.1];
	else
		c.s1 = [0 30];
		c.s2 = [0 20];
		if isinf(f.scale)
			c.s_err = [0 1];
		else
			c.s_err = [0 10];
		end
	end

	im(231, strue(:,:,1), 'True Rays', c.s1), cbar(c.s1)
	ylabel 'Soft Tissue'
	im(234, strue(:,:,2), ' ', c.s2), cbar(c.s2)
	ylabel 'Bone'
	im(232, shat.raw(:,:,1), 'Estimates', c.s1), cbar(c.s1)
	im(235, shat.raw(:,:,2), ' ', c.s2), cbar(c.s2)
	im(233, err(:,:,1), '|Errors|', c.s_err), cbar(c.s_err)
	im(236, err(:,:,2), ' ', c.s_err), cbar(c.s_err)
	savefig(f.fig, 'fig_shat')

prompt
	clear xfbp
end


if ~isvar('shat.err'), printm 'pretty figure of shat errors'
	im clf
	if im
		plot(shat.raw(:,:,1)-strue(:,:,1), shat.raw(:,:,2)-strue(:,:,2), 'y.')
		xlabel 's_1 error', ylabel 's_2 error'
		axis([-20 20 -15 15]), xtick, ytick
		title 'Component line integral errors'
		savefig(f.fig, 'fig_shat_err')
	end

	shat.err.raw(1) = max_percent_diff(shat.raw(:,:,1), strue(:,:,1));
	shat.err.raw(2) = max_percent_diff(shat.raw(:,:,2), strue(:,:,2));
	shat.err.raw(1) = max_percent_diff(shat.raw(:,:,1), strue(:,:,1));
	shat.err.raw(2) = max_percent_diff(shat.raw(:,:,2), strue(:,:,2));
%	printm('sfil err %g%% %g%%', shat.err.fil(1), shat.err.fil(2))
	printm('sraw err %g%% %g%%', shat.err.raw(1), shat.err.raw(2))
prompt
end


%
% FBP reconstruction from processed projections
%
if ~isvar('xfbp'), printm 'xfbp'
	if f.isdirect
		xfbp = shat.raw;
		fbp_title = 'FBP';
	else
		tmp = fbp2(sg, ig);
		fbp_title = ['FBP [1 2 1]/4'];
		f.fbp_kernel = [1 2 1]'; % fixed on 2006-3-3, before it had no '
		f.fbp_kernel = f.fbp_kernel / sum(f.fbp_kernel);
		sino = conv2(shat.raw(:,:,1), f.fbp_kernel, 'same');
		% im([shat.raw(:,:,1); sino])
		xfbp(:,:,1) = fbp2(sino, tmp);
		sino = conv2(shat.raw(:,:,2), f.fbp_kernel, 'same');
		xfbp(:,:,2) = fbp2(sino, tmp);

%		arg = arg_pair('system', 9, 'nx', n.x, 'ny', n.y, ...
%			'nb', n.b, 'na', n.a, 'pixel_size', f.dx, ...
%			'scale', 0, 'support', 'all');
%%			'scale', 0, 'support', ['file ' f.mask]);
%		oxfbp(:,:,1) = fbp_dsc(shat.raw(:,:,1), f.fbp_kernel, arg);
%		oxfbp(:,:,2) = fbp_dsc(shat.raw(:,:,2), f.fbp_kernel, arg);
%		oxfbp = max(oxfbp, 0);
	end
	xfbp = max(xfbp, 0);

	if f.isdirect
		c.err = [0 0.1];
	else
		c.err = [0 1];
	end
	printm('FBP %g %g', Nrms(xfbp, xtrue, 1), Nrms(xfbp, xtrue, 2))
	plot_de_true_recon_err(xtrue, xfbp, fbp_title, c)
	savefig(f.fig, 'fig_true_fbp_err')
prompt
end


%
% examine water correction version
%
if ~isvar('xwater'), printm 'xwater'
	if f.isdirect
		xwater = s_water;
	else
%		xwater(:,:,1) = fbp_dsc(s_water(:,:,1), f.fbp_kernel, arg);
%		xwater(:,:,2) = fbp_dsc(s_water(:,:,2), f.fbp_kernel, arg);
		xwater(:,:,1) = fbp2(s_water(:,:,1), tmp);
		xwater(:,:,2) = fbp2(s_water(:,:,2), tmp);
	end
	xwater = max(xwater, 0);
	dens = sum(xtrue, 3);

	im clf
	im(231, dens, 'True Density', c.bone)
	set(gca, 'position', get(gca, 'position') + [-0.01 -0.25 0 0])
	cbar([0 2])
	text(-0.1*n.x, 2.3*n.y, 'Water')
	text(-0.1*n.x, 2.6*n.y, 'Correction')
	im(232, xwater(:,:,1), f.title1, c.bone), cbar([0 2])
	im(235, xwater(:,:,2), f.title2, c.bone), cbar([0 2])
	im(233, abs(xwater(:,:,1)-dens), '|Error|', c.err), cbar(c.err)
	im(236, abs(xwater(:,:,2)-dens), '|Error|', c.err), cbar(c.err)
	savefig(f.fig, 'fig_water')
prompt
end


% play with color maps
if 0
	err = xfbp - xtrue;
	rgb = uint8(zeros([n.y n.x 3 2]));
	t = err(:,:,1)';
	t = max(t,-1); t = min(t,1);
	blue = (t > 0) .* (255 * t);
	green = (t < 0) .* (-255 * t);
	rgb(:,:,2,1) = green;
	rgb(:,:,3,1) = blue;
	subplot(233)
	h = image(rgb(:,:,:,1)); axis image
%	get(h)
	cbar
prompt
end


%
% build block system object
%
if ~isvar('Gb'), printm 'Gb'
	if f.isdirect
		f.nsubset = 1;
	else
		f.nsubset = n.a/8;
    end
    %clear all;load de_ct_setup_after_FBP.mat;
	Gb = Gblock(G, f.nsubset, 0);

	clear Rml Rwls denom xpwls xpl
prompt
end


%
% see local resolution vs beta
%
f.l2b = 5;	% too blurry at 128
f.l2b = 1;	% pretty good for pwls
f.l2b = 2;
if f.isdirect	% for direct observation case
	f.l2b = -2;
	f.l2b = -inf;
end

if 0, printm 'checking resolution'
	C0 = C2sparse('leak', ig.mask, 8, 0);
	C0 = C0(:,mask(:));
	psf = qpwls_psf(G, C0, 2^f.l2b, ig.mask);
	im clf, subplot(221), fwhm2(psf);
	im(222, psf, 'psf'), cbar

%	tmp = convn(xtrue, psf, 'same');
	tmp = reale(ifft2(fft2(sum(xtrue,3)) .* fft2(fftshift(psf))), 'warn');
	im(223, sum(xtrue,3), 'true')
	im(224, tmp, 'blurred')
	clear C0 tmp
prompt
return
end

% potential function
if ~isvar('f.pot')
	f.pot = 'huber';
	clear Rwls Rml
end


%
% precomputed 'curvatures' and denominator for PWLS
%
if ~isvar('denom.wls'), printm 'denom.wls'
	wi.wls = ymi;
%	denom.wls(:,1) = G' * col(wi.wls(:,:,1) .* gi);
%	denom.wls(:,2) = ;
	denom.wls = zeros(0,2);
%	im clf, im(embed(denom.wls, mask), 'WLS denom'), cbar
prompt
end


%
% penalty for PWLS recon
%
if ~isvar('Rwls'), printm 'Rwls'
	kappa.wls = G' * [col(wi.wls(:,:,1)) col(wi.wls(:,:,2))];
	kappa.wls = sqrt( kappa.wls ./ (G' * ones(n.b*n.a,2)) );
	kappa.wls = embed(kappa.wls, mask);
	im clf, im(211, kappa.wls, 'kappa WLS')
	f.delta_wls = 0.2;	%
	f.delta_wls = 0.05;	% nrms 13.8339 17.6569
	f.delta_wls = 0.1;	% nrms 12.9424 17.5817 (8 iter)
				% nrms 13.601 16.3552 (16 iter)
%	f.nbrs = 8;

	for ii=1:2
%		Rwls{ii} = Rbuild('leak', kappa.wls(:,:,ii), f.nbrs, ...
%			2^f.l2b, f.pot, f.delta_wls, 0);
		Rwls{ii} = Robject(kappa.wls(:,:,ii), ...
			'type_denom', 'matlab', ...
		'beta', 2^f.l2b, 'potential', f.pot, 'delta', f.delta_wls);
	end
	im(212, reshape(Rwls{1}.wt, [n.x n.y 4]), 'wt WLS'), cbar
prompt
	clear C xpwls
end


%
% PWLS reconstruction (FBP initializes the algorithm)
% For now this is individual recon of the two material images.
% It should be generalized to allow for 2x2 covariance matrices...
%
if ~isvar('xpwls'), printm 'xpwls'
%	xinit = xtrue;			% true init
	xinit = xfbp;			% fbp init

	xinit = max(xinit,0);
	im(xinit, 'init')

	if 1, printm 'enforcing object constraints'
		t1 = xinit(:,:,1);
		t2 = xinit(:,:,2);
		t1 = medfilt2(t1, [5 5]); %2D median filtering
		t2 = medfilt2(t2, [3 3]);
		t1 = min(t1, 1.1);
		t2 = min(t2, 2.1);
		t2(t2 < 0.4) = 0;
		xinit(:,:,1) = t1;
		xinit(:,:,2) = t2;
	end, clear t1 t2

	xinit = reshape(xinit, n.x*n.y, 2);	% prepare to iterate
	xinit = xinit(mask,:);

	f.niter = 8+1;
	f.niter = 16+1;
	f.pixmax = 99;

	xpwls = zeros([n.x n.y 2 f.niter]);

	for ll=1:2 %main iteration for pwls
%	profile on
		tmp = pwls_sps_os(xinit(:,ll), ...
			shat.raw(:,:,ll), wi.wls(:,:,ll), Gb, ...
			Rwls{ll}, f.niter, f.pixmax, denom.wls(:,ll), gi);
			xpwls(:,:,ll,:) = ig.embed(tmp);
%	profile report
	end

	im clf, im(xpwls, 'xpwls'), cbar
prompt
end

if 0, printm 'pwls cost'
	cost1 = pwls_cost(squeeze(xpwls(:,:,1,:)), G, spdiag(wi.wls(:,:,1)), ...
		col(strue(:,:,1)), Rwls{1}, ig.mask);
	im clf, plot(0:f.niter-1, cost1, '-o')
prompt
end

if 0
	plot_de_true_recon_err(xtrue, xpwls(:,:,:,end), 'PWLS', c)
	savefig(f.fig, 'fig_true_pwls_err')
prompt
end

%
% pwls vs fbp figure
%
if 0
	figure(1), clf
	im(231, xtrue(:,:,1), 'True', c.soft), cbar([0 1])
	ylabel 'Soft Tissue'
	im(234, xtrue(:,:,2), ' ', c.bone), cbar([0 2])
	ylabel 'Cortical Bone'
	im(232, xfbp(:,:,1), 'FBP', c.soft), cbar([0 1])
	im(235, xfbp(:,:,2), ' ', c.bone), cbar([0 2])
	im(233, xpwls(:,:,1,end), 'PWLS', c.soft), cbar([0 1])
	im(236, xpwls(:,:,2,end), ' ', c.bone), cbar([0 2])
	savefig(f.fig, 'fig_true_fbp_pwls')

	printm('FBP %g %g', Nrms(xfbp, xtrue, 1), Nrms(xfbp, xtrue, 2))
	printm('PWLS %g %g', ...
	Nrms(xpwls(:,:,:,end), xtrue, 1), Nrms(xpwls(:,:,:,end), xtrue, 2))
return
%prompt
end


%
% reproject at 511 keV to form the line integrals that would be used
% for attenuation correction of a PET scan.
% these are among the type of results that could be used for a MIC04 abstract.
%
if 0
	figure(2), clf
	% mass atten. coef at 511
	m_soft = xray_read_atten('soft', 511);
	m_bone = xray_read_atten('bone', 511);

	% true attenuation map projections
	if ~isvar('p_true')
		p1 = G * xtrue(:,:,1);
		p2 = G * xtrue(:,:,2);
		p_true = m_soft * p1 + m_bone * p2;
		im(p_true)
	end

	% estimate attenuation map projections using DE
	if ~isvar('p_de_fbp')
		p1 = G * xfbp(:,:,1);
		p2 = G * xfbp(:,:,2);
		p_de_fbp = m_soft * p1 + m_bone * p2;
		im(p_de_fbp)
	end

	if ~isvar('p_de_pwls')
		p1 = G * xpwls(:,:,1);
		p2 = G * xpwls(:,:,2);
		p_de_pwls = m_soft * p1 + m_bone * p2;
		im(p_de_pwls)
	end

	% estimate attenuation map projections using water correction only
	% and piecewise linear density to attenuation approximation.
	if ~isvar('p_water')
		x = xwater(:,:,2); % use higher energy CT (closer to 511?)
		f_soft = max(min(1.0, 2.0 - x), 0); % soft tissue fraction
		x = (m_soft * f_soft + m_bone * (1-f_soft)) .* x;
		p_water = G * x;
		im(p_water, 'p water'), cbar
	end

	im clf
	clim = [0 4];
	elim = [-0.4 0.4];
	im(241, p_true, 'p true', clim), cbar
	im(242, p_water, 'p water', clim), cbar
	im(243, p_de_fbp, 'p DE FBP', clim), cbar
	im(244, p_de_pwls, 'p DE PWLS', clim), cbar
	im(246, p_water-p_true, 'error', elim), cbar
	title(sprintf('NRMS = %3.1f%%', Nrms(p_true, p_water, 1)))
	im(247, p_de_fbp-p_true, 'error', elim), cbar
	title(sprintf('NRMS = %3.1f%%', Nrms(p_true, p_de_fbp, 1)))
	im(248, p_de_pwls-p_true, 'error', elim), cbar
	title(sprintf('NRMS = %3.1f%%', Nrms(p_true, p_de_pwls, 1)))
%	minmax(p_water-p_true)
%	minmax(p_de-p_true)

return
prompt
end


%
% pwls and fbp errors
%
if 0
	err = [0 0.5];
	im clf
	im(221, xfbp(:,:,1)-xtrue(:,:,1), 'FBP', err), cbar
	ylabel(ftab.mass.type{1})
	im(222, xpwls(:,:,1,end)-xtrue(:,:,1), 'PWLS', err), cbar
	im(223, xfbp(:,:,2)-xtrue(:,:,2), ' ', err), cbar
	ylabel(ftab.mass.type{2})
	im(224, xpwls(:,:,2,end)-xtrue(:,:,2), ' ', err), cbar
	savefig(f.fig, 'fig_fbp_pwls_err')
prompt
	clear err
end


%
% precomputed denominator for PL
%
if ~isvar('denom.ml'), printm 'making denom.ml'
	tic
	[denom.ml wi.ml] = de_pl_denom(G, ymi, ftab.mass.bar, gi);
	printm('denom.ml time %g', toc)
	im clf, im(211, embed(denom.ml, mask), 'ML denom'), cbar
	im(212, wi.ml, 'ML wi'), cbar
prompt
end


%
% penalty for PL recon
%
if ~isvar('Rml'), printm 'making Rml'
	kappa.ml = G' * [col(wi.ml(:,:,1)) col(wi.ml(:,:,2))];
	kappa.ml = sqrt( kappa.ml ./ (G' * ones(n.b*n.a,2)) );
	kappa.ml = embed(kappa.ml, mask);
	im clf, im(211, kappa.ml, 'kappa ML')

	f.l2b_pl = f.l2b - 1;
	f.delta_pl = 0.1;
	f.delta_pl = 0.05;
	f.delta_pl = f.delta_wls;

	for ii=1:2
%		Rml{ii} = Rbuild('leak', kappa.ml(:,:,ii), f.nbrs, ...
%			2^f.l2b_pl, f.pot, f.delta_pl, 0);
		Rml{ii} = Robject(kappa.ml(:,:,ii), ...
			'type_denom', 'matlab', ...
		'beta', 2^f.l2b_pl, 'potential', f.pot, 'delta', f.delta_pl);
	end
	im(212, reshape(Rml{1}.wt, [n.x n.y 4]), 'wt ML'), cbar
prompt
	clear xpl
end


%
% PL reconstruction
%
if ~isvar('xpl'), printm 'xpl'
%	xinit = xpl + 0.1 * rand(size(xpl));
%	xinit = xtrue;			% true init
%	xinit = xfbp;			% fbp init
%	xinit = ones(size(xtrue));	% uniform initial image
	xinit = xpwls(:,:,:,end);	% pwls init!

	if 0, printm 'enforcing object constraints'
		t1 = xinit(:,:,1);
		t2 = xinit(:,:,2);
		t1 = medfilt2(t1, [5 5]);
		t2 = medfilt2(t2, [3 3]);
		t1 = min(t1, 1.1);
		t2 = min(t2, 2.1);
		t2(t2 < 0.4) = 0;
		xinit(:,:,1) = t1;
		xinit(:,:,2) = t2;
	end, clear t1 t2

	xinit = max(xinit,0);
	im(xinit, 'init')

	xtmp = reshape(xinit, n.x*n.y, 2);	% prepare to iterate
	xtmp = xtmp(mask,:);

%	f.niter = 8+1;
	rmi = zeros(size(ymi));		% fix: allow constant...

%	profile on
	xpl = de_pl_osps(xtmp, Gb, ymi, ftab.xray.I, rmi, ...
		Rml, ftab, gtab, ftab.mass.bar, f.niter, ...
		f.pixmax, gi(:), denom.ml);
	xpl = ig.embed(xpl);
%	profile report

	im clf, im(xpl, 'xpl'), cbar

	if 0
		xpl1 = de_pl_osps(xinit, Gb1, ymi, ftab.xray.I, rmi, ...
			Rml, ftab, gtab, ftab.mass.bar, 1, f.niter, ...
			f.pixmax, gi(:), denom.ml);	% 1 subset version
		xpl1 = ig.embed(xpl1);
	end
prompt
	clear xtmp
end

if 0
	printm('fbp error %g', norm(col(xfbp-xtrue)))
	printm('init error %g', norm(col(xinit-xtrue)))
	printm('pl error %g', norm(col(xpl(:,:,:,end)-xtrue)))
end

if 0
	plot_de_true_recon_err(xtrue, xpl(:,:,:,end), 'PL', c)
	savefig(f.fig, 'fig_true_pml_err')
return
prompt
end

if 0
	[t.obj, t.like, t.penal] = ...
		de_pl_obj(xpl, G, ymi, ftab.xray.I, rmi, ftab, Rml, ig.mask);
	if isvar('xpl1')
		t.obj1 = de_pl_obj(xpl1, G, ymi, ftab.xray.I, rmi, ftab, Rml, ig.mask);
	else
		t.obj1 = t.obj;	% lazy
	end
	cost1 = t.obj(1) - t.obj1;
	costn = t.obj(1) - t.obj;
	ii = [0:f.niter-1];
	figure(2)
	im clf, if im, plot(ii, cost1, 'c-o', ii, costn, 'y-x'), end
%	axis([0 16 -6000 0])
	axis tight
	set(gca, 'xtick', [0 4 8 12 16])
%	set(gca, 'ytick', cost1(1:4), 'ygrid', 'on')
	set(gca, 'ytick', [0] * 1000)
	title 'Cost function decrease'
	ylabel '\Psi(x^n) - \Psi(x^0)'
	xlabel 'Iteration n'
	legend(	'# subset = 1', sprintf('# subset = %d', f.nsubset), 1)
%	printc(sprintf([f.fig 'fig_obj_subset%d_vs_1_nx%d'], f.nsubset, n.x))
	for ii=1:3
		hold on
		plot([ii 4*ii], costn(ii+1)*[1 1], 'm--')
		hold off
	end
	figure(1)
end

printm('FBP %g %g', Nrms(xfbp, xtrue, 1), Nrms(xfbp, xtrue, 2))
printm('PWLS %g %g', ...
	Nrms(xpwls(:,:,:,end), xtrue, 1), Nrms(xpwls(:,:,:,end), xtrue, 2))
printm('PL %g %g', ...
	Nrms(xpl(:,:,:,end), xtrue, 1), Nrms(xpl(:,:,:,end), xtrue, 2))

%	error correlation
if 1
	im clf
	ax = [-1 1 -1 1] * 1.5;
	tmp = get(0, 'DefaultLineMarkerSize')
	set(0, 'DefaultLineMarkerSize', 5)
	subplot(131)
	err = reshape(xfbp - xtrue, n.x*n.y,2);
	plot(err(:,1), err(:,2), '.'), title FBP
	xlabel 'soft error', ylabel 'bone error'
	axis(ax), axis square

	subplot(132)
	err = reshape(xpwls(:,:,:,end) - xtrue, n.x*n.y,2);
	plot(err(:,1), err(:,2), '.'), title PWLS
%	xlabel 'soft error', ylabel 'bone error'
	axis(ax), axis square

	subplot(133)
	err = reshape(xpl(:,:,:,end) - xtrue, n.x*n.y,2);
	plot(err(:,1), err(:,2), '.'), title PL
%	xlabel 'soft error', ylabel 'bone error'
	axis(ax), axis square
	set(0, 'DefaultLineMarkerSize', tmp)
	savefig(f.fig, 'fig_fbp_pwls_pml_err')
return
end

% all 4
if 1
	im clf
	im(341, xtrue(:,:,1), 'True', c.soft), cbar([0 1])
	ylabel 'Soft Tissue'
	im(345, xtrue(:,:,2), ' ', c.bone), cbar([0 2])
	ylabel 'Cortical Bone'
	im(342, 'notick', xfbp(:,:,1), 'FBP', c.soft), cbar([0 1])
	im(346, 'notick', xfbp(:,:,2), ' ', c.bone), cbar([0 2])
	im(343, 'notick', xpwls(:,:,1,end), 'PWLS ', c.soft), cbar([0 1])
	im(347, 'notick', xpwls(:,:,2,end), ' ', c.bone), cbar([0 2])
	im(344, 'notick', xpl(:,:,1,end), 'PL ', c.soft), cbar([0 1])
	im(348, 'notick', xpl(:,:,2,end), ' ', c.bone), cbar([0 2])
	savefig(f.fig, 'fig_true_fbp_pwls_pml')
return
end
