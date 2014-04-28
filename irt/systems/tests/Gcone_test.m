% Gcone_test.m
% test the Gcone object
% Copyright 2008-1-1, Jeff Fessler, University of Michigan

printm 'todo: check user_source_zs vs hct2 with file_source_z' 
% todo: sf0

pn = jf_protected_names;

f.sf1 = {'sf1', ...
	'sf1,p:ts,b:st', 'sf1,p:st,b:st', ...
	'sf1,p:ts,b:ts', 'sf1,p:st,b:ts'};
f.sf2 = {'sf2', ...
	'sf2,p:ts,b:st', 'sf2,p:st,b:st', ...
	'sf2,p:ts,b:ts', 'sf2,p:st,b:ts'};
systypes = { ... % lots of variations on system models!
	f.sf1{:}, f.sf2{:}, ...
	'sf1', 'sf2', 'sf3', 'sf4', 'sf5' ...
};
if 0 && exist('dd_ge1_mex') == 3 % UM only
	systypes{end+1} = 'dd1';
	systypes{end+1} = 'dd2'; % todo
end
% systypes = {'nn1', 'pd1'};
%systypes = {'dd2'}; % todo! fails adjoint test!?
%systypes = {'sf2'};
%systypes = f.sf2;
nn = numel(systypes);

%f.class = 'Fatrix';
f.class = 'fatrix2';

% small systems for basic tests
if 0 || ~isvar('A1'), printm 'setup small'
	f.down = 16;

	dfs_list = [0 inf inf]; % arc flat parallel
	dsd_list = [949.075 949.075 inf]; % arc flat parallel

% todo
dfs_list = [0]; % arc
dsd_list = [949.075]; % arc

	for kk = 1:numel(dfs_list)

		cgs = ct_geom('ge2', 'nt', 320, ...
			'source_z0', -20, 'pitch', 0.5, ... % test helix
			'dfs', dfs_list(kk), ... % arc or flat
			'dsd', dsd_list(kk), ... % fan or parallel beam
			'down', f.down);

		for ii=1:nn
			systype = systypes{ii};
			printm('testing type %s dfs=%g dsd=%g', ...
				systype, cgs.dfs, cgs.dsd)

			if streq(systype, 'dd1') || streq(systype, 'dd2')
				f.dy = 1; % DD requires square pixels
				if isinf(cgs.dsd) % no parallel for DD
					continue
				end
			else
				f.dy = 0.7; % stress test SF with non-square
f.dy = 'dx'; % todo
			end

			igs = image_geom('nx', 512, 'ny', 480, 'nz', 416, ...
				'dx', 1, 'dy', f.dy, 'dz', 0.5, ...
				'offset_x', 2.9, 'offset_y', 3.5, ...
				'offset_z', -3.4, ...
				'mask', 'all-but-edge', ... % trick: for hct2
				'down', f.down);

			if im, cgs.plot(igs); end

			args = {cgs, igs, 'type', systype, 'class', f.class};
			if streq(systype, 'sf4') || streq(systype, 'sf5')
				args = {args{:}, 'mexarg', {single(0)}};
			end

			% todo: test both ways for fatrix2, not for Fatrix
			is_zxy = true;
			is_zxy = false;
			args = {args{:}, 'zxy', is_zxy};
			if is_zxy
				to_zxy = @(x) permute(x, [3 1 2]);
			else
				to_zxy = @(x) x;
			end

			A1 = Gcone(args{:}, 'nthread', 1);
			Ac = Gcone(args{:});

			tester_tomo2(A1, to_zxy(igs.mask), 'G2', Ac) % paces
			test_adjoint(A1, 'big', 1, 'tol', 5e-5)
			test_adjoint(Ac, 'big', 1, 'tol', 5e-5)

			if streq(systype, 'dd1') || streq(systype, 'dd2')
				continue % don't bother hct2 test for DD
			end

			% hereafter is hct2 test
			Ah = Gcone(args{:}, 'use_hct2', 1);

			if pn.has_hct2
				if streq(systype, 'nn1') || streq(systype, 'pd1')
					thresh = 3e-2; % big due to rounding in nn1
				else
					thresh = 4e-6;
				end

				if 1 % mex vs hct2
					xs = single(igs.mask);
					xs(round(end/3), round(2*end/3), round(end/4)) = 10;
					xs = to_zxy(xs);
					t1 = A1 * xs;
					t2 = Ah * xs;
					equivs(t1, t2, 'thresh', thresh)

					b1 = A1' * t1;
					b2 = Ah' * t1;
					equivs(b1, b2, 'thresh', thresh)
				end
				if 0 % block
					B1 = Gblock(A1, 2);
					Bh = Gblock(Ah, 2);
					t1 = B1{2} * xs;
					t2 = Bh{2} * xs;
					equivs(t1, t2, 'thresh', thresh)
				end

				tester_tomo2(Ah, to_zxy(igs.mask), ...
					'equiv_thresh', thresh) % paces
				test_adjoint(Ah, 'big', 1, 'tol', 5e-5)
			end
		end % systype
	end % dfs
end % small


if 0 % todo: test sf0 vs dd2
	A0 = Gcone(args{:}, 'type', 'sf0');
	Ad = Gcone(args{:}, 'type', 'dd2');
	x0 = single(igs.mask);
	im(x0)
	cpu etic
	yd = Ad * x0;
	cpu etoc dd
	if 0
		cpu etic
		y0 = A0 * x0; % todo: crashes matlab
		cpu etoc sf0
	end
	im plc 1 3
	im(1, yd), im(2, y0), im(3, y0-yd)
return
end


if ~isvar('x0'), printm 'x0 big'
	f.down = 4;
	igb = image_geom('nx', 512, 'ny', 480, 'nz', 416, ...
		'dx', 1, 'dz', 0.5, ...
'dy', 0.7, ... % stress test
		'offset_x', 12.9, 'offset_y', 3.5, 'offset_z', -3.4, ...
		'down', f.down);
	ell = [3*igb.dx 5*igb.dx -2*igb.dz ...
		igb.dx*igb.nx/3 igb.dy*igb.ny/4 igb.zfov/4 ...
		0 0 10];
	x0 = ellipsoid_im(igb, ell, 'oversample', 2);
%prompt
end


% big systems for accuracy tests
if ~isvar('Ab'), printm 'setup big'
	cgb = ct_geom('ge1', 'nt', 320, ...
		'source_z0', -40, 'pitch', 0.5, ... % test helix
		'down', f.down);
%		'dfs', inf, ... % flat detector
%		'dsd', inf, 'dfs', inf, ... % parallel beam

	clear Ab Ah
	for ii=1:nn
		systype = systypes{ii};
		if streq(systype, 'dd', 2) && isinf(cgb.dsd)
			Ab{ii} = [];
			continue
		end
		Ab{ii} = Gcone(cgb, igb, 'type', systype, ...
				'zxy', is_zxy, 'class', f.class);
	end
end

if ~isvar('ya'), printm 'analytical projections'
	ya = ellipsoid_proj(cgb, ell, 'oversample', 2);
%	im clf, im(ya)
%prompt
end


if ~isvar('yb'), printm 'discrete projections'
	nrmse = @(x,y) norm(y(:)-x(:)) / norm(x(:)) * 100;
	for ii=1:nn
		if ~isempty(Ab{ii})
			cpu etic
			yb{ii} = Ab{ii} * to_zxy(x0);
			f.time(ii) = cpu('etoc');
			printm('%s: time %g nrmse %g %%', ...
				systypes{ii}, f.time(ii), nrmse(ya, yb{ii}))
		else
			f.time(ii) = 0;
			yb{ii} = cgb.zeros;
		end
	end
end

im_toggle(ya(:,end/2,:), yb{1}(:,end/2,:), ...
		yb{1}(:,end/2,:) - ya(:,end/2,:))
im(ya(:,end/2,:) - yb{1}(:,end/2,:))

if 0, printm 'look at error in worst views'
	im('pl', 2, nn)
	for ii=1:numel(systypes)
		err = yb{ii} - ya;
		tmp = reshape(err, [], cgb.na);
		tmp = sum(tmp.^2); % error in each view
		ia = imax(tmp); % worst view
		im(ii, err(:,:,ia)), cbar h
		titlef('%s ia=%d', systypes{ii}, ia)
		im('subplot', ii+nn)
		plot(tmp), axis tight
	end
return
end

if 1, printm 'projection profiles'
	it = cgb.nt;
	it = round(cgb.nt/2); % middle
	it = it + [-2 0 2];
	ia = imin(abs(cgb.ad - 45));
%	ia = ceil(ia/2);
	pro = @(y) col(y(:,it,ia));
	arg = [];
	for ii=1:numel(systypes)
		arg = [arg pro(yb{ii})];
	end
	if im
		clf, plot([arg pro(ya)])
		text(10, 200, sprintf('ia=%d', ia))
		text(10, 400, sprintf('ang=%g', cgb.ad(ia)))
		legend(systypes{:}, 'true')
		axisy(0, 1.2 * max(ya(:)))
		grid
	end
return
end


if 0 % dd1 vs dd2 - they match well
	i_dd1 = strmatch('dd1', systypes);
	i_dd2 = strmatch('dd2', systypes);
	im clf, im(yb{i_dd1} - yb{i_dd2}), cbar
	equivs(yb{i_dd1}, yb{i_dd2}, 'thresh', 2e-5)
return
end

if 1, printm 'show projections and differences'
	im clf, im('pl', 2, 1+nn)
	im(1, x0)
	ia = round([1 cgb.na/4 cgb.na/2 cgb.na]);
	im(nn+2, ya(:,:,ia))

	for ii=1:nn
		tmp = yb{ii};
		im(ii+1, tmp(:,:,ia))
		xlabel(systypes{ii})
		im(ii+2+nn, tmp(:,:,ia)-ya(:,:,ia))
	end
prompt
end

if 1, printm 'show back-projections'
	im clf, im('pl', 1, nn)
	iz = round([1 igb.nz/4 igb.nz/2 igb.nz]);
	iz = 1:2:igb.nz;
	for ii=1:nn
%		tmp = cgb.ones;
		tmp = cgb.zeros; tmp(:,:,20) = 1;
		tmp = Ab{ii}' * tmp;
		im(ii, tmp(:,:,iz))
		xlabel(systypes{ii})
	end
end

if 0
	im pl 1 3
	ia = round([1 cg.na/4 cg.na/2 cg.na]);
	im row 4
	im(1, ya(:,:,ia));
	im(2, yc(:,:,ia));
	im(3, yc(:,:,ia)-ya(:,:,ia));
	im reset
end
