% Cdiff1_test.m
% test Cdiff1 object

% test small size first for adjoint etc.
if 1
	ig = image_geom('nx', 8, 'ny', 6, 'dx', 1);
	list_class = {'fatrix2', 'Fatrix'};
	for ic = 1:numel(list_class)
	 for order = 0:2
		pr order
		arg = {ig.dim, 'order', order, 'class', list_class{ic}};
		if order == 0
			arg = {arg{:}, 'offset', 0};
		else
			arg = {arg{:}, 'offset', [1 2]};
		end

		Cc = Cdiff1(arg{:}, 'type_diff', 'convn');
		C1 = Cdiff1(arg{:}, 'type_diff', 'for1');
		Cf = Cdiff1(arg{:}, 'type_diff', 'imfilter');
		Ci = Cdiff1(arg{:}, 'type_diff', 'ind');
		Cm = Cdiff1(arg{:}, 'type_diff', 'mex');
		Cp = Cdiff1(arg{:}, 'type_diff', 'circshift');
		Cs = Cdiff1(arg{:}, 'type_diff', 'sparse');
		Cz = Cdiff1(arg{:}, 'type_diff', 'spmat');

		Cc_f = Cc(:,:);
		C1_f = C1(:,:);
		Cf_f = Cf(:,:);
		Ci_f = Ci(:,:);
		Cm_f = Cm(:,:);
		Cp_f = Cp(:,:);
		Cs_f = Cs(:,:);
		% Cz already a matrix

		switch list_class{ic}
		case 'Fatrix'
			test_fun = @(C) Fatrix_test_basic(C, true(ig.dim), 'halt', 0);
		case 'fatrix2'
			test_fun = @(C) fatrix2_tests(C);
		otherwise
			fail 'bug'
		end

		test_fun(Cc)
		test_fun(C1)
		test_fun(Cf)
		test_fun(Ci)
		test_fun(Cm)
		test_fun(Cp)
		test_fun(Cs)
%		test_fun(Cz) % Cz is matrix not Fatrix!

		if 1 && order > 0 % trick: Cc,Cf,Cp require Rweights to match
			wt = Rweights(ig.mask, Cc.arg.offset, ...
				'type_wt', 'array', ...
               			'order', order, 'distance_power', 0);

			Ci_f_wt = diag(wt) * Ci_f;

			Cc_f_wt = diag(wt) * Cc_f;
			Cf_f_wt = diag(wt) * Cf_f;
			Cp_f_wt = diag(wt) * Cp_f;

			jf_equal(Ci_f_wt, Cc_f_wt)
			jf_equal(Ci_f_wt, Cf_f_wt)
			jf_equal(Ci_f_wt, Cp_f_wt)

			try
				jf_equal(Ci_f_wt, Cp_f_wt)
			catch
				im plc 2 3
				Cp_f_wt_t = Cp_f_wt';
				Ci_f_wt_t = Ci_f_wt';
				im(1, Ci_f)
				im(4, Cp_f)
				im(2, Ci_f_wt_t)
				im(5, Cp_f_wt_t)
				im(3, Cp_f_wt_t - Ci_f_wt_t)
				im(6, ig.shape(wt))
				fail 'Cp bug'
			end
		else
			jf_equal(Ci_f, Cc_f)
			jf_equal(Ci_f, Cp_f)
		end

%		jf_equal(Ci_f, Cc_f) % see wt'd version above
		jf_equal(Ci_f, C1_f)
%		jf_equal(Ci_f, Cf_f) % see wt'd version above
		jf_equal(Ci_f, Cm_f)
%		jf_equal(Ci_f, Cp_f) % see wt'd version above
		jf_equal(Ci_f, Cs_f)
		jf_equal(Ci_f, Cz)

		% run tests on small cases
		Cdiff1_test1(Cc)
		Cdiff1_test1(C1)
		Cdiff1_test1(Cf)
		Cdiff1_test1(Ci)
		Cdiff1_test1(Cm)
		Cdiff1_test1(Cp)
		Cdiff1_test1(Cs)

	 end
	end
end

% timing test for large size:
if im
%	list in fastest to slowest order (on ire):
	list = {'mex', 'circshift', 'imfilter', 'sparse', 'for1', 'ind', 'convn'};
	ig = image_geom('nx', 2^8, 'ny', 2^8, 'nz', 2^7, 'dx', 1);
	for order = 1:2
		printf ' ' % blank line
		printm('order=%d timing tests: [%d %d %d]', ...
			order, ig.nx, ig.ny, ig.nz)
		if order == 0
			arg = {ig.dim, 'order', order, 'offset', 0};
		else
			arg = {ig.dim, 'order', order, 'offset', [3 2 1]};
		end

		x = ig.xg;
		for ii = 1:numel(list)
			C = Cdiff1(arg{:}, 'type_diff', list{ii});
			if 0 % on ire, warm-up is pointless
				cpu etic
				tmp = C * x; % warm up
				time_warm = cpu('etoc');
			else
				time_warm = 0;
			end

			cpu etic
			tmp = C * x;
			time_forw = cpu('etoc');

			cpu etic
			C' * tmp;
			time_back = cpu('etoc');

			if isfield(C.arg, 'type_diff')
				lab = C.arg.type_diff;
			else
				lab = C.caller;
			end
			printf('Cdiff1 %10s\t%5.3f\t%5.3f\t%5.3f', lab, ...
				time_forw, time_back, time_warm)
		end
	end
end
