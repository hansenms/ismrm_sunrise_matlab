% test_all_mex
% make sure all mex files can execute
% by running the internal "check" of each.

list = {
	@jf_mex
	@dtft_mex
	@exp_xform_mex
	@mri_exp_mult_mex

	@interp1_table_adj_mex
	@interp1_table_mex
	@interp2_table_adj_mex
	@interp2_table_mex
	@interp3_table_adj_mex
	@interp3_table_mex

	@penalty_mex
	@rotmex
	@wtfmex
	@f3d_mex
};

% check for UM-only mex files
if exist('dd_ge1_mex') == 3
	list{end+1} = @dd_ge1_mex;
end
if exist('dd_ge2_mex') == 3
	list{end+1} = @dd_ge2_mex;
end

passed = '';
failed = '';
for ii=1:length(list)
	mex = list{ii};
	try
		mex('check')
		passed = [passed ' ' func2str(mex)];
	catch
		failed = [failed ' ' func2str(mex)];
	end
end

if ~isempty(failed)
	printf(['These mex files failed: ' failed])
	disp 'Sorry, you seem to have mex problems. :-('
	disp 'Probably you are a PC user and Windoze is not supported.'
	disp 'Or (in linux) there may be an annoying gcc library version issue.'
	disp 'Or you may have a path problem.'
	
	if ~isempty(passed)
		printf(['These mex files passed: ' passed])
		printm 'So perhaps some things will still work.'
	end

else
	printm '------------------------------------------------------------'
	printm(['All mex files passed:\n' passed])
	printm '------------------------------------------------------------'
%	printm('All mex files appear to be properly executable.')
end
