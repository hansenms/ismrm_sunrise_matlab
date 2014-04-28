% test_all_mex
% make sure all mex files can execute
% by running the internal "check" of each.

list = {
	'jf_mex', ...
	'dtft_mex', ...
	'exp_xform_mex', ...
	'mri_exp_mult_mex', ...
...
	'interp1_table_adj_mex', ...
	'interp1_table_mex', ...
	'interp2_table_adj_mex', ...
	'interp2_table_mex', ...
	'interp3_table_adj_mex', ...
	'interp3_table_mex', ...
...
	'penalty_mex', ...
	'rotmex', ...
	'wtfmex', ...
	'f3d_mex', ...
};

% check for UM-only mex files
if exist('dd_ge1_mex') == 3
	list{end+1} = 'dd_ge1_mex';
end
if exist('dd_ge2_mex') == 3
	list{end+1} = 'dd_ge2_mex';
end

is_missing = false(numel(list),1);
for ii=1:numel(list)
	mex = list{ii};
%	pr mex
	if exist(mex) ~= 3
		is_missing(ii) = true;
	end
end

if any(is_missing)
	printf(' ')
	printm('These mex files are missing:')
	disp(list{is_missing})
prompt
	list = list{~is_missing};
end

passed = '';
failed = '';
missing = '';
for ii=1:numel(list)
	mex = list{ii};
	if exist(mex) ~= 3
		missing = [missing ' ' mex];
		continue;
	end
	try
		fun = str2func(mex);
		fun('check')
		passed = [passed ' ' mex];
	catch
		failed = [failed ' ' mex];
	end
end

if ~isempty(missing) || ~isempty(failed)

	if ~isempty(missing)
		printf(' ')
		printm(['These mex files are missing: ' missing])
	end

	if ~isempty(failed)
		printf(' ')
		printm(['These mex files failed: ' failed])
	end

	if ~isempty(passed)
		printf(' ')
		printm(['These mex files passed: ' passed])
		printm 'So perhaps some things will still work.'
	end

	disp 'Sorry, you seem to have mex problems. :-('
	disp 'Probably you are a PC user and Windoze is not supported.'
	disp 'Or (in linux) there may be an annoying gcc library version issue.'
	disp 'Or you may have a path problem.'

else
	printm '------------------------------------------------------------'
	printm(['All mex files present and passed:\n' passed])
	printm '------------------------------------------------------------'
%	printm('All mex files appear to be properly executable.')
end
