  function out = ir_fontsize(varargin)
%|function out = ir_fontsize(varargin)
%|
%| set / return default font size for labels, title, etc
%|	'label'		return / set default for label
%|	'title'		return / set default for title
%|	'im_axes'	return / set default for im() axes
%|	'text'		return / set default for text
%|
%| Copyright 2012-12-30, Jeff Fessler, University of Michigan

persistent fs % store preference
if ~isvar('fs') || isempty(fs)
	fs.label = 15;
	fs.title = 15;
	fs.text = 12;
	fs.im_axes = 7;
end

if ~nargin, help(mfilename), pr fs, error(mfilename), end

switch nargin
case 1
	out = fs.(varargin{1});

case 2
	if streq(varargin{1}, 'set')
		val = varargin{2};
		if ischar(val)
			val = str2num(val);
		end
		fs.label = val;
		fs.title = val;
		fs.text = val;
	else
		val = varargin{2};
		if ischar(val)
			val = str2num(val);
		end
		fs.(varargin{1}) = val;
	end

otherwise
	fail('bug')
end
