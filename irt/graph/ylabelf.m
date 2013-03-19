  function h = ylabelf(varargin)
%|function h = ylabelf(varargin)
%| version of ylabel with built-in sprintf
%| also supports default font size from ir_fontsize()

opt = {'fontsize', ir_fontsize('label')};

if 1
	tex = {'interpreter', 'tex'};
	str = varargin{1};
	str = strrep(str, '\', '\\');
	varargin{1} = str;
end

if im
	hh = ylabel(sprintf(varargin{:}), tex{:}, opt{:});
else
	hh = [];
end

if nargout
	h = hh;
end
