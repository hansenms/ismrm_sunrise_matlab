  function h = titlef(varargin)
%|function h = titlef(varargin)
%| version of title with built-in sprintf
%| also supports default font size from ir_fontsize()

opt = {'fontsize', ir_fontsize('title')};
if ir_is_octave
	opt = {opt{:}, 'fontname', 'Helvetica'};
end

if 1
	tex = {'interpreter', 'tex'};
	str = varargin{1};
	str = strrep(str, '\', '\\');
	varargin{1} = str;
end

if isfreemat
	for ii=1:length(varargin)
		if streq(varargin{ii}, 'interpreter') % not supported by freemat
			varargin{ii} = {};
			varargin{ii+1} = {};
		end
	end
end

if im
	hh = title(sprintf(varargin{:}), tex{:}, opt{:});
else
	hh = [];
end

if nargout
	h = hh;
end
