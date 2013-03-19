  function ob = cat(dim, varargin)
%|function ob = cat(dim, varargin)
%|
%| called for ob = cat(dim, ob1, ob2, ...) where any of them is a fatrix2
%|
%| vertcat(ob1,ob2) = [ob1; ob2] = cat(1, ob1, ob2)
%| horzcat(ob1,ob2) = [ob1, ob2] = cat(2, ob1, ob2)
%|
%| dim > 2 is also supported thanks to 'odim'

switch dim
case 1
	ob = fatrix2_vertcat(varargin);
case 2
	ob = fatrix2_horzcat(varargin);
otherwise
	fail('cat(%d, [fatrix2]) is meaningless', dim)
end
