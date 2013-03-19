  function x = ir_project_simplex(y)
%|function x = ir_project_simplex(y)
%|
%| project an n-dim vector y to the simplex Dn
%| Dn = { x : x n-dim, 1 >= x >= 0, sum(x) = 1}
%|
%| (c) Xiaojing Ye
%| xyex19@gmail.com
%|
%| Algorithm is explained as in the linked document
%| http://arxiv.org/abs/1101.6081
%| or
%| http://ufdc.ufl.edu/IR00000353/
%|
%| Jan. 14, 2011.
%|
%| 2012-06-08, modified by JF to include built-in help and example.
%| todo: is there a way to do this without using a loop?

if nargin < 1, help(mfilename), error(mfilename), end
if streq(y, 'test'), ir_project_simplex_test, return, end

m = length(y); bget = false;

s = sort(y,'descend'); tmpsum = 0;

for ii = 1:m-1
	tmpsum = tmpsum + s(ii);
	tmax = (tmpsum - 1)/ii;
	if tmax >= s(ii+1)
		bget = true;
	break
	end
end
    
if ~bget, tmax = (tmpsum + s(m) -1)/m; end;

x = max(y-tmax,0);


% ir_project_simplex_test()
function ir_project_simplex_test

m = 30; % # of trials
rng(0)
zz = 4*(rand(2,m)-0.5);
args = {[0 1], [1 0], '-'};
for ii=1:m
	z = zz(:,ii);
	w = ir_project_simplex(z);
	args = {args{:}, [w(1) z(1)], [w(2) z(2)], '-'};
end
if im
	clf
	plot(args{:})
	axis([-1 1 -1 1]*2), xtick([0 1]), ytick([0 1]), grid
end
