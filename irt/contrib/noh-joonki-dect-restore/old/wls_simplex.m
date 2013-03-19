 function x = wls_simplex(A, y, Wh, x, varargin)
%function x = wls_simplex(A, y, Wh, x, [options])
%| min_x || Wh * (A x - y) ||
%| subject to simplex constraint: 0 <= x <= 1 and sum(x) = 1
%|
%| based on:
%| x = lsqlin(C,d,A,b,Aeq,beq) solves the least-squares
%| (with equality constraints) problem:
%| min_x 0.5*(norm(C*x-d)).^2 subject to A*x <= b and Aeq*x = beq
%| x = lsqlin(C,d,A,b,Aeq,beq,LB,UB) defines a set of lower and upper
%| bounds on the design variables, x, so that LB <= x <= UB.
%|
%| option
%|	'inprodv'	if (row vector) provided, require 1 = inprodv * x
%|
%| Copyright 2006-1-1, Jeff Fessler, University of Michigan

if nargin < 2, help(mfilename), error(mfilename), end

arg.inprodv = [];
arg = vararg_pair(arg, varargin);

n = ncol(A);

if ~isvar('Wh') || isempty(Wh)
	Wh = 1;
end

if ~isvar('x') || isempty(x)
	x = ones(n,1) / n;
else
	x = max(x, 0);
	x = x / sum(x);
end

opt = optimset('largescale', 'off', 'display', 'off');
if isempty(arg.inprodv)
	Aeq = ones(1,n);
	Beq = 1;
else
	Aeq = [ones(1,n); arg.inprodv];
	Beq = [1; 1];
end
x = lsqlin(Wh * A, Wh * y, [], [], Aeq, Beq, zeros(1,n), ones(1,n), x, opt);
