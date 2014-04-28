 function pot_new = potential_shift(pot_orig, vector)
%function pot_new = potential_shift(pot_orig, vector)
%
% Change pot_orig([Cx]_k) to shifted form: pot_new( [Cx]_k - vector_k )
%
% out
%	inline functions:
%	pot.potk(pot, C*x)	potential function value
%	pot.wpot(pot, C*x)	potential 'weights' (aka half-quad. curvatures)
%	pot.dpot(pot, C*x)	potential derivative
%
% Copyright 2005-4-24, Jeff Fessler, The University of Michigan

if nargin == 1 && streq(pot_orig, 'test'), potential_shift_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

pot_new.pot = pot_orig;
pot_new.vector = vector;

pot_new.potk = inline('new.pot.potk(new.pot, t - new.vector)', 'new', 't');
pot_new.wpot = inline('new.pot.wpot(new.pot, t - new.vector)', 'new', 't');
pot_new.dpot = inline('new.pot.dpot(new.pot, t - new.vector)', 'new', 't');


%
% test routine
%
function potential_shift_test

delta = 10;
t = linspace(-4*delta,4*delta,201)';
type = 'hyper3';
pot = potential_func(type, delta);
new = potential_shift(pot, [15]);

pp = pot.potk(pot, t);
pw = pot.wpot(pot, t);
pd = pot.dpot(pot, t);
np = new.potk(new, t);
nw = new.wpot(new, t);
nd = new.dpot(new, t);

if im
	clf
	subplot(311), plot(t, pp, '-', t, np, '--'), title 'potk'
	subplot(312), plot(t, pw, '-', t, nw, '--'), title 'wpot'
	subplot(313), plot(t, pd, '-', t, nd, '--'), title 'dpot'
end
