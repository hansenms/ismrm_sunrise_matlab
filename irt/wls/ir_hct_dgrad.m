  function [grad out com] = ir_hct_dgrad(x, A, yi, wi, varargin)
%|function [grad out com] = ir_hct_dgrad(x, A, yi, wi, varargin)
%|
%| compute gradient of WLS data-fit term using hct binary
%| for UM testing only
%|
%| in
%| option
%|	(many - see below)
%|
%| Jeff Fessler, 2012-06-17

if nargin == 1 && streq(x, 'test'), ir_hct_dgrad_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

arg.clean = true; % remove files after?
arg.dir = ''; % work directory
arg.file_x = 'ir-hct-dgrad-x.fld'; % work file for x
arg.file_grad = 'ir-hct-dgrad-grad.fld'; % work file for grad
arg.file_mask = 'ir-hct-dgrad-mask.fld'; % work file for mask
arg.file_yi = 'ir-hct-dgrad-yi.fld'; % work file for yi
arg.file_wi = 'ir-hct-dgrad-wi.fld'; % work file for wi
arg.nthread = jf('ncore');
arg.chat = 0;

arg = vararg_pair(arg, varargin);
if isempty(arg.dir)
	arg.dir = test_dir;
end

arg.file_x = [arg.dir arg.file_x];
arg.file_grad = [arg.dir arg.file_grad];
arg.file_mask = [arg.dir arg.file_mask];
arg.file_wi = [arg.dir arg.file_wi];
arg.file_yi = [arg.dir arg.file_yi];

pn = jf_protected_names;
if ~pn.has_hct2
	fail 'need hct2'
end

% write files
%tox = @(z) permute(z, [2 3 1]);
fld_write(arg.file_x, x)
fld_write(arg.file_mask, A.imask, 'type', 'byte')
fld_write(arg.file_yi, yi)
fld_write(arg.file_wi, wi)

if exist(arg.file_grad, 'file')
	eval(['!/bin/rm ' arg.file_grad])
end

tmp = pn.hct_arg(A.cg, A.ig);
com = sprintf(['hct2 what iter chat %d sysmod %s nthread %d %s ' ...
	'file_denom_in %s ' ... % trick: to avoid computing denom!
	'file_yi %s file_wi %s ' ...
	'file_mask %s file_init %s file_dgrad %s'], ...
	arg.chat, A.arg.type, arg.nthread, tmp, ...
	arg.file_x, ... % trick
	arg.file_yi, arg.file_wi, ...
	arg.file_mask, arg.file_x, arg.file_grad);

out = os_run(com); % run hct2

grad = fld_read(arg.file_grad);

if arg.clean
	clean = ['!/bin/rm ' arg.file_x ' ' arg.file_grad ' ' arg.file_mask ' ' arg.file_yi ' ' arg.file_wi];
	eval(clean)
end


% ir_hct_dgrad_test()
function ir_hct_dgrad_test
rng(0)
ig = image_geom('nx', 512, 'nz', 32, 'fov', 500, 'down', 4);
tmp = ig.mask;
tmp([1 end],:,:) = false;
tmp(:,[1 end],:) = false;
ig.mask = tmp; % remove 1-pixel border
cg = ct_geom('ge1', 'nt', 16, 'down', 4);
%cg.plot3(ig)

A = Gcone(cg, ig, 'type', 'sf2');
xb = ellipsoid_im(ig, '');
yb = A * xb;
yi = yb + 4 * randn(size(yb));
wi = exp(-yb * (3 / max(yb(:))));
minmax(wi)
sino = @(y) permute(y, [1 3 2]);

x0 = randn(ig.dim);
W = Gdiag(wi);
mat = A' * W * (A * x0 - yi);

printm 'calling hct'
cpu etic
[hct out com] = ir_hct_dgrad(x0, A, yi, wi);
cpu etoc hct2:dgrad:time
printm(['hct2 output: \n ' out])
%printm(com)

err = hct - mat;
clim = [];
im plc 1 3
im(1, mat, clim)
im(2, hct, clim)
im(3, err)
equivs(hct, mat, 'thresh', 2e-6)
