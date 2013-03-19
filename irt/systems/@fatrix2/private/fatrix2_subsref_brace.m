 function out = fatrix2_subsref_brace(ob, iblock)
%function out = fatrix2_subsref_brace(ob, iblock)
% handle subscript references like ob{iblock}
% This will called from ../subsref with ob
% Copyright 2010-12-17, Jeff Fessler, University of Michigan

out = ob;
if isempty(ob.nblock), error 'not a block object', end

% note: user can use subset_starts to select blocks in other orders
% in which case this "iblock" is really "istart"
out.iblock = iblock;
if iblock < 1 || iblock > ob.nblock
	error 'bad block index'
end

odim = ob.odim;
na = odim(end);
ia = out.iblock:ob.nblock:na;
out.odim = [odim(1:end-1) length(ia)];

if ~isempty(out.omask)
	fail 'block object with omask unsupported'
end
out.size = [prod(out.odim) ob.size(2)];
