 function ob = Gblock(ob, nblock, varargin)
%function ob = Gblock(ob, nblock, varargin)

if (nblock > 1) && ...
	(isempty(ob.handle_forw_block) || isempty(ob.handle_back_block))
%	isempty(ob.handle_mtimes_block) &
	fail('The %s of type %s is missing forw/back_block()', ...
		class(ob), ob.caller)
end

if size(ob,1) ~= prod(ob.odim)
	fail 'block fatrix2 needs size(ob,1) == prod(odim)'
end

if nblock == 1
	ob.odim = [ob.odim 1];
end

if nblock > ob.odim(end)
	warn('nblock=%d with dim=[%s]?', nblock, mat2str(ob.odim))
end

ob.nblock = nblock;
ob.caller = sprintf('f:Gblock(%s)', ob.caller);

if ~isempty(ob.handle_block_setup)
	ob = feval(ob.handle_block_setup, ob);
end
