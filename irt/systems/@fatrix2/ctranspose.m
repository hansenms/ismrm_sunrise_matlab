 function ob = ctranspose(ob)
%function ob = ctranspose(ob)
% "ctranspose" method for this class

ob.size = fliplr(ob.size);
[ob.odim ob.idim] = deal(ob.idim, ob.odim);
[ob.omask ob.imask] = deal(ob.imask, ob.omask);

[ob.handle_forw ob.handle_back] = deal(ob.handle_back, ob.handle_forw);
[ob.handle_forw_block ob.handle_back_block] = ...
	deal(ob.handle_back_block, ob.handle_forw_block);

ob.scale = conj(ob.scale);
