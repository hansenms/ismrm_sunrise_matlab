  function out = full(ob)
% function out = full(ob)
%| full(ob) = ob(:,:)
%|
%| note: alternative is ob * eye(np) but this may yield the
%| wrong precision because eye() is double but ob may be single.
%|
%| Copyright 2010-12-05, Jeff Fessler, University of Michigan

out = fatrix2_subsref_colon(ob, ':'); % ob(:,:)
