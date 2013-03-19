  function out = ir_is_octave
%|function out = ir_is_octave
%| determine if this is octave
tmp = version;
out = streq(tmp, '3.6.3');
