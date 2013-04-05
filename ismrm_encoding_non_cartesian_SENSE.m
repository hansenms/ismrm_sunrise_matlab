function outp =  ismrm_encoding_non_cartesian_SENSE(inp,csm,nufft_st,weights,transpose_indicator)

scale = sqrt(prod(prod(nufft_st.Kd))/numel(weights(:)));
if (strcmp(transpose_indicator,'transp')),
    samples = size(nufft_st.om,1);
    coils = numel(inp)/samples;
    inp = reshape(inp,samples,coils);
    outp = (nufft_adj(inp .* repmat(sqrt(weights),[1 coils]),nufft_st)./(sqrt(prod(nufft_st.Kd))))*scale;
    outp = sum(conj(csm) .* outp,3);
    outp = outp(:);
elseif (strcmp(transpose_indicator, 'notransp')),
    outp = repmat(reshape(inp,size(csm,1),size(csm,2)),[1 1 size(csm,3)]) .* csm;
    outp = (nufft(outp,nufft_st)./(sqrt(prod(nufft_st.Kd))))*scale;
    outp = outp .*repmat(sqrt(weights),[1 size(outp,2)]);
    outp = outp(:);
else
    error('Transpose flag not appropriately defined');
end

return

