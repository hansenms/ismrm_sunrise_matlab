function outp =  ismrm_sampling_non_cartesian_SPIRiT(inp,nufft_st,weights,transpose_indicator)

scale = sqrt(prod(prod(nufft_st.Kd))/numel(weights(:)));
if (strcmp(transpose_indicator,'transp')),
    samples = size(nufft_st.om,1);
    coils = numel(inp)/samples;
    inp = reshape(inp,samples,coils);
    outp = (nufft_adj(inp .* repmat(sqrt(weights),[1 coils]),nufft_st)./sqrt(prod(nufft_st.Kd)))*scale;
    %outp = nufft_adj(inp .* repmat(weights,[1 coils]),nufft_st)./sqrt(prod(nufft_st.Nd));
    outp = outp(:);
elseif (strcmp(transpose_indicator, 'notransp')),
    mtx_size = [nufft_st.Nd numel(inp)/prod(nufft_st.Nd)];
    outp = (nufft(reshape(inp,mtx_size),nufft_st)./sqrt(prod(nufft_st.Kd)))*scale;
    outp = outp .*repmat(sqrt(weights),[1 size(outp,2)]);
    outp = outp(:);
else
    error('Transpose flag not appropriately defined');
end

return
