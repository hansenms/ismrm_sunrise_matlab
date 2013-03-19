function outp =  ismrm_system_non_cartesian_SPIRiT(inp,nufft_st,weights,im_kernel,transpose_indicator)

if (strcmp(transpose_indicator,'transp')),
    coils = size(im_kernel,3);
    ele = (numel(weights)*coils);
    inp1 = inp(1:ele);
    inp2 = inp((ele+1):end);
    outp1 = ismrm_sampling_non_cartesian_SPIRiT(inp1,nufft_st,weights,transpose_indicator);
    outp2 = ismrm_image_kernel_SPIRiT(inp2,im_kernel,transpose_indicator);
    outp = outp1(:) + outp2(:);
elseif (strcmp(transpose_indicator, 'notransp')),
    outp1 = ismrm_sampling_non_cartesian_SPIRiT(inp,nufft_st,weights,transpose_indicator);
    outp2 = ismrm_image_kernel_SPIRiT(inp,im_kernel,transpose_indicator);
    outp = [outp1;outp2];
else
    error('Transpose flag not appropriately defined');
end

return
