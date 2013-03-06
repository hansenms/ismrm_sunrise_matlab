function outp =  ismrm_image_kernel_SPIRiT(inp,im_kernel,transpose_indicator)

if (strcmp(transpose_indicator,'transp')),
    outp = squeeze(sum(repmat(reshape(inp,size(im_kernel,1),size(im_kernel,2),size(im_kernel,3)),[1 1 1 size(im_kernel,4)]) .* conj(permute(im_kernel,[1 2 4 3])),3));
    outp = outp(:) - inp(:);
elseif (strcmp(transpose_indicator, 'notransp')),
    outp = squeeze(sum(repmat(reshape(inp,size(im_kernel,1),size(im_kernel,2),size(im_kernel,3)),[1 1 1 size(im_kernel,4)]) .* im_kernel,3));
    outp = outp(:) - inp(:);
else
    error('Transpose flag not appropriately defined');
end

return
