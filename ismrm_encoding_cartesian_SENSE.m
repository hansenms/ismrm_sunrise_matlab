function outp =  ismrm_encoding_cartesian_SENSE(inp,csm,sampling_mask,transpose_indicator)

scale = numel(sampling_mask)/sum(sampling_mask(:) > 0);
if (strcmp(transpose_indicator,'transp')),
    outp = zeros(size(csm));
    outp(repmat(sampling_mask,[1 1 size(csm,3)]) == 1) = inp(:);
    outp = ismrm_transform_kspace_to_image(outp,[1,2])*sqrt(scale);
    outp = sum(conj(csm) .* outp,3);
    outp = outp(:);
elseif (strcmp(transpose_indicator, 'notransp')),
    outp = repmat(reshape(inp,size(csm,1),size(csm,2)),[1 1 size(csm,3)]) .* csm;
    outp = ismrm_transform_image_to_kspace(outp, [1,2])*sqrt(scale);
    outp = outp(repmat(sampling_mask,[1 1 size(csm,3)]) == 1);
else
    error('Transpose flag not appropriately defined');
end

return
