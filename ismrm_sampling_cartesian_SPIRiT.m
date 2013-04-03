function outp =  ismrm_sampling_cartesian_SPIRiT(inp,sampling_mask,transpose_indicator)

scale = numel(sampling_mask)/sum(sampling_mask(:) > 0);
if (strcmp(transpose_indicator,'transp')),
    coils = numel(inp)/sum(sampling_mask(:));
    outp = zeros(size(sampling_mask,1),size(sampling_mask,2),coils);
    outp(repmat(sampling_mask,[1 1 coils]) == 1) = inp(:);
    outp = ismrm_transform_kspace_to_image(outp,[1,2])*sqrt(scale);
    outp = outp(:);
elseif (strcmp(transpose_indicator, 'notransp')),
    coils = numel(inp)/numel(sampling_mask(:));
    outp = ismrm_transform_image_to_kspace(reshape(inp,size(sampling_mask,1),size(sampling_mask,2),coils), [1,2])*sqrt(scale);
    outp = outp(repmat(sampling_mask,[1 1 coils]) == 1);
else
    error('Transpose flag not appropriately defined');
end

return
