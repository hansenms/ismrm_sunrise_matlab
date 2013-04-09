function kernel = ismrm_compute_kspace_unaliasing_coefficients(jer_lookup, kernel_mask)


[kx_offsets, ky_offsets] = ind2sub(size(kernel_mask), find(kernel_mask == 1));

num_source = length(kx_offsets);
num_channel = size(jer_lookup,5);

Rss = zeros(num_source, num_channel, num_source, num_channel);
Rst = zeros(num_source, num_channel, num_channel);


for is2 = 1:num_source,
    for is1 = 1:num_source,
        Rss(is1, :, is2, :) = jer_lookup(kx_offsets(is1), ...
                                         ky_offsets(is1), ...
                                         kx_offsets(is2), ...
                                         ky_offsets(is2), ...
                                         :, :);
    end
end



kx_target = bitshift(size(kernel_mask,1),-1)+1;
ky_target = bitshift(size(kernel_mask,2),-1)+1;


for is1 = 1:num_source,
    Rst(is1, :, :) = jer_lookup(kx_offsets(is1), ...
                                ky_offsets(is1), ...
                                kx_target, ...
                                ky_target, ...
                                :, :);
            
end


num_basis = num_source * num_channel;
Rss = reshape(Rss, [num_basis num_basis]);
Rst = reshape(Rst, [num_basis num_channel]);
w = (Rss + eye(num_basis)*0.01 * trace(Rss) / num_basis ) \ Rst;

kernel = repmat(kernel_mask,[1 1 num_channel num_channel]);
kernel(kernel == 1) = w(:);
