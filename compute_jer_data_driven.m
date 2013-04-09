function jer_lookup = compute_jer_data_driven(cal_data, kernel_shape)

nc = size(cal_data, 3);

jer_lookup = zeros([kernel_shape kernel_shape nc nc]);


for dky = -(kernel_shape(2)-1):(kernel_shape(2)-1),
    for dkx = -(kernel_shape(1)-1):(kernel_shape(1)-1),
        partial_sums = compute_jer_delta(cal_data, kernel_shape, [dkx dky]);
              
        sum_size_y = kernel_shape(2) - abs(dky);
        sum_size_x = kernel_shape(1) - abs(dkx);
        kx_a_min = max(0, -dkx);
        ky_a_min = max(0, -dky);
                
        for iky = 1:sum_size_y,
            for ikx = 1:sum_size_x,
                        
                kxa_index = kx_a_min + ikx;
                kya_index = ky_a_min + iky;
                jer_lookup(kxa_index, kya_index, kxa_index + dkx, kya_index + dky, :, :) = sum(sum(partial_sums(ikx-1 + (1:sum_size_x), iky-1 + (1:sum_size_y),:,:),2),1);
            end
        end
    end
end

return           
            
            
function partial_sums = compute_jer_delta(cal_data, kernel_shape, delta)
    ndata = [size(cal_data,1) size(cal_data,2)] - abs(delta);
    ngroups = 2*(kernel_shape - abs(delta))-1;
    ngroups_2 = bitshift(ngroups,-1);
    
    ndims = numel(ndata);
    
    sum_matrix = cell(ndims);
    
    for dim_index = 1:ndims,
        sum_matrix{dim_index} = zeros(ngroups(dim_index), ndata(dim_index));
        sum_matrix{dim_index}(ngroups_2(dim_index)+1, (ngroups_2(dim_index)+1):(ndata(dim_index)-ngroups_2(dim_index))) = 1;
    
        for index = 1:ngroups_2(dim_index),
            sum_matrix{dim_index}(index, index) = 1;
            sum_matrix{dim_index}(ngroups_2(dim_index)+1+index, ndata(dim_index)-ngroups_2(dim_index)+index) = 1;
        end
    end
    
    kx_a_min = max(0, -delta(1));
    ky_a_min = max(0, -delta(2));
    
    nc = size(cal_data,3);
    partial_sums = zeros([ngroups nc nc]);
    for ic2 = 1:nc,
        sub_data_b = cal_data(kx_a_min + delta(1) + (1:ndata(1)), ky_a_min+delta(2)+(1:ndata(2)), ic2);
        for ic1 = 1:nc,
            sub_data_a = cal_data(kx_a_min + (1:ndata(1)), ky_a_min + (1:ndata(2)), ic1);
    
            partial_sums(:,:,ic1, ic2) = sum_matrix{1} * (conj(sub_data_a) .* ...
                                         sub_data_b) * sum_matrix{2}.';
        end
    end
    
    
    
    
