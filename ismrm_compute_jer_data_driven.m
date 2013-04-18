function jer_lookup = ismrm_compute_jer_data_driven(cal_data, kernel_shape)
%
%   jer_lookup = ismrm_compute_jer_data_driven(cal_data, kernel_shape)
%   
%   Computes a lookup table of joint encoding relationships (JER) using the
%   data driven formulation given in Beatty et al. Proc. ISMRM 2007, p1749.
%   JERs were previous called "correlation values"; we have changed the
%   name to avoid confusion with correlation coefficients, used to relate
%   two random variables.
%
%   This function computes JERs jointly when there is a common delta
%   between the two encoding locations.  This gives some computational
%   reduction.  This function does not take advantage of the symmetry that
%   JER(a,b) = JER(b,a)*.  Also, since it is written in MATLAB with for
%   loops, it's runtime is not fast; it's purpose is mainly pedagogical.
%
%   INPUT:
%       cal_data [kx,ky,coil]   : calibration data (k-space)
%       kernel_shape [2x1]      : kernel shape on a fully sampled grid
%                                 [kx_extent, ky_extent]
%                                 e.g. for acceleration=2, ky_extent=7 
%                                 would use 4 source points along ky; for
%                                 acceleration=4, only 2 source points
%                                 would be used.
%
%   OUTPUT:
%       jer_lookup [kx,ky,coil, kx, ky, coil] : lookup table of all
%       possible relations between kernel sample locations.
%
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%   Philip J. Beatty (philip.beatty@sri.utoronto.ca)
%
assert(nargin==2, '2 arguments needed');
assert(ndims(cal_data)==3, 'cal_data must have 3 dimensions');
assert(length(kernel_shape)==2, 'kernel_shape must be a length 2 vector');

nc = size(cal_data, 3);
cal_shape = [size(cal_data,1) size(cal_data,2)];
    
jer_lookup = zeros([kernel_shape kernel_shape nc nc]);


for dky = -(kernel_shape(2)-1):(kernel_shape(2)-1),
    for dkx = -(kernel_shape(1)-1):(kernel_shape(1)-1),
        partial_sums = compute_partial_sums(cal_data, kernel_shape, [dkx dky]);
              
        nsums = kernel_shape - abs([dkx dky]);
        sum_size = min(nsums, cal_shape - kernel_shape+1);
        kx_a_min = max(0, -dkx);
        ky_a_min = max(0, -dky);
                
        for iky = 1:nsums(2),
            for ikx = 1:nsums(1),
                kxa_index = kx_a_min + ikx;
                kya_index = ky_a_min + iky;
                jer_lookup(kxa_index, kya_index, kxa_index + dkx, kya_index + dky, :, :) = sum(sum(partial_sums(ikx-1 + (1:sum_size(1)), iky-1 + (1:sum_size(2)),:,:),2),1);
            end
        end
    end
end

return           
            
            
function partial_sums = compute_partial_sums(cal_data, kernel_shape, delta)
    cal_shape = [size(cal_data,1) size(cal_data,2)];
    ndata = cal_shape - abs(delta);
    ngroups = min(2*(kernel_shape - abs(delta))-1, cal_shape-abs(delta));
    ngroups_2 = bitshift(ngroups,-1);
    
    ndims = numel(ndata);
    
    sum_matrix = cell(ndims);
    
    for dim_index = 1:ndims,
        matrix_edges = eye(ngroups(dim_index));
        matrix_middle = repmat(matrix_edges(:,ngroups_2(dim_index)+1), [1 ndata(dim_index)-ngroups(dim_index)]);
        sum_matrix{dim_index} = [matrix_edges(:, 1:ngroups_2(dim_index)) matrix_middle matrix_edges(:,(ngroups_2(dim_index)+1):ngroups(dim_index))];
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