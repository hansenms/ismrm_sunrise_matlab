function ccm = ismrm_compute_ccm(csm, noise_matrix)

nx = size(csm, 1);
ny = size(csm, 2);
nc = size(csm, 3);

if( nargin <2 || isempty(noise_matrix) )
    noise_matrix = eye(nc);
end

csm_matrix = reshape(csm, [nx*ny nc]);

relative_ccm = conj(csm_matrix) * noise_matrix;

ccm = relative_ccm ./ repmat(sum(relative_ccm .* csm_matrix, 2), [1 nc]);

ccm = reshape(ccm, [nx ny nc]);