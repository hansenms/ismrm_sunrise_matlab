function na = ismrm_calculate_noise_amplification(unmixing, noise_matrix)

nx = size(unmixing, 1);
ny = size(unmixing, 2);
nc = size(unmixing, 3);

na = zeros([nx ny]);
for ic2 = 1:nc,
    for ic1 = 1:nc,
        na = na + noise_matrix(ic1, ic2) .* unmixing(:,:,ic1) .* conj(unmixing(:,:,ic2));
    end
end

na = sqrt(abs(na));