function gmap = ismrm_calculate_gmap(unmixing, ccm, noise_matrix)

if nargin < 3,
    noise_matrix = [];
end
if isempty(noise_matrix),
    noise_matrix = eye(size(unmixing,3));
end

accel_na = ismrm_calculate_noise_amplification(unmixing, noise_matrix); 
full_na  = ismrm_calculate_noise_amplification(ccm, noise_matrix);
gmap = accel_na ./ full_na;    