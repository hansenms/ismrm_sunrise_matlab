function dmtx = ismrm_calculate_noise_decorrelation_mtx_from_covariance_mtx(noise_covariance_mtx)

dmtx = inv(chol(noise_covariance_mtx, 'lower'));