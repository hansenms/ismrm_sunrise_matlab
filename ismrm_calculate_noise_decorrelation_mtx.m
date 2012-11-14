function dMtx = ismrm_calculate_noise_decorrelation_mtx(noise_samples, noise_dwell_time, sample_dwell_time)
    noise = reshape(noise_samples,numel(noise_samples)/size(noise_samples,length(size(noise_samples))), size(noise_samples,length(size(noise_samples))));
    noise = permute(noise,[2 1]);
    M = size(noise,2);
    dMtx = (1/(M-1))*(noise*noise');
    dMtx = inv(chol(dMtx,'lower'));
    rcvr_noise_bandwidth =0.79; % noise equivalent bandwidth of digital receiver filter
    dMtx = dMtx*sqrt(2*sample_dwell_time/noise_dwell_time*rcvr_noise_bandwidth);
end
