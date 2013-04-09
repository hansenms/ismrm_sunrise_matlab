function cal_data = ismrm_extract_cal_data(data, sp)

    ky_projection = max(sp, [], 1);
    kx_projection = max(sp, [], 2);
    cal_indx = find(kx_projection > 1);
    cal_indy = find(ky_projection > 1);

    kx_cal_bounds = [min(cal_indx):max(cal_indx)];
    ky_cal_bounds = [min(cal_indy):max(cal_indy)];
    cal_data = data(kx_cal_bounds, ky_cal_bounds, :);