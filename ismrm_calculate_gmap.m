function gmap = ismrm_calculate_gmap(unmixing, ccm)

    gmap = sqrt( sum(abs(unmixing).^2,3) ./ sum(abs(ccm).^2,3) );