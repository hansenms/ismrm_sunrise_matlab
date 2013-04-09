function jer_lookup = compute_jer_model_driven(csm, kernel_shape)

nc = size(csm, 3);

jer_lookup = zeros([kernel_shape kernel_shape nc nc]);

nx_2 = bitshift(size(csm,1), -1)+1;
ny_2 = bitshift(size(csm,2), -1)+1;

for ic2 = 1:nc,
    for ic1 = 1:nc,
        lookup = ismrm_transform_image_to_kspace(conj(csm(:,:,ic1)) .* csm(:,:,ic2));
        
        for ikyb = 1:kernel_shape(2),
            for ikxb = 1:kernel_shape(1),
                for ikya = 1:kernel_shape(2),
                    for ikxa = 1:kernel_shape(1),
                        jer_lookup(ikxa, ikya, ikxb, ikyb, ic1, ic2) = lookup(nx_2 + ikxb - ikxa, ny_2 + ikyb-ikya);
                    end
                end
            end
        end
    end
end