%%
%Clean up
close all;
clear all


kernel_mask = ones(3);
kernel_mask(:,2) = 0;
ncalx = 7;
ncaly = 7;
wx = 5;
wy = 7;
nc = 2;

%cal_data = ones([ncalx, ncaly, nc]);
cal_data = complex(rand([ncalx, ncaly, nc]), rand([ncalx, ncaly,nc]));

jer_lookup = compute_jer_data_driven(cal_data, [wx wy]);

%cal_im = ismrm_transform_kspace_to_image(cal_data, [1,2], 2 * size(cal_data));
%jer_lookup = compute_jer_model_driven(cal_im, [wx wy]);

%kernel = ismrm_compute_kspace_unaliasing_coefficients(jer_lookup, kernel_mask);


%%
% compute reference jers

jer_ref = zeros(size(jer_lookup));

for ic2 = 1:nc,
    for ic1 = 1:nc,
        for kyb = 1:wy,
            for kxb = 1:wx,
                for kya = 1:wy,
                    for kxa = 1:wx,
                        jer_ref(kxa, kya, kxb, kyb, ic1, ic2) = compute_jer_dd_reference(cal_data, kxa, kya, kxb, kyb, ic1, ic2, [wx wy]);
                    end
                end
            end
        end
    end
end

diff = jer_lookup - jer_ref;


[y ind] = max(abs(diff(:)))

% [i1,i2, i3, i4, i5, i6] = ind2sub(size(diff), ind)

%diff = reshape(diff, [wx *wy * nc wx * wy * nc]);
%imshow(abs(diff), []);