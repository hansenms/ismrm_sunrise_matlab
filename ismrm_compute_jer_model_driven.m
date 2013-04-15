function jer_lookup = ismrm_compute_jer_model_driven(csm, kernel_shape)
%
%   jer_lookup = ismrm_compute_jer_model_driven(csm, kernel_shape)
%   
%   Computes a lookup table of joint encoding relationships (JER) using the
%   model driven formulation given in Beatty PJ. Reconstruction methods for
%   fast magnetic resonance imaging. PhD thesis, Stanford University, 2006.
%   JERs were previous called "correlation values"; we have changed the
%   name to avoid confusion with correlation coefficients, used to relate
%   two random variables.
%
%
%
%   INPUT:
%       csm [x,y,coil]          : coil sensitivity map (can also be
%                                 weighted coil sensitivity map)
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
nc = size(csm, 3);

jer_lookup = zeros([kernel_shape kernel_shape nc nc]);

nx = size(csm,1);
ny = size(csm,2);
nx_2 = bitshift(nx, -1)+1;
ny_2 = bitshift(ny, -1)+1;

for ic2 = 1:nc,
    for ic1 = 1:nc,
        lookup = ismrm_transform_image_to_kspace(conj(csm(:,:,ic1)) .* csm(:,:,ic2));
        lookup = lookup * sqrt(nx) * sqrt(ny);
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