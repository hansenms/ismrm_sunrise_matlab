function jer_lookup = ismrm_compute_jer_data_driven_reference(cal_data, kernel_shape)
%
%   jer_lookup = ismrm_compute_jer_data_driven_reference(cal_data, kernel_shape)
%   
%   Computes a lookup table of joint encoding relationships (JER) using the
%   data driven formulation given in Beatty et al. Proc. ISMRM 2007, p1749.
%   JERs were previous called "correlation values"; we have changed the
%   name to avoid confusion with correlation coefficients, used to relate
%   two random variables.
%
%   This function computes each JER independently.  This is not very
%   efficient, but it is useful as a reference for testing accelerated computation
%   approaches.
%
%   INPUT:
%       cal_data [kx,ky,coil]   : Calibration data (k-space)
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

assert(nargin==2, '2 arguments needed');
assert(ndims(cal_data)==3, 'cal_data must have 3 dimensions');
assert(length(kernel_shape)==2, 'kernel_shape must be a length 2 vector');

nc = size(cal_data, 3);
wx = kernel_shape(1);
wy = kernel_shape(2);

nfitx = size(data,1) - kernel_shape(1);
nfity = size(data,2) - kernel_shape(2);
jer_lookup = zeros([kernel_shape kernel_shape nc nc]);

for ic2 = 1:nc,
    for ic1 = 1:nc,
        for kyb = 1:wy,
            for kxb = 1:wx,
                for kya = 1:wy,
                    for kxa = 1:wx,
                        jer_lookup(kxa, kya, kxb, kyb, ic1, ic2) = vec(cal_data( kxa+(0:nfitx), kya+(0:nfity), ic1))' * vec(cal_data( kxb+(0:nfitx), kyb+(0:nfity), ic2));                        
                    end
                end
            end
        end
    end
end