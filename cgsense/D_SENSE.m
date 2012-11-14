function D = D_SENSE(fL,lambda,e_args)
%
%  D = simple_SENSE_calcD(fL,lambda,e_args)
% 
%  Calculates the preconditioning matrix D for an iterative reconstruction
%
%  INPUT:
%   fL          : Function handle to regularization matrix
%   lambda      : regularization weight
%   e_args.csm  : coil sensitivities
%
%  Michael Schacht Hansen (michael.schacht.hansen@gmail.com)

D = sum(abs(e_args.csm).^2,length(size(e_args.csm)));

%Make sure there are no zeros
idx = find(D == 0);
if (~isempty(idx)),
   idx_nonzero = find(D > 0);
   D(idx) = min(D(idx_nonzero));
end

D = (D/prod(size(D))) + abs(lambda*fL('H',lambda*fL('N',ones(size(D)),e_args),e_args));

D = sqrt(D) .^ -1;

return