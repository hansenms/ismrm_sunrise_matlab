function b = E_SENSE(direction, x, e_args)
%
%  b = simple_multE(direction, x, e_args)
% 
%  Performs the multiplication
%    b = E * x
%  or
%    b = E^H * x
%
%  INPUT:
%   direction   : 'N' or 'H'
%   x           : data
%   e_args.csm  : coil sensitivities
%   e_args.idx  : sample indices
%
%    Michael Schacht Hansen (michael.schacht.hansen@gmail.com)
%

csm = e_args.csm;
idx = e_args.idx;

dims     = size(csm);
coils    = dims(end); 
dims_vec = [prod(dims)/coils coils];

if (direction == 'N'),
    b = zeros(length(idx),coils);
    x = x(:);
    csm = reshape(csm,dims_vec);
    for c=1:coils,
       tmp = x .* csm(:,c);tmp = reshape(tmp,dims(1:end-1));
       tmp = itok(tmp);
       b(:,c) = tmp(idx);
    end
    clear tmp;
elseif (direction == 'H'),
    b = zeros(dims(1:end-1));
    csm = reshape(csm,dims_vec);
    for c=1:coils,
        tmp = zeros(dims(1:end-1));
        tmp(idx) = x(:,c);
        tmp = ktoi(tmp);
        b = b + reshape(conj(csm(:,c)) .* tmp(:), size(b));
    end
    clear tmp;
else
    error('Invalid direction selection');
end    


return