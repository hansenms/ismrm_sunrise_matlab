function b = L_weight(direction, data, e_args)
%
%  b = L_weight(direction, data, e_args)
% 
%  Performs the multiplication
%    b = L * x
%  or
%    b = L^H * x
%
%  INPUT:
%   direction      : 'N' or 'H' (not used)
%   x              : data
%   e_args.weights : Not used
%
%   Michael Schacht Hansen (michael.schacht.hansen@gmail.com)
%

b = data .* e_args.weights;