function b = L_std_Tikh(direction, data, e_args)
%
%  b = L_std_Tikh(direction, data, e_args)
% 
%  Performs the multiplication
%    b = L * x
%  or
%    b = L^H * x
%
%  INPUT:
%   direction   : 'N' or 'H' (not used in this case)
%   x           : data
%   e_args  :    Not used
%
%   Michael Schacht Hansen (michael.schacht.hansen@gmail.com)
%

b = data;


return