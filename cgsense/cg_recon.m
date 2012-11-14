function [rho,rho_it] = cg_recon(m, fE, e_args, varargin)
% Simple Conjugate gradient reconstruction.
%
%   [rho, rho_it] = cg_recon(m, fE, varargin)
%
%   This function finds the solution to
%     rho = arg min{ ||E*rho-m||2 + lambda*||L(rho-rho_0)||2 }
%   by solving
%     (E^H * E + lambda^2* L^H * L)(rho-rho_0) = E^H * (m - E*rho_0)
%   
%   INPUTS:
%     m    : measured data         
%     fE   : function handle to function implementing multiplication with E.
%   e_args : encoding arguments passed on to function handles
%
%  OPTIONAL PARAMETERS
%   'fD'              : function handle to function which can calculate the
%                       preconditioning matrix D from the avilable parameters.
%   'fL'              : function handle to function implementing multiplication with L.
%   'rho_0'           : initial estimate of rho
%   'lambda'          : regularization factor
%   'limit'           : Residual stopping limit
%   'show_iterations' : '0' or '1'
%   'print_residual'  : '0' or '1'
%
%
%   NOTE:
%     Function handles for E and L must have the following interface:
%     fhandle(direction, data, e_args). Direction is a character indicating
%     normal multiplication 'N' or multiplication with the complex
%     conjugate transpose 'H'.
%     Function hanldes for calculation of preconditioning matrix D have the
%     following interface: fhandle(fL,lambda,e_args).
%   
%
%   Michael Schacht Hansen (michael.schacht.hansen@gmail.com), May 2008
%

%Default input arguments
lambda          = 1;
rho_0           = 0;
fD              = [];
limit           = 1e-5;
maxit           = 100;
show_iterations = 1;
print_residual  = 1;
fL              = @L_default;

%Process the input arguments
for arg=[1:2:length(varargin)],
   if (strcmp(lower(varargin{arg}),'limit')),
       limit = varargin{arg+1};
   elseif (strcmp(lower(varargin{arg}),'iterations')),
       maxit = varargin{arg+1};
   elseif (strcmp(lower(varargin{arg}),'fd')),
       fD = varargin{arg+1};
   elseif (strcmp(lower(varargin{arg}),'fl')),
       fL = varargin{arg+1};
   elseif (strcmp(lower(varargin{arg}),'lambda')),
       lambda = varargin{arg+1};
   elseif (strcmp(lower(varargin{arg}),'rho_0')),
       rho_0 = varargin{arg+1};
   elseif (strcmp(lower(varargin{arg}),'show_iterations')),
       show_iterations = varargin{arg+1};
   elseif (strcmp(lower(varargin{arg}),'print_residual')),
       print_residual = varargin{arg+1};
   end
end

%Form the right hand side of the equation
if (rho_0 ~= 0),
    rhs = fE('H', m - fE('N',rho_0,e_args), e_args);
else
    rhs = fE('H', m, e_args);
end

%Calculate preconditioning matrix
if (isempty(fD)),
    D = ones(size(rhs));
else
    D = fD(fL,lambda,e_args);
end

%Initial residual, we will have a starting guess of rho - rho_0 = 0.
r = D.*rhs;

rr_0 = r(:)'*r(:);
rr = 0;

%Run iteration
rho = 0;
fprintf('Iterating...\n');
for it=[1:maxit],
    rr_1 = rr;
    rr = r(:)'*r(:);
    
    r_it(it) = rr;
    
    if (it == 1),
        p = r;
    else        
        beta = rr/rr_1;
        p =  r + beta*p;    
    end

    %Multiply with system matrix
    %q = D*(E^H * E + lambda^2 * L^H * L)*D*p
    q = D .* (fE('H',fE('N',D .* p,e_args),e_args) + lambda*fL('H',lambda*fL('N',D.*p,e_args),e_args));
    alpha = rr/(p(:)'*q(:));
    rho = rho + alpha*p;
    r = r - alpha*q;

    if (show_iterations == 1),
        imagesc(abs(D .*rho + rho_0)); axis image; drawnow;
    end
    
    if (nargout > 3),
        rho_it(:,:,it) = D .*rho + rho_0;
    end
    
    if (print_residual == 1),
        fprintf('Iteration %d, norm(r) /norm(r_0) = %12.8e\n', it, rr/rr_0);drawnow;
    end
    
    
    if (rr/rr_0 < limit);
       break;
    end
end

rho = D.*rho + rho_0;

return


function b = L_default(direction, x, e_args)
b = 0; %default is no regularization
return