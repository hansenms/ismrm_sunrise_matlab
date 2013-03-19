function [rho,rho_it] = ismrm_cg(m, fE, varargin)
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
%     m    : Measured data         
%     fE   : Function handle to function implementing multiplication with E
%            or its complex conjugate transpose
%            Function protypt is A(x,'transp' or 'notransp')
%
%  OPTIONAL PARAMETERS
%   'D'               : Diagonal preconditioner, vector with diagonal
%                       elements
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
%     fhandle(data, direction). direction is a character string indicating
%     normal multiplication 'notransp' or multiplication with the complex
%     conjugate transpose 'transp'.
%     Function hanldes for calculation of preconditioning matrix D have the
%     following interface: fhandle(x).
%   
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%   Philip Beatty (philip.beatty@sri.utoronto.ca)
%

%Default input arguments
lambda          = 1;
rho_0           = 0;
D               = [];
limit           = 1e-5;
maxit           = 100;
show_iterations = 0;
print_residual  = 1;
fL              = @L_default;

%Process the input arguments
for arg=[1:2:length(varargin)],
   if (strcmp(lower(varargin{arg}),'limit')),
       limit = varargin{arg+1};
   elseif (strcmp(lower(varargin{arg}),'iterations')),
       maxit = varargin{arg+1};
   elseif (strcmp(lower(varargin{arg}),'d')),
       D = varargin{arg+1};
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
    rhs = fE(m - fE(rho_0,'notransp'), 'transp');
else
    rhs = fE(m,'transp');
end

%Calculate preconditioning matrix
if (isempty(D)),
    D = ones(size(rhs));
end

%Initial residual, we will have a starting guess of rho - rho_0 = 0.
r = D.*rhs;

rr_0 = r(:)'*r(:)
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
    q = D .* (fE(fE(D .* p,'notransp'),'transp') + lambda*fL(lambda*fL(D.*p,'notransp'),'transp'));
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


function b = L_default(x, direction)
b = 0; %default is no regularization
return