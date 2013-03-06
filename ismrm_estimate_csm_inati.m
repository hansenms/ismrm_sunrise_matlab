function [sensitivities] = ismrm_estimate_csm_inati(data, method)
% Adaptive combination weights estimation
%
% [sensitivities] = ismrm_estimate_csm_inati(data, method)
%
% Input:
%   data: [Nx,Ny,Nc]
%   2D complex images [Nx,Ny] from Nc coils
%   assumed to be prewhitened
%
%   method: integer
%   0: proposed method
%   1: Walsh et al., MRM 2000;43(5)
%
% Output:
%   sensitivities: [Nx,Ny,Nc]
%   Estimated relative coil sensitivities
%
%  TODO: Add description and reference
%  

if nargin < 2,
    method = 0;
end

[Nx,Ny,Nc] = size(data);

% Initialize the result
sensitivities = zeros(Nx,Ny,Nc);

% Loop over the pixels in the image
% We ignore the edges here for simplicity
% A full implementaion would handle the edges properly.
for x = 4:(Nx-3)
  for y = 4:(Ny-3)

    % Get the data in a box around the current point
    % Use a 7x7 box.  Size is arbitrary, i.e. coil smoothness assumption
    D = reshape(data(x-3:x+3,y-3:y+3,:),[7*7,Nc]);

    % Compute the sensitivity
    if (method == 0)
      % Proposed method

      % The SVD way
      % [U,S,V]   = svd(D);
      % sorted from biggest to smallest
      % u1 = U(:,1); % first left singular vector (Bx*By,1)
      % v1 = V(:,1); % first right singular vector (Nc,1)
      % s1 = S(1,1); % first singular value
      % the above is slow, and since we only care about the first
      % singular vector, we can use the power method to solve
      
      % Power method
      % initialize to the mean of the data
      v1 = transpose(mean(D,1)); 
      v1 = v1/norm(v1);
      % 3 iterations
      for iter = 1:3
        v1 = D'*D*v1; 
        v1 = v1/norm(v1);
      end
      % compute u1 and s1
      u1 = D*v1;
      s1 = norm(u1);
      u1 = u1/s1;
    
      % Compute the phase of the average of u1
      theta = angle(mean(u1));

      % Normalized coil sensitivies (combination weights)
      sensitivities(x,y,:)  = exp(1i*theta)*v1';
      
    else
      % Original Walsh et al. method

      % Compute the covariance over the box (an Nc by Nc matrix)
      R = D'*D;
      % The eigenvalue way
      % [V,D]   = eig(R);
      % sorted from smallest to biggest
      % v1 = V(:,end); % biggest eigenvector (Nc,1)
      % the above is slow, and since we only care about the first
      % eigenvector, we can use the power method to solve
      
      % Power method
      % initialize to the mean of the covariance
      v1 = transpose(mean(R,1)); 
      v1 = v1/norm(v1);
      % 3 iterations
      for iter = 1:3
        v1 = R*v1; 
        v1 = v1/norm(v1);
      end
      
      % Compute the phase of the reference coil at this location
      reference_coil = 1;  % arbitrary
      theta = angle(v1(reference_coil));
      
      % Normalized relative coil sensitivies
      sensitivities(x,y,:)  = exp(1i*theta)*v1';
      
    end

  end % Loop over y
end % Loop over x

end
