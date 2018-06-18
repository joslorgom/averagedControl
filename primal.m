function [zMT] = primal(u, A, B, z0, tout)

% primalApp receives the control vector and returns the primal solutions
% average in the parameter space at time t = T

% Number of time steps
Nt = length(tout);
% Size of parameter space
M = size(A, 3);

% Average of primal solutions
zM = zeros(Nt, 1);
% Solve primal problem for every parameter in nu
for j = 1:M
    % Update matrix A
    Am = A(:, :, j);
    % Update matrix B
    Bm = B(:, :, j);
    % Primal problem zdiff = A*z + B*u
    [tout, zout] = ode45(@(t, z) Am*z + Bm*interp1(tout, u, t), tout, z0);
    % Save final state for parameter nu(j)
    % zend(j, :) = zout(end, :);
    % Update average state
    zM = zM + zout; 
end
% Update average state
zM = zM/M;

zMT = zM(end, :)';
    
end
