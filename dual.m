function [pM] = dual(p0, A, B, tout)

% dualApp receives the adjoint state at time t = T and returns the dual 
% solutions average in the parameter space

% Number of time steps
Nt = length(tout);
% Size of parameter space
M = size(A, 3);

% Adjoint variable average: 1/N*sum_i(B*pout(t, nu(i)))
pM = zeros(Nt, 1);
% Solve adjoint problem forward in time for every parameter in nu
for j = 1:M   
    % Update matrix A
    Am = A(:, :, j);
    % Update matrix B
    Bm = B(:, :, j);
    % Solve adjoint problem forward in time
    [tout, pout] = ode45(@(t, p) Am'*p, tout, -p0); 
    pM = pM + pout*Bm;
end
pM = -pM/M;

% Reverse adjoint variable in time
pM = flipud(pM);

end
