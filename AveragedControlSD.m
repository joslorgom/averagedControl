clc;
close all;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the optimal control for the cost function
% J = beta/2*integral(u^2) + 1/2*|| sum_i( zout(T, nu(i)) ) - ztarget ||^2
% subject to:
% dz/dt = A(nu(i))*z + B(nu(i))*u
% z(0) = z0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vector of parameters
nu = 1:0.5:6;

% Size of state vector
N = 2;

% Initial condition
z0 = ones(N, 1);

% Target at t = T
zt = zeros(N, 1); 

% Parameter beta for the cost function
beta = 1e-3;

% Initial time
T0 = 0;

% Final time
T = 1;

% Number of time steps
Nt = 50;

% Maximum number of iterations
Nmax = 20000000;

% Gradient step
% u = u + d*Du
d = 1;

% Tolerance
tol = 1e-8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System matrices A, B
% zdiff = A*z + B*u
% dz/dt = A*z + B*u
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Length of vector nu (number of parameters)
M = length(nu);

Am = eye(N, N);
for i = 1:N
    for j = 1:N
        if j > i
            Am(i, j) = 1;
        end
    end
end
Am = -Am;

Bm = zeros(N, 1);
Bm(N) = 1;

A = zeros(N, N, M);
B = zeros(N, 1, M);
for j = 1:M
    A(:, :, j) = Am + (nu(j) - 1 )*diag(diag(Am));
    B(:, :, j) = Bm;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time vector
tout = linspace(T0, T, Nt);

% Control at time k
u = zeros(Nt, 1);
% Control at time k-1
ua = u;
% Control update
Du = zeros(Nt, 1);

% Initial error
error = 10;
% Iteration counter
iter = 0;

while (error > tol && iter < Nmax)
    % Update iteration counter
    iter = iter + 1;
    
    % Output average: 1/N*sum_i(zout(t, nu(i)))
    zM = zeros(Nt, N);
    
    % Final state
    zend = zeros(M, N);

    % Solve primal problem for every parameter in nu
    for j = 1:M
        % Update matrix A
        Am = A(:, :, j);
        % Update matrix B
        Bm = B(:, :, j);
        % Primal problem zdiff = A*z + B*u
        [tout, zout] = ode45(@(t, z) Am*z + Bm*interp1(tout, u, t), tout, z0);
        % Save final state for parameter nu(j)
        zend(j, :) = zout(end, :);
        % Update average state
        zM = zM + zout; 
    end
    % Update average state
    zM = zM/M;
    
    % Initial condition for the adjoint problem
    p0 = - (zM(end, :)' - zt);
    
    % Adjoint variable average: 1/N*sum_i(B*pout(t, nu(i)))
    pM = zeros(Nt, 1);
    for j = 1:M   
        % Update matrix A
        Am = A(:, :, j);
        % Update matrix B
        Bm = B(:, :, j);
        % Solve adjoint problem forward in time
        [tout, pout] = ode45(@(t, p) Am'*p, tout, p0); 
        pM = pM + pout*Bm;
    end
    pM = pM/M;

    % Save previous control
    % u(k-1) = u(k)
    ua = u;
    % Reverse adjoint variable
    pM = flipud(pM);
    
    % Control update
    Du = beta*u - pM;
    % Update control
    u = u - d*Du;
    
    % Control update norm
    DU2 = integral(@(t) interp1(tout, Du, t).^2, T0, T);
    U2 = integral(@(t) interp1(tout, u, t).^2, T0, T);
    
    error = sqrt(DU2/U2);

    Jmu = 0.5 * (zM(end, :)' - zt)' * (zM(end, :)' - zt);
    Ju = 0.5 * beta * integral(@(t) interp1(tout, u, t).^2, T0, T);

    J = Jmu + Ju;

    fprintf("Iteration %i - Error %g - Cost %g\n", iter, error, J);
    
end

% Solve state average for the optimal control
% Output average
zM = zeros(Nt, N);
for j = 1:M 
    Am = A(:, :, j);
    Bm = B(:, :, j);
    % Primal problem zdiff = A*z + B*u
    [tout, zout] = ode45(@(t, z) Am*z + Bm*interp1(tout, u, t), tout, z0); 
    zM = zM + zout; 
end
zM = zM/M;

if N == 2
    figure(1)
    plot(tout, zM(:, 1), 'r', 'LineWidth', 2)
    hold on
    plot(tout, zM(:, 2), 'b', 'LineWidth', 2)
    title("Average State")
    legend("X_{av,1}", "X_{av,2}")
    xlabel("t")
    ylabel("x")
    set(gca,'FontSize',14)
    
    figure(2)
    plot(tout, u, 'g', 'LineWidth', 2)
    title("Control")
    xlabel("t")
    ylabel("u")
    set(gca,'FontSize',14)

    figure(3)
    for j = 1:M 
        Am = A(:, :, j);
        Bm = B(:, :, j);
        % Primal problem zdiff = A*z + B*u
        [tout, zout] = ode45(@(t, z) Am*z + Bm*interp1(tout, u, t), tout, z0);
        subplot(1, 2, 1)
        plot(tout, zout(:, 1))
        hold on
        subplot(1, 2, 2)
        plot(tout, zout(:, 2))
        hold on
    end 
    figure(3)
    subplot(1, 2, 1)
    xlabel("t")
    ylabel("x")
    title("X_{i,1}(t)")
    set(gca,'FontSize',14)
    subplot(1, 2, 2)
    xlabel("t")
    ylabel("x")
    title("X_{i,2}(t)")
    set(gca,'FontSize',14)
end

if N == 3
    figure(1)
    plot(tout, zM(:, 1), 'r', 'LineWidth', 2)
    hold on
    plot(tout, zM(:, 2), 'b', 'LineWidth', 2)
    plot(tout, zM(:, 3), 'g', 'LineWidth', 2)
    title("Average State")
    legend("X_{av,1}", "X_{av,2}", "X_{av,3}")
    xlabel("t")
    ylabel("x")
    set(gca,'FontSize',14)
    
    figure(2)
    plot(tout, u, 'g', 'LineWidth', 2)
    title("Control")
    xlabel("t")
    ylabel("u")
    set(gca,'FontSize',14)

    figure(3)
    for j = 1:M 
        Am = A(:, :, j);
        Bm = B(:, :, j);
        % Primal problem zdiff = A*z + B*u
        [tout, zout] = ode45(@(t, z) Am*z + Bm*interp1(tout, u, t), tout, z0);
        subplot(1, 3, 1)
        plot(tout, zout(:, 1))
        hold on
        subplot(1, 3, 2)
        plot(tout, zout(:, 2))
        hold on
        subplot(1, 3, 3)
        plot(tout, zout(:, 3))
        hold on
    end 
    figure(3)
    subplot(1, 3, 1)
    xlabel("t")
    ylabel("x")
    title("X_{i,1}(t)")
    set(gca,'FontSize',14)
    subplot(1, 3, 2)
    xlabel("t")
    ylabel("x")
    title("X_{i,2}(t)")
    set(gca,'FontSize',14)
    subplot(1, 3, 3)
    xlabel("t")
    ylabel("x")
    title("X_{i,3}(t)")
    set(gca,'FontSize',14)
end