tic;

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
N = 3;

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
Nt = 100;

% Maximum number of iterations
Nmax = 100;

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

% Iteration counter
iter = 0;

% Compute b in gradient = A*u-b for CG method
% Free dynamics solution
y = primal(u, A, B, z0, tout);
b = dual(zt - y, A, B, tout);

% Compute A*u
% Zero initial condition solution
zT = primal(u, A, B, zeros(N, 1), tout);
Au = dual(zT, A, B, tout) + beta*u;

% Compute initial gradient and store previous gradient
ga = Au - b; 
g = ga;

% Initial gradient norm
g0L2 = integral(@(t) interp1(tout, g, t).^2, T0, T);
gaL2 = g0L2;
gL2 = g0L2;

% Residual b-A*u
r = -g;
% Residual norm
rn = sqrt( integral(@(t) interp1(tout, r, t).^2, T0, T) );

while (rn > tol && iter <= Nmax)
    
    % A*u
    zT = primal(r, A, B, zeros(N, 1), tout);
    w = dual(zT, A, B, tout) + beta*r;
    
    % alpha
    alpha = integral(@(t) interp1(tout, g, t).^2, T0, T)/...
        integral(@(t) interp1(tout, r, t).*interp1(tout, w, t), T0, T);

    % Update control
    u = u + alpha*r;
    
    % Store previous gradient
    ga = g;
    gaL2 = gL2;
    % Update gradient
    g = g + alpha*w;
    % Compute gradient norm
    gL2 = integral(@(t) interp1(tout, g, t).^2, T0, T);
    
    % gamma
    gamma = gL2/gaL2;
  
    % Update residual
    r = -g + gamma*r;
    
    % Compute residual norm
    rn = sqrt( integral(@(t) interp1(tout, r, t).^2, T0, T) );
    
    % Update iteration counter
    iter = iter + 1;
    
    % Compute functional
    Jmu = 0.5 * (zT - zt)' * (zT - zt);
    Ju = 0.5 * beta * integral(@(t) interp1(tout, u, t).^2, T0, T);
    J = Jmu + Ju;

    % Print status message
    fprintf("Iteration %i - Error %g - Cost %g\n", iter, rn, J);
    
end


% Solve state average for the optimal control
% Output average
zM = zeros(Nt, N);
for j = 1:M 
    % Update matrix A
    Am = A(:, :, j);
    % Update matrix B
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
        % Update matrix A
        Am = A(:, :, j);
        % Update matrix B
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
        % Update matrix A
        Am = A(:, :, j);
        % Update matrix B
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

toc;