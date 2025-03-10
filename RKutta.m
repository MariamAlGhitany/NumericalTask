% Define the ODE
f = @(t, y) -50*y + sin(t);

% Initial condition
y0 = 1;

% Time span
tspan = [0, 1]; % Adjust as needed

% Step size for RK4
h = 0.02; % Adjust as needed

% Classical 4th-order Runge-Kutta (RK4)
[t_rk4, y_rk4] = rungeKutta4(f, tspan, y0, h);

% Adaptive Runge-Kutta (ode45)
[t_ode45, y_ode45] = ode45(f, tspan, y0);

% Plot results
figure;
plot(t_rk4, y_rk4, 'b-', 'LineWidth', 1.5); hold on;
plot(t_ode45, y_ode45, 'r--', 'LineWidth', 1.5);
xlabel('t');
ylabel('y(t)');
legend('RK4', 'ode45');
title('Solution using Runge-Kutta Methods');
grid on;

% RK4 Function
function [t, y] = rungeKutta4(f, tspan, y0, h)
    t = tspan(1):h:tspan(2);
    y = zeros(size(t));
    y(1) = y0;

    for i = 1:length(t)-1
        k1 = f(t(i), y(i));
        k2 = f(t(i) + h/2, y(i) + h/2*k1);
        k3 = f(t(i) + h/2, y(i) + h/2*k2);
        k4 = f(t(i) + h, y(i) + h*k3);
        y(i+1) = y(i) + h/6*(k1 + 2*k2 + 2*k3 + k4);
    end
end
