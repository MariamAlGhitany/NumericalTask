% Adams-Bashforth 2-step method
[t_ab2, y_ab2] = adamsBashforth2(f, tspan, y0, h);

% Adams-Moulton 2-step method
[t_am2, y_am2] = adamsMoulton2(f, tspan, y0, h);

% Plot results
figure;
plot(t_ab2, y_ab2, 'g-', 'LineWidth', 1.5); hold on;
plot(t_am2, y_am2, 'm--', 'LineWidth', 1.5);
xlabel('t');
ylabel('y(t)');
legend('Adams-Bashforth 2-step', 'Adams-Moulton 2-step');
title('Solution using Multi-Step Methods');
grid on;

% Adams-Bashforth 2-step Function
function [t, y] = adamsBashforth2(f, tspan, y0, h)
    t = tspan(1):h:tspan(2);
    y = zeros(size(t));
    y(1) = y0;

    % Use RK4 for the first step
    k1 = f(t(1), y(1));
    k2 = f(t(1) + h/2, y(1) + h/2*k1);
    k3 = f(t(1) + h/2, y(1) + h/2*k2);
    k4 = f(t(1) + h, y(1) + h*k3);
    y(2) = y(1) + h/6*(k1 + 2*k2 + 2*k3 + k4);

    for i = 2:length(t)-1
        y(i+1) = y(i) + h/2*(3*f(t(i), y(i)) - f(t(i-1), y(i-1)));
    end
end

% Adams-Moulton 2-step Function
function [t, y] = adamsMoulton2(f, tspan, y0, h)
    t = tspan(1):h:tspan(2);
    y = zeros(size(t));
    y(1) = y0;

    % Use RK4 for the first step
    k1 = f(t(1), y(1));
    k2 = f(t(1) + h/2, y(1) + h/2*k1);
    k3 = f(t(1) + h/2, y(1) + h/2*k2);
    k4 = f(t(1) + h, y(1) + h*k3);
    y(2) = y(1) + h/6*(k1 + 2*k2 + 2*k3 + k4);

    for i = 2:length(t)-1
        % Predictor (Adams-Bashforth)
        y_pred = y(i) + h/2*(3*f(t(i), y(i)) - f(t(i-1), y(i-1)));
        % Corrector (Adams-Moulton)
        y(i+1) = y(i) + h/2*(f(t(i+1), y_pred) + f(t(i), y(i)));
    end
end
