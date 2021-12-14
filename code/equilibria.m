%% Initialization

clc
clear
provided_params = load("provided_params.mat");

%% Parameters

modified_params = provided_params;
modified_params.rho_C = 0.9;
modified_params.g_T = 10^10;
modified_params.alpha_1 = 0.04;
modified_params.g_C = 2*10^9;
modified_params.tau_C = 7;
modified_params.rho_T = 0.02;
modified_params.alpha_2 = 2.5*10^(-10);

params = modified_params;
rho_C  = params.rho_C;
g_T = params.g_T;
alpha_1 = params.alpha_1;
g_C = params.g_C;
tau_C = params.tau_C;
rho_T = params.rho_T;
alpha_2 = params.alpha_2;

%% Jacobian Definition

syms C T;
J = matlabFunction(jacobian([
    (rho_C * C * T)/(g_T + T) - (alpha_1 * C * T)./(g_C + C) - (1 / tau_C) * C
    rho_T * T - alpha_2 * C * T
]));

%% E_1

C_value = 0;
T_value = 0;

fprintf("E_1:\n\n");
disp([C_value T_value]);

fprintf("Expected Jacobian:\n\n");
disp([
    -1/tau_C 0
    0 rho_T
]);

fprintf("Calculated Jacobian:\n\n");
calculated_jacobian = J(C_value, T_value);
disp(calculated_jacobian);

fprintf("Eigenvalues:\n\n");
[~,D] = eig(calculated_jacobian);
disp(D);

fprintf("Trace:\n\n");
disp(trace(calculated_jacobian));

fprintf("Determinant:\n\n");
disp(det(calculated_jacobian));

fprintf("tr^2 - 4 det:\n\n");
disp(trace(calculated_jacobian)^2 - 4 * det(calculated_jacobian));

%% E_2 and E_3

alpha = alpha_1 / (g_C + rho_T / alpha_2);
discriminant = (alpha * g_T - rho_C + 1/tau_C)^2 - 4 * alpha * g_T / tau_C;
T_star_1 = (-(alpha * g_T - rho_C + 1/tau_C) + sqrt(discriminant))/(2 * alpha);
T_star_2 = (-(alpha * g_T - rho_C + 1/tau_C) - sqrt(discriminant))/(2 * alpha);

fprintf("Discriminant:\n\n");
disp(discriminant);

C_value = rho_T/alpha_2;
T_value = T_star_1;





fprintf("E_2:\n\n");
disp([C_value T_value]);

fprintf("Calculated Jacobian:\n\n");
calculated_jacobian = J(C_value, T_value);
disp(calculated_jacobian);

fprintf("Eigenvalues:\n\n");
[~,D] = eig(calculated_jacobian);
disp(D);

fprintf("Trace:\n\n");
disp(trace(calculated_jacobian));

fprintf("Determinant:\n\n");
disp(det(calculated_jacobian));

fprintf("tr^2 - 4 det:\n\n");
disp(trace(calculated_jacobian)^2 - 4 * det(calculated_jacobian));





C_value = rho_T/alpha_2;
T_value = T_star_2;

fprintf("E_3:\n\n");
disp([C_value T_value]);

fprintf("Calculated Jacobian:\n\n");
calculated_jacobian = J(C_value, T_value);
disp(calculated_jacobian);

fprintf("Eigenvalues:\n\n");
[~,D] = eig(calculated_jacobian);
disp(D);

fprintf("Trace:\n\n");
disp(trace(calculated_jacobian));

fprintf("Determinant:\n\n");
disp(det(calculated_jacobian));

fprintf("tr^2 - 4 det:\n\n");
disp(trace(calculated_jacobian)^2 - 4 * det(calculated_jacobian));

%% Phase Plane

c_span = linspace(0, 3*10^9, 20);
t_span = linspace(0, 8*10^10, 20);
[c,t] = meshgrid(c_span, t_span);
f = paper_model_first(params);

t_0 = 0;
for index = 1:(length(c) * length(t))
    y_prime = f(t_0,[c(index); t(index)]);
    dCdt(index) = y_prime(1);
    dTdt(index) = y_prime(2);
end

figure 
hold on

quiver(c, t, reshape(dCdt, size(c)), reshape(dTdt, size(t)), 'r');

time_span = 0:300;
initial_conditions = [
    0.05*10^9 0.5*10^10
    0.25*10^9 0.5*10^10
    0.50*10^9 0.5*10^10
    0.75*10^9 0.5*10^10
    1.00*10^9 0.5*10^10
    0.05*10^9 1*10^10
    0.25*10^9 1*10^10
    0.50*10^9 1*10^10
    0.75*10^9 1*10^10
    1.00*10^9 1*10^10
    0.05*10^9 2*10^10
    0.25*10^9 2*10^10
    0.50*10^9 2*10^10
    0.75*10^9 2*10^10
    1.00*10^9 2*10^10
    0.05*10^9 3*10^10
    0.25*10^9 3*10^10
    0.50*10^9 3*10^10
    0.75*10^9 3*10^10
    1.00*10^9 3*10^10
    0.05*10^9 4*10^10
    0.25*10^9 4*10^10
    0.50*10^9 4*10^10
    0.75*10^9 4*10^10
    1.00*10^9 4*10^10
    0.05*10^9 5*10^10
    0.25*10^9 5*10^10
    0.50*10^9 5*10^10
    0.75*10^9 5*10^10
    1.00*10^9 5*10^10
];
for index = 1:length(initial_conditions)
    [T,Y] = solve_ode(f, time_span, initial_conditions(index, :));
    plot(Y(:,1),Y(:,2),'b');
end

plot(initial_conditions(:,1), initial_conditions(:,2), 'bo');

title("Figure A1");

xlabel('C');
ylabel('T');

axis tight;

xlim([0 max(c_span)]);
ylim([0 max(t_span)]);

%% Phase Plane with Equilibrium Points

figure 
title("Equilibrium Points");
hold on

quiver(c, t, reshape(dCdt, size(c)), reshape(dTdt, size(t)), 'r');

time_span = 0:240;
initial_conditions = [
    % E1
    0 0
    0 max(t_span)/100
    max(c_span)/100 0
    % E2, E3
    rho_T/alpha_2 T_star_1
    rho_T/alpha_2 T_star_2
    (rho_T/alpha_2 + max(c_span)/100) T_star_1
    (rho_T/alpha_2 + max(c_span)/100) T_star_2
    (rho_T/alpha_2 - max(c_span)/100) T_star_1
    (rho_T/alpha_2 - max(c_span)/100) T_star_2
    rho_T/alpha_2 (T_star_1 + max(t_span)/100)
    rho_T/alpha_2 (T_star_2 + max(t_span)/100)
    (rho_T/alpha_2 + max(c_span)/100) (T_star_1 + max(t_span)/100)
    (rho_T/alpha_2 + max(c_span)/100) (T_star_2 + max(t_span)/100)
    (rho_T/alpha_2 - max(c_span)/100) (T_star_1 + max(t_span)/100)
    (rho_T/alpha_2 - max(c_span)/100) (T_star_2 + max(t_span)/100)
    rho_T/alpha_2 (T_star_1 - max(t_span)/100)
    rho_T/alpha_2 (T_star_2 - max(t_span)/100)
    (rho_T/alpha_2 + max(c_span)/100) (T_star_1 - max(t_span)/100)
    (rho_T/alpha_2 + max(c_span)/100) (T_star_2 - max(t_span)/100)
    (rho_T/alpha_2 - max(c_span)/100) (T_star_1 - max(t_span)/100)
    (rho_T/alpha_2 - max(c_span)/100) (T_star_2 - max(t_span)/100)
];
for index = 1:length(initial_conditions)
    [T,Y] = solve_ode(f, time_span, initial_conditions(index, :));
    plot(Y(:,1),Y(:,2),'b');
end

plot(initial_conditions(:,1), initial_conditions(:,2), 'bo');

xlabel('C');
ylabel('T');

axis tight;

xlim([0 max(c_span)]);
ylim([0 max(t_span)]);

%% Alpha 1 Value for 2 Equilibrium Points Only

alpha_1_eq = (sqrt(rho_C)-sqrt(1/tau_C))^2*(g_C+rho_T/alpha_2)/g_T;