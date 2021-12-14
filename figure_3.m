%% Initialization

clc
clear
provided_params = load("provided_params.mat");

%% Figure 3 Simulations

% Modify Parameters
modified_params = provided_params;
modified_params.rho_C = 0.9;
modified_params.g_T = 10^10;
modified_params.alpha_1 = 0.04;
modified_params.g_C = 2*10^9;
modified_params.tau_C = 7;
modified_params.rho_T = 0.02;
modified_params.alpha_2 = 2.5*10^(-10);

% Conditions
tspan = [0 60];
y_initial_values = [ % [C T]
    5*10^7 3.35*10^10
    8*10^7 3.35*10^10
    1.5*10^8 3.35*10^10
    2*10^8 3.35*10^10
    4*10^8 3.35*10^10
    7*10^8 3.35*10^10
]; 

% Simulation
n = length(y_initial_values);
first_simulation_time = containers.Map;
first_simulation_C = containers.Map;
first_simulation_T = containers.Map;

model = paper_model_first(modified_params);

for index = 1:n
    y_initial = y_initial_values(index,:);
    
    [T, Y] = solve_ode(model, tspan, y_initial);
    
    first_simulation_time(num2str(y_initial_values(index,1))) = T;
    first_simulation_C(num2str(y_initial_values(index,1))) = Y(:,1);
    first_simulation_T(num2str(y_initial_values(index,1))) = Y(:,2);
end

%% Figure 3 Plotting

figure
tiledlayout(2,2)

% a
nexttile
title("Figure 3a");
hold on
ylim([0 2*10^8]);
keys = [string(num2str(5*10^7)) string(num2str(8*10^7)) string(num2str(1.5*10^8))];
plot(first_simulation_time(keys(1)), first_simulation_C(keys(1)), 'b', 'linewidth', 3)
plot(first_simulation_time(keys(2)), first_simulation_C(keys(2)), 'b:', 'linewidth', 3)
plot(first_simulation_time(keys(3)), first_simulation_C(keys(3)), 'b--', 'linewidth', 3)
legend("5*10^7", "8*10^7", "1.5*10^8")

% b
nexttile
title("Figure 3b");
hold on
ylim([0 2.5*10^9]);
keys = [string(num2str(2*10^8)) string(num2str(4*10^8)) string(num2str(7*10^8))];
plot(first_simulation_time(keys(1)), first_simulation_C(keys(1)), 'b', 'linewidth', 3)
plot(first_simulation_time(keys(2)), first_simulation_C(keys(2)), 'b:', 'linewidth', 3)
plot(first_simulation_time(keys(3)), first_simulation_C(keys(3)), 'b--', 'linewidth', 3)
legend("2*10^8", "4*10^8", "7*10^8")

% c
nexttile
title("Figure 3c");
hold on
ylim([0 10*10^10]);
keys = [string(num2str(5*10^7)) string(num2str(8*10^7)) string(num2str(1.5*10^8))];
plot(first_simulation_time(keys(1)), first_simulation_T(keys(1)), 'r', 'linewidth', 3)
plot(first_simulation_time(keys(2)), first_simulation_T(keys(2)), 'r:', 'linewidth', 3)
plot(first_simulation_time(keys(3)), first_simulation_T(keys(3)), 'r--', 'linewidth', 3)
legend("5*10^7", "8*10^7", "1.5*10^8")

% d
nexttile
title("Figure 3d");
hold on
ylim([0 5*10^10]);
keys = [string(num2str(2*10^8)) string(num2str(4*10^8)) string(num2str(7*10^8))];
plot(first_simulation_time(keys(1)), first_simulation_T(keys(1)), 'r', 'linewidth', 3)
plot(first_simulation_time(keys(2)), first_simulation_T(keys(2)), 'r:', 'linewidth', 3)
plot(first_simulation_time(keys(3)), first_simulation_T(keys(3)), 'r--', 'linewidth', 3)
legend("2*10^8", "4*10^8", "7*10^8")