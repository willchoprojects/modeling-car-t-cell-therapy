%% Initialization

clc
clear
provided_params = load("provided_params.mat");

%% Figure 1 Simulations

% Modify Parameters
modified_params = provided_params;
modified_params.rho_C = 0.9;
modified_params.g_T = 10^10;
modified_params.alpha_1 = 0;
modified_params.g_C = 2*10^9;
modified_params.tau_C = 7;
modified_params.rho_T = 0.02;
modified_params.alpha_2 = 2.5*10^(-10);

% Conditions
tspan = [0 60];
y_initial = [8*10^7 3.35*10^10]; % [C T]

alpha_1_values = [0 0.01 0.02 0.03  0.035 0.04 0.1];

% Simulation
n = length(alpha_1_values);
first_simulation_time = containers.Map;
first_simulation_C = containers.Map;
first_simulation_T = containers.Map;

for index = 1:n
    modified_params.alpha_1 = alpha_1_values(index);
    model_with_alpha_1 = paper_model_first(modified_params);
    
    [T, Y] = solve_ode(model_with_alpha_1, tspan, y_initial);
    
    first_simulation_time(num2str(alpha_1_values(index))) = T;
    first_simulation_C(num2str(alpha_1_values(index))) = Y(:,1);
    first_simulation_T(num2str(alpha_1_values(index))) = Y(:,2);
end

%% Figure 1 Plotting

figure
tiledlayout(2,2)

% a
nexttile
title("Figure 1a");
hold on
plot(first_simulation_time("0"), first_simulation_C("0"), 'b', 'linewidth', 3)
plot(first_simulation_time("0.02"), first_simulation_C("0.02"), 'b--', 'linewidth', 3)
plot(first_simulation_time("0.03"), first_simulation_C("0.03"), 'b-.', 'linewidth', 3)
legend("0", "0.02", "0.03")

% b
nexttile
title("Figure 1b");
hold on
plot(first_simulation_time("0.035"), first_simulation_C("0.035"), 'b--', 'linewidth', 3)
plot(first_simulation_time("0.04"), first_simulation_C("0.04"), 'b-.', 'linewidth', 3)
plot(first_simulation_time("0.1"), first_simulation_C("0.1"), 'b:', 'linewidth', 3)
legend("0.035", "0.04", "0.1")

% c
nexttile
title("Figure 1c");
hold on
plot(first_simulation_time("0"), first_simulation_T("0"), 'r', 'linewidth', 3)
plot(first_simulation_time("0.02"), first_simulation_T("0.02"), 'r--', 'linewidth', 3)
plot(first_simulation_time("0.03"), first_simulation_T("0.03"), 'r-.', 'linewidth', 3)
legend("0", "0.02", "0.03")

% d
nexttile
title("Figure 1d");
hold on
plot(first_simulation_time("0.035"), first_simulation_T("0.035"), 'r--', 'linewidth', 3)
plot(first_simulation_time("0.04"), first_simulation_T("0.04"), 'r-.', 'linewidth', 3)
plot(first_simulation_time("0.1"), first_simulation_T("0.1"), 'r:', 'linewidth', 3)
legend("0.035", "0.04", "0.1")