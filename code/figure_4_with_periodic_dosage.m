%% Initialization

clc
clear
provided_params = load("provided_params.mat");


%% Figure 4 Simulations

% Modify Parameters
modified_params = provided_params;
modified_params.rho_C = 0.9;
modified_params.g_T = 10^10;
modified_params.alpha_1 = 0;
modified_params.g_C = 2*10^9;
modified_params.tau_C = 7;
modified_params.rho_T = 0.02;
modified_params.alpha_2 = 2.5*10^(-10);


modified_params.dosage_period = 50;
modified_params.dosage_amount = 2*10^8;

% Conditions
tspan = [0 240];
y_initial_values = [ % [C T]
    4*10^8 3.35*10^10
    8*10^7 6.7*10^9
    4*10^8 6.7*10^9
];

alpha_1_values = [0.04 0.07 0.1];

% Simulation
n = length(alpha_1_values);
second_simulation_time = containers.Map;
second_simulation_C = containers.Map;
second_simulation_T = containers.Map;

for initial_value_index = 1:length(y_initial_values)
    y_initial = y_initial_values(initial_value_index,:);
    
    temp_time = containers.Map;
    temp_C = containers.Map;
    temp_T = containers.Map;
    
    for alpha_1_index = 1:n
        modified_params.alpha_1 = alpha_1_values(alpha_1_index);
        model_with_alpha_1 = paper_model_first_periodic_dosage(modified_params);

        [T, Y] = solve_ode(model_with_alpha_1, tspan, y_initial);

        temp_time(num2str(alpha_1_values(alpha_1_index))) = T;
        temp_C(num2str(alpha_1_values(alpha_1_index))) = Y(:,1);
        temp_T(num2str(alpha_1_values(alpha_1_index))) = Y(:,2);
    end
    
    second_simulation_time(num2str(initial_value_index)) = temp_time;
    second_simulation_C(num2str(initial_value_index)) = temp_C;
    second_simulation_T(num2str(initial_value_index)) = temp_T;
end

%% Figure 4 Plotting

figure
tiledlayout(2,3)

% a
nexttile
title("Figure 4a");
hold on
ylim([0 2*10^9]);
time = second_simulation_time("1");
data = second_simulation_C("1");
plot(time("0.04"), data("0.04"), 'b', 'linewidth', 3)
% plot(time("0.07"), data("0.07"), 'b:', 'linewidth', 3)
plot(time("0.1"), data("0.1"), 'b--', 'linewidth', 3)
legend("0.04", "0.1")

% b
nexttile
title("Figure 4b");
hold on
ylim([0 10*10^8]);
time = second_simulation_time("2");
data = second_simulation_C("2");
plot(time("0.04"), data("0.04"), 'b', 'linewidth', 3)
% plot(time("0.07"), data("0.07"), 'b:', 'linewidth', 3)
plot(time("0.1"), data("0.1"), 'b--', 'linewidth', 3)

% c
nexttile
title("Figure 4c");
hold on
ylim([0 10*10^8]);
time = second_simulation_time("3");
data = second_simulation_C("3");
plot(time("0.04"), data("0.04"), 'b', 'linewidth', 3)
% plot(time("0.07"), data("0.07"), 'b:', 'linewidth', 3)
plot(time("0.1"), data("0.1"), 'b--', 'linewidth', 3)

% d
nexttile
title("Figure 4d");
hold on
ylim([0.001*10^10 10*10^10]);
time = second_simulation_time("1");
data = second_simulation_T("1");
plot(time("0.04"), data("0.04"), 'r', 'linewidth', 3)
% plot(time("0.07"), data("0.07"), 'r:', 'linewidth', 3)
plot(time("0.1"), data("0.1"), 'r--', 'linewidth', 3)
legend("0.04", "0.1")

% e
nexttile
title("Figure 4e");
hold on
ylim([0.001*10^10 2*10^10]);
time = second_simulation_time("2");
data = second_simulation_T("2");
plot(time("0.04"), data("0.04"), 'r', 'linewidth', 3)
% plot(time("0.07"), data("0.07"), 'r:', 'linewidth', 3)
plot(time("0.1"), data("0.1"), 'r--', 'linewidth', 3)

% f
nexttile
title("Figure 4f");
hold on
ylim([0.001*10^10 2*10^10]);
time = second_simulation_time("3");
data = second_simulation_T("3");
plot(time("0.04"), data("0.04"), 'r', 'linewidth', 3)
% plot(time("0.07"), data("0.07"), 'r:', 'linewidth', 3)
plot(time("0.1"), data("0.1"), 'r--', 'linewidth', 3)