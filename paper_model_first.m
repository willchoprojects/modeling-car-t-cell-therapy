function model_function = paper_model_first(params)
    model_function = @(t, y) paper_model_first_parameterized(t, y, params);
end

function dydt = paper_model_first_parameterized(~, y, params)
    rho_C  = params.rho_C;
    g_T = params.g_T;
    alpha_1 = params.alpha_1;
    g_C = params.g_C;
    tau_C = params.tau_C;
    rho_T = params.rho_T;
    alpha_2 = params.alpha_2;
    
    dydt = zeros(2,1);
    
    % dCdt
    dydt(1) = (rho_C .* y(1) .* y(2))./(g_T + y(2)) - (alpha_1 .* y(1) .* y(2))./(g_C + y(1)) - (1 / tau_C) .* y(1);
    % dTdt
    dydt(2) = rho_T * y(2) - alpha_2 * y(1) * y(2);
end