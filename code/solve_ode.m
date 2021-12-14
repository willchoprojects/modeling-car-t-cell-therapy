function [T, Y] = solve_ode(func, t, y_0) 
    [T, Y] = ode15s(func, t, y_0);
end