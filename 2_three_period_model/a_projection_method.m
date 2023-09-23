clc,clear;
% Exogenous variables and parameters
ex.YY=1;
ex.YM=1;
ex.R=1.02;
ex.ss=0.5;
ex.beta=0.9;
ex.gamma=0.5;

%% Grid on state variable in period M
num_aM=200;
aM_min=0.00001;
aM_max=2;
aM_grid=linspace(aM_min,aM_max,num_aM)';
%  dimention of approximate polynominal
dim_app=2;

%% Solve Euler quation between old and middle periods
old_period=1;
% initial value for coefficient (the dimension must be dim_app+1)
coef_ini = [0.1, 0.35, 0.01];
options = optimoptions('fsolve', 'Algorithm', 'levenberg-marquardt', 'MaxFunctionEvaluations', 1000);
%define resid_projection (substitude known variables to a_projection_resid function)
resid=@(coef)a_projection_resid(coef,aM_grid,ex,dim_app,old_period,coef_ini);
% find coefficients that satisfy the zero point of Euler residual
 coef_aO = fsolve(resid, coef_ini, options);

%% Grid on state variable in period Y
num_aY=200;
aY_min=0.00001;
aY_max=2;
aY_grid=linspace(aY_min,aY_max,num_aY)';

%% Solve Euler quation between middle and young
old_period=0;
%define resid_projection (substitude known variables to a_projection_resid function)
resid=@(coef)a_projection_resid(coef,aM_grid,ex,dim_app,old_period,coef_aO);
% find coefficients that satisfy the zero point of Euler residual
 coef_aM = fsolve(resid, coef_ini, options);

%% plot the policy function
%chebyshev polynomial for policy functions
aO_fcn = chebyshev_poly(dim_app,aM_grid,coef_aO);
aM_fcn = chebyshev_poly(dim_app,aY_grid,coef_aM);

figure;
plot(aM_grid, aO_fcn, 'r', 'LineWidth', 2); hold on; % Plotting aO_fcn vs aM_grid in red
plot(aY_grid, aM_fcn, 'b', 'LineWidth', 2); % Plotting aM_fcn vs aY_grid in blue
plot(aM_grid, aM_grid, '--k'); % Plotting the 45-degree line
xlabel('current Assets', 'FontSize', 12);
ylabel('next Assets', 'FontSize', 12);
title('Policy Function Plot', 'FontSize', 14);
legend('Old Period Policy Function', 'Middle Period Policy Function', '45-degree line', 'Location', 'northwest');
grid on;
