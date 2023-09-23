clc,clear;
% Exogenous variables and Parameters
R=1.02;
Y=1;
ss=0.5;
beta=0.9;
gamma=0.5;

% Grid on state variable
num_a=200;
a_min=0.00001;
a_max=1;
a_grid=linspace(a_min,a_max,num_a)';

%  dimention of approximate polynominal
dim_app=2;

% 设置初始系数猜测值(维度必须是dim_app+1)
coef_ini = [0.1, 0.35, 0.01];
% 设置fsolve的优化选项：使用levenberg-marquardt算法，并设置最大函数求值次数为1000
options = optimoptions('fsolve', 'Algorithm', 'levenberg-marquardt', 'MaxFunctionEvaluations', 1000);
%定义resid_projection (将已知变量代入b_projection_resid函数)
resid=@(coef)b_projection_resid(coef,a_grid,beta,gamma,R,Y,ss,dim_app);

% 使用fsolve找到使欧拉方程误差(resid_projection)为零的系数
 coef = fsolve(resid, coef_ini, options);
 
% 使用fsolve得到的系数来计算近似的policy function (saving function)
 for i = 0:dim_app
    XX(:,i+1) = a_grid.^i;
 end
a_prime=XX*coef';

% 绘图部分
% saving function
% The result is a_prime_grid which maps every a_grid(i) to its corresponding a_prime value.
% Plotting the mapping between current assets (a) and next period assets (a_prime)
figure;
plot(a_grid, a_prime, '-b'); % Plot the values with a blue line
xlabel('Current Assets (a)');    % Label for the x-axis
ylabel('Next Period Assets (a'')'); % Label for the y-axis (Note the double '' to escape the single quote in the label)
title('Mapping of Current Assets to Next Period Assets'); % Title for the plot
grid on; % Display grid lines for better readability
% show a 45-degree line
hold on;
plot(a_grid, a_grid, '--r');
% Legend to differentiate between the two lines
legend('Policy Function', '45-degree Line'); 
hold off;

% consumption function
% Consumption (Cy) is determined by total resources (Y + a_grid) minus next period's assets (a_prime_grid)
Cy = Y + a_grid - a_prime;
% Plot the relationship between Consumption (Cy) and Current Assets (a_grid)
figure;
plot(a_grid,Cy);
ylabel('Consumption (Cy)');   % Label for the x-axis
xlabel('Current Assets (a)'); % Label for the y-axis
title('Relationship between Consumption and Current Assets'); % Title for the plot
grid on;

