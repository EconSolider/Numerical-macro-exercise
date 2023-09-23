clc,clear;
% Exogenous variables and Parameters
R=1.02;
Y=1;
ss=0.5;

beta=0.9;
gamma=0.5;

% Grid on state variable
num_a=2000;
a_min=0.00001;
a_max=1;
a_grid=linspace(a_min,a_max,num_a)';

%Define the Euler quation residual
resid=@(a,a_prime) beta*R* ((ss+R*a_prime)/(Y+a-a_prime))^(-gamma) - 1;

% The next part solve the Eulter equation to 
    % get the saving function a_prime=g(a)
a_prime_grid=zeros(num_a,1);
for i = 1:num_a
    a = a_grid(i);
    % Anonymous function for the specific value of 'a'
    func_for_a = @(a_prime) resid(a, a_prime);
    % Starting guess can be the same value of 'a', though different heuristics might be more efficient
    a_prime_guess = a;
    a_prime_grid(i) = fzero(func_for_a, a_prime_guess);  
end

% saving function
% The result is a_prime_grid which maps every a_grid(i) to its corresponding a_prime value.
% Plotting the mapping between current assets (a) and next period assets (a_prime)
figure;
plot(a_grid, a_prime_grid, '-b');
xlabel('Current Assets (a)');  
ylabel('Next Period Assets (a'')');
title('Mapping of Current Assets to Next Period Assets'); 
grid on; % Display grid lines for better readability
% show a 45-degree line
hold on;
plot(a_grid, a_grid, '--r');
% Legend to differentiate between the two lines
legend('Policy Function', '45-degree Line'); 
hold off;

% consumption function
% Consumption (Cy) is determined by total resources (Y + a_grid) minus next period's assets (a_prime_grid)
Cy = Y + a_grid - a_prime_grid;
% Plot the relationship between Consumption (Cy) and Current Assets (a_grid)
figure;
plot(a_grid,Cy);
ylabel('Consumption (Cy)');   % Label for the x-axis
xlabel('Current Assets (a)'); % Label for the y-axis
title('Relationship between Consumption and Current Assets'); % Title for the plot
grid on;


