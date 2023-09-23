clc,clear;
% Exogenous variables and Parameters
R=1.02;
Y=1;
ss=0.5;
beta=0.9;
gamma=0.5;

% Grid on state variable in t+1
num_a=2000;
a_prime_min=0.00001;
a_prime_max=0.999;
a_prime_grid=linspace(a_prime_min,a_prime_max,num_a)';

%consumption and cash-on-hand
Cy=(beta*R*(ss+R*a_prime_grid).^(-gamma)).^(-1/gamma);
a=Cy+a_prime_grid-Y;

% consumption function (policy function)
figure;
plot(a,Cy);
ylabel('Consumption (Cy)');
xlabel('Current Assets (a)');
title('Relationship between Consumption and Current Assets');
grid on;

% saving function (policy function)
figure;
plot(a, a_prime_grid, '-b');
xlabel('Current Assets (a)');
ylabel('Next Period Assets (a'')'); 
title('Mapping of Current Assets to Next Period Assets'); 
grid on;
% show a 45-degree line
hold on;
plot(a, a, '--r');
legend('Policy Function', '45-degree Line'); 
hold off;

