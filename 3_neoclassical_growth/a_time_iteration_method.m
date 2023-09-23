clc,clear;
%% parameters
st.gamma=0.5;
st.delta=1;
st.beta=0.96;
st.alpha=0.4;

%% grid on state variable-- capital k
y_ss=(1/st.beta-1+st.delta)/st.alpha;
k_ss=y_ss^(1/(st.alpha-1));

min_k=0.05;
max_k=0.5;
num_k=200;
k_grid=linspace(min_k,max_k,num_k)';

%% time iteration method

% initial policy function for consumption
c_fcn_0=k_grid;
c_fcn_1=zeros(num_k,1);

% setting for iteration
maxiter = 20; % maximal iteration 
tol  = 1.0e-5;  % 收敛容忍度
options = optimoptions('fsolve','Display','none');

it=1; %iteration counter
diff=1000; % initial value of policy function difference
diff_list = zeros(maxiter,1); 

%substitude known variables into Euler_residual function
while (diff > tol && it < maxiter)
    % find the policy function for consumption
    % to satisfy zero Euler residual
    
    for ii=1:num_k
    capital=k_grid(ii);
    residual=@(c)a_Euler_residual(c,k_grid,capital,st,c_fcn_0);
    c_fcn_1(ii)=fsolve(residual,c_fcn_0(ii),options);
    end
    
    % difference of current and next policy function
    diff = max(abs(c_fcn_1-c_fcn_0));
    diff_list(it)= diff;
    
    % renew of policy function for consumption
    c_fcn_0 = c_fcn_1;
    it=it+1;
    
    % print information in iteration
    fprintf('iteration index: %i \n', it-1);
    fprintf('policy function iteration error: %e\n', diff);
end

c_policy=c_fcn_0;
k_next_policy=k_grid.^st.alpha + (1-st.delta)*k_grid - c_policy;

figure;
plot(k_grid,c_policy);

% Plotting the policy function of next period capital
figure;
plot(k_grid, k_next_policy, '-b'); 
xlabel('Current Assets (k)');
ylabel('Next Period Assets ($k_{t+1}$)');
title('Mapping of Current capital to Next Period capital');
grid on; % Display grid lines for better readability
% show a 45-degree line
hold on;
plot(k_grid,k_grid, '--r');
% Legend to differentiate between the two lines
legend('Policy Function', '45-degree Line'); 
hold off;
