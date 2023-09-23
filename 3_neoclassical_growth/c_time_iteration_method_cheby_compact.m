clc,clear;
%% parameters
st.gamma=2;
st.delta=0.025;
st.beta=0.96;
st.alpha=0.4;

%% grid on state variable-- capital k
y_ss=(1/st.beta-1+st.delta)/st.alpha;
k_ss=y_ss^(1/(st.alpha-1));

min_k=0.5*k_ss;
max_k=1.5*k_ss;
num_k=100;
k_grid=linspace(min_k,max_k,num_k)';

%% time iteration method
% setting for iteration
maxiter = 40; % maximal iteration 
tol  = 1.0e-5;  % 收敛容忍度
% option for minimal searcher
options = optimoptions('fsolve', 'Algorithm', 'levenberg-marquardt', 'MaxFunctionEvaluations', 1000);

it=1; %iteration counter
diff=1000; % initial value of policy function difference
diff_list = zeros(maxiter,1); 

% setting for Chebyshev polynomial approximation
dim_app=1;
coef_0=[0.1,0.3]; % dimension of coef_0 should be (dim_app+1)

% Chebyshev bases
    % this is for the estimation of Chebyshev coefficients
    % and Chebyshev interpolation
num_k=length(k_grid);
XX=zeros(num_k,dim_app);
for ii=0:dim_app
   XX(:,ii+1)=chebyshevT(ii,k_grid);
end

while (diff > tol && it < maxiter)
    % determine the Chebyshev approximate coefficients to meet the
        % zero point of Euler residual
    residual=@(coef) c_Euler_residual_cheby_compact(coef,XX,k_grid,st);
    coef_1  = fsolve(residual,coef_0,options);
    % difference of current and next policy function
    diff = max(abs(coef_1-coef_0));
    diff_list(it)= diff;
    % renew of policy function for consumption
    coef_0 = coef_1;
    it=it+1;
    % print information in iteration
    fprintf('iteration index: %i \n', it-1);
    fprintf('policy function coefficient iteration error: %e\n', diff);
end

%% calculate the policy function
XX=zeros(num_k,dim_app);
for ii=0:dim_app
   XX(:,ii+1)=chebyshevT(ii,k_grid);
end
c_fcn = XX*coef_1';
k_fcn = k_grid.^st.alpha + (1-st.delta)*k_grid - c_fcn;

%% policy function plot
figure;
plot(k_grid, k_fcn, 'r', 'LineWidth', 2); hold on;
plot(k_grid, k_grid, '--k');
xlabel('current capital', 'FontSize', 12);
ylabel('next capital', 'FontSize', 12);
title('Policy Function Plot', 'FontSize', 14);
legend('policy function for next capital', '45-degree line', 'Location', 'northwest');
grid on;


