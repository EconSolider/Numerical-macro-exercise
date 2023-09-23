function residual=c_Euler_residual_cheby_compact(coef,XX,k_grid,st)

%% Description
% inputs
    % coef:     coefficients of Chebyshev approximation (what we want)
    % dim_app:  dimension of Chebyshev approximation
    % k_grid:   descreted state variable (state space)
    % st:       structure for exogenous variables and parameters
% output:
    % residual: Euler residual (minimization target)
    
%% chebyshev approximation for consumption function
c_fcn = XX*coef';

%% next period capital and consumption
k_next = k_grid.^st.alpha + (1-st.delta)*k_grid - c_fcn;
k_next = max(k_next,k_grid(1));
c_fcn_next = interp1(k_grid,c_fcn,k_next,'spline');

% marginal utility
u_prime = c_fcn.^(-st.gamma);
u_prime_next = c_fcn_next.^(-st.gamma);

for ii=1:length(c_fcn)
    if c_fcn(ii)<0
        u_prime=10000;
    end
    if c_fcn_next(ii)<0
        u_prime_next=10000;
    end
end


%% Euler residual
r = st.alpha * k_next .^(st.alpha-1);
residual = u_prime - ...
            st.beta * (r+1-st.delta).* u_prime_next;
end