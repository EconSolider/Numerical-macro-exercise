function residual=a_Euler_residual(c,k_grid,capital,st,c_fcn)

capital_next = capital.^ st.alpha + (1-st.delta) * capital - c;
capital_next = max(k_grid(1),capital_next); % keep kt+1 above 0


% interprete into next period consumption
c_next = interp1(k_grid,c_fcn,capital_next,'spline');

%marginal utility 
u_prime=c.^(-st.gamma);
u_prime_next=c_next.^(-st.gamma);

% Euler residual
r= st.alpha * capital_next.^(st.alpha-1);
residual = u_prime-st.beta*(r+1-st.delta)*u_prime_next;

end