function residual=b_Euler_residual(c_fcn,k_grid,st)

k_prime_grid = k_grid.^st.alpha + (1-st.delta)*k_grid - c_fcn;
% keep k_prime above 0
for ii=1:length(k_prime_grid)
   k_prime_grid(ii)=max(k_grid(1),k_prime_grid(ii)); 
end

c_fcn_next = interp1(k_grid,c_fcn,k_prime_grid);

% Euler residual
r=st.alpha*k_prime_grid.^(st.alpha-1);
residual=c_fcn.^(-st.gamma)- st.beta * (r+1-st.delta) .* c_fcn_next.^(-st.gamma);
end