function [y_grid,dis,P]=dis_AR1_Rouwenhorst(rho,sigma,N)
% Description: discretize the AR1 process by Rouwenhorst method
%----------------------------------------------------------------
% input
% rho: persistence of AR(1)
% sigma: standard diviation of normal innovation in AR(1)
% N: number of grid points
%---------------------------------------------------------------
% choose inner-switching probability p to match persistence rho
p = (1+rho)/2;
% start with states from 0 to n_e-1, 
% scale by alpha to match standard deviation sigma
 e = 0:(N-1);
 alpha = 2 * sigma / sqrt(N-1);
 e = alpha * e;
% obtain Markov transition matrix Pi 
P = Rouwenhorst_P(N,p);
% its stationary distribution
tol=1e-10;
dis=ones(N,1)./N;
diff=1000;
while diff>tol
   dis_new=P'*dis;
   diff=max(max(abs(dis_new-dis)));
   dis=dis_new;
end
% e is log income, get income y and scale so that mean is 1
y_grid=exp(e)/dot(dis,exp(e));
end

function P=Rouwenhorst_P(N,p)
% Description: compute the marcove matrix of rouwenhorst
%--------------------------------------------------------
% input 
% N: number of grid points
% p: parameter of recursion
% output
% P: marcov matrix
%-------------------------------------------------------
% base case Pi_2
Pi=[p 1-p;1-p p];
for n = 3:N
        Pi_old = Pi;
        Pi = zeros(n, n);
        Pi(1:end-1, 1:end-1) = Pi(1:end-1, 1:end-1) + p * Pi_old;
        Pi(1:end-1, 2:end) = Pi(1:end-1, 2:end) + (1 - p) * Pi_old;
        Pi(2:end, 1:end-1) = Pi(2:end, 1:end-1) + (1 - p) * Pi_old;
        Pi(2:end, 2:end) = Pi(2:end, 2:end) + p * Pi_old;
        Pi(2:end-1, :) = Pi(2:end-1, :) / 2;
end
P=Pi;
end