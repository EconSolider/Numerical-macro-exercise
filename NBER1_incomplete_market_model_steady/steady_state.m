function str=steady_state(Pi,a_grid,y,y_dis,r,beta,sigma,tol,iter_max)
%----------------------------------------------------------
% Description: function to calculate steady state of heterogenous
%               agent model
%-----------------------------------------------------------
% Warning: This function need other two functions to operate
%                poicy_ss.m 
%                distribution_ss.m
%------------------------------------------------------------
% Input
% Pi:      Transition matrix of income
% a_grid:  Grid points of assets
% y:       Grid points of income
% y_dis:   Stationary distribution of income process
% r:       Interest rate
% beta:    Discount rate;  
% sigma:   Risk aversion rate
% tol:     Convergence criterioon
% max_iter:Maximun iterate times
%----------------------------------------------------
% Output
% Va(e,a): Deriviative value function
% a'(e,a): Policy function of assets
% c(e,a):  Policy function of consumption
% D:       stationary distribution
% a_i:     Lottery number of a'(e,a)
% a_pi:    Lottery probability of a'(e,a)
% A:       Aggregate assets in steady state
% C:       Aggregate consumption in steady state
%----------------------------------------------------

% calculation of steady state policy function a'(a,e), c(a,e)
[Va,a,c]=policy_ss(Pi,a_grid,y,r,beta,sigma,tol,iter_max);
% calculation of steady state distribution D(a,e)
[D,a_i,a_pi]=distribution_ss(Pi,a,a_grid,y_dis,tol,iter_max);

% aggregate asset and consumption
A=sum(a(:) .* D(:));
C=sum(c(:) .* D(:));

% save inputs
str.Pi=Pi; 
str.a_grid=a_grid; 
str.y=y;       
str.y_dis=y_dis;
str.r=r;        
str.beta=beta;   
str.sigma=sigma; 

% save outputs
str.Va=Va;
str.a=a;
str.c=c;
str.D=D;   
str.a_i=a_i; 
str.a_pi=a_pi;
str.A=A;
str.C=C;

end