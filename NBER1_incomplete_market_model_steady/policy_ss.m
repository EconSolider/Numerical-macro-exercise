function [Va,a,c]=policy_ss(Pi,a_grid,y,r,beta,sigma,tol,iter_max)
%----------------------------------------------------
% Description: Backward iterating the deriviative value function
%               to calculate the policy functions
%---------------------------------------------------
% Input
% Pi:      Transition matrix of income
% a_grid:  Grid points of assets
% y:       Grid points of income
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
%----------------------------------------------------

% initial guess for Va
% cash on hand:coh(e,a),coh(i,j)=(1+r)*a_grid(j) + y(i)
coh=(1+r)*a_grid'+y'; 
% assume MPC is 0.05, initial policy function c(e,a) becomes
c=0.05*coh;
% initial value of diriviated value function
Va= (c.^(-sigma))*(1+r); % Envelop condition

 a_old=zeros(size(coh));
 iter=1;
 diff=1000;
  while diff>tol && iter<iter_max
   % step1: discounting and expectations
   Wa=beta*Pi*Va;
   % Step2: solving for asset policy using FOCs
   c_endo=Wa.^(-1/sigma);
   coh=(1+r)*a_grid'+y';
   at_eap=(c_endo+a_grid'-y')/(1+r);
   
   a=zeros(size(coh));
   for e=1:length(y)
   a(e,:)=interp1(at_eap(e,:),a_grid',a_grid'); 
   % input: samples x,f(x) + interp points x0
    end
   
   % step 3: enforcing the borrowing constraint and backing out consumption
    a=max(a,a_grid(1));
    c=coh-a; 
   %step 4: envelop condition to renew the derivative of value
   %function
   Va=(1+r)*c.^(-sigma);
   
   diff=max(max(abs(a-a_old)));
   a_old=a;
   iter=iter+1;
  end
end