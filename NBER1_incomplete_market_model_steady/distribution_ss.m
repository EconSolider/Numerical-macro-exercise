function [D,a_i,a_pi]=distribution_ss(Pi,a,a_grid,y_dis,tol,iter_max)
%----------------------------------------------------
% Description:a function to calculate the steady state distribution
%             D(e,a) of household asset and income
%---------------------------------------------------
% Input
% Pi:       Transition matrix of income process y(e)
% a:        a(e,a), policy function of t+1 assets
% a_grid:   Gird points of assets
% y_dis:    Stationary distribution of income process
% tol:      Convergence criterion
% iter_max: Maximum of iteration
%----------------------------------------------------
% Output 
% D: stationary distribution
% a_i: Lottery number of a'(e,a)
% a_pi: Lottery probability of a'(e,a)
%---------------------------------------------------
[e_num,a_num]=size(a);
a_i=zeros(size(a));
a_pi=zeros(size(a));
% Using lottery to iterate on distribution
    for e=1:e_num
        for jj=1:a_num
        [a_i(e,jj),a_pi(e,jj)]=get_lottery(a(e,jj),a_grid);
        end
    end
% Initial D, use stationary distribution for e(income),
%           and uniform over a
    D0=y_dis*ones(1,a_num)/a_num;

    % iteration for steady state distribution
diff=1000;
iter=1;
    while diff>tol && iter<iter_max
        D=forward_policy(D0, a_i, a_pi,Pi);
        diff=max(max(abs(D0-D)));
        D0=D;
    end

end

function [a_i, a_pi] = get_lottery(ap, a_grid)
%----------------------------------------------
% Problem: if we start at a gridpoint
%          the asset policy function takes us off the grid
%
% Solution: Assume household follows a "lottery", 
%           going to each of the nearest two gridpoints 
%           with some probability, such that on average, 
%           the household chooses the right asset value
%----------------------------------------------
% Description: a function to get the asset lotteries of households
%-----------------------------------------------
% Input
% ap:      asset policy function a'(e,a).
% a_grid:  grid points of assets
%--------------------------------------------------
% Output
% a_i:  position of each asset policy's lottery  
% a_pi: probability of each asset policy's lottery
%--------------------------------------------------

% Step 1: Find the i such that a' lies between gridpoints 
%         a_i and a_(i+1)
    a_i = find(a_grid <= ap,1,'last') - 1;
    if isempty(a_i)
        a_i = 0;
    end
% Step 2: Obtain lottery probabilities pi
    a_pi = (a_grid(a_i + 2) - ap) / (a_grid(a_i + 2) - a_grid(a_i + 1));
end

function Dnext = forward_policy(D, a_i, a_pi,Pi)
%--------------------------------------------
% Description: function to do forward iteration of stationary
% distribution D(e,a)
%--------------------------------------------
% Input
% Dt:   Current stationary distribution D(e,a)
% a_i:  Lottery values of policy function a'(e,a)
% a_pi: Lottery probability of asset
% Pi:   Transition matrix of income
%----------------------------------------------
% Output
% Dt+1: Next period stationary distribution Dt+1(e',a')
%----------------------------------------------
% Preallocate the resulting matrix
    [m, n] = size(D);
    Dend = zeros(m, n);
  for e = 1:m
    for a = 1:n
     % Send pi(e,a) of the mass to gridpoint i(e,a)
      Dend(e, a_i(e,a)+1) = Dend(e, a_i(e,a)+1) + a_pi(e,a)*D(e,a);
      % Check bounds to ensure we don't exceed matrix dimensions
      if a_i(e, a) + 2 <= n
      % Send 1-pi(e,a) of the mass to gridpoint i(e,a)+1
      Dend(e, a_i(e,a)+2) = Dend(e, a_i(e,a)+2) + (1-a_pi(e,a))*D(e,a);
      end
    end
 end
    % Dend is the "end-of-period distribution" D_end(e,a');
    % Dt+1(e',a')=Pi'*D_end(e,a');
    Dnext=Pi' * Dend;
end