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