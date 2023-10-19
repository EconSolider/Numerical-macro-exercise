function a_grid=DET_discrete(a_min,a_max,a_num)
%----------------------------------------------------
% Description:a function to discrete the asset space by Doube-exponential
%               transformation method. 
%----------------------------------------------------
% Purpose: to get a grid of assets with more points around two sides of
%          assets
%------------------------------------------------------
% input
% a_num: number of grid points of asset space
% a_min: under limitation of assets
% a_max: above limitation of assets
% output
% a_grid: grid points of asset space
%-----------------------------------------------------
u_max=log(1+log(a_max-a_min+1));
u_grid=linspace(0,u_max,a_num)';
a_grid=a_min+exp(exp(u_grid)-1)-1;
end