%% parameters
clc,clear;
close all;
mu     = 3;               % risk aversion (=3 baseline)             
beta   = 0.96;            % subjective discount factor 
delta  = 0.08;            % depreciation
alpha  = 0.36;            % capital's share of income
b      = 3;               % borrowing limit

%% Discrete the process of labor
% tauchen.m to compute transition prob matrix and grid 
% approximate labor endowment shocks with 10 states Markov chain
% log(l_t)=rho*log(l_t-1)+eps_t;  eps_t~ N(0,sig^2)
num_l= 20;         % number of discretized states
rho= 0.6;          % autoregressive coefficient
sigma= 0.4;        % standard diviation of shocks
[P,log_l_grid]=tauchen(rho,0,sigma,num_l,3);
l_grid=exp(log_l_grid);

%% calculate the stationary distribition
%setting for iteration
tol=1e-5;
max_iter=200;
diff=10000;
iter=1;
%initial distribution for iteration
dis_cur=ones(num_l,1)/num_l; 
dis_next=zeros(num_l,1);

% pi{n} = P*pi{n-1} , loop until || pi{n}-pi{n-1} ||<eps
while (diff>tol && iter<max_iter)
dis_next=P'*dis_cur;
diff=max(abs(dis_next-dis_cur));
iter=iter+1;
dis_cur=dis_next; %renew the distribution
end
dis_l=dis_next;
% Total labor supply (expected value of labor)
L=dis_l'*l_grid;

%% compute individual policy function and capital supply Given interest rate
% discrete the space of r
num_R=50;
max_R=(1-beta)/beta-0.001;
min_R=-0.02;
R_grid=linspace(min_R,max_R,num_R);
K=zeros(num_R,1);

% calculate supply curve of capital
for ii=1:num_R
    r=R_grid(ii);
    % self-made function,
    % value function iteration for saving 
    % and stationary state distribution
    [K(ii,1),~]=ayagari_vfi_slow(r,alpha,delta,beta,mu,b,l_grid,P);
end



%% calculate supply and demand curve of capital
K_d=(L*(alpha./(R_grid+delta)).^(1/(1-alpha)))';

figure;
xlabel('capital supply and demand')
    ylabel('Interest rate')
    hold on
    plot(K,R_grid,'r--','LineWidth',2)
    plot(K_d,R_grid,'b','LineWidth',2)
    lg = legend({'capital supply','capital demand'}, ...
    'FontSize',12,'FontName',"Times New Roman");
    hold off
    grid on
    set(gca,'FontSize',12,'FontName','Times New Roman');

 %% COMPUTE K and r in EQ
 %capital supply and capital demand corresponding to r
 K_s=@(r)ayagari_vfi_slow(r,alpha,delta,beta,mu,b,l_grid,P);
 K_d=@(r)(L*(alpha./(r+delta)).^(1/(1-alpha)));
 r0=0.01;
 %使得供需相等的均衡利息,
 EQ_r=fsolve(@(r)(K_s(r)-K_d(r)),r0);
 %均衡下的资本存量，policy function与资产分布
 [EQ_K,EQ_result]=K_s(r);
 
 El=int16(length(l_grid)/2);
 figure('Name','Stationary asset distribution');
 plot(EQ_result.k_grid,EQ_result.stationary_state_dis(:,El))
 
 figure('Name','Policy function');
 plot(EQ_result.k_grid,EQ_result.k_poli_fun(:,El))