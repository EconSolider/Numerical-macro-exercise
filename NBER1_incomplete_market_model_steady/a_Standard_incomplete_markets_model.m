clc,clear; close all;
%% Programs for incomplete markets model for heterogenous agents
%% Discretize the stochastic income space
e_num=7;
[y,y_dis,Pi] = dis_AR1_Rouwenhorst(0.975, 0.7, e_num);
% verify mean and std of income
Ey=dot(y,y_dis)
% standard diviation of innovation
mean_log_y=dot(log(y),y_dis);
std=sqrt(dot(y_dis,(log(y)-mean_log_y).^2)) 
%% Discritize the asset space
a_min=-Ey/(1+0.01/4);
a_max=80;
a_num=2000;
a_grid=DET_discrete(a_min,a_max,a_num);
% histgram of asset space
hist(a_grid);

%% Backward iteration
%parameters
r=0.01/4; %steady state interest rate
beta=1-0.08/4; %discount rate
sigma=1; % risk aversion rate

% initial guess for Va
% cash on hand:coh(e,a),coh(i,j)=(1+r)*a_grid(j) + y(i)
coh=(1+r)*a_grid'+y'; 
% assume MPC is 0.05, initial policy function c(e,a) becomes
c=0.05*coh;
% initial value of diriviated value function
Va= (c.^(-sigma))*(1+r); % Envelop condition

diff=1000;
tol=1e-9;
iter=1;
iter_max=1000;
a_old=zeros(size(coh));
% Description: backward induction
while diff>tol && iter<iter_max
   % step1: discounting and expectations
   Wa=beta*Pi*Va;
   % Step2: solving for asset policy using FOCs
   c_endo=Wa.^(-1/sigma);
   coh=(1+r)*a_grid'+y';
   at_eap=(c_endo+a_grid'-y')/(1+r);
   
   a=zeros(size(coh));
   for e=1:e_num
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

%% 程序封装: Backward iteration
[Va,a,c]=policy_ss(Pi,a_grid,y,r,beta,sigma,tol,iter_max);

%% Figures: policy function c(e,a) given e
figure;
hold on;
for e=1:e_num
 plot(a_grid,c(e,:), 'DisplayName', sprintf('Line y=%.2f', y(e)))
end
hold off;
set(gca, 'FontSize', 12);
hLegend = legend('show');
set(hLegend, 'FontSize', 12);
xlabel('asset (a)', 'FontSize', 14);
ylabel('consumption (c)', 'FontSize', 14);

%% Figures: saving function a'(e,a)-a given e
figure;
hold on;
for e=1:e_num
 plot(a_grid,a(e,:)-a_grid', 'DisplayName', sprintf('Line y=%.2f', y(e)))
end
hold off;
set(gca, 'FontSize', 12);
hLegend = legend('show');
set(hLegend, 'FontSize', 12);
xlabel('asset (a)', 'FontSize', 14);
ylabel('saving', 'FontSize', 14);
h1=yline(0, '--r', 'LineWidth', 1.5);
set(h1, 'HandleVisibility', 'off');

%% Marginal propensity to consume
% Preallocate mpcs array
mpcs = zeros(size(c));
% Symmetric differences away from boundaries
mpcs(:, 2:end-1) = (c(:,3:end) - c(:,1:end-2))./(a_grid(3:end) - a_grid(1:end-2))'/ (1+r);
% Asymmetric first differences at boundaries
mpcs(:,1)  = (c(:,2) - c(:,1)) / (a_grid(2)-a_grid(1))' / (1+r);
mpcs(:,end) = (c(:,end)-c(:,end-1)) / (a_grid(end) - a_grid(end-1))' / (1+r);
% Special case of constrained
mpcs(a == a_grid(1)) = 1;

figure;
hold on;
for e=1:7
 plot(a_grid(1:120),mpcs(e,1:120), 'DisplayName', sprintf('Line y=%.2f', y(e)),'LineWidth',1.5)
end
hold off;
set(gca, 'FontSize', 12);
hLegend = legend('show');
set(hLegend, 'FontSize', 12);
xlabel('asset (a)', 'FontSize', 14);
ylabel('MPC', 'FontSize', 14);

%% forward iteration to obtain the distribution D(e,a)
a_i=zeros(e_num,a_num);
a_pi=zeros(e_num,a_num);
% Using lottery to iterate on distribution
for e=1:e_num
    for jj=1:a_num
    [a_i(e,jj),a_pi(e,jj)]=get_lottery(a(e,jj),a_grid);
    end
end
% as initial D, use stationary distribution for e, 
% and uniform over a
D0=y_dis*ones(1,a_num)/a_num;
% iteration for steady state distribution
diff=1000;
iter=1;
while diff>tol && iter<iter_max
    D=forward_policy(D0, a_i, a_pi,Pi);
    diff=max(max(abs(D0-D)));
    D0=D;
end
% check if D is a distribution
sumD=sum(D,'all')

%% 程序封装: forward iteration
[D,a_i,a_pi]=distribution_ss(Pi,a,a_grid,y_dis,tol,iter_max);

%% figure of steady state distribution of assets
plot(a_grid, cumsum(sum(D,1)),'LineWidth',1.5);
% sum(D, 1):沿着行维度求和
% cumsum:累积和
xlabel('asset','FontSize', 14);
ylabel('Cumulative distribution','FontSize', 14);
set(gca, 'FontSize', 12);
h1=xline(0, '--r', 'LineWidth', 1.5);
set(h1, 'HandleVisibility', 'off');

%% Comparative statics in r
clc,clear;
% grid,steady state distribution, transition matrix of income y
[y,y_dis,Pi] = dis_AR1_Rouwenhorst(0.975, 0.7, 7);
Ey=dot(y,y_dis);

% parameters
r=0.01/4; %steady state interest rate
beta=1-0.08/4; %discount rate
sigma=1; % risk aversion rate
tol=1e-10;
iter_max=1000;

% Given different r,get capital supply curve
r_list=r+linspace(-0.02,0.015,15);
A=zeros(size(r_list));
C=zeros(size(r_list));
for ii=1:length(r_list)
r=r_list(ii);
% grid of assets
a_min=Ey/(1+r);
a_max=80;
a_num=2000;
a_grid=DET_discrete(a_min,a_max,a_num);

[Va,a,c]=policy_ss(Pi,a_grid,y,r,beta,sigma,tol,iter_max);
[D,a_i,a_pi]=distribution_ss(Pi,a,a_grid,y_dis,tol,iter_max);
% aggregate asset A,将矩阵化为向量并点乘分布
A(ii)=sum(a(:) .* D(:));
% aggregate consumption C
C(ii)=sum(c(:) .* D(:));
end

% supply curve of capital
figure;
plot(r_list,A,'LineWidth',1.5);
xlabel('interest rate (r)','FontSize', 14);
ylabel('capital supply (A)','FontSize', 14);
set(gca, 'FontSize', 12);

%% 程序封装: steady_state
Results=steady_state(Pi,a_grid,y,y_dis,r,beta,sigma,tol,iter_max);

%% Calibration
% Total asset (A) is constant times to annual GDP, 
%   in USA,it's about 140%. 
% In the model, A is quarterly data, and quarterly GDP=1. 
%   So, A=1.4*4=5.6; 

% Define the function we want to find the root for
steady_state_diff=@(beta) steady_state(Pi,a_grid,y,y_dis,r,beta,sigma,tol,iter_max).A - 5.6;
% Find the root using fzero
beta_calib = fzero(steady_state_diff, [0.8, 0.995]);
fprintf('Calibrated beta: %.2f\n', beta_calib);

% Generate calibrated result
ss_calib=steady_state(Pi,a_grid,y,y_dis,r,beta_calib,sigma, ...
                        tol,iter_max);
% confirm the calibration is effective.
fprintf('Calibrated A: %.2f\n', ss_calib.A);
% check steady state budget balance C=1+rA
fprintf('1+rA-C: %f\n',1+ss_calib.r*ss_calib.A-ss_calib.C);

%% Steady state in general equilibrium
clc,clear;

% standard settings
[y,y_dis,Pi] = dis_AR1_Rouwenhorst(0.975, 0.7, 7);
r=0.01/4; %steady state interest rate
beta=1-0.08/4; %discount rate
sigma=1; % risk aversion rate
tol=1e-10;
iter_max=1000;

a_min=0;
a_max=80;
a_num=1000;
a_grid=DET_discrete(a_min,a_max,a_num);

% Quarterly bonds/GDP is 560%, and quarterly GDP is 1
B=5.6;
% Labor tax needed to balance steady-state government budget
tau=r*B;
% use y, which had mean 1, for the labor endowment process
e=y; 

% Calibrate beta to be consistent with steady-state asset 
%   market clearing
steady_state_diff=@(beta) steady_state(Pi,a_grid,(1-tau)*e,y_dis,r, ...
                        beta,sigma,tol,iter_max).A - B;
beta_ge = fzero(steady_state_diff, [0.5, 0.995]);

% Compute full steady state with beta_ge
ss_ge=steady_state(Pi,a_grid,(1-tau)*e,y_dis,r, ...
                        beta_ge,sigma,tol,iter_max);

% General equilibrium counter-factuals (vary income risk)
thetas=linspace(0.3,1.2,8)';
rs=zeros(size(thetas));
for ii=1:length(thetas)
   [e_new,e_dis_new,Pi_new]= dis_AR1_Rouwenhorst(0.975, thetas(ii), 7);
   steady_state_fun=@(r) steady_state(Pi_new,a_grid,e_new, ...
                    e_dis_new,r,beta,sigma,tol,iter_max).A - B;
   rs(ii)=fsolve(steady_state_fun,0.01);
end
figure;
plot(thetas,4*rs,'LineWidth',1.5);
xlabel('Standard deviation of income ','FontSize', 14);
ylabel('Equilibrium real interest rare (Annualized)','FontSize', 14);
set(gca, 'FontSize', 12);

%% Expectation functions + Autocorelation of consumption
T=40;
ctilde=ss_ge.c-dot(ss_ge.D,ss_ge.c);
E_ctilde=expectation_functions(ctilde,ss_ge.Pi,ss_ge.a_i,ss_ge.a_pi,T);
Autocov_c = zeros(T,1);
for j = 1:T
    temp = reshape(squeeze(E_ctilde(j,:,:)),[],1);
    Autocov_c(j) = dot(reshape(ss_ge.D,[],1), conj(reshape(ctilde,[],1)) .* temp);
end
Autocorr_c=Autocov_c./ Autocov_c(1);

figure;
plot(0:length(Autocorr_c)-1,Autocorr_c,'LineWidth',1.5);
title('Autocorrelation of consumption');
xlabel('Horizon in quarters');
ylabel('Correlation');
set(gca, 'FontSize', 12);

