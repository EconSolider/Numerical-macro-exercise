%% test
clc,clear;close all;
r=0.02;
alpha=0.35;
delta=0.025;
wage = (1-alpha)*((alpha/(r+delta))^alpha)^(1/(1-alpha));
beta=0.9;

num_l=5;
rho= 0.6;          % autoregressive coefficient
sigma= 0.4;        % standard diviation of shocks
[P,log_l_grid]=tauchen(rho,0,sigma,num_l,3);
l_grid=exp(log_l_grid);

num_k=90;                     % grid size for state and control
max_k=20;                      % maximum value of capital grid  
min_k=-min(3, wage*l_grid(1)/r);  % borrowing constraint
k_grid=linspace(min_k,max_k,num_k)'; % state of assets
mu=3;


u=@(c) c.^(1-mu)/(1-mu);

v = zeros(num_k,num_l); 
Tv=zeros(num_k,num_l);
state_transition_index=zeros(num_k,num_l);

diff=1000;
tol=1e-6;
iter=1;
iter_max=500;

    for ik=1:num_k
    for il=1:num_l
    for ikp=1:num_k
    c(ik,il,ikp)=(1+r)*k_grid(ik)+wage*l_grid(il)-k_grid(ikp);
    end
    end
    end
neg=c<0;
util=u(c);
util(neg)=-1e15;

RHS=zeros(num_k,1);
while diff>tol && iter<iter_max
for ik=1:num_k
    for il=1:num_l
        for ikp=1:num_k    
        RHS(ikp)=util(ik,il,ikp)+beta*v(ikp,:)*P(il,:)';
        end
        [Tv(ik,il),state_transition_index(ik,il)]=max(RHS);
    end
end
        diff=max(max(abs(v-Tv)));
        v=Tv;
        iter=iter+1;
end

k_poli=k_grid(state_transition_index);
hold on;
plot(k_grid,k_poli(:,1)); 
plot(k_grid,k_poli(:,3));
plot(k_grid,k_poli(:,5));
hold off;

