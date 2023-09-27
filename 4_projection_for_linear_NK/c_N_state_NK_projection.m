%% parameters
clc,clear;
rs=3.5/4;
beta=1/(1+rs/100);
sigma=6;
alpha=0.66;
theta=7.66;
omega=0.47;
kappa=(1-alpha)*(1-alpha*beta)/alpha*(1/sigma+omega)/(1+omega*theta);
lambda=0.048/16;
%% process of exogenous variables
rhou=0.6;
rhog=0.8;
sigma_u=1e-5;
sigma_g=1.524;
Ng=20;
Nu=20;
% marcov chain to approximate AR(1) process of gt and ut
[P_g Grid_g]=b_tauchen(rhog,sigma*rs,sigma_g,Ng,3);
[P_u Grid_u]=b_tauchen(rhou,0,sigma_u,Nu,3);

% transition prob matrix of state s=(g,u):
P=kron(P_g,P_u);

% state space s:
Ns=Ng*Nu;
Grid_s=zeros(Ns,2);
% Example for 3by3 case:
% s=[s1 s2 s3 s4 s5 s6]'
% s1=(g1,u1) s2=(g1,u2) s3=(g1,u3) 
% s4=(g2,u1) s5=(g2,u2) s6=(g2,u3)
% s7=(g3,u1) s8=(g3,u2) s9=(g3,u3)
            % In Table form
%       g1         g2           g3
% u1 s(1=0*Nu+1)  s(4=1*Nu+1) s(7=2*Nu+1)
% u2 s(2=0*Nu+2)  s(5=1*Nu+2) s(8=2*Nu+2)
% u3 s(3=0*Nu+3)  s(6=1*Nu+3) s(9=2*Nu+3)
% Let si=(g_ig,u_iu)
% i=(ig-1)*Nu+iu
% therefore, the state space should be
for ig=1:Ng
    for iu=1:Nu
    Grid_s(Nu*(ig-1)+iu,1)=Grid_g(ig);
    Grid_s(Nu*(ig-1)+iu,2)=Grid_u(iu);
    end
end

%% time iteration
tol=1e-5;       % critical value for convergence
max_iter=200;   % maximum iteration
diff=1000;      % initial value of difference of policy functions
iter=1;         % initial counter

% setting for projection
dim_app=2;             %dimention of regression equation

% generate polynominal Matrix 
% colums of X: [1 g u g^2 g*u u^2]
XX=c_gen_PolyMatrix(Grid_s,dim_app);
XX=[ones(length(XX),1) XX];
[Ns,dim_coef]=size(XX);

% initial policy functions
y_coef_cur=zeros(dim_coef,1);
pi_coef_cur=zeros(dim_coef,1);
r_coef_cur=zeros(dim_coef,1);
y_coef_next=zeros(dim_coef,1);
pi_coef_next=zeros(dim_coef,1);
r_coef_next=zeros(dim_coef,1);

% options for fsolve
options = optimoptions('fsolve', 'Algorithm', 'levenberg-marquardt', ...
'MaxFunctionEvaluations', 1000, 'Display','off');

while (diff>tol && iter<max_iter)
    %value of exogenous variable
        g=Grid_s(:,1);
        u=Grid_s(:,2);
    % expectation
        y_cur=XX*y_coef_cur;
        pi_cur=XX*pi_coef_cur;
        ye=P*y_cur;
        pie=P*pi_cur;
    % solve equations without ZLB
        x0=[y_coef_cur r_coef_cur pi_coef_cur];
        eqns=@(x)c_focs(x,XX,ye,pie,g,u,kappa,lambda,beta);
        x=fsolve(eqns,x0,options);
        y_coef_next=x(:,1);
        r_coef_next=x(:,2);
        pi_coef_next=x(:,3);
    
    %difference of policy function
    diff_y=max(abs(y_coef_cur-y_coef_next));
    diff_pi=max(abs(pi_coef_cur-pi_coef_next));
    diff_r=max(abs(r_coef_cur-r_coef_next));
    diff=max([diff_y diff_pi diff_r]);
    
    % renew policy function
    y_coef_cur=y_coef_next;
    pi_coef_cur=pi_coef_next;
    r_coef_cur=r_coef_next;
    iter=iter+1;
end

% policy funxtion for g under u=Eu=0.
base=[ones(Ng,1) Grid_g zeros(Ng,1) Grid_g.^2 zeros(Ng,1) zeros(Ng,1)];
y_poli=base*y_coef_cur;
pi_poli=base*pi_coef_cur;
r_poli=base*r_coef_cur;

figure;
subplot(311);
plot(Grid_g,y_poli,'k-','LineWidth',2.0);
grid on;
ylabel('産出ギャップ, y');

subplot(312);
plot(Grid_g,r_poli,'k-','LineWidth',2.0);
grid on;
ylabel('金利, r');

subplot(313);
plot(Grid_g,pi_poli,'k-','LineWidth',2.0);
grid on;
ylabel('インフレ,\pi');
