function [mean_k,Result]=ayagari_vfi_slow(r,alpha,delta,beta,mu,b,l_grid,P)
%% Settings for iteration

%write wage as a function of interest rate 
wage = (1-alpha)*((alpha/(r+delta))^alpha)^(1/(1-alpha));
% borrowing limit
if r<=0
   phi = b;
else
   phi = min(b, wage*l_grid(1)/r);           
end
% utility function
u=@(c) c.^(1-mu)/(1-mu);

% capital grid (need define in each iteration since it depends on r/phi)
num_k=50;                     % grid size for state and control
max_k=20;                      % maximum value of capital grid  
min_k= -phi;                    % borrowing constraint
k_grid=linspace(min_k,max_k,num_k)'; % state of assets

%  initialize value and policy functions
num_l=length(l_grid);
state_transition_index=zeros(num_k,num_l);
v = zeros(num_k,num_l); 
Tv = zeros(num_k,num_l);
c=zeros(num_k,num_l);

%% value function iteration
tol=1e-5;
diff=10000;
iter=1;
iter_max=800;
while (diff>tol && iter<iter_max)
    for ik=1:num_k
        for il=1:num_l
            %给定资产kt，不同劳动状态lt，计算可实行的消费c
            %c(l,k，k')=(1+r)*k+wl-k'
            c=(1+r)*k_grid(ik)+wage*l_grid(il)-k_grid;
            %识别使得消费小于0的index
            neg_index= c<0;
            c(neg_index)=NaN;
            %计算不同状态下的效用函数
            util(:,il)=u(c);
            %将负的消费对应的效用设为负值
            util(neg_index,il)=-1e10;
        end
        %选择使得值函数最大化的k'的位置
        [Tv(ik,:),state_transition_index(ik,:)]=max(util+beta*(v*P));
    end
    diff=max(max(abs(v-Tv)));
    v=Tv;
    iter=iter+1;
end
k_poli=k_grid(state_transition_index);

%% Stationary distribution of state given policy function of k
iter=1;
diff=1000;
mean_cur=ones(num_k,num_l)/(num_k*num_l);
mean_next=zeros(num_k,num_l);

while (diff>tol && iter<iter_max)
    for ik=1:num_k
        for il=1:num_l
        k_next_index=state_transition_index(ik,il);
            for ipl=1:num_l
            mean_next(k_next_index,ipl)=mean_next(k_next_index,ipl)+P(il,ipl)*mean_cur(ik,il);
            end
        end
    end
    diff=max(max(abs(mean_next-mean_cur)));
    mean_cur=mean_next;
    iter=iter+1;
    mean_next=zeros(num_k,num_l);
end

mean_k=sum(sum(mean_cur.*k_poli));

Result.k_grid=k_grid;
Result.k_poli_fun=k_poli;
Result.stationary_state_dis=mean_cur;
end
