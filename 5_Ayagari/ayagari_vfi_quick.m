function [mean_k,Result]=ayagari_vfi_quick(r,alpha,delta,beta,mu,b,l_grid,P)
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
num_k= 50;                     % grid size for state and control
max_k=20;                      % maximum value of capital grid  
min_k= -phi;                    % borrowing constraint
k_grid=linspace(min_k,max_k,num_k); % state of assets

%  initialize value and policy functions
num_l=length(l_grid);
state_transition_index=zeros(num_k,num_l);

v = zeros(num_k,num_l); 
Tv = zeros(num_k,num_l);
c=zeros(num_k,num_l);

% value function iteration
tol=1e-5;
diff=10000;
iter=1;
iter_max=800;
while (diff>tol && iter<iter_max)
    for ik=1:num_k
        vtemp=-100000*ones(num_k,1);
        for il=1:num_l
        %给定k和l下，valuefuncti    
            for ikp=1:num_k
                %消費(以及使得消费>0的kt+1坐标)
                c=wage*l_grid(il)+(1+r)*k_grid(ik)-k_grid(ikp);
                if c<0
                ikp_max=ikp-1;
                break
                end
                %vt+1期待値の計算
                    vpr=0;
                    for ilp=1:num_l
                    vpr=vpr+P(il,ilp)*v(ikp,ilp);
                    end
                %value functionの計算    
                vtemp(ikp)=u(c)+beta*vpr;
            end %ikp
 [Tv(ik,il),state_transition_index(ik,il)]=max(vtemp(1:ikp_max));
        end %il
    end %ik
    diff=max(max(abs(v-Tv)));
    v=Tv;
    iter=iter+1;
end
k_poli=k_grid(state_transition_index);



% stationary distribution of state given policy function of k
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
