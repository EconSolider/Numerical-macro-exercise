function [PI,a_grid]=tauchen(rho,mu,sigma,N,d)
 %A(t+1)=rho*A(t)+(1-rho)*mu+eps, eps～N(0,sigma^2)
 %rho  : 持续性参数
 %mu   : 随机过程的均值
 %sigma: eps的标准差
 %N    : 离散化空间的size
 %d    : 状态空间上下限对平均的偏离值
 %
 % 返回一个状态转移矩阵和各状态的具体数值
 %
sigma_a=sqrt(sigma^2/(1-rho^2)); %随机状态的标准差
a_min=mu-d*sigma_a;              %离散化状态空间的上下限
a_max=mu+d*sigma_a;
a_grid=linspace(a_min,a_max,N);  %离散化的状态空间
 
    % 为之后求概率分布做准备
    % m(i)=(A(i)+A(i+1))/2
    for ii=1:length(a_grid)-1
    m_list(ii)=(a_grid(ii+1)+a_grid(ii))/2;
    end

    for ii=1:length(a_grid)
        %求PI(i,1)
        PI(ii,1)= normcdf((m_list(1)-(1-rho)*mu-rho*a_grid(ii))/sigma); %#ok<AGROW>
        %求PI(i,N)
        PI(ii,N)= 1-normcdf((m_list(N-1)-(1-rho)*mu-rho*a_grid(ii))/sigma); %#ok<AGROW>
        %求PI(i,k) k=2,...,N-1
        for k=2:length(m_list)
            PI(ii,k)=normcdf((m_list(k)-(1-rho)*mu-rho*a_grid(ii))/sigma)...
           -normcdf((m_list(k-1)-(1-rho)*mu-rho*a_grid(ii))/sigma);
        end 
    end
 a_grid=a_grid';
end