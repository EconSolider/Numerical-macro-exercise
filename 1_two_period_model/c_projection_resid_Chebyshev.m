function resid=c_projection_resid_Chebyshev(coef,a_grid,beta,gamma,R,Y,ss,dim_app)
%输入:
% coef: 近似多项式的系数(求取0点的目标)
% a_grid: 资产a(状态变量的散点)
% beta,gamma,R,Y,ss: 欧拉方程的必要参数
   % gamma:风险回避度，Y:收入，ss:社保收入
% dim_app: 多项式的次数

num_a=length(a_grid);
XX = zeros(num_a, dim_app+1);

% 对于每一个切比雪夫多项式的次数，计算其对应的选点的值。
for i=0:dim_app
   XX(:,i+1)=chebyshevT(i,a_grid);
end
% policy function中的saving function。
a_prime=XX*coef';

%下面准备欧拉residual
Cy=Y+a_grid-a_prime;
Co=ss+R*a_prime;

%计算边际效用
uy_prime=zeros(num_a,1);
uo_prime=zeros(num_a,1);
for i=1:num_a
  if Cy(i)>0    
    uy_prime(i)=Cy(i)^(-gamma);
  else
    uy_prime(i)=10000;
  end
  
  if Co(i)>0
    uo_prime(i)=Co(i)^(-gamma);
  else
    uo_prime(i)=10000;
  end
end

%计算欧拉残差
resid=beta*R*uo_prime./uy_prime-1;

end