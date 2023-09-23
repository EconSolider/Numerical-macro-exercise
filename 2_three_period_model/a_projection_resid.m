function residual=a_projection_resid(coef,grid,str,dim_app,old_period,coef_aO)


%% calculate the saving fucntion / policy function
num=length(grid);
XX=zeros(num,dim_app);

for i=0:dim_app
 XX(:,i+1)=chebyshevT(i,grid);
end

% a_prime is the next period asset (saving function)
a_prime = XX*coef';

%% prepare for Euler residual

% current consumption and marginal utility
if old_period==1
    c_cur = str.YM + str.R * grid - a_prime;
else
    c_cur = str.YY + grid - a_prime;
end

u_prime_cur = c_cur.^(-str.gamma);

% next period consumption and marginal utility
if old_period==1
    c_next = str.R * a_prime + str.ss;
else
    aO = chebyshev_poly(dim_app,grid,coef_aO);
    c_next = str.YM + str.R * a_prime - aO;
end

u_prime_next = c_next.^(-str.gamma);

for i=1:num
   if c_cur(i) < 0
    u_prime_cur(i)=100000;
   end
   if c_next(i) < 0
    u_prime_next(i)=-100000;
   end
end

residual=str.beta * str.R * u_prime_cur ./ u_prime_next -1;

end