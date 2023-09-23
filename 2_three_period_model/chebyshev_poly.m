function output = chebyshev_poly(dim_app,grid,coef)
%% vector of chebyshev basis
col=length(grid);
XX=zeros(col,dim_app);
for i=0:dim_app
 XX(:,i+1)=chebyshevT(i,grid);
end
%% chebyshev polynomial
output= XX * coef';

end