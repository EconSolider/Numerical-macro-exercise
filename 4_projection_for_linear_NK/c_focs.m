function F=c_focs(x,XX,ye,pie,g,u,kappa,lambda,beta)
y=XX*x(:,1);
r=XX*x(:,2);
pi=XX*x(:,3);

   A=r+y-ye-pie-g;
   B=pi-kappa*y-beta*pie-u;
   C=lambda*y+kappa*pi;
   F=[A;B;C];
end