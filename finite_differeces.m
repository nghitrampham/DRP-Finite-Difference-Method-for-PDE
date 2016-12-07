function [x,U, E_2norm] = finite_differeces(a,b,ua,ub,n)
h = (b-a)/n; h1=h*h;
A = sparse(n-1,n-1);
F = zeros(n-1,1);
for i=1:n-2,
  A(i,i) = -2/h1; A(i+1,i) = 1/h1; A(i,i+1)= 1/h1;
end
  A(n-1,n-1) = -2/h1;

for i=1:n-1,
  x(i) = a+i*h;
  
end
  F(1) = F(1) - ua/h1;
  F(n-1) = F(n-1) - ub/h1;
U = A\F;
 

i = 1:n-1; 
u_true(i) = ((ub-ua)/b)*(x(i))+ua;
u_true = u_true';
E = U - u_true;

%using 2_norm
E_2norm = [];
j = 1:length(E);
E_2norm = (h*(sum(abs(E(j)).^2)))^(1/2)
% E_norm = (h^(1/2))*norm(E)
