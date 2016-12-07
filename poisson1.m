function [x,uApprox, E_norm] = poisson1(L,a,b,N)
%INPUT
%u_xx = 0 where x \in [0,L] 
% a,b = interval points 
% N = number of grid points 
%OUTPUT: x(i) is grid point where i is in [1,n-1]
%        uApprox = approximation solution at each grid point

%generating the grid points 
close all 
h = L/N;
h1 = h*h;
i = 1:N-1;
x(i) = i*h

%initilize 
F = zeros(N-1,1);
A = [];
U = [];
% u0 = a, u(L) = u_n = b 

F(1) = -a/h1;
F(end) = -b/h1;

for i = 1:length(F)-1 
    A(i,i) = -2;
    A(i+1,i) = 1;
    A(i, i+1) = 1;
    A(length(F), length(F)) = -2; 
end
A = A/h1;
uApprox = A\F;

%Calculating error
i = 1:N-1; 
u_true(i) = ((b-a)/L)*(x(i))+a;
u_true = u_true'
E = uApprox - u_true;

%using 2_norm
% E_2norm = [];
% j = 1:length(E);
% E_2norm = (h*(sum(abs(E(j)).^2)))^(1/2);
 E_norm = norm(E);
 
 plot(x,u_true,'k-',x,uApprox,'r*')


% loglog(x,E, '-s')



