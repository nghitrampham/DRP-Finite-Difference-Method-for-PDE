function [x,uApprox, E_2norm] = poisson2(L,N)
close all 
%generating the grid points 
h = L/(N+1);
h1 = h*h;
x = h*(0:N+1)';

%initialize 
A = [];
U = [];

F = zeros(N+2,1);
F(1) = 0;
F(2:N+2) = cos(pi*x(1:N+1)/(2*L));
F(2) = (h/2)*F(2);
% F(end) = F(end) -b/h1;
F

A = zeros(N+2,N+2);
A(1,1) = h1;
A(2,1) = -h;  A(2,2) = h;
for i = 3:N+2
    A(i,i-2) =  1;
    A(i,i-1)   = -2;
    A(i,i) =  1; 
end
A
A = A/h1;

uApprox = A\F;

%Calculating error 
u_true = -((2*L/pi)^2)*cos((pi*x(:))/(2*L)) + ((2*L/pi)^2);
E = uApprox - u_true;

%using 2_norm
E_2norm = norm(E);
plot(x,uApprox, 'r*', x,u_true, 'k-')
