function [x,y, uApprox] = poisson2DFD9(f,u0,L,N)
hj = L/(N+1);
x_matrix = 0:hj:L;
y_matrix = 0:hj:L;
[x,y] = meshgrid(x_matrix,y_matrix);

xIndex = 2:N+1;
yIndex = 2:N+1;


F = f(x(xIndex, yIndex), y(xIndex, yIndex))
F = reshape(F, N^2,1 );
[X, Y, uApprox_5p, F_5p ] = poisson2DFD5(f, u0, L, N);
F = F - ((hj^2)*F_5p)/12;
F = reshape(F, N, N);

% create matrix A 
  I = ones(N,1);
  S = spdiags([I 10*I I], [-1 0 1], N, N);
  T = spdiags([-1/2*I I -1/2*I], [-1 0 1], N, N);
  A = kron(T,S)+kron(S,T); 
  A = (1/(6*hj^2))*full(A);
  
% 
  
  
  