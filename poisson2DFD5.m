function [ x,y, uApprox,  uApprox_5p] = poisson2DFD5(f, u0, L, N)

% initializing the variable 
hj = L/(N+1);
x_matrix = 0:hj:L;
y_matrix = 0:hj:L;
[x,y] = meshgrid(x_matrix,y_matrix);

%size of grid point is (n+2) by n+2
xIndex = 2:N+1;
yIndex = 2:N+1;

u_boundary = zeros(N,N);
%left boundary values (keep y fixed and change the index of y)
u_boundary(:,1) = u0(x(xIndex,1), y(yIndex, 1));
u_boundary(:,N) = u0(x(xIndex,N+2), y(yIndex,N+2));
%upper and lower boundary 
u_boundary(1, :) = u_boundary(1,:) + u0(x(1, xIndex), y(1, yIndex));
u_boundary(N, :) = u_boundary(N, :) + u0(x(N+2, xIndex), y(N+2, yIndex));

u_boundary = reshape(u_boundary, N^2, 1);
u_boundary = u_boundary/(hj^2);

%evaluate f(x,y)

F = f(x(xIndex, yIndex), y(xIndex, yIndex)); 
F = reshape(F, N^2, 1);

I = speye(N);
%I = full(I);
e = ones(N, 1);
T = spdiags([e -4*e e],[-1 0 1], N, N);
%T = full(T);
S = spdiags([e e], [-1 1], N, N);
%S = full(S);
A = (kron(I, T) + kron(S, I)); 
%A = full(A);
A = A / (hj^2);

uApprox = A\(F - u_boundary);
 % for poisson 2DFD9 
 uApprox_5p = uApprox;

lower(1, 1:N) = [u0(x(1, xIndex), y(1, yIndex))];
upper(1, 1:N) = [u0(x(N+2, xIndex), y(N+2, yIndex))];
temp_left(1:N+2,1) = [u0(x(1:N+2, 1), y(1:N+2,1))];
temp_right(1:N+2,1) = [u0(x(1:N+2, N+2), y(1:N+2,N+2))];
% add boundary to uApprox 
uApprox = reshape(uApprox, N, N);
uApprox = [ lower; uApprox; upper ];
uApprox = [temp_left, uApprox, temp_right];
 %u2D = reshape(uApprox, size(x));
 surf(x,y, uApprox)
 

 
 
 
 