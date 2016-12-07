function [ x,y,  uApprox_5p] = poisson2DFD5_new(f, u0, L, N)

u_true = @(x,y) sin(pi.*x).*cos(pi.*y);

f = @(x,y) -2*(pi^2)*u_true(x,y);
u0 = u_true;
% compute the index of grid points except the boundary points 
function idx = compute_idx(i,j,N)
    idx = i*(N+2)+j;
end

% initializing the variable 
hj = L/(N+1);
x_matrix = 0:hj:L;
y_matrix = 0:hj:L;
% create mesgrid
[x,y] = meshgrid(x_matrix,y_matrix);

A = speye((N+2)^2,(N+2)^2);
F = sparse((N+2)^2,1);

% loop inside (only loop the points the we want to approximate i.e from
% 2:N+2, not including the boundary points 
for i = 2:N+1
    for j = 2:N+1
        idx = compute_idx(i,j,N); %index of point we approximate
        idxL= compute_idx(i-1,j,N); % left index boundary points
        idxR= compute_idx(i+1,j,N); % right index boundary points
        idxS= compute_idx(i,j-1,N); % south index boundary points 
        idxN= compute_idx(i,j+1,N); % north index boundary points 
        A(idx , idx ) = -4; % putting them in matrix 
        A(idx , idxL) =  1;
        A(idx , idxR) =  1;
        A(idx , idxS) =  1;
        A(idx , idxN) =  1;
        F(idx) = f(x(i,j),y(i,j))*hj*hj;
    end
end

% compute index of left & right boundary points
for j=1:N+2
   % index of left boundary points
   idx = compute_idx(1  ,j,N);
   A(idx,idx) = 1;
   F(idx) = u0(x(1,j),y(1,j));
   % index of right boundary points
   idx = compute_idx(N+2,j,N);
   A(idx,idx) = 1;
   F(idx) = u0(x(N+2,j),y(N+2,j));
end

% compute index of south and north boundary points
for i=2:N+1
    % index of south boundary points
   idx = compute_idx(i,1,N);
   A(idx,idx) = 1;
   F(idx) = u0(x(i,1),y(i,1));
    % index of north boundary points
   idx = compute_idx(i,N+2,N);
   A(idx,idx) = 1;
   F(idx) = u0(x(i,N+2),y(i,N+2));
end
 
uApprox = A\F;

% reshapre the matrix from (N+2)^2 by 1  to N+2 by N+2 matrix 
uApprox_5p = zeros(N+2,N+2);

for i=1:N+2
   for j=1:N+2
       idx = compute_idx(i,j,N);
       uApprox_5p(i,j) = uApprox(idx);
   end
end

figure()
    mesh(x,y, u_true(x,y))
    title('u_{true}')
    
    figure()
    mesh(x,y, uApprox_5p);
    title('u_{approx}');

end 
 
 
 