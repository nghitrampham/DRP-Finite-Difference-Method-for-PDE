function [ x,y,  uApprox_9p] = poisson2DFD9_new(f, u0, L, N)

% compute the index of grid points except the boundary points 
function idx = compute_idx(i,j,N)
    idx = (i-1)*(N+2)+j;
end

% initializing the variable 
hj = L/(N+1);
x_matrix = 0:hj:L;
y_matrix = 0:hj:L;
[x,y] = meshgrid(x_matrix,y_matrix);

A = speye((N+2)^2,(N+2)^2);
F = zeros((N+2)^2,1);

%compute deltaF
deltaF = compute_deltaF5(f, L,N);

% compute left & right boundary points 
for j=1:N+2
    % idx of left boundary points 
   idx = compute_idx(1,j,N);
   A(idx,idx) = 1;
   F(idx) =  u0(x(1,j),y(1,j));
   % idx of right boundary points
   idx = compute_idx(N+2,j,N);
   A(idx,idx) = 1;
   F(idx) =  u0(x(N+2,j),y(N+2,j));
end

% compute index of south and north boundary points
for i=1:N+2
    % index of south boundary points 
   idx = compute_idx(i,1,N);
   A(idx,idx) = 1;
   F(idx) =  u0(x(i,1),y(i,1));
   %index of north boundary points 
   idx = compute_idx(i,N+2,N);
   A(idx,idx) = 1;
   F(idx) = u0(x(i,N+2),y(i,N+2));
end

% loop inside (only loop the points the we want to approximate i.e from
% 2:N+2, not including the boundary points 
for i = 2:N+1
    for j = 2:N+1
        idx = compute_idx(i,j,N); %index of approximated point
        idxL= compute_idx(i-1,j,N); % left boundary index 
        idxR= compute_idx(i+1,j,N); % right boundary index 
        idxS= compute_idx(i,j-1,N); % south boundary index
        idxN= compute_idx(i,j+1,N); % north boundary index
        idx00 = compute_idx(i-1,j-1,N); % lower left boundary index
        idx01 = compute_idx(i-1,j+1,N); % upper left boundary index
        idx10 = compute_idx(i+1,j-1,N); % lower right boundary index
        idx11 = compute_idx(i+1,j+1,N); % upper right boundary index
        % putting them into matrix A 
        A(idx , idx ) = -20; 
        A(idx , idxL) =  4;
        A(idx , idxR) =  4;
        A(idx , idxS) =  4;
        A(idx , idxN) =  4;
        A(idx , idx00) =  1;
        A(idx , idx01) =  1;
        A(idx , idx10) =  1;
        A(idx , idx11) =  1;
        F(idx) = 6*hj*hj*(f(x(i,j),y(i,j))) + ((hj^4)*deltaF(idx))/2;
    end
end
 

uApprox = A\F;
uApprox_9p = reshape(uApprox, N+2, N+2);

disp('testing..');
% reshapre the matrix from (N+2)^2 by 1  to N+2 by N+2 matrix 
% uApprox_9p = zeros(N+2,N+2);
% 
% for i=1:N+2
%    for j=1:N+2
%        idx = compute_idx(i,j,N);
%        uApprox_9p(i,j) = uApprox(idx);
%    end
% end



end 
 
 
 