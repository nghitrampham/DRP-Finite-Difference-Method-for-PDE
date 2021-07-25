function [x ,y, uApproxTf] = heat2D(f,u0, uB, L, Tf, N)

% compute the index of grid points except the boundary points 
function idx = compute_idx(i,j,N)
    idx = (i-1)*(N+2)+j;
end

% initializing the variable, spatial step size( hj = delta_x) 
hj = L/(N+1);
x_matrix = 0:hj:L;
y_matrix = 0:hj:L;
% create mesgrid
[x,y] = meshgrid(x_matrix,y_matrix);

uTM1 = zeros((N+2)^2,1);
uApproxTf = zeros((N+2)^2,1);

% calculate time step size ( k = delta_t)
k = ((hj)^2)/4;
time_step = 0:k:Tf;
coeff = k/((hj)^2);
% initial value at time t = 0 
t = 0; 

% compute index of left & right boundary points at time t = 0
for j=1:N+2
   % index of left boundary points
   idx = compute_idx(1,j,N);
   uTM1(idx) = u0(x(1,j),y(1,j));
   % index of right boundary points
   idx = compute_idx(N+2,j,N);
   uTM1(idx) = u0(x(N+2, j), y(N+2,j));
end

% compute index of south and north boundary points
for i=2:N+1
    % index of south boundary points
   idx = compute_idx(i,1,N);
   uTM1(idx) = u0(x(i,1), y(i,1));
    % index of north boundary points
   idx = compute_idx(i,N+2,N);
   uTM1(idx) = u0(x(i,N+2), y(i,N+2));
end

% loop inside (only loop the points the we want to approximate i.e from
% 2:N+2, not including the boundary points 
for i = 2:N+1
    for j = 2:N+1
        idx = compute_idx(i,j,N); %index of point we approximate
        uTM1(idx) = u0(x(i,j),y(i,j));
    end
end


%Loop for time 
for time = 2:length(time_step)
        %updating the boundary for the next time step
    % compute index of left & right boundary points at time t = 0
    for j=1:N+2
       % index of left boundary points
       idx = compute_idx(1,j,N);
       uApproxTf(idx) = uB(x(1,j),y(1,j),time_step(time));
       % index of right boundary points
       idx = compute_idx(N+2,j,N);
       uApproxTf(idx) = uB(x(N+2,j), y(N+2,j), time_step(time));
    end

    % compute index of south and north boundary points
    for i=2:N+1
        % index of south boundary points
       idx = compute_idx(i,1,N);
       uApproxTf(idx) = uB(x(i,1),y(i,1),time_step(time));
        % index of north boundary points
       idx = compute_idx(i,N+2,N);
       uApproxTf(idx) = uB(x(i, N+2),y(i,N+2), time_step(time));
    end
    for i = 2:N+1
        for j = 2:N+1
            idx = compute_idx(i,j,N); % index of point we approximate
            idxL= compute_idx(i-1,j,N); % left index boundary points
            idxR= compute_idx(i+1,j,N); % right index boundary points
            idxS= compute_idx(i,j-1,N); % south index boundary points 
            idxN= compute_idx(i,j+1,N);
         uApproxTf(idx) = uTM1(idx) + coeff*(uTM1(idxS) + uTM1(idxL) - 4*uTM1(idx) ...
                                + uTM1(idxR) + uTM1(idxN)) + k*f(x(i,j),y(i,j),time_step(time));
       end
    end
    
    uTM1 = uApproxTf;


end

uApproxTf = reshape(uApproxTf, N+2, N+2);

end
 
 