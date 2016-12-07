function [error] = heat2D_error(Tf,N, L)

u_true = @(x,y,t) 
u0 = u_true;
uB = @(x,y,t) 0;
error =[];


for k = 1:5
    n = N;
    n = ((2^k)/2)*n;
    N_val(k) = n;
    
    [x ,y, uApproxTf] = heat2D(f,u0, uB, L, Tf, N);
    figure()
    mesh(X,Y, u_true(X,Y))
    title('u_{true}')
    
    figure()
    mesh(X,Y, uApproxTf);
    title('u_{approx}');
    
    E = uApproxTf' - u_true(X,Y);
    E_norm = norm(reshape(E, (n+2)^2, 1), Inf);
    error =[ error E_norm];
   
end
figure()
loglog(N_val, error, '-s'); hold on;


end 
