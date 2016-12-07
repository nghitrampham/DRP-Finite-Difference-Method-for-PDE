function [N_val matrix_norm] = poisson1_helper(L, a, b, N)
%call the function possion1 and return the matrix_norm which contains the
%norm2 error for each value of N 
close all

matrix_norm =[];
error = [];
error_matrix = [];
E = uApprox - u_true;

for k = 1:5
n = N;
n = ((2^k)/2)*n;
N_val(k) = n;
[x_val, u ,E1_2norm] = poisson1(L,a,b,n);
matrix = E1_2norm;
matrix_norm =[ matrix matrix_norm];
end 

loglog(N_val, matrix_norm, '-s')