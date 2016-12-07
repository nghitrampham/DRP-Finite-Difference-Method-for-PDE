function [matrix_norm1] = poisson2_helper(L, N)
close all

matrix_norm1 =[];
error = [];
error_matrix = [];
for k = 1:5
n = N;
n =((2^k)/2)*n;
N_val1(k) = n;
[x_val1, u1 ,E1_2norm1] = poisson2(L,n);
matrix = E1_2norm1;
matrix_norm1 =[ matrix matrix_norm1];
end 

loglog(N_val1, matrix_norm1, '-s')