function [L,U] = descomposicion_LU_jki(A)
% descomposicion_LU_jki calcula la descomposición LU de una matriz A
% utilizando el algoritmo jki.
%
% entrada:
% A - matriz a descomponer
%
% salida:
% L - matriz triangular inferior
% U - matriz triangular superior

[n,~] = size(A);
L = eye(n); % matriz identidad

for k = 1:n-1
    for j = k+1:n
        for i = k+1:n
            A(i,j) = A(i,j) - A(i,k)*A(k,j)/A(k,k);
        end
    end
end

U = triu(A);
end