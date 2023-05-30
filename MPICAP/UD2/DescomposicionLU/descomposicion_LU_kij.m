function [L,U] = descomposicion_LU_kij(A)
% descomposicion_LU_kij calcula la descomposición LU de una matriz A
% utilizando el algoritmo kij.
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
    for i = k+1:n
        L(i,k) = A(i,k)/A(k,k);
        for j = k+1:n
            A(i,j) = A(i,j) - L(i,k)*A(k,j);
        end
    end
end

U = A;
end