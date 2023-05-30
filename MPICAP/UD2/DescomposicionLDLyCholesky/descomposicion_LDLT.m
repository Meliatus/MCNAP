function [L,D] = descomposicion_LDLT(A)
% descomposicion_LDLT calcula la descomposición LDLT de una matriz A.
%
% entrada:
% A - matriz a descomponer
%
% salida:
% L - matriz triangular inferior
% D - matriz diagonal

[n,~] = size(A);
L = zeros(n,n);
D = zeros(n,n);

for k = 1:n
    suma = 0;
    for p = 1:k-1
        suma = suma + L(k,p)*D(p,p)*L(k,p);
    end
    D(k,k) = A(k,k) - suma;
    for i = k+1:n
        suma = 0;
        for p = 1:k-1
            suma = suma + L(i,p)*D(p,p)*L(k,p);
        end
        L(i,k) = (A(i,k) - suma) / D(k,k);
    end
end

L = L + eye(n);
end