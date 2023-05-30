function L = descomposicion_cholesky(A)
% descomposicion_cholesky calcula la descomposición de Cholesky de una matriz A.
%
% entrada:
% A - matriz a descomponer
%
% salida:
% L - matriz triangular inferior

[n,~] = size(A);
L = zeros(n,n);

for k = 1:n
    suma = 0;
    for p = 1:k-1
        suma = suma + L(k,p)^2;
    end
    L(k,k) = sqrt(A(k,k) - suma);
    for i = k+1:n
        suma = 0;
        for p = 1:k-1
            suma = suma + L(i,p)*L(k,p);
        end
        L(i,k) = (A(i,k) - suma) / L(k,k);
    end
end
end

