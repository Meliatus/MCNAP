function A = descomposicion_LDLT_inplace(A)
% descomposicion_LDLT_inplace calcula la descomposición LDLT de una matriz A sobreescribiéndola sobre ella misma.
%
% entrada:
% A - matriz a descomponer
%
% salida:
% A - matriz descompuesta

[n,~] = size(A);

for k = 1:n
    for i = k+1:n
        A(i,k) = A(i,k) / A(k,k);
        for j = k+1:n
            A(i,j) = A(i,j) - A(i,k)*A(k,j);
        end
    end
end
end