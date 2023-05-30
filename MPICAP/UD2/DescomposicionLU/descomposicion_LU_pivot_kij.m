function [L,U,P] = descomposicion_LU_pivot_kij(A)
[n,~] = size(A);
L = eye(n); % matriz identidad
P = eye(n); % matriz identidad

for k = 1:n-1
    [~,m] = max(abs(A(k:n,k)));
    m = m + k - 1; % ajuste para el subvector
    if m ~= k
        A([k m],:) = A([m k],:); % intercambio de filas
        P([k m],:) = P([m k],:); % intercambio de filas
        if k >= 2
            L([k m],1:k-1) = L([m k],1:k-1); % intercambio de filas
        end
    end
    for i = k+1:n
        L(i,k) = A(i,k)/A(k,k);
        for j = k+1:n
            A(i,j) = A(i,j) - L(i,k)*A(k,j);
        end
    end
end

U = A;
end
