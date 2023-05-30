function x = resuelve_LU_con_pivot(A,b)
% resuelve_LU_con_pivot resuelve un sistema de ecuaciones Ax=b utilizando la descomposición LU con pivotamiento de A.
%
% entrada:
% A - matriz de coeficientes del sistema de ecuaciones
% b - vector de terminos independientes
%
% salida:
% x - vector solución

[n,~] = size(A);
L = eye(n,n);
U = zeros(n,n);
p = 1:n;

% descomposición LU
for k = 1:n-1
    [~,j_max] = max(abs(A(k:n,k)));
    j_max = j_max + k - 1;
    if j_max ~= k
        A([k j_max],:) = A([j_max k],:);
        b([k j_max]) = b([j_max k]);
        p([k j_max]) = p([j_max k]);
    end
    for j = k+1:n
        U(k,j) = A(k,j) / A(k,k);
        A(j,j:n) = A(j,j:n) - U(k,j) * A(k,j:n);
    end
end
U(n,n) = A(n,n);

% resolución del sistema
y = zeros(n,1);
x = zeros(n,1);

y(1) = b(p(1));
for i = 2:n
    suma = 0;
    for j = 1:i-1
        suma = suma + L(i,j) * y(j);
    end
    y(i) = b(p(i)) - suma;
end

x(n) = y(n) / U(n,n);
for i = n-1:-1:1
    suma = 0;
    for j = i+1:n
        suma = suma + U(i,j) * x(j);
    end
    x(i) = (y(i) - suma) / U(i,i);
end
end
