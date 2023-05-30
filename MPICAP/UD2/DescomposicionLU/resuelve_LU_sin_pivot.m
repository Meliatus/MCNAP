function x = resuelve_LU_sin_pivot(A,b)
% resuelve_LU_sin_pivot resuelve un sistema de ecuaciones Ax=b utilizando la descomposición LU sin pivotamiento de A.
%
% entrada:
% A - matriz de coeficientes del sistema de ecuaciones
% b - vector de terminos independientes
%
% salida:
% x - vector solución

[n,~] = size(A);
L = zeros(n,n);
U = zeros(n,n);

% descomposición LU
for k = 1:n
    U(k,k) = A(k,k);
    L(k,k) = 1;
    for j = k+1:n
        L(j,k) = A(j,k) / U(k,k);
        U(k,j) = A(k,j);
    end
    for i = k+1:n
        for j = k+1:n
            A(i,j) = A(i,j) - L(i,k) * U(k,j);
        end
    end
end

% resolución del sistema
y = zeros(n,1);
x = zeros(n,1);

y(1) = b(1);
for i = 2:n
    suma = 0;
    for j = 1:i-1
        suma = suma + L(i,j) * y(j);
    end
    y(i) = b(i) - suma;
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
