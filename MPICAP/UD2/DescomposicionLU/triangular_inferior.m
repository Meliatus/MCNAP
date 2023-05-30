function x = triangular_inferior(A, b)
% triangular_inferior resuelve un sistema de ecuaciones lineales A*x = b
% donde la matriz A es triangular inferior.
%
% entrada:
% A - matriz de coeficientes (triangular inferior)
% b - vector constante
%
% salida:
% x - vector solución

[n,~] = size(A);
x = zeros(n,1);

% resolución hacia adelante
for i = 1:n
    x(i) = (b(i) - A(i,1:i-1)*x(1:i-1))/A(i,i);
end
end