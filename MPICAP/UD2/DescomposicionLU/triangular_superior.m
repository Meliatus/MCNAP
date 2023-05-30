function x = triangular_superior(A, b)
% triangular_superior resuelve un sistema de ecuaciones lineales A*x = b
% donde la matriz A es triangular superior.
%
% entrada:
% A - matriz de coeficientes (triangular superior)
% b - vector constante
%
% salida:
% x - vector solución

[n,~] = size(A);
x = zeros(n,1);

% resolución hacia atrás
for i = n:-1:1
    x(i) = (b(i) - A(i,i+1:n)*x(i+1:n))/A(i,i);
end
end