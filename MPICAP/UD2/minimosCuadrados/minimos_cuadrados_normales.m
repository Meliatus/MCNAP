function x = minimos_cuadrados_normales(A, b)
% Funcion que resuelve el problema de minimos cuadrados utilizando las ecuaciones normales
% A: Matriz de m x n
% b: Vector de m x 1
% x: Vector solución de n x 1

x=A'*A\A'*b;


end
