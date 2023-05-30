function [Q, R] = triangularizacion_ortogonal_householder(A)
% Funcion que triangulariza una matriz mxn utilizando la triangularizacion ortogonal de Householder
% A: Matriz mxn
% Q: Matriz ortogonal
% R: Matriz triangular superior

[m, n] = size(A);
Q = eye(m); % Inicializamos la matriz ortogonal
R = A; % Copiamos A en R

for j = 1:n
    norm_a = norm(R(j:m, j)); % Calculamos la norma de la columna j
    v = R(j:m, j); % Seleccionamos la columna j
    v(1) = v(1) + sign(v(1)) * norm_a; % Calculamos el primer elemento de v
    beta = (v' * v) / 2; % Calculamos beta
    w = v / (v(1) + norm_a); % Calculamos w
    R(j:m, j:n) = R(j:m, j:n) - 2 * w * (w' * R(j:m, j:n)); % Actualizamos R
    Q(j:m, :) = Q(j:m, :) - 2 * w * (w' * Q(j:m, :)); % Actualizamos Q
end

end
