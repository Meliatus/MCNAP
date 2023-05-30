function det = det_LU(A)
% Funcion que calcula el determinante de una matriz A utilizando la descomposicion LU sin pivotamiento
% A: Matriz cuadrada de n x n
% det: determinante de la matriz A

[L,U] = lu(A); % Obtenemos la descomposición LU de A
det = prod(diag(U)); % Calculamos el determinante de U (ya que L es una matriz triangular inferior con diagonal 1)

end
