function C = matriz_producto_vectorial_columnas(A, B)
  [m, n] = size(A);
  [p, q] = size(B);
  if n ~= p
    error('Las dimensiones de las matrices no son compatibles para el producto matriz por matriz.');
  end
  C = zeros(m, q);
  for j = 1:q
    for i = 1:m
      row_i_A = A(i, :);
      col_j_B = B(:, j);
      C(i, j) = row_i_A * col_j_B;
    end
  end
end