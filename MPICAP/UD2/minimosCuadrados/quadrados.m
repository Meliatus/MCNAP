function [x, R, norm_residual] = quadrados(A, b)
  [Q, R] = qr_householder(A);
  y = Q' * b;
  x = R \ y;
  norm_residual = norm(b - A * x);
end

function [Q, R] = qr_householder(A)
  [m, n] = size(A);
  Q = eye(m);
  R = A;
  for j = 1 : n
    x = R(j:m, j);
    v = x;
    v(1) = sign(x(1)) * norm(x) + x(1);
    v = v / norm(v);
    R(j:m, j:n) = R(j:m, j:n) - 2 * v * (v' * R(j:m, j:n));
    Q(:, j:m) = Q(:, j:m) - 2 * (Q(:, j:m) * v) * v';
  end
end
