function y = matriz_vector_denso(A, x)
    [m, n] = size(A);
    y = zeros(m, 1);
    for i = 1:m
        for j = 1:n
            y(i) = y(i) + A(i,j) * x(j);
        end
    end
end