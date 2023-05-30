function y = matriz_vector_disperso(S, IA, JA, x)
    n = length(x);
    y = zeros(n,1);
    for i = 1:n
        for j = IA(i):IA(i+1)-1
            y(i) = y(i) + S(j) * x(JA(j));
        end
    end
end
