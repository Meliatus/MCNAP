function [x, nit, msg] = GaussSeidel(A, b, x0, tol, maxit)

n = length(b);
x = x0;

for k = 1:maxit
    x_old = x;
    for i = 1:n
        s = 0;
        for j = 1:i-1
            s = s + A(i,j)*x(j);
        end
        for j = i+1:n
            s = s + A(i,j)*x_old(j);
        end
        x(i) = (b(i) - s)/A(i,i);
    end
    
    if norm(x - x_old) < tol
        nit = k;
        msg = 'Converged';
        break
    end
end

if norm(x - x_old) >= tol
    nit = maxit;
    msg = 'Did not converge';
end

end
