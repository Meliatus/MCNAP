function [x, nit, msg] = GaussSeidelCombined(A, b, x0, tol, maxit, nit_forward, nit_backward)

for k = 1:nit_forward
    [x, nit_f, msg_f] = GaussSeidel(A, b, x0, tol, maxit);
    if msg_f == 'Converged'
        nit = nit_f;
        msg = msg_f;
        break
    end
end

for k = 1:nit_backward
    [x, nit_b, msg_b] = GaussSeidelBackward(A, b, x0, tol, maxit);
    if msg_b == 'Converged'
        nit = nit_b;
        msg = msg_b;
        break
    end
end

if msg_f == 'Did not converge' && msg_b == 'Did not converge'
    nit = maxit;
    msg = 'Did not converge';
end

end
