function [H,v]=arnoldi_modificado(A,v,k)
n = size(A,1);
m = k+1;
Q = zeros(n,m);
H = zeros(m,m);

v_norm = norm(v);
if v_norm == 0
    return;
end
v(1) = v/v_norm;

for j=1:k
    w(j) = A*v(j);
    for i=1:j
        H(i,j) = w(j)'*v(i);
        w(j) = w(j) - H(i,j)*v(i);
    end
    H(j+1,j) = norm(w(j));
    if H(j+1,j) == 0
        break;
    end
    v(j+1) = w(j)/H(j+1,j);
end
