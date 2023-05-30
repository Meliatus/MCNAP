function [Q,R]=tridiag_qr(T)
n = size(T,1);
Q = eye(n);
R = T;

for i=1:n-1
    x = R(i:i+1,i);
    v = sign(x(1))*norm(x)*[1;0] + x;
    v = v/norm(v);
    R(i:i+1,i:n) = R(i:i+1,i:n) - 2*v*(v'*R(i:i+1,i:n));
    Q(:,i:i+1) = Q(:,i:i+1) - 2*Q(:,i:i+1)*v*v';
end
