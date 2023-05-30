function [Q,R]=hessenberg(A)
[m,n]=size(A);
Q=eye(m);
R=A;

for k=1:n-1
    x=R(k+1:m,k);
    v=sign(x(1))*norm(x)*eye(m-k,1)+x;
    v=v/norm(v);
    R(k+1:m,k:n)=R(k+1:m,k:n)-2*v*(v'*R(k+1:m,k:n));
    R(1:m,k+1:n)=R(1:m,k+1:n)-2*(R(1:m,k+1:n)*v)*v';
    Q(1:m,k+1:m)=Q(1:m,k+1:m)-2*Q(1:m,k+1:m)*v*v';
end
