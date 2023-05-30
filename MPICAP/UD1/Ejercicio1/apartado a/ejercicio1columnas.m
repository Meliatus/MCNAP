A=rand(m,r);
B=rand(r,n);
C=zeros(m,n);
for j=1:n
    for i=1:m
        for k=1:r
            C(i,j)= C(i,j)+A(i,k)*B(k,j);
        end
    end
end
