A=rand(m,r);
B=rand(r,n);
C=zeros(m,n);
for i=1:m
    for j=1:n
        for k=1:r
            C(i,j)= C(i,j)+A(i,k)*B(k,j);
        end
    end
end
