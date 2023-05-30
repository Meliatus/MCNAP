%
%*on entry:
%c----------
% n     = row dimension of A
% x     = real array of length equal to the column dimension of
%         the A matrix.
% na    = integer. The first dimension of arrays a and ja
%         as declared by the calling program.
% ncol  = integer. The number of active columns in array a.
%         (i.e., the number of generalized diagonals in matrix.)
% a, ja = the real and integer arrays of the itpack format
%         (a(i,k),k=1,ncol contains the elements of row i in matrix
%          ja(i,k),k=1,ncol contains their column numbers) 
%
% on return:
%-----------
% y     = real array of length n, containing the product y=y=A*x
%
%-----------------------------------------------------------------------
function y = ellpack_itpack_filas(n, x, na, ncol, a, ja)
    y = zeros(n, 1);
    for i = 1:n
        y(i) = 0;
    end
    for i = 1:n
        for j = 1:ncol
            y(i) = y(i) + a(i,j) * x(ja(i,j));
        end
    end
end