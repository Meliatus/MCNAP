%-----------------------------------------------------------------------
%        A times a vector in Jagged-Diagonal storage format (JAD) 
%----------------------------------------------------------------------- 
% multiplies a matrix by a vector when the original matrix is stored 
% in the jagged diagonal storage format.
%-----------------------------------------------------------------------
%
% on entry:
%----------
% n      = row dimension of A
% x      = real array of length equal to the column dimension of
%         the A matrix.
% jdiag  = integer. The number of jadded-diagonals in the data-structure.
% a      = real array containing the jadded diagonals of A stored
%          in succession (in decreasing lengths) 
% j      = integer array containing the colum indices of the 
%          corresponding elements in a.
% ia     = integer array containing the lengths of the  jagged diagonals
%
% on return:
%-----------
% y      = real array of length n, containing the product y=A*x

function y = jagged_diag(n, x, y, jdiag, a, ja, ia)
    y = zeros(1, n);
    for i=1:n
        y(i) = 0;
    end
    for ii=1:jdiag
        k1 = ia(ii) - 1;
        len = ia(ii + 1) - k1 - 1;
        for j=1:len
            y(j) = y(j) + a(k1 + j) * x(ja(k1 + j));
        end
    end
end
