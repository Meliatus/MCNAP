%-----------------------------------------------------------------------
%        A times a vector in Diagonal storage format (DIA) 
%----------------------------------------------------------------------- 
% multiplies a matrix by a vector when the original matrix is stored 
% in the diagonal storage format.
%-----------------------------------------------------------------------
%
% on entry:
%----------
% n     = row dimension of A
% x     = real array of length equal to the column dimension of
%         the A matrix.
% ndiag  = integer. The first dimension of array adiag as declared in
%         the calling program.
% idiag  = integer. The number of diagonals in the matrix.
% diag   = real array containing the diagonals stored of A.
% idiag  = number of diagonals in matrix.
% diag   = real array of size (ndiag x idiag) containing the diagonals
%          
% ioff   = integer array of length idiag, containing the offsets of the
%   	   diagonals of the matrix:
%          diag(i,k) contains the element a(i,i+ioff(k)) of the matrix.
%
% on return:
%-----------
% y     = real array of length n, containing the product y=A*x

function y = diag_filas(n, x, y, diag, ndiag, idiag, ioff)
    for j = 1:n
        y(j) = 0;
    end
    for j = 1:idiag
        io = ioff(j);
        i1 = max(1, 1-io);
        i2 = min(n, n-io);
        for k = i1:i2
            y(k) = y(k) + diag(k,j) * x(k + io);
        end
    end
end

