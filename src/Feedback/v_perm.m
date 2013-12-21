function [YY] =v_perm(A1, A2, Y)
% permute the left eigenvector matrix to make it have the same order
% of eigenvalues as the right eigenvector matrix
for j=1:length(A1)
   for k=1:length(A2)   
       if(abs(A1(j)-conj(A2(k)))<1e-5)
            YY(:,j) = Y(:,k);
       end
    end
end
