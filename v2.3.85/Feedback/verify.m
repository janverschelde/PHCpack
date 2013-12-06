clear
verify_data % read the input from the output of PHC -> C programe
dim = n+q;  % n and q defined in "verify_data*"
for i=1:length(F)
    fprintf('\n******* Solution No. %d *********\n', i); 
    fprintf('F{%d}=%lf\n', i); disp(F{i})
    fprintf('G{%d}=%lf\n', i); disp(G{i})
    fprintf('H{%d}=%lf\n', i); disp(H{i})
    fprintf('K{%d}=%lf\n', i); disp(K{i})
    CL=[A+B*K{i}*C B*H{i}
          G{i}*C     F{i} ];
    [X,E1]=eig(CL);  % X is the right eigenvectors matrix of close loop
    [Y,E2]=eig(CL'); % Y is the left eigenvectors matrix of close loop
    
    for j=1:dim
        EE1(j)=E1(j,j);
        EE2(j)=E2(j,j);
    end 
  
    YY=v_perm(EE1,EE2, Y);  % Permute Y matrix to make it has the same order
                            % of eigenvalues as the X matrix 
    diff=v_pole(Poles, EE1, dim);% Compare given poles and computed poles
                          
    for k=1:dim   
      fprintf('Eigenvalue ')
      disp(EE1(k))
      fprintf('Difference: %e\n', diff(k))
      y = YY(:,k);
      evcond =1./abs(y'*X(:,k));
      fprintf('The condition number is: %.4e\n\n', evcond)
  end
end
